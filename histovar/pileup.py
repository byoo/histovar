import pysam
import re

class PileupParser:
    def __init__(self, tumor_bam, normal_bam=None, include_indels=False, fasta=None):
        self.tumor_bam = tumor_bam
        self.normal_bam = normal_bam
        self.include_indels = include_indels
        self.fasta = fasta
        self.tumor_samfile = pysam.AlignmentFile(self.tumor_bam, "rb", reference_filename=self.fasta)
        self.normal_samfile = pysam.AlignmentFile(self.normal_bam, "rb", reference_filename=self.fasta) if self.normal_bam else None
        
    def close(self):
        if self.tumor_samfile:
            self.tumor_samfile.close()
        if self.normal_samfile:
            self.normal_samfile.close()

    def _get_bam_chrom(self, samfile, chrom):
        if not samfile or chrom in samfile.references:
            return chrom
        alt1 = chrom.replace('chr', '')
        if alt1 in samfile.references:
            return alt1
        alt2 = f'chr{chrom}'
        if alt2 in samfile.references:
            return alt2
        return chrom

    def _parse_xa_tag(self, read):
        """Parse BWA XA tag to find multi-mapped candidate loci."""
        if not read.has_tag('XA'):
            return []
        
        xa = read.get_tag('XA')
        # Format: chr,pos,CIGAR,NM;
        alt_alignments = xa.strip(';').split(';')
        loci = []
        for alt in alt_alignments:
            parts = alt.split(',')
            if len(parts) >= 4:
                loci.append({
                    "chrom": parts[0],
                    "pos": abs(int(parts[1])),
                    "cigar": parts[2],
                    "nm": int(parts[3])
                })
        return loci

    def extract_features(self, reads, ref_pos):
        """Extract features like depth, AF, MAPQ, strand bias for ML."""
        if not reads:
            return {"depth": 0, "alt_depth": 0, "af": 0.0, "mapq": 0, "strand_bias": 0.0, "pos_bias": 0.0}

        depth = len(reads)
        alt_depth = sum(1 for r in reads if r.get("is_alt", False))
        af = alt_depth / depth if depth > 0 else 0.0
        
        # Mean MAPQ
        mapq = sum(r.get("mapq", 0) for r in reads) / depth if depth > 0 else 0
        
        # Strand bias (forward vs reverse for ALT reads)
        alt_fwd = sum(1 for r in reads if r.get("is_alt", False) and not r.get("is_reverse", False))
        strand_bias = (alt_fwd / alt_depth) if alt_depth > 0 else 0.5
        
        # Position bias (relative position in read)
        pos_bias = sum(r.get("read_pos", 0) for r in reads if r.get("is_alt", False)) / alt_depth if alt_depth > 0 else 0.5
        
        return {
            "depth": depth,
            "alt_depth": alt_depth,
            "af": af,
            "mapq": mapq,
            "strand_bias": strand_bias,
            "pos_bias": pos_bias
        }

    def _read_pileup(self, samfile, chrom, pos, ref, alt):
        """Helper to extract read info from a bam file using pysam."""
        bam_chrom = self._get_bam_chrom(samfile, chrom)
        reads = []
        is_indel_target = '+' in alt or '-' in alt
        for pileupcolumn in samfile.pileup(bam_chrom, pos-1, pos, truncate=True, multiple_iterators=False):
                if pileupcolumn.pos == pos - 1:
                    for pileupread in pileupcolumn.pileups:
                        if pileupread.is_del or pileupread.is_refskip:
                            continue
                        aln = pileupread.alignment
                        if aln.query_sequence is None:
                            continue
                        
                        base = aln.query_sequence[pileupread.query_position]
                        
                        is_alt = False
                        if is_indel_target:
                            if pileupread.indel != 0:
                                if pileupread.indel > 0:
                                    ins_seq = aln.query_sequence[pileupread.query_position + 1 : pileupread.query_position + 1 + pileupread.indel]
                                    read_allele = f"{base}+{ins_seq}"
                                else:
                                    read_allele = f"{base}-{-pileupread.indel}"
                                if read_allele == alt:
                                    is_alt = True
                        else:
                            is_alt = (base == alt and pileupread.indel == 0)

                        reads.append({
                            "is_alt": is_alt,
                            "mapq": aln.mapping_quality,
                            "is_reverse": aln.is_reverse,
                            "read_pos": pileupread.query_position / max(1, len(aln.query_sequence)),
                            "read_id": aln.query_name,
                            "xa_loci": self._parse_xa_tag(aln)
                        })
        return reads

    def fetch_pileup(self, chrom, pos, ref, alt):
        """Fetch pileup from BAM."""
        tumor_reads = self._read_pileup(self.tumor_samfile, chrom, pos, ref, alt)
        tumor_features = self.extract_features(tumor_reads, pos)
        
        normal_features = {"depth": 0, "alt_depth": 0, "af": 0.0, "mapq": 0, "strand_bias": 0.5, "pos_bias": 0.5}
        if self.normal_samfile:
            normal_reads = self._read_pileup(self.normal_samfile, chrom, pos, ref, alt)
            normal_features = self.extract_features(normal_reads, pos)
            somatic_score = (tumor_features["af"] - normal_features["af"]) * tumor_features["depth"]
        else:
            somatic_score = tumor_features["af"] * tumor_features["depth"]
            
        return {
            "tumor": tumor_features,
            "normal": normal_features,
            "somatic_score": somatic_score,
            "tumor_reads": tumor_reads
        }

    def discover_variants(self, chrom, pos, min_alt_reads=2):
        """Scans pileup at a coordinate and identifies reference and alternate alleles."""
        bam_chrom = self._get_bam_chrom(self.tumor_samfile, chrom)
        
        # ULTRA-FAST PATH: Instantly skip coordinates that are strictly homozygous wild-type
        try:
            # count_coverage runs entirely in C (HTSlib) and is exponentially faster than pileup iterations.
            A, C, G, T = self.tumor_samfile.count_coverage(bam_chrom, pos-1, pos, quality_threshold=20)
            base_counts = [A[0], C[0], G[0], T[0]]
            base_counts.sort(reverse=True)
            
            # If the second highest base has fewer than min_alt_reads, it's virtually 100% homozygous.
            # If we aren't specifically looking for indels, we can skip the python pileup entirely!
            if not self.include_indels and base_counts[1] < min_alt_reads:
                return None, []
                
            # Fallback for indels or extremely low coverage
            if sum(base_counts) < min_alt_reads:
                return None, []
        except Exception as e:
            # Catch cases where chrom might not exist in the bam header
            print(f"Warning: skipped {chrom}:{pos} -> {e}")
            return None, []
            
        counts = {}
        for pileupcolumn in self.tumor_samfile.pileup(bam_chrom, pos-1, pos, truncate=True, multiple_iterators=False):
                if pileupcolumn.pos == pos - 1:
                    for pileupread in pileupcolumn.pileups:
                        if pileupread.is_del or pileupread.is_refskip:
                            continue
                        aln = pileupread.alignment
                        if aln.query_sequence is None:
                            continue
                        base = aln.query_sequence[pileupread.query_position]
                        counts[base] = counts.get(base, 0) + 1
                        
                        if self.include_indels and pileupread.indel != 0:
                            if pileupread.indel > 0:
                                ins_seq = aln.query_sequence[pileupread.query_position + 1 : pileupread.query_position + 1 + pileupread.indel]
                                indel_allele = f"{base}+{ins_seq}"
                            else:
                                indel_allele = f"{base}-{-pileupread.indel}"
                            counts[indel_allele] = counts.get(indel_allele, 0) + 1
        
        if not counts:
            return None, []
            
        # Ref is assumed to be the most frequent base (filtering out indels)
        base_counts = {k: v for k, v in counts.items() if '+' not in k and '-' not in k}
        if not base_counts:
            return None, []
            
        sorted_bases = sorted(base_counts.items(), key=lambda x: x[1], reverse=True)
        ref_allele = sorted_bases[0][0]
        
        alt_alleles = []
        for allele, count in counts.items():
            if allele != ref_allele and count >= min_alt_reads:
                alt_alleles.append(allele)
                
        return ref_allele, alt_alleles
