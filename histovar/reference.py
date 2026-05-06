import gzip
import re

H3_GENES = {
    "H3-3A": ["H3F3A", "H3.3A"],
    "H3-3B": ["H3F3B", "H3.3B"],
    "H3C1": ["HIST1H3A"],
    "H3C2": ["HIST1H3B"],
    "H3C3": ["HIST1H3C"],
    "H3C4": ["HIST1H3D"],
    "H3C6": ["HIST1H3E"],
    "H3C7": ["HIST1H3F"],
    "H3C8": ["HIST1H3G"],
    "H3C10": ["HIST1H3H"],
    "H3C11": ["HIST1H3I"],
    "H3C12": ["HIST1H3J"],
    "H3C13": ["HIST2H3C"],
    "H3C14": ["HIST2H3A"],
    "H3C15": ["HIST2H3D"],
}

CODON_TABLE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                 
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}

def translate_codon(codon):
    """Translate a 3-base DNA sequence to an amino acid."""
    if len(codon) != 3:
        return 'X'
    return CODON_TABLE.get(codon.upper(), 'X')

def reverse_complement(seq):
    mapping = str.maketrans('ATCGN', 'TAGCN')
    return seq.upper().translate(mapping)[::-1]

HOTSPOTS = {
    # H3 Hotspots (Note: DNA pos is +1 from legacy nomenclature due to initiator Methionine)
    "R26": {"amino_acid": "R", "pos": 27, "variants": []}, 
    "R27": {"amino_acid": "R", "pos": 27, "variants": []}, # New nomenclature for R26
    
    "K27": {"amino_acid": "K", "pos": 28, "variants": ["M", "I"]},
    "K28": {"amino_acid": "K", "pos": 28, "variants": ["M", "I"]}, # New nomenclature for K27
    
    "S28": {"amino_acid": "S", "pos": 29, "variants": ["C"]},
    "S29": {"amino_acid": "S", "pos": 29, "variants": ["C"]}, # New nomenclature for S28
    
    "G34": {"amino_acid": "G", "pos": 35, "variants": ["R", "V", "W"]},
    "G35": {"amino_acid": "G", "pos": 35, "variants": ["R", "V", "W"]}, # New nomenclature for G34
    
    "K36": {"amino_acid": "K", "pos": 37, "variants": ["M"]},
    "K37": {"amino_acid": "K", "pos": 37, "variants": ["M"]}, # New nomenclature for K36
    
    # H2A Hotspots (e.g. Carcinosarcomas)
    "E92": {"amino_acid": "E", "pos": 92, "variants": ["K"]},
    "E121": {"amino_acid": "E", "pos": 121, "variants": ["Q"]},
    
    # H2B Hotspots (e.g. Carcinosarcomas)
    "G53": {"amino_acid": "G", "pos": 53, "variants": ["D"]},
    "E76": {"amino_acid": "E", "pos": 76, "variants": ["K", "Q"]},
    "E113": {"amino_acid": "E", "pos": 113, "variants": ["K"]},
    
    # H1 Hotspots (Diffuse large B-cell lymphomas)
    # Generic placeholders - specific variants (like S102) can be added when specified
    "H1_GENERIC": {"amino_acid": "X", "pos": 0, "variants": []}
}

class ReferenceBuilder:
    def __init__(self, build="hg38"):
        self.build = build
        self.transcripts = {}

    def is_histone(self, gene_name, whitelist=None):
        """Check if a gene name follows histone nomenclature (H1, H2A, H2B, H3, H4, HIST...) or is in whitelist."""
        # 0. Check explicit whitelist
        if whitelist and gene_name in whitelist:
            return True

        # 1. Check curated H3 map
        for canonical, aliases in H3_GENES.items():
            if gene_name == canonical or gene_name in aliases:
                return True
        
        # 2. Pattern check (H1-H4, HIST, etc.)
        # Matches H1, H2A, H2B, H3, H4, or HIST prefixes
        if re.search(r'^(H[1-4][ABCDEF]?|HIST[1-5]H[1-4]|H[1-4]C[0-9])', gene_name, re.IGNORECASE):
            return True
        return False

    def get_canonical_name(self, gene_name):
        for canonical, aliases in H3_GENES.items():
            if gene_name == canonical or gene_name in aliases:
                return canonical
        return gene_name

    def parse_gtf(self, gtf_path, whitelist=None):
        """Parse GTF to extract exon coordinates for Histone genes."""
        import os, gzip, re
        if not gtf_path or not os.path.exists(gtf_path):
            raise FileNotFoundError(f"GTF file not found: {gtf_path}")

        open_fn = gzip.open if str(gtf_path).endswith('.gz') else open
        cds_blocks = {}
        
        with open_fn(gtf_path, 'rt') as f:
            for line in f:
                if line.startswith('#'): continue
                parts = line.strip().split('\t')
                if len(parts) < 9 or parts[2] != 'CDS': continue
                
                chrom = parts[0]
                if chrom.startswith('NC_0000'):
                    try:
                        num = int(chrom.split('.')[0].split('_')[1])
                        if num == 23: chrom = 'chrX'
                        elif num == 24: chrom = 'chrY'
                        else: chrom = f'chr{num}'
                    except: pass
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                attr_str = parts[8]
                
                m = re.search(r'gene_name "([^"]+)"', attr_str)
                if not m:
                    m = re.search(r'gene "([^"]+)"', attr_str)
                if not m: continue
                gene_name = m.group(1)
                
                if not self.is_histone(gene_name, whitelist=whitelist):
                    continue
                    
                canonical = self.get_canonical_name(gene_name)
                
                if canonical not in cds_blocks:
                    cds_blocks[canonical] = {
                        "chrom": chrom,
                        "strand": strand,
                        "cds_starts": [],
                        "cds_ends": []
                    }
                
                cds_blocks[canonical]["cds_starts"].append(start)
                cds_blocks[canonical]["cds_ends"].append(end)
                
        # Sort the CDS blocks to be strictly left-to-right relative to the genome
        for gene, data in cds_blocks.items():
            starts = sorted(data["cds_starts"])
            ends = sorted(data["cds_ends"])
            self.transcripts[gene] = {
                "chrom": data["chrom"],
                "strand": data["strand"],
                "cds_starts": starts,
                "cds_ends": ends
            }
    def get_genomic_pos_from_cds_offset(self, tx, cds_offset):
        """Map a CDS 0-based offset to a genomic coordinate."""
        strand = tx["strand"]
        starts = tx["cds_starts"]
        ends = tx["cds_ends"]
        
        current_offset = 0
        
        if strand == "+":
            for start, end in zip(starts, ends):
                exon_len = end - start + 1
                if current_offset + exon_len > cds_offset:
                    return start + (cds_offset - current_offset)
                current_offset += exon_len
        else:
            # For reverse strand, we process exons from last to first
            # (assuming they are sorted by genomic coordinate)
            for start, end in reversed(list(zip(starts, ends))):
                exon_len = end - start + 1
                if current_offset + exon_len > cds_offset:
                    return end - (cds_offset - current_offset)
                current_offset += exon_len
                
        return None # Out of bounds

    def get_hotspot_coordinates(self, gene_name, hotspot_id):
        """Map protein hotspot (e.g. K27) to genomic coordinates."""
        gene = self.get_canonical_name(gene_name)
        if gene not in self.transcripts:
            return None
        
        tx = self.transcripts[gene]
        base_offset = (HOTSPOTS[hotspot_id]["pos"] - 1) * 3
        
        codon_coords = []
        for i in range(3):
            pos = self.get_genomic_pos_from_cds_offset(tx, base_offset + i)
            if pos is not None:
                codon_coords.append(pos)
                
        if not codon_coords:
            return None
            
        return {
            "chrom": tx["chrom"], 
            "pos": codon_coords[0], # First base of codon (for compatibility)
            "codon_coords": codon_coords,
            "strand": tx["strand"]
        }
        
    def get_gene_from_pos(self, chrom, pos, pad=500):
        """Map a raw genomic coordinate back to a gene name based on transcript bounds, with padding for read lengths."""
        for gene, tx in self.transcripts.items():
            if tx["chrom"] == chrom:
                start = tx["cds_starts"][0] - pad
                end = tx["cds_ends"][-1] + pad
                if start <= pos <= end:
                    return gene
        return None
        
    def get_all_hotspot_coordinates(self):
        """Returns a list of all potential hotspot genomic coordinates to scan."""
        candidates = []
        seen = set()
        for gene in self.transcripts.keys():
            for hs_id, hs_info in HOTSPOTS.items():
                coords = self.get_hotspot_coordinates(gene, hs_id)
                if coords:
                    key = f'{gene}:{coords["pos"]}'
                    if key not in seen:
                        seen.add(key)
                        candidates.append({
                            "chrom": coords["chrom"],
                            "pos": coords["pos"],
                            "codon_coords": coords["codon_coords"],
                            "strand": coords["strand"],
                            "gene": gene,
                            "hotspot_id": hs_id,
                            "expected_aa": f"{hs_info['amino_acid']}{hs_info['pos']}"
                        })
        return candidates
