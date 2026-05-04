import datetime

class VCFWriter:
    def __init__(self, path):
        self.path = path
        self.file = open(path, 'w')
        self._write_header()
        
    def _write_header(self):
        date_str = datetime.datetime.now().strftime("%Y%m%d")
        self.file.write('##fileformat=VCFv4.2\n')
        self.file.write(f'##fileDate={date_str}\n')
        self.file.write('##source=Histovar_v1.2\n')
        self.file.write('##INFO=<ID=GENES,Number=.,Type=String,Description="Assigned paralog loci">\n')
        self.file.write('##INFO=<ID=AA,Number=1,Type=String,Description="Amino acid change (e.g. K27M)">\n')
        self.file.write('##INFO=<ID=POST,Number=1,Type=Float,Description="Posterior confidence from EM">\n')
        self.file.write('##INFO=<ID=ML,Number=1,Type=Float,Description="Machine Learning artifact probability score">\n')
        self.file.write('##INFO=<ID=TUMOR_AF,Number=1,Type=Float,Description="Tumor allele fraction">\n')
        self.file.write('##INFO=<ID=NORMAL_AF,Number=1,Type=Float,Description="Normal allele fraction">\n')
        self.file.write('##INFO=<ID=SOMATIC_SCORE,Number=1,Type=Float,Description="Somatic score">\n')
        self.file.write('##INFO=<ID=EVENT,Number=1,Type=String,Description="Event ID linking ambiguous multi-mapped records">\n')
        self.file.write('##INFO=<ID=MATE_POS,Number=1,Type=String,Description="Coordinate of the alternate paralog record">\n')
        self.file.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
        self.file.write('##FILTER=<ID=LowQual,Description="Low ML score or Posterior confidence">\n')
        self.file.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')

    def write_record(self, chrom, pos, ref, alt, info_dict, filter_status="PASS"):
        info_strs = []
        for k, v in info_dict.items():
            info_strs.append(f"{k}={v}")
        info_col = ";" .join(info_strs)
        
        row = f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t{filter_status}\t{info_col}\n"
        self.file.write(row)
        
    def close(self):
        self.file.close()

