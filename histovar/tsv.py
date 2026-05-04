class TSVWriter:
    def __init__(self, path):
        self.path = path
        self.file = open(path, 'w')
        self.headers = ["CHROM", "POS", "REF", "ALT", "FILTER", "AA", "GENES", "POST", "ML", "TUMOR_AF", "NORMAL_AF", "SOMATIC_SCORE", "EVENT", "MATE_POS"]
        self.file.write("\t".join(self.headers) + "\n")
        
    def write_record(self, chrom, pos, ref, alt, info_dict, filter_status="PASS"):
        row = [
            str(chrom),
            str(pos),
            str(ref),
            str(alt),
            str(filter_status),
            str(info_dict.get("AA", ".")),
            str(info_dict.get("GENES", ".")),
            str(info_dict.get("POST", ".")),
            str(info_dict.get("ML", ".")),
            str(info_dict.get("TUMOR_AF", ".")),
            str(info_dict.get("NORMAL_AF", ".")),
            str(info_dict.get("SOMATIC_SCORE", ".")),
            str(info_dict.get("EVENT", ".")),
            str(info_dict.get("MATE_POS", "."))
        ]
        self.file.write("\t".join(row) + "\n")
        
    def close(self):
        self.file.close()
