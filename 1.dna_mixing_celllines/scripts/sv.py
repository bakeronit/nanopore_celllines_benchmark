import re, sys
from dataclasses import dataclass
from datetime import datetime
from collections import namedtuple
from typing import List, Any, IO
from cyvcf2 import VCF
import pysam


@dataclass
class SV:
    bnd_pool = dict()
    BP = namedtuple("BP",['chrom','pos','ref', 'strand'])
    Record = namedtuple("Record", ["chrom1","pos","id","ref","alt","sv_type","qual","filter","chrom2","end","sv_len","format","sample_col"])
    output: IO[Any] = sys.stdout

    @classmethod
    def readSv(cls, filename):
        return VCF(filename)

    def write_header(self, header):
        for line in header.split('\n'):
            if not line.startswith('##INFO') and not line == '':
                if line.startswith('#CHROM'):
                    print('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">', file=self.output)
                    print('##INFO=<ID=CHR2,Number=1,Type=String,Description="Mate chromsome for BND SVs">', file=self.output)
                    print('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of structural variation">', file=self.output)
                    print('##INFO=<ID=SVLEN,Number=1,Type=Float,Description="Length of the SV">', file=self.output)
                    print('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR', file=self.output)
                else:
                    print(line, file=self.output)

    def write_record(self, rec):
        print(f"{rec.chrom1}\t{rec.pos}\t{rec.id}\t{rec.ref}\t{rec.alt}\t{rec.qual}\t{rec.filter}\tSVTYPE={rec.sv_type};CHR2={rec.chrom2};END={rec.end};SVLEN={rec.sv_len}\t{":".join(rec.format)}\t{rec.sample_col}", file=self.output)
    
    def get_simple_svtype(self, variant, bp):
        if variant.CHROM != bp.chrom:
            sv_type = "TRA"
        elif bp.strand == "+-":
            sv_type = "DEL"
        elif bp.strand == "-+":
            sv_type = "DUP"
        else:
            sv_type = "INV"

        return sv_type
    

class SAVANA(SV):
        
    def parse_variant(self, variant):
        sv_id = re.match(r"(.+)_[1,2]$", variant.ID).group(1)
        sv_len = variant.INFO.get("SVLEN")
        sv_type = variant.INFO.get('SVTYPE')
        qual = '.' if variant.QUAL is None else variant.QUAL
        mate_id = variant.INFO.get('MATEID')
        sample_col = str(variant).rstrip('\n').split("\t")[-1]
        if sv_type == "INS":
            rec = self.Record(variant.CHROM,variant.POS, sv_id, variant.REF, f'<{sv_type}>', sv_type, qual, variant.FILTERS[0], variant.CHROM, variant.POS, sv_len, variant.FORMAT, sample_col)
            self.write_record(rec)
        elif mate_id is None:
            print(f"WARNING: single BND {id} observed in SAVANA result")
        else:
            if not sv_id in self.bnd_pool:
                self.bnd_pool[sv_id] = self.BP(variant.CHROM, variant.POS, variant.REF, variant.INFO.get('BP_NOTATION'))
            else:
                bp1 = self.bnd_pool[sv_id]
                sv_type = self.get_simple_svtype(variant, bp1)
                rec = self.Record(bp1.chrom, bp1.pos, sv_id, bp1.ref, f'<{sv_type}>', sv_type, qual, variant.FILTERS[0], variant.CHROM, variant.POS, sv_len, variant.FORMAT, sample_col)
                self.write_record(rec)


class SEVERUS(SV):

    def parse_variant_old(self, variant):  ## severus v0.1.2
        sv_id = variant.ID
        sv_len = variant.INFO.get('SVLEN')
        sv_type = variant.INFO['SVTYPE']
        qual = '.' if variant.QUAL is None else variant.QUAL
        sample_col = str(variant).rstrip('\n').split("\t")[-1]
        if sv_type != "BND":
            rec = self.Record(variant.CHROM, variant.POS, sv_id, variant.REF, variant.ALT[0], sv_type, qual, variant.FILTERS[0], variant.INFO['CHR2'], variant.INFO['END'], sv_len, variant.FORMAT, sample_col )
            self.write_record(rec)
        else:
            if sv_id in self.bnd_pool:
                pass
            else:
                self.bnd_pool[sv_id] = 1
                if variant.CHROM != variant.INFO['CHR2']:
                    sv_type = "TRA"
                elif variant.INFO['STRANDS'] == "-+":
                    sv_type = "DUP"
                elif variant.INFO['STRANDS'] == "+-":
                    sv_type = "DEL"
                else:
                    sv_type = "INV"
                rec = self.Record(variant.CHROM, variant.POS, sv_id, variant.REF, f'<{sv_type}>', sv_type, qual, variant.FILTERS[0], variant.INFO['CHR2'], variant.INFO['END'], sv_len, variant.FORMAT, sample_col )
                self.write_record(rec)
    
    def parse_variant(self, variant):  ## severus v1.0.0 
        sv_id = variant.ID
        sv_len = variant.INFO.get('SVLEN')
        sv_type = variant.INFO['SVTYPE']
        qual = '.' if variant.QUAL is None else variant.QUAL
        sample_col = str(variant).rstrip('\n').split("\t")[-1]
        if sv_type != "BND":
            sv_len = int(sv_len)
            pos2 = variant.POS if sv_type=="INS" else variant.POS + sv_len 
            rec = self.Record(variant.CHROM, variant.POS, sv_id, variant.REF, variant.ALT[0], sv_type, qual, variant.FILTERS[0], variant.CHROM, pos2, sv_len, variant.FORMAT, sample_col )
            self.write_record(rec)
        else:
            sv_id = re.match(r"(.+)_[1,2]$", variant.ID).group(1)
            if not sv_id in self.bnd_pool:
                self.bnd_pool[sv_id] = self.BP(variant.CHROM, variant.POS, variant.REF, variant.INFO.get('STRANDS'))
            else:
                bp1 = self.bnd_pool[sv_id]
                sv_type = self.get_simple_svtype(variant, bp1)
                sv_len = 0 if sv_type == 'TRA' else sv_len
                rec = self.Record(bp1.chrom, bp1.pos, sv_id, bp1.ref, f'<{sv_type}>', sv_type, qual, variant.FILTERS[0], variant.CHROM, variant.POS, sv_len, variant.FORMAT, sample_col )
                self.write_record(rec)

class DELLY(SV):
    """
    Need to make END equals to POS2 for BND and fix header. so jasmine can read it.
    """
    def parse_variant(self, variant):
        sv_type = variant.INFO['SVTYPE']
        qual = '.' if variant.QUAL is None else variant.QUAL
        sample_col = str(variant).rstrip('\n').split("\t")[-2]
        if sv_type != "BND":
            sv_len = variant.INFO['SVLEN'] if sv_type == 'INS' else variant.INFO['END'] - variant.POS
            rec = self.Record(variant.CHROM, variant.POS, variant.ID, variant.REF, variant.ALT[0], sv_type, qual, variant.FILTERS[0], variant.CHROM, variant.INFO['END'], sv_len, variant.FORMAT, sample_col )
        else:
            sv_type = "TRA"
            rec = self.Record(variant.CHROM, variant.POS, variant.ID, variant.REF, variant.ALT[0], sv_type, qual, variant.FILTERS[0], variant.INFO['CHR2'], variant.INFO['POS2'], 0, variant.FORMAT, sample_col )
        self.write_record(rec)

@dataclass
class NANOMONSV(SV):
    reference: str = ''

    def __post_init__(self):
        self.ref_tb = pysam.FastaFile(self.reference)

    @classmethod
    def readSv(cls, filename):
        with open(filename, 'rt') as fh:
            for line in fh:
                line = line.strip()
                if line.startswith("Chr_1"):
                    header = line.split("\t")
                else:
                    yield({index:element for index, element in zip(header, line.split("\t"))})

    def write_header(self):
        header = '##fileformat=VCFv4.2\n'\
                 f'##fileDate={datetime.today().strftime("%Y%m%d")}\n'\
                 f'##source=nanomonsv-0.7.1\n'\
                 f'##reference={self.reference}\n'

        for (tchr, tlen) in zip(self.ref_tb.references, self.ref_tb.lengths):
            header += f'##config=<ID={tchr},length={tlen}>\n'
        
        header += '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n'\
                  '##INFO=<ID=CHR2,Number=1,Type=String,Description="Mate chromsome for BND SVs">\n'\
                  '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of structural variation">\n'\
                  '##INFO=<ID=SVLEN,Number=1,Type=Float,Description="Length of the SV">\n'\
                  '##FORMAT=<ID=TR,Number=1,Type=Integer,Description="The number of reads around the breakpoints">\n'\
                  '##FORMAT=<ID=VR,Number=1,Type=Integer,Description="The number of variant supporting reads determined in the validation realignment step">\n'\
                  '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR'
        
        print(header, file=self.output)

    def parse_variant(self, variant):
        match variant['SV_Type']:
            case "Deletion":
                sv_type = "DEL"
            case "Duplication":
                sv_type = "DUP"
            case "Inversion":
                sv_type = "INV"
            case "Insertion":
                sv_type = "INS"
            case "Translocation":
                sv_type = "TRA"
        
        sample_col = f"{variant["Checked_Read_Num_Tumor"]}:{variant["Supporting_Read_Num_Tumor"]}"
        format_col = ['TR','VR']
        chrom1 = variant['Chr_1']
        pos1 = int(variant['Pos_1'])
        ref = self.ref_tb.fetch(chrom1, pos1 - 1, pos1)
        alt = variant['Inserted_Seq'] if sv_type == "INS" else f'<{sv_type}>'
        end = int(variant['Pos_2']) - 1
        sv_len = len(variant['Inserted_Seq']) if sv_type== "INS" else end - pos1 
        rec = self.Record(chrom1, pos1, variant['SV_ID'], ref, alt, sv_type, '.', variant['Is_Filter'], variant['Chr_2'], end, sv_len, format_col, sample_col )

        self.write_record(rec)

@dataclass
class QSV(SV):
    reference: str = ''

    def __post_init__(self):
        self.ref_tb = pysam.FastaFile(self.reference)

    @classmethod
    def readSv(cls, filename):
        with open(filename, 'rt') as fh:
            for line in fh:
                line = line.strip()
                if line.startswith("#"):
                    continue
                elif line.startswith("analysis_id"):
                    header = line.split("\t")
                else:
                    yield({index:element for index, element in zip(header, line.split("\t"))})

    def write_header(self):
        header = '##fileformat=VCFv4.2\n'\
                 f'##fileDate={datetime.today().strftime("%Y%m%d")}\n'\
                 f'##source=qsv\n'\
                 f'##reference={self.reference}\n'

        for (tchr, tlen) in zip(self.ref_tb.references, self.ref_tb.lengths):
            header += f'##config=<ID={tchr},length={tlen}>\n'
        
        header += '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n'\
                  '##INFO=<ID=CHR2,Number=1,Type=String,Description="Mate chromsome for BND SVs">\n'\
                  '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of structural variation">\n'\
                  '##INFO=<ID=SVLEN,Number=1,Type=Float,Description="Length of the SV">\n'\
                  '##FORMAT=<ID=TUMOUR_DP,Number=1,Type=Integer,Description="The number of reads in tumour">\n'\
                  '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR'
        
        print(header, file=self.output)
    
    def parse_variant(self, variant):
        chrom1 = 'chr' + variant['chr_from'].replace('chr','')
        chrom2 = 'chr' + variant['chr_to'].replace('chr','')
        chrom1 = 'chrX' if chrom1 == 'chr23' else chrom1
        chrom2 = 'chrX' if chrom2 == 'chr23' else chrom2
        chrom1 = 'chrY' if chrom1 == 'chr24' else chrom1
        chrom2 = 'chrY' if chrom2 == 'chr24' else chrom2
        
        if not variant['orientation_category'] == "2":
            pos1 = int(variant['chr_from_bkpt'])
            pos2 = int(variant['chr_to_bkpt'])
        else:
            pos2 = int(variant['chr_from_bkpt'])
            pos1 = int(variant['chr_to_bkpt'])
        
        sv_id = variant['sv_id']
        ref = "N"
        sv_len = pos2 - pos1
        match variant['annotation']:
            case "CTX":
                sv_type = "TRA"
                sv_len = 0
            case "DEL/ITX":
                sv_type = "DEL"
            case "INV/ITX":
                sv_type = "INV"
            case "DUP/INS/ITX":
                sv_type = "DUP"
            case _:
                sv_type = "NA"

        alt = f"<{sv_type}>"
        format_col = ["TUMOUR_DP"]
        sample_col = f"{variant['number_of_reads']}"

        #if variant['category'] == '1' or variant['category'] == '2':
        if True: # test using all events
            rec = self.Record(chrom1, pos1, sv_id, ref, alt, sv_type, '.', 'PASS', chrom2, pos2, sv_len, format_col, sample_col )
            self.write_record(rec)
        

class GRIDSS(SV):

    def parse_variant(self, variant):
        sv_id = variant.INFO.get("EVENT")
        sv_type = variant.INFO.get("EVENTTYPE")
        sv_type = 'TRA' if sv_type == "BND" else sv_type
        if sv_type != "SGL":
            #match_alt = re.match(r".*[\[,\]](chr(?:[0-9]+|X|Y)):(\d+).+", variant.ALT[0])
            match_alt = re.match(r".*[\[,\]](chr(?:[^:]+)):(\d+).+", variant.ALT[0])
            chrom2 = match_alt.group(1)
            pos2 = int(match_alt.group(2))
        qual = "."
        sample_col = str(variant).strip('\n').split('\t')[-1]

        if not sv_id in self.bnd_pool and sv_type != "SGL":
            self.bnd_pool[sv_id] = self.BP(variant.CHROM, variant.POS, variant.REF, 'unknown')
            sv_len = 0 if sv_type == 'TRA' else pos2 - variant.POS   # JZ: I am having issue getting SVLEN for INS type from GRIDSS.
            rec = self.Record(variant.CHROM, variant.POS, sv_id, variant.REF, f'<{sv_type}>', sv_type, qual, variant.FILTERS[0], chrom2, pos2, sv_len, variant.FORMAT, sample_col)
            self.write_record(rec)

class LUMPY(SV):
    
    def parse_variant(self, variant):
        sv_id = variant.ID.replace("_1","").replace("_2","")
        sv_type = variant.INFO.get("SVTYPE")
        sample_col = str(variant).rstrip('\n').split("\t")[-2]
        qual = '.' if variant.QUAL is None else variant.QUAL
        if sv_type == "BND":
            if not sv_id in self.bnd_pool:
                self.bnd_pool[sv_id] = self.BP(variant.CHROM, variant.POS, variant.REF, variant.INFO.get("STRANDS"))
            else:
                bp1 = self.bnd_pool[sv_id]
                sv_type = self.get_simple_svtype(variant, bp1)
                sv_len = 0 if sv_type == "TRA" else variant.POS - bp1.pos
                rec = self.Record(bp1.chrom, bp1.pos, sv_id, bp1.ref, f'<{sv_type}>', sv_type, qual, variant.FILTERS[0], variant.CHROM, variant.POS, sv_len, variant.FORMAT, sample_col)
                self.write_record(rec)
        else:
            sv_len = variant.INFO.get("SVLEN")
            pos2 = variant.INFO.get("END")
            rec = self.Record(variant.CHROM, variant.POS, sv_id, variant.REF, f'<{sv_type}>', sv_type, qual, variant.FILTERS[0], variant.CHROM, pos2, sv_len, variant.FORMAT, sample_col)
            self.write_record(rec)

# This only works for the published colo829 SV golden standard vcf file. TODO: should remove this and put it into a new file.
class COLO829(SV):

    def parse_variant(self, variant):
        sv_id = re.match(r"(.+)_[1,2]$", variant.ID).group(1)
        sv_type = "TRA" if variant.INFO['SVTYPE'] == "BND" else variant.INFO['SVTYPE']
        sv_len = 0 if variant.INFO.get('SVLEN') is None else int(variant.INFO['SVLEN'])
        qual = '.' if variant.QUAL is None else variant.QUAL
        sample_col = str(variant).rstrip('\n').split("\t")[-1]
        chrom = "chr" + variant.CHROM.replace("chr","")
        if variant.ID == "truthset_49_2": ## this is a bug in this vcf file
            pass
        elif sv_type == "INS":
            rec = self.Record(chrom,variant.POS, sv_id, variant.REF, f'<{sv_type}>', sv_type, qual, variant.FILTERS[0], chrom, variant.POS ,sv_len, variant.FORMAT, sample_col)
            self.write_record(rec)
        elif not sv_id in self.bnd_pool:
            self.bnd_pool[sv_id] = self.BP(chrom, variant.POS, variant.REF, None)
        else:
            bp1 = self.bnd_pool[sv_id]
            rec = self.Record(bp1.chrom, bp1.pos, sv_id, bp1.ref, f'<{sv_type}>', sv_type, qual, variant.FILTERS[0], chrom, variant.POS, sv_len, variant.FORMAT, sample_col)
            self.write_record(rec)