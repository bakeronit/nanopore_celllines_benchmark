import sys, argparse

import sv

parser = argparse.ArgumentParser(
    prog = "simple_sv_type",
    description = "Formating the SV events called by different tools to a simple four-type SV VCF format that SV merge tool jasmine can use"
)
parser.add_argument("filename")
parser.add_argument("-t", "--type", help="The type of the filename provided, can be delly, savana, nanomonsv, severus and qsv")
parser.add_argument("-o","--output",nargs='?', type=argparse.FileType('w'),default=sys.stdout, help="The output filename")
parser.add_argument("-r","--reference", nargs="?", help="The path to reference fasta, required for convert nanomonsv, and qsv output")
args = parser.parse_args()



type = args.type.lower()

if type in ['nanomonsv','qsv'] and args.reference == "" :
    raise "ERROR: reference fasta file is need for converting nanomonsv and qsv output for header information"
    sys.exit()

match args.type:
    case "delly":
        sv_parser = sv.DELLY()
    case "savana":
        sv_parser = sv.SAVANA()
    case "severus":
        sv_parser = sv.SEVERUS()
    case "nanomonsv":
        sv_parser = sv.NANOMONSV(reference=args.reference)
    case "qsv":
        sv_parser = sv.QSV(reference=args.reference)
    case "gridss":
        sv_parser = sv.GRIDSS()
    case "lumpy":
        sv_parser = sv.LUMPY()
    case "colo829":
        sv_parser = sv.COLO829()
    case _:
        raise "ERROR: currently only support delly, savana, severus (v1^), nanomonsv (txt output), and qsv (dcc) format"

file_in = sv_parser.readSv(args.filename)
sv_parser.output = args.output

if not args.type in ["nanomonsv","qsv"]:
    sv_parser.write_header(file_in.raw_header)
else:
    sv_parser.write_header()

for variant in file_in:
    sv_parser.parse_variant(variant)
