import argparse
import re

parser = argparse.ArgumentParser(description='Merge CCS reports.')
parser.add_argument('ccs_report', metavar='R', type=str, nargs='+', help='CCS report(s)')
args = parser.parse_args()

d = {}

for ccs_report in args.ccs_report:
    file = open(ccs_report, "r")
    for line in file:
        if len(line) > 1 and 'Exclusive ZMW counts for' not in line:
            a = line.rstrip().split(":")

            k = a[0].rstrip()
            v = float(re.sub(" ", "", re.sub(" \(.*$", "", a[1])))

            if k not in d:
                d[k] = 0.0;

            d[k] = d[k] + v

print(f'ZMWs input          (A)  : {d["ZMWs input          (A)"]}')
print(f'ZMWs generating CCS (B)  : {d["ZMWs generating CCS (B)"]} ({(100.0*d["ZMWs generating CCS (B)"]/d["ZMWs input          (A)"]):.2f}%)')
print(f'ZMWs filtered       (C)  : {d["ZMWs filtered       (C)"]} ({(100.0*d["ZMWs filtered       (C)"]/d["ZMWs input          (A)"]):.2f}%)')
print(f'')
print(f'Exclusive ZMW counts for (C):')
print(f'Median length filter     : {d["Median length filter"]} ({(100.0*d["Median length filter"]/d["ZMWs filtered       (C)"]):.2f}%)')
print(f'Below SNR threshold      : {d["Below SNR threshold"]} ({(100.0*d["Below SNR threshold"]/d["ZMWs filtered       (C)"]):.2f}%)')
print(f'Lacking full passes      : {d["Lacking full passes"]} ({(100.0*d["Lacking full passes"]/d["ZMWs filtered       (C)"]):.2f}%)')
print(f'Heteroduplex insertions  : {d["Heteroduplex insertions"]} ({(100.0*d["Heteroduplex insertions"]/d["ZMWs filtered       (C)"]):.2f}%)')
print(f'Coverage drops           : {d["Coverage drops"]} ({(100.0*d["Coverage drops"]/d["ZMWs filtered       (C)"]):.2f}%)')
print(f'Insufficient draft cov   : {d["Insufficient draft cov"]} ({(100.0*d["Insufficient draft cov"]/d["ZMWs filtered       (C)"]):.2f}%)')
print(f'Draft too different      : {d["Draft too different"]} ({(100.0*d["Draft too different"]/d["ZMWs filtered       (C)"]):.2f}%)')
print(f'Draft generation error   : {d["Draft generation error"]} ({(100.0*d["Draft generation error"]/d["ZMWs filtered       (C)"]):.2f}%)')
print(f'Draft above --max-length : {d["Draft above --max-length"]} ({(100.0*d["Draft above --max-length"]/d["ZMWs filtered       (C)"]):.2f}%)')
print(f'Draft below --min-length : {d["Draft below --min-length"]} ({(100.0*d["Draft below --min-length"]/d["ZMWs filtered       (C)"]):.2f}%)')
print(f'Reads failed polishing   : {d["Reads failed polishing"]} ({(100.0*d["Reads failed polishing"]/d["ZMWs filtered       (C)"]):.2f}%)')
print(f'Empty coverage windows   : {d["Empty coverage windows"]} ({(100.0*d["Empty coverage windows"]/d["ZMWs filtered       (C)"]):.2f}%)')
print(f'CCS did not converge     : {d["CCS did not converge"]} ({(100.0*d["CCS did not converge"]/d["ZMWs filtered       (C)"]):.2f}%)')
print(f'CCS below minimum RQ     : {d["CCS below minimum RQ"]} ({(100.0*d["CCS below minimum RQ"]/d["ZMWs filtered       (C)"]):.2f}%)')
print(f'Unknown error            : {d["Unknown error"]} ({(100.0*d["Unknown error"]/d["ZMWs filtered       (C)"]):.2f}%)')