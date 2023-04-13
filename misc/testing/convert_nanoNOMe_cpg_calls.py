#!/usr/bin/env python3

'''
This script convert CpG calls file from nanoNOMe output to modbed format
nanoNOMe results can be found at https://zenodo.org/record/3969567#.ZCztE-zMI52
'''

import sys


def write_record(fh, c, prid, rs, re, m, u):
    rm = [str(x-rs) for x in m]
    ru = [str(x-rs) for x in u]
    fh.write(f'{c}\t{rs}\t{re}\t{prid}\t{",".join(rm)}\t{",".join(ru)}\n')


m = []  # meth positions
u = []  # unmeth positions
rs = 0  # read start, assume smallest
re = 0  # read end, assume biggest
prid = ''  # read id
c = ''
# with open('../test/200731_GM12878_nanoNOMe_CTCF_motifs_10kb_singleread_cpg_calls_first_100k_lines.txt') as fin, open('../test/nanoNOMe_test.modbed', 'w') as fout:
with open('200731_GM12878_nanoNOMe_CTCF_motifs_10kb_singleread_cpg_calls.txt') as fin, open('200731_GM12878_nanoNOMe_CTCF_motifs_10kb_singleread_cpg_calls.modbed', 'w') as fout:
    next(fin)
    for line in fin:
        t = line.strip().split()
        rid = t[3]
        s = int(t[1])
        if rid != prid:
            # a new read record
            if prid != '':
                write_record(fout, c, prid, rs, re, m, u)
                m = []
                u = []
            rs = s - 2
            if t[4] == '1':
                m.append(s)
            elif t[4] == '0':
                u.append(s)
        else:
            c = t[0]
            re = s + 2
            if t[4] == '1':
                m.append(s)
            elif t[4] == '0':
                u.append(s)
        prid = rid
    write_record(fout, c, prid, rs, re, m, u)  # last record
