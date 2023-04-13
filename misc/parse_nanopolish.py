#!/usr/bin/env python3

'''
this script parses methylation calls from nanopolish to modbed format
it requires the bam file to get the mapping position of reads
'''

import sys
import pysam


def get_reads_mapping(bamfile):
    d = {}
    bam = pysam.AlignmentFile(bamfile, 'rb', check_sq=False)
    for read in bam.fetch(until_eof=True):
        if not (read.is_supplementary or read.is_secondary or read.is_unmapped):
            d[read.query_name] = [read.reference_name,
                                  read.reference_start, read.reference_end]
    return d


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print(f'''
usage: python3 {sys.argv[0]} bamfile nanopolish_methy_call_file
        ''')
        sys.exit(1)
    d = get_reads_mapping(sys.argv[1])
    outfile = f'{sys.argv[2].split(".")[0]}.modbed'
    m = {}  # key: read id, value: [ [meth positions], [unmeth positions]]
    with open(sys.argv[2]) as fin:
        for line in fin:
            t = line.strip().split('\t')
            lik = float(5)
            pos = int(t[2])
            rid = t[4]
            if rid not in m:
                m[rid] = [[], []]
            if lik >= 0:
                m[rid][0].append(pos)
            else:
                m[rid][1].append(pos)
    with open(outfile, 'w') as fout:
        for k in m:
            if not k in d:
                continue
            chrom, start, end = d[k]
            ms, um = m[k]
            mss = [str(x-start) for x in ms]
            ums = [str(x-start) for x in um]
            fout.write(
                f'{chrom}\t{start}\t{end}\t{k}\t{",".join(mss)}\t{",".join(ums)}\n')
