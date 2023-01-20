import sys
import os
import gzip
import pysam


def process_read(read, cutoff):
    '''
    convert bam file with Ml/Mm tags to bed file with methylation information in format: chr, start, end, methylated position array, unmethylated position array.
    One line per read.
    '''
    if not (read.is_supplementary or read.is_secondary or read.is_unmapped):
        line = []
        chrom = read.reference_name
        start = read.reference_start
        end = read.reference_end
        align = read.get_aligned_pairs(matches_only=True)
        alignd = {}
        for x in align:
            # aligned dict, key: query pos in read, value: ref pos
            alignd[x[0]] = x[1]
        modbase_key = ('C', 1, 'm') if read.is_reverse else ('C', 0, 'm')
        if modbase_key not in read.modified_bases:
            return []
        modbase_list = read.modified_bases[modbase_key]
        modbase_methy_string = ''
        modbase_unmet_string = ''
        modbase_methy_list = []
        modbase_unmet_list = []
        for j in modbase_list:
            if j[0] in alignd:
                if j[1]/255. >= cutoff:  # methylated base
                    if read.is_reverse:
                        modbase_methy_list.append(str(start - alignd[j[0]]))
                    else:
                        modbase_methy_list.append(str(alignd[j[0]] - start))
                else:
                    if read.is_reverse:
                        modbase_unmet_list.append(str(start - alignd[j[0]]))
                    else:
                        modbase_unmet_list.append(str(alignd[j[0]] - start))
        if len(modbase_methy_list):
            modbase_methy_string = ','.join(modbase_methy_list)
        if len(modbase_unmet_list):
            modbase_unmet_string = ','.join(modbase_unmet_list)
        line = [chrom, str(start), str(end), modbase_methy_string,
                modbase_unmet_string]
        return line
    else:
        return []


def bam2mod(bamfile, outfile, cutoff=0.5):
    bam = pysam.AlignmentFile(bamfile, 'rb', check_sq=False)
    if os.path.exists(bamfile+'.bai'):
        num_reads = bam.count()  # this needs index
        print(f'[info] total reads: {num_reads}', file=sys.stderr)
    outf = '{}.modbed'.format(outfile)
    print(f'[info] writing file {outf}', file=sys.stderr)
    with open(outf, "w") as out:
        for read in bam.fetch(until_eof=True):  # this makes bam index optional
            items = process_read(read, cutoff)
            if len(items):
                if items[3] or items[4]:  # for output, need either has modified base or unmodified base
                    line = '\t'.join(items)
                    out.write(line+'\n')


rct = {
    'A': 'T',
    'C': 'G',
    'T': 'A',
    'G': 'C',
}


def xopen(fn):
    if not os.path.exists(fn):
        print(f'[error]: {fn} not exist, please check', file=sys.stderr)
        sys.exit(2)
    if fn.endswith('.gz'):
        return gzip.open(fn, 'rt')
    else:
        return open(fn, 'r')


def read_fa(fasta_file):
    print(f'[info] reading file {fasta_file}...', file=sys.stderr)
    d = {}
    k = ''
    s = ''
    with xopen(fasta_file) as fin:
        for line in fin:
            if line.startswith('>'):
                if s:
                    d[k] = s
                k = line.strip().lstrip('>')
                s = ''
            else:
                s += line.strip()
        d[k] = s
    return d


def addbg(bedfile, fasta_file, output, base):
    fa = read_fa(fasta_file)
    outf = '{}.modbed'.format(output)
    print(f'[info] writing file {outf}', file=sys.stderr)
    with xopen(bedfile) as fin, open(outf, 'w') as out:
        for line in fin:
            t = line.strip().split('\t')
            chrom = t[0]
            start = int(t[1])
            end = int(t[2])
            # the fiber-seq data should not count first and last one as told
            bs = t[-1].split(',')[1:-1]
            # bs = t[-1].split(',')
            bs2 = [int(x) for x in bs]
            # contains start position in genome for each methylated base
            bs3 = [start+x for x in bs2]
            s = fa[chrom]  # the sequence
            c1 = []
            c2 = []
            for x in range(start, end+1):
                if x in bs3:
                    # a methylated base
                    if s[x].upper() == base:
                        # + strand
                        c1.append(str(x-start))
                    elif s[x].upper() == rct[base]:
                        # - strand
                        c1.append(str(start-x))
                else:
                    # an unmethylated base but need check seq base
                    if s[x].upper() == base:
                        # + strand
                        c2.append(str(x-start))
                    elif s[x].upper() == rct[base]:
                        # - strand
                        c2.append(str(start-x))
            out.write(
                f"{chrom}\t{start}\t{end}\t{','.join(c1)}\t{','.join(c2)}\n")


def main():
    pass


if __name__ == "__main__":
    main()
