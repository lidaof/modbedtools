import sys
import os
import gzip
import pysam

# SAM/BAM base modification codes -> human readable abbreviation.
# Single letter codes and their ChEBI numeric equivalents both map here, per the
# SAMtags specification (https://samtools.github.io/hts-specs/SAMtags.pdf).
# Used to name the output file for each (base, modification) found in the input,
# so e.g. 5mC and 5hmC (both C mods) go to separate files instead of being mixed.
MOD_CODE_NAMES = {
    # cytosine
    'm': '5mC', 'h': '5hmC', 'f': '5fC', 'c': '5caC', 'C': 'modC',
    27551: '5mC', 76792: '5hmC', 76794: '5fC', 76793: '5caC',
    # thymine / uracil
    'g': '5hmU', 'e': '5fU', 'b': '5caU', 'T': 'modT', 'U': 'modU',
    16964: '5hmU', 80961: '5fU', 17477: '5caU',
    # adenine
    'a': '6mA', 'A': 'modA', 28871: '6mA',
    # guanine
    'o': '8oxoG', 'G': 'modG', 44605: '8oxoG',
    # any base
    'n': 'Xao', 'N': 'modN', 18107: 'Xao',
}


def mod_label(base, mod_code):
    '''
    human readable label for a (base, modification code) pair, used in output file names.
    falls back to base+code (e.g. 'C21839') for modifications not in the spec table.
    '''
    if mod_code in MOD_CODE_NAMES:
        return MOD_CODE_NAMES[mod_code]
    return f'{base}{mod_code}'


def process_read(read, cutoff, cpg):
    '''
    convert bam file with Ml/Mm tags to bed file with methylation information in format: chr, start, end, name, score, strand, methylated position array, unmethylated position array.
    One line per read.
    if cpg mode is on: for pacbio bam, it assumes C in both strands of an CpG has same methylation level, both C will show at bp level vis
    Returns a dict keyed by (base, modification code), e.g. ('C', 'm') for 5mC,
    ('C', 'h') for 5hmC, ('A', 'a') for 6mA, so different modifications on the same
    base are kept separate instead of being mixed together.
    '''
    if not (read.is_supplementary or read.is_secondary or read.is_unmapped):
        chrom = read.reference_name
        start = read.reference_start
        end = read.reference_end
        name = read.query_name
        align = read.get_aligned_pairs(matches_only=True)
        alignd = {}
        for x in align:
            # aligned dict, key: query pos in read, value: ref pos
            alignd[x[0]] = x[1]
        # modbase_key = ('C', 1, 'm') if read.is_reverse else ('C', 0, 'm')
        # if base == 'A':
        #     modbase_key = ('A', 1, 'a') if read.is_reverse else ('A', 0, 'a')
        # if modbase_key not in read.modified_bases:
        #     return []
        modified = read.modified_bases
        # modified_bases is None when the MM/ML tags fail to parse (malformed tags)
        # and an empty dict when the read simply carries no modifications; skip both.
        if not modified:
            return {}
        modbase_keys = list(modified.keys())
        mod_dict = {} # key: (base, modification code); base can be A or C, T is combined into A. value: the output list
        '''
        modified_bases keys are (canonical_base, strand, modification_code) tuples,
        the modification_code (e.g. 'm'=5mC, 'h'=5hmC, 'a'=6mA) can be a single
        letter or a ChEBI integer:
        >>> a[0].modified_bases.keys()
        dict_keys([('A', 1, 'a'), ('C', 1, 'm'), ('T', 0, 'a')])
        >>> a[2].modified_bases.keys()
        dict_keys([('T', 0, 'a'), ('C', 1, 'm'), ('A', 1, 'a')])
        '''
        for modbase_key in modbase_keys:
            real_base = modbase_key[0]
            mod_code = modbase_key[2]
            if real_base == 'T':
                base = 'A'
            else:
                base = real_base
            modbase_list = modified[modbase_key]
            modbase_methy_string = '.'
            modbase_unmet_string = '.'
            modbase_methy_list = []
            modbase_unmet_list = []
            strand = '-' if read.is_reverse else '+'
            for j in modbase_list:
                if j[0] in alignd:
                    if j[1]/255. >= cutoff:  # methylated base
                        if read.is_reverse:
                            if cpg:
                                modbase_methy_list.append(
                                    str(alignd[j[0]] - start - 1))
                            if real_base == 'T': # use the A on the other strand, add - sign
                                modbase_methy_list.append(str(-(start - alignd[j[0]])))
                            else:
                                modbase_methy_list.append(str(start - alignd[j[0]]))
                        else:
                            if real_base == 'T':
                                modbase_methy_list.append(str(-(alignd[j[0]] - start)))
                            else:
                                modbase_methy_list.append(str(alignd[j[0]] - start))
                            if cpg:
                                modbase_methy_list.append(
                                    str(-(alignd[j[0]] - start+1)))
                    else:
                        if read.is_reverse:
                            if cpg:
                                modbase_unmet_list.append(
                                    str(alignd[j[0]] - start - 1))
                            if real_base == 'T':
                                modbase_unmet_list.append(str(-(start - alignd[j[0]])))
                            else:
                                modbase_unmet_list.append(str(start - alignd[j[0]]))
                        else:
                            if real_base == 'T':
                                modbase_unmet_list.append(str(-(alignd[j[0]] - start)))
                            else:
                                modbase_unmet_list.append(str(alignd[j[0]] - start))
                            if cpg:
                                modbase_unmet_list.append(
                                    str(-(alignd[j[0]] - start+1)))
            if len(modbase_methy_list):
                modbase_methy_string = ','.join(modbase_methy_list)
            if len(modbase_unmet_list):
                modbase_unmet_string = ','.join(modbase_unmet_list)
            mod_key = (base, mod_code)
            if mod_key in mod_dict:
                # combine info from different modified base keys of the same
                # (base, modification), e.g. 6mA reported on both A and T strands
                existing = mod_dict[mod_key]
                # existing: [chrom, str(start), str(end), name, '0', strand, modbase_methy_string,
                #             modbase_unmet_string]
                if existing[6] == '.':
                    existing[6] = modbase_methy_string
                elif modbase_methy_string != '.':
                    existing[6] = existing[6] + ',' + modbase_methy_string
                if existing[7] == '.':
                    existing[7] = modbase_unmet_string
                elif modbase_unmet_string != '.':
                    existing[7] = existing[7] + ',' + modbase_unmet_string
                mod_dict[mod_key] = existing
            else:
                mod_dict[mod_key] = [chrom, str(start), str(end), name, '0', strand, modbase_methy_string,
                    modbase_unmet_string]
        return mod_dict
    else:
        return {}


def bam2mod(bamfile, outfile, cutoff=0.5, cpg=False, reference=None):
    # remove 'rb' mode for auto-detect
    # For CRAM files, reference_filename is required
    print(f'[info] reading file {bamfile}', file=sys.stderr)
    if reference:
        assert os.path.exists(reference), f'Reference file {reference} does not exist'
        bam = pysam.AlignmentFile(bamfile, reference_filename=reference, check_sq=False)
    else:
        bam = pysam.AlignmentFile(bamfile, check_sq=False)
    hasIndex = False
    # Check for both BAM (.bai) and CRAM (.crai) index files
    if os.path.exists(bamfile+'.bai') or os.path.exists(bamfile+'.crai'):
        hasIndex = True
        # num_reads = bam.count()  # this needs index
        # print(f'[info] total reads: {num_reads}', file=sys.stderr)
    cpgtag = '.cpg' if cpg else ''
    fhs = {} # file handles keyed by output file name, opened on the fly per (base, modification)
    # this makes bam index optional
    for read in bam.fetch(until_eof=(not hasIndex)):
        # process_read returns a dict keyed by (base, modification code); the base
        # and modification are discovered dynamically from the reads' MM/ML tags, so
        # any base/modification present (5mC, 5hmC, 6mA, ...) gets its own output file.
        items = process_read(read, cutoff, cpg)
        for (base, mod_code), line_items in items.items():
            # need either a modified or an unmodified base with an aligned position;
            # skip empty entries so we don't create spurious 0-byte output files
            if not line_items or (line_items[6] == '.' and line_items[7] == '.'):
                continue
            outf_base = f'{outfile}{cpgtag}.{mod_label(base, mod_code)}.modbed'
            if outf_base not in fhs:
                print(f'[info] writing file {outf_base}', file=sys.stderr)
                fhs[outf_base] = open(outf_base, 'w')
            out = fhs[outf_base]
            line = '\t'.join(line_items)
            out.write(line+'\n')

    for fh in fhs.values():
        fh.close()
        


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
    sbase = base.lower()
    rc_base = rct[base]
    rc_sbase = rct[base].lower()
    fa = read_fa(fasta_file)
    outf = '{}.modbed'.format(output)
    print(f'[info] writing file {outf}', file=sys.stderr)
    with xopen(bedfile) as fin, open(outf, 'w') as out:
        for line in fin:
            t = line.strip().split('\t')
            # the fiber-seq data should not count first and last one as told
            bs = t[-1].split(',')[1:-1]
            if not len(bs):
                continue
            chrom = t[0]
            start = int(t[1])
            end = int(t[2])
            # bs = t[-1].split(',')
            bs2 = [int(x) for x in bs]
            # contains start position in genome for each methylated base
            bs3 = [start+x for x in bs2]
            s = fa[chrom]  # the sequence
            c1 = []
            c2 = []
            c1str = '.'
            c2str = '.'
            for x in range(start, end+1):
                if x >= len(s):  # Check if x exceeds the chromosome length, from @mitomac @github
                    continue  # Skip if it does
                if x in bs3:
                    # a methylated base
                    if s[x] == base or s[x] == sbase:
                        # + strand
                        c1.append(str(x-start))
                    elif s[x] == rc_base or s[x] == rc_sbase:
                        # - strand
                        c1.append(str(start-x))
                else:
                    # an unmethylated base but need check seq base
                    if s[x] == base or s[x] == sbase:
                        # + strand
                        c2.append(str(x-start))
                    elif s[x] == rc_base or s[x] == rc_sbase:
                        # - strand
                        c2.append(str(start-x))
            if len(c1):
                c1str = ','.join(c1)
            if len(c2):
                c2str = ','.join(c2)
            if c1str != '.' or c2str != '.':
                out.write(
                    f"{chrom}\t{start}\t{end}\t{t[3]}\t{t[4]}\t+\t{c1str}\t{c2str}\n")


def main():
    pass


if __name__ == "__main__":
    main()
