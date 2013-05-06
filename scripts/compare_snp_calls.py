""" Compare SNP calls from vcf file with matching samples from a ped formatted
file.

'Statistical test should include allele frquency in normal population'
"""
import argparse

from itertools import izip


def compare_snp_calls(gt_call, vcf_call):
    """ Compare the calls from both data sources.

    'replace with proper statistical test when we have real data'
    """
    if vcf_call == 'NA':
        print(vcf_call)
        return 'NA'

    if not set('ACGT ').issuperset(gt_call):
        return 'NA'

    allele_freq, ref_allele, alt_allele = vcf_call.split(':')

    if gt_call == '{0} {1}'.format(ref_allele, ref_allele):
        gt_call = 0

    elif gt_call == '{0} {1}'.format(alt_allele, alt_allele):
        gt_call = 1

    elif ref_allele not in gt_call and alt_allele not in gt_call:
        gt_call = 0.5

    # 'If alleles do not match, e.g. genotyping gives C/T and SNP calling A/G!'
    else:
        print('')
        print(gt_call)
        print('{0} {1}'.format(ref_allele, ref_allele))
        print('{0} {1}'.format(alt_allele, alt_allele))
        print('')
        return 'NA'

    diff = gt_call - float(allele_freq)
    return diff


def ped_line_reader(line):
    """ Since lines in ped files can be very long, read one tab seperated
    field at a time.
    """
    s = []
    for c in line:
        if c in '\t\n':
            yield ''.join(s)
            s = []
            continue

        s.append(c)


def read_ped(ped_file, sample):
    """ Parse ped file, extract genotypes for given sample.
    """
    genotypes = {}
    with open(ped_file) as fh:
        markers = fh.next().split('\t')[5:]
        for line in fh:
            pline = ped_line_reader(line)

            # Skip 'FAMILY' columns
            pline.next()

            sample_id = pline.next()
            if sample_id != sample:
                continue

            # Skip 'FATHER', 'MOTHER', 'SEX', 'AFFECTION_STATUS' columns
            for i in range(4):
                pline.next()

            for marker, call in izip(markers, pline):
                genotypes[marker] = call

    return genotypes


def read_snp_map(bed_file):
    """ Parse out relevant information from .bed file
    """
    snp_map = {}
    with open(bed_file) as fh:
        for line in fh:
            chrom, pos, _, sample_id, _ = line.split('\t')
            snp_map[sample_id] = '{0} {1}'.format(chrom, pos)

    return snp_map


def read_varscan_vcf(vcf_file, min_depth):
    """ Parse out relevant information from .vcf file
    """
    vcf = {}
    with open(vcf_file) as fh:
        for line in fh:
            row = line.split('\t')
            chrom = row[0]
            bp = row[1]
            ref = row[3]
            alt = row[4]
            info = row[7]

            read_depth = int(info.partition('DP=')[-1].partition(';')[0])
            allele_freq = float(info.partition('AF1=')[-1].partition(';')[0])

            if read_depth >= min_depth:
                vcf['{0} {1}'.format(chrom, bp)] = '{0}:{1}:{2}'.format(allele_freq, ref, alt)

            else:
                vcf['{0} {1}'.format(chrom, bp)] = 'NA'

    return vcf


def main(args):
    snp_map = read_snp_map(args.bed_file)
    vcf = read_varscan_vcf(args.vcf_file, args.min_depth)
    genotypes = read_ped(args.ped_file, args.sample)

    is_callable = 0
    uncallable = 0
    identical = 0
    totaldiff = 0.0

    for marker in snp_map.keys():
        marker_pos = snp_map[marker]

        if marker_pos not in vcf or marker not in genotypes:
            uncallable += 1
            continue

        result = compare_snp_calls(genotypes[marker], vcf[marker_pos])

        if result == 'NA':
            uncallable += 1

        else:
            if result < 0.1:
                identical += 1

            totaldiff += result
            is_callable += 1

    different = is_callable - identical
    print('Sample {0}: SAME: {1} DIFFERENT: {2} UNCALLABLE: {3}'
          '\n'.format(args.sample, identical, different, uncallable))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('bed_file')
    parser.add_argument('vcf_file')
    parser.add_argument('ped_file')
    parser.add_argument('sample')
    parser.add_argument('--min-depth', type=int, default=6)

    args = parser.parse_args()

    main(args)
