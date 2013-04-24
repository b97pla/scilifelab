""" Compare SNP calls from vcf file with matching samples from a ped formatted
file.

'Statistical test should include allele frquency in normal population'
"""
import argparse

from itertools import izip


def compare_snp_calls(vcf_call, gt_call):
    if vcf_call == 'NA' or not set(['ACGT ']).issuperset(vcf_call):
        return 'NA'

    # TODO: Implement


def ped_line_reader(line):
    s = []
    for c in line:
        if c == '\t':
            yield ''.join(s)
            s = []
            continue

        s.append(c)

    yield ''.join(s)


def read_ped(ped_file, sample):
    genotypes = {}
    with open(ped_file) as fh:
        markers = fh.next().split('\t')[5:]
        for line in fh:
            pline = ped_line_reader(line)
            pline.next()
            sample_id = pline.next()
            if sample_id != sample:
                continue

            for i in range(4):
                pline.next()

            for marker, call in izip(markers, pline):
                genotypes[marker] = call


def read_snp_map(bed_file):
    snp_map = {}
    with open(bed_file) as fh:
        for line in fh:
            chrom, pos, _, sample_id, _ = line.split('\t')
            snp_map[sample_id] = '{0} {1}'.format(chrom, pos)

    return snp_map


def read_varscan_vcf(vcf_file, min_depth):
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

        result = compare_snp_calls(vcf[marker_pos], genotypes[marker])

        if result == 'NA':
            uncallable += 1

        else:
            if result < 0.1:
                identical += 1

            totaldiff += result
            is_callable += 1

    different = is_callable - identical
    print('Sample {0}: SAME: {1} DIFFERENT: {1} UNCALLABLE: {2}'
          '\n'.format(args.sample, identical, different, uncallable))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('bed_file')
    parser.add_argument('vcf_file')
    parser.add_argument('ped_file')

    args = parser.parse_args()

    main(args)
