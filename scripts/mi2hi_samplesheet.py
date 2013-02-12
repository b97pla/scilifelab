"""Convert a MiSeq samplesheet to a valid HiSeq samplesheet
"""
import argparse

from scilifelab.illumina.miseq import MiSeqSampleSheet
from scilifelab.illumina.hiseq import HiSeqSampleSheet


def main(miseq_samplesheet, hiseq_samplesheet):
    m_samplesheet = MiSeqSampleSheet(miseq_samplesheet)
    h_samplesheet = HiSeqSampleSheet(m_samplesheet.to_hiseq())

    if hiseq_samplesheet is None:
        hiseq_samplesheet = miseq_samplesheet + ".out"

    h_samplesheet.write(hiseq_samplesheet)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("miseq_samplesheet", \
        help="Samplesheet in MiSeq Samplesheet format")
    parser.add_argument("-o", "--out", \
        default=None, \
        help="Output Samplesheet in HiSeq format")

    args = parser.parse_args()

    main(args.miseq_samplesheet, args.out)
