"""Convert a MiSeq samplesheet to a valid HiSeq samplesheet
"""
import argparse

from scilifelab.illumina.miseq import MiSeqSampleSheet
from scilifelab.illumina.hiseq import HiSeqSampleSheet

import pprint
pp = pprint.PrettyPrinter(indent=2)

def main(miseq_samplesheet):
    m_samplesheet = MiSeqSampleSheet(miseq_samplesheet)

    pp.pprint(m_samplesheet.samples)
    print
    pp.pprint(m_samplesheet.to_hiseq())

    h_samplesheet = HiSeqSampleSheet("scripts/C1DLBACXX.csv")
    h_style_ss = HiSeqSampleSheet(m_samplesheet.to_hiseq())

    print
    pp.pprint(h_samplesheet[:2])
    print
    pp.pprint(h_style_ss)

    h_style_ss.write("hiseqstyle.csv")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("miseq_samplesheet", \
        help="Samplesheet in MiSeq Samplesheet format")

    args = parser.parse_args()

    main(args.miseq_samplesheet)
