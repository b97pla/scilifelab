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

    print
    pp.pprint(h_samplesheet[:2])



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("miseq_samplesheet", \
        help="Samplesheet in MiSeq Samplesheet format")

    args = parser.parse_args()

    main(args.miseq_samplesheet)
