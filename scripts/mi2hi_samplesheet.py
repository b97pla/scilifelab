"""Convert a MiSeq samplesheet to a valid HiSeq samplesheet
"""
import argparse

from scilifelab.illumina.miseq import MiSeqSampleSheet


def main(miseq_samplesheet):
	samplesheet = MiSeqSampleSheet(miseq_samplesheet)
	print(samplesheet.samples)
	print(samplesheet.to_hiseq())


if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	parser.add_argument("miseq_samplesheet", \
		help="Samplesheet in MiSeq Samplesheet format")

	args = parser.parse_args()

	main(args.miseq_samplesheet)
