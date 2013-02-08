"""Convert a MiSeq samplesheet to a valid HiSeq samplesheet
"""
import argparse

from scilifelab.illumina.miseq import MiSeqSampleSheet


def main():
	pass


if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	parser.add_argument("miseq_samplesheet", \
		help="Samplesheet in MiSeq Samplesheet format")

	args = parser.parse_args()

	main(args.miseq_samplesheet)
