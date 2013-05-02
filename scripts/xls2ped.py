""" Convert Excel document (.xls) delivered from MAF to a ped.txt file.
"""
import argparse
import csv

import xlrd


def main(xls_file, out_file, sheet_name):
    with xlrd.open_workbook(xls_file) as workbook:
        worksheet = workbook.sheet_by_name(sheet_name)

        with open(out_file, 'w') as fh:
            c = csv.writer(fh, delimiter='\t')
            for row in range(worksheet.nrows):
                c.writerow(worksheet.row_values(row))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('xls_file')
    parser.add_argument('out_file')
    parser.add_argument('--sheet_name', default='HaploView_ped_0')

    args = parser.parse_args()

    main(args.xls_file, args.out_file, args.sheet_name)
