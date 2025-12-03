import pandas as pd
from pysam import AlignmentFile
import numpy as np
import sys


def main():
    base_dic = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    chrom = sys.argv[1]
    pos = int(sys.argv[2])
    ref = sys.argv[3]
    alt = sys.argv[4]
    bam_file = sys.argv[5]
    bam = AlignmentFile(bam_file, 'rb')
    coverage = bam.count_coverage(chrom, pos - 1, pos)
    reads_a = coverage[0][0]
    reads_c = coverage[1][0]
    reads_g = coverage[2][0]
    reads_t = coverage[3][0]
    ratio = coverage[base_dic[alt]][0] / (reads_a + reads_c + reads_g + reads_t)
    print(f'{chrom}-{pos}-{ref}-{alt}支持reads数目: A:{reads_a}, C:{reads_c}, G:{reads_g}, T:{reads_t}')
    print(f'{chrom}-{pos}-{ref}-{alt}突变ratio: {ratio}')
    bam.close()


if __name__ == '__main__':
    main()

