"""运行
"""
from argparse import ArgumentParser
from .utils import run


def main():
    """输入参数
    """
    parser = ArgumentParser()
    parser.add_argument('-mode', help='mode')
    parser.add_argument('-gc_fix', help='gc fix result file list')
    parser.add_argument('-control', help='control file info')
    parser.add_argument('-bam', help='bam info')
    parser.add_argument('-out_dir', help='out dir')
    parser.add_argument('-out', help='out file')
    parser.add_argument('-cpu', help='cpu count', type=int, default=1)
    parser.add_argument('-bed', help='target bed or gc bed')
    parser.add_argument('-gender', help='gender info')
    parser.add_argument('-win', help='win size', type=int, default=50)
    parser.add_argument('-slide', help='slide size', type=int, default=25)
    parser.add_argument('-win_min', help='win min size', type=int, default=40)
    parser.add_argument('-fasta', help='fasta file')
    parser.add_argument('-conf', help='conf file')
    parser.add_argument('-genome', help='genome version', default='GRCh37')
    parsed_args = parser.parse_args()
    run(parsed_args)


if __name__ == '__main__':
    main()
