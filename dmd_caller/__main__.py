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
    parser.add_argument('-conf', help='config file')
    parser.add_argument('-cpu', help='cpu count', type=int, default=1)
    parser.add_argument('-bed', help='config bed')
    parser.add_argument('-gender', help='gender info')
    parser.add_argument('-gene', help='gene list')
    parsed_args = parser.parse_args()
    run(parsed_args)


if __name__ == '__main__':
    main()
