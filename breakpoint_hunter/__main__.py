"""精确断点deletion hunter
"""
import os
from multiprocessing import Pool, cpu_count
from collections import namedtuple
from argparse import ArgumentParser
import pandas as pd
import yaml
from .breakpoint_hunter import BreakpointHunter


ROOT = os.path.dirname(os.path.abspath(__file__))
CONF_DIR = os.path.join(ROOT, 'config')
Param = namedtuple('Param', ['work_dir', 'sample', 'bam', 'chrom', 'bp_left', 'bp_right', 'cnv_name', 'thres_normal',
                             'thres_hom_del', 'thres_support_reads', 'thres_mean_reads', 'extend', 'min_len'])


def read_conf(conf_file) -> dict:
    """读取配置文件

    :param conf_file: 配置文件路径
    :type conf_file: str
    :return: 配置文件内容字典
    :rtype: dict
    """
    with open(conf_file, 'r', encoding='utf-8') as y:
        conf_dic = yaml.load(y, Loader=yaml.FullLoader)
    return conf_dic


def main():
    """参数解析
    """
    parser = ArgumentParser()
    parser.add_argument('-extend', help='breakpoint extend region', type=int, default=350)
    parser.add_argument('-min_len', help='del cnv min length', type=int, default=300)
    parser.add_argument('-out', help='output result file', required=True)
    parser.add_argument('-work_dir', help='work dir', required=True)
    parser.add_argument('-cpu', help='cpu number', type=int, default=1)
    parser.add_argument('-info', help='sample info')
    parser.add_argument('-genome', help='genome build version', type=str, default='GRCh37')
    parser.add_argument('-conf', help='config file')
    parsed_args = parser.parse_args()
    genome = 'GRCh38' if '38' in parsed_args.genome else 'GRCh37'
    cpu = min(cpu_count(), parsed_args.cpu)
    work_dir = parsed_args.work_dir
    if not os.path.exists(work_dir):
        os.makedirs(work_dir)
    conf_dic = read_conf(parsed_args.conf) if parsed_args.conf else read_conf(f'{CONF_DIR}/conf.yaml')
    cnv_info = pd.read_csv(os.path.join(ROOT, conf_dic[f'info_{genome}']), sep='\t')
    info = pd.read_csv(parsed_args.info, sep='\t', header=None, names=['sample', 'bam'])
    jobs = []
    for i in info.itertuples(False):
        for j in cnv_info.itertuples(False):
            jobs.append(Param(work_dir, i.sample, i.bam, j.chrom, j.bp_left, j.bp_right, j.cnv_name, j.thres_normal, j.thres_hom_del, j.thres_support_reads, j.thres_mean_reads, parsed_args.extend, parsed_args.min_len))
    with Pool(cpu) as pool:
        result_list = pool.map(BreakpointHunter.run, jobs)
    df = pd.concat(result_list, axis=0)
    df.sort_values(by=['sample', 'cnv_name'], inplace=True)
    df.to_csv(parsed_args.out, sep='\t', index=False)


if __name__ == '__main__':
    main()
