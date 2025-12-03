from collections import defaultdict
from argparse import ArgumentParser
from multiprocessing.dummy import Process
from multiprocessing import Queue
from itertools import product
from functools import partial
import portion
import pysam
import glob
import re
import os


import logging
logger = logging.getLogger(__name__)
formatter = logging.Formatter(
    '[%(asctime)s-%(filename)s:%(lineno)s-%(levelname)s] %(message)s', '%Y-%m-%d %H:%M:%S'
)
handler = logging.StreamHandler()
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.DEBUG)


class GenomicRegion:
    """
    基因组区域
    """
    def __init__(self):
        self._chromosomes = defaultdict(portion.empty)

    def union(self, contig, start, end, extend=0):
        """
        合入一个基因组区域
        :param contig: contig, 染色体编号
        :param start: 起始位置
        :param end: 终止位置
        :param extend: 尝试进行左右进行延伸，延伸后若与其他区域有重叠，则执行延伸，否则不进行延伸
        """
        target_region = portion.closed(start, end)
        extend_region = target_region.replace(lower=lambda x: x - extend)
        if self._chromosomes[contig].overlaps(extend_region):
            target_region = extend_region
        extend_region = target_region.replace(upper=lambda x: x + extend)
        if self._chromosomes[contig].overlaps(extend_region):
            target_region = extend_region
        self._chromosomes[contig] |= target_region

    def __iter__(self):
        """
        迭代所有的基因组区域
        """
        for contig, regions in self._chromosomes.items():
            for region in regions:
                yield contig, region.lower, region.upper

    def __len__(self):
        """
        计算所有区域的个数
        """
        return sum(len(regions) for regions in self._chromosomes.values())


def get_mate_regions(job_queue, region_queue, reads_length):
    """
    搜索给定区域内所有mate落入的所有区域
    :param job_queue: 待查找的所有区域队列
    :param region_queue: 所有mate落入的区域队列
    :param reads_length: reads 长度，即左右延伸的长度
    """
    mate_regions = GenomicRegion()
    # 查找队列内所有待查找区域
    for filename, contig, start, end in iter(job_queue.get, 'END'):
        # reads所在区域
        mate_regions.union(contig, start, end)
        # 读取文件
        bam = pysam.AlignmentFile(filename, 'rb')
        try:
            # 查找所有的目标区域
            for record in bam.fetch(contig, start, end):
                mate_regions.union(
                    record.next_reference_name,
                    record.next_reference_start,
                    record.next_reference_start + 1,
                    reads_length
                )
        except ValueError:
            continue
    # 结果放入队列
    for contig, start, end in mate_regions:
        region_queue.put((contig, start, end))
    region_queue.put('END')


def get_target_regions(filenames, regions, processes, reads_length=0):
    """
    获得所有目标区域及其mate所在的所有基因组区域
    :param filenames: 所有的文件
    :param regions: 所有的区域
    :param processes: 线程数
    :param reads_length: reads长度
    :returns 不重叠的所有目标区域
    """
    # merge regions
    target_regions = GenomicRegion()
    for contig, start, end in regions:
        target_regions.union(contig, start, end, reads_length)

    # create queue
    jobs = Queue()
    for filename, (contig, start, end) in product(filenames, target_regions):
        jobs.put((filename, contig, start, end))

    # start workers
    mate_regions = Queue()
    worker = partial(get_mate_regions, jobs, mate_regions, reads_length)
    for _ in range(processes):
        Process(target=worker).start()
        # 放入标记技术区域搜索线程
        jobs.put('END')

    # 从结果队列中去重并进行合并
    for _ in range(processes):
        for region in iter(mate_regions.get, 'END'):
            target_regions.union(*region, reads_length)

    logger.debug('found %d regions ...', len(target_regions))

    for contig, start, end in target_regions:
        yield contig, start, end


def get_records(job_queue, record_queue):
    """
    从给定的文件和区域中读取reads
    :param job_queue: 文件和区域队列
    :param record_queue: 结果队列
    """
    for filename, contig, start, end in iter(job_queue.get, 'END'):
        bam = pysam.AlignmentFile(filename, 'rb')
        try:
            for record in bam.fetch(contig, start, end):
                record_queue.put(record.to_string())
        except ValueError:
            continue
    record_queue.put('END')


def run(filenames, regions, processes, reads_length=0):
    logger.debug('run in %s regions ...', len(regions))
    # 打开所有文件，确保文件都存在
    all_bam = [pysam.AlignmentFile(filename, 'rb') for filename in filenames]
    # 获取所有目标区及mate落入的区域
    target_regions = get_target_regions(filenames, regions, processes, reads_length)
    # create queue
    jobs = Queue()
    for filename, (contig, start, end) in product(filenames, target_regions):
        jobs.put((filename, contig, start, end))

    # 提取待提取区域的reads
    records = Queue()
    worker = partial(get_records, jobs, records)
    for _ in range(processes):
        Process(target=worker).start()
        # 放入结束标记，结束提取线程
        jobs.put('END')

    # 输出结果
    # header
    print(all_bam[0].header, end='')
    # records
    for _ in range(processes):
        for record in iter(records.get, 'END'):
            print(record)


def main():
    parser = ArgumentParser()
    parser.add_argument('bam', metavar='<in.bam>')
    parser.add_argument('region', metavar='<region>')
    parser.add_argument('-r', action='store_true', help='find bam in subdirectories recursively')
    parser.add_argument('-e', type=int, default=0, help='extend given region')
    parser.add_argument('-p', type=int, default=1, help='processes number')
    parser.add_argument('-l', type=int, default=0, help='reads length')
    parsed_args = parser.parse_args()
    filenames = glob.glob(parsed_args.bam)
    # 解析所有目标区
    regions = list()
    if parsed_args.region.endswith('.bed'):
        if os.path.isfile(parsed_args.region):
            with open(parsed_args.region) as f:
                for i, line in enumerate(f, start=1):
                    try:
                        contig, start, end = line.strip().split('\t')
                        start, end = sorted(map(int, (start, end)))
                        start, end = start - parsed_args.e, end + parsed_args.e
                        regions.append((contig, start, end))
                    except ValueError:
                        raise ValueError("Can't parse region from line {} from file \"{}\"".format(
                            i, parsed_args.region
                        ))
        else:
            raise ValueError("Can't read bed file: \"{}\"".format(parsed_args.region))
    else:
        try:
            contig, start, end = re.findall(r'(\w+):(\d+)-(\d+)', parsed_args.region).pop(0)
            start, end = sorted(map(int, (start, end)))
            start, end = start - parsed_args.e, end + parsed_args.e
            regions.append((contig, start, end))
        except IndexError:
            raise ValueError("Can't parse region \"{}\"".format(parsed_args.region))
    run(filenames, regions, parsed_args.p, parsed_args.l)
    logger.debug('finished')


if __name__ == '__main__':
    main()

