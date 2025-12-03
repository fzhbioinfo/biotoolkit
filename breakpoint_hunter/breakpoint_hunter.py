"""
精确断点cnv(大片段del)检测
"""
from pysam import AlignmentFile
import numpy as np
import pandas as pd


class BreakpointHunter:
    """精确断点reads识别
    """
    def __init__(self, bam):
        self.bam = AlignmentFile(bam, 'rb')

    def close(self):
        """关闭bam
        """
        self.bam.close()

    def cal_mean_depth(self, chromosome: str, start: int, stop: int):
        """计算平均深度

        :param chromosome: 染色体
        :type chromosome: str
        :param start: 起始坐标
        :type start: int
        :param stop: 终止坐标
        :type stop: int
        :return: 平均深度
        :rtype: float
        """
        assert start < stop
        return np.sum(self.bam.count_coverage(chromosome, start, stop)) / (stop - start)

    @staticmethod
    def read_filter(read):
        """reads过滤

        :param read: read segment
        :type read: pysam.AlignedSegment
        :return: True|False
        :rtype: bool
        """
        return read.is_duplicate or read.is_qcfail or read.is_unmapped or read.is_secondary or read.is_supplementary or read.mate_is_unmapped

    def find_del_support_reads(self, chromosome: str, break_point_1: int, break_point_2: int, extend: int=350, min_length: int=300):
        """搜索支持reads

        :param chromosome: 染色体
        :type chromosome: str
        :param break_point_1: 左侧断点坐标
        :type break_point_1: int
        :param break_point_2: 右侧断点坐标
        :type break_point_2: int
        :param extend: 断点位置单侧搜索区域延伸长度, defaults to 350
        :type extend: int, optional
        :param min_length: 左右断点最小长度, 随意性设置, defaults to 300
        :type min_length: int, optional
        :return: 支持del的reads数目, reads比对信息
        :rtype: (int, pysam.AlignedSegment)
        """
        assert break_point_1 < break_point_2 - min_length
        support_reads = []
        support_reads_segment = []
        left_candidate_reads = {}
        right_candidate_reads = {}
        left_region_reads = self.bam.fetch(chromosome, break_point_1 - extend, break_point_1)
        right_region_reads = self.bam.fetch(chromosome, break_point_2 - 1, break_point_2 + extend)
        reads_bp1 = self.bam.fetch(chromosome, break_point_1 - 1, break_point_1)
        reads_bp2 = self.bam.fetch(chromosome, break_point_2 - 1, break_point_2)
        for read in left_region_reads:
            try:
                mate_read = self.bam.mate(read)
            except ValueError:
                continue
            if read.query_name in support_reads or self.read_filter(read) or self.read_filter(mate_read) or mate_read.reference_name != chromosome:
                continue
            aligned_end = read.reference_end
            mate_aligned_start = mate_read.reference_start + 1
            mate_aligned_end = mate_read.reference_end
            if aligned_end <= break_point_1 and mate_aligned_start >= break_point_2:
                if (break_point_2 - break_point_1) / (mate_aligned_start - aligned_end) >= 0.9:
                    support_reads.append(read.query_name)
                    support_reads_segment.extend([read, mate_read])
            elif mate_aligned_end > aligned_end and mate_aligned_end <= break_point_1 and 'S' in mate_read.cigarstring:
                left_candidate_reads[read.query_name] = [read, mate_read]
        for read in right_region_reads:
            try:
                mate_read = self.bam.mate(read)
            except ValueError:
                continue
            if read.query_name in support_reads or self.read_filter(read) or self.read_filter(mate_read) or mate_read.reference_name != chromosome:
                continue
            aligned_start = read.reference_start + 1
            mate_aligned_start = mate_read.reference_start + 1
            mate_aligned_end = mate_read.reference_end
            if aligned_start >= break_point_2 and mate_aligned_end <= break_point_1:
                if (break_point_2 - break_point_1) / (aligned_start - mate_aligned_end) >= 0.9:
                    support_reads.append(read.query_name)
                    support_reads_segment.extend([read, mate_read])
            elif mate_aligned_start < aligned_start and mate_aligned_start >= break_point_2 and 'S' in mate_read.cigarstring:
                right_candidate_reads[read.query_name] = [read, mate_read]
        for read in reads_bp1:
            if not ((read.is_secondary or read.is_supplementary) and ('S' in read.cigarstring or 'H' in read.cigarstring)):
                continue
            if read.query_name in support_reads:
                support_reads_segment.append(read)
            elif read.query_name in right_candidate_reads:
                support_reads.append(read.query_name)
                support_reads_segment.append(read)
                support_reads_segment.extend(right_candidate_reads[read.query_name])
        for read in reads_bp2:
            if not ((read.is_secondary or read.is_supplementary) and ('S' in read.cigarstring or 'H' in read.cigarstring)):
                continue
            if read.query_name in support_reads:
                support_reads_segment.append(read)
            elif read.query_name in left_candidate_reads:
                support_reads.append(read.query_name)
                support_reads_segment.append(read)
                support_reads_segment.extend(left_candidate_reads[read.query_name])
        return support_reads, support_reads_segment


    def write_bam(self, reads, out_file):
        """输出bam

        :param reads: 输出的reads信息
        :type reads: pysam.AlignedSegment
        :param out_file: 输出bam文件名
        :type out_file: str
        """
        with AlignmentFile(out_file, 'wb', template=self.bam) as out:
            for read in reads:
                out.write(read)

    @staticmethod
    def write_fq_name(fq_names, out_file):
        """输出fq中reads名

        :param fq_names: reads名列表
        :type fq_names: list
        :param out_file: 输出文件名
        :type out_file: str
        """
        with open(out_file, 'w', encoding='utf-8') as out:
            for fq_name in fq_names:
                out.write(f'{fq_name}\n')
    
    @staticmethod
    def ratio_to_copy_number(ratio: float, support_reads: int, mean_reads: int, thres_normal: float, thres_hom_del: float, thres_support_reads: int, thres_mean_reads: int):
        """ratio值转化为拷贝数

        :param ratio: 根据ratio值判断拷贝数
        :type ratio: float
        :param support_reads: 支持cnv的reads数目
        :type support_reads: int
        :param mean_reads: 断点侧翼平均reads数目
        :type mean_reads: int
        :param thres_normal: 正常拷贝数ratio阈值
        :type thres_normal: float
        :param thres_hom_del: 纯合缺失ratio阈值
        :type thres_hom_del: float
        :param thres_support_reads: 支持cnv的reads数目阈值
        :type thres_support_reads: int
        :param thres_mean_reads: 断点侧翼平均reads数目阈值
        :type thres_mean_reads: int
        :return: 拷贝数分析结果
        :rtype: str
        """
        assert thres_normal <= thres_hom_del
        if ratio is None or mean_reads < thres_mean_reads:
            copy_num = 'unknown'
        else:
            if ratio < thres_normal and support_reads < thres_support_reads:
                copy_num = '2'
            elif ratio >= thres_hom_del and support_reads >= thres_support_reads:
                copy_num = '0'
            elif ratio >= thres_normal and support_reads >= thres_support_reads:
                copy_num = '1'
            else:
                copy_num = 'unknown'
        return copy_num
    
    @classmethod
    def run(cls, args) -> pd.DataFrame:
        """运行

        :param args: 输入参数
        :type args: namedtuple|argparse.Namespace
        :return: 支持cnv的reads统计结果
        :rtype: pd.DataFrame
        """
        work_dir, sample, bam, chromosome, bp_left, bp_right, cnv_name, thres_normal, thres_hom_del, thres_support_reads, thres_mean_reads, extend, min_len = args.work_dir, args.sample, args.bam, args.chrom, args.bp_left, args.bp_right, args.cnv_name, args.thres_normal, args.thres_hom_del, args.thres_support_reads, args.thres_mean_reads, args.extend, args.min_len
        out_reads_name = f'{work_dir}/{sample}.{cnv_name}.reads.list'
        out_bam = f'{work_dir}/{sample}.{cnv_name}.bam'
        bh = cls(bam)
        support_reads, support_reads_segment = bh.find_del_support_reads(chromosome, bp_left, bp_right, extend, min_len)
        support_reads_count = len(set(support_reads))
        mean_reads_count = (bh.cal_mean_depth(chromosome, bp_left - extend, bp_left - 1) + bh.cal_mean_depth(chromosome, bp_right, bp_right + extend)) * 0.5
        try:
            ratio = support_reads_count / float(mean_reads_count)
        except ZeroDivisionError:
            ratio = None
        copy_num = cls.ratio_to_copy_number(ratio, support_reads_count, mean_reads_count, thres_normal, thres_hom_del, thres_support_reads, thres_mean_reads)
        df = pd.DataFrame([{'sample': sample, 'mean_reads_count': mean_reads_count, 'support_reads_count': support_reads_count, 'ratio': ratio, 'cnv_name': cnv_name, 'copy_num': copy_num}])
        bh.write_fq_name(set(support_reads), out_reads_name)
        bh.write_bam(support_reads_segment, out_bam)
        bh.close()
        return df
