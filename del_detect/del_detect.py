from argparse import ArgumentParser
import subprocess
from io import StringIO
import pandas as pd
from interval3 import Interval


class DelDetect:
    def __init__(self, info):
        self.info = pd.read_csv(info, sep='\t')

    def detected(self, df):
        result = []
        for cnv in self.info.itertuples(False):
            size_cnv = cnv.stop - cnv.start
            dic = cnv._asdict()
            for row in df.itertuples(False):
                if cnv.chrom != row.CHROM:
                    continue
                overlap = Interval(cnv.start, cnv.stop) & Interval(row.POS, row.END)
                size_overlap = overlap.upper_bound - overlap.lower_bound
                size_del = row.END - row.POS
                overlap_rate_aim = size_overlap / size_cnv
                overlap_rate_self = size_overlap / size_del
                # 重叠50%以上信息保留用于后续统计
                if overlap_rate_aim >= 0.5 and overlap_rate_self >= 0.5:
                    dic['sample_id'] = row.SAMPLE
                    dic['pe_support_count'] = row.PE
                    dic['sr_support_count'] = row.SR
                    dic['overlap_rate_aim'] = overlap_rate_aim
                    dic['overlap_rate_self'] = overlap_rate_self
                    result.append(dic)
        result_df = pd.DataFrame(result)
        columns = self.info.columns.tolist()[:-3] + ['sample_id', 'pe_support_count', 'sr_support_count', 'overlap_rate_aim', 'overlap_rate_self']
        if result_df.empty:
            result_df = pd.DataFrame(columns=columns)
        else:
            result_df = result_df[columns]
        return result_df
    
    @staticmethod
    def yes_or_no(cnv, pe, sr, rate_aim, rate_self, dic):
        return True if pe >= dic[cnv]['pe_threshold'] and sr >= dic[cnv]['sr_threshold'] and rate_aim >= dic[cnv]['overlap_rate'] and rate_self >= dic[cnv]['overlap_rate'] else False
    
    def filter_cnv(self, df):
        info_dic = self.info.set_index('cnv_name').to_dict(orient='index')
        if not df.empty:
            df = df[df.apply(lambda x: self.yes_or_no(x.cnv_name, x.pe_support_count, x.sr_support_count, x.overlap_rate_aim, x.overlap_rate_self, info_dic), axis=1)]
        return df
        
    @staticmethod
    def read_vcf(vcf_file):
        grep = 'zgrep' if vcf_file.endswith('.gz') else 'grep'
        df = pd.read_csv(StringIO(subprocess.getoutput(f'{grep} -v ^## {vcf_file}')), sep='\t', dtype=str)
        df['SAMPLE'] = df.columns[-1]
        df['SVTYPE'] = df['INFO'].str.extract('SVTYPE=(.*?);')
        df['END'] = df['INFO'].str.extract('END=(.*?);')
        df['PE'] = df['INFO'].str.extract(';PE=(.*?);')
        df['SR'] = df['INFO'].str.extract('SR=(.*?)$')
        df = df[(df['SVTYPE'] == 'DEL') & (~df['#CHROM'].str.contains('chrM'))].copy()
        df.rename(columns={'#CHROM': 'CHROM'}, inplace=True)
        df['POS'] = df['POS'].astype(int)
        df['END'] = df['END'].astype(int)
        df.fillna(0, inplace=True)
        df['PE'] = df['PE'].astype(int)
        df['SR'] = df['SR'].astype(int)
        columns = ['SAMPLE', 'CHROM', 'POS', 'END', 'SVTYPE', 'PE', 'SR']
        return df[columns]

    @classmethod
    def run(cls, args):
        dd = cls(args.info)
        df = dd.read_vcf(args.vcf)
        df = dd.detected(df)
        df.to_csv(f'{args.out}.stat.tsv', index=False, sep='\t')
        df = dd.filter_cnv(df)
        df.to_csv(args.out, index=False, sep='\t')


def main():
    parser = ArgumentParser()
    parser.add_argument("-vcf", help="lumpy vcf")
    parser.add_argument("-info", help="del cnv info")
    parser.add_argument("-out", help="output file")
    parsed_args = parser.parse_args()
    DelDetect.run(parsed_args)


if __name__ == "__main__":
    main()
