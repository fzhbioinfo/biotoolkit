# -*- coding:utf-8 -*-
import os
import pandas as pd
from itertools import count
from argparse import ArgumentParser


def family_role(relation):
    if relation == '女方':
        role = 'Mother'
    elif relation == '男方':
        role = 'Father'
    elif '儿子' in relation or '女儿' in relation or '孩子' in relation:
        role = 'Proband'
    elif '女方' in relation:
        role = 'Maternal'
    elif '男方' in relation:
        role = 'Paternal'
    elif relation == '胚胎':
        role = 'Embryo'
    else:
        raise Exception('undefined family role!')
    return role


def member_gender(relation):
    if relation == '女方':
        gender = 'F'
    elif relation == '男方':
        gender = 'M'
    elif '儿子' in relation or '父亲' in relation or '弟' in relation or '哥' in relation:
        gender = 'M'
    elif '女儿' in relation or '母亲' in relation or '姐' in relation or '妹' in relation:
        gender = 'F'
    else:
        gender = 'unknown'
    return gender


def generate_info(parsed_args):
    df = pd.read_excel(parsed_args.info, dtype={'barcode': str})
    df['家系编号'] = df['家系编号'].fillna(method='ffill')
    df = df[(df['分析进度'] == '已下机') & (df['备注'].isna())].copy()
    df['sample'] = df.apply(lambda x: family_role(x['家系关系']), axis=1)
    df['gender'] = df.apply(lambda x: member_gender(x['家系关系']), axis=1)
    for family in df['家系编号'].unique():
        info = list()
        family_df = df.loc[df['家系编号'] == family].copy()
        family_df.reset_index(drop=True, inplace=True)
        embryo = count(1)
        for index in family_df.index:
            sample = family_df.loc[index, 'sample']
            relation = family_df.loc[index, '家系关系']
            gender = family_df.loc[index, 'gender']
            if sample == 'Embryo':
                sample = f'{sample}{next(embryo)}'
            sampleID = family_df.loc[index, '样本编号']
            barcode = family_df.loc[index, 'barcode']
            barcode = barcode.split('-')
            chips = family_df.loc[index, '上机信息'].split('/')
            lanes = family_df.loc[index, 'lane'].replace('L', '').replace('0', '').split('/')
            for i, chip in enumerate(chips):
                lane = lanes[i].split('-')
                for b in range(int(barcode[0]), int(barcode[-1])+1):
                    for l in range(int(lane[0]), int(lane[-1])+1):
                        fq1 = os.path.join(parsed_args.data_dir, chip, f'L0{l}', f'{chip}_L0{l}_{b}_1.fq.gz')
                        info.append({'sample': sample, 'sampleID': sampleID, 'fq1': fq1, 'relation': relation, 'gender': gender})
        info_df = pd.DataFrame(info)
        info_df.to_csv(f'{os.path.join(parsed_args.out_dir, f"sample.info.{family}.txt")}', sep='\t', index=False)


def main():
    parser = ArgumentParser()
    parser.add_argument('-data_dir', help='fq data dir', required=True)
    parser.add_argument('-out_dir', help='out dir', required=True)
    parser.add_argument('-info', help='sample info', required=True)
    parsed_args = parser.parse_args()
    generate_info(parsed_args)


if __name__ == '__main__':
    main()
