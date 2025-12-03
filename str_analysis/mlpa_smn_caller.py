from sklearn.linear_model import LinearRegression
from collections import defaultdict
from scipy.signal import find_peaks
from argparse import ArgumentParser
from interval3 import Interval
from Bio.SeqIO import AbiIO
import pandas as pd
import numpy as np
import logging
import sys
import os
import glob
import configparser
from PIL import Image
from io import BytesIO
from matplotlib import gridspec
import matplotlib.pyplot as plt


logger = logging.getLogger()
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setFormatter(formatter)
logger.addHandler(handler)


# 已知的探针/内标(ROX500)信息和先验的分布/阈值
standard_size = np.array([35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500])
probe_name = ['Ref_1', 'Ref_7', 'Ref_2', 'SMN2_7', 'Ref_8', 'SMN2_8', 'SMN1_7', 'Ref_9', 'SMN1_8', 'Ref_3', 'Ref_4', 'Ref_5', 'Ref_6']
probe_name_order = ['Ref_1', 'Ref_2', 'Ref_3', 'SMN1_7', 'SMN1_8', 'Ref_4', 'Ref_5', 'Ref_6', 'SMN2_7', 'SMN2_8', 'Ref_7', 'Ref_8', 'Ref_9']
probe_num = len(probe_name)
probe_name_ref = [n for n in probe_name if n.startswith('Ref')]
probe_name_test = [n for n in probe_name_order if n.startswith('SMN')]
probe_size = np.array([109, 118, 122, 132, 138, 148, 162, 170, 178, 188, 198, 209, 218])
size_deviation = {'Ref_1': Interval(-4, 4, closed=False), 'Ref_7': Interval(-4, 2, closed=False), 'Ref_2': Interval(-2, 4, closed=False),
                  'SMN2_7': Interval(-4, 3, closed=False), 'Ref_8': Interval(-3, 4, closed=False), 'SMN2_8': Interval(-4, 4, closed=False),
                  'SMN1_7': Interval(-4, 4, closed=False), 'Ref_9': Interval(-4, 3, closed=False), 'SMN1_8': Interval(-3, 4, closed=False),
                  'Ref_3': Interval(-4, 4, closed=False), 'Ref_4': Interval(-4, 4, closed=False), 'Ref_5': Interval(-4, 4, closed=False),
                  'Ref_6': Interval(-4, 4, closed=False)}
r_square_threshold = 0.997
probe_detect_fail_count_threshold = 3
loss = 0.7
gain = 1.3
ref_peaks_num = 21
test_peaks_num = 20
interference_peak_ratio = 0.25
peak_distance = 40
peak_height_min = 80
standard_data = 4
dye = 'PET'
dye_color = 'red'


def global_params(standard):
    # 不同的内标设置不同的参数
    global standard_size, probe_size, size_deviation, ref_peaks_num, test_peaks_num, interference_peak_ratio, peak_distance, peak_height_min, standard_data, dye, dye_color
    logger.info(f'内标为{standard}')
    if standard == 'LIZ500':
        standard_size = np.array([35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350, 400, 450, 490, 500])
        probe_size = np.array([109, 118, 122, 132, 138, 148, 162, 170, 178, 188, 198, 209, 218])
        size_deviation = {'Ref_1': Interval(-4, 4, closed=False), 'Ref_7': Interval(-4, 2, closed=False), 'Ref_2': Interval(-2, 4, closed=False),
                  'SMN2_7': Interval(-4, 3, closed=False), 'Ref_8': Interval(-3, 4, closed=False), 'SMN2_8': Interval(-4, 4, closed=False),
                  'SMN1_7': Interval(-4, 4, closed=False), 'Ref_9': Interval(-4, 3, closed=False), 'SMN1_8': Interval(-3, 4, closed=False),
                  'Ref_3': Interval(-4, 4, closed=False), 'Ref_4': Interval(-4, 4, closed=False), 'Ref_5': Interval(-4, 4, closed=False),
                  'Ref_6': Interval(-4, 4, closed=False)}
        ref_peaks_num = 21
        test_peaks_num = 20
        interference_peak_ratio = 0.25
        peak_distance = 40
        peak_height_min = 80
        standard_data = 105
        dye = 'LIZ'
        dye_color = 'orange'
    elif standard == 'ROX500':
        pass
    else:
        raise Exception('内标必须为ROX500或LIZ500')


def capillary_electrophoresis_signal(fsa):
    # ABI信号读取
    abi = AbiIO.AbiIterator(open(fsa, 'rb'))
    record = next(abi)
    df = pd.DataFrame()
    df['FAM'] = record.annotations['abif_raw'][f'DATA1']
    df[f'{dye}'] = record.annotations['abif_raw'][f'DATA{standard_data}']
    return df


def peak_calling(df):
    # 峰图检测
    ref_kwargs = {'height': peak_height_min, 'distance': peak_distance}
    test_kwargs = {'height': peak_height_min, 'distance': peak_distance}
    # 动态阈值
    peaks_ref, properties_ref = find_peaks(df[f'{dye}'].values, **ref_kwargs)
    if peaks_ref.size > ref_peaks_num:
        height = sorted(properties_ref['peak_heights'], reverse=True)
        peaks_ref, properties_ref = find_peaks(df[f'{dye}'].values, height=height[ref_peaks_num] + 1, distance=peak_distance)
        # 过滤杂峰
        index = np.argmin(properties_ref['peak_heights'][-standard_size.size:]) - standard_size.size
        if np.min(properties_ref['peak_heights'][-standard_size.size:]) / np.median(properties_ref['peak_heights'][-standard_size.size:]) < interference_peak_ratio:
            peaks_ref = np.delete(peaks_ref, index)
            properties_ref['peak_heights'] = np.delete(properties_ref['peak_heights'], index)
    peaks_test, properties_test = find_peaks(df['FAM'].values, **test_kwargs)
    if peaks_test.size > test_peaks_num:
        height = sorted(properties_test['peak_heights'], reverse=True)
        peaks_test, properties_test = find_peaks(df['FAM'].values, height=height[test_peaks_num] + 1, distance=peak_distance)
    return peaks_ref, properties_ref, peaks_test, properties_test


def fit_sizing(peaks_ref, peaks_test):
    # 内标片段大小拟合,右对齐
    lr = LinearRegression()
    peaks = peaks_ref[-standard_size.size:].reshape(-1, 1)
    lr.fit(peaks, standard_size)
    score = lr.score(peaks, standard_size)
    size_fit = lr.predict(peaks)
    size_test = lr.predict(peaks_test.reshape(-1, 1))
    return peaks, np.round(score, 4), np.round(size_fit, 1), np.round(size_test, 1)


def linear_regression(actual_size):
    # 线性回归测试探针片段大小
    x = probe_size[actual_size != 0].reshape(-1, 1)
    y = actual_size[actual_size != 0]
    lr = LinearRegression()
    lr.fit(x, y)
    score = lr.score(x, y)
    y_ = lr.predict(x)
    return x, y_, np.round(score, 4)


def determine_test_probe(peaks_test, properties_test, size_test):
    # 测试探针识别,杂峰过滤
    probe_detect = {p: {'size': 0, 'pos': 0, 'height': 0} for p in probe_name}
    for i, size in enumerate(size_test):
        deviation = size - probe_size
        absolute_deviation = np.abs(deviation)
        index = np.argmin(absolute_deviation)
        abs_deviation = absolute_deviation[index]
        probe = probe_name[index]
        if deviation[index] in size_deviation[probe]:
            if 'abs_deviation' in probe_detect[probe].keys():
                if properties_test['peak_heights'][i] > probe_detect[probe]['height']:
                    probe_detect[probe]['abs_deviation'] = abs_deviation
                    probe_detect[probe]['size'] = size
                    probe_detect[probe]['pos'] = peaks_test[i]
                    probe_detect[probe]['height'] = properties_test['peak_heights'][i]
            else:
                probe_detect[probe]['abs_deviation'] = abs_deviation
                probe_detect[probe]['size'] = size
                probe_detect[probe]['pos'] = peaks_test[i]
                probe_detect[probe]['height'] = properties_test['peak_heights'][i]
    return probe_detect


def plot_signal(sample, df, peaks_ref, properties_ref, peaks_test, properties_test):
    # 信号峰画图
    fig = plt.figure(figsize=(16, 8))
    gs = gridspec.GridSpec(2, 1)
    ax0 = fig.add_subplot(gs[0, :])
    ax0.set_title(sample, {'fontweight': 'bold'})
    ax0.plot(df.index.values, df[f'{dye}'].values, color=f'{dye_color}', label=f'{dye}', linewidth=0.8)
    ax0.scatter(peaks_ref, properties_ref['peak_heights'], marker='*', color='black', label='peak')
    ax0.set_xlabel('pos')
    ax0.set_ylabel('intensity')
    ax0.legend()
    ax1 = fig.add_subplot(gs[1, :])
    ax1.plot(df.index.values, df['FAM'].values, color='blue', label='FAM', linewidth=0.8)
    ax1.scatter(peaks_test, properties_test['peak_heights'], marker='*', color='black', label='peak')
    ax1.set_xlabel('pos')
    ax1.set_ylabel('intensity')
    ax1.legend()


def plot_detect(sample, df, peaks_ref, properties_ref, peaks_test, properties_test, size_test, probe_detect, peaks, score, size_fit, actual_size, probe_size_fit, actual_size_fit, score_fit):
    # 峰图检测结果画图
    fig = plt.figure(figsize=(16, 16))
    gs = gridspec.GridSpec(3, 2)
    ax0 = fig.add_subplot(gs[0, :])
    ax0.set_title(sample, {'fontweight': 'bold'})
    ax0.plot(df.index.values, df[f'{dye}'].values, color=f'{dye_color}', label=f'{dye}', linewidth=0.8)
    ax0.scatter(peaks_ref, properties_ref['peak_heights'], marker='*', color='black', label='peak')
    ax0.set_xlabel('pos')
    ax0.set_ylabel('intensity')
    ax0.legend()
    for i in range(1, standard_size.size + 1):
        ax0.text(peaks_ref[-i], properties_ref['peak_heights'][-i], str(standard_size[-i]), horizontalalignment='center', verticalalignment='bottom')
    ax1 = fig.add_subplot(gs[1, :])
    ax1.plot(df.index.values[:peaks_ref[-1]], df['FAM'].values[:peaks_ref[-1]], color='blue', label='FAM', linewidth=0.8)
    ax1.scatter(peaks_test, properties_test['peak_heights'], marker='*', color='black', label='peak')
    ax1.set_xlabel('pos')
    ax1.set_ylabel('intensity')
    ax1.legend()
    for i in range(peaks_test.size):
        ax1.text(peaks_test[i], properties_test['peak_heights'][i], str(size_test[i]), horizontalalignment='center', verticalalignment='bottom')
    for probe in probe_detect:
        if probe_detect[probe]['size'] == 0:
            continue
        ax1.text(probe_detect[probe]['pos'], probe_detect[probe]['height'] / 2, probe, fontweight='bold', rotation='vertical', horizontalalignment='center', verticalalignment='center', color='green')
    ax2 = fig.add_subplot(gs[2, 0])
    ax2.scatter(peaks, standard_size, marker='*', label='theoretical')
    ax2.plot(peaks, size_fit, linestyle='--', label='fit')
    ax2.set_xlabel('pos')
    ax2.set_ylabel('size')
    ax2.legend(loc='upper left')
    ax2.text(np.median(peaks), np.median(size_fit), f'$R^2$={score}', horizontalalignment='left', verticalalignment='top')
    ax3 = fig.add_subplot(gs[2, 1])
    ax3.scatter(probe_size, actual_size, marker='*', label='actual')
    ax3.plot(probe_size_fit, actual_size_fit, linestyle='--', label='fit')
    ax3.text(np.median(probe_size_fit), np.median(actual_size_fit), f'$R^2$={score_fit}', horizontalalignment='left', verticalalignment='top')
    for i in range(probe_num):
        ax3.text(probe_size[i], actual_size[i], probe_name[i], horizontalalignment='right', verticalalignment='bottom')
    ax3.set_xlabel('theoretical size')
    ax3.set_ylabel('actual size')
    ax3.legend(loc='upper left')


def cal_dosage_quotient(df, all_as_ref):
    # 计算DQ
    if all_as_ref:
        df['is_ref'] = 'yes'
        logger.info('使用批次样本作为参考')
    df_ = df.loc[df['QC'] == 'pass'].copy()
    if np.sum(df_['is_ref'] == 'yes') < 3:
        logger.info('分析失败!可用的参考样本必须>=3例!退出...')
        sys.exit(0)
    result = list()
    for row in df_.itertuples(index=False):
        logger.info(f'{row.sample}计算DQ')
        dic = row._asdict()
        for probe in probe_name:
            for ref in probe_name_ref:
                if probe == ref:
                    continue
                dic[f'{probe}.{ref}.DQ'] = dic[f'{probe}_height'] / dic[f'{ref}_height']
        result.append(dic)
    df_dq = pd.DataFrame(result)
    df_dq.index = df_dq['sample'].tolist()
    ref_samples = df_dq.loc[df_dq['is_ref'] == 'yes', 'sample'].values
    dq_dic = defaultdict(dict)
    for i in df_dq.index:
        logger.info(f'{i}参考DQ校正')
        for j in probe_name:
            dq_ij_list = list()
            mad_ij_list = list()
            for h in ref_samples:
                if i == h:
                    continue
                dq_ihjz = np.array([df_dq.loc[i, f'{j}.{z}.DQ'] / df_dq.loc[h, f'{j}.{z}.DQ'] for z in probe_name_ref if j != z and df_dq.loc[h, f'{j}.{z}.DQ'] != 0])
                if len(dq_ihjz) == 0:
                    continue
                dq_ihj = np.median(dq_ihjz)
                mad_ihj = np.median(np.abs(dq_ihjz - dq_ihj))
                dq_ij_list.append(dq_ihj)
                mad_ij_list.append(mad_ihj)
            dq_ij = np.mean(dq_ij_list)
            mad_ij = np.mean(mad_ij_list)
            var_ij = np.mean(np.square(np.array(dq_ij_list) - dq_ij))
            r = 1.96 * np.sqrt(np.square(1.4826 * mad_ij) + np.square(var_ij))
            dq_dic[i][j] = {'DQ': dq_ij, 'low': np.round(dq_ij - r, 4), 'up': np.round(dq_ij + r, 4)}
    dq_dic_ = defaultdict(dict)
    for probe in probe_name:
        dq = [dq_dic[sample][probe]['DQ'] for sample in ref_samples]
        dq_mean = np.mean(dq)
        r_ = 1.96 * np.mean(np.square(np.array(dq) - dq_mean))
        dq_dic_[probe]['low'] = np.round(dq_mean - r_, 4)
        dq_dic_[probe]['up'] = np.round(dq_mean + r_, 4)
    return dq_dic, dq_dic_


def plot_ratio(dq_dic, dq_dic_, prefix):
    # ratio画图
    imgs = list()
    for sample in dq_dic:
        logger.info(f'{sample}画Ratio图')
        fig, ax = plt.subplots(figsize=(12, 6))
        ax.set_title(sample, {'fontweight': 'bold'})
        ax.set_ylabel('ratio')
        ax.plot(range(1, probe_num + 1), [loss] * probe_num, color='red')
        ax.plot(range(1, probe_num + 1), [gain] * probe_num, color='red')
        ax.scatter(range(1, probe_num + 1), [dq_dic[sample][p]['DQ'] for p in probe_name_order], color='red')
        for j, probe in enumerate(probe_name_order):
            i = j + 1
            background = 'blue' if probe.startswith('SMN') else 'darkviolet'
            ax.fill_between([i - 0.48, i + 0.48], 0, 3, color=background, alpha=0.5)
            ax.fill([i - 0.25, i + 0.25, i + 0.25, i - 0.25], [dq_dic_[probe]['low'], dq_dic_[probe]['low'], dq_dic_[probe]['up'], dq_dic_[probe]['up']], color='green', alpha=0.5)
            if dq_dic[sample][probe]['low'] > gain or dq_dic[sample][probe]['up'] < loss:
                color = 'red'
            elif (dq_dic[sample][probe]['low'] < loss < dq_dic[sample][probe]['up']) or (dq_dic[sample][probe]['low'] < gain < dq_dic[sample][probe]['up']):
                color = 'orange'
            else:
                color = 'black'
            ax.plot([i - 0.25, i + 0.25], [dq_dic[sample][probe]['low'], dq_dic[sample][probe]['low']], color=color)
            ax.plot([i - 0.25, i + 0.25], [dq_dic[sample][probe]['up'], dq_dic[sample][probe]['up']], color=color)
            ax.plot([i, i], [dq_dic[sample][probe]['low'], dq_dic[sample][probe]['up']], color=color)
        ax.set_xticks(range(1, probe_num + 1), labels=probe_name_order)
        buf = BytesIO()
        plt.savefig(buf, format='jpg', dpi=200)
        buf.seek(0)
        img = Image.open(buf)
        imgs.append(img)
        plt.close()
    imgs[0].save(f'{prefix}.ratio.pdf', "PDF", resolution=300.0, save_all=True, append_images=imgs[1:])


def smn_caller(dq_dic, dq_dic_):
    # SMN拷贝数分析
    ref_sample_qc = 'pass'
    abnormal_ref_probe_thres = 2
    abnormal_ref_probe_interval_thres = 8
    high_overlap = 0.6
    high_overlap_ref_probe_thres = 5
    hom_del_ratio = 0.2
    het_del_ratio_low = 0.35
    het_del_ratio_up = 0.65
    normal_ratio_low = 0.8
    normal_ratio_up = 1.3
    for probe in dq_dic_:
        if dq_dic_[probe]['low'] < loss or dq_dic_[probe]['up'] > gain:
            ref_sample_qc = 'ref_sample_fail'
            break
    result = list()
    for sample in dq_dic:
        logger.info(f'{sample}拷贝数分析')
        smn_dic = {'sample': sample}
        qc = list()
        abnormal_ref_probe_count = 0
        abnormal_ref_probe_interval_count = 0
        high_overlap_ref_probe_count = 0
        for probe in probe_name_order:
            if probe.startswith('Ref'):
                if dq_dic[sample][probe]['DQ'] < loss or dq_dic[sample][probe]['DQ'] > gain:
                    abnormal_ref_probe_count += 1
                if dq_dic[sample][probe]['low'] > gain or dq_dic[sample][probe]['up'] < loss or (dq_dic[sample][probe]['low'] < loss < dq_dic[sample][probe]['up']) or (dq_dic[sample][probe]['low'] < gain < dq_dic[sample][probe]['up']):
                    abnormal_ref_probe_interval_count += 1
                overlap = Interval(dq_dic[sample][probe]['low'], dq_dic[sample][probe]['up']) & Interval(loss, gain)
                overlap_len = overlap.upper_bound - overlap.lower_bound
                dq_interval_len = dq_dic[sample][probe]['up'] - dq_dic[sample][probe]['low']
                if dq_interval_len != 0 and overlap_len / dq_interval_len >= high_overlap:
                    high_overlap_ref_probe_count += 1
            else:
                if dq_dic[sample][probe]['DQ'] >= normal_ratio_up:
                    smn_dic[probe] = '>2'
                elif normal_ratio_low < dq_dic[sample][probe]['DQ'] < normal_ratio_up:
                    smn_dic[probe] = '2'
                elif het_del_ratio_up <= dq_dic[sample][probe]['DQ'] <= normal_ratio_low:
                    smn_dic[probe] = '1.5'
                elif het_del_ratio_low < dq_dic[sample][probe]['DQ'] < het_del_ratio_up:
                    smn_dic[probe] = '1'
                elif hom_del_ratio <= dq_dic[sample][probe]['DQ'] <= het_del_ratio_low:
                    smn_dic[probe] = '0.5'
                else:
                    smn_dic[probe] = '0'
                smn_dic[f'{probe}_DQ'] = dq_dic[sample][probe]['DQ']
                smn_dic[f'{probe}_DQ_low'] = dq_dic[sample][probe]['low']
                smn_dic[f'{probe}_DQ_up'] = dq_dic[sample][probe]['up']
        if ref_sample_qc != 'pass':
            qc.append(ref_sample_qc)
        # DQ在[0.7, 1.3]外的参考探针数目>=2
        if abnormal_ref_probe_count >= abnormal_ref_probe_thres:
            qc.append('ref_probe_fail')
        # 参考探针DQ区间异常数目>=8且与[0.7, 1.3]overlap比例超过0.6的数目<=5
        elif abnormal_ref_probe_interval_count >= abnormal_ref_probe_interval_thres and high_overlap_ref_probe_count <= high_overlap_ref_probe_thres:
            qc.append('ref_probe_fail')
        # DQ在[0.7, 1.3]外的参考探针数目>=1且与[0.7, 1.3]overlap比例超过0.6的数目<=6
        elif abnormal_ref_probe_count >= abnormal_ref_probe_thres - 1 and high_overlap_ref_probe_count <= high_overlap_ref_probe_thres + 1:
            qc.append('ref_probe_fail')
        smn_dic['qc'] = ';'.join(qc) if len(qc) else 'pass'
        if len(qc):
            for smn in probe_name_test:
                smn_dic[smn] = '.'
        result.append(smn_dic)
    return pd.DataFrame(result)


def run(args):
    # 工作流程
    global_params(args.standard)
    info = pd.read_csv(args.info, sep='\t')
    result = list()
    imgs = list()
    for row in info.itertuples(index=False):
        logger.info(f'{row.sample}峰图识别')
        dic = {'sample': row.sample, 'is_ref': row.is_ref}
        df = capillary_electrophoresis_signal(row.fsa)
        peaks_ref, properties_ref, peaks_test, properties_test = peak_calling(df)
        qc = list()
        if peaks_ref.size < standard_size.size:
            logger.info(f'{row.sample}内参峰识别失败')
            qc.append('standard_detect_fail')
            plot_signal(row.sample, df, peaks_ref, properties_ref, peaks_test, properties_test)
        else:
            peaks, score, size_fit, size_test = fit_sizing(peaks_ref, peaks_test)
            probe_detect = determine_test_probe(peaks_test, properties_test, size_test)
            actual_size = np.array([probe_detect[probe]['size'] for probe in probe_name])
            probe_size_, actual_size_, score_ = linear_regression(actual_size)
            undetected_probe = list(filter(lambda x: probe_detect[x]['size'] == 0, probe_detect.keys()))
            for probe in probe_name:
                dic[f'{probe}_size'] = probe_detect[probe]['size']
                dic[f'{probe}_height'] = probe_detect[probe]['height']
            dic['r_square_sizing'] = score
            dic['r_square_detect'] = score_
            if score < r_square_threshold:
                logger.info(f'{row.sample}内参峰sizing拟合失败')
                qc.append('standard_fit_sizing_fail')
            if len(set(undetected_probe) & set(probe_name_ref)) > 0 or len(set(undetected_probe) & set(probe_name_test)) > probe_detect_fail_count_threshold:
                logger.info(f'{row.sample}探针峰识别失败')
                qc.append('probe_detect_fail')
            if score_ < r_square_threshold:
                logger.info(f'{row.sample}探针峰detect拟合失败')
                qc.append('probe_fit_sizing_fail')
            plot_detect(row.sample, df, peaks_ref, properties_ref, peaks_test, properties_test, size_test, probe_detect, peaks, score, size_fit, actual_size, probe_size_, actual_size_, score_)
        dic['QC'] = ';'.join(qc) if len(qc) else 'pass'
        result.append(dic)
        buf = BytesIO()
        plt.savefig(buf, format='jpg', dpi=200)
        buf.seek(0)
        img = Image.open(buf)
        imgs.append(img)
        plt.close()
    result_df = pd.DataFrame(result)
    logger.info('输出峰图检测结果')
    result_df.to_excel(f'{args.prefix}.peak.xlsx', index=False)
    imgs[0].save(f'{args.prefix}.peak.pdf', "PDF", resolution=300.0, save_all=True, append_images=imgs[1:])
    dq_dic, ref_sample_dq_dic = cal_dosage_quotient(result_df, args.all_as_ref)
    plot_ratio(dq_dic, ref_sample_dq_dic, args.prefix)
    smn = smn_caller(dq_dic, ref_sample_dq_dic)
    smn_ = pd.merge(result_df, smn, on=['sample'], how='left').fillna('.')
    smn_.loc[smn_['qc'].str.contains('fail'), 'QC'] = smn_.loc[smn_['qc'].str.contains('fail'), 'qc']
    columns = ['sample', 'QC', 'is_ref', 'r_square_sizing', 'r_square_detect'] + probe_name_test + [f'{p}_{q}' for p in probe_name_test for q in ['DQ', 'DQ_low', 'DQ_up']]
    smn_.to_excel(f'{args.prefix}.smn.xlsx', index=False, columns=columns)


def prepare(work_dir):
    # 生成样本路径信息
    logger.info('读取样本信息sample_info.xlsx')
    df = pd.read_excel(os.path.join(work_dir, 'sample_info.xlsx'), dtype=str)
    info = list()
    for i in range(df.shape[0]):
        sample = df.loc[i, 'sample']
        is_ref = df.loc[i, 'is_ref'].lower()
        fsa = os.path.join(work_dir, f'{sample}.fsa')
        ab1 = os.path.join(work_dir, f'{sample}.ab1')
        file_path = fsa if os.path.exists(fsa) else ab1
        if not os.path.exists(ab1) and not os.path.exists(fsa):
            logger.warning(f'未找到{sample}.fsa或{sample}.ab1 !!! 跳过{sample}')
            continue
        info.append({'sample': sample, 'is_ref': is_ref, 'fsa': file_path})
    logger.info('生成样本路径信息')
    pd.DataFrame(info).drop_duplicates(subset=['sample']).to_csv(os.path.join(work_dir, 'sample_info.txt'), sep='\t', index=False)


def generate_sample_info(work_dir):
    # 自动生成样本信息,需要手动设置参考样本
    sample_info = os.path.join(work_dir, 'sample_info.xlsx')
    if not os.path.exists(sample_info):
        logger.info('无样本信息表格,自动生成中')
        samples = [os.path.splitext(os.path.basename(fsa))[0] for fsa in glob.glob(f'{os.path.join(work_dir, "*.fsa")}') + glob.glob(f'{os.path.join(work_dir, "*.ab1")}')]
        df = pd.DataFrame(set(samples), columns=['sample'])
        df['is_ref'] = 'no'
        df.to_excel(sample_info, index=False)
        logger.info('请手动设置样本表格中的参考样本,在下次运行程序前需确保参考样本数目>=3')
        logger.info('样本信息表格生成完成,退出!')
        sys.exit(0)


def main():
    parser = ArgumentParser()
    parser.add_argument('-prefix', help='output file prefix')
    parser.add_argument('-info', help='sample info')
    parser.add_argument('-standard', help='standard ROX500,LIZ500', default='ROX500')
    parser.add_argument('-all_as_ref', help='all sample as ref', action='store_true')
    args = parser.parse_args()
    logger.info('分析开始!')
    if args.prefix is None and args.info is None:
        config = configparser.ConfigParser()
        config.read(f'{os.path.dirname(os.path.realpath(__file__))}/../config.ini')
        work_dir = config.get('Config', 'work_dir')
        try:
            standard = config.get('Config', 'standard')
            args.standard = standard
        except:
            pass
        batch = os.path.basename(work_dir)
        args.prefix = os.path.join(work_dir, batch)
        generate_sample_info(work_dir)
        prepare(work_dir)
        args.info = os.path.join(work_dir, 'sample_info.txt')
    run(args)
    logger.info('分析完成!')


if __name__ == '__main__':
    main()
