from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy.signal import find_peaks
from argparse import ArgumentParser
from Bio.SeqIO import AbiIO
import numpy as np
import logging
import os
import glob
import configparser
import pandas as pd
from pyecharts import options
from collections import defaultdict
from pyecharts.charts import Line, Scatter


logger = logging.getLogger()
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
handler.setFormatter(formatter)
logger.addHandler(handler)


standard_size = np.array([20, 30, 40, 60, 80, 100, 114, 120, 140, 160, 180, 200, 214, 220, 240, 250, 260, 280, 300, 314, 320, 340, 360, 380, 400, 414, 
                          420,440, 460, 480, 500, 514, 520, 540, 560, 580, 600, 614, 620, 640,660, 680, 700, 714, 720, 740, 760, 780, 800, 820, 840, 
                          850, 860, 880, 900, 920, 940, 960, 980, 1000, 1020, 1040, 1060, 1080, 1100, 1120, 1160, 1200])
slope_ratio = {20: 0.70, 30: 0.80, 40: 0.90, 60: 0.90, 80: 0.90, 100: 0.90, 114: 0.90, 120: 0.90, 140: 0.90, 160: 0.90, 180: 0.90, 200: 0.90, 214: 0.90, 220: 0.90, 
               240: 0.90, 250: 0.85, 260: 0.90, 280: 0.90, 300: 0.80, 314: 0.75, 320: 0.90, 340: 0.90, 360: 0.90, 380: 0.90, 400: 0.90, 414: 0.90, 420: 0.90, 
               440: 0.90, 460: 0.90, 480: 0.90, 500: 0.90, 514: 0.90, 520: 0.90, 540: 0.90, 560: 0.90, 580: 0.90, 600: 0.90, 614: 0.90, 620: 0.90, 640: 0.90, 
               660: 0.90, 680: 0.90, 700: 0.90, 714: 0.90, 720: 0.90, 740: 0.90, 760: 0.90, 780: 0.90, 800: 0.90, 820: 0.85, 840: 0.85, 850: 0.85, 860: 0.85, 
               880: 0.85, 900: 0.85, 920: 0.85, 940: 0.85, 960: 0.85, 980: 0.85, 1000: 0.80, 1020: 0.80, 1040: 0.80, 1060: 0.75, 1080: 0.75, 1100: 0.65, 
               1120: 0.60}
start = 900


def capillary_electrophoresis_signal(fsa):
    # ABI信号读取
    abi = AbiIO.AbiIterator(open(fsa, 'rb'))
    record = next(abi)
    liz = record.annotations['abif_raw'][f'DATA105']
    fam = record.annotations['abif_raw'][f'DATA1']
    return liz, fam


def cal_slope(x1, y1, x2, y2):
    return (y2 - y1) / (x2 - x1)


def peak_filter(peaks, properties):
    # prominences filter
    prominences_ratio = 0.58
    prominences_pass = properties['prominences'] / properties['peak_heights'] >= prominences_ratio
    peaks = peaks[prominences_pass]
    properties['peak_heights'] = properties['peak_heights'][prominences_pass]
    # 去除末端杂峰
    drop_first_times = 8
    while drop_first_times:
        peaks_heights = defaultdict(dict)
        ratio_dic = dict()
        height_dic = dict()
        peak_dic = dict()
        index = 0
        x2, y2 = peaks[-1], standard_size[-1]
        x1, y1 = peaks[-2], standard_size[-2]
        slope = cal_slope(x1, y1, x2, y2)
        for i, size in enumerate(reversed(standard_size)):
            for j in range(index + 1, peaks.size + 1):
                if i == 0 or i == 1:
                    peaks_heights[size]['peak'] = peaks[-j]
                    peaks_heights[size]['height'] = properties['peak_heights'][-j]
                    height_dic[size] = properties['peak_heights'][-j]
                    peak_dic[size] = peaks[-j]
                    index += 1
                    break
                x0, y0 = peaks[-j], size
                k = cal_slope(x0, y0, x1, y1)
                ratio = np.min([k, slope]) / np.max([k, slope])
                if ratio >= slope_ratio[size]:
                    peaks_heights[size]['peak'] = peaks[-j]
                    peaks_heights[size]['height'] = properties['peak_heights'][-j]
                    height_dic[size] = properties['peak_heights'][-j]
                    peak_dic[size] = peaks[-j]
                    index += 1
                    slope = k
                    x1 = x0
                    y1 = y0
                    ratio_dic[size] = ratio
                    break
                index += 1
        if len(peaks_heights) < 30:
            peaks = np.delete(peaks, -1)
            properties['peak_heights'] = np.delete(properties['peak_heights'], -1)
            drop_first_times -= 1
        else:
            break
    peaks = np.array([peaks_heights[key]['peak'] for key in peaks_heights][::-1])
    heights = np.array([peaks_heights[key]['height'] for key in peaks_heights][::-1])
    properties['peak_heights'] = heights
    return peaks, properties, ratio_dic, height_dic, peak_dic


def peak_calling(liz, fam):
    liz = liz[start:]
    # 峰图检测
    peak_distance = 10
    peak_distance_standard = 23
    low_initial_standard = 100
    low_initial = 30
    up_limit = 20000
    peaks_num = standard_size.size + 9
    withd_max = 40
    # 内标动态阈值
    peaks, properties = find_peaks(liz, height=(low_initial_standard, up_limit), distance=peak_distance_standard, width=(None, withd_max))
    if peaks.size > peaks_num:
        height = sorted(properties['peak_heights'], reverse=True)
        peaks, properties = find_peaks(liz, height=(height[peaks_num] + 1, up_limit), distance=peak_distance_standard, width=(None, withd_max))
    peaks = peaks + start
    if peaks.size >= standard_size.size:
        peaks, properties, ratio_dic, height_dic, peak_dic = peak_filter(peaks, properties)
    else:
        ratio_dic = dict()
        height_dic = dict()
        peak_dic = dict()
    peaks_, properties_ = find_peaks(fam, height=(low_initial, None), distance=peak_distance)
    return peaks, properties, peaks_, properties_, ratio_dic, height_dic, peak_dic


def lowess_fit(x, y):
    # lowess 拟合
    y_ = lowess(y, x, frac=0.1, return_sorted=False)
    u = np.sum(np.square(y - y_))
    v = np.sum(np.square(y - np.mean(y)))
    r_square = np.round(1 - u / v, 4)
    return y_, r_square


def linear_predict(x1, y1, x2, y2, x0):
    # 线性预测
    k = (y2 - y1) / (x2 - x1)
    b = y1 - k * x1
    y0 = k * x0 + b
    return np.round(y0, 1)


def find_two_nearest(x, x0):
    # 2近邻
    x1_index, x2_index = np.argsort(np.abs(x - x0))[:2]
    return x1_index, x2_index


def make_markpoint(peaks, peaks_, heights_):
    # 标注
    markpoint = list()
    for i, peak in enumerate(peaks_):
        if peak < peaks[0]:
            continue
        p1, p2 = find_two_nearest(peaks, peak)
        size = linear_predict(peaks[p1], standard_size[-len(peaks):][p1], peaks[p2], standard_size[-len(peaks):][p2], peak)
        markpoint.append(([peak, heights_[i]], size))
    return markpoint


def plot_signal(liz, fam, title, peaks, properties, out):
    # 信号画图
    line = Line(init_opts=options.InitOpts(bg_color='white'))
    line.set_global_opts(datazoom_opts=options.DataZoomOpts(is_show=True), 
                        toolbox_opts=options.ToolboxOpts(is_show=True),
                        xaxis_opts=options.AxisOpts(name='pos', name_location='center', name_gap=15, type_='value'), 
                        yaxis_opts=options.AxisOpts(name='intensity', name_location='center', name_gap=50, min_=0, type_='value'),
                        title_opts=options.TitleOpts(title=title),
                        tooltip_opts=options.TooltipOpts(position='top', is_always_show_content=True)
                        )
    line.add_xaxis(range(len(liz)))
    line.add_yaxis('LIZ', liz, color='orange', is_symbol_show=False, 
                   markpoint_opts=options.MarkPointOpts(data=[options.MarkPointItem(coord=[str(j), properties['peak_heights'][i]], value=f"({j},{properties['peak_heights'][i]})") for i, j in enumerate(peaks)], 
                                                        symbol_size=10, label_opts=options.LabelOpts(is_show=False))
                    )
    line.add_yaxis('FAM', fam, color='blue', is_symbol_show=False)
    line.render(path=out, template_name='simple_chart.html')


def plot_perfect(liz, fam, title, peaks, properties, peaks_, properties_, out):
    y_, r_square = lowess_fit(peaks, standard_size[-len(peaks):])
    markpoint = make_markpoint(peaks, peaks_, properties_['peak_heights'])
    # 完整版图
    scatter = Scatter(init_opts=options.InitOpts(bg_color='white'))
    scatter.set_global_opts(datazoom_opts=options.DataZoomOpts(is_show=True),
                            toolbox_opts=options.ToolboxOpts(is_show=True),
                            xaxis_opts=options.AxisOpts(name='pos', name_location='center', name_gap=15),
                            yaxis_opts=options.AxisOpts(name='intensity', name_location='center', name_gap=50, min_=0),
                            title_opts=options.TitleOpts(title=title, subtitle=f'R^2={r_square}'),
                            tooltip_opts=options.TooltipOpts(position='top', is_always_show_content=True)
                           )
    scatter.extend_axis(yaxis=options.AxisOpts(type_='value', position='right', name='size', name_location='center', name_gap=35))
    scatter.add_xaxis(xaxis_data=peaks.astype(str))
    scatter.add_yaxis(series_name='size', y_axis=standard_size.tolist()[-len(peaks):], color='green', symbol='triangle', yaxis_index=1, symbol_size=5, 
                      label_opts=options.LabelOpts(font_size=5, position='top'))
    line = Line()
    line.add_xaxis(range(len(liz)))
    line.add_yaxis('LIZ', liz, color='orange', is_symbol_show=False, 
                   markpoint_opts=options.MarkPointOpts(data=[options.MarkPointItem(coord=[str(j), properties['peak_heights'][i]], value=f"({j},{properties['peak_heights'][i]})") for i, j in enumerate(peaks)], 
                                                        symbol_size=10, label_opts=options.LabelOpts(is_show=False))
                    )
    line.add_yaxis('FAM', fam, color='blue', is_symbol_show=False, 
                   markpoint_opts=options.MarkPointOpts(data=[options.MarkPointItem(coord=[str(coord[0]), coord[1]], value=value) for coord, value in markpoint], 
                                                        symbol_size=5, label_opts=options.LabelOpts(is_show=False))
                    )
    line.add_xaxis(peaks.astype(str))
    line.add_yaxis('fit', y_, is_symbol_show=False, yaxis_index=1, color='green')
    scatter.overlap(line).render(path=out, template_name='simple_chart.html')


def run(args):
    ratio_list = list()
    height_list = list()
    peak_list = list()
    for fsa in glob.glob(os.path.join(args.work_dir, '*.fsa')) + glob.glob(os.path.join(args.work_dir, '*.ab1')):
        sample = os.path.splitext(os.path.basename(fsa))[0]
        logger.info(f'{sample}画图')
        liz, fam = capillary_electrophoresis_signal(fsa)
        baseline = np.median(fam)
        if baseline < 0:
            fam = np.array(fam) - baseline
        peaks, properties, peaks_, properties_, ratio_dic, height_dic, peak_dic = peak_calling(liz, fam)
        ratio_dic['sample'] = sample
        ratio_list.append(ratio_dic)
        height_dic['sample'] = sample
        height_list.append(height_dic)
        peak_dic['sample'] = sample
        peak_list.append(peak_dic)
        if peaks.size < standard_size.size - 1:
            plot_signal(liz, fam, sample, peaks, properties, os.path.join(args.work_dir, f'{sample}.html'))
        else:
            plot_perfect(liz, fam, sample, peaks, properties, peaks_, properties_, os.path.join(args.work_dir, f'{sample}.html'))
    pd.DataFrame(ratio_list).to_excel(f'{args.work_dir}/size_slope_ratio.xlsx', index=False)
    pd.DataFrame(height_list).to_excel(f'{args.work_dir}/size_height.xlsx', index=False)
    pd.DataFrame(peak_list).to_excel(f'{args.work_dir}/size_peak.xlsx', index=False)


def main():
    parser = ArgumentParser()
    parser.add_argument('-work_dir', help='work dir')
    args = parser.parse_args()
    logger.info('分析开始!')
    if args.work_dir is None:
        config = configparser.ConfigParser()
        config.read(f'{os.path.dirname(os.path.realpath(__file__))}/../config.ini')
        args.work_dir = config.get('Config', 'work_dir')
    run(args)
    logger.info('分析完成!')


if __name__ == '__main__':
    main()
