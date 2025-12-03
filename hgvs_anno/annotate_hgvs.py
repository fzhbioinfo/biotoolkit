# -*- coding:utf-8 -*-
from vcf.parser import _Info as VcfInfo, field_counts as vcf_field_counts
from hgvs.exceptions import HGVSError, HGVSParseError, HGVSUsageError
from multiprocessing import Queue, Process, cpu_count
from bioutils.assemblies import make_name_ac_map
from hgvs.assemblymapper import AssemblyMapper
from hgvs.dataproviders.uta import connect
from hgvs.normalizer import Normalizer
from collections import namedtuple
from operator import itemgetter
from hgvs.parser import Parser
from itertools import count
import pandas as pd
import optparse
import logging
import queue
import yaml
import vcf
import re


logger = logging.getLogger()
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
handler = logging.StreamHandler()
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.WARNING)


def yaml_read(yaml_file) -> dict:
    # 读取配置文件
    with open(yaml_file, 'r') as y:
        yaml_dic = yaml.load(y, Loader=yaml.FullLoader)
    return yaml_dic


def remove_trans_acc(trans) -> list:
    # 去除RefSeq转录本版本号: NM_007294.3 -> NM_007294
    trans_no_acc = [re.sub(r'\..*?$', '', tran) for tran in trans]
    return trans_no_acc


def get_sorted_trans_acc(trans) -> list:
    # RefSeq转录本列表若同一转录本按版本号降序排列
    trans = [re.findall(r'(.*?)\.(.*?)$', tran)[0] for tran in trans]
    trans_acc = list()
    for tran in trans:
        try:
            trans_acc.append((tran[0], int(tran[1])))
        except:
            pass
    trans_acc = list(map(lambda x: (x[0], int(x[1])), trans_acc))
    trans_acc.sort(key=itemgetter(0, 1), reverse=True)
    trans_acc_sorted = list(map(lambda x: '.'.join((x[0], str(x[1]))), trans_acc))
    return trans_acc_sorted


def select_trans(trans_related, trans_provided, how) -> list:
    # UTA中相关的转录本与所提供的选择最新的版本, 若未提供或无交集选择查询UTA得到的第一个转录本
    if len(trans_related) == 0:
        return []
    if len(trans_provided) == 0 and how == 'first':
        return [trans_related[0]]
    elif len(trans_provided) == 0 and how != 'first':
        return trans_related
    trans_related_no_acc = remove_trans_acc(trans_related)
    trans_no_acc_overlap = set(trans_related_no_acc) & set(trans_provided)
    if len(trans_no_acc_overlap) == 0 and how == 'first':
        return [trans_related[0]]
    elif len(trans_no_acc_overlap) == 0 and how != 'first':
        return trans_related
    trans_filter = [tran for tran_overlap in trans_no_acc_overlap for tran in trans_related if tran.startswith(tran_overlap + '.')]
    trans_filter = get_sorted_trans_acc(trans_filter)
    trans_selected = list()
    for tran_overlap in trans_no_acc_overlap:
        for tran in trans_filter:
            if tran.startswith(tran_overlap + '.'):
                trans_selected.append(tran)
                break
    return trans_selected


def bgi_anno_select_trans(trans_related, trans_provided, how) -> list:
    # UTA中相关的转录本与所提供的选择最新的版本, 若未提供或无交集选择查询UTA得到的第一个转录本
    if len(trans_related) == 0:
        return []
    if len(trans_provided) == 0 and how == 'first':
        return [trans_related[0]]
    elif len(trans_provided) == 0 and how != 'first':
        return trans_related
    trans_selected = list(set(trans_related) & set(trans_provided))
    # trans_related_no_acc = remove_trans_acc(trans_related)
    # trans_no_acc_overlap = set(trans_related_no_acc) & set(trans_provided)
    # if len(trans_no_acc_overlap) == 0 and how == 'first':
    #     return [trans_related[0]]
    # elif len(trans_no_acc_overlap) == 0 and how != 'first':
    #     return []
    # trans_filter = [tran for tran_overlap in trans_no_acc_overlap for tran in trans_related if tran.startswith(tran_overlap + '.')]
    # trans_filter = get_sorted_trans_acc(trans_filter)
    # trans_selected = list()
    # for tran_overlap in trans_no_acc_overlap:
    #     for tran in trans_filter:
    #         if tran.startswith(tran_overlap + '.'):
    #             trans_selected.append(tran)
    #             break
    return trans_selected


# def generate_chrome_dic(annotation) -> dict:
#     # 染色体对应关系字典, GRCh37: chr1 or 1 -> NC_000001.10, GRCh38: chr1 or 1 -> NC_000001.11
#     chrome_dic = dict()
#     if annotation == 'GRCh37':
#         nc_suffix = ['01.10', '02.11', '03.11', '04.11', '05.9', '06.11', '07.13', '08.10', '09.11', '10.10', '11.9', '12.11', '13.10',
#                      '14.8', '15.9', '16.9', '17.10', '18.9', '19.9', '20.10', '21.8', '22.10', '23.10', '24.9']
#     else:
#         nc_suffix = ['01.11', '02.12', '03.12', '04.12', '05.10', '06.12', '07.14', '08.11', '09.12', '10.11', '11.10', '12.12', '13.11',
#                      '14.9', '15.10', '16.10', '17.11', '18.10', '19.10', '20.11', '21.9', '22.11', '23.11', '24.10']
#     nc = ['NC_0000' + i for i in nc_suffix] + ['NC_012920.1']
#     chromes = [str(j) for j in range(1, 23)] + ['X', 'Y', 'MT']
#     for k, chrome in enumerate(chromes):
#         chrome_dic[chrome] = nc[k]
#     chromes_chr = ['chr' + m if m != 'MT' else 'chrM_NC_012920.1' for m in chromes]
#     for n, chrome in enumerate(chromes_chr):
#         chrome_dic[chrome] = nc[n]
#     chrome_dic['chrMT'] = 'NC_012920.1'
#     return chrome_dic


def generate_chrome_dic(annotation) -> dict:
    # 染色体对应关系字典, GRCh37: chr1 or 1 -> NC_000001.10, GRCh38: chr1 or 1 -> NC_000001.11
    chrome_dic = make_name_ac_map(annotation)
    chromes = [str(j) for j in range(1, 23)] + ['X', 'Y']
    for chrome in chromes:
        chrome_dic['chr'+chrome] = chrome_dic[chrome]
    if 'MT' not in chrome_dic:
        chrome_dic['MT'] = 'NC_012920.1'
    chrome_dic['chrMT'] = chrome_dic['MT']
    chrome_dic['chrM_NC_012920.1'] = chrome_dic['MT']
    return chrome_dic


def annotator(annotation):
    # 连接UTA, 创建Mapper、Parser、Normalizer
    hdp = connect()
    am = AssemblyMapper(hdp, assembly_name=annotation, alt_aln_method='splign', replace_reference=True)
    hp = Parser()
    hn3 = Normalizer(hdp, shuffle_direction=3)
    hn5 = Normalizer(hdp, shuffle_direction=5)
    return am, hp, hn3, hn5, hdp


def get_transcript_strand(opts, hdp, g, tran):
    if opts.rule_3 == 'transcript':
        nc = str(g.split(':')[0])
        tx_exons = hdp.get_tx_exons(tran, nc, 'splign')
        strand = tx_exons[0][4]
        if strand == 1:
            return 3
        elif strand == -1:
            return 5
        else:
            return 3
    else:
        return 3


def judge_var_type(ref, call):
    if ref == '.':
        var_type = 'ins'
    elif call == '.':
        var_type = 'del'
    elif ref != '.' and call != '.' and len(ref) == 1 and len(call) == 1 and ref != call:
        var_type = 'snv'
    elif ref != '.' and call != '.' and (len(ref) > 1 or len(call) > 1) and ref != call:
        var_type = 'delins'
    else:
        var_type = 'ref'
    return var_type


def generate_g(record_parser, chrome_dic):
    # vcf record 转成 g
    chrome = chrome_dic[record_parser.chrome]
    if record_parser.chrome.__contains__('M'):
        dna = 'm'
    else:
        dna = 'g'
    if record_parser.var_type == 'snv':
        g = chrome + ':' + dna + '.' + str(record_parser.stop) + record_parser.ref + '>' + record_parser.call
    elif record_parser.var_type == 'del':
        start = record_parser.start + 1
        if start == record_parser.stop:
            g = chrome + ':' + dna + '.' + str(record_parser.stop) + 'del'
        else:
            g = chrome + ':' + dna + '.' + str(start) + '_' + str(record_parser.stop) + 'del'
    elif record_parser.var_type == 'delins':
        start = record_parser.start + 1
        g = chrome + ':' + dna + '.' + str(start) + '_' + str(record_parser.stop) + 'delins' + record_parser.call
    elif record_parser.var_type == 'ref':
        if len(record_parser.ref) == 1 and len(record_parser.call) == 1:
            g = chrome + ':' + dna + '.' + str(record_parser.stop) + record_parser.ref + '>' + record_parser.call
        else:
            g = '.'
    else:
        g = chrome + ':' + dna + '.' + str(record_parser.start) + '_' + str(record_parser.start+1) + 'ins' + record_parser.call
    return g


def create_annotate_result(file_out):
    # 创建结果文件句柄
    fp = open(file_out, 'w')
    fp.write('#Chr' + '\t')
    fp.write('Start' + '\t')
    fp.write('Stop' + '\t')
    fp.write('Ref' + '\t')
    fp.write('Call' + '\t')
    fp.write('VarType' + '\t')
    fp.write('cHGVS' + '\t')
    fp.write('cHGVS Normalise' + '\t')
    fp.write('pHGVS' + '\t')
    fp.write('pHGVS Normalise' + '\n')
    return fp


def write_annotate_result(fp, record_annotation):
    # 写结果
    fp.write(record_annotation.chrome + '\t')
    fp.write(str(record_annotation.start) + '\t')
    fp.write(str(record_annotation.stop) + '\t')
    fp.write(record_annotation.ref + '\t')
    fp.write(record_annotation.call + '\t')
    fp.write(record_annotation.var_type + '\t')
    fp.write(record_annotation.cHGVS + '\t')
    fp.write(record_annotation.cHGVS_Normalise + '\t')
    fp.write(record_annotation.pHGVS + '\t')
    fp.write(record_annotation.pHGVS_Normalise + '\n')


RecordAnnotation = namedtuple('RecordAnnotation', ['chrome', 'start', 'stop', 'ref', 'call', 'var_type', 'cHGVS', 'cHGVS_Normalise',
                                                   'pHGVS', 'pHGVS_Normalise'])
VariantRecord = namedtuple('VariantRecord', ['chrome', 'start', 'stop', 'ref', 'call', 'var_type'])


def serial_annotate(opts, trans_provided_no_acc):
    # 串行注释, 生成vcf格式
    am, hp, hn3, hn5, hdp = annotator(opts.annotation)
    chrome_dic = generate_chrome_dic(opts.annotation)
    vcf_reader = vcf.Reader(filename=opts.file_in)
    vcf_reader.infos['HGVS'] = VcfInfo('HGVS', vcf_field_counts['A'], 'String', 'VCF record alleles in HGVS syntax', version=None, source=None)
    vcf_reader.infos['HGVS_Normalise'] = VcfInfo('HGVS_Normalise', vcf_field_counts['A'], 'String', 'VCF record alleles in HGVS syntax (Normalised)', version=None, source=None)
    vcf_writer = vcf.Writer(open(opts.file_out, 'w'), vcf_reader)
    for record in vcf_reader:
        chrome = str(record.CHROM)
        start = record.affected_start
        stop = record.affected_end
        record_hgvs_list = list()
        record_hgvs_normalise_list = list()
        for alt in record.ALT:
            hgvs_list = list()
            hgvs_normalise_list = list()
            if record.is_snp:
                var_type = 'snv'
                ref = record.REF
                call = str(alt)
            else:
                if len(record.REF) == 1 and len(str(alt)) > 1:
                    var_type = 'ins'
                    ref = '.'
                    call = str(alt)[1:]
                elif len(record.REF) > 1 and len(str(alt)) == 1:
                    var_type = 'del'
                    ref = record.REF[1:]
                    call = '.'
                else:
                    var_type = 'delins'
                    if record.REF[0] == str(alt)[0]:
                        ref = record.REF[1:]
                        call = str(alt)[1:]
                    else:
                        ref = record.REF
                        call = str(alt)
                        start = record.affected_start - 1
            record_parser = VariantRecord(chrome, start, stop, ref, call, var_type)
            g = generate_g(record_parser, chrome_dic)
            try:
                g_parser = hp.parse_hgvs_variant(g)
                g_normalise_3 = hn3.normalize(g_parser)
                g_normalise_5 = hn5.normalize(g_parser)
                trans_related = am.relevant_transcripts(g_parser)
            except (HGVSParseError, HGVSError, HGVSUsageError) as e:
                error = str(e)
                logging.error('{chrome} {start} {stop} {ref} {call} {g} annotate error. {error}.'.format(**locals()))
                record_hgvs_list.append('.|.|.')
                record_hgvs_normalise_list.append('.|.|.')
                continue
            trans = select_trans(trans_related, trans_provided_no_acc, opts.how)
            if len(trans) == 0:
                logging.warning('{chrome} {start} {stop} {ref} {call} {g} no related transcripts in UTA.'.format(**locals()))
                record_hgvs_list.append(g+'|.|.')
                record_hgvs_normalise_list.append(str(g_normalise_3)+'|.|.')
                continue
            for tran in trans:
                try:
                    t = am.g_to_t(g_parser, tran)
                    strand = get_transcript_strand(opts, hdp, g, tran)
                    if strand == 3:
                        g_normalise = g_normalise_3
                    else:
                        g_normalise = g_normalise_5
                    t_normalise = am.g_to_t(g_normalise, tran)
                    p = am.t_to_p(t)
                    p_normalise = am.t_to_p(t_normalise)
                    hgvs_ = '|'.join([g, str(t), str(p)])
                    hgvs_normalise = '|'.join([str(g_normalise), str(t_normalise), str(p_normalise)])
                except (HGVSError, HGVSUsageError, NotImplementedError, IndexError) as e:
                    error = str(e)
                    logging.error('{chrome} {start} {stop} {ref} {call} {tran} {g} annotate error. {error}.'.format(**locals()))
                    hgvs_ = '|'.join([g, '.', '.'])
                    hgvs_normalise = '|'.join([str(g_normalise_3), '.', '.'])
                hgvs_list.append(hgvs_)
                hgvs_normalise_list.append(hgvs_normalise)
            hgvs_alt = '/'.join(hgvs_list)
            hgvs_normalise_alt = '/'.join(hgvs_normalise_list)
            record_hgvs_list.append(hgvs_alt)
            record_hgvs_normalise_list.append(hgvs_normalise_alt)
        record_hgvs = ','.join(record_hgvs_list)
        record_hgvs_normalise = ','.join(record_hgvs_normalise_list)
        record.add_info('HGVS', record_hgvs)
        record.add_info('HGVS_Normalise', record_hgvs_normalise)
        vcf_writer.write_record(record)
    vcf_writer.close()


def serial_annotate_to_bed(opts, trans_provided_no_acc):
    # 串行注释, 生成bed格式
    fp = create_annotate_result(opts.file_out)
    am, hp, hn3, hn5, hdp = annotator(opts.annotation)
    chrome_dic = generate_chrome_dic(opts.annotation)
    vcf_reader = vcf.Reader(filename=opts.file_in)
    for record in vcf_reader:
        chrome = str(record.CHROM)
        start = record.affected_start
        stop = record.affected_end
        for alt in record.ALT:
            if record.is_snp:
                var_type = 'snv'
                ref = record.REF
                call = str(alt)
            else:
                if len(record.REF) == 1 and len(str(alt)) > 1:
                    var_type = 'ins'
                    ref = '.'
                    call = str(alt)[1:]
                elif len(record.REF) > 1 and len(str(alt)) == 1:
                    var_type = 'del'
                    ref = record.REF[1:]
                    call = '.'
                else:
                    var_type = 'delins'
                    if record.REF[0] == str(alt)[0]:
                        ref = record.REF[1:]
                        call = str(alt)[1:]
                    else:
                        ref = record.REF
                        call = str(alt)
                        start = record.affected_start - 1
            record_parser = VariantRecord(chrome, start, stop, ref, call, var_type)
            g = generate_g(record_parser, chrome_dic)
            try:
                g_parser = hp.parse_hgvs_variant(g)
                g_normalise_3 = hn3.normalize(g_parser)
                g_normalise_5 = hn5.normalize(g_parser)
                trans_related = am.relevant_transcripts(g_parser)
            except (HGVSParseError, HGVSError, HGVSUsageError) as e:
                error = str(e)
                logging.error('{chrome} {start} {stop} {ref} {call} {g} annotate error. {error}.'.format(**locals()))
                continue
            trans = select_trans(trans_related, trans_provided_no_acc, opts.how)
            if len(trans) == 0:
                logging.warning('{chrome} {start} {stop} {ref} {call} {g} no related transcripts in UTA.'.format(**locals()))
            for tran in trans:
                try:
                    t = am.g_to_t(g_parser, tran)
                    strand = get_transcript_strand(opts, hdp, g, tran)
                    if strand == 3:
                        g_normalise = g_normalise_3
                    else:
                        g_normalise = g_normalise_5
                    t_normalise = am.g_to_t(g_normalise, tran)
                    p = am.t_to_p(t)
                    p_normalise = am.t_to_p(t_normalise)
                    record_annotation = RecordAnnotation(chrome, start, stop, ref, call, var_type, str(t), str(t_normalise), str(p), str(p_normalise))
                    write_annotate_result(fp, record_annotation)
                except (HGVSError, HGVSUsageError, NotImplementedError, IndexError) as e:
                    error = str(e)
                    logging.error('{chrome} {start} {stop} {ref} {call} {tran} {g} annotate error. {error}.'.format(**locals()))
    fp.close()


def process_record(records, results, opts, trans_provided):
    # 一个进程创建uta连接， 从records队列里获取数据, 注释结果返回到results队列, 适用输出为vcf格式
    am, hp, hn3, hn5, hdp = annotator(opts.annotation)
    while True:
        try:
            record = records.get(False)
        except queue.Empty:
            continue
        if record != 'END':
            record_id, record_infos = record[0], record[1]
            record_hgvs_list = list()
            record_hgvs_normalise_list = list()
            for record_info in record_infos:
                record_parser = record_info[0]
                g = record_info[1]
                chrome, start, stop, ref, call, var_type = record_parser.chrome, record_parser.start, record_parser.stop, record_parser.ref, record_parser.call, record_parser.var_type
                hgvs_list = list()
                hgvs_normalise_list = list()
                try:
                    g_parser = hp.parse_hgvs_variant(g)
                    g_normalise_3 = hn3.normalize(g_parser)
                    g_normalise_5 = hn5.normalize(g_parser)
                    trans_related = am.relevant_transcripts(g_parser)
                except (HGVSParseError, HGVSError, HGVSUsageError) as e:
                    error = str(e)
                    logging.error('{chrome} {start} {stop} {ref} {call} {g} annotate error. {error}.'.format(**locals()))
                    record_hgvs_list.append('.|.|.')
                    record_hgvs_normalise_list.append('.|.|.')
                    continue
                trans = select_trans(trans_related, trans_provided, opts.how)
                if len(trans) == 0:
                    logging.warning('{chrome} {start} {stop} {ref} {call} {g} no related transcripts in UTA.'.format(**locals()))
                    record_hgvs_list.append(g + '|.|.')
                    record_hgvs_normalise_list.append(str(g_normalise_3) + '|.|.')
                    continue
                for tran in trans:
                    try:
                        t = am.g_to_t(g_parser, tran)
                        strand = get_transcript_strand(opts, hdp, g, tran)
                        if strand == 3:
                            g_normalise = g_normalise_3
                        else:
                            g_normalise = g_normalise_5
                        t_normalise = am.g_to_t(g_normalise, tran)
                        p = am.t_to_p(t)
                        p_normalise = am.t_to_p(t_normalise)
                        hgvs_ = '|'.join([g, str(t), str(p)])
                        hgvs_normalise = '|'.join([str(g_normalise), str(t_normalise), str(p_normalise)])
                    except (HGVSError, HGVSUsageError, NotImplementedError, IndexError) as e:
                        error = str(e)
                        logging.error('{chrome} {start} {stop} {ref} {call} {tran} {g} annotate error. {error}.'.format(**locals()))
                        hgvs_ = '|'.join([g, '.', '.'])
                        hgvs_normalise = '|'.join([str(g_normalise_3), '.', '.'])
                    hgvs_list.append(hgvs_)
                    hgvs_normalise_list.append(hgvs_normalise)
                hgvs_alt = '/'.join(hgvs_list)
                hgvs_normalise_alt = '/'.join(hgvs_normalise_list)
                record_hgvs_list.append(hgvs_alt)
                record_hgvs_normalise_list.append(hgvs_normalise_alt)
            record_hgvs = ','.join(record_hgvs_list)
            record_hgvs_normalise = ','.join(record_hgvs_normalise_list)
            results.put((record_id, record_hgvs, record_hgvs_normalise))
        else:
            records.put('END')
            break


def process_record_to_bed(records, results, opts, trans_provided):
    # 一个进程创建uta连接， 从records队列里获取数据, 注释结果返回到results队列, 适用输出为bed格式
    am, hp, hn3, hn5, hdp = annotator(opts.annotation)
    while True:
        try:
            record = records.get(False)
        except queue.Empty:
            continue
        if record != 'END':
            record_parser = record[0]
            g = record[1]
            chrome, start, stop, ref, call, var_type = record_parser.chrome, record_parser.start, record_parser.stop, record_parser.ref, record_parser.call, record_parser.var_type
            try:
                g_parser = hp.parse_hgvs_variant(g)
                g_normalise_3 = hn3.normalize(g_parser)
                g_normalise_5 = hn5.normalize(g_parser)
                trans_related = am.relevant_transcripts(g_parser)
            except (HGVSParseError, HGVSError, HGVSUsageError) as e:
                error = str(e)
                logging.error('{chrome} {start} {stop} {ref} {call} {g} annotate error. {error}.'.format(**locals()))
                continue
            trans = select_trans(trans_related, trans_provided, opts.how)
            if len(trans) == 0:
                logging.warning('{chrome} {start} {stop} {ref} {call} {g} no related transcripts in UTA.'.format(**locals()))
                continue
            for tran in trans:
                try:
                    t = am.g_to_t(g_parser, tran)
                    strand = get_transcript_strand(opts, hdp, g, tran)
                    if strand == 3:
                        g_normalise = g_normalise_3
                    else:
                        g_normalise = g_normalise_5
                    t_normalise = am.g_to_t(g_normalise, tran)
                    p = am.t_to_p(t)
                    p_normalise = am.t_to_p(t_normalise)
                    results.put(RecordAnnotation(chrome, start, stop, ref, call, var_type, str(t), str(t_normalise), str(p), str(p_normalise)))
                except (HGVSError, HGVSUsageError, NotImplementedError, IndexError) as e:
                    error = str(e)
                    logging.error('{chrome} {start} {stop} {ref} {call} {tran} {g} annotate error. {error}.'.format(**locals()))
        else:
            records.put('END')
            break


def parallel_annotate(opts, trans_provided_no_acc, process_num):
    # 并行注释
    chrome_dic = generate_chrome_dic(opts.annotation)
    # 创建队列, 初始化
    records, results = Queue(100*process_num), Queue()
    input_finished = False
    output_finished = False
    wait_records = dict()
    records_id = count()
    processes = list()
    # 开启多个进程监听队列, 注释
    for i in range(process_num):
        p = Process(target=process_record, args=(records, results, opts, trans_provided_no_acc))
        processes.append(p)
        p.start()
    # 读取vcf信息, 写入新的vcf
    vcf_reader = vcf.Reader(filename=opts.file_in)
    vcf_reader.infos['HGVS'] = VcfInfo('HGVS', vcf_field_counts['A'], 'String', 'VCF record alleles in HGVS syntax', version=None, source=None)
    vcf_reader.infos['HGVS_Normalise'] = VcfInfo('HGVS_Normalise', vcf_field_counts['A'], 'String', 'VCF record alleles in HGVS syntax (Normalised)', version=None, source=None)
    vcf_writer = vcf.Writer(open(opts.file_out, 'w'), vcf_reader)
    while True:
        while not records.full() and not input_finished:
            try:
                record = next(vcf_reader)
                chrome = str(record.CHROM)
                start = record.affected_start
                stop = record.affected_end
                record_id = next(records_id)
                wait_records[record_id] = record
                record_infos = list()
                for alt in record.ALT:
                    if record.is_snp:
                        var_type = 'snv'
                        ref = record.REF
                        call = str(alt)
                    else:
                        if len(record.REF) == 1 and len(str(alt)) > 1:
                            var_type = 'ins'
                            ref = '.'
                            call = str(alt)[1:]
                        elif len(record.REF) > 1 and len(str(alt)) == 1:
                            var_type = 'del'
                            ref = record.REF[1:]
                            call = '.'
                        else:
                            var_type = 'delins'
                            if record.REF[0] == str(alt)[0]:
                                ref = record.REF[1:]
                                call = str(alt)[1:]
                            else:
                                ref = record.REF
                                call = str(alt)
                                start = record.affected_start - 1
                    record_parser = VariantRecord(chrome, start, stop, ref, call, var_type)
                    g = generate_g(record_parser, chrome_dic)
                    record_infos.append((record_parser, g))
                records.put((record_id, record_infos))
            except StopIteration:
                input_finished = True
                records.put('END')
                break
        processes_status = list()
        for p in processes:
            processes_status.append(p.is_alive())
        if True not in processes_status:
            results.put('END')
        while True:
            try:
                result = results.get(False)
            except queue.Empty:
                break
            if result != 'END':
                record_id, record_hgvs, record_hgvs_normalise = result[0], result[1], result[2]
                record_write = wait_records.pop(record_id)
                record_write.add_info('HGVS', record_hgvs)
                record_write.add_info('HGVS_Normalise', record_hgvs_normalise)
                vcf_writer.write_record(record_write)
            else:
                output_finished = True
                break
        if output_finished:
            break
    vcf_writer.close()


def parallel_annotate_to_bed(opts, trans_provided_no_acc, process_num):
    # 并行注释, 输出bed格式结果
    chrome_dic = generate_chrome_dic(opts.annotation)
    fp = create_annotate_result(opts.file_out)
    records, results = Queue(100*process_num), Queue()
    input_finished = False
    output_finished = False
    processes = list()
    for i in range(process_num):
        p = Process(target=process_record_to_bed, args=(records, results, opts, trans_provided_no_acc))
        processes.append(p)
        p.start()
    vcf_reader = vcf.Reader(filename=opts.file_in)
    while True:
        while not records.full() and not input_finished:
            try:
                record = next(vcf_reader)
                chrome = str(record.CHROM)
                start = record.affected_start
                stop = record.affected_end
                for alt in record.ALT:
                    if record.is_snp:
                        var_type = 'snv'
                        ref = record.REF
                        call = str(alt)
                    else:
                        if len(record.REF) == 1 and len(str(alt)) > 1:
                            var_type = 'ins'
                            ref = '.'
                            call = str(alt)[1:]
                        elif len(record.REF) > 1 and len(str(alt)) == 1:
                            var_type = 'del'
                            ref = record.REF[1:]
                            call = '.'
                        else:
                            var_type = 'delins'
                            if record.REF[0] == str(alt)[0]:
                                ref = record.REF[1:]
                                call = str(alt)[1:]
                            else:
                                ref = record.REF
                                call = str(alt)
                                start = record.affected_start - 1
                    record_parser = VariantRecord(chrome, start, stop, ref, call, var_type)
                    g = generate_g(record_parser, chrome_dic)
                    records.put((record_parser, g))
            except StopIteration:
                input_finished = True
                records.put('END')
                break
        processes_status = list()
        for p in processes:
            processes_status.append(p.is_alive())
        if True not in processes_status:
            results.put('END')
        while True:
            try:
                result = results.get(False)
            except queue.Empty:
                break
            if result != 'END':
                write_annotate_result(fp, result)
            else:
                output_finished = True
                break
        if output_finished:
            break
    fp.close()


def bgianno_serial_annotate(opts, trans_provided):
    # 串行注释, 生成bed格式
    fp = create_annotate_result(opts.file_out)
    am, hp, hn3, hn5, hdp = annotator(opts.annotation)
    chrome_dic = generate_chrome_dic(opts.annotation)
    if opts.file_format == 'excel':
        df = pd.read_excel(opts.file_in, skiprows=opts.skip_rows)
    else:
        df = pd.read_csv(opts.file_in, sep='\t', low_memory=False, skiprows=opts.skip_rows, dtype={'#Chr': str})
    df.rename(columns={'#Chr': 'Chr'}, inplace=True)
    for record in df.itertuples():
        chrome = getattr(record, 'Chr')
        start = getattr(record, 'Start')
        stop = getattr(record, 'Stop')
        ref = getattr(record, 'Ref')
        call = getattr(record, 'Call')
        var_type = judge_var_type(ref, call)
        record_parser = VariantRecord(chrome, start, stop, ref, call, var_type)
        g = generate_g(record_parser, chrome_dic)
        try:
            g_parser = hp.parse_hgvs_variant(g)
            g_normalise_3 = hn3.normalize(g_parser)
            g_normalise_5 = hn5.normalize(g_parser)
            trans_related = am.relevant_transcripts(g_parser)
        except (HGVSParseError, HGVSError, HGVSUsageError) as e:
            error = str(e)
            logging.error('{chrome} {start} {stop} {ref} {call} {g} annotate error. {error}.'.format(**locals()))
            continue
        trans = bgi_anno_select_trans(trans_related, trans_provided, opts.how)
        if len(trans) == 0:
            logging.warning('{chrome} {start} {stop} {ref} {call} {g} no related transcripts in UTA.'.format(**locals()))
        for tran in trans:
            try:
                t = am.g_to_t(g_parser, tran)
                strand = get_transcript_strand(opts, hdp, g, tran)
                if strand == 3:
                    g_normalise = g_normalise_3
                else:
                    g_normalise = g_normalise_5
                t_normalise = am.g_to_t(g_normalise, tran)
                p = am.t_to_p(t)
                p_normalise = am.t_to_p(t_normalise)
                record_annotation = RecordAnnotation(chrome, start, stop, ref, call, var_type, str(t), str(t_normalise), str(p), str(p_normalise))
                write_annotate_result(fp, record_annotation)
            except (HGVSError, HGVSUsageError, NotImplementedError, IndexError) as e:
                error = str(e)
                logging.error('{chrome} {start} {stop} {ref} {call} {tran} {g} annotate error. {error}.'.format(**locals()))
    fp.close()


def bgi_anno_process_record(records, results, opts, trans_provided):
    # 一个进程创建uta连接， 从records队列里获取数据, 注释结果返回到results队列, 适用输出为bed格式
    am, hp, hn3, hn5, hdp = annotator(opts.annotation)
    while True:
        try:
            record = records.get(False)
        except queue.Empty:
            continue
        if record != 'END':
            record_parser = record[0]
            g = record[1]
            chrome, start, stop, ref, call, var_type = record_parser.chrome, record_parser.start, record_parser.stop, record_parser.ref, record_parser.call, record_parser.var_type
            try:
                g_parser = hp.parse_hgvs_variant(g)
                g_normalise_3 = hn3.normalize(g_parser)
                g_normalise_5 = hn5.normalize(g_parser)
                trans_related = am.relevant_transcripts(g_parser)
            except (HGVSParseError, HGVSError, HGVSUsageError) as e:
                error = str(e)
                logging.error('{chrome} {start} {stop} {ref} {call} {g} annotate error. {error}.'.format(**locals()))
                continue
            trans = bgi_anno_select_trans(trans_related, trans_provided, opts.how)
            if len(trans) == 0:
                logging.warning('{chrome} {start} {stop} {ref} {call} {g} no related transcripts in UTA.'.format(**locals()))
                continue
            for tran in trans:
                try:
                    t = am.g_to_t(g_parser, tran)
                    strand = get_transcript_strand(opts, hdp, g, tran)
                    if strand == 3:
                        g_normalise = g_normalise_3
                    else:
                        g_normalise = g_normalise_5
                    t_normalise = am.g_to_t(g_normalise, tran)
                    p = am.t_to_p(t)
                    p_normalise = am.t_to_p(t_normalise)
                    results.put(RecordAnnotation(chrome, start, stop, ref, call, var_type, str(t), str(t_normalise), str(p), str(p_normalise)))
                except (HGVSError, HGVSUsageError, NotImplementedError, IndexError) as e:
                    error = str(e)
                    logging.error('{chrome} {start} {stop} {ref} {call} {tran} {g} annotate error. {error}.'.format(**locals()))
        else:
            records.put('END')
            break


def bgianno_parallel_annotate(opts, trans_provided, process_num):
    # 并行注释, 输出bed格式结果
    chrome_dic = generate_chrome_dic(opts.annotation)
    fp = create_annotate_result(opts.file_out)
    records, results = Queue(100*process_num), Queue()
    input_finished = False
    output_finished = False
    processes = list()
    for i in range(process_num):
        p = Process(target=bgi_anno_process_record, args=(records, results, opts, trans_provided))
        processes.append(p)
        p.start()
    if opts.file_format == 'excel':
        df = pd.read_excel(opts.file_in, skiprows=opts.skip_rows)
    else:
        df = pd.read_csv(opts.file_in, sep='\t', low_memory=False, skiprows=opts.skip_rows, dtype={'#Chr': str})
    df.rename(columns={'#Chr': 'Chr'}, inplace=True)
    while True:
        while not records.full() and not input_finished:
            try:
                record = next(df.itertuples())
                chrome = getattr(record, 'Chr')
                start = getattr(record, 'Start')
                stop = getattr(record, 'Stop')
                ref = getattr(record, 'Ref')
                call = getattr(record, 'Call')
                var_type = judge_var_type(ref, call)
                record_parser = VariantRecord(chrome, start, stop, ref, call, var_type)
                g = generate_g(record_parser, chrome_dic)
                records.put((record_parser, g))
            except StopIteration:
                input_finished = True
                records.put('END')
                break
        processes_status = list()
        for p in processes:
            processes_status.append(p.is_alive())
        if True not in processes_status:
            results.put('END')
        while True:
            try:
                result = results.get(False)
            except queue.Empty:
                break
            if result != 'END':
                write_annotate_result(fp, result)
            else:
                output_finished = True
                break
        if output_finished:
            break
    fp.close()


def main():
    parser = optparse.OptionParser()
    parser.add_option('-i', dest='file_in', help='input file', default=None, metavar='file')
    parser.add_option('-o', dest='file_out', help='output file', default=None, metavar='file')
    parser.add_option('--it', dest='input_type', help='input file type, vcf or bed', default='vcf', metavar='string')
    parser.add_option('--ot', dest='out_type', help='output file type, vcf or bed', default='vcf', metavar='string')
    parser.add_option('-c', dest='config', help='config file', default=None, metavar='file')
    parser.add_option('-a', dest='annotation', help='annotation', default='GRCh37', metavar='string')
    parser.add_option('-p', dest='processes', help='process num', default=1, type=int, metavar='int')
    parser.add_option('--how', dest='how', help='how to select trans', default='all', metavar='string')
    parser.add_option('--rule_3', dest='rule_3', help='3\' rule', default='transcript', metavar='string')
    parser.add_option('--format', dest='file_format', help='file format', default='excel', metavar='string')
    parser.add_option('-k', dest='skip_rows', help='skip rows', default=0, type=int, metavar='int')
    (opts, args) = parser.parse_args()
    if opts.config is None:
        trans_provided = list()
        trans_provided_no_acc = list()
    else:
        config_dic = yaml_read(opts.config)
        trans_df = pd.read_csv(config_dic['Transcript'], sep='\t', header=None)
        trans_provided = trans_df[0].values.tolist()
        trans_provided_no_acc = remove_trans_acc(trans_provided)
    process_num = min(opts.processes, cpu_count())
    if opts.input_type == 'vcf':
        if opts.out_type == 'vcf':
            serial_annotate(opts, trans_provided_no_acc) if process_num == 1 else parallel_annotate(opts, trans_provided_no_acc, process_num)
        else:
            serial_annotate_to_bed(opts, trans_provided_no_acc) if process_num == 1 else parallel_annotate_to_bed(opts, trans_provided_no_acc, process_num)
    else:
        bgianno_serial_annotate(opts, trans_provided) if process_num == 1 else bgianno_parallel_annotate(opts, trans_provided, process_num)


if __name__ == '__main__':
    main()
