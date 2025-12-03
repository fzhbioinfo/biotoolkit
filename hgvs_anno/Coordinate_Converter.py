# -*- coding:utf-8 -*-
from hgvs.exceptions import HGVSError, HGVSParseError, HGVSUsageError
from multiprocessing import Queue, Process, cpu_count
from bioutils.assemblies import make_name_ac_map
from hgvs.assemblymapper import AssemblyMapper
from hgvs.dataproviders.uta import connect
from hgvs.normalizer import Normalizer
from collections import namedtuple
from operator import itemgetter
from hgvs.parser import Parser
import pandas as pd
import optparse
import logging
import queue
import re


logger = logging.getLogger()
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
handler = logging.StreamHandler()
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.WARNING)


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
    # UTA中相关的转录本与所提供的选择最新的版本, 若未提供或无交集返回空列表
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
        return []
    trans_filter = [tran for tran_overlap in trans_no_acc_overlap for tran in trans_related if tran.startswith(tran_overlap + '.')]
    trans_filter = get_sorted_trans_acc(trans_filter)
    trans_selected = list()
    for tran_overlap in trans_no_acc_overlap:
        for tran in trans_filter:
            if tran.startswith(tran_overlap + '.'):
                trans_selected.append(tran)
                break
    return trans_selected


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
    hn = Normalizer(hdp)
    return am, hp, hn, hdp


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
    fp.write('hgvs-cHGVS' + '\t')
    fp.write('hgvs-cHGVS_Normalise' + '\t')
    fp.write('hgvs-pHGVS' + '\t')
    fp.write('hgvs-pHGVS_Normalise' + '\n')
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


def serial_annotate_to_bed(opts, trans_provided_no_acc):
    # 串行注释, 生成bed格式
    fp = create_annotate_result(opts.file_out)
    am, hp, hn, hdp = annotator(opts.annotation)
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
            g_normalise = hn.normalize(g_parser)
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
                t_normalise = am.g_to_t(g_normalise, tran)
                p = am.t_to_p(t)
                p_normalise = am.t_to_p(t_normalise)
                record_annotation = RecordAnnotation(chrome, start, stop, ref, call, var_type, str(t), str(t_normalise), str(p), str(p_normalise))
                write_annotate_result(fp, record_annotation)
            except (HGVSError, HGVSUsageError, NotImplementedError, IndexError) as e:
                error = str(e)
                logging.error('{chrome} {start} {stop} {ref} {call} {tran} {g} annotate error. {error}.'.format(**locals()))
    fp.close()


def process_record_to_bed(records, results, annotation, trans_provided, how):
    # 一个进程创建uta连接， 从records队列里获取数据, 注释结果返回到results队列, 适用输出为bed格式
    am, hp, hn, hdp = annotator(annotation)
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
                g_normalise = hn.normalize(g_parser)
                trans_related = am.relevant_transcripts(g_parser)
            except (HGVSParseError, HGVSError, HGVSUsageError) as e:
                error = str(e)
                logging.error('{chrome} {start} {stop} {ref} {call} {g} annotate error. {error}.'.format(**locals()))
                continue
            trans = select_trans(trans_related, trans_provided, how)
            if len(trans) == 0:
                logging.warning('{chrome} {start} {stop} {ref} {call} {g} no related transcripts in UTA.'.format(**locals()))
                continue
            for tran in trans:
                try:
                    t = am.g_to_t(g_parser, tran)
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


def parallel_annotate_to_bed(opts, trans_provided_no_acc, process_num):
    # 并行注释, 输出bed格式结果
    chrome_dic = generate_chrome_dic(opts.annotation)
    fp = create_annotate_result(opts.file_out)
    records, results = Queue(100*process_num), Queue()
    input_finished = False
    output_finished = False
    processes = list()
    for i in range(process_num):
        p = Process(target=process_record_to_bed, args=(records, results, opts.annotation, trans_provided_no_acc, opts.how))
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


def t_annotate(opts):
    am, hp, hn, hdp = annotator(opts.annotation)
    if opts.file_format == 'excel':
        df = pd.read_excel(opts.file_in, header=None)
    else:
        df = pd.read_csv(opts.file_in, sep='\t', header=None)
    df.columns = ['Input Variant']
    df['Input Variant'] = df['Input Variant'].str.replace(' ', '', regex=True)
    df['Normalise'] = '.'
    df['g.'] = '.'
    df['#Chr'] = '.'
    df['Start'] = '.'
    df['Stop'] = '.'
    df['Ref'] = '.'
    df['Call'] = '.'
    df['MuType'] = '.'
    df['State'] = '.'
    for i in range(df.shape[0]):
        try:
            var_c = hp.parse_hgvs_variant(df.loc[i, 'Input Variant'])
            df.loc[i, 'Normalise'] = str(hn.normalize(var_c))
            var_g = am.t_to_g(var_c)
        except (HGVSError, HGVSParseError, HGVSUsageError, NotImplementedError, IndexError) as e:
            df.loc[i, 'State'] = str(e)
            continue
        if var_g.posedit.edit.type == 'dup':
            df.loc[i, '#Chr'] = var_g.ac
            df.loc[i, 'Start'] = var_g.posedit.pos.end.base
            df.loc[i, 'Stop'] = var_g.posedit.pos.end.base
            df.loc[i, 'Call'] = var_g.posedit.edit.ref
            df.loc[i, 'MuType'] = var_g.posedit.edit.type
            df.loc[i, 'State'] = 'success'
            df.loc[i, 'g.'] = str(var_g)
        elif var_g.posedit.edit.type == 'ins':
            df.loc[i, '#Chr'] = var_g.ac
            df.loc[i, 'Start'] = var_g.posedit.pos.start.base
            df.loc[i, 'Stop'] = var_g.posedit.pos.start.base
            df.loc[i, 'Call'] = var_g.posedit.edit.alt
            df.loc[i, 'MuType'] = var_g.posedit.edit.type
            df.loc[i, 'State'] = 'success'
            df.loc[i, 'g.'] = str(var_g)
        elif var_g.posedit.edit.type == 'sub':
            df.loc[i, '#Chr'] = var_g.ac
            df.loc[i, 'Start'] = var_g.posedit.pos.start.base - 1
            df.loc[i, 'Stop'] = var_g.posedit.pos.end.base
            df.loc[i, 'Ref'] = var_g.posedit.edit.ref
            df.loc[i, 'Call'] = var_g.posedit.edit.alt
            df.loc[i, 'MuType'] = var_g.posedit.edit.type
            df.loc[i, 'State'] = 'success'
            df.loc[i, 'g.'] = str(var_g)
        elif var_g.posedit.edit.type.startswith('del'):
            df.loc[i, '#Chr'] = var_g.ac
            df.loc[i, 'Start'] = var_g.posedit.pos.start.base - 1
            df.loc[i, 'Stop'] = var_g.posedit.pos.end.base
            df.loc[i, 'Ref'] = var_g.posedit.edit.ref
            df.loc[i, 'Call'] = var_g.posedit.edit.alt
            df.loc[i, 'MuType'] = var_g.posedit.edit.type
            df.loc[i, 'State'] = 'success'
            df.loc[i, 'g.'] = str(var_g)
        else:
            df.loc[i, '#Chr'] = var_g.ac
            df.loc[i, 'Start'] = var_g.posedit.pos.start.base - 1
            df.loc[i, 'Stop'] = var_g.posedit.pos.end.base
            df.loc[i, 'Ref'] = var_g.posedit.edit.ref
            df.loc[i, 'Call'] = var_g.posedit.edit.alt
            df.loc[i, 'MuType'] = var_g.posedit.edit.type
            df.loc[i, 'State'] = 'success'
            df.loc[i, 'g.'] = str(var_g)
    df.fillna(value={'Call': '.'}, inplace=True)
    df['#Chr'] = df['#Chr'].str.extract(r'NC_0+(\d+)\.\d+')
    df.fillna(value={'#Chr': '.'}, inplace=True)
    df.replace({'#Chr': {'23': 'X', '24': 'Y'}}, inplace=True)
    df.to_csv(opts.file_out, sep='\t', index=False)


def g_annotate(opts):
    am, hp, hn, hdp = annotator(opts.annotation)
    if opts.file_format == 'excel':
        df = pd.read_excel(opts.file_in, header=None)
    else:
        df = pd.read_csv(opts.file_in, sep='\t', header=None)
    df.columns = ['Input Variant']
    df['Input Variant'] = df['Input Variant'].str.replace(' ', '', regex=True)
    df['#Chr'] = '.'
    df['Start'] = '.'
    df['Stop'] = '.'
    df['Ref'] = '.'
    df['Call'] = '.'
    df['MuType'] = '.'
    df['State'] = '.'
    for i in range(df.shape[0]):
        try:
            var_g = hp.parse_hgvs_variant(df.loc[i, 'Input Variant'])
        except (HGVSError, HGVSParseError, HGVSUsageError, NotImplementedError, IndexError) as e:
            df.loc[i, 'State'] = str(e)
            continue
        if var_g.posedit.edit.type == 'dup':
            df.loc[i, '#Chr'] = var_g.ac
            df.loc[i, 'Start'] = var_g.posedit.pos.end.base
            df.loc[i, 'Stop'] = var_g.posedit.pos.end.base
            df.loc[i, 'Call'] = var_g.posedit.edit.ref
            df.loc[i, 'MuType'] = var_g.posedit.edit.type
            df.loc[i, 'State'] = 'success'
        elif var_g.posedit.edit.type == 'ins':
            df.loc[i, '#Chr'] = var_g.ac
            df.loc[i, 'Start'] = var_g.posedit.pos.start.base
            df.loc[i, 'Stop'] = var_g.posedit.pos.start.base
            df.loc[i, 'Call'] = var_g.posedit.edit.alt
            df.loc[i, 'MuType'] = var_g.posedit.edit.type
            df.loc[i, 'State'] = 'success'
        elif var_g.posedit.edit.type == 'sub':
            df.loc[i, '#Chr'] = var_g.ac
            df.loc[i, 'Start'] = var_g.posedit.pos.start.base - 1
            df.loc[i, 'Stop'] = var_g.posedit.pos.end.base
            df.loc[i, 'Ref'] = var_g.posedit.edit.ref
            df.loc[i, 'Call'] = var_g.posedit.edit.alt
            df.loc[i, 'MuType'] = var_g.posedit.edit.type
            df.loc[i, 'State'] = 'success'
        elif var_g.posedit.edit.type.startswith('del'):
            df.loc[i, '#Chr'] = var_g.ac
            df.loc[i, 'Start'] = var_g.posedit.pos.start.base - 1
            df.loc[i, 'Stop'] = var_g.posedit.pos.end.base
            df.loc[i, 'Ref'] = var_g.posedit.edit.ref
            df.loc[i, 'Call'] = var_g.posedit.edit.alt
            df.loc[i, 'MuType'] = var_g.posedit.edit.type
            df.loc[i, 'State'] = 'success'
        else:
            df.loc[i, '#Chr'] = var_g.ac
            df.loc[i, 'Start'] = var_g.posedit.pos.start.base - 1
            df.loc[i, 'Stop'] = var_g.posedit.pos.end.base
            df.loc[i, 'Ref'] = var_g.posedit.edit.ref
            df.loc[i, 'Call'] = var_g.posedit.edit.alt
            df.loc[i, 'MuType'] = var_g.posedit.edit.type
            df.loc[i, 'State'] = 'success'
    df.fillna(value={'Call': '.'}, inplace=True)
    df['#Chr'] = df['#Chr'].str.extract(r'NC_0+(\d+)\.\d+')
    df.fillna(value={'#Chr': '.'}, inplace=True)
    df.replace({'#Chr': {'23': 'X', '24': 'Y'}}, inplace=True)
    df.to_csv(opts.file_out, sep='\t', index=False)


def gene_c_annotate(opts):
    am, hp, hn, hdp = annotator(opts.annotation)
    if opts.file_format == 'excel':
        df = pd.read_excel(opts.file_in, header=None)
    else:
        df = pd.read_csv(opts.file_in, sep='\t', header=None)
    df.columns = ['gene', 'c.']
    df.replace(' ', '', regex=True, inplace=True)
    annotate_res = list()
    for i in range(df.shape[0]):
        gene = df.loc[i, 'gene']
        c = df.loc[i, 'c.']
        tx_gene = pd.DataFrame.from_records(hdp.get_tx_for_gene(gene))
        if not tx_gene.empty:
            tx_gene_trans = tx_gene[(tx_gene[5] == 'splign') & (tx_gene[4].str.contains('NC_'))].drop_duplicates(subset=[3])[3].tolist()
            if len(tx_gene_trans) == 0:
                annotate_list = [gene, c, '.', '.', '.', '.', '.', '.', '.', '.', 'transcript not find']
                annotate_res.append(annotate_list)
                continue
            for tran in tx_gene_trans:
                try:
                    var_c = hp.parse_hgvs_variant(tran+':'+c)
                    var_g = am.t_to_g(var_c)
                except (HGVSError, HGVSParseError, HGVSUsageError, NotImplementedError, IndexError) as e:
                    annotate_list = [gene, c, tran, '.', '.', '.', '.', '.', '.', '.', str(e)]
                    annotate_res.append(annotate_list)
                    continue
                if var_g.posedit.edit.type == 'dup':
                    chrome = var_g.ac
                    start = var_g.posedit.pos.end.base
                    stop = var_g.posedit.pos.end.base
                    ref = '.'
                    call = var_g.posedit.edit.ref
                    mutype = var_g.posedit.edit.type
                elif var_g.posedit.edit.type == 'ins':
                    chrome = var_g.ac
                    start = var_g.posedit.pos.start.base
                    stop = var_g.posedit.pos.start.base
                    ref = '.'
                    call = var_g.posedit.edit.alt
                    mutype = var_g.posedit.edit.type
                elif var_g.posedit.edit.type == 'sub':
                    chrome = var_g.ac
                    start = var_g.posedit.pos.start.base - 1
                    stop = var_g.posedit.pos.end.base
                    ref = var_g.posedit.edit.ref
                    call = var_g.posedit.edit.alt
                    mutype = var_g.posedit.edit.type
                elif var_g.posedit.edit.type.startswith('del'):
                    chrome = var_g.ac
                    start = var_g.posedit.pos.start.base - 1
                    stop = var_g.posedit.pos.end.base
                    ref = var_g.posedit.edit.ref
                    call = var_g.posedit.edit.alt
                    mutype = var_g.posedit.edit.type
                else:
                    chrome = var_g.ac
                    start = var_g.posedit.pos.start.base - 1
                    stop = var_g.posedit.pos.end.base
                    ref = var_g.posedit.edit.ref
                    call = var_g.posedit.edit.alt
                    mutype = var_g.posedit.edit.type
                annotate_list = [gene, c, tran, str(var_g), chrome, start, stop, ref, call, mutype, 'success']
                annotate_res.append(annotate_list)
        else:
            annotate_list = [gene, c, '.', '.', '.', '.', '.', '.', '.', '.', 'transcript not find! Please check your gene name!']
            annotate_res.append(annotate_list)
    annotate_cols = ['gene', 'c.', 'Transcript', 'g.', '#Chr', 'Start', 'Stop', 'Ref', 'Call', 'MuType', 'State']
    annotate_df = pd.DataFrame(annotate_res)
    annotate_df.columns = annotate_cols
    annotate_df['#Chr'] = annotate_df['#Chr'].str.extract(r'NC_0+(\d+)\.\d+')
    annotate_df.fillna(value={'#Chr': '.'}, inplace=True)
    annotate_df.replace({'#Chr': {'23': 'X', '24': 'Y'}}, inplace=True)
    annotate_df.to_csv(opts.file_out, sep='\t', index=False)


def main():
    parser = optparse.OptionParser()
    parser.add_option('-i', dest='file_in', help='input file', default=None, metavar='file')
    parser.add_option('-o', dest='file_out', help='output file', default=None, metavar='file')
    parser.add_option('-t', dest='input_type', help='input file type', default=None, metavar='string')
    parser.add_option('--NM', dest='nm_list', help='Transcript list file', default=None, metavar='file')
    parser.add_option('-a', dest='annotation', help='annotation, default=GRCh37', default='GRCh37', metavar='string')
    parser.add_option('-p', dest='processes', help='process num, default=1', default=1, type=int, metavar='int')
    parser.add_option('-k', dest='skip_rows', help='skip rows num, default=0', default=0, type=int, metavar='int')
    parser.add_option('--how', dest='how', help='how to select trans, default=all', default='all', metavar='string')
    parser.add_option('--format', dest='file_format', help='input file format, tsv or excel', default='tsv', metavar='string')
    (opts, args) = parser.parse_args()
    if opts.nm_list is None:
        trans_provided_no_acc = list()
    else:
        trans_df = pd.read_csv(opts.nm_list, sep='\t', header=None)
        trans_provided = trans_df[0].values.tolist()
        trans_provided_no_acc = remove_trans_acc(trans_provided)
    if opts.input_type == 'c':
        t_annotate(opts)
    elif opts.input_type == 'g':
        g_annotate(opts)
    elif opts.input_type == 'gene_c':
        gene_c_annotate(opts)
    else:
        process_num = min(opts.processes, cpu_count())
        serial_annotate_to_bed(opts, trans_provided_no_acc) if process_num == 1 else parallel_annotate_to_bed(opts, trans_provided_no_acc, process_num)


if __name__ == '__main__':
    main()
