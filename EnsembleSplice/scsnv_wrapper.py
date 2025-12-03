# -*- coding:utf-8 -*-
from vcf.parser import _Info as VcfInfo, field_counts as vcf_field_counts
from multiprocessing import Queue, Process, cpu_count
from argparse import ArgumentParser
from collections import namedtuple
from itertools import count
import pandas as pd
import queue
import tabix
import vcf
import os


ROOT = os.path.abspath(os.path.dirname(__file__))
ANNOTATION = os.path.join(ROOT, 'annotation')
VariantRecord = namedtuple('VariantRecord', ['chromosome', 'pos', 'ref', 'alt'])


class ScSNV:
    def __init__(self, annotation):
        self.annotation = tabix.open(os.path.join(ANNOTATION, 'dbscSNV1.1_' + annotation + '.tsv.gz'))

    def scsnv_score(self, record):
        score = record.alt + '|.|.'
        if not (record.ref in ['A', 'T', 'C', 'G', 'a', 't', 'c', 'g'] and record.alt in ['A', 'T', 'C', 'G', 'a', 't', 'c', 'g']):
            return score
        records_query = self.annotation.query(record.chromosome, record.pos - 1, record.pos)
        if records_query:
            for record_query in records_query:
                if [record.chromosome, str(record.pos), record.ref, record.alt] == record_query[0: 4]:
                    score = '|'.join([record.alt, record_query[4], record_query[5]])
                    break
        return score


def annotation_pseudo_vcf(parsed_args):
    if parsed_args.format_in == 'vcf-4cols':
        df = pd.read_csv(parsed_args.file_in, sep='\t', dtype={'#CHROM': str})
        df['#chr'] = df['#CHROM']
        df['pos'] = df['POS']
        df['ref'] = df['REF']
        df['alt'] = df['ALT']
    else:
        df = pd.read_csv(parsed_args.file_in, sep='\t', dtype={'#Chr': str})
        df['#chr'] = df['#Chr']
        df['pos'] = df['Stop']
        df['ref'] = df['Ref']
        df['alt'] = df['Call']
    db = pd.read_csv(os.path.join(ANNOTATION, 'dbscSNV1.1_' + parsed_args.annotation + '.tsv.gz'),
                     sep='\t', dtype={'#chr': str, 'ada_score': str, 'rf_score': str})
    df_score = pd.merge(df, db, on=['#chr', 'pos', 'ref', 'alt'], how='left')
    df_score.fillna('.', inplace=True)
    df_score['dbscSNV'] = df_score['ada_score'] + '|' + df_score['rf_score']
    df_score.drop(columns=['#chr', 'pos', 'ref', 'alt', 'ada_score', 'rf_score'], inplace=True)
    df_score.to_csv(parsed_args.file_out, sep='\t', index=False)


def score_vcf(records, results, annotation):
    scsnv = ScSNV(annotation)
    while True:
        try:
            record = records.get(False)
        except queue.Empty:
            continue
        if record != 'END':
            record_id, record_infos, record_score_list = record[0], record[1], list()
            for record_info in record_infos:
                record_score_list.append(scsnv.scsnv_score(record_info))
            results.put((record_id, ','.join(record_score_list)))
        else:
            records.put('END')
            break


def annotation_vcf(parsed_args, process_num):
    records, results = Queue(100 * process_num), Queue()
    input_finished = False
    output_finished = False
    wait_records = dict()
    processes = list()
    records_id = count()
    for i in range(process_num):
        p = Process(target=score_vcf, args=(records, results, parsed_args.annotation))
        processes.append(p)
        p.start()
    vcf_reader = vcf.Reader(filename=parsed_args.file_in)
    vcf_reader.infos['dbscSNV'] = VcfInfo('dbscSNV', vcf_field_counts['A'], 'String',
                                          'dbscSNV Score for VCF record alleles, Format: ALLELE|ada_score|rf_score', version=None, source=None)
    vcf_writer = vcf.Writer(open(parsed_args.file_out, 'w'), vcf_reader)
    while True:
        while not records.full() and not input_finished:
            try:
                record = next(vcf_reader)
                record_id = next(records_id)
                wait_records[record_id] = record
                record_infos = list()
                chromosome = str(record.CHROM)
                pos = record.POS
                ref = record.REF
                for alt in record.ALT:
                    record_infos.append(VariantRecord(chromosome, pos, ref, str(alt)))
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
                record_id, record_score = result[0], result[1]
                record_write = wait_records.pop(record_id)
                record_write.add_info('dbscSNV', record_score)
                vcf_writer.write_record(record_write)
            else:
                output_finished = True
                break
        if output_finished:
            break
    vcf_writer.close()


def main():
    parser = ArgumentParser()
    parser.add_argument('-a', help='Genome hg19 or hg38', default='hg19', dest='annotation')
    parser.add_argument('-i', help='Input file', required=True, dest='file_in')
    parser.add_argument('-o', help='Output file', required=True, dest='file_out')
    parser.add_argument('--format_in', help='Input file format: vcf,vcf-4cols,bgianno', default='vcf-4cols', dest='format_in')
    parser.add_argument('-p', help='Process number', type=int, default=1, dest='processes')
    parsed_args = parser.parse_args()
    if parsed_args.annotation not in ['hg19', 'hg38']:
        raise Exception('Genome not recognized! Must be hg19 or hg38')
    if parsed_args.format_in not in ['vcf', 'vcf-4cols', 'bgianno']:
        raise Exception('Input file format not recognized! Must be vcf, vcf-4cols or bgianno')
    process_num = min(cpu_count(), parsed_args.processes)
    if parsed_args.format_in != 'vcf':
        annotation_pseudo_vcf(parsed_args)
    else:
        annotation_vcf(parsed_args, process_num)


if __name__ == '__main__':
    main()
