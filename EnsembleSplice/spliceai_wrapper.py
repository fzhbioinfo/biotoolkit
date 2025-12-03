# -*- coding:utf-8 -*-
from vcf.parser import _Info as VcfInfo, field_counts as vcf_field_counts
from multiprocessing import Queue, Process, cpu_count
from argparse import ArgumentParser
from collections import namedtuple
from itertools import count
import queue
import tabix
import vcf
import os


ROOT = os.path.abspath(os.path.dirname(__file__))
ANNOTATION = os.path.join(ROOT, 'annotation')
VariantRecord = namedtuple('VariantRecord', ['chromosome', 'pos', 'ref', 'alt'])


class SpliceAI:
    def __init__(self, annotation):
        self.annotation_snv = tabix.open(os.path.join(ANNOTATION, 'spliceai_scores.raw.snv.' + annotation + '.vcf.gz'))
        self.annotation_indel = tabix.open(os.path.join(ANNOTATION, 'spliceai_scores.raw.indel.' + annotation + '.vcf.gz'))

    def spliceai_score(self, record):
        score = record.alt + '|.|.|.|.|.|.|.|.|.'
        if SpliceAI.mutation_type(record.ref, record.alt) == 'snv':
            records_query = self.annotation_snv.query(record.chromosome, record.pos - 1, record.pos)
        else:
            records_query = self.annotation_indel.query(record.chromosome, record.pos - 1, record.pos)
        if records_query:
            score_list = list()
            for record_query in records_query:
                if [record.chromosome, str(record.pos), record.ref, record.alt] == [record_query[0], record_query[1], record_query[3], record_query[4]]:
                    score_list.append(record_query[7].replace('SpliceAI=', ''))
            if score_list:
                score = '/'.join(score_list)
            else:
                score = record.alt + '|.|.|.|.|.|.|.|.|.'
        return score

    @staticmethod
    def mutation_type(ref, alt):
        if ref in ['A', 'T', 'C', 'G', 'a', 't', 'c', 'g'] and alt in ['A', 'T', 'C', 'G', 'a', 't', 'c', 'g']:
            return 'snv'
        else:
            return 'indel'


def score_vcf(records, results, annotation):
    spliceai = SpliceAI(annotation)
    while True:
        try:
            record = records.get(False)
        except queue.Empty:
            continue
        if record != 'END':
            record_id, record_infos, record_score_list = record[0], record[1], list()
            for record_info in record_infos:
                record_score_list.append(spliceai.spliceai_score(record_info))
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
    vcf_reader.infos['SpliceAI'] = VcfInfo('SpliceAI', vcf_field_counts['A'], 'String',
                                           'SpliceAIv1.3 variant annotation. These include delta scores (DS) and delta positions (DP) for '
                                           'acceptor gain (AG), acceptor loss (AL), donor gain (DG), and donor loss (DL). '
                                           'Format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL', version=None, source=None)
    vcf_writer = vcf.Writer(open(parsed_args.file_out, 'w'), vcf_reader)
    while True:
        while not records.full() and not input_finished:
            try:
                record = next(vcf_reader)
                record_id = next(records_id)
                wait_records[record_id] = record
                record_infos = list()
                chromosome = str(record.CHROM).replace('chr', '')
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
                record_write.add_info('SpliceAI', record_score)
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
    parser.add_argument('--format_in', help='Input file format: vcf,vcf-4cols,bgianno', default='vcf', dest='format_in')
    parser.add_argument('-p', help='Process number', type=int, default=1, dest='processes')
    parsed_args = parser.parse_args()
    if parsed_args.annotation not in ['hg19', 'hg38']:
        raise Exception('Genome not recognized! Must be hg19 or hg38')
    if parsed_args.format_in not in ['vcf', 'vcf-4cols', 'bgianno']:
        raise Exception('Input file format not recognized! Must be vcf, vcf-4cols or bgianno')
    process_num = min(cpu_count(), parsed_args.processes)
    if parsed_args.format_in != 'vcf':
        raise Exception('Only support vcf Input file format!')
    else:
        annotation_vcf(parsed_args, process_num)


if __name__ == '__main__':
    main()
