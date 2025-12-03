# -*- coding:utf-8 -*-
from vcf.parser import _Info as VcfInfo, field_counts as vcf_field_counts
from multiprocessing import Queue, Process, Pool, cpu_count
from maxentpy.maxent import load_matrix5, load_matrix3
from maxentpy.maxent_fast import score5, score3
from functools import partial, reduce
from argparse import ArgumentParser
from collections import namedtuple
from itertools import count
import pandas as pd
import pyfaidx
import queue
import vcf
import os


ROOT = os.path.abspath(os.path.dirname(__file__))
ANNOTATION = os.path.join(ROOT, 'annotation')
VariantRecord = namedtuple('VariantRecord', ['chromosome', 'pos', 'ref', 'alt'])
VariantAnnotation = namedtuple('VariantAnnotation', ['location', 'exon_nearby', 'exon_start', 'exon_end', 'strand', 'gene'])


class MaxEntScan:
    def __init__(self, annotation, reference):
        self.annotation = pd.read_csv(os.path.join(ANNOTATION, annotation + '.exon.gff.tsv'), sep='\t', header=None, dtype={0: str})
        self.matrix5 = load_matrix5()
        self.matrix3 = load_matrix3()
        self.reference = pyfaidx.Fasta(reference)

    def mes_score(self, record):
        if not (record.ref in ['A', 'T', 'C', 'G', 'a', 't', 'c', 'g'] and record.alt in ['A', 'T', 'C', 'G', 'a', 't', 'c', 'g']):
            return record.alt + '|.|.|.|.'
        record_anno = MaxEntScan.variant_annotation(self.annotation, record, 0)
        seq_ss, seq_ss_mut, ss = '', '', ''
        if record_anno.location == 'exon':
            if record_anno.strand == '+':
                if record.pos + 2 >= record_anno.exon_end:
                    seq_ss, seq_ss_mut, ss = MaxEntScan.forward_exon_ss5(self.reference, record, record_anno.exon_end)
                elif record.pos - 2 <= record_anno.exon_start:
                    seq_ss, seq_ss_mut, ss = MaxEntScan.forward_exon_ss3(self.reference, record, record_anno.exon_start)
            else:
                if record.pos - 2 <= record_anno.exon_start:
                    seq_ss, seq_ss_mut, ss = MaxEntScan.reverse_exon_ss5(self.reference, record, record_anno.exon_start)
                elif record.pos + 2 >= record_anno.exon_end:
                    seq_ss, seq_ss_mut, ss = MaxEntScan.reverse_exon_ss3(self.reference, record, record_anno.exon_end)
            return '|'.join([record.alt, record_anno.gene, *MaxEntScan.score_ss_seq(self.matrix5, self.matrix3, seq_ss, seq_ss_mut, ss)])
        else:
            record_anno = MaxEntScan.variant_annotation(self.annotation, record, -6)
            if record_anno.exon_nearby:
                if record_anno.strand == '+':
                    seq_ss, seq_ss_mut, ss = MaxEntScan.forward_intron_ss5(self.reference, record, record_anno.exon_end)
                    return '|'.join([record.alt, record_anno.gene, *MaxEntScan.score_ss_seq(self.matrix5, self.matrix3, seq_ss, seq_ss_mut, ss)])
            record_anno = MaxEntScan.variant_annotation(self.annotation, record, 20)
            if record_anno.exon_nearby:
                if record_anno.strand == '+':
                    seq_ss, seq_ss_mut, ss = MaxEntScan.forward_intron_ss3(self.reference, record, record_anno.exon_start)
                    return '|'.join([record.alt, record_anno.gene, *MaxEntScan.score_ss_seq(self.matrix5, self.matrix3, seq_ss, seq_ss_mut, ss)])
            record_anno = MaxEntScan.variant_annotation(self.annotation, record, 6)
            if record_anno.exon_nearby:
                if record_anno.strand == '-':
                    seq_ss, seq_ss_mut, ss = MaxEntScan.reverse_intron_ss5(self.reference, record, record_anno.exon_start)
                    return '|'.join([record.alt, record_anno.gene, *MaxEntScan.score_ss_seq(self.matrix5, self.matrix3, seq_ss, seq_ss_mut, ss)])
            record_anno = MaxEntScan.variant_annotation(self.annotation, record, -20)
            if record_anno.exon_nearby:
                if record_anno.strand == '-':
                    seq_ss, seq_ss_mut, ss = MaxEntScan.reverse_intron_ss3(self.reference, record, record_anno.exon_end)
                    return '|'.join([record.alt, record_anno.gene, *MaxEntScan.score_ss_seq(self.matrix5, self.matrix3, seq_ss, seq_ss_mut, ss)])
            return '|'.join([record.alt, '.', *MaxEntScan.score_ss_seq(self.matrix5, self.matrix3, seq_ss, seq_ss_mut, ss)])

    @staticmethod
    def variant_annotation(annotation, record, offset):
        exon = annotation[(annotation[0] == record.chromosome) & (annotation[1] <= record.pos + offset) & (record.pos + offset <= annotation[2])]
        if exon.empty:
            location = 'intron'
            exon_nearby = False
            exon_start = None
            exon_end = None
            strand = None
            gene = None
        else:
            location = 'exon'
            exon_nearby = True
            exon_start = exon[1].values[0]
            exon_end = exon[2].values[0]
            strand = exon[3].values[0]
            gene = exon[4].values[0]
        return VariantAnnotation(location, exon_nearby, exon_start, exon_end, strand, gene)

    @staticmethod
    def forward_exon_ss5(reference, record, exon_end):
        seq_ss5 = str(reference.get_seq(record.chromosome, exon_end - 2, exon_end + 6))
        seq_ss5 = seq_ss5[0:3].lower() + seq_ss5[3:5] + seq_ss5[5:].lower()
        seq_ss5_mut = seq_ss5[0:record.pos - exon_end + 2] + record.alt + seq_ss5[record.pos - exon_end + 3:]
        return seq_ss5, seq_ss5_mut, 'ss5'

    @staticmethod
    def reverse_exon_ss5(reference, record, exon_start):
        seq_ss5 = str(reference.get_seq(record.chromosome, exon_start - 6, exon_start + 2).reverse.complement)
        seq_ss5 = seq_ss5[0:3].lower() + seq_ss5[3:5] + seq_ss5[5:].lower()
        seq_ss5_mut = seq_ss5[0:exon_start - record.pos + 2] + MaxEntScan.base_complement(record.alt) + seq_ss5[exon_start - record.pos + 3:]
        return seq_ss5, seq_ss5_mut, 'ss5'

    @staticmethod
    def forward_exon_ss3(reference, record, exon_start):
        seq_ss3 = str(reference.get_seq(record.chromosome, exon_start - 20, exon_start + 2))
        seq_ss3 = seq_ss3[0:18].lower() + seq_ss3[18:20] + seq_ss3[20:].lower()
        seq_ss3_mut = seq_ss3[0:20] + seq_ss3[20:record.pos - exon_start + 20] + record.alt + seq_ss3[record.pos - exon_start + 21:]
        return seq_ss3, seq_ss3_mut, 'ss3'

    @staticmethod
    def reverse_exon_ss3(reference, record, exon_end):
        seq_ss3 = str(reference.get_seq(record.chromosome, exon_end - 2, exon_end + 20).reverse.complement)
        seq_ss3 = seq_ss3[0:18].lower() + seq_ss3[18:20] + seq_ss3[20:].lower()
        seq_ss3_mut = seq_ss3[0:20] + seq_ss3[20:exon_end - record.pos + 20] + MaxEntScan.base_complement(record.alt) + seq_ss3[exon_end - record.pos + 21:]
        return seq_ss3, seq_ss3_mut, 'ss3'

    @staticmethod
    def forward_intron_ss5(reference, record, exon_end):
        seq_ss5 = str(reference.get_seq(record.chromosome, exon_end - 2, exon_end + 6))
        seq_ss5 = seq_ss5[0:3].lower() + seq_ss5[3:5] + seq_ss5[5:].lower()
        seq_ss5_mut = seq_ss5[0:3] + seq_ss5[3:record.pos - exon_end + 2] + record.alt + seq_ss5[record.pos - exon_end + 3:]
        return seq_ss5, seq_ss5_mut, 'ss5'

    @staticmethod
    def reverse_intron_ss5(reference, record, exon_start):
        seq_ss5 = str(reference.get_seq(record.chromosome, exon_start - 6, exon_start + 2).reverse.complement)
        seq_ss5 = seq_ss5[0:3].lower() + seq_ss5[3:5] + seq_ss5[5:].lower()
        seq_ss5_mut = seq_ss5[0:3] + seq_ss5[3:exon_start - record.pos + 2] + MaxEntScan.base_complement(record.alt) + seq_ss5[exon_start - record.pos + 3:]
        return seq_ss5, seq_ss5_mut, 'ss5'

    @staticmethod
    def forward_intron_ss3(reference, record, exon_start):
        seq_ss3 = str(reference.get_seq(record.chromosome, exon_start - 20, exon_start + 2))
        seq_ss3 = seq_ss3[0:18].lower() + seq_ss3[18:20] + seq_ss3[20:].lower()
        seq_ss3_mut = seq_ss3[0:record.pos - exon_start + 20] + record.alt + seq_ss3[record.pos - exon_start + 21:]
        return seq_ss3, seq_ss3_mut, 'ss3'

    @staticmethod
    def reverse_intron_ss3(reference, record, exon_end):
        seq_ss3 = str(reference.get_seq(record.chromosome, exon_end - 2, exon_end + 20).reverse.complement)
        seq_ss3 = seq_ss3[0:18].lower() + seq_ss3[18:20] + seq_ss3[20:].lower()
        seq_ss3_mut = seq_ss3[0:exon_end - record.pos + 20] + MaxEntScan.base_complement(record.alt) + seq_ss3[exon_end - record.pos + 21:]
        return seq_ss3, seq_ss3_mut, 'ss3'

    @staticmethod
    def score_ss_seq(matrix5, matrix3, seq_ss, seq_ss_mut, ss):
        if 'N' in seq_ss or 'N' in seq_ss_mut or 'n' in seq_ss or 'n' in seq_ss_mut:
            return '.', '.', '.'
        if ss == 'ss5':
            score_wt = score5(seq_ss, matrix=matrix5)
            score_mut = score5(seq_ss_mut, matrix=matrix5)
            return str(round(score_wt, 3)), str(round(score_mut, 3)), str(round((score_mut - score_wt) / score_wt, 3))
        elif ss == 'ss3':
            score_wt = score3(seq_ss, matrix=matrix3)
            score_mut = score3(seq_ss_mut, matrix=matrix3)
            return str(round(score_wt, 3)), str(round(score_mut, 3)), str(round((score_mut - score_wt) / score_wt, 3))
        else:
            return '.', '.', '.'

    @staticmethod
    def base_complement(dna_base):
        complement_rule = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'a': 'T', 't': 'A', 'c': 'G', 'g': 'C'}
        return complement_rule[dna_base]


def split_df(df, split_num):
    df.reset_index(drop=True, inplace=True)
    df_list = list()
    step = round(df.shape[0]/split_num)
    for i in range(split_num):
        if i == 0:
            df_list.append(df.loc[0: step-1])
        elif i == split_num-1:
            df_list.append(df.loc[step*i:])
        else:
            df_list.append(df.loc[step*i:step*(i+1)-1])
    return df_list


def score_pseudo_vcf(annotation, reference, format_in, df):
    mes = MaxEntScan(annotation, reference)
    if format_in == 'vcf-4cols':
        df['MaxEntScan'] = df.apply(lambda x: mes.mes_score(VariantRecord(x['#CHROM'], x['POS'], x['REF'], x['ALT'])), axis=1)
    else:
        df['MaxEntScan'] = df.apply(lambda x: mes.mes_score(VariantRecord(x['#Chr'], x['Stop'], x['Ref'], x['Call'])), axis=1)
    return df


def annotation_pseudo_vcf(parsed_args, process_num):
    if parsed_args.format_in == 'vcf-4cols':
        df = pd.read_csv(parsed_args.file_in, sep='\t', dtype={'#CHROM': str})
    else:
        df = pd.read_csv(parsed_args.file_in, sep='\t', dtype={'#Chr': str})
    df_list = split_df(df, process_num)
    score_partial = partial(score_pseudo_vcf, parsed_args.annotation, parsed_args.reference, parsed_args.format_in)
    with Pool(process_num) as pool:
        df_list_score = pool.map(score_partial, df_list)
    df_score = reduce(lambda x, y: x.append(y), df_list_score)
    df_score.to_csv(parsed_args.file_out, sep='\t', index=False)


def score_vcf(records, results, annotation, reference):
    mes = MaxEntScan(annotation, reference)
    while True:
        try:
            record = records.get(False)
        except queue.Empty:
            continue
        if record != 'END':
            record_id, record_infos, record_score_list = record[0], record[1], list()
            for record_info in record_infos:
                record_score_list.append(mes.mes_score(record_info))
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
        p = Process(target=score_vcf, args=(records, results, parsed_args.annotation, parsed_args.reference))
        processes.append(p)
        p.start()
    vcf_reader = vcf.Reader(filename=parsed_args.file_in)
    vcf_reader.infos['MaxEntScan'] = VcfInfo('MaxEntScan', vcf_field_counts['A'], 'String',
                                             'MaxEntScan Score for VCF record alleles, related_score = (mut_score-wild_score)/wild_score Format: ALLELE|SYMBOL|wild_score|mut_score|related_score',
                                             version=None, source=None)
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
                record_write.add_info('MaxEntScan', record_score)
                vcf_writer.write_record(record_write)
            else:
                output_finished = True
                break
        if output_finished:
            break
    vcf_writer.close()


def main():
    parser = ArgumentParser()
    parser.add_argument('-r', help='Reference genome fasta file', required=True, dest='reference')
    parser.add_argument('-a', help='Genome annotation: GRCh37 or GRCh38', default='GRCh37', dest='annotation')
    parser.add_argument('-i', help='Input file', required=True, dest='file_in')
    parser.add_argument('-o', help='Output file', required=True, dest='file_out')
    parser.add_argument('--format_in', help='Input file format: vcf,vcf-4cols,bgianno', default='vcf-4cols', dest='format_in')
    parser.add_argument('-p', help='Process number', type=int, default=1, dest='processes')
    parsed_args = parser.parse_args()
    if parsed_args.annotation not in ['GRCh37', 'GRCh38']:
        raise Exception('Genome not recognized! Must be GRCh37 or GRCh38')
    if parsed_args.format_in not in ['vcf', 'vcf-4cols', 'bgianno']:
        raise Exception('Input file format not recognized! Must be vcf, vcf-4cols or bgianno')
    process_num = min(cpu_count(), parsed_args.processes)
    if parsed_args.format_in != 'vcf':
        annotation_pseudo_vcf(parsed_args, process_num)
    else:
        annotation_vcf(parsed_args, process_num)


if __name__ == '__main__':
    main()
