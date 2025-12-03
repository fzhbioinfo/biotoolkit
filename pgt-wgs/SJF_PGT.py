# -*- coding:utf-8 -*-
import os
import re
import sys
import time
import glob
import yaml
import logging
import subprocess
import pandas as pd
from argparse import ArgumentParser


logger = logging.getLogger()
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
logger.setLevel(logging.INFO)

# 配置文件目录和每步分析流程目录
ROOT = os.path.abspath(os.path.dirname(__file__))
CONFIG = os.path.join(ROOT, 'etc')
PIPELINE = os.path.join(ROOT, 'script')

# 染色体和突变类型,后续相应步骤进行拆分
chromosomes = [f'chr{str(i)}' for i in range(1, 23)] + ['chrX', 'chrY', 'chrM_NC_012920.1']
mut_type = ['SNP', 'INDEL']


class CreateJob:
    def __init__(self, flow, work_dir, info):
        # 流程步骤配置
        self.flow = flow
        self.flow.index = self.flow['Name'].tolist()
        # 分析的工作目录
        self.work_dir = work_dir
        # 样本信息
        self.info = info
        self.info['sampleID'] = self.info['sample'] + '_' + self.info['sampleID']

    def make(self, gene):
        # 创建每个分析步骤的脚本,根据配置中的SplitBy进行拆分
        # 每组家系有特定的目录结构
        # 特定的步骤有特定的输入参数: fq1,fq2,chromosome,mut_type,gene
        job_conf_dic = self.flow.to_dict('index')
        for step in job_conf_dic.keys():
            pipe = os.path.join(PIPELINE, f'{step}.sh')
            for sampleID in self.info['sampleID'].unique():
                sample_info = self.info.loc[self.info['sampleID'] == sampleID]
                sample = sample_info['sample'].values[0]
                if sample in ['Mother', 'Father'] or sample.startswith('Embryo'):
                    work_dir = self.work_dir
                    script_dir = os.path.join(self.work_dir, 'shell')
                else:
                    work_dir = os.path.join(self.work_dir, sampleID)
                    script_dir = os.path.join(work_dir, 'shell')
                if not os.path.exists(script_dir):
                    os.makedirs(script_dir)
                if job_conf_dic[step]['Type'] == 'single':
                    if job_conf_dic[step]['SplitBy'] == 'barcode':
                        for index in sample_info.index:
                            fq1 = sample_info.loc[index, 'fq1']
                            fq2 = fq1.replace('_1.fq.gz', '_2.fq.gz')
                            fq_name = os.path.basename(fq1).replace('_1.fq.gz', '')
                            with open(os.path.join(script_dir, f'step{job_conf_dic[step]["Order"]}-{step}-{sample}-{fq_name}.sh'), 'w') as fp:
                                content = '''#!/bin/bash
sh {} {} {} {} {} {}'''.format(pipe, ROOT, work_dir, sample, fq1, fq2)
                                fp.write(content)
                    elif job_conf_dic[step]['SplitBy'] == 'chromosome':
                        for chromosome in chromosomes:
                            with open(os.path.join(script_dir, f'step{job_conf_dic[step]["Order"]}-{step}-{sample}-{chromosome}.sh'), 'w') as fp:
                                content = '''#!/bin/bash
sh {} {} {} {} {}'''.format(pipe, ROOT, work_dir, sample, chromosome)
                                fp.write(content)
                    elif job_conf_dic[step]['SplitBy'] == 'mut_type':
                        for mt in mut_type:
                            with open(os.path.join(script_dir, f'step{job_conf_dic[step]["Order"]}-{step}-{sample}-{mt}.sh'), 'w') as fp:
                                content = '''#!/bin/bash
sh {} {} {} {} {}'''.format(pipe, ROOT, work_dir, sample, mt)
                                fp.write(content)
                    else:
                        with open(os.path.join(script_dir, f'step{job_conf_dic[step]["Order"]}-{step}-{sample}.sh'), 'w') as fp:
                            if step in ['Stat']:
                                content = '''#!/bin/bash
sh {} {} {} {} {}'''.format(pipe, ROOT, work_dir, sample, gene)
                            else:
                                content = '''#!/bin/bash
sh {} {} {} {}'''.format(pipe, ROOT, work_dir, sample)
                            fp.write(content)
                else:
                    if sample in ['Mother', 'Father'] or sample.startswith('Embryo'):
                        continue
                    else:
                        if job_conf_dic[step]['SplitBy'] == 'chromosome':
                            for chromosome in chromosomes:
                                with open(os.path.join(script_dir, f'step{job_conf_dic[step]["Order"]}-{step}-{chromosome}.sh'), 'w') as fp:
                                    content = '''#!/bin/bash
sh {} {} {} {}'''.format(pipe, ROOT, work_dir, chromosome)
                                    fp.write(content)
                        elif job_conf_dic[step]['SplitBy'] == 'mut_type':
                            for mt in mut_type:
                                with open(os.path.join(script_dir, f'step{job_conf_dic[step]["Order"]}-{step}-{mt}.sh'), 'w') as fp:
                                    content = '''#!/bin/bash
sh {} {} {} {}'''.format(pipe, ROOT, work_dir, mt)
                                    fp.write(content)
                        else:
                            with open(os.path.join(script_dir, f'step{job_conf_dic[step]["Order"]}-{step}.sh'), 'w') as fp:
                                if step in ['Anno', 'Report']:
                                    content = '''#!/bin/bash
sh {} {} {} {}'''.format(pipe, ROOT, work_dir, gene)
                                else:
                                    content = '''#!/bin/bash
sh {} {} {}'''.format(pipe, ROOT, work_dir)
                                fp.write(content)

    def link(self):
        # 将Mother,Father,Embryo*链接到每组家系目录下面,优化了原始流程的重复分析
        for sampleID in self.info['sampleID'].unique():
            sample_info = self.info.loc[self.info['sampleID'] == sampleID]
            sample = sample_info['sample'].values[0]
            if sample in ['Mother', 'Father'] or sample.startswith('Embryo'):
                sample_dir = os.path.join(self.work_dir, sample)
            else:
                sample_dir = os.path.join(self.work_dir, sampleID, sample)
            if not os.path.exists(sample_dir):
                os.makedirs(sample_dir)
        sample_uni = self.info[['sample', 'sampleID']].drop_duplicates()
        parents = sample_uni.loc[sample_uni['sample'].isin(['Mother', 'Father']) | sample_uni['sample'].str.contains('Embryo')]['sample']
        grandparents = sample_uni.loc[~(sample_uni['sample'].isin(['Mother', 'Father']) | sample_uni['sample'].str.contains('Embryo'))]['sampleID']
        # 已链接的就跳过
        for sampleID in grandparents:
            for sample in parents:
                try:
                    os.symlink(f'{os.path.join(self.work_dir, sample)}', f'{os.path.join(self.work_dir, sampleID, sample)}')
                except FileExistsError:
                    continue

    @classmethod
    def create_job(cls, parsed_args):
        # 根据输入创建目录和脚本
        flow = pd.read_csv(parsed_args.step, sep='\t')
        work_dir = parsed_args.work_dir
        info = pd.read_csv(parsed_args.info, sep='\t')
        logger.info('Start Create Scripts!')
        cj = cls(flow, work_dir, info)
        cj.make(parsed_args.gene)
        cj.link()
        logger.info('Create Scripts Finished!')


class SubJobFrame:
    def __init__(self, flow, work_dir, info):
        # 同CreateJob,此处可改成类继承
        self.flow = flow
        self.flow.index = self.flow['Name'].tolist()
        self.work_dir = work_dir
        self.info = info
        self.info['sampleID'] = self.info['sample'] + '_' + self.info['sampleID']

    @staticmethod
    def initialize_job(job_graph, job_conf_dic, job, step):
        # 每个任务的初始化信息
        job_graph[job] = {}
        if os.path.exists(job + '.complete'):
            job_graph[job]['Status'] = 'complete'
        else:
            job_graph[job]['Status'] = 'incomplete'
        job_graph[job]['Name'] = step
        job_graph[job]['Type'] = job_conf_dic[step]['Type']
        job_graph[job]['Resources'] = job_conf_dic[step]['Resources']
        job_graph[job]['JobID'] = ''
        job_graph[job]['Children'] = []
        job_graph[job]['Parents'] = []
        return job_graph

    @staticmethod
    def find_parents(job_graph, job_conf_dic):
        # 根据配置文件和目录结构寻找每个任务的父任务
        # 使用了glob匹配,依赖特定目录结构和特定分析步骤名,很难写成万能通用的,此处修改应慎重
        for job in job_graph.keys():
            job_dir = os.path.dirname(job)
            job_name = os.path.basename(job)
            step = job_graph[job]['Name']
            job_type = job_conf_dic[step]['Type']
            parents = job_conf_dic[step]['Parents'].split(',')
            if job_type == 'single':
                for parent in parents:
                    job_name_split = job_name.split('-')
                    if parent != 'None':
                        order = f'step{job_conf_dic[parent]["Order"]}'
                        job_name_split[0], job_name_split[1] = order, parent
                        if parent == 'SortSam':
                            job_graph[job]['Parents'].extend(glob.glob(f'{job_dir}/{"-".join(job_name_split[0:3])}-*.sh'))
                        elif parent == 'BaseRecalibrator' and job_conf_dic[step]['SplitBy'] == 'None':
                            job_name_split[2] = job_name_split[2].replace('.sh', '')
                            job_graph[job]['Parents'].extend(glob.glob(f'{job_dir}/{"-".join(job_name_split[0:3])}-*.sh'))
                        else:
                            job_graph[job]['Parents'].append(os.path.join(job_dir, '-'.join(job_name_split)))
            else:
                for parent in parents:
                    job_name_split = job_name.split('-')
                    order = f'step{job_conf_dic[parent]["Order"]}'
                    job_name_split[0], job_name_split[1] = order, parent
                    if parent == 'HaplotypeCaller':
                        job_name_split.insert(2, '*')
                        job_graph[job]['Parents'].extend(glob.glob(f'{job_dir}/{"-".join(job_name_split)}'))
                        job_graph[job]['Parents'].extend(glob.glob(f'{os.path.dirname(os.path.dirname(job_dir))}/shell/{"-".join(job_name_split)}'))
                    elif parent == 'Stat':
                        job_graph[job]['Parents'].extend(glob.glob(f'{job_dir}/{"-".join(job_name_split[0:2])}-*.sh'))
                        job_graph[job]['Parents'].extend(glob.glob(f'{os.path.dirname(os.path.dirname(job_dir))}/shell/{"-".join(job_name_split[0:2])}-*.sh'))
                    elif parent == 'CombineGVCFs':
                        job_graph[job]['Parents'].append(os.path.join(job_dir, '-'.join(job_name_split)))
                    else:
                        job_graph[job]['Parents'].extend(glob.glob(f'{job_dir}/{"-".join(job_name_split[0:2])}*.sh'))
        return job_graph

    @staticmethod
    def find_children(job_graph):
        # 根据每个任务的父任务就能找到每个任务的所有子任务
        for job in job_graph.keys():
            for parent in job_graph[job]['Parents']:
                job_graph[parent]['Children'].append(job)
        return job_graph

    def generate_job_graph(self) -> dict:
        # 所有任务的关系字典
        job_graph = dict()
        job_conf_dic = self.flow.to_dict('index')
        for step in job_conf_dic.keys():
            for sampleID in self.info['sampleID'].unique():
                sample_info = self.info.loc[self.info['sampleID'] == sampleID]
                sample = sample_info['sample'].values[0]
                if sample in ['Mother', 'Father'] or sample.startswith('Embryo'):
                    script_dir = os.path.join(self.work_dir, 'shell')
                else:
                    script_dir = os.path.join(self.work_dir, sampleID, 'shell')
                if job_conf_dic[step]['Type'] == 'single':
                    if job_conf_dic[step]['SplitBy'] == 'barcode':
                        for index in sample_info.index:
                            fq_name = os.path.basename(sample_info.loc[index, 'fq1']).replace('_1.fq.gz', '')
                            job = os.path.join(script_dir, f'step{job_conf_dic[step]["Order"]}-{step}-{sample}-{fq_name}.sh')
                            job_graph = self.initialize_job(job_graph, job_conf_dic, job, step)
                    elif job_conf_dic[step]['SplitBy'] == 'chromosome':
                        for chromosome in chromosomes:
                            job = os.path.join(script_dir, f'step{job_conf_dic[step]["Order"]}-{step}-{sample}-{chromosome}.sh')
                            job_graph = self.initialize_job(job_graph, job_conf_dic, job, step)
                    else:
                        job = os.path.join(script_dir, f'step{job_conf_dic[step]["Order"]}-{step}-{sample}.sh')
                        job_graph = self.initialize_job(job_graph, job_conf_dic, job, step)
                else:
                    if sample in ['Mother', 'Father'] or sample.startswith('Embryo'):
                        continue
                    else:
                        if job_conf_dic[step]['SplitBy'] == 'chromosome':
                            for chromosome in chromosomes:
                                job = os.path.join(script_dir, f'step{job_conf_dic[step]["Order"]}-{step}-{chromosome}.sh')
                                job_graph = self.initialize_job(job_graph, job_conf_dic, job, step)
                        elif job_conf_dic[step]['SplitBy'] == 'mut_type':
                            for mt in mut_type:
                                job = os.path.join(script_dir, f'step{job_conf_dic[step]["Order"]}-{step}-{mt}.sh')
                                job_graph = self.initialize_job(job_graph, job_conf_dic, job, step)
                        else:
                            job = os.path.join(script_dir, f'step{job_conf_dic[step]["Order"]}-{step}.sh')
                            job_graph = self.initialize_job(job_graph, job_conf_dic, job, step)
        job_graph = self.find_parents(job_graph, job_conf_dic)
        job_graph = self.find_children(job_graph)
        return job_graph

    @staticmethod
    def job_num_in_sge():
        # sge中任务数
        command = "qstat | grep `whoami` |wc -l"
        status, output = subprocess.getstatusoutput(command)
        return status, output

    @staticmethod
    def job_id_in_sge(command):
        # 获取任务id
        status, output = subprocess.getstatusoutput(command)
        try:
            job_id = re.findall(r"Your job (\d+) ", output)[0]
        except IndexError:
            job_id = ''
        return status, job_id

    @staticmethod
    def job_status_in_sge(job_id):
        # 分析任务在sge中状态
        command = "qstat | grep " + "\"" + job_id + " " + "\""
        status, output = subprocess.getstatusoutput(command)
        return status, output

    @staticmethod
    def parents_status(job_graph, job):
        # 判断任务的全部父任务是否完成
        if len(job_graph[job]['Parents']) == 0:
            return 'complete'
        status_list = [job_graph[parent]['Status'] for parent in job_graph[job]['Parents']]
        if 'incomplete' not in status_list:
            return 'complete'
        else:
            return 'incomplete'

    @staticmethod
    def kill_job(job_graph, jobs):
        # 杀除任务
        for job in jobs:
            if os.path.exists(job + '.complete'):
                continue
            if job_graph[job]['JobID'] != '':
                _ = subprocess.getoutput('qdel ' + job_graph[job]['JobID'])
        logger.info('Running and pending jobs were killed!')

    def submit(self, job_graph, job):
        # 任务投递
        status, job_num = self.job_num_in_sge()
        if status != 0:
            logger.info('qstat error!')
            sys.exit(1)
        while int(job_num) >= 4000:
            time.sleep(600)
            status, job_num = self.job_num_in_sge()
            if status != 0:
                logger.info('qstat error!')
                sys.exit(1)
        job_path = os.path.dirname(job)
        command = "qsub -wd " + job_path + " " + job_graph[job]['Resources'] + " " + job
        times = 3
        status, job_id = self.job_id_in_sge(command)
        while times:
            if status == 0:
                break
            time.sleep(10)
            status, job_id = self.job_id_in_sge(command)
            times -= 1
        return status, job_id

    @classmethod
    def work_flow(cls, parsed_args):
        # 读取流程、工作路径、样本编号
        flow = pd.read_csv(parsed_args.step, sep='\t')
        work_dir = parsed_args.work_dir
        info = pd.read_csv(parsed_args.info, sep='\t')
        # 生成流程图
        sjf = cls(flow, work_dir, info)
        job_graph = sjf.generate_job_graph()
        with open(os.path.join(parsed_args.work_dir, 'job.yaml'), 'w') as f:
            yaml.dump(job_graph, f)
        logger.info('All Jobs Graph Created!')
        # 按流程图依赖顺序投递并监控任务状态
        jobs = list(filter(lambda x: len(job_graph[x]['Parents']) == 0, job_graph.keys()))
        while len(jobs) > 0:
            jobs_add = []
            jobs_remove = []
            for job in jobs:
                if os.path.exists(job + '.complete'):
                    job_graph[job]['Status'] = 'complete'
                if job_graph[job]['Status'] == 'incomplete':
                    if job_graph[job]['JobID'] == '':
                        if sjf.parents_status(job_graph, job) == 'complete':
                            status, job_id = sjf.submit(job_graph, job)
                            if status != 0:
                                logger.error(job + ' Submit Failed!')
                                sys.exit(1)
                            job_graph[job]['JobID'] = job_id
                            logger.info(job + ' Submit Success! JobID is ' + job_id)
                    else:
                        status, output = sjf.job_status_in_sge(job_graph[job]['JobID'])
                        if status != 0 and output == '' and not os.path.exists(job + '.complete'):
                            time.sleep(60)
                            if not os.path.exists(job + '.complete'):
                                logger.error(job + ' Run Failed!')
                                sjf.kill_job(job_graph, jobs)
                                sys.exit(1)
                else:
                    jobs_remove.append(job)
                    jobs_add += job_graph[job]['Children']
                    logger.info(job + ' Finished!')
            if len(jobs_add) == 0:
                time.sleep(60)
            jobs = list(set(list(set(jobs) - set(jobs_remove)) + jobs_add))
        logger.info('All Jobs Finished!')


def auto_remove(args):
    # 如果发现了新的胚胎就自动删除已完成的家系步骤的任务的.complete标识
    if not os.path.exists(f'{args.work_dir}/input.list'):
        logger.info(f'{args.work_dir}/input.list not exists and copy from {args.info}')
    else:
        info_new = pd.read_csv(args.info, sep='\t')
        info_new.drop_duplicates(subset=['sampleID'], inplace=True)
        info_new = info_new[info_new['sample'].str.contains('Embryo')]
        info_old = pd.read_csv(f'{args.work_dir}/input.list', sep='\t')
        info_old.drop_duplicates(subset=['sampleID'], inplace=True)
        info_old = info_old[info_old['sample'].str.contains('Embryo')]
        if info_old.shape[0] < info_new.shape[0]:
            steps = pd.read_csv(args.step, sep='\t')
            steps = steps[steps['Type'] == 'family']
            logger.info('New Embryo Available! Remove completed family jobs!')
            for step in steps['Name'].tolist():
                _ = subprocess.getoutput(f'rm {args.work_dir}/*/shell/*{step}*sh.complete')
                logger.info(f'rm {args.work_dir}/*/shell/*{step}*sh.complete')


def main():
    parser = ArgumentParser()
    parser.add_argument('-step', help='pipeline all steps', default=os.path.join(CONFIG, 'allsteps.tsv'))
    parser.add_argument('-work_dir', help='work dir', required=True)
    parser.add_argument('-info', help='sample info', required=True)
    parser.add_argument('-gene', help='detect gene', required=True)
    parser.add_argument('-create_only', help='only create scripts', action='store_true')
    parser.add_argument('-submit_only', help='only submit scripts', action='store_true')
    parsed_args = parser.parse_args()
    if not os.path.exists(parsed_args.work_dir):
        os.makedirs(parsed_args.work_dir)
    handler = logging.FileHandler(f'{parsed_args.work_dir}/auto.log')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    # 如果有新的胚胎样本则删除批次任务的.complete文件
    auto_remove(parsed_args)
    _ = subprocess.getoutput(f'cp {parsed_args.info} {parsed_args.work_dir}/input.list')
    if not parsed_args.submit_only:
        CreateJob.create_job(parsed_args)
    if parsed_args.create_only:
        logger.info('Only Create Scripts and Exit!')
        sys.exit(0)
    SubJobFrame.work_flow(parsed_args)


if __name__ == '__main__':
    main()
