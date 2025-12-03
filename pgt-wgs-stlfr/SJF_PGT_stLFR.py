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

# 配置文件目录和每步流程目录
ROOT = os.path.abspath(os.path.dirname(__file__))
CONFIG = os.path.join(ROOT, 'etc')
PIPELINE = os.path.join(ROOT, 'script')

# 按染色体拆分
chromosomes = pd.read_csv(f'{CONFIG}/chr.list', sep='\t', header=None)
chromosomes = chromosomes[0].tolist()


class SubJobFrame:
    def __init__(self, flow, work_dir, info):
        # 类实例化,读取流程步骤
        self.flow = pd.read_csv(flow, sep='\t')
        # 分析目录和样本信息
        self.work_dir = work_dir
        self.info = pd.read_csv(info, sep='\t')
        # 父母缺少任何一人
        if 'Mother' not in self.info['sample'].values or 'Father' not in self.info['sample'].values:
            self.flow = self.flow.loc[~self.flow['Name'].str.contains('Stat')].copy()
            if 'yes' in self.info['is_lfr'].values:
                self.flow = self.flow.loc[self.flow['Type'] != 'family'].copy() if 'no' in self.info['is_lfr'].values else self.flow.loc[(self.flow['Type'] != 'family') & (self.flow['is_lfr'] == 'yes')].copy()
            else:
                self.flow = self.flow.loc[(self.flow['Type'] != 'family') & (self.flow['is_lfr'] != 'yes')].copy()
        # 只有LFR数据则只分析LFR相关步骤和家系步骤
        elif 'yes' in self.info['is_lfr'].values and 'no' not in self.info['is_lfr'].values:
            self.flow = self.flow.loc[(self.flow['Type'] == 'family') | (self.flow['is_lfr'] == 'yes')].copy()
        self.flow.index = self.flow['Name'].tolist()
        # 流程配置转换成字典,去除不存在的父任务
        self.job_conf_dic = self.flow.to_dict('index')
        for step in self.job_conf_dic:
            if 'None' in self.job_conf_dic[step]['Parents']:
                continue
            self.job_conf_dic[step]['Parents'] = ','.join(set(self.job_conf_dic[step]['Parents'].split(',')) & set(self.job_conf_dic.keys()))

    def make(self):
        # 任务字典
        job_graph = dict()
        # 脚本的创建的目录
        script_dir = os.path.join(self.work_dir, 'shell')
        if not os.path.exists(script_dir):
            os.makedirs(script_dir)
        # 遍历流程每个步骤
        for step in self.job_conf_dic.keys():
            pipe = os.path.join(PIPELINE, f'{step}.sh')
            # 是否为单样本步骤
            if self.job_conf_dic[step]['Type'] == 'single':
                # 遍历每个样本
                for sample in self.info['sample'].unique():
                    # 样本的输入信息
                    sample_info = self.info.loc[self.info['sample'] == sample]
                    # LFR相关步骤
                    if self.job_conf_dic[step]['is_lfr'] == 'yes':
                        # 如果不是LFR样本则此样本跳过该步骤
                        if sample_info['is_lfr'].values[0] != 'yes':
                            continue
                        # LFR样本第一步需要输入fq信息
                        if self.job_conf_dic[step]['Parents'] == 'None':
                            fq1 = ','.join(sample_info['fq1'].tolist())
                            fq2 = fq1.replace('_1.fq.gz', '_2.fq.gz')
                            job = f'{script_dir}/{step}-{sample}.sh'
                            job_graph = self.initialize_job(job_graph, job, step)
                            with open(job, 'w') as fp:
                                content = f'sh {pipe} {ROOT} {self.work_dir} {sample} {fq1} {fq2}'
                                fp.write(content)
                        else:
                            # 按染色体进行拆分
                            if self.job_conf_dic[step]['SplitBy'] == 'chromosome':
                                for chromosome in chromosomes:
                                    job = f'{script_dir}/{step}-{sample}-{chromosome}.sh'
                                    job_graph = self.initialize_job(job_graph, job, step)
                                    with open(job, 'w') as fp:
                                        content = f'sh {pipe} {ROOT} {self.work_dir} {sample} {chromosome}'
                                        fp.write(content)
                            else:
                                # 不用拆分的LFR单样本步骤
                                job = f'{script_dir}/{step}-{sample}.sh'
                                job_graph = self.initialize_job(job_graph, job, step)
                                with open(job, 'w') as fp:
                                    content = f'sh {pipe} {ROOT} {self.work_dir} {sample}'
                                    fp.write(content)
                    # 非LFR相关步骤
                    else:
                        # 如果是LFR样本则此样本跳过该步骤
                        if sample_info['is_lfr'].values[0] == 'yes':
                            continue
                        # 按barcode维度进行拆分
                        if self.job_conf_dic[step]['SplitBy'] == 'barcode':
                            for index in sample_info.index:
                                fq1 = sample_info.loc[index, 'fq1']
                                fq2 = fq1.replace('_1.fq.gz', '_2.fq.gz')
                                fq_name = os.path.basename(fq1).replace('_1.fq.gz', '')
                                job = f'{script_dir}/{step}-{sample}-{fq_name}.sh'
                                job_graph = self.initialize_job(job_graph, job, step)
                                with open(job, 'w') as fp:
                                    content = f'sh {pipe} {ROOT} {self.work_dir} {sample} {fq1} {fq2}'
                                    fp.write(content)
                        # 按染色体维度拆分
                        elif self.job_conf_dic[step]['SplitBy'] == 'chromosome':
                            for chromosome in chromosomes:
                                job = f'{script_dir}/{step}-{sample}-{chromosome}.sh'
                                job_graph = self.initialize_job(job_graph, job, step)
                                with open(job, 'w') as fp:
                                    content = f'sh {pipe} {ROOT} {self.work_dir} {sample} {chromosome}'
                                    fp.write(content)
                        # 不用拆分的非LFR单样本步骤
                        else:
                            job = f'{script_dir}/{step}-{sample}.sh'
                            job_graph = self.initialize_job(job_graph, job, step)
                            with open(job, 'w') as fp:
                                content = f'sh {pipe} {ROOT} {self.work_dir} {sample}'
                                fp.write(content)
            # 家系分析步骤
            else:
                # 按染色体维度拆分
                if self.job_conf_dic[step]['SplitBy'] == 'chromosome':
                    for chromosome in chromosomes:
                        job = f'{script_dir}/{step}-{chromosome}.sh'
                        job_graph = self.initialize_job(job_graph, job, step)
                        with open(job, 'w') as fp:
                            content = f'sh {pipe} {ROOT} {self.work_dir} {chromosome}'
                            fp.write(content)
                # 不用拆分
                else:
                    job = f'{script_dir}/{step}.sh'
                    job_graph = self.initialize_job(job_graph, job, step)
                    with open(job, 'w') as fp:
                        content = f'sh {pipe} {ROOT} {self.work_dir}'
                        fp.write(content)
        return job_graph

    def initialize_job(self, job_graph, job, step):
        # 每个任务的初始化信息
        job_graph[job] = {}
        if os.path.exists(job + '.complete'):
            job_graph[job]['Status'] = 'complete'
        else:
            job_graph[job]['Status'] = 'incomplete'
        job_graph[job]['Name'] = step
        job_graph[job]['Type'] = self.job_conf_dic[step]['Type']
        job_graph[job]['is_lfr'] = self.job_conf_dic[step]['is_lfr']
        job_graph[job]['Resources'] = self.job_conf_dic[step]['Resources']
        job_graph[job]['JobID'] = ''
        job_graph[job]['Children'] = []
        job_graph[job]['Parents'] = []
        return job_graph

    def find_parents(self, job_graph):
        # 根据配置文件和目录结构寻找每个任务的父任务,SplitBy: None,chromosome,barcode
        script_dir = os.path.join(self.work_dir, 'shell')
        for job in job_graph.keys():
            step = job_graph[job]['Name']
            # 单样本任务,其父任务只会是single类型
            if job_graph[job]['Type'] == 'single':
                for parent in self.job_conf_dic[step]['Parents'].split(','):
                    job_name_split = os.path.basename(job).split('-')
                    if parent == 'None':
                        continue
                    # 根据step的SplitBy进行分情况讨论
                    if self.job_conf_dic[step]['SplitBy'] == 'None':
                        # 根据parent的SplitBy进行分情况讨论
                        if self.job_conf_dic[parent]['SplitBy'] == 'None':
                            job_name_split[0] = parent
                            # 非LFR步骤依赖LFR步骤
                            if self.job_conf_dic[step]['is_lfr'] != 'yes' and self.job_conf_dic[parent]['is_lfr'] == 'yes':
                                job_graph[job]['Parents'].extend(glob.glob(f'{script_dir}/{job_name_split[0]}-*.sh'))
                            else:
                                job_graph[job]['Parents'].append(f'{script_dir}/{"-".join(job_name_split)}')
                        elif self.job_conf_dic[parent]['SplitBy'] == 'chromosome':
                            job_name_split[0], job_name_split[1] = parent, job_name_split[1].replace('.sh', '')
                            job_graph[job]['Parents'].extend(glob.glob(f'{script_dir}/{"-".join(job_name_split)}-*.sh'))
                        else:
                            raise Exception('unexpect condition')
                    elif self.job_conf_dic[step]['SplitBy'] == 'chromosome':
                        # 根据parent的SplitBy进行分情况讨论
                        if self.job_conf_dic[parent]['SplitBy'] == 'chromosome':
                            job_name_split[0] = parent
                            job_graph[job]['Parents'].append(f'{script_dir}/{"-".join(job_name_split)}')
                        elif self.job_conf_dic[parent]['SplitBy'] == 'barcode':
                            job_name_split[0] = parent
                            job_graph[job]['Parents'].extend(glob.glob(f'{script_dir}/{"-".join(job_name_split[0:2])}-*.sh'))
                        else:
                            job_name_split[0] = parent
                            job_graph[job]['Parents'].append(f'{script_dir}/{"-".join(job_name_split[0:2])}.sh')
                    else:
                        if self.job_conf_dic[parent]['SplitBy'] == 'barcode':
                            job_name_split[0] = parent
                            job_graph[job]['Parents'].append(f'{script_dir}/{"-".join(job_name_split)}')
                        else:
                            raise Exception('unexpect condition')
            # 家系任务,SplitBy只有None,chromosome
            else:
                for parent in self.job_conf_dic[step]['Parents'].split(','):
                    job_name_split = os.path.basename(job).split('-')
                    # parent为单样本任务
                    if self.job_conf_dic[parent]['Type'] == 'single':
                        if self.job_conf_dic[step]['SplitBy'] == 'None' and self.job_conf_dic[parent]['SplitBy'] == 'None':
                            job_name_split[0] = parent
                            job_graph[job]['Parents'].extend(glob.glob(f'{script_dir}/{"-".join(job_name_split)}-*.sh'))
                        elif self.job_conf_dic[step]['SplitBy'] == 'chromosome' and self.job_conf_dic[parent]['SplitBy'] == 'chromosome':
                            job_name_split[0] = parent
                            job_name_split.insert(1, '*')
                            job_graph[job]['Parents'].extend(glob.glob(f'{script_dir}/{"-".join(job_name_split)}'))
                        else:
                            raise Exception('unexpect condition')
                    # parent也为家系任务
                    else:
                        if self.job_conf_dic[step]['SplitBy'] == 'None':
                            if self.job_conf_dic[parent]['SplitBy'] == 'None':
                                job_name_split[0] = parent
                                job_graph[job]['Parents'].append(f'{script_dir}/{"-".join(job_name_split)}.sh')
                            else:
                                job_name_split[0] = parent
                                job_graph[job]['Parents'].extend(glob.glob(f'{script_dir}/{"-".join(job_name_split)}-*.sh'))
                        else:
                            if self.job_conf_dic[parent]['SplitBy'] == 'chromosome':
                                job_name_split[0] = parent
                                job_graph[job]['Parents'].append(f'{script_dir}/{"-".join(job_name_split)}')
                            else:
                                raise Exception('unexpect condition')
        return job_graph

    @staticmethod
    def find_children(job_graph):
        # 根据每个任务的父任务就能找到每个任务的所有子任务
        for job in job_graph.keys():
            for parent in job_graph[job]['Parents']:
                job_graph[parent]['Children'].append(job)
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
    def work_flow(cls, args):
        # 生成流程图
        sjf = cls(args.step, args.work_dir, args.info)
        logger.info('Start Create Scripts!')
        job_graph = sjf.make()
        logger.info('Create Scripts Finished!')
        job_graph = sjf.find_parents(job_graph)
        job_graph = sjf.find_children(job_graph)
        with open(os.path.join(args.work_dir, 'job.yaml'), 'w') as f:
            yaml.dump(job_graph, f)
        logger.info('All Jobs Graph Created!')
        if args.create_only:
            logger.info('Only Create Scripts and Exit!')
            sys.exit(0)
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
        info_new.drop_duplicates(subset=['sample'], inplace=True)
        info_new = info_new[info_new['sample'].str.contains('Embryo')]
        info_old = pd.read_csv(f'{args.work_dir}/input.list', sep='\t')
        info_old.drop_duplicates(subset=['sample'], inplace=True)
        info_old = info_old[info_old['sample'].str.contains('Embryo')]
        if info_old.shape[0] < info_new.shape[0]:
            steps = pd.read_csv(args.step, sep='\t')
            steps = steps[(steps['Type'] == 'family') & (steps['is_lfr'] != 'yes')]
            logger.info('New Embryo Available! Remove completed family jobs!')
            for step in steps['Name'].tolist():
                _ = subprocess.getoutput(f'rm {args.work_dir}/shell/*{step}*sh.complete')
                logger.info(f'rm {args.work_dir}/shell/*{step}*sh.complete')


def main():
    parser = ArgumentParser()
    parser.add_argument('-step', help='pipeline all steps', default=os.path.join(CONFIG, 'allsteps.tsv'))
    parser.add_argument('-work_dir', help='work dir', required=True)
    parser.add_argument('-info', help='sample info', required=True)
    parser.add_argument('-create_only', help='only create scripts', action='store_true')
    parsed_args = parser.parse_args()
    if not os.path.exists(parsed_args.work_dir):
        os.makedirs(parsed_args.work_dir)
    handler = logging.FileHandler(f'{parsed_args.work_dir}/auto.log')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    # 如果有新的胚胎样本则删除批次任务的.complete文件
    auto_remove(parsed_args)
    _ = subprocess.getoutput(f'cp {parsed_args.info} {parsed_args.work_dir}/input.list')
    # 生成基因list
    info = pd.read_csv(parsed_args.info, sep='\t')
    info = info.loc[info['gene'] != '.']
    genes = list()
    for gene in info['gene'].unique():
        genes.extend(gene.split(','))
    pd.DataFrame(set(genes)).to_csv(f'{parsed_args.work_dir}/gene.list', sep='\t', header=None, index=False)
    SubJobFrame.work_flow(parsed_args)


if __name__ == '__main__':
    main()
