"""自动投递shell脚本到队列

Usage:
    python auto_qsub.py <shell_file>

Example:
    python auto_qsub.py shell.txt

shell_file:
    shell脚本文件，每行一个shell脚本
"""
import sys
import time
import subprocess
import pandas as pd


def job_num_in_sge():
    # SGE任务数目
    command = "qstat | grep `whoami` |wc -l"
    status, output = subprocess.getstatusoutput(command)
    return status, output

def main():
    shell = pd.read_csv(sys.argv[1], sep="\t", header=None)[0].tolist()
    qsub = f"qsub -cwd -l vf={sys.argv[2]}G,p={sys.argv[3]},h='!tj-compute-44-4&!tj-compute-31-6' -P P23Z15000N0119 -S /bin/bash -q b2c_rd1.q"
    for i in shell:
        _, job_num = job_num_in_sge()
        while int(job_num) > 1000:
            time.sleep(100)
            _, job_num = job_num_in_sge()
        subprocess.getoutput(f'{qsub} {i}')


if __name__ == '__main__':
    main()
