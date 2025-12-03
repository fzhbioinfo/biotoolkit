#!/bin/bash
export UTA_DB_URL=postgresql://uta_admin:bgi_2020@10.2.1.4:6543/uta/uta_20190920
export HGVS_SEQREPO_DIR=/ifs7/B2C_RD_P2/SHARE/seqrepo/2020-04-13
/ifs7/B2C_RD_P2/USER/fangzhonghai/software/anaconda3/envs/hgvs/bin/python Coordinate_Converter.py -i g.txt -o g.Coordinate_Converter.txt -t g
