#!/bin/bash
#$ -l h_rt=0:20:00
#$ -cwd
#$ -l redhat_release=rhel7
#$ -l m_mem_free=100G

module load python/3.7
python /home/eha56862/code/DLS_cluster_conversion/mib2hdf_watch_convert.py e02 2017 em16952-1 1 512 'processing/raw'
