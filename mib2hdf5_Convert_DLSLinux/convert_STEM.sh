#!/bin/bash
#$ -l h_rt=1:00:00
#$ -cwd
#$ -l redhat_release=rhel7
#$ -l m_mem_free=100G

module load python/3.7
python /dls/science/groups/e02/code/python_conversion_quadMedipix/working_versions/v1/development/mib2hdf_watch04_np_TEM.py e02 2019 cm22979-3
