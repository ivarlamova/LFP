# -*- coding: utf-8 -*-
"""
Control of project
"""
from __future__ import division, unicode_literals, absolute_import 

"""General parameters"""
PATH2POTENTIALS = '/home/irina/potcars/'
pmgkey = "AWqKPyV8EmTRlf1t" # please get your own key from materials project for pymatgen
PATH2DATABASE        = '/home/aksenov/Data/CEStorage/'
PATH2JMOL     = 'jmol'
PATH2PHONOPY  = 'phonopy'
PATH2NEBMAKE = '~/Simulation_wrapper/vts/nebmake.pl' # http://theory.cm.utexas.edu/vasp/scripts.html
geo_folder           = './' # duplicate structures using .geo format files in this folder



"""Cluster parameters"""
PATH2PROJECT = '' # path to project on cluster relative to home folder
DEFAULT_CLUSTER = 'pardus'
PBS_PROCS = True # if true than #PBS -l procs="+str(number_cores) is used
WALLTIME_LIMIT = True

#description of user clusters:
CLUSTERS = {}

CLUSTERS['pardus'] = {
'address':'i.varlamova@10.30.17.12',
'vasp_com':'mpirun  vasp_std',
'schedule':'PBS',
'walltime':'48:00:00',
'corenum':16,
'pythonpath':'/usr/lib64/python2.7/site-packages/numpy',
# 'modules':'module load Compilers/Intel/psxe_2015.6; module load MPI/intel/5.1.3.258/intel; module load QCh/VASP/5.4.1p1/psxe2015.6; module load ScriptLang/python/2.7\n',
'modules':'module load Compilers/Intel/psxe_2017.4; module load MPI/hpcx-openmpi/1.9.5/Intel/psxe_2017.4; module load QCh/VASP/5.4.4/psxe2017.4; module load ScriptLang/python/2.7\n',
}

CLUSTERS['cee'] = {
'address':'i.varlamova@10.30.16.168',
'vasp_com':'mpirun vasp_std',
'homepath':'/home/i.varlamova/',
'schedule':'SLURM',
'corenum':16,
'modules':'module load Compiler/Intel/17u8; module load Q-Ch/VASP/5.4.4_OPT; module load ScriptLang/python/3.6i_2018u3; \nulimit -s unlimited\n',
}

CLUSTERS['1'] = {
'address':'username@ip',
'vasp_com':'prun /opt/vasp/bin/vasp',
'homepath':'/home/username/',
'schedule':'SLURM',
'corenum':16,
'pythonpath':'/usr/lib64/python2.7/site-packages/numpy'
}

CLUSTERS['2'] = {
'address':'username@ip',
'vasp_com':'mpirun  vasp_std',
'homepath':'/home/username/',
'schedule':'PBS',
'walltime':'72:00:00',
'corenum':16,
'pythonpath':'/usr/lib64/python2.7/site-packages/numpy',
'modules':'module load Compilers/Intel/psxe_2015.6; module load MPI/intel/5.1.3.258/intel; module load QCh/VASP/5.4.1p1/psxe2015.6; module load ScriptLang/python/2.7',
}




"""Other parameters for developers, can be removed"""
NEW_BATCH = True   # for testing new batch system based on set sequences 
# RAMDISK              = '/mnt/ramdisk/'
RAMDISK              = None
EXCLUDE_NODES  = 1
CIF2CELL = True 
path_to_wrapper      = '/home/aksenov/Simulation_wrapper/'
path_to_paper        = '/home/aksenov/Research/CEStorage/aksenov_report/'
path_to_images       = path_to_paper+'/fig/'


"""List of manually added calculations:"""
#Deprecated, now not really needed
MANUALLY_ADDED = [# calc name, calc folder, calc des  
    ( 'Li111'        ,"Li",        "2 Li"                                  ),
    ( 'Rb111'        ,"Rb/bcc",        "2 Rb"                                  ),
    ]

