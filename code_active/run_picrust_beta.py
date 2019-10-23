#!/usr/bin/python

import shutil
import urllib2
from contextlib import closing
import gzip
import time
import os
import getopt
import sys
from subprocess import call

import numpy as np
import pandas as pd

from itertools import izip
from StringIO import StringIO
from collections import defaultdict

from lxml import etree
from skbio.parse.sequences import parse_fastq, parse_fasta

from joblib import Parallel, delayed  
import multiprocessing
import subprocess

import qiime_default_reference as qdr


def main():
   data_path = '/data/sw1/Dropbox/stm_microbiome/data_active'
   folder_name = str(sys.argv[1])
   print os.path.join(data_path,folder_name,'beta_table.biom')

   otu_table = os.path.join(data_path,folder_name,'beta_table.biom')
   pred_metagenome = os.path.join(data_path,folder_name,'beta_predicted_metagenome.biom')
   pred_cogs = os.path.join(data_path,folder_name,'beta_predicted_cogs.biom')
   beta_div_path = os.path.join(data_path,folder_name,'beta_diversity')

   call(["!predict_metagenomes.py","-i",_otu_table,"-o",_pred_metagenome])
   call(["!predict_metagenomes.py","--type_of_prediction","cog","-i",_otu_table,"-o",_pred_cogs])
   call(["!beta_diversity.py","-i",_otu_table,"-o",_beta_div_path,"-t",_tree])

   otu_table = os.path.join(data_path,folder_name,'beta_ppd_table.biom')
   pred_metagenome = os.path.join(data_path,folder_name,'beta_ppd_predicted_metagenome.biom')
   pred_cogs = os.path.join(data_path,folder_name,'beta_ppd_predicted_cogs.biom')
   beta_div_path = os.path.join(data_path,folder_name,'beta_diversity')

   call(["!predict_metagenomes.py","-i",_otu_table,"-o",_pred_metagenome])
   call(["!predict_metagenomes.py","--type_of_prediction","cog","-i",_otu_table,"-o",_pred_cogs])
   call(["!beta_diversity.py","-i",_otu_table,"-o",_beta_div_path,"-t",_tree])

   otu_table = os.path.join(data_path,folder_name,'beta_repl_table.biom')
   pred_metagenome = os.path.join(data_path,folder_name,'beta_repl_predicted_metagenome.biom')
   pred_cogs = os.path.join(data_path,folder_name,'beta_repl_predicted_cogs.biom')
   beta_div_path = os.path.join(data_path,folder_name,'beta_diversity')

   call(["!predict_metagenomes.py","-i",_otu_table,"-o",_pred_metagenome])
   call(["!predict_metagenomes.py","--type_of_prediction","cog","-i",_otu_table,"-o",_pred_cogs])
   call(["!beta_diversity.py","-i",_otu_table,"-o",_beta_div_path,"-t",_tree])

   otu_table = os.path.join(data_path,folder_name,'beta_table_permuted.biom')
   pred_metagenome = os.path.join(data_path,folder_name,'beta_permuted_predicted_metagenome.biom')
   pred_cogs = os.path.join(data_path,folder_name,'beta_permuted_predicted_cogs.biom')
   beta_div_path = os.path.join(data_path,folder_name,'beta_diversity')

   call(["!predict_metagenomes.py","-i",_otu_table,"-o",_pred_metagenome])
   call(["!predict_metagenomes.py","--type_of_prediction","cog","-i",_otu_table,"-o",_pred_cogs])

   otu_table = os.path.join(data_path,folder_name,'beta_ppd_table_permuted.biom')
   pred_metagenome = os.path.join(data_path,folder_name,'beta_ppd_permuted_predicted_metagenome.biom')
   pred_cogs = os.path.join(data_path,folder_name,'beta_ppd_permuted_predicted_cogs.biom')
   beta_div_path = os.path.join(data_path,folder_name,'beta_diversity')

   call(["!predict_metagenomes.py","-i",_otu_table,"-o",_pred_metagenome])
   call(["!predict_metagenomes.py","--type_of_prediction","cog","-i",_otu_table,"-o",_pred_cogs])

if __name__ == "__main__":
   main()
