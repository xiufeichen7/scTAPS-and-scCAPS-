import sys
import pysam
import os
import re
import numpy as np
import itertools
import importlib.util
import pandas as pd
import warnings
import argparse
import time
import glob
import shutil
import gzip

warnings.filterwarnings('ignore')




# functions
def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

# paramters
parser = argparse.ArgumentParser(
    description="Add....",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
parser.add_argument("-r", help="read meth....", type=str, dest="read_meth")
parser.add_argument("-s", help="start base....", type=int, dest="start_base")
parser.add_argument("-e", help="end base....", type=int, dest="end_base")
parser.add_argument("-c", help="clip overlap....", type=bool, default=True, dest="clip_overlap")
args = parser.parse_args()


# get the output name
read_meth =  args.read_meth
start_base =  args.start_base 
end_base =  args.end_base 

read_filter = read_meth.replace("read.meth.txt.gz", "filter.read.meth.txt.gz")
all_sta = read_meth.replace("read.meth.txt.gz", "filter.all.sta.txt.gz")
meth_sta = read_meth.replace("read.meth.txt.gz", "filter.meth.sta.txt.gz")


print("""summarize methylation ...""", read_meth)
final_df = pd.read_csv(read_meth, compression='gzip', sep="\t")
final_df = final_df[(final_df['readpos']>start_base) & (final_df['readpos']<end_base)]
final_df = final_df.sort_values(by=['query_name', 'flag'], ascending=[True, True])

if args.clip_overlap:
    print("""before clipping""", len(final_df))
    final_df = final_df.drop_duplicates(subset=['query_name', 'chr', 'pos'], keep='first')
    print("""after clipping""", len(final_df))
    final_df.to_csv(read_filter, sep='\t', index=False)

table_pos_context = final_df.groupby(['chr','pos','context'],as_index=False).count()
table_pos_context =table_pos_context.pivot_table(index=['chr','pos'],
                                                columns='context',
                                                values='query_name',
                                                aggfunc=sum,
                                                fill_value=0).reset_index()
if all(table_pos_context.columns != "CA"):
    table_pos_context['CA'] = 0

if all(table_pos_context.columns != "TG"):
    table_pos_context['TG'] = 0 

table_pos_context['mC'] = table_pos_context['CA']+table_pos_context['TG']
table_pos_context['aC'] = table_pos_context['CA']+table_pos_context['TG'] + table_pos_context['CG']                             
table_pos_context['ratio'] = round(table_pos_context['mC'] / table_pos_context['aC'] *100,2)


table_pos_context.to_csv(all_sta, sep='\t', index=False)
table_pos_context[['chr','pos','mC','aC','ratio']].to_csv(meth_sta, sep='\t', index=False)
