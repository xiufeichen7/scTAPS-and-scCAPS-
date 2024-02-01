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
from plotnine import *
import time
from multiprocessing import Pool, Manager

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

def trans_query_by_cigar(read):
    """ 
    Return converted query_alignment_sequence by CIGAR.
    
    Args:
        Read: from samfile.fetch
    
    Returns:
        Convert seq: converted sequence (str)
    
    """
    cigartuples = pd.DataFrame(read.cigartuples, columns=["ops","len"])
    # divide query_sequence into match / insertion / deletion blocks; generate relative_position for each block
    # 0/M match; 1/I for indels; 2/D for del (todo: whether other POS need to be included)
    cigartuples.loc[cigartuples['ops']!=2, 'pos'] = list(itertools.accumulate(cigartuples[cigartuples['ops']!=2]['len']))  
    if 2 in list(cigartuples['ops']): cigartuples.loc[cigartuples['ops']==2, 'pos'] = cigartuples.loc[list(cigartuples['ops']).index(2) - 1, 'pos'] 
    # query_alignment_sequence exclude flanking bases that were soft clipped (None if not present)
    query_seq = read.query_alignment_sequence
    query_qual = list(read.query_alignment_qualities)
    # extract block sequence
    cigartuples['start'] = cigartuples.pos.shift(1, fill_value=0).astype('int')
    cigartuples['end'] = cigartuples.pos.shift(0, fill_value=list(cigartuples.pos)[-1]).astype('int')
    cigartuples['seq'] = cigartuples.apply(lambda x: query_seq[int(x['start']): int(x['end'])], axis=1)
    cigartuples['quals'] = cigartuples.apply(lambda x: str(query_qual[int(x['start']): int(x['end'])]), axis=1)   
    # replace insertion block as ""
    if 1 in list(cigartuples['ops']): cigartuples.loc[list(cigartuples['ops']).index(1),'seq'] = ""
    if 1 in list(cigartuples['ops']): cigartuples.loc[list(cigartuples['ops']).index(1),'quals'] = ""
    # replace deletion block as "D..D"
    if 2 in list(cigartuples['ops']): cigartuples.loc[list(cigartuples['ops']).index(2),'seq'] = cigartuples.loc[list(cigartuples['ops']).index(2),'len'] * "D"
    if 2 in list(cigartuples['ops']): cigartuples.loc[list(cigartuples['ops']).index(2),'quals'] = " "+ cigartuples.loc[list(cigartuples['ops']).index(2),'len'] * "0, "
    convert_seq = cigartuples['seq'].str.cat(sep='')
    convert_quals = re.sub(",[^0-9]*,",", ",cigartuples['quals'].str.strip('[]').str.cat(sep=', ')).split(", ")
    return(convert_seq, convert_quals)

def read_cg_var(read):
    """
    Return the read context in CG position
    
    Args:
        Read: from samfile.fetch
        
    Returns:
        Convert seq: sequece in CG pos (str)
        
    """
    querytrans = trans_query_by_cigar(read)[0]  # get converted query_alignment_sequence by CIGAR.
    iter = re.finditer("CG", read.get_reference_sequence().upper())  # find CG in reference sequence 
    # get the methylation context in CG position
    res = [[read.query_name, read.flag, read.cigarstring, m.start(), read.reference_name, m.start()+read.reference_start, querytrans[m.start():m.end()].upper()] for m in iter]
    return pd.DataFrame(res,columns=['query_name','flag', 'cigar', 'readpos','chr','pos','context'])

def bam_to_csv(bam, save_name, mq):
    """
    Return the read context in CG position
    Args:
        Read: bam
        
    Returns:       
    """
    pysam.index(bam)
    samfile=pysam.AlignmentFile(bam,"rb" )
    final_df = pd.DataFrame()
    for read in samfile.fetch():
        if read.mapping_quality > mq and str(read.flag) in ["83","99","147","163"]:
            final_df = pd.concat([read_cg_var(read), final_df])
    final_df.to_csv(save_name, sep='\t', index=False, header=False, compression="gzip")

def mk_temp_file(input_file_name, temp_dir, reads_num_mk):
    temp_file_list = []     # store the name of temp_file
    file_index = 1
    smp_name = os.path.basename(input_file_name)

    temp_file_name = os.path.join(temp_dir, f'{smp_name[:-4]}.temp.{file_index}.bam')
    temp_file_list.append(temp_file_name)
    # store whole input bam
    save_mk = pysam.set_verbosity(0)        # help suppressing HTSlib's message
    bam_input_all = pysam.AlignmentFile(input_file_name, "rb")
    pysam.set_verbosity(save_mk)
    temp_bam_header = bam_input_all.header  # record the header
    bam_temp = pysam.AlignmentFile(temp_file_name, "wb", header = temp_bam_header)
    all_bam_reads = bam_input_all.fetch(until_eof = True)
    last_read_name = ''                 # to store last read name, make sure paired read will be in the same file
    count = 1
    # read by read, write into a new temp file, 50000 reads per file
    for read in all_bam_reads:
        if count > reads_num_mk and read.query_name != last_read_name:
            # close this temp file
            count = 1
            bam_temp.close()
            # create the next temp file
            file_index += 1
            temp_file_name = os.path.join(temp_dir, f'{smp_name[:-4]}.temp.{file_index}.bam')
            temp_file_list.append(temp_file_name)
            bam_temp = pysam.AlignmentFile(temp_file_name, "wb", header = temp_bam_header)
            bam_temp.write(read)
        else:
            count += 1
            bam_temp.write(read)
            last_read_name = read.query_name
    # close the last bam_temp
    bam_temp.close()
    return temp_file_list, file_index


# paramters
parser = argparse.ArgumentParser(
    description="Add....",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
)
parser.add_argument("-b", help="bam....", type=str, dest="bam_path")
parser.add_argument("-t", help="core....", type=int, default=8, dest="core_num")
parser.add_argument("-q", help="mapping quality....", type=int, default=10, dest="mq")
parser.add_argument("-n", help="read per file ....", type=int, default=30000, dest="per_n")

args = parser.parse_args()


# get the output name
bam =  args.bam_path #"scTAPS_i7_8_i5_2_S58_0010.temp.6.bam"#
core_num = args.core_num
per_n = args.per_n
mq = args.mq

save_name1 = bam.replace("bam", "read.meth.txt.gz")
save_name2 = bam.replace("bam", "all.sta.txt.gz")
save_name3 = bam.replace("bam", "meth.sta.txt.gz")

# processing
print("start:", time.strftime("%Y-%m-%d %H:%M:%S"))
pysam.index(bam)

print("""processing bam ... """)
temp_dir=os.path.basename(bam).replace("bam", "call_temp")
os.makedirs(temp_dir)

mk_temp_result = mk_temp_file(bam, temp_dir, per_n)
temp_file = mk_temp_result[0]
temp_out_file = []
temp_out_filter_file = []
for ii in range(len(temp_file)):        
    temp_out_file.append(temp_file[ii][:-4] + '.read.meth.txt.gz')

po = Pool(core_num)
results = []
for ii in range(len(temp_file)):
    result=po.apply_async(bam_to_csv, (temp_file[ii], temp_out_file[ii], mq))
    results.append(result)

po.close()
po.join() 


for result in results:
    print(result.get())

print("""write output ...""")

import gzip
import shutil

# output file name
header = b"query_name\tflag\tcigar\treadpos\tchr\tpos\tcontext\n"
with gzip.open(save_name1, 'wb') as outfile:
    outfile.write(header)
    for infile_name in temp_out_file:
        with gzip.open(infile_name, 'rb') as infile:
            shutil.copyfileobj(infile, outfile)

shutil.rmtree(temp_dir)

print("""summarize methylation ...""")
final_df = pd.read_csv(save_name1, compression='gzip', sep="\t")



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

table_pos_context.to_csv(save_name2, sep='\t', index=False)
table_pos_context[['chr','pos','mC','aC','ratio']].to_csv(save_name3, sep='\t', index=False)

print("""finished ...""")
print("start:", time.strftime("%Y-%m-%d %H:%M:%S"))
