import vcfpy
import pandas as pd
import numpy as np

import argparse
import re
from glob import glob
import os

ALPHABET = pd.Series([0,1,2,3,4], index=["A","C","T","G","N"])




def grab_vcf(vcffile,name=""):
    
    reader = vcfpy.Reader.from_path(vcffile)

    
    #Records:
    #Ref, Position, alt, depth, AF, SB, DP4[-2:]
    #More than one ALT - tough

    out_records = []

    
    for record in reader:

        rec = [record.CHROM, record.POS, record.ALT[0].value, record.INFO["AF"],record.INFO["DP"], record.INFO["SB"], *record.INFO["DP4"][-2:]]

        

        out_records.append(rec)

    
    df = pd.DataFrame(out_records, columns=["chrom","pos", "alt", "af", "depth", "sb", "alt_fw","alt_rv"])

    df["name"] = name



    

    return df






def tally_pos(df_list):

    if type(df_list) is pd.core.frame.DataFrame:
        pos = df_list.pos.drop_duplicates().sort_values().values
    else:
        pos = sorted(set(np.concatenate([df.pos.values for df in df_list])))

    return pos




def read_vcfs(filelist=None, path=None):
    
    vcf_files = []

    if path is not None:
        vcf_files.extend(glob("{}/*.vcf*".format(path)))

    if filelist is not None:
        vcf_files.extend(filelist)
    
    print(vcf_files)
    re_sub = re.compile("\\.vcf.*$")
    name_list = [re.sub(re_sub, "", os.path.basename(filename)) for filename in vcf_files]


    print(vcf_files)
    print(name_list)


    vcf_df = pd.concat([grab_vcf(filename,name=name) for name,filename in zip(name_list,vcf_files)])



    return vcf_df


def filter_vcf(vcf_df, strand_min_count=1, strand_match_min=5, min_frac=0.02, min_depth=30,SB=60):

    vcf_df_trim = vcf_df[(vcf_df.alt_fw > 0) & (vcf_df.alt_rv > 0) & \
            ((vcf_df.alt_fw >= strand_match_min) | (vcf_df.alt_rv >= strand_match_min)) & \
            (vcf_df.sb < SB)] 

    return vcf_df_trim



def vcf_matrix(vcf_df, template, pos_index, col_index):

    template = np.copy(template[:])

    template[:,4] = 1

    template[pos_index[vcf_df.pos].values, col_index[vcf_df.alt].values] = vcf_df.af.values

    template_sum = np.sum(template[:,:4],axis=1)

    template[:,4] = template[:,4] - template_sum



    return template

def vcf_matrices(vcf_df):

    pos = np.array(sorted(tally_pos(vcf_df)))

    pos_index = pd.Series(np.arange(len(pos)), index=pos)

    matrix_template = np.zeros(5*len(pos)).reshape(len(pos),5)

    col_index = pd.Series([0,1,2,3,4], index=["A","C","T","G","X"])

    g = vcf_df.groupby("name")

    vcf_matrix_dict = {}
    for name, group in g:

        vm = vcf_matrix(group, matrix_template, pos_index, col_index)

        vcf_matrix_dict[name] = vm[:]


    return vcf_matrix_dict





















if __name__ == "__main__":
    parser = argparse.ArgumentParser()


    parser.add_argument("-d", dest="vcfdir", type=str, help="Folder containing VCF files")
    parser.add_argument("-v", dest="vcffiles", type=str, help="Single or list of vcf files (separated by a comma)")
    parser.add_argument("-o", dest="output", type=str, help="Folder containing VCF files")





