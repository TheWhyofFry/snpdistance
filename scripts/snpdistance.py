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
    
    re_sub = re.compile("\\.vcf.*$")
    name_list = [re.sub(re_sub, "", os.path.basename(filename)) for filename in vcf_files]




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

def badpos(vcf_folder, pos_set, min_depth=100):

    reg_ext = re.compile("\\.depth.+")
    filelist = glob("{vcf}/*depth.gz".format(vcf=vcf_folder))


    pos_set = set(pos_set) 
    basename = [re.sub(reg_ext,"",os.path.basename(filename)) for filename in filelist]

    depth_df_list = []

    for base, filename in zip(basename, filelist):

        df = pd.read_table(filename, sep="\t", names=["ref","pos","depth"])
        df = df[(df.depth < min_depth)]

        df["sample"] = base

        depth_df_list.append(df)
    
    return pd.concat(depth_df_list)

                



    










def dist(m1,m2,badpositions=None):
    if badpos is None:
        return np.sum(np.abs(m1-m2))
    else:
        return np.sum(np.abs(m1-m2))
        



# Will output a "melted" DF
# SAMPLE1, SAMPLE2, Distance
def all_vs_all(vcf_matrix_dict):

    keys = list(vcf_matrix_dict.keys())

    n_solutions = (len(keys)**2 + len(keys))/2

    out_df_list = []
    
    n_sites = len(vcf_matrix_dict[keys[0]])


    for i, key_start in enumerate(keys[:-1]):
        distances = [(key_start, key_end, dist(vcf_matrix_dict[key_start], vcf_matrix_dict[key_end]), i+j+i,) \
                     for j,key_end in enumerate(keys[i+1:])]

        distances_rev = [(key_end, key_start, score, comp) for key_start, key_end, score, comp in distances]

        out_df_list.append(pd.DataFrame.from_records(distances + distances_rev,columns=["Sample1","Sample2","L1norm","comp"]))

        n_solutions -= 2
    
    out_df = pd.concat(out_df_list)
    
    out_df["n_sites"] = n_sites

    return out_df
            




    




















if __name__ == "__main__":
    parser = argparse.ArgumentParser()


    parser.add_argument("-p", dest="vcfpath", default="",type=str, help="Folder containing VCF files")
    parser.add_argument("-v", dest="vcffiles", default="", type=str, help="Single or list of vcf files (separated by a comma)")
    parser.add_argument("-o", dest="output", type=str, help="Folder containing VCF files")

    args = parser.parse_args()

    if args.vcffiles != "":
        filelist = args.vcffiles.split(",")
    else:
        filelist = None

    if args.vcfpath != "":
        vcf_path = args.vcfpath
    else:
        vcf_path = None


    vcf_df = read_vcfs(filelist=filelist, path=vcf_path)

    vcf_m = vcf_matrices(vcf_df)

    vcf_m_filt = filter_vcf(vcf_m)

    all_vs_all = all_vs_all(vcf_m)


    all_vs_all.to_csv(args.output)









