"""

SNP distances 


1. Iterate through BAMs
2. Do the Viterbi thing
3. Lofreq calls
4. All vs All comparison

Tailored for SARS-CoV-2 right now

"""

import pandas as pd 
import os
from glob import glob

def get_input(inputdir):

	bam_files = glob("%s/*.bam"%inputdir)
	basenames = list(map(os.path.basename, bam_files))

	bam_files_ref = pd.Series(bam_files, index=basenames)
	delim_df = pd.DataFrame.from_records([(basename, basename.count("_"),) for basename in basenames], columns=["basename","delim_count"])

	out_df_list = []

	for c, df in delim_df.groupby("delim_count"):
		if c == 0:
			out_df_list.append(pd.DataFrame(dict(samplename=df.basename.apply(lambda x:x.replace(".bam","")).values,
							     filename=bam_files_ref[df.basename.values].values)))
			continue

		for i in range(c):
			s = pd.Series(map(lambda x:"_".join(x.split("_")[:(i+1)]), df.basename))
			if len(s.drop_duplicates()) == len(s):
				out_df_list.append(pd.DataFrame(dict(samplename=s.values, filename=bam_files_ref[df.basename.values].values)))
				break





	return pd.concat(out_df_list).set_index("samplename",drop=False)



sample_df = get_input(config["inputdir"])

sample_df.to_csv("./sampledf.csv")




rule all:
	input:
		l1norm="output/l1norm.csv"



rule viterbi:
	input:
		BAM=lambda wildcards: sample_df.loc[wildcards.sample].filename
	output:
		BAM=temp("output/bam/{sample}.viterbi.bam"),
		BAI=temp("output/bam/{sample}.viterbi.bam.bai"),
	
	params:
		ref=config["reference"]
	threads: 2

	shell:
		"lofreq viterbi -f {params.ref} {input.BAM} | samtools sort -@ {threads} -o {output.BAM} && samtools index {output.BAM}"

rule depth:
	input:
		BAM="output/bam/{sample}.viterbi.bam",
		BAI="output/bam/{sample}.viterbi.bam.bai"


	output:
		"output/vcf/{sample}.depth.gz"
	shell:
		"samtools depth {input.BAM} | gzip -c -9 - > {output}"



rule vcf:
	input:
		BAM="output/bam/{sample}.viterbi.bam",
		BAI="output/bam/{sample}.viterbi.bam.bai",
	output:
		"output/vcf/{sample}.vcf"
	params:
		ref=config["reference"]
	threads: 4
	shell:
		"lofreq call-parallel --pp-threads {threads} -f {params.ref} {input} -o {output}"



rule lnorm:
	input:
		expand("output/vcf/{sample}.{what}", sample=sample_df["samplename"],what=["vcf","depth.gz"])
	output:
		"output/l1norm.csv"
	params:
		path="output/vcf"
	shell:
		"python3 scripts/snpdistance.py -p {params.path} -o {output}"







