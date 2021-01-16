# SNPDistance

Create L1-norm distances between BAMs/VCFs from the same reference.  It specifically uses `lofreq` to get a quantitative representation of all the alleles at a site in a sample. While we do output position depth to discover sites with low/no coverage, it does not utilize it yet, i.e. when calculating the L1-norm distance between two samples where either have no information, it does not take this into account.  We should fix this in the near future.  There are also some internal filtering parameters that are fixed, but arguments to the `snpdistance` script should be added soon. We may include a "reference" L1-norm distance. . .for calibration?


# Requirements


General:

- lofreq
- samtools

Python:

- pandas
- numpy
- vcfpy



# Workflow

1. Input BAMs are run through lofreq's viterbi algorithm to realign reads
2. SNPs are called with `lofreq`
3. The VCFs are processed with `scripts/snpdistance.py`
4. Output is:
- output/vcf: contains all the VCFs and positional depths
- output/l1norm.csv: A csv containing the all-vs-all (redundant) distance comparisons. 




Example of the output:
```
   Sample1            Sample2            L1norm               comp  n_sites
1  K002191            KPCOVID-0299       40.44569             1     175
1  KPCOVID-0299       K002191            40.44569             1     175
```

- `Sample1/Sample2` are the samples compared with each other
- `L1norm` is the L1norm distance
- `comp` is a field used if you would like to remove redundancy (only A v B and not B v A in the data frame)
- `n_sites` shows how many sites (positions) were used to calculate the L1norm score


# Running

General:

`snakemake -j {threads} -p -k --config inputdir={inputdir} reference=assets/references/2697049.fasta # This is the reference of the BAMs`

Where `threads` is the maximum threads needed, `{inputdir}` is the path to the BAMs. The reference can be modified as needed (make sure the FASTA header's accession number matches the one in the BAMs. Working on a way to automate this.


