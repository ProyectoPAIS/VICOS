# VICOS

## Software

The whole pipeline runs using Docker (https://docs.docker.com/engine/install/), there is no need to install nothing else.
Then, the *vicos* script/wrapper can be used to run the pipeline more conveniently
```shell script
# it must be donwloaded in the folder where all the data is stored
wget https://raw.githubusercontent.com/ProyectoPAIS/VICOS/main/vicos -O vicos
chmod +x vicos
```

## Whole Pipeline

To run the pipeline from begining to end:
```shell script
./vicos minority_analysis.py complete -i ./bams_folder -o ./output
```
This assumes that: 
- All bam files are in bams_folder and that folder is INSIDE the current directory
- There is enough space to process everything (usually much less than the space occupied by the bam files)
- Writing permissions over the current directory

## General usage

```shell script
# Every command has this structure:
./vicos minority_analysis.py [command] [options]
# using -h as an argument you can check the rest of the options
./vicos minority_analysis.py [command] -h
```

All scripts bellow assumes that the used data is INSIDE the current directory. If that is not the case, you have
to add a volume mapping. For example if you have the bam files in /something/bams and you are in /some/other/dir
when you execute the docker command you have to add "-v /something/bams:/something/bams" like this:

```shell script
./vicos minority_analysis.py [command] [options]
```

## Example step by step 
```shell script

# download and process required files (reference sequence and annotation)
./vicos minority_analysis.py download
# by default creates a folder named 'data' where the reference is stored. 

# build the vcfs from BAM files. Each bam input should have GATK mandatory fields (RG ID SM PL LB ) -> http://broadinstitute.github.io/picard/command-line-overview.html#AddOrReplaceReadGroups
./vicos minority_analysis.py bam2vcf --bams_folder ./coinfection -o ./vcfs
# this creates one vcf per sample and stores information for each position

# merge variants: combines all samples and variant positions from all vcf files
mkdir results
./vicos minority_analysis.py merge_vcfs --vcfs_dir ./vcfs -o ./results/combined.vcf
./vicos vcffixer.py ./results/combined.vcf -o ./results/combined_fixed.vcf


# filter and process low freq vars
# by default
# min_allele_depth=10   Minimun allele read depth
# min_coverage=0.8      max %N to discard a position 
# min_freq=0.2          minimun minority variant frequency
./vicos minority_analysis.py iSNVs --vcf ./results/combined_fixed.vcf --out ./results/variants.json
# the output, in this case ./results/variants.json, has all the information of the FILTERED low frequency variants, 
# and some distribution stats

# comparative analysis: coinfection candidates are detected by analyzing low_frequency variant counts in each sample.
# by default deviation_lowfreq(default 2) is used (samples with more than mean + 2*STD low_frequency variants are classified as coinfection candidates)  
# This assumes that most samples will NOT be a coinfection. If that is not the case, min_lowfreq can be used, where you 
# have to set the number of low_frequency mutations a sample has to have to be a coinfection candidate
./vicos  minority_analysis.py candidates --data ./results/variants.json --out_dir ./results/report
# Output files are:
# - candidates_freqs.png : low frequency variant counts per sample
# - min_variants_count_per_sample.png : low frequency variant counts per sample
# - min_variants_count_per_sample_with_0.png : low frequency variant counts per sample
# - candidates_summary.csv:
#   - pos
#   - sample_consesnsus
#   - dataset_consensus
#   - exlusive_minoritary
#   - exlusive_majoritary
#   - depth
#   - min_mut
#   - min_freq
# - variants_list.csv:
#   - pos
#   - sample_consesnsus
#   - dataset_consensus
#   - exlusive_minoritary
#   - exlusive_majoritary
#   - depth
#   - min_mut
#   - min_freq
# - one CSV file per candidate with the following fields:
#   - pos
#   - sample_consesnsus
#   - dataset_consensus
#   - exlusive_minoritary
#   - exlusive_majoritary
#   - depth
#   - min_mut
#   - min_freq
```

## Programs used by the docker image

* GATK 4.1.8.0
* BCFtools
* Samtools

