# SARS-CoV-2 co-infections

## Software

```shell script
conda create --name covid19 python=3.8
conda activate covid19
conda install -c conda-forge biopython pandas matplotlib numpy seaborn
```

The whole pipeline runs using Docker and a basic Python installation, there is no need to install nothing else.

## Steps

```shell script
# download required files
python minority_analysis.py download

# get vcf files from BAM files
python minority_analysis.py bam2vcf --bams_folder ./coinfection -o ./vcfs
 
# merge variants
python minority_analysis.py merge_vcfs --vcfs_dir ./vcfs -o ./results/combined.vcf

# process minority variants
mkdir results
python minority_analysis.py vcf2minconsensus --vcf ./results/combined.vcf --out ./results/variants.json 

# comparative analysis
python minority_analysis.py report --data ./results/variants.json --out_dir ./results/
```
 
# Used Programs

* GATK 4.1.8.0
* BCFtools
* Mafft
* ...
