# SARS-CoV-2 co-infections

## Software
The whole pipeline is runned using Docker and a basic Python instalation, there is no need to install nothing else.

## Steps

```shell script
# download required files
python minority_analysis.py download

# get combined vcf
python minority_analysis.py bam2vcf --bams_folder ./coinfection > combined.vcf
# created the file  

# process minority variants
python minority_analysis.py bam2vcf --vcf combined.vcf


```


 
## Example







## Anex
### How to get the mapping
There are many ways to build the bam file


### Used programs

* GATK 4.1.8.0
* BCFtools
* Mafft
* ...
