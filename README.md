# SARS-CoV-2 co-infections

## Software

```shell script
git clone https://github.com/ProyectoPAIS/CoV-2-co-infections.git
cd CoV-2-co-infections
```

The whole pipeline runs using Docker and a basic Python installation, there is no need to install nothing else.

## Steps

```shell script
# download required files
docker run -u $(id -u ${USER}):$(id -g ${USER}) --rm -w /out -v $PWD:/out sndg/covminanalysis minority_analysis.py download


# get vcf files from BAM files
docker run -u $(id -u ${USER}):$(id -g ${USER}) --rm -w /out -v $PWD:/out sndg/covminanalysis minority_analysis.py bam2vcf --bams_folder ./coinfection -o ./vcfs
 
# merge variants
mkdir results
docker run -u $(id -u ${USER}):$(id -g ${USER}) --rm -w /out -v $PWD:/out sndg/covminanalysis minority_analysis.py merge_vcfs --vcfs_dir ./vcfs -o ./results/combined.vcf

# process low freq vars
docker run -u $(id -u ${USER}):$(id -g ${USER}) --rm -w /out -v $PWD:/out sndg/covminanalysis minority_analysis.py minor_freq_vars --vcf ./results/combined.vcf --out ./results/variants.json 

# comparative analysis
docker run -u $(id -u ${USER}):$(id -g ${USER}) --rm -w /out -v $PWD:/out sndg/covminanalysis minority_analysis.py report --data ./results/variants.json --out_dir ./results/
```
 
# Used Programs

* GATK 4.1.8.0
* BCFtools
* Mafft
* ...
