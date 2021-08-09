import os
import sys

if __name__ == '__main__':
    import argparse
    import subprocess as sp
    from glob import glob
    import gzip
    from collections import defaultdict

    parser = argparse.ArgumentParser(description='Process minory variants')

    subparsers = parser.add_subparsers(help='commands', description='valid subcommands', dest='command', required=True)

    cmd = subparsers.add_parser('download', help='Downloads reference and sample data')
    # cmd.add_argument('-o', '--output_folder', default="./data")

    cmd = subparsers.add_parser('bam2vcf', help='variant calling pipeline')
    cmd.add_argument('-bams', '--bamfiles_dir', required=True)
    cmd.add_argument('-ref', '--reference', default="data/MN996528.fna", help='fasta file. Can be gziped')
    cmd.add_argument('-o', '--output', default="results/")



    cmd = subparsers.add_parser('vcf2minconsensus', help='genotyping for multi sample')
    cmd.add_argument('-mc', '--min_coverage', default=0.8,type=float, help='max %N to discard a position. Between 0.0-1.0 Default 0.8')
    cmd.add_argument('-mf', '--min_freq', default=0.2,type=float, help='minimun minority variant frequence. Between 0.0-1.0 Default 0.8')
    cmd.add_argument('-ocut', '--outlayers_cutoff', default=0.2,type=float, help='minimum minority variant frequency. Between 0.0-1.0 Default 0.8')

    cmd.add_argument('--vcf', required=True,help="Multi Sample VCF. GT and AD fields are mandatory")

    args = parser.parse_args()
    if args.command == 'download':
        outfolder = "data" # args.output_folder
        if not os.path.exists(outfolder):
            os.makedirs(outfolder)
        if not os.path.exists(outfolder):
            sys.stderr.write(f"'{outfolder}' could not be created")
            sys.exit(1)
        cmd = f'wget -O {outfolder}/MN996528.fna "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=1802633808&conwithfeat=on&withparts=on&hide-cdd=on"'
        sp.run(cmd)
        cmd = f"docker run -v $PWD:/out -w /out staphb/samtools samtools dict -o {outfolder}/MN996528.dict {outfolder}/MN996528.fna  "
        # makeblastdb -dbtype nucl -in /ref/MN996528.fna && bwa index MN996528.fna && samtools faidx MN996528.fna
        sp.run(cmd)

    if args.command == 'bam2vcf':
        outfolder = args.output
        if not os.path.exists(outfolder):
            os.makedirs(outfolder)
        if not os.path.exists(outfolder):
            sys.stderr.write(f"'{outfolder}' could not be created")
            sys.exit(1)

        for bam_file in glob(args.bamfiles_dir + "/*.bam"):
            sample = bam_file.split("/")[-1].split(".bam")[0]
            f"""docker run -v $PWD:/out -w /out broadinstitute/gatk:4.1.0.0 \
             java -jar /gatk/gatk-package-4.1.0.0-local.jar  HaplotypeCaller -ERC GVCF -R {args.reference_fasta} \
	            -ploidy 2 -I {bam_file} --output-mode EMIT_ALL_SITES -O {args.output}/{sample}.g.vcf.gz"""

    if args.command == 'vcf2minconsensus':
        """MN996528.1	1879	.	A	G	14552.79	.	AC=2;AF=5.556e-03;AN=360;BaseQRankSum=2.16;DP=113297;ExcessHet=0.0061;FS=0.838;InbreedingCoeff=0.8880;MLEAC=2;MLEAF=5.556e-03;MQ=59.99;MQRankSum=0.00;QD=28.76;ReadPosRankSum=1.03;SOR=0.597	GT:AD:DP:GQ:PGT:PID:PL:PS	0/0:717,0:717:99:.:.:0,120,1800	0/0:166,0:166:99:.:.:0,120,1800	0/0:395,0:395:99:.:.:0,120,1800	0/0:354,0:354:99:.:.:0,120,1800	0/0:464,0:464:99:.:.:0,120,1800	0/0:363,0:363:99:.:.:0,120,1800	0/0:464,0:464:99:.:.:0,120,1800	0/0:396,0:396:99:.:.:0,120,1800	0/0:288,0:288:99:.:.:0,120,1800	0/0:371,0:371:99:.:.:0,120,1800	0/0:347,0:347:99:.:.:0,120,1800	0/0:422,0:422:99:.:.:0,120,1800	0/0:517,0:517:99:.:.:0,120,1800	0/0:392,0:392:99:.:.:0,120,1800	0/0:392,0:392:99:.:.:0,120,1800	0/0:465,0:465:99:.:.:0,120,1800
        """
        if not os.path.exists(args.vcf):
            sys.stderr.write(f"'{args.vcf}' does not exists\n")
            sys.exit(1)
        h = gzip.open(args.vcf,"rt") if args.vcf.endswith(".gz") else open(args.vcf)
        try:
            samples = None
            variants = defaultdict(list)
            for line in h:
                if line.startswith("#CHROM"):
                    samples = line.split()[9:]
                elif not line.startswith("#"):
                    vec = line.split()
                    pos = int(vec [1])
                    ref = vec [3]
                    alts = {i:x for i,x in enumerate(vec[4].split(","))}
                    format_f = vec[8]#GT:AD:DP:GQ:PGT:PID:PL:PS
                    gt_index = [i for i,x in enumerate(format_f.split(":")) if x == "GT"]
                    if not gt_index:
                        sys.stderr.write("not all pos/samples have a GT field")
                        sys.exit(2)
                    gt_index = gt_index[0]
                    ad_index = [i for i,x in enumerate(format_f.split(":")) if x == "AD"]
                    if not ad_index:
                        sys.stderr.write("not all pos/samples have a AD field")
                        sys.exit(2)
                    ad_index = ad_index[0]
                    for idx,sample in enumerate(samples):
                        gt = vec[9+idx].split(":")[gt_index]

                        ad = vec[9+idx].split(":")[ad_index]
                        if len(set(gt.split("/"))) > 1:
                            variants[pos].append( [sample,gt,ad] )
                            break
        finally:
            h.close()
        variants = dict(variants)
        