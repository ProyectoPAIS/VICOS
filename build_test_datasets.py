#!/usr/bin/env python3

import subprocess as sp
import os


# /mnt/data2/projects/covid

def e(cmd):
    if os.environ.get("verbose"):
        print(cmd)
    sp.run(cmd, shell=True)


if __name__ == '__main__':
    import argparse
    import random
    from tqdm import tqdm

    parser = argparse.ArgumentParser(description='Creates in silico mixtures of samples')
    parser.add_argument('-v', '--verbose', action='store_true')

    parser.add_argument('-f1-1', '--fastq1_1', required=True)
    parser.add_argument('-f1-2', '--fastq1_2', required=True)
    parser.add_argument('-f2-1', '--fastq2_1', required=True)
    parser.add_argument('-f2-2', '--fastq2_2', required=True)

    parser.add_argument('-o', '--output_folder', default="./")
    parser.add_argument('-p', '--proportions', nargs="+",
                        default=["10", "20", "30", "40", "50", "60", "70", "80", "90"])

    args = parser.parse_args()

    if args.verbose:
        os.environ["verbose"] = "y"

    assert os.path.exists(args.fastq1_1), f"'{args.fastq1_1}' does not exits"
    assert os.path.exists(args.fastq1_2), f"'{args.fastq1_2}' does not exits"
    assert os.path.exists(args.fastq2_1), f"'{args.fastq2_1}' does not exits"
    assert os.path.exists(args.fastq2_2), f"'{args.fastq2_2}' does not exits"

    if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder)

    assert os.path.exists(args.output_folder), f"'{args.output_folder}' could not be created"

    for proportion in tqdm(args.proportions):
        p1 = float(proportion) / 100
        p2 = 1 - p1
        seed = random.randint(1,9999)
        e(f"seqtk sample -s{seed} {args.fastq1_1} {p1} > {args.output_folder}/s{proportion}.R1.fq.gz ")
        e(f"seqtk sample -s{seed} {args.fastq2_1} {p1} >> {args.output_folder}/s{proportion}.R1.fq.gz ")

        e(f"seqtk sample -s{seed} {args.fastq1_2} {p2} > {args.output_folder}/s{proportion}.R2.fq.gz ")
        e(f"seqtk sample -s{seed} {args.fastq2_2} {p2} >> {args.output_folder}/s{proportion}.R2.fq.gz ")
