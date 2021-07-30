import os
import sys

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Process minory variants')

    subparsers = parser.add_subparsers(help='commands', description='valid subcommands', dest='command', required=True)

    cmd = subparsers.add_parser('bam2vcf', help='variant calling pipeline')
    cmd.add_argument('-ref', '--reference', required=True, help='fasta file. Can be gziped')
    cmd.add_argument('-b', '--bam', required=True)

    cmd = subparsers.add_parser('bam2vcf', help='variant calling pipeline')
    cmd.add_argument('-ref', '--reference', required=True, help='fasta file. Can be gziped')
    cmd.add_argument('-b', '--bam', required=True)

    cmd = subparsers.add_parser('vcf2minconsensus', help='genotyping for multi sample')
    cmd.add_argument('-ref', '--reference', required=True, help='fasta file. Can be gziped')
    cmd.add_argument('--vcf', required=True)

    args = parser.parse_args()
    if args.command == 'vc':
        pass
    if args.command == 'gt':
        pass