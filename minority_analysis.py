#!/usr/bin/env python3

import os
import sys
import json
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import Bio.SeqIO as bpio
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from tqdm import tqdm
from itertools import groupby
from Bio.SeqUtils import seq1


def e(cmd):
    if os.environ.get("verbose"):
        print(cmd)
    sp.run(cmd, shell=True)


def merge_vcfs(args):
    outfolder = os.path.dirname(args.output)
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    if not os.path.exists(outfolder):
        sys.stderr.write(f"'{outfolder}' could not be created")
        sys.exit(1)
    if not os.path.exists(args.reference):
        sys.stderr.write(f"'{args.reference}' does not exists")
        sys.exit(1)
    vcf_files = []
    for x in glob(args.vcfs_dir + "/*vcf*"):
        if x.endswith(".vcf") or x.endswith(".vcf.gz") or x.endswith(".gvcf") or x.endswith(".gvcf.gz"):
            vcf_files.append(x)
    if not vcf_files:
        raise FileNotFoundError(f'no .vcf or .vcf.gz files where found at {args.vcfs_dir}')
    vcfs = " ".join([f"--variant {x}" for x in vcf_files])
    # docker run -u $(id -u ${{USER}}):$(id -g ${{USER}})  -v $PWD:/out -w /out broadinstitute/gatk:4.2.2.0 \
    cmdx = f"""/gatk/gatk CombineGVCFs -R {args.reference} {vcfs} -O {args.output}.raw.gz"""
    e(cmdx)
    # docker run -u $(id -u ${{USER}}):$(id -g ${{USER}})  -v $PWD:/out -w /out broadinstitute/gatk:4.2.2.0 \
    cmdx = f"""/gatk/gatk GenotypeGVCFs \
                    -R "{args.reference}" -ploidy 2 \
                    -V "{args.output}.raw.gz" \
                    -O "{args.output}.unann" 
        """
    e(cmdx)
    cmdx = f"""java -jar /opt/snpEff/snpEff.jar ann covid19 "{args.output}.unann"  > "{args.output}" """
    e(cmdx)


def download(args):
    outfolder = args.output_folder
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    if not os.path.exists(outfolder):
        sys.stderr.write(f"'{outfolder}' could not be created")
        sys.exit(1)
    cmdx = f'wget -O {outfolder}/MN996528.fna "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=1802633808&conwithfeat=on&withparts=on&hide-cdd=on"'
    sp.run(cmdx, shell=True)
    # docker run -u $(id -u ${{USER}}):$(id -g ${{USER}}) -v $DIND:/out -w /out staphb/samtools \
    cmdx = f"""samtools dict -o {outfolder}/MN996528.dict {outfolder}/MN996528.fna"""
    # makeblastdb -dbtype nucl -in /ref/MN996528.fna && bwa index MN996528.fna && samtools faidx MN996528.fna
    e(cmdx)
    # docker run -u $(id -u ${{USER}}):$(id -g ${{USER}}) -v $DIND:/out -w /out staphb/samtools \
    cmdx = f"""samtools faidx {outfolder}/MN996528.fna"""
    e(cmdx)


def variant_filter(args, lineage_data={}):
    """MN996528.1	1879	.	A	G	14552.79	.	AC=2;AF=5.556e-03;AN=360;BaseQRankSum=2.16;DP=113297;ExcessHet=0.0061;FS=0.838;InbreedingCoeff=0.8880;MLEAC=2;MLEAF=5.556e-03;MQ=59.99;MQRankSum=0.00;QD=28.76;ReadPosRankSum=1.03;SOR=0.597	GT:AD:DP:GQ:PGT:PID:PL:PS	0/0:717,0:717:99:.:.:0,120,1800	0/0:166,0:166:99:.:.:0,120,1800	0/0:395,0:395:99:.:.:0,120,1800	0/0:354,0:354:99:.:.:0,120,1800	0/0:464,0:464:99:.:.:0,120,1800	0/0:363,0:363:99:.:.:0,120,1800	0/0:464,0:464:99:.:.:0,120,1800	0/0:396,0:396:99:.:.:0,120,1800	0/0:288,0:288:99:.:.:0,120,1800	0/0:371,0:371:99:.:.:0,120,1800	0/0:347,0:347:99:.:.:0,120,1800	0/0:422,0:422:99:.:.:0,120,1800	0/0:517,0:517:99:.:.:0,120,1800	0/0:392,0:392:99:.:.:0,120,1800	0/0:392,0:392:99:.:.:0,120,1800	0/0:465,0:465:99:.:.:0,120,1800
        """
    if not os.path.exists(args.vcf):
        sys.stderr.write(f"'{args.vcf}' does not exists\n")
        sys.exit(1)

    print(f"Running lowfreq variant detection with:")
    print(f'- Minimun allele read depth: {args.min_allele_depth}')
    print(f'- max %N to discard a position: {args.min_coverage}')
    print(f'- minimun minority variant frequency: {args.min_freq}')
    print("----------------")

    lineage_variants_count = defaultdict(lambda: 0)
    for lineajes in lineage_data.values():
        for lin, _ in lineajes:
            lineage_variants_count[lin] += 1

    for k, v in lineage_variants_count.items():
        if v <= 8:
            for var, lin_freqs in lineage_data.items():
                lineage_data[var] = [(x, y) for x, y in lin_freqs if x != k]
    lineage_variants_count = {k: v for k, v in lineage_variants_count.items() if v >= 8}

    h = gzip.open(args.vcf, "rt") if args.vcf.endswith(".gz") else open(args.vcf)
    number_of_variable_sites = 0
    number_of_mutations = 0
    ns_per_sample = defaultdict(lambda: 0)
    sample_lineages = defaultdict(dict)
    try:
        samples = None
        variants = {}
        for line in tqdm(h):
            if line.startswith("#CHROM"):
                samples = line.split()[9:]


            elif not line.startswith("#"):

                # ;ANN=G|missense_variant|MODERATE|S|Gene_21562_25383|transcript|QHR63260.2|protein_coding|1/1|c.3815A>G|p.Tyr1272Cys|3815/3822|3815/3822|1272/1273||

                number_of_variable_sites += 1

                vec = line.split()
                ann = vec[7].split(";ANN=")[1].split(";")[0].split(",")[0].split("|")
                gene = ann[3]
                gene_nt = ann[9][2:]
                gene_aa = ann[10][2:]

                # MN996528.1      3539    .       GCTA    G       2837.64 .       AC=1;AF=4.281e-04;AN=2336;BaseQRankSum=1.57;DP=317701;
                # ExcessHet=3.0103;FS=1.823;InbreedingCoeff=0.5465;MLEAC=1;MLEAF=4.281e-04;MQ=59.97;MQRankSum=0.00;QD=9.21;
                # ReadPosRankSum=-4.450e-01;SOR=0.792;ANN=
                # G|disruptive_inframe_deletion|MODERATE|ORF1a|Gene_265_13467|transcript|QHR63259.1|protein_coding|1/1|c.3278_3280delCTA|p.Thr1093del|3278/13203|3278/13203|1093/4400||WARNING_TRANSCRIPT_NO_STOP_CODON&INFO_REALIGN_3_PRIME

                # MN996528.1      3545    .       AATG    A       1742.27 .       AC=3;AF=1.290e-03;AN=2326;BaseQRankSum=0.888;DP=317477;
                # ExcessHet=3.0159;FS=0.000;InbreedingCoeff=0.5466;MLEAC=3;MLEAF=1.290e-03;MQ=59.97;MQRankSum=0.00;QD=5.66;
                # ReadPosRankSum=-3.580e-01;SOR=0.552;ANN=
                # A|disruptive_inframe_deletion|MODERATE|ORF1a|Gene_265_13467|transcript|QHR63259.1|protein_coding|1/1|c.3281_3283delATG|p.Asn1094_Gly1095delinsArg|3281/13203|3281/13203|1094/4400||WARNING_TRANSCRIPT_NO_STOP_CODON

                if ("del" in gene_aa) and ("disruptive_inframe_deletion" in ann[1]):
                    gene_aa = "".join([x for x in gene_aa.split("del")[0].replace("_", "/") if x.isdigit() or x == "/"])
                    if "/" not in gene_aa:
                        gene_aa = gene_aa + "/" + gene_aa
                    gene_aa = "DEL" + gene_aa
                if ("stop_gained" in ann[1]) or ("missense_variant" in ann[1]):
                    gene_aa = seq1(gene_aa[:3]) + "".join([x for x in gene_aa[3:] if x.isdigit()]) + seq1(
                        "".join([x for x in gene_aa[3:] if not x.isdigit()]))

                lineage_variant_key = gene + ":" + gene_aa
                variant_lineages = []
                if lineage_variant_key in lineage_data:
                    variant_lineages = lineage_data[lineage_variant_key]


                pos = int(vec[1])

                # if pos < 28000:
                #     continue

                ref = vec[3]
                alts = {i: x for i, x in enumerate(vec[4].split(","))}

                gt_options = {k + 1: v for k, v in alts.items()}
                gt_options[0] = ref
                gt_options[99] = "N"
                format_f = vec[8]  # GT:AD:DP:GQ:PGT:PID:PL:PS
                gt_index = [i for i, x in enumerate(format_f.split(":")) if x == "GT"]
                if not gt_index:
                    sys.stderr.write("not all pos/samples have a GT field")
                    sys.exit(2)
                gt_index = gt_index[0]
                ad_index = [i for i, x in enumerate(format_f.split(":")) if x == "AD"]
                if not ad_index:
                    sys.stderr.write("not all pos/samples have a AD field")
                    sys.exit(2)
                dp_index = [i for i, x in enumerate(format_f.split(":")) if x == "DP"]
                if not dp_index:
                    sys.stderr.write("not all pos/samples have a GT field")
                    sys.exit(2)

                ad_index = ad_index[0]
                dp_index = dp_index[0]
                number_of_mutations += len(set(vec[4].split(",")) - set(["*", "N"]))
                low_freq = False
                pos_data = {}
                for idx, sample in enumerate(samples):
                    gt = vec[9 + idx].split(":")[gt_index]
                    if gt == "./.":
                        ns_per_sample[sample] += 1
                    # dp = int(vec[9 + idx].split(":")[dp_index])

                    ad = vec[9 + idx].split(":")[ad_index]
                    # (28280 == pos) and (gt.replace("|","/") in [("0/1")])
                    ads = [int(x) for x in ad.split(",") if x != "."]
                    if ads:
                        min_ad = sorted(ads)[-2] # biggest second
                        gt_vec = [int(x) if x != "." else 99 for x in gt.replace("|", "/").split("/")]

                        pos_data[sample] = [{gt_options[gt_num]:
                                                 ([int(x) for x in ad.split(",")[gt_num]] if gt_num != 99 and (
                                                         min_ad >= args.min_allele_depth) else "?")
                                             for gt_num in gt_vec},
                                            {gt_options[i]: int(ad_num) for i, ad_num in enumerate(ad.split(","))}, ref,
                                            (gene, gene_nt, gene_aa)
                                            ]

                        if ads and min_ad >= args.min_allele_depth and len(set(gt.replace("|", "/").split("/"))) > 1:
                            low_freq = True

                        if variant_lineages:
                            sample_lineages[sample][lineage_variant_key] = variant_lineages

                if low_freq:
                    variants[pos] = pos_data

    finally:
        h.close()

    variants = dict(variants)
    ns_per_sample = dict(ns_per_sample)
    sample_lineages2 = {}
    for sample, lineage_sample_data in sample_lineages.items():
        lineages_count = defaultdict(lambda: 0)
        for key, variant_lineages in lineage_sample_data.items():
            for lineage, freq in variant_lineages:
                lineages_count[lineage] += 1
        lineages_count = dict(lineages_count)

        sample_lineages2[sample] = sorted(
            [(y[0], y[1], round(y[1] * 1.0 / lineage_variants_count[y[0]], 2)) for y in lineages_count.items()
             ], key=lambda x: x[2], reverse=True)
        sample_lineages2[sample] = [y for y in sample_lineages2[sample] if
                                    round(y[1] * 1.0 / lineage_variants_count[y[0]], 2) > 0.8]

    badqualitysamples = {s: v for s, v in ns_per_sample.items() if v > args.badq_strain_ns_threshold}
    print(f"Number of samples: {len(samples)}")
    print(f'Number of bad quality samples: {len(badqualitysamples)}')
    print(f'number of variable sites: {number_of_variable_sites}')
    print(f'number of mutations: {number_of_mutations}')

    excluded_positions = []
    low_freq_freq = []
    high_freq_freq = []
    # low_freq_pos = []
    low_freq_muts = []
    entries_data = {}

    for pos, samples_variants in variants.items():
        ns = 0
        for sample, (gts, ads, ref, ann) in samples_variants.items():
            is_n = (1 if "N" in gts else 0)
            ns += is_n

        if (1 - (1.0 * ns / len(samples))) < args.min_coverage:
            excluded_positions.append(pos)
            #     if not is_n:
            #         low_freq_muts.append(low_freq_mut)
            #         low_freq_freq += min_freqs
            #         high_freq_freq += consensus_freqs
            # else:

    print(f'excluded positions( Ns count greater than threashold):{len(excluded_positions)}')
    discarded_low_depth = defaultdict(list)
    variant_samples = defaultdict(list)
    for pos, samples_variants in variants.items():
        variant_samples_pos = defaultdict(list)
        if pos not in excluded_positions:
            min_freqs = []
            consensus_freqs = []
            entries = {}
            valid_min_variant = False
            for sample, (gts, ads, ref, (gene, gene_nt, gene_aa)) in samples_variants.items():

                if len(gts) > 1:
                    depth = sum(ads.values())
                    freqs = [(k, 1 * v / depth) for k, v in sorted(ads.items(), key=lambda x: x[1])]
                    min_variant = freqs[-2]
                    dp = sum(ads.values())
                    if min_variant[1] >= args.min_freq:
                        if (dp >= args.min_allele_depth):
                            min_freqs.append(min_variant[1])
                            consensus_variant = freqs[-1]
                            consensus_freqs.append(min_variant[1])
                            low_freq_mut = [pos, min_freqs]
                            variant_samples_pos[f'{pos}_{ref}_{min_variant[0]}'].append(sample)
                            valid_min_variant = True
                            low_freq_freq.append(freqs[-2][1])
                            high_freq_freq.append(freqs[-1][1])
                        else:
                            discarded_low_depth[sample].append([pos, ads])
                    else:
                        consensus_variant = list(gts.items())[0][0]
                        min_variant = [""]
                        high_freq_freq.append(freqs[-1][1])
                else:
                    consensus_variant = list(gts.items())[0][0]
                    min_variant = [""]
                lineage_key = gene + ":" + gene_aa

                if (sample in sample_lineages) and (lineage_key in sample_lineages[sample]):
                    lineages = [(lineage, lfreq) for lineage, lfreq in sample_lineages[sample][lineage_key]
                                if lineage in [x[0] for x in sample_lineages2[sample]]]
                else:
                    lineages = []

                entries[sample] = [consensus_variant, min_variant, gts, ads, ref, (gene, gene_nt, gene_aa), lineages]
            if valid_min_variant:
                entries_data[pos] = entries
                for k, v in variant_samples_pos.items():
                    variant_samples[k] = v

    print(f'low freq positions:{len(entries_data)}')
    print(f'low freq mutations:{len(variant_samples)}')

    variant_samples = dict(variant_samples)
    entries_data = dict(entries_data)

    with open(f"{args.out}", "w") as h:
        data = {"variant_samples": variant_samples, "entries_data": entries_data,
                "discarded_low_depth": discarded_low_depth,
                "excluded_positions": excluded_positions, "low_freq_freq": low_freq_freq,
                "high_freq_freq": high_freq_freq, "badquality_samples": badqualitysamples,
                "sample_lineages": sample_lineages2}
        json.dump(data, h)


def aln(h, output, refseq=None, included_samples=None):
    # if hasattr(vcf_file, "read"):
    #     h = vcf_file
    # else:
    #     h = open(vcf_file)
    try:
        base_idx = 0
        for line in h:
            if line.startswith("#"):
                if line.startswith("#CHROM"):
                    samples = [x.strip() for x in line.split()[9:]]

                    seqmap = {s: "" for s in samples}
                continue
            break

        for line in h:
            _, pos, _, ref, alts = line.split()[:5]
            pos = int(pos)
            alts = [ref] + alts.split(",")
            gts = [x[0] for x in line.split()[9:]]
            gts = ["N" if x[0] == "." else alts[int(x[0])] for x in gts]
            pos_size = max([len(x) for x in alts])
            for i, s in enumerate(samples):
                if not included_samples or s in included_samples:
                    subseq = refseq[base_idx:pos] + gts[i].ljust(pos_size, "-")
                    seqmap[s] += subseq

            sizes = {}
            for s in samples:
                if not included_samples or s in included_samples:
                    sizes[s] = len(seqmap[s])
            assert len(set(sizes.values())) == 1, [base_idx, set(sizes.values()),
                                                   json.dumps({k: [x[0] for x in v] for k, v in
                                                               groupby(sorted(sizes.items(), key=lambda x: x[1]),
                                                                       lambda x: x[1])})]

            base_idx = pos + len(ref)
    finally:
        h.close()

    for s in samples:
        if not samples or s in included_samples:
            seqmap[s] += refseq[base_idx:]

    if hasattr(output, "write"):
        h = output
    else:
        h = open(output, "w")
    try:
        for k, v in seqmap.items():
            bpio.write(SeqRecord(id=k, name="", description="", seq=Seq(v)), h, "fasta")
    finally:
        h.close()


def comparative_analysis(json_file, output_dir, deviation_lowfreq=1, min_lowfreq=None, min_depth=10):
    assert os.path.exists(json_file), f'"{json_file}" does not exists'
    with open(json_file) as h:
        data = json.load(h)

    min_per_sample = defaultdict(list)
    pos_data = {}

    variant_samples = defaultdict(dict)
    min_variant_samples = defaultdict(list)
    ann_variants = {}

    for pos, sample_data in data["entries_data"].items():
        pos_data[pos] = {"consensus": defaultdict(list), "mins": defaultdict(list)}

        for sample, (
                consensus_variant, min_variant, gts, ads, ref, (gene, gene_nt, gene_aa),
                lineages) in sample_data.items():
            min_var_freq = 0
            for allele in ads.keys():
                if (allele in gts) and ads[allele]:
                    variant_samples[f'{pos}_{ref}_{allele}'][sample] = len(gts)

            if min_variant and min_variant[0] and (ads[min_variant[0]] >= min_depth):
                min_var_mut = f'{pos}_{ref}_{min_variant[0]}'
                min_variant_samples[min_var_mut].append(sample)
                ann_variants[min_var_mut] = (gene + ":" + gene_aa) if ref != min_variant[0] else ""

                min_var_freq = min_variant[1]
                min_per_sample[sample].append(min_var_mut)

                pos_data[pos]["mins"][min_var_mut].append(min_var_freq)

            pos_data[pos]["consensus"][consensus_variant[0]].append(1 - min_var_freq)
    min_per_sample = dict(min_per_sample)

    df = pd.DataFrame([{"sample": k, "count": len(v)} for k, v in min_per_sample.items()])

    if min_lowfreq:
        sys.stderr.write(f"using min_lowfreq method to select candidates. cutoff:  {min_lowfreq} \n")
        cutoff = min_lowfreq
    else:

        median = np.mean(df["count"])
        deviation = np.std(df["count"])
        cutoff = median + deviation_lowfreq * deviation
        sys.stderr.write(
            f"deviation_lowfreq method to select candidates: mean {median:.2f} deviation {deviation:.2f} cutoff {cutoff:.2f}\n")

    candidates = list(df[df["count"] > cutoff]["sample"])

    ax = plt.subplot()
    ax.axvline(cutoff, 0, max(df["count"]), color="red")
    plt.xlabel("Low freq variants count")
    plt.ylabel("Samples count")

    ax.hist(data=df, x="count", bins=30)
    ax.plot(df["count"], [0.01] * len(df["count"]), '|', color='k')
    plt.savefig(f'{output_dir}/distplot.png')
    plt.close()

    dfs = {}

    for c in candidates:
        cdata = []
        for min_var in sorted(min_per_sample[c], key=lambda x: int(x.split("_")[0])):
            pos = min_var.split("_")[0]
            aln_consensus = defaultdict(lambda: 0)

            pos_muts = []
            sample_ads = {}
            for sample, (consensus_variant, min_variant, gts, ads, ref, (gene, gene_nt, gene_aa), lineages
                         ) in data["entries_data"][pos].items():

                if sample == c:
                    consensus = consensus_variant[0]
                    min_freq = min_variant[1]
                    min_mut = min_variant[0]
                    sample_ads = ads
                else:
                    pos_muts.append(consensus_variant[0])
                    if min_variant and min_variant[0]:
                        pos_muts.append(min_variant[0])
                aln_consensus[consensus_variant[0]] += 1

            exclusive_minority = bool(set([min_var.split("_")[1]]) - set(pos_muts))
            exclusive_consensus = bool(set([consensus]) - set(pos_muts))
            if "n" in aln_consensus:
                del aln_consensus["N"]

            aln_consensus_str = " ".join([f"{k}:{v}" for k, v in aln_consensus.items()])

            r = {"pos": int(pos), "ref": ref, "sample_consensus": consensus,
                 "exclusive_consensus": exclusive_consensus,
                 "exclusive_min": exclusive_minority,
                 "depth": sum([int(x) for x in sample_ads.values()]), "depth_min": sample_ads[min_mut],
                 "allele_min": min_mut, "freq_min": np.round(min_freq, 2),
                 "dataset_consensus": aln_consensus_str,
                 "gene": gene if gene_aa else "", "gene_nt": gene_nt if gene_aa else "", "gene_aa": gene_aa,
                 "lineages": " ".join([x[0] + "|" + str(round(x[1], 2)) for x in lineages])
                 }
            if sample_ads and (sample_ads[min_mut] >= min_depth):
                cdata.append(r)

        dfs[c] = pd.DataFrame(cdata,
                              columns=["pos", "ref", "gene", "gene_nt", "gene_aa", "sample_consensus",
                                       "dataset_consensus",
                                       "exclusive_min", "exclusive_consensus",
                                       "depth", "allele_min", "depth_min", "freq_min",
                                       "lineages"])

    # variant_samples[f'{pos}_{allele}'].append(sample)
    with open(f'{output_dir}/candidates_summary.csv', "w") as h:
        columns = ["sample", "variants", "mean_freq", "mean_depth", "exclusive_consensus", "exclusive_min",
                   "bad_quality"] # , "lineages"
        h.write("\t".join(columns) + "\n")

        for sample_name, df in dfs.items():
            sample_summary = {
                "sample": sample_name, "variants": len(df),
                "mean_freq": round(df.freq_min.mean(), 2), "mean_depth": round(df.depth_min.mean(), 2),
                "exclusive_consensus": len(df[df.exclusive_consensus]),
                "exclusive_min": len(df[df.exclusive_min]),
                #"lineages": " ".join([f'{x[0]}|{x[1]}|{x[2]}' for x in data["sample_lineages"][sample_name]]),

                "bad_quality": data["badquality_samples"][sample_name]
                if sample_name in data["badquality_samples"] else 0
            }
            df.to_csv(f'{output_dir}/report_{sample_name}.csv', index=False)
            h.write(("\t".join([str(sample_summary[c]) for c in columns]) + "\n"))

    plt.figure(figsize=(15, 10))
    plt.xlabel("Samples")
    plt.ylabel("Min Variants Freqs")
    plt.xticks(rotation=90)
    plt.boxplot(x=[list(dfs[c].freq_min) for c in candidates], labels=candidates)
    plt.savefig(f'{output_dir}/candidates_freqs.png')
    plt.close()

    print(f"Report Complete: {len(dfs)} candidate/s were processed")

    plt.figure(figsize=(15, 10))
    plt.xlabel("Variants")
    plt.ylabel("Samples Count")
    plt.xticks(rotation=90)
    plt.barpolt([len(v) for v in min_variant_samples.values()])
    plt.savefig(f'{output_dir}/variants_hist.png')
    plt.savefig(f'{output_dir}/variants_hist.eps',format="eps")
    plt.close()

    with open(f'{output_dir}/variants_list.csv', "w") as h:
        columns = ["variant", "ann" ,"consensus", "lowfreq", "in_candidate", "candidate_list"]
        h.write("\t".join(columns) + "\n")
        for var_str, samples in sorted(variant_samples.items(), key=lambda x: int(x[0].split("_")[0])):
            in_candidates = set(samples) & set(candidates)
            row = {"variant": var_str,"ann":ann_variants[var_str] if var_str in ann_variants else "",
                   "consensus": len([k for k, v in samples.items() if v == 1]),
                   "lowfreq": len([k for k, v in samples.items() if v != 1]),
                   "in_candidate": len(in_candidates),
                   "candidate_list": ",".join(in_candidates)}
            h.write(("\t".join([str(row[c]) for c in columns]) + "\n"))

    not_candidate_sample_variants = []

    for pos, sample_data in data["entries_data"].items():

        not_candidate_sample_variant = False
        vnoncandidates = []
        vcandidates = []
        freqs = []
        depths = []
        min_var_mut = ""
        for sample, (
                consensus_variant, min_variant, gts, ads, ref, (gene, gene_nt, gene_aa),
                lineages) in sample_data.items():

            if min_variant and min_variant[0]:
                min_var_mut = f'{pos}_{ref}_{min_variant[0]}'
                depths.append(ads[min_variant[0]])
                min_variant_samples[min_var_mut].append(sample)
                min_var_freq = min_variant[1]
                freqs.append(min_var_freq)
                if sample not in candidates:
                    not_candidate_sample_variant = True
                    vnoncandidates.append(sample)
                else:
                    vcandidates.append(sample)
        if not_candidate_sample_variant and min_var_mut:
            not_candidate_sample_variants.append({"pos": pos, "key": min_var_mut,  "depths": depths,
                                                  "freq": freqs, "candidates": vcandidates,
                                                  "non_candidates": vnoncandidates
                                                     , "gene": gene, "gene_nt": gene_nt, "gene_aa": gene_aa})

    pd.DataFrame(not_candidate_sample_variants).to_csv(f'{output_dir}/non_candidate_variants_list.csv')

    return candidates


if __name__ == '__main__':
    import argparse
    import subprocess as sp
    from glob import glob
    import gzip
    from collections import defaultdict

    parser = argparse.ArgumentParser(description='Process minory variants')

    subparsers = parser.add_subparsers(help='commands', description='valid subcommands', dest='command', required=True)

    cmd = subparsers.add_parser('download', help='Downloads reference and sample data')
    cmd.add_argument('-o', '--output_folder', default="./data")
    cmd.add_argument('-v', '--verbose', action='store_true')

    cmd = subparsers.add_parser('bam2vcf', help='calls haplotype caller for each bam sample')
    cmd.add_argument('-bams', '--bams_folder', required=True, help="directory were bam files are located")
    cmd.add_argument('-ref', '--reference', default="data/MN996528.fna",
                     help='fasta file. Can be gziped. Default "data/MN996528.fna"')
    cmd.add_argument('-o', '--output', default="results/", help="output dir. Default resutls/")
    cmd.add_argument('-v', '--verbose', action='store_true')

    cmd = subparsers.add_parser('merge_vcfs', help='joins all haplotype calls in one gvcf')
    cmd.add_argument('-vcfs', '--vcfs_dir', required=True, help="directory were bam files are located")
    cmd.add_argument('-ref', '--reference', default="data/MN996528.fna",
                     help='fasta file. Can be gziped. Default "data/MN996528.fna"')
    cmd.add_argument('-o', '--output',
                     default='results/variants.vcf.gz', help='output file. Default "results/variants.vcf.gz"')
    cmd.add_argument('-v', '--verbose', action='store_true')

    cmd = subparsers.add_parser('minor_freq_vars', help='gets a list of minor freq variants')
    cmd.add_argument('-mad', '--min_allele_depth', default=10, type=int,
                     help='Minimun allele read depth. Default 10')
    cmd.add_argument('-mc', '--min_coverage', default=0.8, type=float,
                     help='max percentaje N to discard a position. Between 0.0-1.0 Default 0.8'
                     )
    cmd.add_argument('-mf', '--min_freq', default=0.2, type=float,
                     help='minimun minority variant frequency. Between 0.0-1.0 Default 0.8')

    cmd.add_argument('--badq_strain_ns_threshold', default=1000, type=int,
                     help='sets the threshold (Ns) to tag a sample as a bad quality one')

    cmd.add_argument('--lineage_json', default='data/lineage_variants.json',
                     help='JSON file created from outbreak.info API')
    cmd.add_argument('--vcf', required=True, help="Multi Sample VCF. GT and AD fields are mandatory")
    cmd.add_argument('--out', default="results/data.json", help="Output data")
    cmd.add_argument('-v', '--verbose', action='store_true')

    cmd = subparsers.add_parser('report', help='comparative analysis between samples')
    cmd.add_argument('--data', required=True,
                     help='JSON file created by "minor_freq_vars"')

    cmd.add_argument('--min_lowfreq', default=None, type=int,
                     help='Instead of using a standard deviation to classify a sample as a coinfection "candidate"'
                          'we use the number of minority variants, as a hard limit. Should be used in a known set of samples '
                          'or if the fact that all samples are candidates is known beforehand. '
                          'If min_lowfreq is active, deviation_lowfreq value is ignored.'
                     )
    cmd.add_argument('--deviation_lowfreq', default=1,
                     help='Coinfection candidates are determined if the number of minority variants are greater than'
                          'N (this parameter) deviations from the mean. Default 1. If min_lowfreq is active, this parameter is ignored '
                     )

    cmd.add_argument('-mad', '--min_allele_depth', default=10, type=int,
                     help='Minimun allele read depth. Default 10')

    cmd.add_argument('--out_dir', default="./results")
    cmd.add_argument('-v', '--verbose', action='store_true')

    cmd = subparsers.add_parser('minconsensus', help='creates sequences using minor frequency variants')
    cmd.add_argument('--data', required=True,
                     help='JSON file created by minor_freq_vars')
    cmd.add_argument('--vcf', required=True, help="Multi Sample VCF. GT and AD fields are mandatory")
    cmd.add_argument('--out_dir', default="./results")
    cmd.add_argument('-ref', '--reference', default="data/MN996528.fna",
                     help='fasta file. Can be gziped. Default "data/MN996528.fna"')
    cmd.add_argument('-v', '--verbose', action='store_true')

    args = parser.parse_args()

    if args.verbose:
        os.environ["verbose"] = "y"

    if args.command == 'download':
        download(args)

    if args.command == 'bam2vcf':

        if not os.path.exists(args.reference):
            sys.stderr.write(f"'{args.reference}' does not exists")
            sys.exit(1)

        outfolder = args.output
        if not os.path.exists(outfolder):
            os.makedirs(outfolder)
        if not os.path.exists(outfolder):
            sys.stderr.write(f"'{outfolder}' could not be created")
            sys.exit(1)

        for bam_file in glob(args.bams_folder + "/*.bam"):
            sample = bam_file.split("/")[-1].split(".bam")[0]
            e(f"samtools index  {bam_file}")
            cmd = f"""/gatk/gatk  HaplotypeCaller -ERC GVCF -R {args.reference} \
                -ploidy 2 -I {bam_file} --output-mode EMIT_ALL_CONFIDENT_SITES -O {args.output}/{sample}.g.vcf.gz"""
            e(cmd)

    if args.command == 'merge_vcfs':
        merge_vcfs(args)

    if args.command == 'minor_freq_vars':
        if os.path.exists(args.lineage_json):
            with open(args.lineage_json) as h:
                lineage_data = json.load(h)

        else:
            lineage_data = {}
        variant_filter(args, lineage_data)

    if args.command == 'report':
        assert os.path.exists(args.out_dir), f'"{args.out_dir}" could not be created'
        candidates = comparative_analysis(args.data, args.out_dir, args.deviation_lowfreq, args.min_lowfreq,
                                          min_depth=args.min_allele_depth)
        # --sequences_dir ./results/sequences/

    if args.command == 'minconsensus':
        with open("data/combined.vcf") as h, open("./min_seqs.fasta", "w") as output:
            aln(h, output, refseq=str(bpio.read("data/MN996528.fna", "fasta").seq),
                included_samples=candidates)
