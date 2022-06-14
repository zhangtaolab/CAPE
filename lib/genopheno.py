from Bio import pairwise2
import re
from tqdm import tqdm
from itertools import chain
import numpy as np
import pandas as pd
from scipy import stats


def load_fasta(seqfile):
    seqinfo = {}
    with open(seqfile, "r") as infile:
        for line in infile:
            if line.startswith(">"):
                info = line.split("|")
                name = info[0][1:]
                if name == "REF":
                    sample_num = 0
                else:
                    sample_num = int(info[1].split(":")[1])
            else:
                seq = line.rstrip()
                if name != "REF":
                    if seq == seqinfo["REF"]["seq"]:
                        refid = name
                seqinfo[name] = {}
                seqinfo[name]["seq"] = seq
                seqinfo[name]["num"] = sample_num
    infile.close()
    return seqinfo, refid


def parse_alignment(alignment):
    aligninfo = {}
    refseq = alignment[0]
    altseq = alignment[1]
    indels = re.compile(r'-+')
    inspos = {}
    for i in range(len(refseq.replace("-", ""))):
        inspos[i] = 0
    for m in indels.finditer(refseq):
        start = m.span()[0]
        end = m.span()[1]
        for j in range(start+1, len(inspos)):
            inspos[j] += end - start
        aligninfo[start-inspos[start]] = {}
        aligninfo[start-inspos[start]]["ref"] = refseq[start-inspos[start]]
        aligninfo[start-inspos[start]]["alt"] = altseq[start-inspos[start]:end]
    for m in indels.finditer(altseq):
        start = m.span()[0]
        end = m.span()[1]
        aligninfo[start-inspos[start]] = {}
        aligninfo[start-inspos[start]]["ref"] = refseq[start-inspos[start]:end]
        aligninfo[start-inspos[start]]["alt"] = altseq[start-inspos[start]]
    for i in range(len(refseq)):
        refbase = refseq[i]
        altbase = altseq[i]
        if refbase != altbase:
            if refbase != "-" and altbase != "-":
                aligninfo[i-inspos[i]] = {}
                aligninfo[i-inspos[i]]["ref"] = refbase
                aligninfo[i-inspos[i]]["alt"] = altbase
    return aligninfo


def pairwise_alignment(seqfile):
    seqinfo, refid = load_fasta(seqfile)
    mutinfo = {}
    refseq = seqinfo["REF"]["seq"]
    total_num = sum([seqinfo[x]["num"] for x in seqinfo])
    count = 0
    for sample in tqdm(seqinfo, desc="Finding mutations"):
        if sample == "REF":
            continue
#         if count >= 5:
#             break
        altseq = seqinfo[sample]["seq"]
        num = seqinfo[sample]["num"]
        ratio = round(num / total_num, 4)
        alignments = pairwise2.align.globalms(refseq, altseq, 2, -1, -1.5, -.5)
        # print(sample, ratio, alignments[0], sep="\n")
        mutinfo[sample] = {}
        mutinfo[sample]["ratio"] = ratio
        mutinfo[sample]["alignment"] = parse_alignment(alignments[0])
        count += 1
    return mutinfo, refid


def mut2pos(seqfile):
    mutinfo, refid = pairwise_alignment(seqfile)
    vcfinfo = {}
    for sample in mutinfo:
        for pos in mutinfo[sample]["alignment"]:
            refbase = mutinfo[sample]["alignment"][pos]["ref"]
            altbase = mutinfo[sample]["alignment"][pos]["alt"]
            if altbase in ["a", "c", "g", "t", "n"]:
                altbase = altbase.upper()
                homozygous = 1
            else:
                homozygous = 0
            if pos not in vcfinfo:
                vcfinfo[pos] = {}
                vcfinfo[pos]["ref"] = refbase
                vcfinfo[pos]["alt"] = {}
            if altbase not in vcfinfo[pos]["alt"]:
                vcfinfo[pos]["alt"][altbase] = {}
            ratio = mutinfo[sample]["ratio"]
            vcfinfo[pos]["alt"][altbase][sample] = [ratio, homozygous]
    return vcfinfo, refid


def load_phenodata(phenodata):
    gid_info = {}
    with open(phenodata, "r") as infile:
        for line in infile:
            if line.startswith("Genotype_ID"):
                continue
            info = line.rstrip().split("\t")
            sample = info[0]
            if len(info) > 1:
                values = list(map(float, [x for x in info[1].split(", ")]))
                if len(values) > 1:
                    gid_info[sample] = values
    infile.close()
    return gid_info


def link_genopheno(genoinfo, seqfile, phenodata):
    posinfo, refid = mut2pos(seqfile)
    phenoinfo = load_phenodata(phenodata)
    startpos = genoinfo.start
    outfile = seqfile.replace(".fasta", "_geno_pheno.txt")
    outf = open(outfile, "w")
    print("name", "pos", "ref", "alt", "value", "avg", "sd", sep="\t", file=outf)
    for pos in sorted(posinfo):
        pos_abs = pos + startpos
        ref = posinfo[pos]["ref"]
        flag = 0
        for alt in posinfo[pos]["alt"]:
            input_lst = [phenoinfo[x] for x in posinfo[pos]["alt"][alt] if x in phenoinfo]
            values = list(chain(*input_lst))
            name = str(pos_abs) + "_" + ref + "/" + alt
            if values:
                flag = 1
                avg_value = round(np.average(values), 4)
                sd = round(np.std(values), 4)
                for value in values:
                    print(name, pos_abs, ref, alt, value, avg_value, sd, sep="\t", file=outf)
        if flag:
            ref_values = phenoinfo[refid]
            ref_avg = round(np.average(ref_values), 4)
            ref_sd = round(np.std(ref_values), 4)
            ref_name = str(pos_abs) + "_" + ref + "/" + ref
            for value in ref_values:
                print(ref_name, pos_abs, ref, ref, value, ref_avg, ref_sd, sep="\t", file=outf)
    outf.close()
    return outfile


def output_genopheno(genoinfo, seqfile, phenodata, outfile = "", startpos = 0):
    infile = link_genopheno(genoinfo, seqfile, phenodata)
    geno_pheno = pd.read_table(infile)
    chrom = genoinfo.chrom
    if not outfile:
        outfile = infile.replace(".txt", ".bed")
    outf = open(outfile, "w")
    for pos in pd.unique(geno_pheno.pos):
        value_lst = []
        ref = pd.unique(geno_pheno[geno_pheno.pos==pos].ref)
        for alt in pd.unique(geno_pheno[geno_pheno.pos==pos].alt):
            value_lst.append(geno_pheno[(geno_pheno.pos==pos) & (geno_pheno.alt==alt)].value.tolist())
        if len(value_lst) > 1:
            kruskal = stats.kruskal(*value_lst)
            statistic = kruskal[0]
            pvalue1 = kruskal[1]
            pvalue2 = -np.log10(pvalue1)
            # print(statistic, pvalue, pvalue2)
            real_pos = pos + startpos
            print(chrom, real_pos, real_pos+len(ref), statistic, pvalue1, pvalue2, sep="\t", file=outf)
    outf.close()
    return outfile

