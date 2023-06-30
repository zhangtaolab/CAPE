import os
import random
import pandas as pd
import numpy as np
from scipy import stats
from pybedtools import BedTool

from lib import misc


class Geneinfo():
    # Provide Gene infomation
    def __init__(self):
        self.gene = "Gene"
        self.alias = "NA"
        self.chrom = "chrom"
        self.start = 0
        self.end = 1000
        self.strand = "+"
        self.binsize = 10
        self.step = 10


def kmeans_like_diff(wt, control):
    # Calculate the k-means-like phenotype difference between mutants and WT
    ## Do not consider the length of mutations
    y_avg = np.average(wt)
    diff = 0
    for x in control:
        diff += abs(x - y_avg)
    diff_score = diff / len(control)
    return diff_score


def kmeans_like_diff2(wt, control, binsize):
    # Calculate the k-means-like phenotype difference between mutants and WT
    ## Consider the influence of length of mutations
    y_avg = np.average(wt)
    diff = 0
    for x in control:
        diff += abs(float(x[4]) - y_avg) * ((int(x[2]) - int(x[1]))/binsize)
    diff_score = diff / len(control)
    return diff_score


def openchromatin_scores(geneinfo, bedfile, peakfile = "", outdir = "./", samplename = "openchromatin"):
    """
    Generate the open chromatin feature in specific bins.
    (Alternative data: ATAC-seq, DNase-seq, MNase-seq)

    Mandatory parameters:
    1. geneinfo - A class that defines the information of target gene
    2. bedfile - Open chromatin values in bedGraph format
    3. peakfile - Enrichment regions called from open chromatin data in BED format

    Alternative parameters:
    1. outdir - Output directory for saving the scores file (bedGraph format)
    """

    # Get gene info
    gene = geneinfo.gene
    genename = geneinfo.alias
    if genename == "NA":
        gene_alias = gene
    else:
        gene_alias = genename
    chromosome = geneinfo.chrom
    binstart = geneinfo.start
    binstop = geneinfo.end
    binsize = geneinfo.binsize
    step = geneinfo.step
    
    # Check output directory
    misc.check_outdir(outdir)
    if not os.path.exists(bedfile):
        smooth_openchromatin = {}
        return smooth_openchromatin

    # Convert BigWig file to bedGraph file
    # Load bedGraph file as bed file
    if peakfile:
        oc_peak = BedTool(peakfile)
    oc_score = BedTool(bedfile)

    # Calculate scores
    oc_info = []
    overlap_list = []
    posinfo = {}
    for i, pos in enumerate(range(binstart, binstop, step)):
        posinfo[i] = pos
        binbed = BedTool("\t".join([chromosome, str(pos), str(pos+binsize)])+"\n", 
                         from_string=True)
        score_in_bin = oc_score.intersect(binbed)
        if peakfile:
            peak_in_bin = oc_peak.intersect(binbed)
            overlap = [(int(str(x).split()[2])-int(str(x).split()[1]))/binsize for x in peak_in_bin]
            if overlap:
                overlap = 1
            else:
                overlap = 0.5
        else:
            overlap = 0.5
        sclst = [float(str(x).split()[3]) for x in score_in_bin]
        if sum(sclst):
            score = np.average(sclst)
        else:
            score = 0
        oc_info.append(score)
        overlap_list.append(overlap)

    # Smooth the scores
    max_score = max(oc_info)
    smooth_openchromatin = {}
    if max_score:
        oc_info = [x*overlap_list[i]/max_score for i,x in enumerate(oc_info)]
    else:
        return smooth_openchromatin
    outf = open(outdir + "/" + gene_alias + "/" + samplename + ".bedGraph", "w")
    for i in posinfo:
        pos = posinfo[i]
        score = oc_info[i]
        smooth_openchromatin[pos] = score
        print(chromosome, pos, pos+binsize, score, sep="\t", file=outf)
    outf.close()

    return smooth_openchromatin


def ptm_scores(geneinfo, bedfile, ocname, outdir = "./", samplename = "PTM", minratio = 0.2):
    """
    Generate the histone modification feature in specific bins.
    (Alternative data: ChIP-seq)

    Mandatory parameters:
    1. geneinfo - A class that defines the information of target gene
    2. bedfile - histone modification values in bedGraph format

    Alternative parameters:
    1. outdir - Output directory for saving the scores file (bedGraph format)
    """

    # Get gene info
    gene = geneinfo.gene
    genename = geneinfo.alias
    if genename == "NA":
        gene_alias = gene
    else:
        gene_alias = genename
    chromosome = geneinfo.chrom
    binstart = geneinfo.start
    binstop = geneinfo.end
    binsize = geneinfo.binsize
    step = geneinfo.step

    # Check output directory
    misc.check_outdir(outdir)

    # Load bedGraph file as bed file
    ptm_score = BedTool(bedfile)

    # Calculate scores
    ptm_info = []
    posinfo = {}
    for i, pos in enumerate(range(binstart, binstop, step)):
        posinfo[i] = pos
        binbed = BedTool("\t".join([chromosome, str(pos), str(pos+binsize)])+"\n", 
                         from_string=True)
        ptm_in_bin = ptm_score.intersect(binbed)
        score = np.average([float(str(x).split()[3]) for x in ptm_in_bin])
        if pd.isna(score):
            score = 0
        ptm_info.append(score)

    max_score = max(ptm_info)
    ptm_info = [x/max_score for x in ptm_info]
    
    # Get ratios from open chromatin results
    ocfile = outdir + "/" + gene_alias + "/" + ocname + ".bedGraph"
    oc_scores = BedTool(ocfile)
    oc_ratios = {}
    for interval in oc_scores:
        chrom, start, end, score = str(interval).rstrip().split("\t")
        if float(score) > minratio:
            oc_ratios[int(start)] = float(score)
        else:
            oc_ratios[int(start)] = minratio

    # Smooth the scores
    smooth_ptm = misc.smooth_scores_fill2(ptm_info, posinfo)
    outf = open(outdir + "/" + gene_alias + "/" + samplename + ".bedGraph", "w")
    ptm_scores = [(1-smooth_ptm[x])*oc_ratios[x] for x in smooth_ptm]
    max_score2 = max(ptm_scores)
    for pos in smooth_ptm:
#         score = smooth_ptm[pos]
        score = (1 - smooth_ptm[pos]) * oc_ratios[pos] / max_score2
        smooth_ptm[pos] = score
        print(chromosome, pos, pos+binsize, score, sep="\t", file=outf)
    outf.close()

    return smooth_ptm


def merge_reps(geneinfo, feature_list, outdir = "./", samplename = "merged"):
    """
    Merge the NGS feature in specific bins.
    (Alternative data: DNase-seq, ATAC-seq, ChIP-seq)

    Mandatory parameters:
    1. geneinfo - A class that defines the information of target gene
    2. feature_list - A list contains features need to be merged

    Alternative parameters:
    1. outdir - Output directory for saving the scores file (bedGraph format)
    """

    # Get gene info
    gene = geneinfo.gene
    genename = geneinfo.alias
    if genename == "NA":
        gene_alias = gene
    else:
        gene_alias = genename
    chromosome = geneinfo.chrom
    binsize = geneinfo.binsize

    # Check output directory
    misc.check_outdir(outdir)

    # Merge the feature scores
    scorelist = []
    cnt = 0
    for feature in feature_list:
        if not feature:
            continue
        if cnt:
            for i, pos in enumerate(feature):
                scorelist[i] += feature[pos]
        else:
            for pos in feature:
                scorelist.append(feature[pos])
        cnt += 1
    scores_merge = {}
    if not scorelist:
        return scores_merge
    outf = open(outdir + "/" + gene_alias + "/" + samplename + ".bedGraph", "w")
    for i, pos in enumerate(feature_list[0]):
        score = scorelist[i] / max(scorelist)
        scores_merge[pos] = score
        print(chromosome, pos, pos+binsize, score, sep="\t", file=outf)
    outf.close()

    return scores_merge


def motif_scores(geneinfo, bedfile, outdir = "./", flanking = 3):
    """
    Generate the TF motifs feature in specific bins.
    (Alternative data: Motif sites calculated by FIMO with PlantTFDB/JASPAR PWM files)

    Mandatory parameters:
    1. geneinfo - A class that defines the information of target gene
    2. bedfile - Motif positions in BED format

    Alternative parameters:
    1. outdir - Output directory for saving the scores file (bedGraph format)
    """

    # Get gene info
    gene = geneinfo.gene
    genename = geneinfo.alias
    if genename == "NA":
        gene_alias = gene
    else:
        gene_alias = genename
    chromosome = geneinfo.chrom
    binstart = geneinfo.start
    binstop = geneinfo.end
    strand = geneinfo.strand
    binsize = geneinfo.binsize
    step = geneinfo.step

    # Check output directory
    misc.check_outdir(outdir)

    # Load bed file
    tf_motif = BedTool(bedfile)

    # Calculate scores
    motif_density = []
    motif_info = []
    posinfo = {}
    count = 0
    for i, pos in enumerate(range(binstart, binstop, step)):
        posinfo[i] = pos
        binbed = BedTool("\t".join([chromosome, str(pos), str(pos+binsize)])+"\n", 
                         from_string=True)
        motif_in_bin = tf_motif.intersect(binbed)
        motif_lens = [int(str(x).split()[1]) for x in motif_in_bin]+[int(str(x).split()[2]) for x in motif_in_bin]
        if motif_lens:
            motif_density.append((max(motif_lens) - min(motif_lens)) / binsize)
        else:
            motif_density.append(0)
    bincount = len(motif_density)
    for i, score in enumerate(motif_density):
        if i > flanking:
            if i+flanking+1 > bincount:
                density = sum(motif_density[i-flanking:]) / (flanking+bincount-i)
            else:
                density = sum(motif_density[i-flanking:i+flanking+1]) / (2*flanking+1)
        else:
            density = sum(motif_density[:i+flanking+1]) / (flanking+i+1)
        motif_info.append(density)

    # Smooth the scores
    smooth_motif = misc.smooth_scores_fill2(motif_info, posinfo, minratio=1)
    max_score = max(smooth_motif.values())
    outf = open(outdir + "/" + gene_alias + "/motifs.bedGraph", "w")
    for pos in smooth_motif:
        score = smooth_motif[pos] / max_score
        smooth_motif[pos] = score
        print(chromosome, pos, pos+binsize, score, sep="\t", file=outf)
    outf.close()

    return smooth_motif


def cns_scores(geneinfo, bedfile, outdir = "./"):
    """
    Generate the conservation feature in specific bins.
    (Alternative data: Phastcons scores)

    Mandatory parameters:
    1. geneinfo - A class that defines the information of target gene
    2. bedfile - Conservation scores in BED format

    Alternative parameters:
    1. outdir - Output directory for saving the scores file (bedGraph format)
    """

    # Get gene info
    gene = geneinfo.gene
    genename = geneinfo.alias
    if genename == "NA":
        gene_alias = gene
    else:
        gene_alias = genename
    chromosome = geneinfo.chrom
    binstart = geneinfo.start
    binstop = geneinfo.end
    binsize = geneinfo.binsize
    step = geneinfo.step

    # Check output directory
    misc.check_outdir(outdir)
    if not os.path.exists(bedfile):
        smooth_cns = {}
        return smooth_cns

    # Load bed file
    cns = BedTool(bedfile)

    # Calculate scores
    cns_info = []
    posinfo = {}
    for i, pos in enumerate(range(binstart, binstop, step)):
        posinfo[i] = pos
        binbed = BedTool("\t".join([chromosome, str(pos), str(pos+binsize)])+"\n", 
                         from_string=True)
        cns_in_bin = cns.intersect(binbed)
        sclst = [float(str(x).split()[3]) for x in cns_in_bin]
        if sum(sclst):
            score = np.average(sclst)
        else:
            score = 0
        cns_info.append(score)

    # Smooth the scores
    smooth_cns = misc.smooth_scores2(cns_info, posinfo)
    outf = open(outdir + "/" + gene_alias + "/CNS.bedGraph", "w")
    for pos in smooth_cns:
        score = smooth_cns[pos]
        print(chromosome, pos, pos+binsize, score, sep="\t", file=outf)
    outf.close()

    return smooth_cns


def genopheno_scores(geneinfo, bedfile, outdir = "./"):
    """
    Generate the genotype and phenotype relationship feature in specific bins.
    (Alternative data: SNPs&Indels and Phenotype data)

    Mandatory parameters:
    1. geneinfo - A class that defines the information of target gene
    2. bedfile - genotype and phenotype relationship scores in BED format

    Alternative parameters:
    1. outdir - Output directory for saving the scores file (bedGraph format)
    """

    # Get gene info
    gene = geneinfo.gene
    genename = geneinfo.alias
    if genename == "NA":
        gene_alias = gene
    else:
        gene_alias = genename
    chromosome = geneinfo.chrom
    binstart = geneinfo.start
    binstop = geneinfo.end
    binsize = geneinfo.binsize
    step = geneinfo.step

    # Check output directory
    misc.check_outdir(outdir)

    # Load bed file
    genopheno = BedTool(bedfile)

    # Calculate scores
    genopheno_info = []
    posinfo = {}
    for i, pos in enumerate(range(binstart, binstop, step)):
        posinfo[i] = pos
        binbed = BedTool("\t".join([chromosome, str(pos), str(pos+binsize)])+"\n", 
                         from_string=True)
        genopheno_in_bin = genopheno.intersect(binbed)
        values = [float(str(x).split()[-1]) for x in genopheno_in_bin]
        # values = [x if x <= highest else highest for x in values]
        if values:
            score = sum(values)
        else:
            score = 0
        genopheno_info.append(score)

    # Smooth the scores
    smooth_genopheno = misc.smooth_scores_fill2(genopheno_info, posinfo)
    max_score = max(smooth_genopheno.values())
    outf = open(outdir + "/" + gene_alias + "/genopheno.bedGraph", "w")
    for pos in smooth_genopheno:
        score = smooth_genopheno[pos] / max_score
        smooth_genopheno[pos] = score
        print(chromosome, pos, pos+binsize, score, sep="\t", file=outf)
    outf.close()

    return smooth_genopheno


def aggregate_scores(geneinfo, scorelist, weightlist, outdir = "./"):
    """
    Generate the aggregate score in specific bins.

    Mandatory parameters:
    1. geneinfo - A class that defines the information of target gene
    2. scorelist - A list host multiple feature scores from different data
    3. weightlist - A list contains different weights assigned to different features
       (Should have the same order and numbers as scorelist)

    Alternative parameters:
    1. outdir - Output directory for saving the scores file (bedGraph format)
    """

    # Get gene info
    gene = geneinfo.gene
    genename = geneinfo.alias
    if genename == "NA":
        gene_alias = gene
    else:
        gene_alias = genename
    chromosome = geneinfo.chrom
    binsize = geneinfo.binsize

    # Check output directory
    misc.check_outdir(outdir)

    # Calculate scores
    outf = open(outdir + "/" + gene_alias + "/aggregate.bedGraph", "w")
    ziplist = zip(scorelist, weightlist)
    aggregate_info = {}
    total = sum(weightlist)
    for item in ziplist:
        scorelist = item[0]
        weight = item[1]
        for pos in scorelist:
            aggregate = scorelist[pos] * weight / total
            if pos in aggregate_info:
                aggregate_info[pos] += aggregate
            else:
                aggregate_info[pos] = aggregate
    if aggregate_info:
        max_score = max(aggregate_info.values())
    else:
        print(gene_alias)
        return aggregate_info
    if not max_score:
        return aggregate_info
    for pos in aggregate_info:
        score = aggregate_info[pos] / max_score
        print(chromosome, pos, pos+binsize, score, sep="\t", file=outf)
    outf.close()

    return aggregate_info


def phenodata_scores(geneinfo, bedfile, method = "kmeans1", outdir = "./", randbg = 0.02):
    """
    Calcuate the average phenodata value from multiple samples in specific bins.

    Mandatory parameters:
    1. geneinfo - A class that defines the information of target gene
    2. phenodata - Phenotype data of mutants in BED format
    chrom start end samplename avg_value

    Alternative parameters:
    1. method - Methods used for calculating phenotype difference between WT and mutants
    ["ratio", "stdev", "utest", "kmeans1", "kmeans2"]
    2. outdir - Output directory for saving the scores file (bedGraph format)
    """

    # Get gene info
    gene = geneinfo.gene
    genename = geneinfo.alias
    if genename == "NA":
        gene_alias = gene
    else:
        gene_alias = genename
    chromosome = geneinfo.chrom
    binstart = geneinfo.start
    binstop = geneinfo.end
    binsize = geneinfo.binsize
    step = geneinfo.step

    # Check output directory
    misc.check_outdir(outdir)

    # Load bed file
    phenodata = BedTool(bedfile)

    # Calculate scores
    methods = ["ratio", "stdev", "utest", "kmeans1", "kmeans2"]
    phenoinfo = []
    posinfo = {}
    for i, pos in enumerate(range(binstart, binstop, step)):
        posinfo[i] = pos
        binbed = BedTool("\t".join([chromosome, str(pos), str(pos+binsize)])+"\n", 
                         from_string=True)
        pheno_in_bin = phenodata.intersect(binbed)
        mutant_phenos = [float(str(x).split()[4]) for x in pheno_in_bin if str(x).split()[3] != "WT"]
        wt_phenos = [float(str(x).split()[4]) for x in pheno_in_bin if str(x).split()[3] == "WT"]
        if mutant_phenos:
            if method == methods[0]:
                score = np.average(mutant_phenos) / np.average(wt_phenos)
            elif method == methods[1]:
                score = np.std(wt_phenos + mutant_phenos)
            elif method == methods[2]:
                mannwhitneyu = stats.mannwhitneyu(wt_phenos, mutant_phenos)
                score = -np.log10(mannwhitneyu[1])
            elif method == methods[3]:
                score = kmeans_like_diff(wt_phenos, mutant_phenos)
            elif method == methods[4]:
                mutant_phenos = [str(x).split()[:5] for x in pheno_in_bin if str(x).split()[3] != "WT"]
                score = kmeans_like_diff2(wt_phenos, mutant_phenos, binsize)
            else:
                print("Cannot find this method. Available methods are:", methods)
        else:
            score = 0
        phenoinfo.append(score)

    # Smooth the scores
    max_score = max(phenoinfo)
    random.seed(81)
    phenoinfo = [max(x/max_score+random.uniform(-randbg, randbg), 0) if x else x for x in phenoinfo]
    # Output raw scores of phenotypes
    outraw = open(outdir + "/" + gene_alias + "/phenoscores_" + method + "_raw.bedGraph", "w")
    for i in posinfo:
        pos = posinfo[i]
        score = phenoinfo[i]
        print(chromosome, pos, pos+binsize, score, sep="\t", file=outraw)
    outraw.close()
    # Output smooth and gap-filled scores of phenotypes
    smooth_phenos = misc.smooth_scores1(phenoinfo, posinfo)
    max_score = max(smooth_phenos.values())
    min_score = min([x for x in smooth_phenos.values() if x])
    outf = open(outdir + "/" + gene_alias + "/phenoscores_" + method + ".bedGraph", "w")
    for pos in smooth_phenos:
        if smooth_phenos[pos]:
            score = (smooth_phenos[pos] - min_score) / (max_score - min_score)
        else:
            score = 0
        smooth_phenos[pos] = score
        print(chromosome, pos, pos+binsize, score, sep="\t", file=outf)
    outf.close()

    return smooth_phenos


def define_key_regions(geneinfo, aggregate, phenodata, threshold = 0, outdir = "./"):
    """
    Define the key regions of the target site.

    Mandatory parameters:
    1. geneinfo - A class that defines the information of target gene
    2. aggregate - Aggregate scores
    3. phenotypes - Phenotype scores
       (Should have the same order and numbers as scorelist)

    Alternative parameters:
    1. threshold - Bin with score above the threshold is defined as a key region
       (Default: average of aggregate scores)
    2. outdir - Output directory for saving the scores file (bedGraph format)

    Outputs:
    1. plot_scores - Phenotype and aggregate scores for R/ggplot2
    2. key_regions - Key regions in the target site
    3. stats - Statistics of Pearson correlation and differential significance
    """

    # Get gene info
    gene = geneinfo.gene
    genename = geneinfo.alias
    if genename == "NA":
        gene_alias = gene
    else:
        gene_alias = genename
    chromosome = geneinfo.chrom
    binsize = geneinfo.binsize

    # Check output directory
    misc.check_outdir(outdir)

    # Define the cutoff
    if threshold:
        cutoff = threshold
        cutoff_dev = 0
    else:
        cutoff_dev = np.std(list(aggregate.values()))
        cutoff = np.average(list(aggregate.values()))

    # Classify key regions and other regions
    key_regions = []
    aggregate_all = []
    phenotype_all = []
    for pos in aggregate:
        score = aggregate[pos]
        if score >= cutoff:
            aggregate_all.append(score)
            key_regions.append([pos, score])
        else:
            aggregate_all.append(score)

    # Output key regions info
    merged_regions = misc.merge_regions(key_regions, geneinfo)
    raw_file = outdir + "/" + gene_alias + "/key_regions_raw.bed"
    outregion1 = open(raw_file, "w")
    merged_file = outdir + "/" + gene_alias + "/key_regions_merged.bed"
    outregion2 = open(merged_file, "w")
    for region in key_regions:
        pos, score = region
        print(chromosome, pos, pos+binsize, score, sep="\t", file=outregion1)
    for lst in merged_regions:
        print("\t".join(list(map(str, lst))), file=outregion2)
    outregion1.close()
    outregion2.close()

    if os.path.exists(phenodata):
        # Calculate statistical values
        outf = open(outdir + "/" + gene_alias + "/plot_scores.txt", "w")
        print("sample", "group", "ratio", "difference", sep="\t", file=outf)
        outstat = open(outdir + "/" + gene_alias + "/statistics.txt", "w")
        # Cutoff of key regions definition
        print("Cutoff for defining key regions: %s" % cutoff, file=outstat)
        print("Cutoff deviation: %s" % cutoff_dev, file=outstat)
        # Calculate difference
        pheno_all = new_stats(geneinfo, phenodata, outdir = outdir)
        mean_ratio = np.average([x[1] for x in pheno_all])
        min_ratio = min(([x[1] for x in pheno_all]))
        max_ratio = max(([x[1] for x in pheno_all]))
        high_edited = []
        high_edited2 = []
        low_edited = []
        low_edited2 = []
        for scores in pheno_all:
            diff = scores[0]
            ratio = scores[1]
            sample = scores[2]
            if ratio > mean_ratio:
                high_edited.append((diff))
                high_edited2.append((diff-min_ratio)/(max_ratio-min_ratio))
                print(sample, "high", ratio, diff, sep="\t", file=outf)
            else:
                low_edited.append(diff)
                low_edited2.append((diff-min_ratio)/(max_ratio-min_ratio))
                print(sample, "low", ratio, diff, sep="\t", file=outf)
        outf.close()
        
        phe_high = np.average(high_edited2)
        phe_low = np.average(low_edited2)
        phe_ratio = phe_high / phe_low
        phe_pvalue = stats.mannwhitneyu(low_edited, high_edited)
        phe_pvalue2 = stats.ks_2samp(low_edited, high_edited, alternative="greater")
        phe_pvalue3 = stats.f_oneway(low_edited, high_edited)
        print("Phenotype differential ratio:", phe_ratio)
        print("Phenotype significance (U test):", phe_pvalue[1])
        print("Phenotype significance (KS test):", phe_pvalue2[1])
        print("Phenotype significance (ANOVA):", phe_pvalue3[1])
        print("Phenotype differential ratio:", phe_ratio, file=outstat)
        print("Phenotype significance (U test):", phe_pvalue[1], file=outstat)
        print("Phenotype significance (KS test):", phe_pvalue2[1], file=outstat)
        print("Phenotype significance (ANOVA):", phe_pvalue3[1], file=outstat)
        outstat.close()
    else:
        print("No Phenotype data detected, output key regions.")

    return key_regions


def new_stats(geneinfo, phenodata, outdir = "./", side="both"):

    pheno_bed = BedTool(phenodata)
    gene = geneinfo.gene
    genename = geneinfo.alias
    if genename == "NA":
        gene_alias = gene
    else:
        gene_alias = genename
    key_regions = outdir + "/" + gene_alias + "/key_regions_merged.bed"
    region_bed = BedTool(key_regions)
    regionlens = sum([int(str(x).split()[2])-int(str(x).split()[1]) for x in region_bed])
    intersect = pheno_bed.intersect(region_bed, wao=True)
    sample_values = {}
    for interval in intersect:
        info = str(interval).rstrip().split("\t")
        sample = info[3]
        pheno = float(info[4])
        if sample == "WT":
            wt_value = pheno
            continue
        length = int(info[-1])
        if sample not in sample_values:
            if side == "none":
                phenoscore = abs(pheno - wt_value)
            else:
                phenoscore = pheno - wt_value
            sample_values[sample] = [phenoscore, 0]
        sample_values[sample][1] += length / regionlens
    max_ratio = max([x[1] for x in sample_values.values()])
    mean_pheno = np.average([x[0] for x in sample_values.values()])
    if mean_pheno < 0:
        for s in sample_values:
            sample_values[s][0] *= -1
    scores_list = sorted([(sample_values[s][0], sample_values[s][1]/max_ratio, s) for s in sample_values], 
                         key=lambda x:x[1], reverse=True)

    return scores_list

