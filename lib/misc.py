import os
import numpy as np
import pandas as pd
from scipy import stats
import pyBigWig
from pybedtools import BedTool


def check_outdir(path):
    dirpath = os.path.abspath(os.path.dirname(path))
    if not os.path.exists(dirpath):
        print("Create directory:", dirpath)
        os.makedirs(dirpath)


def split_region(geneinfo):
    # Position info split by binsize and step
    step = geneinfo.step
    binstart = geneinfo.start
    binstop = geneinfo.end
    posinfo = {}
    for i, pos in enumerate(range(binstart, binstop, step)):
        posinfo[i] = pos

    return posinfo


def get_chrom_sizes(file):
    chrlens = {}
    with open(file) as infile:
        for line in infile:
            if line.startswith("#"):
                continue
            info = line.rstrip().split("\t")
            if len(info) == 2:
                chrom = info[0]
                length = info[1]
            else:
                chrom = info[0]
                length = info[2]
            chrlens[chrom] = int(length)
    return chrlens


def bigwig2bedGraph(bwfile, geneinfo, chrlens, outfile, ext = 50):
    # Convert the bigwig file from Deeptools/Popera to bedGraph file in single-base-pair resolution
    # Suitable for DNase-seq/ATAC-seq/MNase-seq/ChIP-seq
    bwin = pyBigWig.open(bwfile)
    chrom = geneinfo.chrom
    start = geneinfo.start
    end = geneinfo.end
    chrom_len = chrlens[chrom]
    check_outdir(outfile)
    outf = open(outfile, "w")
    for i in range(max(1, start-ext), min(end+ext, chrom_len)):
        try:
            value = bwin.values(chrom, i, i+1)[0]
            if np.isnan(value):
                value = 0
            print(chrom, i, i+1, value, sep="\t", file=outf)
        except:
            continue
    outf.close()


def fimo_filter(gfffile, matrixinfo, geneinfo, outfile, pcut = 1e-5, qcut = 1):
    # Filter FIMO results with p-value or q-value cutoff
    motif_family = {}
    # Matrix from JASPAR
    if matrixinfo.startswith("JASPAR"):
        with open(matrixinfo, "r") as infile:
            for line in infile:
                if line.startswith("MOTIF"):
                    info = line.rstrip().split()
                    motif_id = info[1]
                    motif_name = info[2]
                    motif_family[motif_id] = motif_name
    # Matrix from PlantTFDB
    else:
        with open(matrixinfo, "r") as infile:
            for line in infile:
                if line.startswith("#"):
                    continue
                info = line.rstrip().split()
                genename = info[0]
                family = info[1]
                motif_family[genename] = family
    # Get gene info
    chrom = geneinfo.chrom
    begin = geneinfo.start
    # Output filtered motifs
    motif_list = []
    with open(gfffile, "r") as infile:
        for line in infile:
            if line.startswith("#"):
                continue
            info = line.rstrip().split("\t")
            start = int(info[3]) + begin
            end = int(info[4]) + begin
            strand = info[6]
            desc = info[8].split(";")
            motif_id = desc[0].split("=")[1]
            motif_name = motif_family[motif_id]
            pvalue = float(desc[2].split("=")[1])
            qvalue = float(desc[3].split("= ")[1])
            if pvalue <= pcut and qvalue <= qcut:
                motif_list.append([chrom, start, end, motif_name, ".", strand, pvalue, qvalue])
    outf = open(outfile, "w")
    for lst in sorted(motif_list):
        print("\t".join(list(map(str, lst))), file=outf)
    outf.close()


def smooth_scores_fill2(info, posinfo, minscore=0.01, minratio=0.5):
    """
    Make the discrete score values smoothly (fill zero scores).

    Mandatory parameters:
    1. info - A list contains scores in different bins
    2. posinfo - Position information of each bin

    """

    # In case original score info be modified
    new_info = info.copy()
    minval = max(min([x for x in new_info if x]), minscore)*minratio
    zerocnt = 0
    flag = 0
    for i in posinfo:
        pos = posinfo[i]
        score = new_info[i]
        if i == 0:
            if score == 0:
                flag = 1
                zerocnt += 1
                zerostart = i
                continue
        elif i == len(new_info)-1:
            if score == 0:
                for j in range(zerostart+1, len(new_info), 1):
                    score1 = new_info[j-1]
                    score2 = new_info[j]
                    if score1 == new_info[zerostart]:
                        new_info[j] = np.average([score1*minratio, minval])
                    else:
                        new_info[j] = np.average([score1, minval])
        else:
            if score == 0:
                zerocnt += 1
                zerostart = i
                continue
            else:
                if flag:
                    for j in range(i, zerostart, -1):
                        new_info[j-1] = np.average([score*minratio, minval])
                    flag = 0
                else:
                    if zerocnt:
                        right = int(score*zerocnt/(score+new_info[zerostart]))
                        left = zerocnt - right
                        for j in range(zerostart+1, zerostart+left+1, 1):
                            score1 = new_info[j-1]
                            score2 = new_info[j]
                            if score1 == new_info[zerostart]:
                                new_info[j] = np.average([score1*minratio, minval])
                            else:
                                new_info[j] = np.average([score1, minval])
                        for k in range(i-1, zerostart+left, -1):
                            score1 = new_info[k]
                            score2 = new_info[k+1]
                            if score2 == score:
                                new_info[k] = np.average([minval, score2*minratio])
                            else:
                                new_info[k] = np.average([minval, score2])
                zerostart = i
                zerocnt = 0
    # smooth scores
    smooth_info = smooth_scores2(new_info, posinfo, keep_tails=False)

    return smooth_info


def smooth_scores_fill(info, posinfo):
    """
    Make the discrete score values smoothly (fill zero scores).

    Mandatory parameters:
    1. info - A list contains scores in different bins
    2. posinfo - Position information of each bin

    """

    new_info = info.copy()
    smooth_info = {}
    nonzero = [x for x in new_info if x]
    if sum(nonzero):
        minscore = min(nonzero)
        maxscore = max(new_info)
    else:
        return smooth_info
    # Set the minimum fill score
    if minscore / maxscore > 0.1:
        bottom = 0.1
    else:
        bottom = minscore
    for i in range(len(new_info)):
        if i:
            score0 = new_info[i-1]
            score1 = new_info[i]
            if not score1:
                for j in range(i+1, len(new_info)):
                    score2 = new_info[j]
                    if score2:
                        break
                    if j == len(new_info)-1 and score2 == 0:
                        score2 = bottom
                ranges = j - i
                diff1 = abs(score0 - bottom)
                diff2 = abs(score2 - bottom)
                total = diff1 + diff2
                if total:
                    mid = int(ranges * diff1 / total)
                else:
                    mid = 0
                # print(i, j, ranges, mid, score0, score2, diff1, diff2, sep="\t")
                if ranges > 1:
                    for k in range(mid):
                        new_info[i+k] = score0 - diff1 * (k+1)/(mid+1)
                    for k in range(mid+1, ranges):
                        new_info[i+k] = bottom + diff2 * (k-mid)/(ranges-mid)
                new_info[i+mid] = bottom
        else:
            score = new_info[i]
            if score:
                pass
            else:
                new_info[i] = bottom

    smooth_info = smooth_scores2(new_info, posinfo, keep_tails=False)

    return smooth_info


def smooth_scores1(info, posinfo, keep_tails=True):
    """
    Make the discrete score values smoothly.
    (Remove missing values between two scores)

    Mandatory parameters:
    1. info - A list contains scores in different bins
    2. posinfo - Position information of each bin

    """

    # In case original score info be modified
    new_info = info.copy()
    score_num = len(new_info)
    # Fill gap between two scores
    for i in range(score_num):
        score = new_info[i]
        if i == 0:
            tmp_score = score
            tmp_idx = i
        else:
            if score and tmp_score:
                interval = i - tmp_idx
                if interval > 1:
                    for n, j in enumerate(range(tmp_idx+1, i)):
                        new_info[j] = tmp_score + (score - tmp_score) * n / (i - tmp_idx)
                tmp_score = score
                tmp_idx = i
    smooth_info = {}
    if max(new_info):
        new_info = [x/max(new_info) for x in new_info]
    else:
        return smooth_info
    # Smooth the scores
    smooth_info = smooth_scores2(new_info, posinfo, keep_tails=keep_tails)

    return smooth_info


def smooth_scores2(info, posinfo, keep_tails=False):
    """
    Make the discrete score values smoothly.

    Mandatory parameters:
    1. info - A list contains scores in different bins
    2. posinfo - Position information of each bin

    Alternative parameters:
    1. keep_tails - Whether or not to keep the missing values in the two tails

    """

    # In case original score info be modified
    new_info = info.copy()
    smooth_info = {}
    if not max(new_info):
        return smooth_info
    score_num = len(new_info)
    begin = 0
    end = score_num
    # Find the two tails
    for i in range(end):
        if i:
            begin_avg = np.average(new_info[:i])
        else:
            begin_avg = new_info[i]
        if i == end-1:
            end_avg = new_info[i]
        else:
            end_avg = np.average(new_info[i:])
        if begin_avg == 0:
            begin = i
        if end_avg == 0:
            end = i
            break
    # Get average value in adjacent scores
    if not keep_tails:
        for i in range(begin, 0, -1):
            if i:
                if begin == score_num-1:
                    score = new_info[i]
                else:
                    score = (new_info[i-1] + new_info[i] + new_info[i+1]) / 3
            else:
                score = (new_info[i] + new_info[i+1]) / 2
            new_info[i] = score
        for i in range(end, score_num):
            if i < score_num - 1:
                score = (new_info[i-1] + new_info[i] + new_info[i+1]) / 3
            else:
                score = (new_info[i-1] + new_info[i]) / 2
            new_info[i] = score
    for i in range(begin, end):
        if i == begin:
            if begin == score_num-1:
                score = new_info[i]
            else:
                score = (new_info[i] + new_info[i+1]) / 2
        elif i == end - 1:
            score = (new_info[i-1] + new_info[i]) / 2
        else:
            score = (new_info[i-1] + new_info[i] + new_info[i+1]) / 3
        new_info[i] = score
    for i in posinfo:
        # provide real positions for smoothed scores
        pos = posinfo[i]
        smooth_info[pos] = new_info[i] / max(new_info)

    return smooth_info


def merge_regions(regions, geneinfo, minlen = 2, mindist = 1):
    # Filter and merge key regions
    chromosome = geneinfo.chrom
    binsize = geneinfo.binsize
    merged = {}
    for pos, score in regions:
        start = pos
        end = pos + binsize
        if not merged:
            tmppos = start
            merged[tmppos] = [end, [score]]
            tmp_end = end
            continue
        if start - tmp_end <= binsize * mindist:
            merged[tmppos][0] = end
            merged[tmppos][1].append(score)
        else:
            merged[start] = [end, [score]]
            tmppos = start
        tmp_end = end

    merged_regions = []
    for pos in merged:
        start = pos
        end = merged[pos][0]
        if end - start >= binsize * minlen:
            score = np.average(merged[pos][1])
            merged_regions.append([chromosome, start, end, score])

    return merged_regions


def calc_importance(phenotypes, scorelist, namelist, geneinfo, outdir="./", side="none"):
    # Calculate the correlation between phenodata and scores from different features
    ziplist = zip(scorelist, namelist)
    gene = geneinfo.gene
    genename = geneinfo.alias
    if genename == "NA":
        gene_alias = gene
    else:
        gene_alias = genename
    sample_scores = {}
    for item in ziplist:
        scores = item[0]
        name = item[1]
        score_bed = BedTool("\n".join(["\t".join(map(str, [geneinfo.chrom, x, x+geneinfo.binsize, scores[x]])) 
                                       for x in scores]), 
                            from_string=True)
        pheno_bed = BedTool(phenotypes)
        intersect = pheno_bed.intersect(score_bed, wo=True)
        fscores = {}
        for interval in intersect:
            info = str(interval).rstrip().split("\t")
            sample = info[3]
            if sample == "WT":
                wt_value = float(info[4])
                continue
            ratio = int(info[-1]) / geneinfo.binsize
            if sample in fscores:
                fscores[sample]["feature"] += float(info[-2]) * ratio
            else:
                fscores[sample] = {}
                if side == "none":
                    fscores[sample]["pheno"] = abs(float(info[4]) - wt_value)
                else:
                    fscores[sample]["pheno"] = float(info[4]) - wt_value
                fscores[sample]["feature"] = float(info[-2]) * ratio
        min_score = min([fscores[x]["feature"] for x in fscores])
        max_score = max([fscores[x]["feature"] for x in fscores])
        avg_pheno = np.average([fscores[x]["pheno"] for x in fscores])
        if avg_pheno < 0:
            for s in fscores:
                fscores[s]["pheno"] *= -1
        min_pheno = min([fscores[x]["pheno"] for x in fscores])
        max_pheno = max([fscores[x]["pheno"] for x in fscores])
        feature_scores = []
        pheno_scores = []
        for s in fscores:
            score1 = (fscores[s]["feature"]-min_score)/(max_score-min_score)
            feature_scores.append(score1)
            if side == "none":
                score2 = fscores[s]["pheno"]
            else:
                score2 = (fscores[s]["pheno"]-min_pheno)/(max_pheno-min_pheno)
            pheno_scores.append(score2)
            if s not in sample_scores:
                sample_scores[s] = {}
                sample_scores[s]["pheno"] = score2
            sample_scores[s][name] = score1
        pearson = stats.pearsonr(feature_scores, pheno_scores)
        print(name, "Pearson correlation:", pearson[0])
        
    outfile = outdir + "/" + gene_alias + "/scores_by_sample.txt"
    df = pd.DataFrame(sample_scores).T
    df.index.name = "sample"
    df.to_csv(outfile, sep="\t")

