import os
import numpy as np
from pybedtools import BedTool


def get_scores(geneinfo, scorefile, regionfile):
    score_bed = BedTool(scorefile)
    region_bed = BedTool(regionfile)
    scores = {}
    for interval in score_bed.intersect(region_bed, wo=True):
        info = str(interval).rstrip().split("\t")
        chrom = info[0]
        score = float(info[3])
        region_start = info[5]
        region_end = info[6]
        name = "_".join([chrom, region_start, region_end])
        if name not in scores:
            scores[name] = []
        scores[name].append(score)
    return scores


def get_cores(geneinfo, scores):
    binsize = geneinfo.binsize
    cores = ""
    for region in scores:
        values = scores[region]
        if len(values) >= 5:
            cutoff = np.average(values)
        else:
            cutoff = 0
        chrom, rstart, rend = region.split("_")
        for i, score in enumerate(values):
            if score >= cutoff:
                start = int(rstart) + int(binsize * i)
                end = start + binsize
                cores += "\t".join([chrom, str(start), str(end)]) + "\n"
    core_regions = BedTool(cores, from_string=True).merge()
    return core_regions


def output_cores(geneinfo, scorefile, regionfile, minlen = 2, outfile = ""):
    scores = get_scores(geneinfo, scorefile, regionfile)
    cores = get_cores(geneinfo, scores)
    binsize = geneinfo.binsize
    core_regions = []
    if not outfile:
        outfile = regionfile.replace("key_regions_merged", "core_regions")
    outf = open(outfile, "w")
    for interval in cores:
        info = str(interval).rstrip().split("\t")
        chrom = info[0]
        start = int(info[1])
        end = int(info[2])
        if end - start >= binsize * minlen:
            core_regions.append([chrom, start, end])
            print(chrom, start, end, sep="\t", file=outf)
    outf.close()
    return core_regions

