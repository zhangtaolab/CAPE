#! /usr/bin/env python

##################################################
# CRISPR-Cas12a promoter editing (CAPE)          #
# Script: Single Mode                            #
##################################################

import os
import sys
import shutil
from glob import glob
import configparser
from time import time
from multiprocessing import Pool
from pybedtools import BedTool, cleanup

from lib import misc
from lib.features import *
from lib.cores import output_cores


class Features_info():
    # Provide Gene infomation
    def __init__(self):
        self.geneinfo = Geneinfo()
        self.feature = "feature"
        self.workdir = "results"
        self.outname = "name"
        self.slop = 200
        self.config = {}
        self.chrlens = {}


def get_gene_info(gene_file):

    genes_info = {}
    with open(gene_file) as infile:
        for line in infile:
            if line.startswith("#") or line.startswith("\n"):
                continue
            info = line.rstrip().split("\t")
            chrom = info[0]
            start = int(info[1])
            end = int(info[2])
            gene_name = info[3]
            strand = info[5]
            genes_info[gene_name] = [chrom, start, end, strand]
            break

    print("Genes infomation loaded.\n")

    return genes_info


def generate_regions(geneinfo, workdir, gene, chrlens):
    chrom = geneinfo.chrom
    start = geneinfo.start
    end = geneinfo.end
    strand = geneinfo.strand
    chrom_len = chrlens[chrom]
    outfile = workdir + "/" + gene + "/analysis_region.bed"
    misc.check_outdir(outfile)
    if os.path.exists(outfile):
        return outfile
    outf = open(outfile, "w")
    print(chrom, max(0, start), min(end, chrom_len), gene, '.', strand, 
          sep="\t", file=outf)
    outf.close()

    return outfile


def generate_features(Features_info):

    geneinfo = Features_info.geneinfo
    feature_file = Features_info.feature
    workdir = Features_info.workdir
    outname = Features_info.outname
    slop = Features_info.slop
    chrlens = Features_info.chrlens

    chrom = geneinfo.chrom
    start = geneinfo.start
    end = geneinfo.end
    gene = geneinfo.gene

    if "peak" in outname:
        outfile = workdir + "/" + gene + "/" + outname + "_raw.bed"
    else:
        outfile = workdir + "/" + gene + "/" + outname + "_raw.bedGraph"
    if os.path.exists(outfile):
        return outfile
    if feature_file.endswith(".bw") or feature_file.endswith(".bigwig"):
        misc.bigwig2bedGraph(feature_file, geneinfo, chrlens, outfile, ext = slop)
    else:
        target_bed = BedTool("\t".join([chrom, str(max(0, start-slop)), str(end+slop)]), from_string = True)
        feature_bed = BedTool(feature_file)
        feature_bed.intersect(target_bed, wa=True).moveto(outfile)
    
    cleanup()

    return outfile


def generate_features_from_large(inputfile, genes_info, upstream, slop, workdir, feature):
    
    basemap = {}
    existed = set()
    num = 0
    genelens = len(genes_info)
    for gene in genes_info:
        outfile = os.path.join(workdir, gene, feature+"_raw.bedGraph")
        if os.path.exists(outfile):
            filesize = os.path.getsize(outfile)
            if filesize > 10:
                existed.add(gene)
        chrom, start, end, strand = genes_info[gene][:4]
        if chrom not in basemap:
            basemap[chrom] = {}
        if strand == "+":
            for i in range(max(0, start-upstream-slop), start+slop+1):
                if i in basemap[chrom]:
                    basemap[chrom][i].append(gene)
                else:
                    basemap[chrom][i] = [gene]
        else:
            for i in range(end-slop, end+upstream+slop+1):
                if i in basemap[chrom]:
                    basemap[chrom][i].append(gene)
                else:
                    basemap[chrom][i] = [gene]
        print("%s / %s genes processed, %s existed genes." % (num, genelens, len(existed)), 
              end="\r")
        num += 1
    print("Load genes completed.", " "*30)
    genenums = []
    total_num = max(1, genelens-len(existed))
    outf = {}
    split = 500
    kept = split * 0.9
    tmp_cnt = 0
    tmp_mod = 0
    num = 0
    with open(inputfile) as infile:
        for line in infile:
            chrom, start, end = line.rstrip().split("\t")[:3]
            if chrom not in basemap:
                continue
            if feature == "CNS":
                s = int(start)
            else:
                s = int((int(start) + int(end)) / 2)
            if s in basemap[chrom]:
                genes = basemap[chrom][s]
                for gene in genes:
                    if gene in existed:
                        continue
                    else:
                        outfile = os.path.join(workdir, gene, feature+"_raw.bedGraph")
                        if gene not in outf:
                            outf[gene] = open(outfile, "w")
                            # try:
                            #     outf[gene] = open(outfile, "w")
                            # except:
                            #     opened = len(outf)
                            #     print("Processing %s, %s genes opened." % (gene, opened))
                            #     outf[gene] = open(outfile, "w")
                    print(line.rstrip(), file=outf[gene])
                    if gene not in genenums:
                        genenums.append(gene)
                cnt = len(genenums)
                remain = cnt % split
                mod = cnt // split
            if mod - tmp_mod > 0:
                st = max(0, int(split * (mod - 1) - kept - 1))
                ed = int(split * mod - kept)
                # print("#"*100+"\n", tmp_cnt, cnt, tmp_mod, mod, st, ed, genes, sep=", ")
                for j in genenums[st:ed]:
                    outf[j].close()
            if tmp_cnt != cnt:
                pct = round(cnt * 100 / total_num, 2)
                print(pct, "%", " output.", end="\r")
            tmp_cnt = cnt
            tmp_mod = mod
    print("All files output.")

    for gene in outf:
        outf[gene].close()

    return cnt


def run_analysis(feature_info):
    
    workdir = feature_info.workdir
    geneinfo = feature_info.geneinfo
    gene = feature_info.geneinfo.gene

    # Check if calculated
    # check = os.path.join(workdir, gene, "key_regions_merged.bed")
    # if os.path.exists(check):
    #     return (gene, 0)
    check = os.path.join(workdir, gene, "aggregate.bedGraph")
    if os.path.exists(check):
        filesize = os.path.getsize(check)
        if filesize > 10:
            return (gene, 0)

    # Open chromatin
    ocscores = glob(os.path.join(workdir, gene, "OCscores*_raw.bedGraph"))
    ocpeaks = glob(os.path.join(workdir, gene, "OCpeaks*_raw.bed"))
    # Calculate scores
    ocscorelist = []
    for idx, ocscorefile in enumerate(ocscores):
        if idx + 1 > len(ocpeaks):
            ocpeakfile = ""
        else:
            ocpeakfile = ocpeaks[idx]
        if len(ocscores) > 1:
            ocname = os.path.basename(ocscorefile).split("_raw")[0]
        else:
            ocname = "OCscores"
        scores_oc1 = openchromatin_scores(geneinfo, ocscorefile, ocpeakfile, 
                                          samplename = ocname, outdir = workdir)
        ocscorelist.append(scores_oc1)
    if len(ocscores) > 1:
        scores_oc = merge_reps(geneinfo, ocscorelist, samplename = "OCscores", outdir = workdir)
    else:
        scores_oc = scores_oc1

    # Histone modification
    ptmfiles = glob(os.path.join(workdir, gene, "PTM*_raw.bedGraph"))
    # Calculate scores
    ptmscorelist = []
    for ptmscorefile in ptmfiles:
        if len(ptmfiles) > 1:
            ptmname = os.path.basename(ptmscorefile).split("_raw")[0]
        else:
            ptmname = "PTMscores"
        scores_ptm1 = ptm_scores(geneinfo, ptmscorefile, ocname="OCscores",
                                 samplename = ptmname, outdir = workdir)
        ptmscorelist.append(scores_ptm1)
    if len(ptmfiles) > 1:
        scores_ptm = merge_reps(geneinfo, ptmscorelist, samplename = "PTMscores", outdir = workdir)
    else:
        scores_ptm = scores_ptm1

    # TF motifs
    motiffile = os.path.join(workdir, gene, "motifs_raw.bedGraph")
    # Calculate scores
    scores_motif = motif_scores(geneinfo, motiffile, outdir = workdir)

    # Conserved sequences
    cnsfile = os.path.join(workdir, gene, "CNS_raw.bedGraph")
    # Calculate scores
    scores_cns = cns_scores(geneinfo, cnsfile, outdir = workdir)

    # Genotype versus Phenotype (MBKbase)
    genopheno = os.path.join(workdir, gene, "genopheno_raw.bedGraph")
    # Calculate scores
    if os.path.exists(genopheno):
        scores_genopheno = genopheno_scores(geneinfo, genopheno, outdir = workdir) 
    else: 
        scores_genopheno = {}
    
    # Aggregate scores
    if scores_genopheno:
        scorelist = [scores_oc, scores_motif, scores_cns, scores_ptm, scores_genopheno]
        weightlist = [0.25, 0.2, 0.3, 0.1, 0.05]
    else:
        scorelist = [scores_oc, scores_motif, scores_cns, scores_ptm]
        weightlist = [0.25, 0.2, 0.3, 0.1]
    scores_aggregate = aggregate_scores(geneinfo, scorelist, weightlist, outdir = workdir)

    # Load phenodata from CRISPR-edited results
    phenodata = os.path.join(workdir, gene, "phenoscores_raw.bedGraph")
    # Calculate scores
    if os.path.exists(phenodata):
        scores_phenodata = phenodata_scores(geneinfo, phenodata, method = "kmeans2", 
                                            outdir = workdir)
    else:
        scores_phenodata = {}

    # Find the feature importance
    if scores_phenodata:
        namelist = ["DHS", "H3K27ac", "TF motif", "CNS", "GenoPheno", "Aggregate"]
        misc.calc_importance(phenodata, scorelist+[scores_aggregate], 
                             namelist, geneinfo, side="both", outdir = workdir)

    # Define key regions
    key_regions = define_key_regions(geneinfo, scores_aggregate, phenodata, 
                                     outdir = workdir)

    # Get the core of key regions
    scorefile = os.path.join(workdir, gene, "aggregate.bedGraph")
    regionfile = os.path.join(workdir, gene, "key_regions_merged.bed")
    core_regions = output_cores(geneinfo, scorefile, regionfile)
    
    cleanup()

    return (gene, 1)


def check_options(config):

    print("# Using the following options:")
    if config["General"]["workdir"]:
        config["General"]["workdir"] = os.path.abspath(config["General"]["workdir"])
    else:
        config["General"]["workdir"] = "results"
    misc.check_outdir(config["General"]["workdir"])
    for section in config.sections():
        for param in config.options(section):
            values = config[section][param]
            if section == "Features":
                if "," in values:
                    values = values.split(",")
                    for file in values:
                        if file and not os.path.exists(file):
                            print("# Error, cannot find the %s: %s" % (param, file))
                            sys.exit(1)
                else:
                    file = values
                    if file and not os.path.exists(file):
                        print("# Error, cannot find the %s: %s" % (param, file))
                        sys.exit(1)
            print("%s: %s" % (param, values))
    if int(config["General"]["threads"]) > os.cpu_count():
        config["General"]["threads"] = os.cpu_count()
    if int(config["General"]["slop"]) > 5e4:
        config["General"]["slop"] = 5e4
    if int(config["General"]["upstream"]) > 1e4:
        config["General"]["upstream"] = 1e4
    if int(config["General"]["binsize"]) > int(config["General"]["upstream"]) / 2:
        config["General"]["binsize"] = int(config["General"]["upstream"]) / 2
    if int(config["General"]["step"]) > int(config["General"]["binsize"]):
        config["General"]["step"] = int(config["General"]["binsize"])
    if config["Genes"]["gene_file"]:
        print("\n# Using Single mode.\n")
    
    return config


def main():

    # Load configs
    config = configparser.ConfigParser()
    if len(sys.argv) == 1:
        config_file = "config.ini"
    elif len(sys.argv) == 2:
        config_file = sys.argv[1]
    else:
        print("Usage:\n    python single.py [configfile]\n")
        sys.exit(1)
    config.read(config_file)

    config = check_options(config)
    workdir = config["General"]["workdir"]
    threads = int(config["General"]["threads"])
    slop = int(config["General"]["slop"])
    upstream = int(config["General"]["upstream"])
    binsize = int(config["General"]["binsize"])
    step = int(config["General"]["step"])
    gene_file = config["Genes"]["gene_file"]
    chrom_sizes = config["Genes"]["chrom_sizes"]

    # Load genes
    if gene_file:
        genes_info = get_gene_info(gene_file)
    else:
        print("No gene annotation file found, stop!")
        sys.exit(1)
    
    # Define the input numbers of multiprocessing list
    inputnum = 512
    if inputnum < threads:
        inputnum = threads
    else:
        roundnum = (inputnum // threads) * threads
        inputnum = int(max(roundnum, threads*4))

    # Load chromosome sizes
    chrlens = misc.get_chrom_sizes(chrom_sizes)

    # Define features information
    feature_map = {"ocfiles":"OCscores", "ocpeaks":"OCpeaks", "ptmfiles":"PTM", 
                   "motifs":"motifs", "cnss":"CNS", "genopheno":"genopheno", 
                   "phenodata":"phenoscores"}
    for item in config["Features"]:
        feature_files = config["Features"][item]
        if not feature_files:
            continue
        filelist = feature_files.split(",")
        count = 1
        for file in filelist:
            feature_infos = []
            num = 1
            for gene in genes_info:
                chrom = genes_info[gene][0]
                start = genes_info[gene][1]
                end = genes_info[gene][2]
                strand = genes_info[gene][3]
                feature_info = Features_info()
                feature_info.workdir = workdir
                feature_info.slop = slop
                feature_info.config = config
                feature_info.idx = num
                feature_info.geneinfo = Geneinfo()
                feature_info.geneinfo.gene = gene
                feature_info.geneinfo.chrom = chrom
                feature_info.geneinfo.strand = strand
                if strand == "+":
                    feature_info.geneinfo.start = start - upstream
                    feature_info.geneinfo.end = start - 1
                else:
                    feature_info.geneinfo.start = end
                    feature_info.geneinfo.end = end + upstream - 1
                feature_info.geneinfo.binsize = binsize
                feature_info.geneinfo.step = step
                num += 1
                # Output analyzed gene regions
                generate_regions(feature_info.geneinfo, workdir, gene, chrlens)
                feature_info.feature = file
                if len(filelist) > 1:
                    outname = feature_map[item] + "_" + str(count)
                else:
                    outname = feature_map[item]
                feature_info.outname = outname
                feature_info.chrlens = chrlens
                feature_infos.append(feature_info)
            count += 1
            # Generate features file
            time_st = time()
            file_suffix = file.split(".")[-1].lower()
            filesize = os.path.getsize(file)
            if file_suffix in ["bed", "bedgraph", "txt"] and filesize > 1e8:
                results = generate_features_from_large(file, genes_info, upstream, slop, 
                                                       workdir, outname)
            else:
                # Multiprocessing
                results = generate_features(feature_infos[0])
            time_ed = time()
            time_elapse = round(time_ed - time_st)
            print("Generate %s features files finished.\nUsing %ss" % (outname, time_elapse))
    
    # Perform analysis
    time_st = time()
    result = run_analysis(feature_infos[0])
    if result[1]:
        time_total = round(time() - time_st, 2)
        print("\nGene analysis finished using %ss. %s\n" % (time_total, " "*30))

    print("All the processes completed.", " "*10)



if __name__ == '__main__':

    try:
        main()

    except KeyboardInterrupt:
        sys.stderr.write("User interrupt\n")
        sys.exit(0)

