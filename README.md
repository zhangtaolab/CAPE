# CAPE

The computational pipeline of CAPE (<u>C</u>RISPR-Cas12<u>a</u> <u>p</u>romoter <u>e</u>diting)


## Prerequisition

1. Python >= 3.5
2. Open chromatin data (profiles in BigWig format, peaks in BED format)
3. TF binding motifs (identified by FIMO, matrix files are from PlantTFDB or JARSPR)
4. Sequence conservation (Scores are from PhastCons/mVISTA, or manually calculate with PHAST package)
5. Genome annotation file (in GFF3 format) and chromosome sizes file
6. (Optional) H3K27ac histone modification profile (BigWig format), genomic variations and phenotypes from rice3K/RFGB/MBKBase/etc.

## Install

```bash
# Install CAPE dependencies
git clone https://github.com/zhangtaolab/CAPE.git
cd CAPE
pip install -r requirements.txt

# Run test for single gene
cd test/single
python ../../single.py config.ini
```

### Run the pipeline for single gene

```bash
# Modify the config.ini file
[General]
workdir = results
binsize = 10
step = 10
upstream = 2000
slop = 200
withutr = 0
threads = 16

[Features]
ocfiles = Rice_leaf_DNase.bw,Rice_callus_DNase.bw
ocpeaks = TIGR7_DHSs.bed
ptmfiles = rice_H3K27ac.bw
motifs = genome_wide_motifs_JASPAR.bed
cnss = Osj_PhastCons.bedGraph
genopheno = 
phenodata = 

[Genes]
gene_file = gene.bed
gff_file = 
chrom_sizes = osativa_7.chrom.sizes
```

```bash
# Run the pipeline
python single.py config.ini
```

### Run the pipeline for whole genome genes

```bash
# Modify the config.ini file
[General]
workdir = results
binsize = 10
step = 10
upstream = 2000
slop = 200
withutr = 0
threads = 16

[Features]
ocfiles = Rice_leaf_DNase.bw,Rice_callus_DNase.bw
ocpeaks = TIGR7_DHSs.bed
ptmfiles = rice_H3K27ac.bw
motifs = genome_wide_motifs_JASPAR.bed
cnss = Osj_PhastCons.bedGraph
genopheno = 
phenodata = 

[Genes]
gene_file = 
gff_file = TIGR7_all.gff3
chrom_sizes = osativa_7.chrom.sizes
```

```bash
# Run the pipeline
python batch.py config.ini
```

## Output

All output files are stored in the workdir defined in the config.ini file.  
A folder will be created for each gene analyzed.  
In the output gene folder, several files are generated:  
1. analysis_region.bed (File records the analyzed regions in the genome for this gene)
2. OCpeaks_*_raw.bed (Open chromatin regions overlap with the analysis region)
3. OCscores*.bedGraph (Open chromatin scores for the analysis region, raw means raw scores from BigWig file, others are normalized in range 0 to 1)
4. motifs*.bedGraph (Raw file contains motifs identified in the analysis region, another file is the normalized motifs scores)
5. CNS*.bedGraph (Raw file contains raw conserved score in the analysis region, another file is the normalized CNS scores)
6. PTM*.bedGraph (H3K27ac profile for the analysis region, scores from BigWig file, others are normalized in range 0 to 1)
7. aggregate.bedGraph (The aggregate scores (AS) calculated from all above features)
8. key_regions_*.bed (Merged file means merged key regions when two key regions are adjacent)
9. core_regions.bed (Core regions which have high AS within the key regions)  
( **Optional:** if CRISPR edited phenotype data are provided, also export the statistical analysis results. )
10. phenoscores_*.bedGraph (phenotype scores, measured by kmeans-like method)
11. scores_by_sample.txt (Features scores and aggregate scores for each CRISPR edited sample)
12. plot_scores.txt (Comparison between phenotype difference and estimated scores)
13. statistics.txt (Cutoff for defining key regions and significance analysis)

