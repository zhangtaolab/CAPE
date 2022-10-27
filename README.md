# CAPE

The computational pipeline of CAPE (<ins>C</ins>RISPR-Cas12<ins>a</ins> <ins>p</ins>romoter <ins>e</ins>diting)


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

## Input (Feature data processing)

The instruction of how to generate feature data for calculation:  
1. Open chromatin data:  
(1) Raw sequencing data (from DNase-seq/ATAC-seq/MNase-seq) first align to reference genome by BWA/Bowtie2;  
(2) Call peaks from the alignment using Macs2/Genrich/F-seq2/Popera;  
(3) Generate profiles from the alignment (BigWig format, using DeepTools/F-seq2/Popera).  
2. TF binding motifs:  
(1) Download the TF PFM data from database (PlantTFDB/JASPAR/CisBP);  
(2) Find the occurrences of TF motifs in the genome by FIMO;  
(3) Merge results of all TF motifs (BED format, TFs from the same family can be merged into one).  
3. Sequence conservation:  
(1) Pre-calculated sequence conservation of plant genomes can be retrieved from PlantRegMap database;  
(2) If no existed result for the target genome, calculate conservation scores using multiple close related genomes with PHAST/mVISTA.  
4. H3K27ac histone modification:  
(1) Raw sequencing data (from ChIP-seq) first align to reference genome by BWA/Bowtie2;  
(2) Generate profiles from the alignment (BigWig format, using DeepTools).  
5. Relationships between genomic variations and phenotypes (GenoPheno):  
(1) Get the genotype data from public database, in FASTA format (for rice, using rice3K/RFGB/MBKBase/etc);  
(2) Get the corresponding phenotype data from public database.   
    (two column tab format, first column is Genotype_ID, second is Phenotype_Values separated by comma)  
6. Genome annotation file (BED/GFF3 format) is required for getting the promoter of target gene.  
7. Chromosome sizes file is required for converting input file format.  
   (two column tab format, first column is chromosome name, second is chromosome length)  

\* Note that H3K27ac and GenoPheno data are optional for analysis.

## Output

All output files are stored in the workdir defined in the config.ini file.  
A folder will be created for each gene analyzed.  
In the output gene folder, several files are generated:  
1. analysis_region.bed (File records the analyzed regions in the genome for this gene)
2. OCpeaks_*_raw.bed (Open chromatin regions overlap with the analysis region)
3. OCscores*.bedGraph (Open chromatin scores for the analysis region, suffix 'raw' means raw scores from BigWig file, others are normalized in range 0 to 1)
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

