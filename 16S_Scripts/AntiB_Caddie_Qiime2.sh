#!/bin/bash

cd /Users/lab/Desktop/Loomis_Caddies/Seq_Data

#IMPORT DATA; completed on 04/23/24 
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --input-path Seqs/ \
  --output-path Qiime2_outputs/CaddieD-sequencesFINAL.qza

cd /Users/lab/Desktop/Loomis_Caddies/Seq_Data/Qiime2_outputs

#VISUALIZE DATA IMPORT; completed on 4/16/24 
qiime demux summarize \
 --i-data CaddieD-sequencesFINAL.qza \
 --o-visualization CaddieD-sequencesFINAL.qzv

#FILTER BY Q SCORE IF USING DEBLUR (we aren't though, so note how is commented (#) out). 
#qiime quality-filter q-score \
# --i-demux CaddieD-sequencesR1.qza \
# --o-filtered-sequences CaddieD-filtered.qza \
# --o-filter-stats CaddieD-filter-stats.qza

#VISUALIZE FILTER STATS
#qiime metadata tabulate \
#  --m-input-file CaddieD-filter-stats.qza \
#  --o-visualization CaddieD-filter-stats.qzv

#PAIRED END DATA (Forward and Reverse) and DADA2; completed on 4/16/24 
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs CaddieD-sequencesFINAL.qza \
  --p-trunc-len-r 140 \
  --p-trunc-len-f 145 \
  --p-max-ee-f 2 \
  --p-max-ee-r 2 \
  --p-pooling-method "pseudo" \
  --p-chimera-method "consensus" \
  --o-representative-sequences CaddieD-paired-asv-sequencesFINAL.qza \
  --o-table CaddieD-paired-asv-tableFINAL.qza \
  --o-denoising-stats CaddieD-paired-denoising-statsFINAL.qza

##VISUALIZE PAIRED END DATA; completed on 4/16/24 
qiime metadata tabulate \
  --m-input-file CaddieD-paired-denoising-statsFINAL.qza \
  --o-visualization CaddieD-paired-denoising-statsFINAL.qzv

qiime feature-table tabulate-seqs \
  --i-data CaddieD-paired-asv-sequencesFINAL.qza \
  --o-visualization CaddieD-paired-asv-sequencesFINAL.qzv

qiime feature-table summarize \
  --i-table CaddieD-paired-asv-tableFINAL.qza \
  --m-sample-metadata-file 16S_CaddieD_metadata.tsv \
  --o-visualization CaddieD-paired-asv-tableFINAL.qzv

#SILVA 138; TRAIN -- make sure its for 515-806
#qiime feature-classifier fit-classifier-naive-bayes \
#  --i-reference-reads silva-138-99-seqs-515-806.qza \
#  --i-reference-taxonomy silva-138-99-tax-515-806.qza \
#  --o-classifier silva-138-99-515F-806R-nb-classifierFINAL.qza

#CLASSIFY AGAINST SILVA, PAIRED
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-99-515F-806R-nb-classifierFINAL.qza \
 --i-reads CaddieD-paired-asv-sequencesFINAL.qza \
  --output-dir CaddieD-paired-taxonomyFINAL

#VISUALIZE TAXONOMY PAIRED
qiime metadata tabulate \
  --m-input-file CaddieD-paired-taxonomyFINAL/classification.qza \
  --o-visualization CaddieD-paired-taxonomyFINAL.qzv

#BUILD PHYLOGENY PAIRED 
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences CaddieD-paired-asv-sequencesFINAL.qza \
  --o-alignment CaddieD-aligned-paired-seqsFINAL.qza \
  --o-masked-alignment CaddieD-masked-aligned-paired-seqsFINAL.qza \
  --o-tree CaddieD-unrooted-tree-paired-seqsFINAL.qza \
  --o-rooted-tree CaddieD-rooted-tree-paired-seqsFINAL.qza

#EXPORT PAIRED
qiime tools export \
  --input-path CaddieD-paired-asv-tableFINAL.qza \
  --output-path CaddieD_paired_dadatableFINAL

qiime tools export \
  --input-path CaddieD-paired-taxonomyFINAL/classification.qza \
  --output-path CaddieD_paired_taxaFINAL

#EXPORT PHYLOGENY
qiime tools export \
  --input-path CaddieD-rooted-tree-paired-seqsFINAL.qza \
  --output-path CaddieD-rooted-tree-paired-seqsFINAL

