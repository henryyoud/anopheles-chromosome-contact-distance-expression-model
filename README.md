# Anopheles Chorosome Contact Distance Expression Model

This R script models the impact of 3D chromosomal contact, linear genomic distance, and phylogenetic distance on pairiwse gene expression correlation. This script requires processed differntial expression RNAseq data.

Linear genomic distance and gene annotation were computed using the file "Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.7.gff3.gz" taken from VectorBase.

Chromosome contact data were extracted using Hi-C data taken from Lukyanchikova et al. (2022).

Phylogenetic distance p-values calculated using Mash (https://github.com/marbl/Mash).
