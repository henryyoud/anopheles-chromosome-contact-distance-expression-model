# This script is used to model the effects of:
# linear genomic distance,
# three-dimensional chromosomal contact,
# and gene relatedness
# on pairwise gene expression correlation

#### Setup ####

library(ggplot2)
library(WGCNA)
library(ggpubr)
library(strawr)
library(data.table)
library(dplyr)
library(tidyr)
library(plyr)
library(effectsize)
library(boot)

setwd("D:/OneDrive - Lancaster University/PR2/data/") # Change as appropriate

#### Process Expression Data ####

# Load the datasets
datasets <- setNames(nm = c('Bouake', 'Busia', 'Tiassale'))

load.expression.table <- function(dataset, expression.threshold){
  filepath <- paste('normalisedCounts/normCounts', dataset, '.tsv', sep = '')
  expression.table <- read.table(filepath, sep = '\t', header = T, row.names = 1)
  expression.table[rowSums(expression.table) >= expression.threshold, ]
}

# Since these are normalised expression data, I don't know how to set a read count threshold equivalent to what 
# I set for tsetse. Using 1000 removes around 2000 genes per gambiae dataset, which is probably about right
data.tables <- lapply(datasets, load.expression.table, 100)

# Take logs of the data
log.data.tables <- lapply(data.tables, function(M) log2(M + 1))

# Use the adjacency function from wgcna package to quickly calculate correlation between all pairs of genes
cor.matrices <- lapply(log.data.tables, function(M) adjacency(t(M), power = 1))

# Straighten out the matrices into lower triangle vectors 
cor.lower.tri <- lapply(cor.matrices, function(M) M[lower.tri(M)])

#### Distance Calculations ####

# Anopheles gambiae (PEST) gene annotation
agam.gff <- read.table(gzfile('Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.12.gff3.gz', 'r'), sep = '\t', quote = '') # Read file
colnames(agam.gff) <- c('chromosome', 'annotation', 'type', 'start', 'end', 'score', 'strand', 'frame', 'attribute') # Rename columns
agam.genes <- subset(agam.gff, grepl('gene', type)) # Subset just genes
agam.gene.names <- sub('^ID=', '', sub(';.*', '', agam.genes$attribute)) # Extract AGAP gene names
rownames(agam.genes) <- agam.gene.names # Assign gene names to row names
agam.genes$name = agam.gene.names # Assign gene names to name column

# We now want to create a matrix of distances between genes. For any two genes that share a scaffold, we take 
# the absolute value of the difference between their start points. 
# Create a matrix where all distances are assumed to be the maximum value (double the size of the largest 
# observed gene position)
# This arbitrary maximum value is retained if genes are on different chromosomes
agam.genes.by.chrom <- split(agam.genes, as.character(agam.genes$chromosome)) # Split by chromosome
agam.max.dist <- 2*max(agam.gff$end)
agam.dist.matrix <- matrix(agam.max.dist, nrow(agam.genes), nrow(agam.genes), dimnames = list(agam.gene.names, agam.gene.names))
agam.within.chrom.distances <- lapply(agam.genes.by.chrom, function(X) apply(X[, 'start', drop = F], 1, function(x) abs(x['start'] - X$start)))
for (wcd in agam.within.chrom.distances){
  agam.dist.matrix[colnames(wcd), colnames(wcd)] <- wcd
}

# Get a separate distance matrix for each dataset, with genes subsetted and ordered
dist.matrices <- c(lapply(cor.matrices, function(M) agam.dist.matrix[rownames(M), colnames(M)]))

# Straighten out the matrices into lower triangle vectors 
dist.lower.tri <- lapply(dist.matrices, function(M) M[lower.tri(M)])


#### Extract and Process Contact Data ####

# Hi-C data from Lukyanchikova et al. (20202) - check resolution and length of chromosomes
# Anopheles coluzzi 
strawr::readHicBpResolutions("AcolNg_V4.hic") # 5000 is minimum resolution
strawr::readHicChroms("AcolNg_V4.hic") # chromosome lengths

# Have to loop straw function as it only lets you read data in specified regions, not whole genome

# All chromosome combinations used in for loop
var1 = c("X", "X", "X", "X", "X", 
         "2R", "2R", "2R", "2R", 
         "2L", "2L", "2L", 
         "3R", "3R", 
         "3L")

var2 = c("X", "2R", "2L", "3R", "3L", 
         "2R", "2L", "3R", "3L", 
         "2L", "3R", "3L", 
         "3R", "3L", 
         "3L")

chromosomes = data.frame(var1, var2)

# Create simplified straw function (KR normalisation and expected contact frequency used)
straw_all = function(chrom1, chrom2, data){
  output = strawr::straw(norm = "KR", fname = data, 
                         chr1loc = chrom1, chr2loc = chrom2, 
                         unit = "BP", binsize = 5000, matrix = "expected")
}

# Empty list and name vector to rename list for loop
chromosome_contact = list()
list_names = c()

# Run loop to read contact data for each chromosome pair and produce contact value for each gene combination
# As minimum resolution is 5000 BP, if genes overlap contact bins then the average contact value across the whole gene length is taken
# Change Hi-C file in first line of loop for different data
for(i in 1:nrow(chromosomes)){
  
  tmp = straw_all(chromosomes$var1[i], chromosomes$var2[i], "AcolNg_V4.hic") # Produce data frame for chromosome combo
  tmp$chromosome_1_end = tmp[,1] + 4999 # Make end window for first chromosome
  tmp$chromosome_2_end = tmp[,2] + 4999 # Make end window for second chromosome
  tmp = tmp[,c(1,4,2,5,3)] # Reorder columns
  colnames(tmp) = c("chromosome_1_start","chromosome_1_end","chromosome_2_start","chromosome_2_end","count") # Rename columns
  
  
  chromosome_1_name = chromosomes[i,1] # Extract name of first chromosome
  chromosome_2_name = chromosomes[i,2] # Extract name of second chromosome
  tmp_genes_chromosome_1 = agam.genes %>% filter(chromosome == chromosome_1_name) # Subset genes for first chromosome
  tmp_genes_chromosome_1$gene_id = 1:nrow(tmp_genes_chromosome_1) # Assign gene ID for first chromosome
  colnames(tmp_genes_chromosome_1) = c("chromosome","annotation","type", "chromosome_1_start","chromosome_1_end","score","strand","frame","attribute","name","gene_id") # Rename columns
  tmp_genes_chromosome_2 = agam.genes %>% filter(chromosome == chromosome_2_name) # Subset genes for second chromosome
  tmp_genes_chromosome_2$gene_id = 1:nrow(tmp_genes_chromosome_2) # Assign gene ID for second chromosome
  colnames(tmp_genes_chromosome_2) = c("chromosome","annotation","type", "chromosome_2_start","chromosome_2_end","score","strand","frame","attribute","name","gene_id") # Rename columns
  
  tmp_genes_chromosome_1 = data.table(tmp_genes_chromosome_1) # Convert gene data frame to data table
  setkey(tmp_genes_chromosome_1, chromosome_1_start, chromosome_1_end) # Set key for first chromosome genes
  tmp_genes_chromosome_2 = data.table(tmp_genes_chromosome_2) # Convert gene data frame to data table
  setkey(tmp_genes_chromosome_2, chromosome_2_start, chromosome_2_end) # Set key for second chromosome genes
  
  
  # Assign window IDs
  tmp = tmp[order(tmp[,1]),] # Order by first chromosome
  tmp$chromosome_1_window_id = 1:nrow(tmp) # Assign window ID
  tmp = tmp[order(tmp[,3]),] # Order by second chromosome
  tmp$chromosome_2_window_id = 1:nrow(tmp) # Assign window ID
  
  
  # Assign first chromosome genes
  tmp_chromosome_1 = tmp[,c(1,2,5,6)] # Split into first chromosome only
  tmp_chromosome_1 = tmp_chromosome_1[order(tmp_chromosome_1[,1]),] # Order by first chromosome
  tmp_chromosome_1 = data.table(tmp_chromosome_1) # Convert contact data frame to data table
  setkey(tmp_chromosome_1, chromosome_1_start, chromosome_1_end) # Set key for first chromosome
  
  ov_chromosome_1 = foverlaps(tmp_chromosome_1, tmp_genes_chromosome_1, type = "any", which  = TRUE) # Check overlap of windows
  ov_chromosome_1 = merge(ov_chromosome_1, tmp_genes_chromosome_1[,c("gene_id", "name")], by.x = "yid", by.y = "gene_id", all = TRUE) # Merge to assign gene AGAP names
  tmp_chromosome_1 = merge(tmp_chromosome_1, ov_chromosome_1[,c("xid", "name")], by.x = "chromosome_1_window_id", by.y = "xid", all = TRUE) # Merge with contact data
  
  
  # Assign second chromosome genes
  tmp_chromosome_2 = tmp[,c(3,4,5,7)] # Split into second chromosome only
  tmp_chromosome_2 = tmp_chromosome_2[order(tmp_chromosome_2[,1]),] # Order by second chromosome
  tmp_chromosome_2 = data.table(tmp_chromosome_2) # Convert contact data frame to data table
  setkey(tmp_chromosome_2, chromosome_2_start, chromosome_2_end) # Set key for second chromosome
  
  ov_chromosome_2 = foverlaps(tmp_chromosome_2, tmp_genes_chromosome_2, type = "any", which  = TRUE) # Check overlap of windows
  ov_chromosome_2 = merge(ov_chromosome_2, tmp_genes_chromosome_2[,c("gene_id", "name")], by.x = "yid", by.y = "gene_id", all = TRUE) # Merge to assign gene AGAP names
  tmp_chromosome_2 = merge(tmp_chromosome_2, ov_chromosome_2[,c("xid", "name")], by.x = "chromosome_2_window_id", by.y = "xid", all = TRUE) # Merge with contact data
  
  
  # Merge both tables back together
  tmp = merge(tmp, tmp_chromosome_1, by = "chromosome_1_window_id")
  tmp = merge(tmp, tmp_chromosome_2, by = "chromosome_2_window_id")
  tmp = na.omit(tmp) # Remove NAs
  tmp = tmp[,c(10,11,15)] # Select only gene names and contact value
  tmp$id = cumsum(!duplicated(tmp[2:3])) # Make unique ID for each gene combination
  # Therefore any duplicated ID is where a gene is spanning across multiple contact bins
  # Group by unique ID and take contact average 
  tmp = tmp %>%
    group_by(id) %>%
    mutate(mean_count = mean(count.y)) %>%
    distinct(id, .keep_all = TRUE)
  
  tmp = tmp[,c(2,3,5)] # Select only gene names and mean contact value
  
  # Rename columns
  names(tmp)[names(tmp) == "name.x"] <- paste("chromosome_1", chromosome_1_name, "gene_name", sep = "_")
  names(tmp)[names(tmp) == "name.y"] <- paste("chromosome_2", chromosome_2_name, "gene_name", sep = "_")
  names(tmp)[names(tmp) == "mean_count"] <- "contact_count"
  
  list_names = c(list_names, paste(chromosomes$var1[i], chromosomes$var2[i], sep = "_")) # Rename list element
  chromosome_contact[[length(chromosome_contact) + 1]] <- tmp # Append data frame to list
  
  
  print(paste(i, "of", nrow(chromosomes), "chromosome pairings complete.", sep = " ")) # Print progress
}

# Rename list elements
names(chromosome_contact) <- list_names


# Append all chromosome combination tables in list to a dataframe
ColNames = c("chromosome_1_gene","chromosome_2_gene","contact_count") # Rename data frames in list so bind_rows works well
all_contact_df = lapply(chromosome_contact, setNames, ColNames) # Apply names
all_contact_df = bind_rows(all_contact_df) # Bind all data frames in list

# DO NOT USE IF LOOKING AT EXPECTED INSTEAD OF OBSERVED IN HI-C DATA
# Check which genes are missing
# Less genes are in contact matrix compared to distance matrix as not all observed genes have contact data from Hi-C file
unique_contact_genes = unique(all_contact_df$chromosome_1_gene) # Get all genes in contact data
unique_contact_genes = agam.gene.names[!agam.gene.names %in% unique_contact_genes] # Check which genes are missing vs Anopheles gambiae annotation file
missing_contact_genes = expand.grid(unique_contact_genes,unique_contact_genes) # Get all combinations of missing genes
missing_contact_genes$contact_count = 0 # Set these contact interactions as zero
colnames(missing_contact_genes) <- ColNames # Rename for binding to all_contact_df
all_contact_df = rbind(all_contact_df, missing_contact_genes) # Bind together

# Make matrix 
all_contact_df = all_contact_df[order(all_contact_df$chromosome_1_gene),] # Order by first chromosome
contact_matrix = pivot_wider(all_contact_df, names_from = chromosome_1_gene, values_from = contact_count)
contact_matrix = contact_matrix[order(contact_matrix$chromosome_2_gene),] # Order by second chromosome
matrix_row_names = contact_matrix$chromosome_2_gene # Row names from first column
contact_matrix = contact_matrix[,-1] # Remove first column
rownames(contact_matrix) <- matrix_row_names # Rename rows - ignore warning
contact_matrix = as.matrix(contact_matrix) # Set as matrix instead of tibble from pivot

# NA means there is no contact value between these genes
# Set NA as zero
contact_matrix[is.na(contact_matrix)] = 0

# Get a separate contact matrix for each dataset, with genes subsetted and ordered
contact.matrices <- c(lapply(cor.matrices, function(M) contact_matrix[rownames(M), colnames(M)]))

# Straighten out the matrices into lower triangle vectors 
contact.lower.tri <- lapply(contact.matrices, function(M) M[lower.tri(M)])

#### Gene Family Relatedness (InterPro approach) ####
agam_gene_family = read.csv("GeneByLocusTag_Summary.csv") # Table for each gene with InterPro IDs produced using VetorBase
agam_gene_family = agam_gene_family[,c(1,9)] # Select relevant columns
interpro_ids = ldply(strsplit(agam_gene_family$Interpro.ID, ";", fixed = TRUE), rbind) # Split InterPro IDs using ; delim
colnames(interpro_ids) = paste("interpro_id", colnames(interpro_ids), sep = "_")
agam_gene_family = cbind(agam_gene_family[,1], interpro_ids) # Bind
colnames(agam_gene_family)[1] = "gene_name" # Rename column

# Check missing genes and set InterPro ID to NA
unique_genes = unique(agam_gene_family$gene_name) # Get all genes in gene family data
missing_genes = (agam.gene.names[!agam.gene.names %in% unique_genes]) # Check which genes are missing vs Anopheles gambiae annotation file
missing_genes_df = data.frame(
  gene_name = missing_genes,
  interpro_id_1 = NA,
  interpro_id_2 = NA,
  interpro_id_3 = NA,
  interpro_id_4 = NA,
  interpro_id_5 = NA,
  interpro_id_6 = NA,
  interpro_id_7 = NA,
  interpro_id_8 = NA,
  interpro_id_9 = NA,
  interpro_id_10 = NA,
  interpro_id_11 = NA,
  interpro_id_12 = NA,
  interpro_id_13 = NA,
  interpro_id_14 = NA,
  interpro_id_15 = NA,
  interpro_id_16 = NA,
  interpro_id_17 = NA,
  interpro_id_18 = NA,
  interpro_id_19 = NA,
  interpro_id_20 = NA
) # Max InterPro IDs is 20

agam_gene_family = rbind(agam_gene_family, missing_genes_df)

# Set N/A character as NA class
agam_gene_family[agam_gene_family == "N/A" ] <- NA

# Calculate pairwise total InterPro domain matches between genes
agam_gene_family_matrix = crossprod(table(cbind(unlist(agam_gene_family[-1]), agam_gene_family[1])))

# Calculate pairwise total InterPro domain matches between genes as a percentage
# to avoid bias of some genes having more InterPro domains
# DON'T USE - this doesn't work well
agam_gene_family_sum = agam_gene_family
agam_gene_family_sum$total_domains = 20 - rowSums(is.na(agam_gene_family_sum))
agam_gene_family_sum = agam_gene_family_sum[,c(1,22)]
agam_gene_family_sum_matrix = outer(agam_gene_family_sum$total_domains, agam_gene_family_sum$total_domains, FUN = "+") # Total domains
colnames(agam_gene_family_sum_matrix) = agam_gene_family_sum$gene_name # Rename columns
rownames(agam_gene_family_sum_matrix) = agam_gene_family_sum$gene_name # Rename rows
agam_gene_family_percent_matrix = agam_gene_family_matrix/agam_gene_family_sum_matrix # Calculate percent domains shared

# Calculate pairwise yes/no any InterPro domain/family/homologous superfamily match between genes
# DON'T USE - this doesn't work well
agam_gene_family_matrix_binary = apply(agam_gene_family_matrix, 2, function(m) as.integer(m >= 1))
row.names(agam_gene_family_matrix_binary) <- colnames(agam_gene_family_matrix_binary)

# Get a separate InterPro domain match matrix for each dataset, with genes subsetted and ordered
gene.fam.matrices.interpro <- c(lapply(cor.matrices, function(M) agam_gene_family_matrix[rownames(M), colnames(M)]))

# Straighten out the matrices into lower triangle vectors 
gene.fam.interpro.lower.tri <- lapply(gene.fam.matrices.interpro, function(M) M[lower.tri(M)])


#### Gene Family Relatedness (Mash approach) ####

# This requires the use of mash (Ondov et al., 2016), using the Anopheles gambiae PEST reference fasta sequences
# Manually change k-mer size file due to computation limitations can't do for loop
# k-mer size = 12 is best

matrix_nrow = nrow(read.csv("output_12.csv")) - 1 # -1 as top row contains number of sequences
agam_gene_mash_matrix = read.table("output_12.csv", sep = "\t", col.names = 1:matrix_nrow, fill = TRUE) # Convert to table to process data
agam_gene_mash_matrix = agam_gene_mash_matrix[-1,] # Remove top row as contains number of sequences
agam_gene_mash_matrix = agam_gene_mash_matrix[-nrow(agam_gene_mash_matrix),] # Remove bottom row
row.names(agam_gene_mash_matrix) = agam_gene_mash_matrix[,1] # Set row names from first column
agam_gene_mash_matrix = agam_gene_mash_matrix[,-1] # Remove first column
agam_gene_mash_matrix = cbind(agam_gene_mash_matrix, data.frame(X13022 = rep(NA, 13022),
                                                            X13023 = rep(NA, 13022))) # Add two columns to end which are missing due to formatting
colnames(agam_gene_mash_matrix) = row.names(agam_gene_mash_matrix) # Set column names as row names
agam_gene_mash_matrix = as.matrix(agam_gene_mash_matrix) # Coerce to matrix
diag(agam_gene_mash_matrix) <- 0 # Set diagonal p-value as 0 as no phylogenetic distance between same gene
agam_gene_mash_matrix[nrow(agam_gene_mash_matrix), nrow(agam_gene_mash_matrix) - 1] = 1 # This cell missing, replace as 1

mash_missing_genes = agam.gene.names[!agam.gene.names %in% row.names(agam_gene_mash_matrix)] # Check which genes are missing vs Anopheles gambiae annotation file
mash_missing_genes_df = data.frame(matrix(1, nrow = length(mash_missing_genes), ncol = length(mash_missing_genes))) # Create dataframe
mash_missing_genes_df = cbind(as.data.frame(mash_missing_genes), mash_missing_genes_df) # Append missing genes into dataframe
rownames(mash_missing_genes_df) <- mash_missing_genes_df[,1] # Rename rows
mash_missing_genes_df = mash_missing_genes_df[,-1] # Remove first row
colnames(mash_missing_genes_df) <- rownames(mash_missing_genes_df) # Rename columns
mash_missing_genes_df = as.matrix(mash_missing_genes_df) # Coerce to matrix

agam_gene_mash_matrix = rbind.fill.matrix(agam_gene_mash_matrix, mash_missing_genes_df)
row.names(agam_gene_mash_matrix) <- colnames(agam_gene_mash_matrix)
agam_gene_mash_matrix = agam_gene_mash_matrix[,order(colnames(agam_gene_mash_matrix))]
agam_gene_mash_matrix = agam_gene_mash_matrix[order(row.names(agam_gene_mash_matrix)),]
agam_gene_mash_matrix[is.na(agam_gene_mash_matrix)] <- 1


# Get a separate mash phylogenetic distance matrix for each dataset, with genes subsetted and ordered
gene.fam.matrices.mash <- c(lapply(cor.matrices, function(M) agam_gene_mash_matrix[rownames(M), colnames(M)]))

# Straighten out the matrices into lower triangle vectors 
gene.fam.mash.lower.tri <- lapply(gene.fam.matrices.mash, function(M) M[lower.tri(M)])

# Combine data
bouake_mash = data.frame(expression = cor.lower.tri[["Bouake"]], gene_family = gene.fam.mash.lower.tri[["Bouake"]])
busia_mash = data.frame(expression = cor.lower.tri[["Busia"]], gene_family = gene.fam.mash.lower.tri[["Busia"]])
tiassale_mash = data.frame(expression = cor.lower.tri[["Tiassale"]], gene_family = gene.fam.mash.lower.tri[["Tiassale"]])

# Take average for each phylogenetic distance point as can't compute otherwise due to computational limitations
bouake_mash = bouake_mash %>%
  group_by(gene_family) %>%
  summarise(expression = mean(expression))

busia_mash = busia_mash %>%
  group_by(gene_family) %>%
  summarise(expression = mean(expression))

tiassale_mash = tiassale_mash %>%
  group_by(gene_family) %>%
  summarise(expression = mean(expression))

# Bouake
bouake_k_mer_glm = glm(formula = expression ~ gene_family, data = bouake_mash, family = Gamma(link = log))
summary(bouake_k_mer_glm)
with(summary(bouake_k_mer_glm), 1 - deviance/null.deviance) # Deviance explained by model (D-squared)

# Busia
busia_k_mer_glm = glm(formula = expression ~ gene_family, data = busia_mash, family = Gamma(link = log))
summary(busia_k_mer_glm)
with(summary(busia_k_mer_glm), 1 - deviance/null.deviance) # Deviance explained by model (D-squared)

# Tiassale
tiassale_k_mer_glm = glm(formula = expression ~ gene_family, data = tiassale_mash, family = Gamma(link = log))
summary(tiassale_k_mer_glm)
with(summary(tiassale_k_mer_glm), 1 - deviance/null.deviance) # Deviance explained by model (D-squared)

# Make note of D-squared for each model:
# 12: bouake = 0.06100829, busia = 0.04265849, tiassale = 0.08552134, all sig
# 13: bouake = 0.06911114, busia = 0.04158308, tiassale = 0.07652623, all sig
# 14: bouake = 0.05926568, busia = 0.03662493, tiassale = 0.06007802, all sig
# 15: bouake = 0.0525267, busia = 0.02999048, tiassale = 0.05178486, all sig
# 20: bouake = 0.02916834, busia = 0.01477157, tiassale = 0.008289245, all sig
# 25: bouake = 0.006902908, busia = 0.000004518568, tiassale = 0.0006590827, bouake ONLY sig
# 30: bouake = 0.009434342, busia = 0.0006091245, tiassale = 0.00001241168, bouake ONLY sig


# Make data frame for figure
k_mer_plot_df = data.frame(bouake = c(0.06100829, 0.06911114, 0.05926568, 0.0525267, 0.02916834, 0.006902908, 0.009434342),
                           busia = c(0.04265849, 0.04158308, 0.03662493, 0.02999048, 0.01477157, 0.000004518568, 0.0006091245),
                           tiassale = c(0.08552134, 0.07652623, 0.06007802, 0.05178486, 0.008289245, 0.0006590827, 0.00001241168)
                           )

k_mer_plot_df$mean = rowMeans(k_mer_plot_df)
k_mer_plot_df$k_mer_size = as.factor(c("12", "13", "14", "15", "20", "25", "30"))
k_mer_plot_df = k_mer_plot_df %>%
  gather(key = "variable", value = "value", -k_mer_size)

# Plot
ggplot(k_mer_plot_df, aes(x = k_mer_size, y = value, group = variable, color = variable)) +
  geom_line(size = 1) + geom_point() +
  labs(x = "K-mer size", y = "D-squared", color = "Dataset") +
  scale_color_manual(breaks = c("bouake", "busia", "tiassale", "mean"), values = c("#7CAE00","#00BFC4","#C77CFF","#F8766D"), 
                     labels = c("Bouaké", "Busia", "Tiassalé", "Mean"))



#### Combine distance, contact, and gene relatedness data ####
# Bouake
bouake_df = data.frame(
  contact = contact.lower.tri[["Bouake"]],
  distance = dist.lower.tri[["Bouake"]],
  expression = cor.lower.tri[["Bouake"]],
  gene_family_interpro = gene.fam.interpro.lower.tri[["Bouake"]],
  gene_family_mash = gene.fam.mash.lower.tri[["Bouake"]]
)
# Busia
busia_df = data.frame(
  contact = contact.lower.tri[["Busia"]],
  distance = dist.lower.tri[["Busia"]],
  expression = cor.lower.tri[["Busia"]],
  gene_family_interpro = gene.fam.interpro.lower.tri[["Busia"]],
  gene_family_mash = gene.fam.mash.lower.tri[["Busia"]]
)

# Tiassale
tiassale_df = data.frame(
  contact = contact.lower.tri[["Tiassale"]],
  distance = dist.lower.tri[["Tiassale"]],
  expression = cor.lower.tri[["Tiassale"]],
  gene_family_interpro = gene.fam.interpro.lower.tri[["Tiassale"]],
  gene_family_mash = gene.fam.mash.lower.tri[["Tiassale"]]
)

# Assign if genes pairs are on separate chromosomes or not
# This is done using the max distance value of agam.max.dist
bouake_df$chromosome[bouake_df$distance == agam.max.dist] <- "diff"
bouake_df$chromosome[bouake_df$distance != agam.max.dist] <- "same"

busia_df$chromosome[busia_df$distance == agam.max.dist] <- "diff"
busia_df$chromosome[busia_df$distance != agam.max.dist] <- "same"

tiassale_df$chromosome[tiassale_df$distance == agam.max.dist] <- "diff"
tiassale_df$chromosome[tiassale_df$distance != agam.max.dist] <- "same"

bouake_df_same = bouake_df[bouake_df$chromosome == "same",]
busia_df_same = busia_df[busia_df$chromosome == "same",]
tiassale_df_same = tiassale_df[tiassale_df$chromosome == "same",]

# Set Interpro domains >10 as 10
bouake_df_same$gene_family_interpro[bouake_df_same$gene_family_interpro > 10] <- 10
busia_df_same$gene_family_interpro[busia_df_same$gene_family_interpro > 10] <- 10
tiassale_df_same$gene_family_interpro[tiassale_df_same$gene_family_interpro > 10] <- 10

#### Figures ####

set.seed(123) # Set seed because sampling for graphs to ensure reproducibility

# Check distribution of response variable (expression) for model family
ggplot(bouake_df_same, aes(x = expression)) +
  geom_histogram(bins = 500)

ggplot(busia_df_same, aes(x = expression)) +
  geom_histogram(bins = 500)

ggplot(tiassale_df_same, aes(x = expression)) +
  geom_histogram(bins = 500)

# Poisson distribution, but not integers/counts, therefore gamma distribution

# Distance expression - use windowed average data (wa.cor and wa.log.dist) instead from expression_distance_correlation.r

# Bouake
tmp_bouake = data.frame(expression = wa.cor$Bouake, log10_distance = wa.log.dist$Bouake)
tmp_bouake = tmp_bouake[!tmp_bouake$log10_distance == log10(agam.max.dist + 1),]
ggplot(tmp_bouake, aes(x = log10_distance, y = expression)) +
  geom_point(size = 1, alpha = 1/5) +
  stat_smooth(method = "loess", color = "#F8766D", se = F) +
  labs(x = 'log10(Genomic Distance)', y ='Expression Correlation Coefficient', title = "Bouaké") # Plot

# Busia
tmp_busia = data.frame(expression = wa.cor$Busia, log10_distance = wa.log.dist$Busia)
tmp_busia = tmp_busia[!tmp_busia$log10_distance == log10(agam.max.dist + 1),]
ggplot(tmp_busia, aes(x = log10_distance, y = expression)) +
  geom_point(size = 1, alpha = 1/5) +
  stat_smooth(method = "loess", color = "#F8766D", se = F) +
  labs(x = 'log10(Genomic Distance)', y ='Expression Correlation Coefficient', title = "Busia") # Plot

# Tiassale
tmp_tiassale = data.frame(expression = wa.cor$Tiassale, log10_distance = wa.log.dist$Tiassale)
tmp_tiassale = tmp_tiassale[!tmp_tiassale$log10_distance == log10(agam.max.dist + 1),]
ggplot(tmp_tiassale, aes(x = log10_distance, y = expression)) +
  geom_point(size = 1, alpha = 1/5) +
  stat_smooth(method = "loess", color = "#F8766D", se = F) +
  labs(x = 'log10(Genomic Distance)', y ='Expression Correlation Coefficient', title = "Tiassalé") # Plot

# Contact expression

# Bouake
tmp_bouake = bouake_df_same[,c(1,3)] # Subset
tmp_bouake$log10_contact = log10(tmp_bouake$contact + 1) # Take log10

# Either take sample or group

tmp_bouake = sample_n(tmp_bouake, 1000000, replace = FALSE) # Take 1,000,000 sample for plotting due to computational limitations

tmp_bouake = tmp_bouake %>%
  group_by(log10_contact) %>%
  summarise(expression = mean(expression)) # Take average for each unique log10 contact frequency value

ggplot(tmp_bouake, aes(x = log10_contact, y = expression)) +
  geom_point(size = 1, alpha = 1/5) + 
  geom_smooth(se = F, color = "#F8766D", method = "loess") +
  scale_x_reverse() +
  labs(x = 'log10(Contact Frequency)', y ='Expression Correlation Coefficient', title = "Bouaké") # Plot
  
# Busia
tmp_busia = busia_df_same[,c(1,3)] # Subset
tmp_busia$log10_contact = log10(tmp_busia$contact + 1) # Take log10

# Either take sample or group

tmp_busia = sample_n(tmp_busia, 1000000, replace = FALSE) # Take 1,000,000 sample for plotting due to computational limitations

tmp_busia = tmp_busia %>%
  group_by(log10_contact) %>%
  summarise(expression = mean(expression)) # Take average for each unique log10 contact frequency value

ggplot(tmp_busia, aes(x = log10_contact, y = expression)) +
  geom_point(size = 1, alpha = 1/5) +
  geom_smooth(se = F, color = "#F8766D", method = "loess") +
  scale_x_reverse() +
  labs(x = 'log10(Contact Frequency)', y ='Expression Correlation Coefficient', title = "Busia") # Plot

# Tiassale
tmp_tiassale = tiassale_df_same[,c(1,3)] # Subset
tmp_tiassale$log10_contact = log10(tmp_tiassale$contact + 1) # Take log10

# Either take sample or group

tmp_tiassale = sample_n(tmp_tiassale, 1000000, replace = FALSE) # Take 1,000,000 sample for plotting due to computational limitations

tmp_tiassale = tmp_tiassale %>%
  group_by(log10_contact) %>%
  summarise(expression = mean(expression)) # Take average for each unique log10 contact frequency value

ggplot(tmp_tiassale, aes(x = log10_contact, y = expression)) +
  geom_point(size = 1, alpha = 1/5) +
  geom_smooth(se = F, color = "#F8766D", method = "loess") +
  scale_x_reverse() +
  labs(x = 'log10(Contact Frequency)', y ='Expression Correlation Coefficient', title = "Tiassalé") # Plot

# Distance contact

# Using Bouake as it has the most data points
# Doesn't matter which one is plotted as expression data is only variable changing
tmp_bouake = bouake_df_same[,1:2] # Subset
tmp_bouake$log10_distance = log10(tmp_bouake$distance + 1) # Take log10
tmp_bouake$log10_contact = log10(tmp_bouake$contact + 1) # Take log10
tmp_bouake = tmp_bouake[tmp_bouake$contact > 0,] # Remove values with no contact
tmp_bouake = tmp_bouake[tmp_bouake$log10_distance > 3,] # Truncate
tmp_bouake = sample_n(tmp_bouake, 1000000, replace = FALSE) # Take 1,000,000 sample for plotting due to computational limitations

# Plot
ggplot(tmp_bouake, aes(x = log10_distance, y = log10_contact)) +
  geom_point(size = 1, alpha = 1/5) +
  geom_smooth(se = F, color = "#F8766D", method = "loess") +
  labs(x = 'log10(Genomic Distance)', y ='log10(Contact Frequency)')

# Gene relatedness expression (InterPro)
# Merge to one data frame for long ggplot format for groups
tmp_bouake = bouake_df_same[,3:4]
tmp_bouake$dataset = "Bouaké"
tmp_busia = busia_df_same[,3:4]
tmp_busia$dataset = "Busia"
tmp_tiassale = tiassale_df_same[,3:4]
tmp_tiassale$dataset = "Tiassalé"

tmp_interpro_df = rbind(tmp_bouake, tmp_busia, tmp_tiassale)


# Plot
ggplot(tmp_interpro_df, aes(as.factor(gene_family_interpro), expression, fill = dataset)) +
  geom_boxplot(width = 0.5) +
  labs(x = 'Number of shared InterPro domains/families/homologous superfamilies', y ='Expression Correlation Coefficient', fill = "Dataset")

# Gene relatedness expression (mash)

# Bouake
tmp_bouake = bouake_df_same[,c(3,5)] # Subset

# Either take sample or group

tmp_bouake = sample_n(tmp_bouake, 1000000, replace = FALSE) # Take 1,000,000 sample for plotting due to computational limitations

tmp_bouake = tmp_bouake %>%
  group_by(gene_family_mash) %>%
  summarise(expression = mean(expression)) # Take average for each unique log10 contact frequency value

ggplot(tmp_bouake, aes(x = gene_family_mash, y = expression)) +
  geom_point(size = 1, alpha = 1/5) +
  geom_smooth(se = F, color = "#F8766D", method = "loess") +
  labs(x = 'Phylogenetic distance (p-value)', y ='Expression Correlation Coefficient', title = "Bouaké") # Plot

# Busia
tmp_busia = busia_df_same[,c(3,5)] # Subset

# Either take sample or group

tmp_busia = sample_n(tmp_busia, 1000000, replace = FALSE) # Take 1,000,000 sample for plotting due to computational limitations

tmp_busia = tmp_busia %>%
  group_by(gene_family_mash) %>%
  summarise(expression = mean(expression)) # Take average for each unique log10 contact frequency value

ggplot(tmp_busia, aes(x = gene_family_mash, y = expression)) +
  geom_point(size = 1, alpha = 1/5) +
  geom_smooth(se = F, color = "#F8766D", method = "loess") +
  labs(x = 'Phylogenetic distance (p-value)', y ='Expression Correlation Coefficient', title = "Busia") # Plot

# Tiassale
tmp_tiassale = tiassale_df_same[,c(3,5)] # Subset

# Either take sample or group

tmp_tiassale = sample_n(tmp_tiassale, 1000000, replace = FALSE) # Take 1,000,000 sample for plotting due to computational limitations

tmp_tiassale = tmp_tiassale %>%
  group_by(gene_family_mash) %>%
  summarise(expression = mean(expression)) # Take average for each unique log10 contact frequency value

ggplot(tmp_tiassale, aes(x = gene_family_mash, y = expression)) +
  geom_point(size = 1, alpha = 1/5) +
  geom_smooth(se = F, color = "#F8766D", method = "loess") +
  labs(x = 'Phylogenetic distance (p-value)', y ='Expression Correlation Coefficient', title = "Tiassalé") # Plot

# InterPro vs mash

# Percentage of data points with 0 InterPro matches
(sum(bouake_df_same$gene_family_interpro == 0) / nrow(bouake_df_same)) * 100 # 99.23%
(sum(busia_df_same$gene_family_interpro == 0) / nrow(busia_df_same)) * 100 # 99.22%
(sum(tiassale_df_same$gene_family_interpro == 0) / nrow(tiassale_df_same)) * 100 # 99.21%

# Percent of data points with p-value 1 Mash
(sum(bouake_df_same$gene_family_mash == 1) / nrow(bouake_df_same)) * 100 # 80.49%
(sum(busia_df_same$gene_family_mash == 1) / nrow(busia_df_same)) * 100 # 80.24%
(sum(tiassale_df_same$gene_family_mash == 1) / nrow(tiassale_df_same)) * 100 # 79.91%

# Using Bouake as it has the most data points
# Doesn't matter which one is plotted as expression data is only variable changing
tmp_bouake = bouake_df_same[,4:5] # Subset
tmp_bouake = tmp_bouake[tmp_bouake$gene_family_interpro > 0,] # Remove values with no interpro matches

tmp_bouake = tmp_bouake %>%
  group_by(gene_family_mash) %>%
  summarise(gene_family_interpro = mean(gene_family_interpro)) # Take average for each unique log10 contact frequency value

# Separate out maximum mash p-value so linear model can be computed without large value
tmp_bouake_0 = tmp_bouake[tmp_bouake$gene_family_mash == 1,] 
tmp_bouake = tmp_bouake[tmp_bouake$gene_family_mash < 1,]

# Linear model to extract R-squared
interpro_mash_lm = lm(gene_family_interpro ~ gene_family_mash, data = tmp_bouake)

# Plot
ggplot(tmp_bouake, aes(x = gene_family_mash, y = gene_family_interpro)) +
  geom_point(size = 1, alpha = 1/5) +
  geom_smooth(se = F, color = "#F8766D", method = "lm", size = 1) +
  scale_x_reverse() +
  geom_point(data = tmp_bouake_0, aes(gene_family_mash, gene_family_interpro)) +
  labs(x = 'Phylogenetic distance (p-value)', y ='Number of shared InterPro domains/families/homologous superfamilies')

#### Models ####

# Bouake
bouake_glm_contact = glm(formula = expression ~ contact, data = bouake_df_same, family = Gamma(link = log)) # Contact model
with(summary(bouake_glm_contact), 1 - deviance/null.deviance) # Deviance explained by model (D-squared)

bouake_glm_distance = glm(formula = expression ~ distance, data = bouake_df_same, family = Gamma(link = log)) # Distance model
with(summary(bouake_glm_distance), 1 - deviance/null.deviance) # Deviance explained by model (D-squared)

bouake_glm_gene_relatedness = glm(formula = expression ~ gene_family_mash, data = bouake_df_same, family = Gamma(link = log)) # Gene relatedness model
with(summary(bouake_glm_gene_relatedness), 1 - deviance/null.deviance) # Deviance explained by model (D-squared)

bouake_glm_final = glm(formula = expression ~ contact * distance * gene_family_mash, data = bouake_df_same, family = Gamma(link = log)) # Combined model
with(summary(bouake_glm_final), 1 - deviance/null.deviance) # Deviance explained by model (D-squared)

effectsize(bouake_glm_final, type = "omega") # Effect size

glm.diag.plots(bouake_glm_final) # Diagnostics

# Busia
busia_glm_contact = glm(formula = expression ~ contact, data = busia_df_same, family = Gamma(link = log)) # Contact model
with(summary(busia_glm_contact), 1 - deviance/null.deviance) # Deviance explained by model (D-squared)

busia_glm_distance = glm(formula = expression ~ distance, data = busia_df_same, family = Gamma(link = log)) # Distance model
with(summary(busia_glm_distance), 1 - deviance/null.deviance) # Deviance explained by model (D-squared)

busia_glm_gene_relatedness = glm(formula = expression ~ gene_family_mash, data = busia_df_same, family = Gamma(link = log)) # Gene relatedness model
with(summary(busia_glm_gene_relatedness), 1 - deviance/null.deviance) # Deviance explained by model (D-squared)

busia_glm_final = glm(formula = expression ~ contact * distance * gene_family_mash, data = busia_df_same, family = Gamma(link = log)) # Combined model
with(summary(busia_glm_final), 1 - deviance/null.deviance) # Deviance explained by model (D-squared)

effectsize(busia_glm_final, type = "omega") # Effect size

glm.diag.plots(busia_glm_final) # Diagnostics


# Tiassale
tiassale_glm_contact = glm(formula = expression ~ contact, data = tiassale_df_same, family = Gamma(link = log)) # Contact model
with(summary(tiassale_glm_contact), 1 - deviance/null.deviance) # Deviance explained by model (D-squared)

tiassale_glm_distance = glm(formula = expression ~ distance, data = tiassale_df_same, family = Gamma(link = log)) # Distance model
with(summary(tiassale_glm_distance), 1 - deviance/null.deviance) # Deviance explained by model (D-squared)

tiassale_glm_gene_relatedness = glm(formula = expression ~ gene_family_mash, data = tiassale_df_same, family = Gamma(link = log)) # Gene relatedness model
with(summary(tiassale_glm_gene_relatedness), 1 - deviance/null.deviance) # Deviance explained by model (D-squared)

tiassale_glm_final = glm(formula = expression ~ contact * distance * gene_family_mash, data = tiassale_df_same, family = Gamma(link = log)) # Combined model
with(summary(tiassale_glm_final), 1 - deviance/null.deviance) # Deviance explained by model (D-squared)

effectsize(tiassale_glm_final, type = "omega") # Effect size

glm.diag.plots(tiassale_glm_final) # Diagnostics
