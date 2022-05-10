library(strawr)
library(ggplot2)
library(dplyr)
library(data.table)
library(tidyr)

setwd("D:/OneDrive - Lancaster University/PR2/data/") # Change as appropriate

# Hi-C data from Lukyanchikova et al. (20202) - check resolution and length of chromosomes
# Anopheles coluzzi 
strawr::readHicBpResolutions("AcolNg_V4.hic") # 5000 is minimum resolution
strawr::readHicChroms("AcolNg_V4.hic") # chromosome lengths


# Anopheles gambiae (PEST) gene annotation
agam.gff <- read.table(gzfile('Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.12.gff3.gz', 'r'), sep = '\t', quote = '') # Read file
colnames(agam.gff) <- c('chromosome', 'annotation', 'type', 'start', 'end', 'score', 'strand', 'frame', 'attribute') # Rename columns
agam.genes <- subset(agam.gff, grepl('gene', type)) # Subset just genes
agam.gene.names <- sub('^ID=', '', sub(';.*', '', agam.genes$attribute)) # Extract AGAP gene names
rownames(agam.genes) <- agam.gene.names # Assign gene names to row names
agam.genes$name = agam.gene.names # Assign gene names to name column


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

# Create simplified straw function (KR normalisation used)
straw_all = function(chrom1, chrom2, data){
  output = strawr::straw(norm = "KR", fname = data, 
                       chr1loc = chrom1, chr2loc = chrom2, 
                       unit = "BP", binsize = 5000, matrix = "observed")
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

# Check which genes are missing
# Less genes are in contact matrix compared to distance matrix as not all genes have contact data from Hi-C file
unique_contact_genes = unique(all_contact_df$chromosome_1_gene) # Get all genes in contact data
unique_contact_genes = agam.gene.names[!agam.gene.names %in% unique_contact_genes] # Check which genes are missing vs Anopheles gambiae annotation file
missing_contact_genes = expand.grid(unique_contact_genes,unique_contact_genes) # Get all combinations of missing genes
missing_contact_genes$contact_count = 0 # Set these contact interactions as zero
colnames(missing_contact_genes) <- ColNames # Rename for binding to all_contact_df
all_contact_df = rbind(all_contact_df, missing_contact_genes) # Bind together
all_contact_df = all_contact_df[order(all_contact_df$chromosome_1_gene),] # Order by first chromosome

# Make matrix 
contact_matrix = pivot_wider(all_contact_df, names_from = chromosome_1_gene, values_from = contact_count)
contact_matrix = contact_matrix[order(contact_matrix$chromosome_2_gene),] # Order by second chromosome
matrix_row_names = contact_matrix$chromosome_2_gene # Row names from first column
contact_matrix = contact_matrix[,-1] # Remove first column
rownames(contact_matrix) <- matrix_row_names # Rename rows - ignore warning
contact_matrix = as.matrix(contact_matrix) # Set as matrix instead of tibble from pivot

# NA means there is no contact value between these genes
# Set NA as zero
contact_matrix[is.na(contact_matrix)] = 0

# We now want to create a matrix of distances between genes. For any two genes that share a scaffold, we take 
# the absolute value of the difference between their start points. 
# Create a matrix where all distances are assumed to be the maximum value (double the size of the largest 
# observed gene position)
agam.genes.by.chrom <- split(agam.genes, as.character(agam.genes$chromosome)) # Split by chromosome
agam.max.dist <- 2*max(agam.gff$end)
agam.dist.matrix <- matrix(agam.max.dist, nrow(agam.genes), nrow(agam.genes), dimnames = list(sort(agam.gene.names), sort(agam.gene.names))) # Ordered
agam.within.chrom.distances <- lapply(agam.genes.by.chrom, function(X) apply(X[, 'start', drop = F], 1, function(x) abs(x['start'] - X$start)))
for (wcd in agam.within.chrom.distances){
  agam.dist.matrix[colnames(wcd), colnames(wcd)] <- wcd
}

# Straighten out contact and distance matrices into lower triangle vectors
contact_vec = contact_matrix[lower.tri(contact_matrix)]
distance_vec = agam.dist.matrix[lower.tri(agam.dist.matrix)]

# Take log for visualisation
contact_vec_log = log10(contact_vec + 1)
distance_vec_log = log10(distance_vec + 1)

# Combine distance and contact data
contact_distance_df = data.frame(contact_vec_log, distance_vec_log)

# Set NA as zero again
contact_distance_df[is.na(contact_distance_df)] = 0

# Subset
plot_df = contact_distance_df # Randomly sample fewer than total observations if required, otherwise leave as is
plot_df = plot_df[which(plot_df$distance_vec < 8),] # Remove < 8 log10(distance)
plot_df = plot_df[which(plot_df$contact_vec > 0),] # Remove 0 values for ease of interpretation of plot

# Plot to check data
ggplot(plot_df, aes(x = distance_vec_log, y = contact_vec_log)) +
  geom_point(size = 0.75, alpha = 1/10) +
  xlab("Log10 (Linear Genomic Distance)") + ylab("Log10 (Chromosomal Contact Value)")
