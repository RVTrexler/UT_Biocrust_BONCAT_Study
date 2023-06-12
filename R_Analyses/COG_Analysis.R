# BONCAT Biocrust Project - Analysis of COG data from IMG metagenomes
# Version 6.0
# Ryan V. Trexler
# 12/01/2022

# Note: COG Annotation data downloaded here: https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/


# Clear environment
remove(list =ls())

# Load libraries
library(tidyverse)
library(vegan)
library(ggpubr)


setwd("~/Documents/PSU_Documents/Couradeau_Lab/Biocrust_Boncat_MG/metagenome_analysis/IMG_Data")


#### 1.1. Import metadata ####

# Import metadata
metadata <- read.table(file = "Metadata_key.txt", sep = '\t', header = T)
# Note: metadata_key.txt was modified from previous version to make an additional column for a Descriptive_Sample_ID column

# The following samples were not available from IMG: SYTO_8 , SYTO_18
# Filter these from the metadata file
metadata <- metadata %>% filter(Sample != "SYTO_8") %>% filter(Sample != "SYTO_18")

#### 1.2. Import cog data ####

cog_dir='~/Documents/PSU_Documents/Couradeau_Lab/Biocrust_Boncat_MG/metagenome_analysis/IMG_Data/COG_Data/'

# Initialize compiled table data frame
cog.table <- data.frame()
infile <- NULL

for (infile in list.files(cog_dir, pattern = ".txt"))
{
	raw.inp <- data.frame()
	
	# Read in sample data
	raw.inp <- read.table(file = paste(cog_dir,infile,sep=""), sep = '\t', header = T, quote = "", fill = T, comment.char = "")

	# Add sample name column to data
	colnames(raw.inp) <- c("COG_ID","COG_Name", str_sub(infile, end=-5))
	
	# Add to final cog table
	if(dim(cog.table)[1] == 0){
	  cog.table <- raw.inp
	}else{
	  cog.table <- full_join(cog.table, raw.inp, by = c("COG_ID","COG_Name"))
	}
}

# Replace NAs with zeros
cog.table <- cog.table %>% mutate(across(where(is.integer), replace_na, 0))


## Merge cog.table with COG annotation data
# Import annotations
cog.annotations <- read_tsv(file = "COG_ID_Annotation.txt") %>% mutate(Pathway = replace_na(.$Pathway,"Not_Assigned"))
cog.categories <- read_tsv(file = "COG_Categories_Expanded.txt")

# Add annotations based on COG_ID
cog.table <- cog.table %>% left_join(cog.annotations, by = "COG_ID") 
# Add COG Category annotations based on COG_Category_ID
cog.table <- cog.table %>% left_join(cog.categories, by = "COG_Category_ID") %>% relocate(COG_ID, COG_Name, COG_Category_ID, COG_Category_1_Annotation, COG_Category_2_Annotation, COG_Category_3_Annotation, COG_Category_4_Annotation, Annotation, Gene, Pathway)

#write_excel_csv(cog.table,"manuscript_data_output/COG_Table.csv")


#### 2.0 NMDS ####

cog.table.nmds <- metaMDS(cog.table.pcoa,k=2) # Perform NMDS on Bray Curtis distances

cog.nmds.scores <- as.data.frame(scores(cog.table.nmds, "sites")) #Using the scores function from vegan to extract the site scores and convert to a data.frame
cog.nmds.scores <- cog.nmds.scores %>% rownames_to_column(., var = "Sample") %>% left_join(metadata, by = "Sample")

# Get NMDS Group (Treatment) centroids
centroids <- aggregate(cbind(NMDS1, NMDS2) ~ Treatment, data = cog.nmds.scores, FUN = mean) %>% rename(Cent_NMDS1 = NMDS1) %>% rename(Cent_NMDS2 = NMDS2)
# Get NMDS dispersions (distances from group centroid)
distances <- cog.nmds.scores %>% left_join(centroids, by = "Treatment") %>% mutate(Distance = sqrt(((NMDS1-Cent_NMDS1)^2 + (NMDS2-Cent_NMDS2)^2)))

# Visualize by Treatment
#dev.new()
ggplot(cog.nmds.scores, aes(x=NMDS1,y=NMDS2)) + geom_point(aes(shape=Timepoint, color=Fraction,fill=Fraction),size=4) + labs(title="Bray-Curtis NMDS of COGs") + theme(plot.title = element_text(hjust = 0.5)) + theme(panel.background = element_rect(fill = "white", colour = "grey50"), legend.position = c(0.8,0.7), legend.background = element_rect(color = "black", fill = "grey90", size = 0.5, linetype = "solid")) + geom_text(x=0.25, y=0.0, label= paste("Stress: ", as.character(round(cog.table.nmds$stress,4))))

# Permanova of bray-curtis distances
cog.table.bray.permanova <- adonis(cog.table.bray ~ Fraction*Timepoint, dat=metadata, perm=1e3)
cog.table.bray.permanova$aov.tab


# Beta-dispersion of NMDS distances
# Sort it so that treatment is ordered correctly
distances$Treatment <- factor(distances$Treatment, levels = c("Dry", "4h_SYTO_plus", "4h_BC_minus", "4h_BC_plus", "21h_SYTO_plus", "21h_BC_minus", "21h_BC_plus"))

# Plot boxplots of distance to group centroid
ggplot(distances, aes(x=Treatment, y=Distance, fill=Fraction)) + geom_boxplot() + geom_point(position=position_dodge(width=0.75), aes(group=Fraction))

# ANOVA of beta-dispersion (distances from group centroid)
distances.anova <- aov(Distance ~ Treatment, data = distances)
summary(distances.anova)
TukeyHSD(distances.anova)

distances.sum <-  distances %>% group_by(Treatment) %>% summarize(n = n(), mean = mean(Distance), median = median(Distance), sd = sd(Distance), se = sd/sqrt(n()))

#write_excel_csv(distances.sum,"manuscript_data_output/COG_BetaDisp_summary.csv")

#### 3.0 Plot Pathway abundances ####

## Pathway
cog.table.pathway <- cog.table %>% pivot_longer(-c("COG_ID", "COG_Name", "COG_Category_ID","COG_Category_1_Annotation", "COG_Category_2_Annotation", "COG_Category_3_Annotation", "COG_Category_4_Annotation", "Annotation", "Gene", "Pathway"), names_to = "Sample", values_to = "COG_Count") %>% group_by(Sample, Pathway) %>% summarize(COG_Count_Pathway = sum(COG_Count))
cog.table.pathway <- cog.table.pathway %>% group_by(Sample) %>% mutate(COG_Count_Pathway_RA = COG_Count_Pathway/sum(COG_Count_Pathway))
ggplot(cog.table.pathway, aes(x=Sample,y=COG_Count_Pathway_RA, fill = Pathway)) + geom_bar(stat='identity') + theme_bw()

## Category
cog.table.category <- cog.table %>% pivot_longer(-c("COG_ID", "COG_Name", "COG_Category_ID","COG_Category_1_Annotation", "COG_Category_2_Annotation", "COG_Category_3_Annotation", "COG_Category_4_Annotation", "Annotation", "Gene", "Pathway"), names_to = "Sample", values_to = "COG_Count") %>% group_by(Sample, COG_Category_1_Annotation) %>% summarize(COG_Count_Category = sum(COG_Count))
# Calcuate category relative abundances
cog.table.category <- cog.table.category %>% group_by(Sample) %>% mutate(COG_Count_Category_RA = COG_Count_Category/sum(COG_Count_Category))
ggplot(cog.table.category, aes(x=Sample,y=COG_Count_Category_RA, fill = COG_Category_1_Annotation)) + geom_bar(stat='identity') + theme_bw()


#### 4.0 Differential abundance of COGS (at COG_ID level) between active and inactive fractions: DeSeq2 ####

# Transpose the cog table, format properly
cog.table.diffabund <- t(cog.table %>% select(-c("COG_Name", "COG_Category_ID","COG_Category_1_Annotation", "COG_Category_2_Annotation", "COG_Category_3_Annotation", "COG_Category_4_Annotation", "Annotation", "Gene", "Pathway")))
colnames(cog.table.diffabund) <- cog.table.diffabund[1,1:ncol(cog.table.diffabund)]
cog.table.diffabund <- cog.table.diffabund[-1,]
cog.table.diffabund <- rownames_to_column(as.data.frame(cog.table.diffabund)) %>% rename(Sample = "rowname")

# Get diffabund table for 4hrs
cog.table.diffabund.4h <- cog.table.diffabund %>% left_join(metadata, by = 'Sample') %>% filter(Timepoint == "4h") %>% filter(Fraction == "BC_plus" | Fraction == "BC_minus") %>% select(-c("Fraction", "Desc_SampleID", "Bio_Rep",	"Tech_Rep", "Treatment", "Timepoint")) %>% relocate(Sample) %>% as.tibble() %>% t()
colnames(cog.table.diffabund.4h) <- cog.table.diffabund.4h[1,]
cog.table.diffabund.4h <- cog.table.diffabund.4h[2:nrow(cog.table.diffabund.4h),]
cog.table.diffabund.4h <- cog.table.diffabund.4h %>% as.data.frame() %>% rownames_to_column() %>% mutate(COG_ID = rowname) %>% select(-rowname) %>% relocate(COG_ID)
cog.table.diffabund.4h <- cog.table.diffabund.4h %>% mutate_at(vars(-("COG_ID")), as.integer) %>% mutate_at("COG_ID", as.factor)

# Get diffabund table for 21hrs
cog.table.diffabund.21h <- cog.table.diffabund %>% left_join(metadata, by = 'Sample') %>% filter(Timepoint == "21h") %>% filter(Fraction == "BC_plus" | Fraction == "BC_minus") %>% select(-c("Fraction", "Desc_SampleID", "Bio_Rep",	"Tech_Rep", "Treatment", "Timepoint")) %>% relocate(Sample) %>% as.tibble() %>% t()
colnames(cog.table.diffabund.21h) <- cog.table.diffabund.21h[1,]
cog.table.diffabund.21h <- cog.table.diffabund.21h[2:nrow(cog.table.diffabund.21h),]
cog.table.diffabund.21h <- cog.table.diffabund.21h %>% as.data.frame() %>% rownames_to_column() %>% mutate(COG_ID = rowname) %>% select(-rowname) %>% relocate(COG_ID)
cog.table.diffabund.21h <- cog.table.diffabund.21h %>% mutate_at(vars(-("COG_ID")), as.integer) %>% mutate_at("COG_ID", as.factor)


# Get metadata table for deseq2
deseq.md.4h <- metadata %>% filter(Timepoint=="4h") %>% filter(Fraction != "SYTO_plus")
deseq.md.21h <- metadata %>% filter(Timepoint=="21h") %>% filter(Fraction != "SYTO_plus")


#library(DESeq2)
# Note: loading DESeq2 will mess up some tidyverse functions!
# Create DESeq2 dataset
dds.4h <- DESeqDataSetFromMatrix(countData=cog.table.diffabund.4h, colData=deseq.md.4h, design = ~ Fraction, tidy=T)
dds.21h <- DESeqDataSetFromMatrix(countData=cog.table.diffabund.21h, colData=deseq.md.21h, design = ~ Fraction, tidy=T)


# Run DESeq function
dds.4h <- DESeq(dds.4h)
dds.21h <- DESeq(dds.21h)

# Contrast Fractions
res.4h.fraction <- results(dds.4h, contrast=c("Fraction", "BC_plus", "BC_minus"), tidy = TRUE)
res.21h.fraction <- results(dds.21h, contrast=c("Fraction", "BC_plus", "BC_minus"), tidy = TRUE)


# Look at results
#res.SYTO4h <- results(dds)
summary(res.4h.fraction)
res.4h.fraction <- res.4h.fraction[order(res.4h.fraction$padj),]
summary(res.21h.fraction)
res.21h.fraction <- res.21h.fraction[order(res.21h.fraction$padj),]


# Fold change bar plot
res.4h.boncat.plot <- res.4h.fraction %>% filter(padj <= 0.05) %>% mutate(COG_ID = row) %>% select(-row) %>% left_join(., cog.annotations, by = 'COG_ID') %>% left_join(., cog.categories, by = 'COG_Category_ID')
res.21h.boncat.plot <- res.21h.fraction %>% filter(padj <= 0.05) %>% mutate(COG_ID = row) %>% select(-row) %>% left_join(., cog.annotations, by = 'COG_ID') %>% left_join(., cog.categories, by = 'COG_Category_ID')

res.4h.boncat.plot <- res.4h.boncat.plot %>% unite(ID_Label, COG_ID, Gene, sep = "_", remove = FALSE)
res.21h.boncat.plot <- res.21h.boncat.plot %>% unite(ID_Label, COG_ID, Gene, sep = "_", remove = FALSE)


## 4 hrs
#dev.new()
p1 <- ggplot(res.4h.boncat.plot, aes(x = reorder(COG_ID, log2FoldChange), y = log2FoldChange, fill = Pathway)) + geom_bar(stat='identity') + theme_bw() + theme(axis.title.y = element_blank()) + coord_flip() #+ theme(legend.position = "none")
p2 <- ggplot(res.4h.boncat.plot, aes(x = reorder(COG_ID, log2FoldChange), y = padj)) + geom_point() + theme_bw() + theme(axis.text.y=element_blank(), axis.title.y = element_blank()) + coord_flip() + theme(legend.position = "none")
leg <- get_legend(p1)
ggarrange(p1, p2, as_ggplot(leg), ncol=3, nrow=1, legend="none", widths = c(1, 0.25, .5))


# Just plot it without pvalues and with legend
# Generate custom color pallete
pal <- c("#2bd0c0", "#3643e4", "#f972f2", "#ccc59b", "#8d8c22", "#b0f9c5", "#966223", "#ebfd9a", "#f9dcf8", "#637b5e", "#b51dbf", "#573ec2", "#ef6703", "#3ea128", "#c36b56", "#0b1d6c", "#ce144a")
ggplot(res.4h.boncat.plot %>% filter(log2FoldChange >= 2.0 | log2FoldChange <= -1.0), aes(x = reorder(COG_ID, log2FoldChange), y = log2FoldChange, fill = Pathway)) + geom_bar(stat='identity') + theme_bw() + theme(axis.title.y = element_blank()) + coord_flip() #+ theme(legend.position = "none")
ggplot(res.4h.boncat.plot %>% filter(log2FoldChange >= 1.5 | log2FoldChange <= -1.0) %>% filter(COG_Category_1_Annotation != "Function unknown") %>% filter(COG_Category_1_Annotation != "General function prediction only"), aes(x = reorder(ID_Label, log2FoldChange), y = log2FoldChange, fill = COG_Category_1_Annotation)) + geom_bar(stat='identity') + scale_fill_manual(values = pal) + theme_bw() + theme(axis.title.y = element_blank()) + coord_flip() #+ theme(legend.position = "none")

#write_csv(res.4h.boncat.plot,"~/Documents/PSU_Documents/Couradeau_Lab/Biocrust_Boncat_MG/metagenome_analysis/IMG_Data/deseq_cog_id_4hr.csv")

## 21 hrs
p1 <- ggplot(res.21h.boncat.plot, aes(x = reorder(COG_ID, log2FoldChange), y = log2FoldChange, fill = Pathway)) + geom_bar(stat='identity') + theme_bw() + theme(axis.title.y = element_blank()) + coord_flip() #+ theme(legend.position = "none")
p2 <- ggplot(res.21h.boncat.plot, aes(x = reorder(COG_ID, log2FoldChange), y = padj)) + geom_point() + theme_bw() + theme(axis.text.y=element_blank(), axis.title.y = element_blank()) + coord_flip() + theme(legend.position = "none")
leg <- get_legend(p1)
ggarrange(p1, p2, as_ggplot(leg), ncol=3, nrow=1, legend="none", widths = c(1, 0.25, .5))

# Just plot it without pvalues and with legend

# Generate custom color pallete
pal <- c("#3643e4", "#ccc59b", "#ebfd9a", "#e64323", "#b51dbf", "#ce144a")
ggplot(res.21h.boncat.plot, aes(x = reorder(COG_ID, log2FoldChange), y = log2FoldChange, fill = Pathway)) + geom_bar(stat='identity') + theme_bw() + theme(axis.title.y = element_blank()) + coord_flip() #+ theme(legend.position = "none")
ggplot(res.21h.boncat.plot, aes(x = reorder(ID_Label, log2FoldChange), y = log2FoldChange, fill = COG_Category_1_Annotation)) + geom_bar(stat='identity') + scale_fill_manual(values = pal) + theme_bw() + theme(axis.title.y = element_blank()) + coord_flip() #+ theme(legend.position = "none")

#write_excel_csv(res.4h.boncat.plot,"manuscript_data_output/DESeq2_COG_ID_4h.csv")
#write_excel_csv(res.21h.boncat.plot,"manuscript_data_output/DESeq2_COG_ID_21h.csv")


#### 5.0 Differential abundance of COG Pathways (at COG Pathway level) between active and inactive fractions: DeSeq2 ####

cog.table.pathway <- cog.table.pathway %>% select(-COG_Count_Pathway_RA) %>% pivot_wider(names_from = Sample, values_from = COG_Count_Pathway)


# Transpose the cog table, format properly
cog.table.pathway.diffabund <- t(cog.table.pathway) %>% as.data.frame() #%>% select(-c("COG_Name", "COG_Category_ID", "Annotation", "Gene", "Pathway")))
colnames(cog.table.pathway.diffabund) <- cog.table.pathway.diffabund[1,1:ncol(cog.table.pathway.diffabund)]
cog.table.pathway.diffabund <- cog.table.pathway.diffabund[-1,]
cog.table.pathway.diffabund <- add_rownames(cog.table.pathway.diffabund, var = "rowname") %>% rename(Sample = "rowname")

# Get diffabund table for 4hrs
cog.table.pathway.diffabund.4h <- cog.table.pathway.diffabund %>% left_join(metadata, by = 'Sample') %>% filter(Timepoint == "4h") %>% filter(Fraction == "BC_plus" | Fraction == "BC_minus") %>% select(-c("Fraction", "Desc_SampleID", "Bio_Rep",	"Tech_Rep", "Treatment", "Timepoint")) %>% relocate(Sample) %>% as.tibble() %>% t()
colnames(cog.table.pathway.diffabund.4h) <- cog.table.pathway.diffabund.4h[1,]
cog.table.pathway.diffabund.4h <- cog.table.pathway.diffabund.4h[2:nrow(cog.table.pathway.diffabund.4h),]
cog.table.pathway.diffabund.4h <- cog.table.pathway.diffabund.4h %>% as.data.frame() %>% rownames_to_column() %>% mutate(COG_ID = rowname) %>% select(-rowname) %>% relocate(COG_ID)
cog.table.pathway.diffabund.4h <- cog.table.pathway.diffabund.4h %>% mutate_at(vars(-("COG_ID")), as.integer) %>% mutate_at("COG_ID", as.factor)
cog.table.pathway.diffabund.4h <- cog.table.pathway.diffabund.4h[1:nrow(cog.table.pathway.diffabund.4h)-1,]


# Get diffabund table for 21hrs
cog.table.pathway.diffabund.21h <- cog.table.pathway.diffabund %>% left_join(metadata, by = 'Sample') %>% filter(Timepoint == "21h") %>% filter(Fraction == "BC_plus" | Fraction == "BC_minus") %>% select(-c("Fraction", "Desc_SampleID", "Bio_Rep",	"Tech_Rep", "Treatment", "Timepoint")) %>% relocate(Sample) %>% as.tibble() %>% t()
colnames(cog.table.pathway.diffabund.21h) <- cog.table.pathway.diffabund.21h[1,]
cog.table.pathway.diffabund.21h <- cog.table.pathway.diffabund.21h[2:nrow(cog.table.pathway.diffabund.21h),]
cog.table.pathway.diffabund.21h <- cog.table.pathway.diffabund.21h %>% as.data.frame() %>% rownames_to_column() %>% mutate(COG_ID = rowname) %>% select(-rowname) %>% relocate(COG_ID)
cog.table.pathway.diffabund.21h <- cog.table.pathway.diffabund.21h %>% mutate_at(vars(-("COG_ID")), as.integer) %>% mutate_at("COG_ID", as.factor)


# Get metadata table for deseq2
deseq.md.4h <- metadata %>% filter(Timepoint=="4h") %>% filter(Fraction != "SYTO_plus")
deseq.md.21h <- metadata %>% filter(Timepoint=="21h") %>% filter(Fraction != "SYTO_plus")


#library(DESeq2)
# Note: loading DESeq2 will mess up some tidyverse functions!
# Create DESeq2 dataset
dds.pathway.4h <- DESeqDataSetFromMatrix(countData=cog.table.pathway.diffabund.4h, colData=deseq.md.4h, design = ~ Fraction, tidy=T)
dds.pathway.21h <- DESeqDataSetFromMatrix(countData=cog.table.pathway.diffabund.21h, colData=deseq.md.21h, design = ~ Fraction, tidy=T)


# Run DESeq function
dds.pathway.4h <- DESeq(dds.pathway.4h)
dds.pathway.21h <- DESeq(dds.pathway.21h)

# Contrast Fractions
res.pathway.4h.fraction <- results(dds.pathway.4h, contrast=c("Fraction", "BC_plus", "BC_minus"), tidy = TRUE)
res.pathway.21h.fraction <- results(dds.pathway.21h, contrast=c("Fraction", "BC_plus", "BC_minus"), tidy = TRUE)

# Look at results
#res.SYTO4h <- results(dds)
summary(res.pathway.4h.fraction)
res.pathway.4h.fraction <- res.pathway.4h.fraction[order(res.pathway.4h.fraction$padj),]
summary(res.pathway.21h.fraction)
res.pathway.21h.fraction <- res.pathway.21h.fraction[order(res.pathway.21h.fraction$padj),]


# Fold change bar plot
res.pathway.4h.boncat.plot <- res.pathway.4h.fraction %>% filter(padj <= 0.05) %>% mutate(COG_Pathway = row) %>% select(-row) #%>% left_join(., cog.annotations, by = 'COG_ID') %>% left_join(., cog.categories, by = 'COG_Category_ID')
res.pathway.21h.boncat.plot <- res.pathway.21h.fraction %>% filter(padj <= 0.05) %>% mutate(COG_Pathway = row) %>% select(-row) #%>% left_join(., cog.annotations, by = 'COG_ID') %>% left_join(., cog.categories, by = 'COG_Category_ID')

#dev.new()
p1 <- ggplot(res.pathway.4h.boncat.plot, aes(x = reorder(COG_Pathway, log2FoldChange), y = log2FoldChange)) + geom_bar(stat='identity') + theme_bw() + theme(axis.title.y = element_blank()) + coord_flip() + theme(legend.position = "none")
p2 <- ggplot(res.pathway.4h.boncat.plot, aes(x = reorder(COG_Pathway, log2FoldChange), y = padj)) + geom_point() + theme_bw() + theme(axis.text.y=element_blank(), axis.title.y = element_blank()) + coord_flip() + theme(legend.position = "none")
ggarrange(p1, p2, ncol=2, nrow=1, legend="none", widths = c(1, 0.25))

p1 <- ggplot(res.pathway.21h.boncat.plot, aes(x = reorder(COG_Pathway, log2FoldChange), y = log2FoldChange)) + geom_bar(stat='identity') + theme_bw() + theme(axis.title.y = element_blank()) + coord_flip() + theme(legend.position = "none")
p2 <- ggplot(res.pathway.21h.boncat.plot, aes(x = reorder(COG_Pathway, log2FoldChange), y = padj)) + geom_point() + theme_bw() + theme(axis.text.y=element_blank(), axis.title.y = element_blank()) + coord_flip() + theme(legend.position = "none")
ggarrange(p1, p2, ncol=2, nrow=1, legend="none", widths = c(1, 0.25))

#write_excel_csv(res.pathway.4h.boncat.plot,"manuscript_data_output/DESeq2_COG_Pathway_4h.csv")
#write_excel_csv(res.pathway.21h.boncat.plot,"manuscript_data_output/DESeq2_COG_Pathway_21h.csv")


#### 5.0 Differential abundance of COG Categories (at COG Category level) between active and inactive fractions: DeSeq2 ####

cog.table.category <- cog.table.category %>% select(-COG_Count_Category_RA) %>% pivot_wider(names_from = Sample, values_from = COG_Count_Category)

# Transpose the cog table, format properly
cog.table.category.diffabund <- t(cog.table.category) %>% as.data.frame() #%>% select(-c("COG_Name", "COG_Category_ID", "Annotation", "Gene", "category")))
colnames(cog.table.category.diffabund) <- cog.table.category.diffabund[1,1:ncol(cog.table.category.diffabund)]
cog.table.category.diffabund <- cog.table.category.diffabund[-1,]
cog.table.category.diffabund <- add_rownames(cog.table.category.diffabund, var = "rowname") %>% rename(Sample = "rowname")


# Get diffabund table for 4hrs
cog.table.category.diffabund.4h <- cog.table.category.diffabund %>% left_join(metadata, by = 'Sample') %>% filter(Timepoint == "4h") %>% filter(Fraction == "BC_plus" | Fraction == "BC_minus") %>% select(-c("Fraction", "Desc_SampleID", "Bio_Rep",	"Tech_Rep", "Treatment", "Timepoint")) %>% relocate(Sample) %>% as.tibble() %>% t()
colnames(cog.table.category.diffabund.4h) <- cog.table.category.diffabund.4h[1,]
cog.table.category.diffabund.4h <- cog.table.category.diffabund.4h[2:nrow(cog.table.category.diffabund.4h),]
cog.table.category.diffabund.4h <- cog.table.category.diffabund.4h %>% as.data.frame() %>% rownames_to_column() %>% mutate(COG_ID = rowname) %>% select(-rowname) %>% relocate(COG_ID)
cog.table.category.diffabund.4h <- cog.table.category.diffabund.4h %>% mutate_at(vars(-("COG_ID")), as.integer) %>% mutate_at("COG_ID", as.factor)
cog.table.category.diffabund.4h <- cog.table.category.diffabund.4h[1:nrow(cog.table.category.diffabund.4h)-1,]


# Get diffabund table for 21hrs
cog.table.category.diffabund.21h <- cog.table.category.diffabund %>% left_join(metadata, by = 'Sample') %>% filter(Timepoint == "21h") %>% filter(Fraction == "BC_plus" | Fraction == "BC_minus") %>% select(-c("Fraction", "Desc_SampleID", "Bio_Rep",	"Tech_Rep", "Treatment", "Timepoint")) %>% relocate(Sample) %>% as.tibble() %>% t()
colnames(cog.table.category.diffabund.21h) <- cog.table.category.diffabund.21h[1,]
cog.table.category.diffabund.21h <- cog.table.category.diffabund.21h[2:nrow(cog.table.category.diffabund.21h),]
cog.table.category.diffabund.21h <- cog.table.category.diffabund.21h %>% as.data.frame() %>% rownames_to_column() %>% mutate(COG_ID = rowname) %>% select(-rowname) %>% relocate(COG_ID)
cog.table.category.diffabund.21h <- cog.table.category.diffabund.21h %>% mutate_at(vars(-("COG_ID")), as.integer) %>% mutate_at("COG_ID", as.factor)


# Get metadata table for deseq2
deseq.md.4h <- metadata %>% filter(Timepoint=="4h") %>% filter(Fraction != "SYTO_plus")
deseq.md.21h <- metadata %>% filter(Timepoint=="21h") %>% filter(Fraction != "SYTO_plus")


#library(DESeq2)
# Note: loading DESeq2 will mess up some tidyverse functions!
# Create DESeq2 dataset
dds.category.4h <- DESeqDataSetFromMatrix(countData=cog.table.category.diffabund.4h, colData=deseq.md.4h, design = ~ Fraction, tidy=T)
dds.category.21h <- DESeqDataSetFromMatrix(countData=cog.table.category.diffabund.21h, colData=deseq.md.21h, design = ~ Fraction, tidy=T)


# Run DESeq function
dds.category.4h <- DESeq(dds.category.4h)
dds.category.21h <- DESeq(dds.category.21h)

# Contrast Fractions
res.category.4h.fraction <- results(dds.category.4h, contrast=c("Fraction", "BC_plus", "BC_minus"), tidy = TRUE)
res.category.21h.fraction <- results(dds.category.21h, contrast=c("Fraction", "BC_plus", "BC_minus"), tidy = TRUE)

# Look at results
#res.SYTO4h <- results(dds)
summary(res.category.4h.fraction)
res.category.4h.fraction <- res.category.4h.fraction[order(res.category.4h.fraction$padj),]
summary(res.category.21h.fraction)
res.category.21h.fraction <- res.category.21h.fraction[order(res.category.21h.fraction$padj),]


# Fold change bar plot
res.category.4h.boncat.plot <- res.category.4h.fraction %>% filter(padj <= 0.05) %>% mutate(COG_Category_1_Annotation = row) %>% select(-row) #%>% left_join(., cog.categories, by = 'COG_Category_1_Annotation') #%>% left_join(., cog.annotations, by = 'COG_Category_ID')
res.category.21h.boncat.plot <- res.category.21h.fraction %>% filter(padj <= 0.05) %>% mutate(COG_Category_1_Annotation = row) %>% select(-row) #%>% left_join(., cog.annotations, by = 'COG_ID') %>% left_join(., cog.categories, by = 'COG_Category_ID')

#dev.new()
p1 <- ggplot(res.category.4h.boncat.plot, aes(x = reorder(COG_Category_1_Annotation, log2FoldChange), y = log2FoldChange)) + geom_bar(stat='identity') + theme_bw() + ylim(c(-1,1))  + theme(axis.title.y = element_blank()) + coord_flip() + theme(legend.position = "none")
p2 <- ggplot(res.category.4h.boncat.plot, aes(x = reorder(COG_Category_1_Annotation, log2FoldChange), y = padj)) + geom_point() + theme_bw() + ylim(c(0,0.05)) + theme(axis.text.y=element_blank(), axis.title.y = element_blank()) + coord_flip() + theme(legend.position = "none")
ggarrange(p1, p2, ncol=2, nrow=1, legend="none", widths = c(1, 0.25))

p1 <- ggplot(res.category.21h.boncat.plot, aes(x = reorder(COG_Category_1_Annotation, log2FoldChange), y = log2FoldChange)) + geom_bar(stat='identity') + theme_bw() + ylim(c(-1,1)) + theme(axis.title.y = element_blank()) + coord_flip() + theme(legend.position = "none")
p2 <- ggplot(res.category.21h.boncat.plot, aes(x = reorder(COG_Category_1_Annotation, log2FoldChange), y = padj)) + geom_point() + theme_bw() + ylim(c(0,0.05)) + theme(axis.text.y=element_blank(), axis.title.y = element_blank()) + coord_flip() + theme(legend.position = "none")
ggarrange(p1, p2, ncol=2, nrow=1, legend="none", widths = c(1, 0.25))


#write_excel_csv(res.category.4h.boncat.plot,"manuscript_data_output/DESeq2_COG_Category_4h.csv")
#write_excel_csv(res.category.21h.boncat.plot,"manuscript_data_output/DESeq2_COG_Category_21h.csv")
