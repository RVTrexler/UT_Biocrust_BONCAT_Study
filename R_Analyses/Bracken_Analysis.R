# BONCAT Biocrust Project - Analysis of bracken output data of metagenomes
# Version 3.0
# Ryan V. Trexler
# 12/06/2022

# Clear environment
remove(list =ls())

# Load libraries
library(tidyverse)
library(vegan)
library(microshades)
library(phyloseq)
library(speedyseq)
library(cowplot)
library(patchwork)
library(ggpubr)


setwd("~/Documents/PSU_Documents/Couradeau_Lab/Biocrust_Boncat_MG/Metagenome_Analysis/kraken_bracken")


#### 1.1. Import metadata ####

# Import metadata
metadata <- read.table(file = "Metadata_key.txt", sep = '\t', header = T)
# Note: metadata_key.txt was modified from previous version to make an additional column for a Descriptive_Sample_ID column

#### 1.2. Import and format bracken data ####

bracken.table <- read_tsv(file = "bracken_combined_v2.txt") %>% filter(Domain %in% c("Bacteria", "Archaea"))

# Select only the 'num' data (number of reads per taxa), and delete "_num" from the column names.
bracken.table <- bracken.table %>% select(-ends_with("_frac")) %>% rename_at(.vars = vars(ends_with("_num")),.funs = funs(sub("_num$", "", .)))

# Concatenate species and taxid to make unique
bracken.table <- bracken.table %>% unite(species_id, Species:tax_id, sep = "_", remove = F) #unite("z", x:y, na.rm = TRUE, remove = FALSE)


#### 2.0 NMDS ####

# Format the data
bracken.table.frac <- bracken.table %>% select(-c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "tax_level", "tax_id"))
bracken.table.frac <- t(bracken.table.frac)
colnames(bracken.table.frac) <- bracken.table.frac[1,1:ncol(bracken.table.frac)]
bracken.table.frac <- bracken.table.frac[-1,]
class(bracken.table.frac) <- "numeric"

# Convert to fraction
bracken.table.frac <- bracken.table.frac/rowSums(bracken.table.frac)*100

# Calculate Bray-Curtis distance between samples
bracken.table.bray <- vegdist(bracken.table.frac)

bracken.table.nmds <- metaMDS(bracken.table.frac, k=2) # Perform NMDS on Bray Curtis distances
# Get NMDS scores
bracken.nmds.scores <- as.data.frame(scores(bracken.table.nmds, "sites")) #Using the scores function from vegan to extract the site scores and convert to a data.frame
bracken.nmds.scores <- bracken.nmds.scores %>% rownames_to_column(., var = "Sample") %>% left_join(metadata, by = "Sample")

# Get NMDS Group (Treatment) centroids
centroids <- aggregate(cbind(NMDS1, NMDS2) ~ Treatment, data = bracken.nmds.scores, FUN = mean) %>% rename(Cent_NMDS1 = NMDS1) %>% rename(Cent_NMDS2 = NMDS2)
# Get NMDS dispersions (distances from group centroid)
distances <- bracken.nmds.scores %>% left_join(centroids, by = "Treatment") %>% mutate(Distance = sqrt(((NMDS1-Cent_NMDS1)^2 + (NMDS2-Cent_NMDS2)^2)))

# Visualize by Fraction
#dev.new()
ggplot(bracken.nmds.scores, aes(x=NMDS1,y=NMDS2)) + geom_point(aes(shape=Timepoint, color=Fraction,fill=Fraction),size=4) + labs(title="NMDS of Bracken taxonomy classifications") + theme(plot.title = element_text(hjust = 0.5)) + theme(panel.background = element_rect(fill = "white", colour = "grey50"), legend.position = c(0.75,0.7), legend.background = element_rect(color = "black", fill = "grey90", linewidth = 0.5, linetype = "solid")) + geom_text(x=-0.3, y=0.0, label= paste("Stress: ", as.character(round(bracken.table.nmds$stress,4))))
# To add group centroids: # + geom_point(data = centroids %>% rename(NMDS1 = Cent_NMDS1) %>% rename(NMDS2 = Cent_NMDS2) , color="black")

# Permanova of bray-curtis distances
bracken.table.bray.permanova <- adonis(bracken.table.bray ~ Fraction*Timepoint, dat=metadata, perm=1e3)
bracken.table.bray.permanova$aov.tab


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

#write_excel_csv(distances.sum,"manuscript_data_output/Taxa_BetaDisp_summary.csv")


#### 3.0 Phyloseq Analysis - 100% Stacked Bar Plots at Phylum level ####

# Make Phyloseq 'esv' table object
PhySq_ESV <- bracken.table %>% select(-Domain, -Phylum, -Class, -Order, -Family, -Genus, -Species, -tax_id, -tax_level) %>% as.data.frame()
rownames(PhySq_ESV) <- PhySq_ESV[,1]
PhySq_ESV[,1] <- NULL
PhySq_ESV <- PhySq_ESV/rowSums(PhySq_ESV)*100
PhySq_ESV <- PhySq_ESV %>% otu_table(taxa_are_rows=T)

# Make Phyloseq metadata object

# Make factor levels correct order
metadata$Treatment <- factor(metadata$Treatment, levels=c("Dry", "4h_SYTO_plus", "4h_BC_minus", "4h_BC_plus", "21h_SYTO_plus", "21h_BC_minus", "21h_BC_plus"))
#  Make phyloseq object
PhySq_metadat <- metadata
rownames(PhySq_metadat) <- PhySq_metadat[,1]
PhySq_metadat[,1] <- NULL
PhySq_metadat <- sample_data(PhySq_metadat,errorIfNULL=TRUE)

# Make Phyloseq taxonomy object
PhySq_taxo <- bracken.table %>% select(species_id, Domain, Phylum, Class, Family, Genus, Species) %>% as.data.frame()
rownames(PhySq_taxo) <- PhySq_taxo[,1]
PhySq_taxo[,1] <- NULL
PhySq_taxo <- PhySq_taxo %>% as.matrix() %>% tax_table(errorIfNULL=TRUE)

# Create PhyloSeq object
PhySq.16S <- phyloseq(PhySq_taxo, PhySq_ESV, PhySq_metadat)

# ggplot Rank Abundance Curve - Select entire dataset, top x. Note currently, the RA is not divided correctly
#RACurve.all <- sort(taxa_sums(PhySq.16S), TRUE) %>% '/'(nsamples(PhySq.16S)) %>% as.data.frame() %>% rownames_to_column(var="OTUID") %>% 'colnames<-'(c("OTUID","RA")) %>% top_n(40,wt=RA)



# Agglomerate and normalize the phyloseq object, and melt to a data frame
mdf_prep <- prep_mdf(PhySq.16S)

# Generate a color object for the specified data
color_objs_GP <- create_color_dfs(mdf_prep,selected_groups = c("Cyanobacteria", "Proteobacteria", "Firmicutes", "Actinobacteria", "Bacteroidota") , cvd = TRUE)

# Extract
mdf_GP <- color_objs_GP$mdf
cdf_GP <- color_objs_GP$cdf

#write_excel_csv(mdf_GP,"manuscript_data_output/Figure_1_data.csv")

# Plot
plot <- plot_microshades(mdf_GP, cdf_GP)

plot_1 <- plot + scale_y_continuous(labels = scales::percent, expand = expansion(0)) +
  theme(legend.key.size = unit(0.2, "cm"), text=element_text(size=10)) +
  theme(axis.text.x = element_text(size= 6)) 

#plot_1 



plot_2 <- plot + scale_y_continuous(labels = scales::percent, expand = expansion(0)) +
  theme(legend.key.size = unit(0.2, "cm"), text=element_text(size=10)) +
  theme(axis.text.x = element_text(size= 6)) +
  facet_wrap(~Treatment, scales = "free_x", nrow = 1) +
  theme (strip.text.x = element_text(size = 6))

#plot_2


GP_legend <-custom_legend(mdf_GP, cdf_GP)

plot_3 <- plot + scale_y_continuous(labels = scales::percent, expand = expansion(0)) +
  theme(legend.position = "none")  +
  theme(axis.text.x = element_text(size= 6)) +
  facet_wrap(~Treatment, scales = "free_x", nrow = 1) +
  theme(axis.text.x = element_text(size= 6)) + 
  theme(plot.margin = margin(6,20,6,6))

# Final microshades plot
plot_grid(plot_3, GP_legend,  rel_widths = c(1, .25))


#### 4.0 Aplha diversity plots ####

# Make sure DESeq2 is not loaded....
bracken.table.div <- bracken.table.frac
bracken.table.div <- merge(as.data.frame(diversity(bracken.table.div,"shannon")), as.data.frame(rowSums(bracken.table.div > 0)), by="row.names") %>% rename(Sample = Row.names, Shannon_Div = 'diversity(bracken.table.div, "shannon")', ESV_Richness = 'rowSums(bracken.table.div > 0)') %>% mutate(Pielous_evenness = Shannon_Div/log(ESV_Richness))


# Merge diversity with metadata
metadata.div <- left_join(metadata, bracken.table.div, by ="Sample")

# Collapse technical replicates to means
#metadata.div <- metadata.div %>% group_by(Timepoint, Fraction, Bio_Rep) %>% summarize(Shannon_Div = mean(Shannon_Div), ESV_Richness = mean(ESV_Richness), Pielous_evenness = mean(Pielous_evenness))

### ANOVA Statistics ###
# ESV Richness
div.anova.er <- aov(ESV_Richness ~ Fraction * Timepoint, data = metadata.div)
summary(div.anova.er)
TukeyHSD(div.anova.er)

qqnorm(metadata.div$ESV_Richness, pch = 1, frame = FALSE)
qqline(metadata.div$ESV_Richness, col = "steelblue", lwd = 2)
shapiro.test(div.anova.er$residuals)

# Shannon Diversity
div.anova.sd <- aov(Shannon_Div ~ Fraction * Timepoint, data = metadata.div)
summary(div.anova.sd)
TukeyHSD(div.anova.sd)

qqnorm(metadata.div$Shannon_Div, pch = 1, frame = FALSE)
qqline(metadata.div$Shannon_Div, col = "steelblue", lwd = 2)
shapiro.test(div.anova.sd$residuals)

# Pielou's evenness
div.anova.pe <- aov(Pielous_evenness ~ Fraction * Timepoint, data = metadata.div)
summary(div.anova.pe)
TukeyHSD(div.anova.pe)

qqnorm(metadata.div$Pielous_evenness, pch = 1, frame = FALSE)
qqline(metadata.div$Pielous_evenness, col = "steelblue", lwd = 2)
hist(div.anova.pe$residuals)
shapiro.test(div.anova.pe$residuals)


# Sort it so that timepoint is ordered correctly
metadata.div$Timepoint <- factor(metadata.div$Timepoint,levels = c("Dry", "4h", "21h"))
metadata.div$Fraction <- factor(metadata.div$Fraction,levels = c("Dry", "SYTO_plus", "BC_minus", "BC_plus"))

# Plot calculated diversity via boxplots
p1 <- ggplot(metadata.div, aes(x=Timepoint, y=ESV_Richness, fill=Fraction)) + geom_boxplot() + geom_point(position=position_dodge(width=0.75), aes(group=Fraction)) #+ ggtitle("Total MG Reads") + theme(plot.title = element_text(size=10))
p2 <- ggplot(metadata.div, aes(x=Timepoint, y=Shannon_Div, fill=Fraction)) + geom_boxplot() + geom_point(position=position_dodge(width=0.75), aes(group=Fraction)) #+ ggtitle("Total MG Reads") + theme(plot.title = element_text(size=10))
p3 <- ggplot(metadata.div, aes(x=Timepoint, y=Pielous_evenness, fill=Fraction)) + geom_boxplot() + geom_point(position=position_dodge(width=0.75), aes(group=Fraction)) #+ ggtitle("Total MG Reads") + theme(plot.title = element_text(size=10))
#dev.new()
#grid.arrange(p1,p2,p3,ncol=3,nrow=1,top=textGrob("Alpha Diversity: Timepoint x Fraction ",gp=gpar(fontsize=20)))
ggarrange(p1, p2, p3, ncol=3, nrow=1, legend="none", widths = c(1,1,1))

#write_excel_csv(metadata.div,"manuscript_data_output/Figure_2_data.csv")
metadata.div.sum <- metadata.div %>% group_by(Timepoint,Fraction) %>% summarize(Shannon_Div_mean = mean(Shannon_Div), Shannon_Div_se = sd(Shannon_Div)/(sqrt(n())), ESV_Richness_mean = mean(ESV_Richness), ESV_Richness_se = sd(ESV_Richness)/(sqrt(n())), Pielous_evenness_mean = mean(Pielous_evenness), Pielous_evenness_se = sd(Pielous_evenness)/(sqrt(n())))
#write_excel_csv(metadata.div.sum,"manuscript_data_output/Figure_2_summarized_data.csv")


#### 5.0 DESeq2 on taxa at species level ####

braken.table.sp <- bracken.table %>% select(-c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "tax_id", "tax_level")) %>% t()
colnames(braken.table.sp) <- braken.table.sp[1,1:ncol(braken.table.sp)]
braken.table.sp <- braken.table.sp[-1,]
braken.table.sp <- rownames_to_column(as.data.frame(braken.table.sp)) %>% rename(Sample = "rowname")

# Get diffabund table for 4hrs
bracken.table.sp.diffabund.4h <- braken.table.sp %>% left_join(metadata, by = 'Sample') %>% filter(Timepoint == "4h") %>% filter(Fraction == "BC_plus" | Fraction == "BC_minus") %>% select(-c("Fraction", "Desc_SampleID", "Bio_Rep",	"Tech_Rep", "Treatment", "Timepoint")) %>% relocate(Sample) %>% as.tibble() %>% t()
colnames(bracken.table.sp.diffabund.4h) <- bracken.table.sp.diffabund.4h[1,]
bracken.table.sp.diffabund.4h <- bracken.table.sp.diffabund.4h[2:nrow(bracken.table.sp.diffabund.4h),]
bracken.table.sp.diffabund.4h <- bracken.table.sp.diffabund.4h %>% as.data.frame() %>% rownames_to_column() %>% mutate(species_id = rowname) %>% select(-rowname) %>% relocate(species_id)
bracken.table.sp.diffabund.4h <- bracken.table.sp.diffabund.4h %>% mutate_at(vars(-("species_id")), as.integer) %>% mutate_at("species_id", as.factor)

# Get diffabund table for 21hrs
bracken.table.sp.diffabund.21h <- braken.table.sp %>% left_join(metadata, by = 'Sample') %>% filter(Timepoint == "21h") %>% filter(Fraction == "BC_plus" | Fraction == "BC_minus") %>% select(-c("Fraction", "Desc_SampleID", "Bio_Rep",	"Tech_Rep", "Treatment", "Timepoint")) %>% relocate(Sample) %>% as.tibble() %>% t()
colnames(bracken.table.sp.diffabund.21h) <- bracken.table.sp.diffabund.21h[1,]
bracken.table.sp.diffabund.21h <- bracken.table.sp.diffabund.21h[2:nrow(bracken.table.sp.diffabund.21h),]
bracken.table.sp.diffabund.21h <- bracken.table.sp.diffabund.21h %>% as.data.frame() %>% rownames_to_column() %>% mutate(species_id = rowname) %>% select(-rowname) %>% relocate(species_id)
bracken.table.sp.diffabund.21h <- bracken.table.sp.diffabund.21h %>% mutate_at(vars(-("species_id")), as.integer) %>% mutate_at("species_id", as.factor)


# Get metadata table for deseq2
deseq.md.4h <- metadata %>% filter(Timepoint=="4h") %>% filter(Fraction != "SYTO_plus")
deseq.md.21h <- metadata %>% filter(Timepoint=="21h") %>% filter(Fraction != "SYTO_plus")


#library(DESeq2)
# Note: loading DESeq2 will mess up some tidyverse functions!
# Create DESeq2 dataset
dds.4h <- DESeqDataSetFromMatrix(countData=bracken.table.sp.diffabund.4h, colData=deseq.md.4h, design = ~ Fraction, tidy=T)
dds.21h <- DESeqDataSetFromMatrix(countData=bracken.table.sp.diffabund.21h, colData=deseq.md.21h, design = ~ Fraction, tidy=T)


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
res.4h.boncat.plot <- res.4h.fraction %>% filter(padj <= 0.05) %>% mutate(species_id = row) %>% select(-row) %>% left_join(bracken.table %>% select(c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "species_id", "tax_id", "tax_level")), by = 'species_id')
res.21h.boncat.plot <- res.21h.fraction %>% filter(padj <= 0.05) %>% mutate(species_id = row) %>% select(-row) %>% left_join(bracken.table %>% select(c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "species_id", "tax_id", "tax_level")), by = 'species_id')


#dev.new()
p1 <- ggplot(res.4h.boncat.plot, aes(x = reorder(species_id, log2FoldChange), y = log2FoldChange, fill = Phylum)) + geom_bar(stat='identity') + theme_bw() + theme(axis.title.y = element_blank()) + coord_flip() #+ theme(legend.position = "none")
p2 <- ggplot(res.4h.boncat.plot, aes(x = reorder(species_id, log2FoldChange), y = padj)) + geom_point() + theme_bw() + theme(axis.text.y=element_blank(), axis.title.y = element_blank()) + coord_flip() + theme(legend.position = "none")
leg <- get_legend(p1) %>% as_ggplot()
plot <- ggarrange(p1, p2, leg, ncol=3, nrow=1, legend="none", widths = c(2, 0.25, 1))
plot
annotate_figure(plot, top = text_grob("Species-level Differential Abundance - 4hrs", face = "bold", size = 16))

p1 <- ggplot(res.21h.boncat.plot, aes(x = reorder(species_id, log2FoldChange), y = log2FoldChange, fill = Phylum)) + geom_bar(stat='identity') + theme_bw() + theme(axis.title.y = element_blank()) + coord_flip() #+ theme(legend.position = "none")
p2 <- ggplot(res.21h.boncat.plot, aes(x = reorder(species_id, log2FoldChange), y = padj)) + geom_point() + theme_bw() + theme(axis.text.y=element_blank(), axis.title.y = element_blank()) + coord_flip() + theme(legend.position = "none")
leg <- get_legend(p1) %>% as_ggplot()
plot <- ggarrange(p1, p2, leg, ncol=3, nrow=1, legend="none", widths = c(2, 0.25, 1))
plot
annotate_figure(plot, top = text_grob("Species-level Differential Abundance - 21hrs", face = "bold", size = 16))

#write_excel_csv(res.4h.boncat.plot,"manuscript_data_output/DESeq2_Species_4hr.csv")
#write_excel_csv(res.21h.boncat.plot,"manuscript_data_output/DESeq2_Species_21hr.csv")

#### 6.0 DESeq2 on taxa at genus level ####

# Summarize to genus level
braken.table.gen <- bracken.table %>% pivot_longer(-c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "species_id", "tax_id", "tax_level"), names_to = "Sample", values_to = "Species_Count") %>% group_by(Sample, Genus) %>% summarize(Genus_Count = sum(Species_Count))
# Pivot it back to wide
braken.table.gen <- braken.table.gen %>% pivot_wider(names_from = Sample, values_from = Genus_Count) %>% t()
# Remove last column because it is NA (Genera that were not assigned)
braken.table.gen <- braken.table.gen[,1:ncol(braken.table.gen)-1]
colnames(braken.table.gen) <- braken.table.gen[1,1:ncol(braken.table.gen)] 
braken.table.gen <- braken.table.gen[-1,]
braken.table.gen <- rownames_to_column(as.data.frame(braken.table.gen)) %>% rename(Sample = "rowname")

# Get diffabund table for 4hrs
bracken.table.gen.diffabund.4h <- braken.table.gen %>% left_join(metadata, by = 'Sample') %>% filter(Timepoint == "4h") %>% filter(Fraction == "BC_plus" | Fraction == "BC_minus") %>% select(-c("Fraction", "Desc_SampleID", "Bio_Rep",	"Tech_Rep", "Treatment", "Timepoint")) %>% relocate(Sample) %>% as.tibble() %>% t()
colnames(bracken.table.gen.diffabund.4h) <- bracken.table.gen.diffabund.4h[1,]
bracken.table.gen.diffabund.4h <- bracken.table.gen.diffabund.4h[2:nrow(bracken.table.gen.diffabund.4h),]
bracken.table.gen.diffabund.4h <- bracken.table.gen.diffabund.4h %>% as.data.frame() %>% rownames_to_column() %>% mutate(species_id = rowname) %>% select(-rowname) %>% relocate(species_id)
bracken.table.gen.diffabund.4h <- bracken.table.gen.diffabund.4h %>% mutate_at(vars(-("species_id")), as.integer) %>% mutate_at("species_id", as.factor)

# Get diffabund table for 21hrs
bracken.table.gen.diffabund.21h <- braken.table.gen %>% left_join(metadata, by = 'Sample') %>% filter(Timepoint == "21h") %>% filter(Fraction == "BC_plus" | Fraction == "BC_minus") %>% select(-c("Fraction", "Desc_SampleID", "Bio_Rep",	"Tech_Rep", "Treatment", "Timepoint")) %>% relocate(Sample) %>% as.tibble() %>% t()
colnames(bracken.table.gen.diffabund.21h) <- bracken.table.gen.diffabund.21h[1,]
bracken.table.gen.diffabund.21h <- bracken.table.gen.diffabund.21h[2:nrow(bracken.table.gen.diffabund.21h),]
bracken.table.gen.diffabund.21h <- bracken.table.gen.diffabund.21h %>% as.data.frame() %>% rownames_to_column() %>% mutate(species_id = rowname) %>% select(-rowname) %>% relocate(species_id)
bracken.table.gen.diffabund.21h <- bracken.table.gen.diffabund.21h %>% mutate_at(vars(-("species_id")), as.integer) %>% mutate_at("species_id", as.factor)


# Get metadata table for deseq2
deseq.md.4h <- metadata %>% filter(Timepoint=="4h") %>% filter(Fraction != "SYTO_plus")
deseq.md.21h <- metadata %>% filter(Timepoint=="21h") %>% filter(Fraction != "SYTO_plus")


library(DESeq2)
# Note: loading DESeq2 will mess up some tidyverse functions!
# Create DESeq2 dataset
dds.4h <- DESeqDataSetFromMatrix(countData=bracken.table.gen.diffabund.4h, colData=deseq.md.4h, design = ~ Fraction, tidy=T)
dds.21h <- DESeqDataSetFromMatrix(countData=bracken.table.gen.diffabund.21h, colData=deseq.md.21h, design = ~ Fraction, tidy=T)


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
tax <- bracken.table %>% select(c("Domain", "Phylum", "Class", "Order", "Family", "Genus")) %>% distinct()
res.4h.boncat.plot <- res.4h.fraction %>% filter(padj <= 0.05) %>% mutate(Genus = row) %>% select(-row) %>% left_join(tax, by = 'Genus') %>% arrange(desc(log2FoldChange)) #%>% mutate(Genus=factor(Genus, levels=Genus))
res.21h.boncat.plot <- res.21h.fraction %>% filter(padj <= 0.05) %>% mutate(Genus = row) %>% select(-row) %>% left_join(tax, by = 'Genus') %>% arrange(desc(log2FoldChange)) #%>% mutate(Genus=factor(Genus, levels=Genus))


#dev.new()
p1 <- ggplot(res.4h.boncat.plot, aes(x = reorder(Genus, log2FoldChange), y = log2FoldChange, fill = Phylum)) + geom_bar(stat='identity') + theme_bw() + theme(axis.title.y = element_blank()) + coord_flip() #+ theme(legend.position = "none")
p2 <- ggplot(res.4h.boncat.plot, aes(x = reorder(Genus, log2FoldChange), y = padj)) + geom_point() + theme_bw() + theme(axis.text.y=element_blank(), axis.title.y = element_blank()) + coord_flip() + theme(legend.position = "none")
leg <- get_legend(p1) %>% as_ggplot()
plot <- ggarrange(p1, p2, leg, ncol=3, nrow=1, legend="none", widths = c(2, 0.25, 1))
plot
annotate_figure(plot, top = text_grob("Genus-level Differential Abundance - 4hrs", face = "bold", size = 16))

p1 <- ggplot(res.21h.boncat.plot, aes(x = reorder(Genus, log2FoldChange), y = log2FoldChange, fill = Phylum)) + geom_bar(stat='identity') + theme_bw() + theme(axis.title.y = element_blank()) + coord_flip() #+ theme(legend.position = "none")
p2 <- ggplot(res.21h.boncat.plot, aes(x = reorder(Genus, log2FoldChange), y = padj)) + geom_point() + theme_bw() + theme(axis.text.y=element_blank(), axis.title.y = element_blank()) + coord_flip() + theme(legend.position = "none")
leg <- get_legend(p1) %>% as_ggplot()
plot <- ggarrange(p1, p2, leg, ncol=3, nrow=1, legend="none", widths = c(2, 0.25, 1))
plot
annotate_figure(plot, top = text_grob("Genus-level Differential Abundance - 21hrs", face = "bold", size = 16))

#write_excel_csv(res.4h.boncat.plot,"manuscript_data_output/DESeq2_Genus_4hr.csv")
#write_excel_csv(res.21h.boncat.plot,"manuscript_data_output/DESeq2_Genus_21hr.csv")


#### 7.0 DESeq2 on taxa at family level ####

# Summarize to family level
braken.table.fam <- bracken.table %>% pivot_longer(-c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "species_id", "tax_id", "tax_level"), names_to = "Sample", values_to = "Species_Count") %>% group_by(Sample, Family) %>% summarize(Family_Count = sum(Species_Count))
# Pivot it back to wide
braken.table.fam <- braken.table.fam %>% pivot_wider(names_from = Sample, values_from = Family_Count) %>% t()
# Remove last column because it is NA (Genera that were not assigned)
braken.table.fam <- braken.table.fam[,1:ncol(braken.table.fam)-1]
colnames(braken.table.fam) <- braken.table.fam[1,1:ncol(braken.table.fam)] 
braken.table.fam <- braken.table.fam[-1,]
braken.table.fam <- rownames_to_column(as.data.frame(braken.table.fam)) %>% rename(Sample = "rowname")

# Get diffabund table for 4hrs
bracken.table.fam.diffabund.4h <- braken.table.fam %>% left_join(metadata, by = 'Sample') %>% filter(Timepoint == "4h") %>% filter(Fraction == "BC_plus" | Fraction == "BC_minus") %>% select(-c("Fraction", "Desc_SampleID", "Bio_Rep",	"Tech_Rep", "Treatment", "Timepoint")) %>% relocate(Sample) %>% as.tibble() %>% t()
colnames(bracken.table.fam.diffabund.4h) <- bracken.table.fam.diffabund.4h[1,]
bracken.table.fam.diffabund.4h <- bracken.table.fam.diffabund.4h[2:nrow(bracken.table.fam.diffabund.4h),]
bracken.table.fam.diffabund.4h <- bracken.table.fam.diffabund.4h %>% as.data.frame() %>% rownames_to_column() %>% mutate(species_id = rowname) %>% select(-rowname) %>% relocate(species_id)
bracken.table.fam.diffabund.4h <- bracken.table.fam.diffabund.4h %>% mutate_at(vars(-("species_id")), as.integer) %>% mutate_at("species_id", as.factor)

# Get diffabund table for 21hrs
bracken.table.fam.diffabund.21h <- braken.table.fam %>% left_join(metadata, by = 'Sample') %>% filter(Timepoint == "21h") %>% filter(Fraction == "BC_plus" | Fraction == "BC_minus") %>% select(-c("Fraction", "Desc_SampleID", "Bio_Rep",	"Tech_Rep", "Treatment", "Timepoint")) %>% relocate(Sample) %>% as.tibble() %>% t()
colnames(bracken.table.fam.diffabund.21h) <- bracken.table.fam.diffabund.21h[1,]
bracken.table.fam.diffabund.21h <- bracken.table.fam.diffabund.21h[2:nrow(bracken.table.fam.diffabund.21h),]
bracken.table.fam.diffabund.21h <- bracken.table.fam.diffabund.21h %>% as.data.frame() %>% rownames_to_column() %>% mutate(species_id = rowname) %>% select(-rowname) %>% relocate(species_id)
bracken.table.fam.diffabund.21h <- bracken.table.fam.diffabund.21h %>% mutate_at(vars(-("species_id")), as.integer) %>% mutate_at("species_id", as.factor)


# Get metadata table for deseq2
deseq.md.4h <- metadata %>% filter(Timepoint=="4h") %>% filter(Fraction != "SYTO_plus")
deseq.md.21h <- metadata %>% filter(Timepoint=="21h") %>% filter(Fraction != "SYTO_plus")


library(DESeq2)
# Note: loading DESeq2 will mess up some tidyverse functions!
# Create DESeq2 dataset
dds.4h <- DESeqDataSetFromMatrix(countData=bracken.table.fam.diffabund.4h, colData=deseq.md.4h, design = ~ Fraction, tidy=T)
dds.21h <- DESeqDataSetFromMatrix(countData=bracken.table.fam.diffabund.21h, colData=deseq.md.21h, design = ~ Fraction, tidy=T)


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
tax <- bracken.table %>% select(c("Domain", "Phylum", "Class", "Order", "Family")) %>% distinct()
res.4h.boncat.plot <- res.4h.fraction %>% filter(padj <= 0.05) %>% mutate(Family = row) %>% select(-row) %>% left_join(tax, by = 'Family') %>% arrange(desc(log2FoldChange)) #%>% mutate(Genus=factor(Genus, levels=Genus))
res.21h.boncat.plot <- res.21h.fraction %>% filter(padj <= 0.05) %>% mutate(Family = row) %>% select(-row) %>% left_join(tax, by = 'Family') %>% arrange(desc(log2FoldChange)) #%>% mutate(Genus=factor(Genus, levels=Genus))

# Plot only ones with L2FC >= |1|
res.4h.boncat.plot <- res.4h.boncat.plot %>% filter(log2FoldChange >= 1 | log2FoldChange <= -1)

#dev.new()
p1 <- ggplot(res.4h.boncat.plot, aes(x = reorder(Family, log2FoldChange), y = log2FoldChange, fill = Phylum)) + geom_bar(stat='identity') + theme_bw() + theme(axis.title.y = element_blank()) + coord_flip() #+ theme(legend.position = "none")
p2 <- ggplot(res.4h.boncat.plot, aes(x = reorder(Family, log2FoldChange), y = padj)) + geom_point() + theme_bw() + theme(axis.text.y=element_blank(), axis.title.y = element_blank()) + coord_flip() + theme(legend.position = "none")
leg <- get_legend(p1) %>% as_ggplot()
plot <- ggarrange(p1, p2, leg, ncol=3, nrow=1, legend="none", widths = c(2, 0.25, 1))
plot
annotate_figure(plot, top = text_grob("Family-level Differential Abundance - 4hrs", face = "bold", size = 16))


p1 <- ggplot(res.21h.boncat.plot, aes(x = reorder(Family, log2FoldChange), y = log2FoldChange, fill = Phylum)) + geom_bar(stat='identity') + theme_bw() + theme(axis.title.y = element_blank()) + coord_flip() #+ theme(legend.position = "none")
p2 <- ggplot(res.21h.boncat.plot, aes(x = reorder(Family, log2FoldChange), y = padj)) + geom_point() + theme_bw() + theme(axis.text.y=element_blank(), axis.title.y = element_blank()) + coord_flip() + theme(legend.position = "none")
leg <- get_legend(p1) %>% as_ggplot()
plot <- ggarrange(p1, p2, leg, ncol=3, nrow=1, legend="none", widths = c(2, 0.25, 1))
plot
annotate_figure(plot, top = text_grob("Family-level Differential Abundance - 21hrs", face = "bold", size = 16))

#write_excel_csv(res.4h.boncat.plot,"manuscript_data_output/DESeq2_Family_4hr.csv")
#write_excel_csv(res.21h.boncat.plot,"manuscript_data_output/DESeq2_Family_21hr.csv")


#### 8.0 DESeq2 on taxa at phylum level ####

# Summarize to phylum level
braken.table.phyl <- bracken.table %>% pivot_longer(-c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "species_id", "tax_id", "tax_level"), names_to = "Sample", values_to = "Species_Count") %>% group_by(Sample, Phylum) %>% summarize(Phylum_Count = sum(Species_Count))
# Pivot it back to wide
braken.table.phyl <- braken.table.phyl %>% pivot_wider(names_from = Sample, values_from = Phylum_Count) %>% t()
# Remove last column because it is NA (Genera that were not assigned)
braken.table.phyl <- braken.table.phyl[,1:ncol(braken.table.phyl)-1]
colnames(braken.table.phyl) <- braken.table.phyl[1,1:ncol(braken.table.phyl)] 
braken.table.phyl <- braken.table.phyl[-1,]
braken.table.phyl <- rownames_to_column(as.data.frame(braken.table.phyl)) %>% rename(Sample = "rowname")

# Get diffabund table for 4hrs
bracken.table.phyl.diffabund.4h <- braken.table.phyl %>% left_join(metadata, by = 'Sample') %>% filter(Timepoint == "4h") %>% filter(Fraction == "BC_plus" | Fraction == "BC_minus") %>% select(-c("Fraction", "Desc_SampleID", "Bio_Rep",	"Tech_Rep", "Treatment", "Timepoint")) %>% relocate(Sample) %>% as.tibble() %>% t()
colnames(bracken.table.phyl.diffabund.4h) <- bracken.table.phyl.diffabund.4h[1,]
bracken.table.phyl.diffabund.4h <- bracken.table.phyl.diffabund.4h[2:nrow(bracken.table.phyl.diffabund.4h),]
bracken.table.phyl.diffabund.4h <- bracken.table.phyl.diffabund.4h %>% as.data.frame() %>% rownames_to_column() %>% mutate(species_id = rowname) %>% select(-rowname) %>% relocate(species_id)
bracken.table.phyl.diffabund.4h <- bracken.table.phyl.diffabund.4h %>% mutate_at(vars(-("species_id")), as.integer) %>% mutate_at("species_id", as.factor)

# Get diffabund table for 4hrs
bracken.table.phyl.diffabund.21h <- braken.table.phyl %>% left_join(metadata, by = 'Sample') %>% filter(Timepoint == "21h") %>% filter(Fraction == "BC_plus" | Fraction == "BC_minus") %>% select(-c("Fraction", "Desc_SampleID", "Bio_Rep",	"Tech_Rep", "Treatment", "Timepoint")) %>% relocate(Sample) %>% as.tibble() %>% t()
colnames(bracken.table.phyl.diffabund.21h) <- bracken.table.phyl.diffabund.21h[1,]
bracken.table.phyl.diffabund.21h <- bracken.table.phyl.diffabund.21h[2:nrow(bracken.table.phyl.diffabund.21h),]
bracken.table.phyl.diffabund.21h <- bracken.table.phyl.diffabund.21h %>% as.data.frame() %>% rownames_to_column() %>% mutate(species_id = rowname) %>% select(-rowname) %>% relocate(species_id)
bracken.table.phyl.diffabund.21h <- bracken.table.phyl.diffabund.21h %>% mutate_at(vars(-("species_id")), as.integer) %>% mutate_at("species_id", as.factor)


# Get metadata table for deseq2
deseq.md.4h <- metadata %>% filter(Timepoint=="4h") %>% filter(Fraction != "SYTO_plus")
deseq.md.21h <- metadata %>% filter(Timepoint=="21h") %>% filter(Fraction != "SYTO_plus")


library(DESeq2)
# Note: loading DESeq2 will mess up some tidyverse functions!
# Create DESeq2 dataset
dds.4h <- DESeqDataSetFromMatrix(countData=bracken.table.phyl.diffabund.4h, colData=deseq.md.4h, design = ~ Fraction, tidy=T)
dds.21h <- DESeqDataSetFromMatrix(countData=bracken.table.phyl.diffabund.21h, colData=deseq.md.21h, design = ~ Fraction, tidy=T)


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
tax <- bracken.table %>% select(c("Domain", "Phylum")) %>% distinct()
res.4h.boncat.plot <- res.4h.fraction %>% filter(padj <= 0.05) %>% mutate(Phylum = row) %>% select(-row) %>% left_join(tax, by = 'Phylum') %>% arrange(desc(log2FoldChange)) #%>% mutate(Genus=factor(Genus, levels=Genus))
res.21h.boncat.plot <- res.21h.fraction %>% filter(padj <= 0.05) %>% mutate(Phylum = row) %>% select(-row) %>% left_join(tax, by = 'Phylum') %>% arrange(desc(log2FoldChange)) #%>% mutate(Genus=factor(Genus, levels=Genus))


#dev.new()
p1 <- ggplot(res.4h.boncat.plot, aes(x = reorder(Phylum, log2FoldChange), y = log2FoldChange, fill = Phylum)) + geom_bar(stat='identity') + theme_bw() + theme(axis.title.y = element_blank()) + coord_flip() #+ theme(legend.position = "none")
p2 <- ggplot(res.4h.boncat.plot, aes(x = reorder(Phylum, log2FoldChange), y = padj)) + geom_point() + theme_bw() + theme(axis.text.y=element_blank(), axis.title.y = element_blank()) + coord_flip() + theme(legend.position = "none")
leg <- get_legend(p1) %>% as_ggplot()
plot <- ggarrange(p1, p2, leg, ncol=3, nrow=1, legend="none", widths = c(2, 0.25, 1))
plot
annotate_figure(plot, top = text_grob("Phylum-level Differential Abundance - 4hrs", face = "bold", size = 16))


p1 <- ggplot(res.21h.boncat.plot, aes(x = reorder(Phylum, log2FoldChange), y = log2FoldChange, fill = Phylum)) + geom_bar(stat='identity') + theme_bw() + theme(axis.title.y = element_blank()) + coord_flip() #+ theme(legend.position = "none")
p2 <- ggplot(res.21h.boncat.plot, aes(x = reorder(Phylum, log2FoldChange), y = padj)) + geom_point() + theme_bw() + theme(axis.text.y=element_blank(), axis.title.y = element_blank()) + coord_flip() + theme(legend.position = "none")
leg <- get_legend(p1) %>% as_ggplot()
plot <- ggarrange(p1, p2, leg, ncol=3, nrow=1, legend="none", widths = c(2, 0.25, 1))
plot
annotate_figure(plot, top = text_grob("Phylum-level Differential Abundance - 21hrs", face = "bold", size = 16))

#write_excel_csv(res.4h.boncat.plot,"manuscript_data_output/DESeq2_Phylum_4hr.csv")
#write_excel_csv(res.21h.boncat.plot,"manuscript_data_output/DESeq2_Phylum_21hr.csv")
