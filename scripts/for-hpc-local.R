## ============================================================
## Libraries
## ============================================================


library(plyr)
library(ggplot2)
library(phyloseq)
library(ggthemes)
library(vegan)
library(pheatmap)
library(pheatmap)
library(RColorBrewer)
library(data.table)

## ============================================================
## Heatmap of only significant WT vs OIL MAGs
## ============================================================

setwd('~/Box Sync/Oilcane/MAG-based/')

# Loading into a phyloseq object with QC filtering 
annot <- read.csv('all-mags-derep-renamed.gff.annot2', sep='\t', header=FALSE)
cog_list <- read.csv2(file="cog-20.def.tab.txt", sep="\t", header=FALSE)
annot2 <- merge(annot, cog_list, by.x = "V3", by.y = "V1", all.x = TRUE)
rownames(annot2) <- annot2$V1
annot2$V1 <- NULL
annot2 <- annot2[, colSums(is.na(annot2)) < nrow(annot2)]
colnames(annot2) <- c("COG", "genome", "COG_letters", "description", "gene", "descr2", 
                      "cog_num", "alphnum")
annot3 <- tax_table(as.matrix(annot2))
foo <- read.csv('summary-count.txt')
count <- foo
count <- count[!grepl("^__", rownames(count)),]
cn <- colnames(count)
length(cn)
cn2 <- cn[2:66]
cn3 <- strtrim(cn2, 9)
cn[2:66] <- cn3
cn
colnames(count) <- cn
rownames(count) <- count$Feature
count <- count[!grepl("^__", rownames(count)),]
count$Feature <- NULL
count_m <- as.matrix(count)
otu <- otu_table(count_m, taxa_are_rows = TRUE)
meta1 <- read.csv(file="../metadata_oilcane2.csv", header=TRUE, row.names = 1)
meta2 <- read.csv(file="../jgi-oilcane.csv", header=TRUE)
meta3 <- read.csv(file="../sample_key.txt", header=FALSE, sep="\t")
meta1$X
meta2$Seq.Project.Name
meta2$Seq.Project.Name <- gsub(" metagenome", "", meta2$Seq.Project.Name)
meta2$Seq.Project.Name
meta_merge1 <- merge(meta2, meta1, by.x = c("Seq.Project.Name"), by.y = c("sample_id"))
meta <- merge(meta_merge1, meta3, by.x=c("Seq.Project.Name"), by.y =  c("V2"))
head(meta)
head(meta)
rownames(meta) <- meta$V3
meta_phy <- sample_data(meta)
phy <- phyloseq(otu, meta_phy, annot3)
#phy_filter = prune_samples(sample_sums(phy)>=10000, phy)
phy_filter2 = phy
phy_filter2
phy_filter3 <- subset_samples(phy_filter2, Original.material..Soil..Roots..Stalks..Leaves. == "Root-associated soils")
phy_filter5 <- subset_samples(phy_filter3, Orig_name != "17T")
phy_filter3 <- phy_filter5

signif_genomes_list <- c("3300056388_1021", "3300056388_1321", "3300056388_133", "3300056388_1390",
                         "3300056388_1539", "3300056388_1866", "3300056388_1875", "3300056388_1894",
                         "3300056388_1898", "3300056388_1985", "3300056388_2063", "3300056388_2089",
                         "3300056388_2120", "3300056388_229", "3300056388_2568", "3300056388_2580",
                         "3300056388_2701", "3300056388_2924", "3300056388_3025", "3300056388_387",
                         "3300056388_390", "3300056388_398", "3300056388_638", "3300056388_874",
                         "3300056388_882")
phy_filter4 <- subset_taxa(phy_filter3, genome %in% c("3300056388_1021", "3300056388_1321", "3300056388_133", "3300056388_1390",
                                                      "3300056388_1539", "3300056388_1866", "3300056388_1875", "3300056388_1894",
                                                      "3300056388_1898", "3300056388_1985", "3300056388_2063", "3300056388_2089",
                                                      "3300056388_2120", "3300056388_229", "3300056388_2568", "3300056388_2580",
                                                      "3300056388_2701", "3300056388_2924", "3300056388_3025", "3300056388_387",
                                                      "3300056388_390", "3300056388_398", "3300056388_638", "3300056388_874",
                                                      "3300056388_882"))


ps <- psmelt(phy_filter4)
head(ps)
bin_info <- read.csv('exported_img_data.tsv', sep="\t")
ps_merged <- merge(ps, bin_info, by.x = "genome", by.y = "bin_oid")

x <- scan("./list-of-all-hkgs-in-mags.uniq.list2", what="", sep="\n")
hkg_phy <- prune_taxa(x, phy_filter2)
phy_melt_hkg = psmelt(hkg_phy)
f_hkg <- ddply(phy_melt_hkg, .(Sample), summarise, SUM=sum(Abundance))
ps_merged2 <- merge(ps_merged, f_hkg, by.x = "Sample", by.y = "Sample")
ps_merged2$norm <- ps_merged2$Abundance/ps_merged2$SUM

mat <- read.csv('matrix_data3.csv')
rownames(mat) <- mat$X
mat$X <- NULL
mat2 <- mat[,-c(4,5,6,22, 23, 24)]
matrix_data4 <- subset(mat2, rownames(mat) %in% signif_genomes_list )
cn_updated_man <- c("WT (CP88-1762) GH", "WT (CP88-1762) GH", "WT (CP88-1762) GH", "1565 GH", "1565 GH", "1565 GH", "1566 GH", "1566 GH", "1566 GH", 
                    "1569 GH", "1569 GH", "1569 GH", "1580 GH", "1580 GH", "1580 GH", "WT (CP88-1762) F", 
                    "WT (CP88-1762) F", "WT (CP88-1762) F", "1565 F", 
                    "1565 F", "1565 F", "1566 F", "1566 F", "1566 F", "1569 F", "1569 F", "1569 F", 
                    "1580 F", "1580 F", "1580 F")
colnames(matrix_data4) <- cn_updated_man
rn <- rownames(matrix_data4)
df1 <- as.data.frame(rn)
mag_ann <- read.csv(file="exported_img_data.tsv", sep="\t")
df2 <- merge(df1, mag_ann, by.x = "rn", by.y = "bin_oid")
annotation_r <-as.data.frame(df2$T1)
rownames(annotation_r) <- rn
colnames(annotation_r) <- "T1"

t1_colors <- c(  
  "#E69F00", # orange
  "#56B4E9", # sky blue
  "#009E73", # bluish green
  "#F0E442", # yellow
  "#0072B2", # blue
  "#D55E00", # vermillion
  "#CC79A7", # reddish purple
  "#999999", # gray
  "#000000", # black
  "#8B4513", # saddle brown (extra)
  "#800080"  # purple (extra)
  
)

list_t1 <- unique(annotation_r$T1)
taxon_col <- list(T1=setNames(t1_colors, list_t1))

out <- pheatmap(matrix_data4, color = brewer.pal(10, 'RdYlGn'), cluster_cols = FALSE,  clustering_distance_rows = 'euclidean', clustering_method = 'average', fontsize_row = 3, annotation_row = annotation_r, annotation_colors=taxon_col)


## ============================================================
## WT:Oilcane Ratio for each MAG (Figure 4)
## ============================================================


library(tidyr)

f <- ddply(ps_merged2, .(Sample, genome), summarise, MEAN=mean(Abundance), SE = sd(Abundance)/sqrt(length(Abundance)))
filtered_list <- unique(subset(f, MEAN > 5)$genome) 
ps_merged3 <- subset(ps_merged2, genome %in% filtered_list)

# Finding only significant genomes between WT and Other accessions
ps_WT_all <- ps_merged3 %>% mutate(is_WT = ifelse(Orig_name == "WT (CP88-1762)", "WT (CP88-1762)", "Oilcane"))

# Calculates normalized Coverage per Genome in Each Sample
ps_WT_all_T1 <- subset(ps_WT_all, Sampling.time == "T1")
ps_WT_norm <- ddply(ps_WT_all_T1, .(Sample, genome, is_WT), summarise, MEAN=mean(norm), SE=sd(norm)/sqrt(length(norm))) #This makes sense, you should be comparing the 'norm' coverage


# Calculates normalized Coverage per Genome in Each Sample
ps_WT_all_T2 <- subset(ps_WT_all, Sampling.time == "T2")
ps_WT_norm <- ddply(ps_WT_all_T2, .(Sample, genome, is_WT), summarise, MEAN=mean(norm), SE=sd(norm)/sqrt(length(norm))) #This makes sense, you should be comparing the 'norm' coverage


t1_colors <- brewer.pal(n = 10, name = "Set3")   # Good for taxonomy

# This is the summary of the genomes
ps_WT_subset <- subset(ps_WT_all, genome %in%  signif_genomes_list)
f <- ddply(ps_WT_subset, .(is_WT, Sample, genome, Sampling.time,  T1), summarise, MEAN=mean(norm), SE=sd(norm)/sqrt(length(norm)))
f2 <- ddply(f, .(is_WT, genome, T1, Sampling.time), summarise, MEAN2=mean(MEAN), SE2=mean(SE))
f3 <- f2
unique_genomes <- unique(f3$genome)

mags <- mag_ann
library(tidyverse)
wide_data <- f3 %>% select(is_WT, Sampling.time, genome, MEAN2) %>%
  pivot_wider(names_from = is_WT, values_from = MEAN2)
wide_data$ratio <- log2(wide_data$`WT (CP88-1762)`/wide_data$Oilcane)
wide_data_T1 <- subset(wide_data, Sampling.time == "T1")
ratio_sorted <- wide_data_T1[order(-wide_data$ratio),]
ratio_sorted$genome <- factor(ratio_sorted$genome, levels = ratio_sorted$genome[order(-ratio_sorted$ratio)],)
ratio_sorted2 <- merge(wide_data, mags, by.x = "genome", by.y = "bin_oid")

t1_colors <- c(
  
  "#E69F00", # orange
  "#56B4E9", # sky blue
  "#009E73", # bluish green
  "#F0E442", # yellow
  "#0072B2", # blue
  "#D55E00", # vermillion
  "#CC79A7", # reddish purple
  "#999999", # gray
  "#000000", # black
  "#8B4513", # saddle brown (extra)
  "#800080"  # purple (extra)
  
  )

list_t1 <- unique(ratio_sorted2$T1)
taxon_col <- setNames(t1_colors, list_t1)

ss <- subset(ratio_sorted2, Sampling.time == "T2")
ord <- order(ss$ratio)
ord_gen <- ss[ord,]$genome

#This produces the main image
ratio_sorted2$genome <- factor(ratio_sorted2$genome, levels=ord_gen)
p = ggplot(ratio_sorted2, aes(x=genome, y=ratio, fill=T1, shape=Sampling.time))
p = p+geom_point(size=5) +theme_bw()+theme(axis.text.x = element_text(angle = 90))+ylab("Log2 Ratio of WT:Oilcane") + xlab("Accession")
p = p+theme(strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5), axis.text.x = element_text(vjust=0.5, hjust=1, angle=90, size=14))+scale_fill_manual(values=taxon_col)
p = p+ theme(axis.title.x = element_text(size = 18),
             axis.title.y = element_text(size = 18)) +scale_shape_manual(values = c(21, 22))+scale_fill_manual(values=taxon_col)
p
ggsave(p, file="/Users/adina/Library/CloudStorage/GoogleDrive-adina.chuang@gmail.com/My Drive/Work Drive/Oilcane/Paper Figures and Tables/ratio1.pdf", width = 10)


#this produces the legend
ratio_sorted2$genome <- factor(ratio_sorted2$genome, levels=ord_gen)
p = ggplot(ratio_sorted2, aes(x=genome, y=ratio, fill=T1, shape=Sampling.time))
p = p+geom_point(size=5, aes(color=T1)) +theme_bw()+theme(axis.text.x = element_text(angle = 90))+ylab("Log2 Ratio of WT:Oilcane") + xlab("Accession")
p = p+theme(strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5), axis.text.x = element_text(vjust=0.5, hjust=1, angle=90, size=14))+scale_fill_manual(values=taxon_col)
p = p+ theme(axis.title.x = element_text(size = 18),
             axis.title.y = element_text(size = 18), legend.text=element_text(size=18)) +scale_shape_manual(values = c(21, 22))+scale_color_manual(values=taxon_col)
p
ggsave(p, file="/Users/adina/Library/CloudStorage/GoogleDrive-adina.chuang@gmail.com/My Drive/Work Drive/Oilcane/Paper Figures and Tables/ratio2.pdf", width = 10)

subset(ratio_sorted2, ratio > 0)$genome
more2 <- unique(subset(ratio_sorted2, ratio > 0)$genome)
less2 <- unique(subset(ratio_sorted2, ratio < 0)$genome)
unique(subset(ratio_sorted2, genome %in% more2 & Sampling.time == "T1"))
unique(subset(ratio_sorted2, genome %in% less2 & Sampling.time == "T1"))


## ============================================================
## WT:Oilcane Ratio for each MAG (Greenhouse Vs Field) (Supp Figure 1)
## ============================================================

head(f3)
wide <- f3 %>%
  pivot_wider(
    id_cols = c(is_WT, genome, T1),     # drop taxon if not present
    names_from = Sampling.time,
    values_from = c(MEAN2, SE2)
  )

wide2 <- wide %>%
  mutate(
    ratio_T2_to_T1  = MEAN2_T2 / MEAN2_T1,
    log2FC_T2_vs_T1 = log2(ratio_T2_to_T1)
  )


wide2$genome <- factor(wide2$genome, levels=ord_gen)

p = ggplot(wide2, aes(x=genome, y=log2FC_T2_vs_T1, color = T1, fill=T1, shape=is_WT))
p = p+geom_point(size=5) #+theme_bw()+theme(axis.text.x = element_text(angle = 90))+ylab("Log2 Ratio of WT:Oilcane") + xlab("Accession")
p = p +ylab("Log2 Ratio of Abundances in Field:Greenhouse") + xlab("MAG ID")+ theme_bw()
p = p+theme(strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5), axis.text.x = element_text(vjust=0.5, hjust=1, angle=90, size=14))+scale_color_manual(values=taxon_col)
p 

ggsave(p, file="/Users/adina/Library/CloudStorage/GoogleDrive-adina.chuang@gmail.com/My Drive/Work Drive/Oilcane/Paper Figures and Tables/ratio_t1_t2.pdf", width = 10)

##============================================================
## WT:Oilcane Ratio for each MAG (Greenhouse Vs Field) Faceted by Accession (Supp Figure 3)
## ============================================================

ps_WT_subset <- subset(ps_WT_all, genome %in%  signif_genomes_list)
f <- ddply(ps_WT_subset, .(Orig_name, Sample, genome, Sampling.time,  T1), summarise, MEAN=mean(norm), SE=sd(norm)/sqrt(length(norm)))
f2 <- ddply(f, .(Orig_name, genome, T1, Sampling.time), summarise, MEAN2=mean(MEAN), SE2=mean(SE))
wide <- f2 %>%
  pivot_wider(
    id_cols = c(Orig_name, genome, T1),     # drop taxon if not present
    names_from = Sampling.time,
    values_from = c(MEAN2, SE2)
  )
wide2 <- wide %>%
  mutate(
    ratio_T2_to_T1  = MEAN2_T2 / MEAN2_T1,
    log2FC_T2_vs_T1 = log2(ratio_T2_to_T1)
  )

p = ggplot(wide2, aes(x=genome, y=log2FC_T2_vs_T1, color=T1))
p = p+geom_point(size=5) #+theme_bw()+theme(axis.text.x = element_text(angle = 90))+ylab("Log2 Ratio of WT:Oilcane") + xlab("Accession")
p = p + theme_bw()+theme(axis.text.x = element_text(angle = 90))+ylab("Log2 Ratio of Abundances in Field:Greenhouse") + xlab("MAG")+facet_wrap(~Orig_name, nrow=5)+scale_color_manual(values=taxon_col)
p
ggsave(p, file="/Users/adina/Library/CloudStorage/GoogleDrive-adina.chuang@gmail.com/My Drive/Work Drive/Oilcane/Paper Figures and Tables/ratio_t1_t2_per_acc.pdf", width = 10)


## ============================================================
## Abundances of oilcane- or WT-enriched metagenome-assembled genomes (Supp Figure 2)
## ============================================================

# This is the summary of the genomes

# #Figure in Paper - WT Enriched MAG Abundances
ps_WT_subset <- subset(ps_WT_all, genome %in%  c(more2, less2))

oil <- less2
wt <- more2

ps_WT_subset$enrich <- ifelse(ps_WT_subset$genome %in% more2, "WT", "OIL")

ps_targets <- subset(ps_WT_all, genome %in% c(oil, wt))


plant_col <- c(
  "#E69F00", # Orange
  "#56B4E9", # Sky Blue
  "#009E73", # Bluish Green
  "#0072B2", # Blue
  "#D55E00")  # Vermillion
p_list <- c("WT (CP88-1762)", "1565", "1566", "1569", "1580")
plant_col2 <- setNames(plant_col, p_list)


f <- ddply(ps_WT_subset, .(genome, Sampling.time, T1, Sample, Orig_name, enrich), summarise, MEAN=mean(norm)) #Avg coverage of a genome
f2 <- ddply(f, .(genome, Sampling.time, Orig_name, T1, enrich), summarise, MEAN2=mean(MEAN), SE2=sd(MEAN)/sqrt(length(MEAN))) #Average of the observation of genomes
limits<-aes(ymin=MEAN2-SE2, ymax=MEAN2+SE2)
ggplot(f2, aes(x = Orig_name, y = MEAN2, shape=Sampling.time, color=Orig_name)) +
  geom_point(size=5)+theme_minimal()+geom_errorbar(limits, width=0)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_color_manual(values=plant_col2)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14))+ theme(strip.text = element_text(size = 20),
        legend.title = element_text(size = 14))+facet_wrap(enrich~genome, nrow=5, scales="free_y")+xlab('MAG ID')+ylab('Abundance (normalized)')
ggsave(file="/Users/adina/Library/CloudStorage/GoogleDrive-adina.chuang@gmail.com/My Drive/Work Drive/Oilcane/Paper Figures and Tables/per-genome.pdf", height=20, width = 30)


#Checking the abundance distribution in replicates
f <- ddply(ps_WT_subset, .(genome, Sample, Orig_name, Sampling.time), summarise, MEAN=mean(norm), SE=sd(norm)/sqrt(length(norm))) #Avg coverage of a genome
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
sample_ids <- unique(ps$Sample)
items <- sample_ids
random_colors <- rgb(runif(length(items)), runif(length(items)), runif(length(items)))
names(random_colors) <- sample_ids
ggplot(f, aes(x = Orig_name, y = MEAN, shape=Sampling.time, color=Orig_name)) +
  geom_point(size=3)+geom_errorbar(limits, width=0)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))+facet_wrap(~genome, nrow=5, scales="free_y")+xlab('WT-enriched MAG')+ylab('Abundance (normalized)')



#Figure in Paper - WT Enriched MAG Abundances

plant_col <- c(
  "#E69F00", # Orange
  "#56B4E9", # Sky Blue
  "#009E73", # Bluish Green
  "#0072B2", # Blue
  "#D55E00")  # Vermillion
p_list <- c("WT (CP88-1762)", "1565", "1566", "1569", "1580")
plant_col2 <- setNames(plant_col, p_list)


ps_WT_subset <- subset(ps_WT_all, genome %in%  less2)
f <- ddply(ps_WT_subset, .(genome, Sampling.time, T1, Sample, Orig_name), summarise, MEAN=mean(norm)) #Avg coverage of a genome
f2 <- ddply(f, .(genome, Sampling.time, Orig_name, T1), summarise, MEAN2=mean(MEAN), SE2=sd(MEAN)/sqrt(length(MEAN))) #Average of the observation of genomes 
limits<-aes(ymin=MEAN2-SE2, ymax=MEAN2+SE2)
f2$Orig_name <- factor(f2$Orig_name, levels= p_list)
f2$genome <- factor(f2$genome, levels = ord_gen[1:19])
                       
ggplot(f2, aes(x = Orig_name, y = MEAN2, shape=Sampling.time, color=T1)) +
  geom_point(size=5)+theme_minimal()+geom_errorbar(limits, width=0)+
  theme(strip.text.x = element_text(angle = 90, hjust = 1, size=12)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))+facet_wrap(~genome, nrow=2, scales="free_y")+xlab('WT-enriched MAG')+ylab('Abundance')+
  scale_color_manual(values=taxon_col)
ggsave(file='less.png', dpi=300, width=30)
# 


t1_colors <- c(
  
  "#E69F00", # orange
  "#56B4E9", # sky blue
  "#009E73", # bluish green
  "#F0E442", # yellow
  "#0072B2", # blue
  "#D55E00", # vermillion
  "#CC79A7", # reddish purple
  "#999999", # gray
  "#000000", # black
  "#8B4513", # saddle brown (extra)
  "#800080"  # purple (extra)
  
)

list_t1 <- unique(ratio_sorted2$T1)
taxon_col <- setNames(t1_colors, list_t1)
ps_WT_subset <- subset(ps_WT_all, genome %in%  more2)
f <- ddply(ps_WT_subset, .(genome, Sampling.time, T1, Sample, Orig_name), summarise, MEAN=mean(norm)) #Avg coverage of a genome
f2 <- ddply(f, .(genome, Sampling.time, Orig_name, T1), summarise, MEAN2=mean(MEAN), SE2=sd(MEAN)/sqrt(length(MEAN))) #Average of the observation of genomes 
limits<-aes(ymin=MEAN2-SE2, ymax=MEAN2+SE2)
f2$Orig_name <- factor(f2$Orig_name, levels= p_list)
f2$genome <- factor(f2$genome, levels = ord_gen[20:25])

ggplot(f2, aes(x = Orig_name, y = MEAN2, shape=Sampling.time, color=T1)) +
  geom_point(size=5)+theme_minimal()+geom_errorbar(limits, width=0)+
  theme(strip.text.x = element_text(angle = 90, hjust = 1, size=12)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))+facet_wrap(~genome, nrow=2, scales="free_y")+xlab('WT-enriched MAG')+ylab('Abundance')+
  scale_color_manual(values=taxon_col)


# #Checking the abundance distribution in replicates
f <- ddply(ps_WT_subset, .(genome, Sample, Orig_name, Sampling.time), summarise, MEAN=mean(norm), SE=sd(norm)/sqrt(length(norm))) #Avg coverage of a genome
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
sample_ids <- unique(ps$Sample)
items <- sample_ids
random_colors <- rgb(runif(length(items)), runif(length(items)), runif(length(items)))
names(random_colors) <- sample_ids
ggplot(f, aes(x = Orig_name, y = MEAN, shape=Sampling.time, color=Orig_name)) +
  geom_point(size=3)+geom_errorbar(limits, width=0)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 16),
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))+facet_wrap(~genome, nrow=3, scales="free_y")+xlab('WT-enriched MAG')+ylab('Abundance (normalized)')


# 
# 
# 
# #Figure in Paper - WT Enriched Functions
# # All the COGs split into their V5 categories in WT vs Others

ps_WT_subset <- subset(ps_WT_all, genome %in%  more2)
ps_WT_sig <- ps_WT_subset
cog1 <- read.csv(file = "cog-20.def.tab.txt", header=FALSE, sep="\t")
cog2 <- read.csv(file = "fun-20.tab.txt", header=FALSE, sep="\t")
cog_merge <- merge(cog1, cog2, by.x = "V2", by.y = "V1")
ps_cog <- merge(ps_WT_sig, cog_merge, by.x = "COG", by.y = "V1")
f <- ddply(ps_cog, .(V3.y, genome, Sample, Orig_name), summarise, SUM=sum(norm))
f2 <- ddply(f, .(V3.y, genome, Orig_name), summarize, MEAN = mean(SUM), SE = sd(SUM)/sqrt(length(SUM)))
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
p = ggplot(f2, aes(x=Orig_name, y=MEAN, shape=Orig_name, color=genome))
#p = p+geom_bar(stat="identity", position=position_dodge2(preserve="single", width = 0.8)) +theme_bw()+theme(axis.text.x = element_text(angle = 90))
p = p+geom_point(size=2) +theme_bw()+theme(axis.text.x = element_text(angle = 90))+geom_errorbar(limits, width=0)
p+facet_grid(col=vars(V3.y), scales="free_x")+
  theme(axis.title.x = element_text(size=10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 10),
        strip.text.x = element_text(size=10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10), strip.text = element_text(angle = 90) )+ylab('Abundance (norm)')+ xlab('Functional Genes')

ps_WT_subset <- subset(ps_WT_all, genome %in%  c(more2, less2))
ps_WT_sig <- ps_WT_subset
cog1 <- read.csv(file = "cog-20.def.tab.txt", header=FALSE, sep="\t")
cog2 <- read.csv(file = "fun-20.tab.txt", header=FALSE, sep="\t")
cog_merge <- merge(cog1, cog2, by.x = "V2", by.y = "V1")
ps_cog <- merge(ps_WT_sig, cog_merge, by.x = "COG", by.y = "V1")
f <- ddply(ps_cog, .(V3.y, genome, Sample, Orig_name), summarise, SUM=sum(norm))
f2 <- ddply(f, .(V3.y, genome, Orig_name), summarize, MEAN = mean(SUM), SE = sd(SUM)/sqrt(length(SUM)))
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
p = ggplot(f2, aes(x=Orig_name, y=MEAN, shape=Orig_name, color=genome))
#p = p+geom_bar(stat="identity", position=position_dodge2(preserve="single", width = 0.8)) +theme_bw()+theme(axis.text.x = element_text(angle = 90))
p = p+geom_point(size=2) +theme_bw()+theme(axis.text.x = element_text(angle = 90))+geom_errorbar(limits, width=0)
p+facet_grid(col=vars(V3.y), scales="free_x")+
  theme(axis.title.x = element_text(size=10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 10),
        strip.text.x = element_text(size=10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10), strip.text = element_text(angle = 90) )+ylab('Abundance (norm)')+ xlab('Functional Genes')


# f2_hunt <- subset(f2, genome == "3300056388_2568")
# p = ggplot(f2_hunt, aes(x=Orig_name, y=MEAN, shape=Orig_name, color=genome))        
# #p = p+geom_bar(stat="identity", position=position_dodge2(preserve="single", width = 0.8)) +theme_bw()+theme(axis.text.x = element_text(angle = 90))
# p = p+geom_point(size=2) +theme_bw()+theme(axis.text.x = element_text(angle = 90))+geom_errorbar(limits, width=0)
# p+facet_grid(col=vars(V3.y), scales="free_x")+
#   theme(axis.title.x = element_text(size=10),
#         axis.text.y = element_text(size = 10),
#         axis.title = element_text(size = 10),
#         strip.text.x = element_text(size=10),
#         legend.text = element_text(size = 10),
#         legend.title = element_text(size = 10), strip.text = element_text(angle = 90) )+ylab('Abundance (norm)')+ xlab('Functional Genes')
# 
# 
# 
# ps_WT_subset <- subset(ps_WT_all, genome %in%  more2)
# ps_WT_sig <- ps_WT_subset
# cog1 <- read.csv(file = "cog-20.def.tab.txt", header=FALSE, sep="\t")
# cog2 <- read.csv(file = "fun-20.tab.txt", header=FALSE, sep="\t")
# cog_merge <- merge(cog1, cog2, by.x = "V2", by.y = "V1")
# ps_cog <- merge(ps_WT_sig, cog_merge, by.x = "COG", by.y = "V1")
# f <- ddply(ps_cog, .(V3.y, genome, Sample, Orig_name), summarise, SUM=sum(norm))
# f2 <- ddply(f, .(V3.y, genome, Orig_name), summarize, MEAN = mean(SUM), SE = sd(SUM)/sqrt(length(SUM)))
# f3 <- ddply(f2, .(V3.y, Orig_name), summarize, MEAN2 = mean(MEAN) , SE2 = sd(MEAN)/sqrt(length(MEAN)))
# limits<-aes(ymin=MEAN2-SE2, ymax=MEAN2+SE2)
# p = ggplot(f3, aes(x=Orig_name, y=MEAN2, color=Orig_name))        
# #p = p+geom_bar(stat="identity", position=position_dodge2(preserve="single", width = 0.8)) +theme_bw()+theme(axis.text.x = element_text(angle = 90))
# p = p+geom_point(size=2) +theme_bw()+theme(axis.text.x = element_text(angle = 90))+geom_errorbar(limits, width=0)
# p+facet_grid(col=vars(V3.y), scales="free_x")+
#   theme(axis.title.x = element_text(size=10),
#         axis.text.y = element_text(size = 10),
#         axis.title = element_text(size = 10),
#         strip.text.x = element_text(size=10),
#         legend.text = element_text(size = 10),
#         legend.title = element_text(size = 10), strip.text = element_text(angle = 90) )+ylab('Abundance (norm)')+ xlab('Functional Genes')
# 
# 
# # 
# ps_WT_subset <- subset(ps_WT_all, genome %in%  less2)
# ps_WT_sig <- ps_WT_subset
# cog1 <- read.csv(file = "cog-20.def.tab.txt", header=FALSE, sep="\t")
# cog2 <- read.csv(file = "fun-20.tab.txt", header=FALSE, sep="\t")
# cog_merge <- merge(cog1, cog2, by.x = "V2", by.y = "V1")
# ps_cog <- merge(ps_WT_sig, cog_merge, by.x = "COG", by.y = "V1")
# f <- ddply(ps_cog, .(V3.y, genome, Sample, is_WT), summarise, SUM=sum(norm))
# f2 <- ddply(f, .(V3.y, genome, is_WT), summarize, MEAN = mean(SUM), SE = sd(SUM)/sqrt(length(SUM)))
# limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
# p = ggplot(f2, aes(x=is_WT, y=MEAN, shape=is_WT, color=genome))        
# #p = p+geom_bar(stat="identity", position=position_dodge2(preserve="single", width = 0.8)) +theme_bw()+theme(axis.text.x = element_text(angle = 90))
# p = p+geom_point(size=2) +theme_bw()+theme(axis.text.x = element_text(angle = 90))+geom_errorbar(limits, width=0)
# p+facet_grid(col=vars(V3.y), scales="free_x")+
#   theme(axis.title.x = element_text(size=10),
#         axis.text.y = element_text(size = 10),
#         axis.title = element_text(size = 10),
#         strip.text.x = element_text(size=10),
#         legend.position = "none",
#         legend.title = element_text(size = 10), strip.text = element_text(angle = 90) )+ylab('Abundance (norm)')+ xlab('Functional Genes')
# 
# ps_WT_subset <- subset(ps_WT_all, genome %in%  more2)
# ps_WT_sig <- ps_WT_subset
# cog1 <- read.csv(file = "cog-20.def.tab.txt", header=FALSE, sep="\t")
# cog2 <- read.csv(file = "fun-20.tab.txt", header=FALSE, sep="\t")
# cog_merge <- merge(cog1, cog2, by.x = "V2", by.y = "V1")
# ps_cog <- merge(ps_WT_sig, cog_merge, by.x = "COG", by.y = "V1")
# ps_cog2 <- subset(ps_cog, V5 != "")
# f <- ddply(ps_cog2, .(COG, V3.y, genome, Sample,Orig_name), summarise, SUM=sum(norm)) #sums COGs over ORFs
# f2 <- ddply(f, .(COG, V3.y, genome, Orig_name), summarize, MEAN = mean(SUM), SE = sd(SUM)/sqrt(length(SUM))) #averages over the samples
# limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
# p = ggplot(f2, aes(x=COG, y=MEAN, fill=genome))        
# p = p+geom_bar(stat="identity") +theme_bw()+theme(axis.text.x = element_text(angle = 90))
# p+facet_grid(row=vars(Orig_name), col=vars(V3.y), scales="free_x")+
#   theme(axis.text.x = element_text(size = 5, vjust=0.5),
#         axis.text.y = element_text(size = 15),
#         axis.title = element_text(size = 15),
#         legend.text = element_text(size = 15),
#         legend.title = element_text(size = 15), strip.text = element_text(angle = 90))
# 
# 
# ps_WT_subset <- subset(ps_WT_all, genome %in%  more2)
# ps_WT_sig <- ps_WT_subset
# cog1 <- read.csv(file = "cog-20.def.tab.txt", header=FALSE, sep="\t")
# cog2 <- read.csv(file = "fun-20.tab.txt", header=FALSE, sep="\t")
# cog_merge <- merge(cog1, cog2, by.x = "V2", by.y = "V1")
# ps_cog <- merge(ps_WT_sig, cog_merge, by.x = "COG", by.y = "V1")
# ps_cog2 <- subset(ps_cog, V5 != "")
# f <- ddply(ps_cog2, .(COG, V3.y, genome, T1, Sample,Orig_name), summarise, SUM=sum(norm)) #sums COGs over ORFs
# f2 <- ddply(f, .(COG, V3.y, genome, Orig_name, T1), summarize, MEAN = mean(SUM), SE = sd(SUM)/sqrt(length(SUM))) #averages over the samples
# limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
# p = ggplot(f2, aes(x=genome, y=MEAN, fill=T1))        
# p = p+geom_bar(stat="identity") +theme_bw()+theme(axis.text.x = element_text(angle = 90))
# p+facet_grid(row=vars(Orig_name), col=vars(V3.y), scales="free_x")+
#   theme(axis.text.x = element_text(size = 14, vjust=0.5),
#         axis.text.y = element_text(size = 15),
#         axis.title = element_text(size = 15),
#         legend.text = element_text(size = 15),
#         legend.title = element_text(size = 15), strip.text = element_text(angle = 90))
# 
# 
# ps_WT_subset <- subset(ps_WT_all, genome %in%  less2)
# ps_WT_sig <- ps_WT_subset
# cog1 <- read.csv(file = "cog-20.def.tab.txt", header=FALSE, sep="\t")
# cog2 <- read.csv(file = "fun-20.tab.txt", header=FALSE, sep="\t")
# cog_merge <- merge(cog1, cog2, by.x = "V2", by.y = "V1")
# ps_cog <- merge(ps_WT_sig, cog_merge, by.x = "COG", by.y = "V1")
# ps_cog2 <- subset(ps_cog, V5 != "")
# f <- ddply(ps_cog2, .(COG, V3.y, genome, T1, Sample,Orig_name), summarise, SUM=sum(norm)) #sums COGs over ORFs
# f2 <- ddply(f, .(COG, V3.y, genome, Orig_name, T1), summarize, MEAN = mean(SUM), SE = sd(SUM)/sqrt(length(SUM))) #averages over the samples
# limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
# p = ggplot(f2, aes(x=genome, y=MEAN, fill=T1))        
# p = p+geom_bar(stat="identity") +theme_bw()+theme(axis.text.x = element_text(angle = 90))
# p+facet_grid(row=vars(Orig_name), col=vars(V3.y), scales="free_x")+
#   theme(axis.text.x = element_text(size = 14, vjust=0.5),
#         axis.text.y = element_text(size = 15),
#         axis.title = element_text(size = 15),
#         legend.position = "none",
#         legend.title = element_text(size = 15), strip.text = element_text(angle = 90))


#targets <- subset(ps_WT_all, genome %in% c(oil, wt))
targets <- subset(ps_WT_all, genome %in% c(more2, less2))
#targets <- subset(ps_WT_all, genome %in% c(less2))
#targets <- subset(ps_WT_all, genome %in% c(more2))
wt <- more2
oil <- less2
targets$group <- ifelse(targets$genome %in% oil, "oil", "wt")
#targets$group <- ifelse(targets$genome %in% less2, "less", "more")

#f <- ddply(targets, .(genome,  T1, Sample, Orig_name), summarise, MEAN=mean(norm)) #Avg coverage of a genome, average of all ORFs
#f2 <- ddply(f, .(genome, Orig_name, T1), summarise, MEAN2=mean(MEAN), SE2=sd(MEAN)/sqrt(length(MEAN))) #Average of the observation of genomes 
#limits<-aes(ymin=MEAN2-SE2, ymax=MEAN2+SE2)
f <- ddply(targets, .(genome,  T1, Sample, Orig_name, group), summarise, MEAN=mean(norm)) #Avg coverage of a genome, average of all ORFs

# Plot of genome in each accession
# ggplot(f2, aes(x = genome, y = MEAN2, color=T1)) +
#   geom_point()+theme_minimal()+geom_errorbar(limits, width=0)+
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
#         axis.text.y = element_text(size = 14),
#         axis.title = element_text(size = 14),
#         legend.text = element_text(size = 14),
#         strip.text.x = element_text(angle = 90),
#         legend.title = element_text(size = 14))+facet_wrap(~Orig_name, nrow=1, scales="free_x")+ylab("Coverage (normalized by HKG)")+xlab('Genome')
# 
f_wt <- subset(f, group == "wt")
f_oil <- subset(f, group == "oil")

library(data.table)
sample_sums_wt <- as.data.table(f_wt)[, .(sum_MEAN_wt = sum(MEAN, na.rm = TRUE)), by = Sample]
wt <- as.data.frame(sample_sums_wt)
sample_sums_oil <- as.data.table(f_oil)[, .(sum_MEAN_oil = sum(MEAN, na.rm = TRUE)), by = Sample]
oil <- as.data.frame(sample_sums_oil)



plant_col <- c(
  "#E69F00", # Orange
  "#56B4E9", # Sky Blue
  "#009E73", # Bluish Green
  "#0072B2", # Blue
  "#D55E00")  # Vermillion
p_list <- c("WT (CP88-1762)", "1565", "1566", "1569", "1580")
plant_col2 <- setNames(plant_col, p_list)

combo <- merge(wt, oil, by.x ="Sample", by.y = "Sample")
ggplot(combo, aes(x= sum_MEAN_wt, y = sum_MEAN_oil)) + geom_point()
meta1 <- read.csv(file="~/Box Sync/Oilcane/metadata_oilcane2.csv", header=TRUE, row.names = 1)
meta2 <- read.csv(file="~/Box Sync/Oilcane/jgi-oilcane.csv", header=TRUE)
meta3 <- read.csv(file="~/Box Sync/Oilcane/sample_key.txt", header=FALSE, sep="\t")
meta1$X
meta2$Seq.Project.Name
meta2$Seq.Project.Name <- gsub(" metagenome", "", meta2$Seq.Project.Name)
meta2$Seq.Project.Name
meta_merge1 <- merge(meta2, meta1, by.x = c("Seq.Project.Name"), by.y = c("sample_id"))
meta <- merge(meta_merge1, meta3, by.x=c("Seq.Project.Name"), by.y =  c("V2"))
head(meta)
combo2 <- merge(combo, meta, by.x = "Sample", by.y = "V3")
combo3 <- subset(combo2, sum_MEAN_wt < 0.0050)
combo_s <- subset(combo2, Orig_name == "WT (CP88-1762)")
combo2b <- subset(combo2, Sampling.time == "T2")
ggplot(combo2b, aes(x= sum_MEAN_wt, y = sum_MEAN_oil, color = Orig_name, shape=Sampling.time)) + geom_point(size=6)+
  geom_smooth(method = "lm", se = TRUE, color = "black") +xlab('Average abundance (norm) of WT')+
  scale_color_manual(values = plant_col2) +theme_bw()+xlab('Average abundance (norm) of WT-enriched MAGs')+ylab('Average abundance (norm) of oilcane-enriched MAGs')
model <- lm(sum_MEAN_oil ~ sum_MEAN_wt, data = combo2b)
summary(model)
ggsave(file="/Users/adina/Library/CloudStorage/GoogleDrive-adina.chuang@gmail.com/My Drive/Work Drive/Oilcane/Paper Figures and Tables/t1-model.pdf", height=10, width = 10)

library(lme4)
model <- lmer(sum_MEAN_oil ~ sum_MEAN_wt + (1 | Sampling.time), data = combo2)
ggplot(combo3, aes(x= sum_MEAN_wt, y = sum_MEAN_oil, color=Orig_name)) + geom_point(size=6)+
  geom_smooth(method = "lm", se = TRUE, color = "black") +xlab('Average abundance (norm) of WT-enriched MAGs')+ylab('Average abundance (norm) of oilcane-enriched MAGs')
model <- lm(sum_MEAN_oil ~ sum_MEAN_wt, data = combo3)
summary(model)


library(lmerTest)
model <- lmer(sum_MEAN_oil ~ sum_MEAN_wt + (1 | Sampling.time), data = combo2)
summary(model)

ggplot(combo2, aes(x= sum_MEAN_wt, y = sum_MEAN_oil, shape=Sampling.time, color=Orig_name)) + geom_point(size=6)+
  geom_smooth(method = "lm", se=FALSE, linetype='dotted') +xlab('Average abundance (norm) of WT-enriched MAGs')+ylab('Average abundance (norm) of oilcane-enriched MAGs')
model <- lm(sum_MEAN_oil ~ sum_MEAN_wt, data = combo2)
summary(model)


list_t1 <- unique(ratio_sorted2$T1)
taxon_col <- setNames(t1_colors, list_t1)

oil <- less2
wt <- more2
targets <- subset(ps_WT_all, genome %in% c(oil, wt))
targets$group <- ifelse(targets$genome %in% wt, "WT", "oilcane")

f <- ddply(targets, .(genome,  T1, Sample, Orig_name, group), summarise, MEAN=mean(norm)) #Avg coverage of a genome
f2 <- ddply(f, .(genome, Orig_name, T1,  group), summarise, MEAN2=mean(MEAN), SE2=sd(MEAN)/sqrt(length(MEAN))) #Average of the observation of genomes 
list_t1 <- unique(f2$T1)
taxon_col <- setNames(t1_colors, list_t1)
limits<-aes(ymin=MEAN2-SE2, ymax=MEAN2+SE2)
ggplot(f2, aes(x = T1, y = MEAN2, color=T1)) +
  geom_point(size=3)+theme_minimal()+geom_errorbar(limits, width=0)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        strip.text.x = element_text(angle = 90),
        legend.title = element_text(size = 14))+facet_wrap(group~Orig_name, nrow=2, scales="free_x")+ylab("Coverage (normalized by HKG)")+xlab('Genome')+
  scale_color_manual(values=taxon_col)

f2$plant <- ifelse(f2$Orig_name != "WT (CP88-1762)", "OIL", "WT")
f3 <- ddply(f2, .(plant, T1, genome, group), summarise, MEAN3=mean(MEAN2), SE3 =sd(MEAN2)/sqrt(length(MEAN2)) ) #Average of the observation of genomes 
list_t1 <- unique(f3$T1)
limits<-aes(ymin=MEAN3-SE3, ymax=MEAN3+SE3)
taxon_col <- setNames(t1_colors, list_t1)
ggplot(f3, aes(x = T1, y = MEAN3, color=plant)) +
  geom_point(size=5)+theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 20),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        strip.text.x = element_text(angle = 90),
        legend.title = element_text(size = 20))+ylab("Mean Normalized Coverage")+xlab('Genome')+facet_wrap(~group, scale="free_x")+
  scale_color_manual(values=c("red", "black"))

# there are going to be variations of abundances of each ORF, which we accoutn for in teh average genome abundance
# we should go back and overlay all functions with that average abundane
# first we need to name each genome_sample pairing
targets$name2 <- paste0(targets$Sample, "_", targets$genome) 
f$name2 <- paste0(f$Sample, "_", f$genome) 

targets2 <- targets %>% select(Sample, genome, OTU,  Orig_name, Sampling.time, name2, COG)
merged_abund <- merge(targets2, f, by = "name2") #so now there is avg mean assigned to each abundance, MEAN
targets_cog <- merged_abund %>% select(Sample.y, genome.y, OTU, MEAN, Orig_name.y, T1, COG)
#This is the averaged abundance of genoem overlaid on each function and then used to make the functional plots
#targets_cog2 <- subset(targets_cog, genome.y %in% more2)
cog1 <- read.csv(file = "cog-20.def.tab.txt", header=FALSE, sep="\t")
cog2 <- read.csv(file = "fun-20.tab.txt", header=FALSE, sep="\t")
cog_merge <- merge(cog1, cog2, by.x = "V2", by.y = "V1")
ps_cog <- merge(targets_cog, cog_merge, by.x = "COG", by.y = "V1")
ps_cog2 <- subset(ps_cog, V3.y != "")
f <- ddply(ps_cog2, .(COG, V3.y, genome.y, T1, Sample.y,Orig_name.y), summarise, SUM=sum(MEAN)) #sums COGs over ORFs
f2 <- ddply(f, .(V3.y, genome.y, Orig_name.y, T1, Sample.y), summarize, SUM2 = sum(SUM)) 
f2b <- subset(f2, genome.y != "3300056388_325")
f3 <- ddply(f2b, .(V3.y, Orig_name.y, genome.y, T1), summarize, MEAN = mean(SUM2) , SE = sd(SUM2)/sqrt(length(SUM2)))
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
p = ggplot(f3, aes(x=Orig_name.y, y=MEAN, color=T1))        
p = p+geom_point(size=3) +theme_bw()+theme(axis.text.x = element_text(angle = 90))+geom_errorbar(limits, width=0)
p+facet_grid(col=vars(V3.y), scales="free_x")+
  theme(axis.title.x = element_text(size=15),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),            
        axis.title = element_text(size = 15),
        strip.text.x = element_text(size=15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15), strip.text = element_text(angle = 90) )+ylab('Abundance')+ xlab('Functional Genes')



t1_colors <- c(
  
  "#E69F00", # orange
  "#56B4E9", # sky blue
  "#009E73", # bluish green
  "#F0E442", # yellow
  "#0072B2", # blue
  "#D55E00", # vermillion
  "#CC79A7", # reddish purple
  "#999999", # gray
  "#000000", # black
  "#8B4513", # saddle brown (extra)
  "#800080"  # purple (extra)
  
)
list_t1 <- unique(f3$T1.y)
taxon_col <- setNames(t1_colors, list_t1)

# This is all functions unfiltered by what is unique
p = ggplot(f3, aes(x=Orig_name.y, y=MEAN, fill=T1))        
p = p+geom_bar(stat="identity", position="stack") +theme_bw()+theme(axis.text.x = element_text(angle = 90))
p+facet_grid(col=vars(V3.y), scales="free_x")+
  theme(axis.title.x = element_text(size=15),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.title = element_text(size = 15),
        strip.text.x = element_text(size=15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15), strip.text = element_text(angle = 90) )+ylab('Abundance')+ xlab('Functional Genes')+
        scale_fill_manual(values=taxon_col)

list_of_func <- c("Amino acid transport and metabolism", "Carbohydrate transport and metabolism", 
"Coenzyme transport and metabolism", "Defense mechanisms", 
"Energy production and conversion",  "Inorganic ion transport and metabolism", 
"Lipid transport and metabolism", 
"Secondary metabolites biosynthesis, transport and catabolism", 
"Signal transduction mechanisms", "Transcription", "Translation, ribosomal structure and biogenesis"
)

f4 <- subset(f3, V3.y %in% list_of_func)



df_rel <- f4 %>% select(V3.y, genome.y, Orig_name.y, T1, MEAN) %>%
  group_by(Orig_name.y) %>%
  mutate(
    MEAN_REL = MEAN / sum(MEAN)
  ) %>%
  ungroup()

# This is all functions unfiltered by what is unique
p = ggplot(f4, aes(x=Orig_name.y, y=MEAN, fill=T1))        
p = p+geom_bar(stat="identity", position="stack") +theme_bw()+theme(axis.text.x = element_text(angle = 90))
p+facet_grid(col=vars(V3.y), scales="free_x")+
  theme(axis.title.x = element_text(size=15),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.title = element_text(size = 15),
        strip.text.x = element_text(size=15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15), strip.text = element_text(angle = 90) )+ylab('Abundance (norm)')+ xlab('Functional Genes')+
        scale_fill_manual(values=taxon_col)
ggsave(file="/Users/adina/Library/CloudStorage/GoogleDrive-adina.chuang@gmail.com/My Drive/Work Drive/Oilcane/Paper Figures and Tables/function.pdf", height=10, width = 13)


# This is all functions unfiltered by what is unique
p = ggplot(df_rel, aes(x=Orig_name.y, y=MEAN_REL, fill=T1))        
p = p+geom_bar(stat="identity") +theme_bw()+theme(axis.text.x = element_text(angle = 90))
p+facet_grid(col=vars(V3.y), scales="free_x")+
  theme(axis.title.x = element_text(size=15),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.title = element_text(size = 15),
        strip.text.x = element_text(size=15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15), strip.text = element_text(angle = 90) )+ylab('Relative Abundance')+ xlab('Functional Genes')+
        scale_fill_manual(values=taxon_col)
ggsave(file="/Users/adina/Library/CloudStorage/GoogleDrive-adina.chuang@gmail.com/My Drive/Work Drive/Oilcane/Paper Figures and Tables/function_rel.pdf", height=10, width = 13)



 a <- subset(df_rel, genome.y == "3300056388_1021" & Orig_name.y == "1565")
 a <- subset(df_rel, Orig_name.y == "1566")
# Comparing Time Points and Functions
 
targets_cog <- merged_abund %>% select(Sample.y, genome.y, OTU, MEAN, Orig_name.y, T1, COG, Sampling.time)

#Oilcane samples consistently exhibited higher abundances of genes involved in energy production and conversion, lipid transport and metabolism, secondary metabolite biosynthesis, transport and catabolism and signal‑transduction mechanisms. 

list_of_func <- c(
"Energy production and conversion",  
"Lipid transport and metabolism", 
"Secondary metabolites biosynthesis, transport and catabolism", 
"Signal transduction mechanisms")





#This is the averaged abundance of genoem overlaid on each function and then used to make the functional plots
cog1 <- read.csv(file = "cog-20.def.tab.txt", header=FALSE, sep="\t")
cog2 <- read.csv(file = "fun-20.tab.txt", header=FALSE, sep="\t")
cog_merge <- merge(cog1, cog2, by.x = "V2", by.y = "V1")
ps_cog <- merge(targets_cog, cog_merge, by.x = "COG", by.y = "V1")
ps_cog2 <- subset(ps_cog, V3.y != "")
f <- ddply(ps_cog2, .(COG, V3.y, genome.y, T1, Sampling.time, Sample.y,Orig_name.y), summarise, SUM=sum(MEAN)) #sums COGs over ORFs
f2 <- ddply(f, .(V3.y, genome.y, Orig_name.y, T1, Sample.y, Sampling.time), summarize, SUM2 = sum(SUM)) 
f3 <- ddply(f2, .(V3.y, Orig_name.y, genome.y, T1, Sampling.time), summarize, MEAN = mean(SUM2) , SE = sd(SUM2)/sqrt(length(SUM2)))
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
p = ggplot(f3, aes(x=Orig_name.y, y=MEAN, fill=T1))        
p = p+geom_bar(stat="identity") +theme_bw()+theme(axis.text.x = element_text(angle = 90))
p+facet_grid(col=vars(V3.y), scales="free_x")+
  theme(axis.title.x = element_text(size=15),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),            
        axis.title = element_text(size = 15),
        strip.text.x = element_text(size=15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15), strip.text = element_text(angle = 90) )+ylab('Abundance')+ xlab('Functional Genes')





t1_colors <- c(
  
  "#E69F00", # orange
  "#56B4E9", # sky blue
  "#009E73", # bluish green
  "#F0E442", # yellow
  "#0072B2", # blue
  "#D55E00", # vermillion
  "#CC79A7", # reddish purple
  "#999999", # gray
  "#000000", # black
  "#8B4513", # saddle brown (extra)
  "#800080"  # purple (extra)
  
)
list_t1 <- unique(f3$T1.y)
taxon_col <- setNames(t1_colors, list_t1)

f3b <- subset(f3, Sampling.time == "T2")
# This is all functions unfiltered by what is unique
f4 <- subset(f3, V3.y %in% list_of_func)
p = ggplot(f4, aes(x=Orig_name.y, y=MEAN, fill=T1))        
p = p+geom_bar(stat="identity", position="stack") +theme_bw()+theme(axis.text.x = element_text(angle = 90))
p+facet_grid(V3.y~Sampling.time)+ scale_fill_manual(values=taxon_col)+
  theme(axis.title.x = element_text(size=15),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.title = element_text(size = 15),
        strip.text.x = element_text(size=15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15), strip.text = element_text(angle = 90) )+ylab('Abundance')
ggsave(file="/Users/adina/Library/CloudStorage/GoogleDrive-adina.chuang@gmail.com/My Drive/Work Drive/Oilcane/Paper Figures and Tables/function_focus.pdf", height=10, width = 13)


####################################################################################################
head(ps_cog2)
func_chosen = list_of_func[1]
ps_cog_select <- subset(ps_cog2, V3.y %in% func_chosen)
ps_cog_select$group <- ifelse(ps_cog_select$genome.y %in% more2, "WT", "OIL")
ps_cog_select$plant <- ifelse(ps_cog_select$Orig_name.y == "WT (CP88-1762)", "WT", "OIL")
f2 <- ddply(ps_cog_select, .(V3.y, V5, T1,  plant, genome.y, Sample.y), summarize, SUM = sum(MEAN)) 
f3 <- ddply(f2, .(T1, V5,  V3.y, plant, Sample.y), summarize, SUM2 = sum(SUM)) 
f4 <- ddply(f3, .(T1, V5, plant), summarize, MEAN = mean(SUM2), SE = sd(SUM2)/sqrt(length(SUM2)))
f4b <- subset(f4, V5 != "")
#limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)

t1_colors <- c(
  
  "#E69F00", # orange
  "#56B4E9", # sky blue
  "#009E73", # bluish green
  "#F0E442", # yellow
  "#0072B2", # blue
  "#D55E00", # vermillion
  "#CC79A7", # reddish purple
  "#999999", # gray
  "#000000", # black
  "#8B4513", # saddle brown (extra)
  "#800080"  # purple (extra)
  
)
list_t1 <- unique(f4b$T1)
taxon_col <- setNames(t1_colors, list_t1)

#p = ggplot(f4b, aes(x=V5, y=MEAN, fill = T1))      
#p = p+geom_bar(stat="identity", position="stack") +theme_bw()+theme(axis.text.x = element_text(angle = 90))
#p+ theme(axis.title.x = element_text(size=15),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15, vjust=0.5, hjust=1),
        axis.title = element_text(size = 15),
        strip.text.x = element_text(size=15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))+ylab('Abundance')+ xlab('Functional Genes')+facet_wrap(~plant)


df_mirror <- f4b %>%
  mutate(
    MEAN_signed = if_else(plant == "WT", -MEAN, MEAN),
    SE_signed   = if_else(plant == "WT", -SE,   SE)
  ) %>%
  # (optional) order genes by total abundance (OIL) so the most abundant appear at top
  group_by(V5) %>%
  mutate(total_OIL = sum(MEAN[plant == "OIL"], na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    V5 = fct_reorder(V5, total_OIL, .desc = TRUE),
    # (optional) control phylum legend order by overall abundance
    T1 = fct_reorder(T1, MEAN, .fun = sum, .desc = TRUE)
  )


p <- ggplot(df_mirror, aes(x = V5, y = MEAN_signed, fill = T1)) +
  geom_col(width = 0.85, position = "stack", color = NA) +
  # Zero line for visual reference
  geom_hline(yintercept = 0, color = "grey25") +
  coord_flip() + scale_fill_manual(values=taxon_col)+
  # Show absolute values on axis while keeping the mirror effect
  scale_y_continuous(labels = function(x) number(abs(x), accuracy = 0.0001),
                     expand = expansion(mult = c(0.05, 0.05))) +
  labs(x = "Functional Genes", y = "Abundance (mean)",
       fill = "Phylum",
       title = func_chosen,
       subtitle = "OIL (right, positive) vs WT (left, negative)") +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 9),
    legend.position = "right"
  )

p



###################################################################################################
head(f3)

df_rel <- f4 %>% select(V3.y, genome.y, Orig_name.y, T1.y, MEAN) %>%
  group_by(Orig_name.y) %>%
  mutate(
    MEAN_REL = MEAN / sum(MEAN)
  ) %>%
  ungroup()

# This is all functions unfiltered by what is unique
p = ggplot(f4, aes(x=Orig_name.y, y=MEAN, fill=T1.y))        
p = p+geom_bar(stat="identity", position="stack") +theme_bw()+theme(axis.text.x = element_text(angle = 90))
p+facet_grid(col=vars(V3.y), scales="free_x")+
  theme(axis.title.x = element_text(size=15),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.title = element_text(size = 15),
        strip.text.x = element_text(size=15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15), strip.text = element_text(angle = 90) )+ylab('Abundance (norm)')+ xlab('Functional Genes')+
        scale_fill_manual(values=taxon_col)


# Abundance of interesting MAGS
mag_descr <- unique(merged_abund %>% select(Sample.y, genome.y, MEAN, Orig_name.y, T1.y))
head(mag_ann)
mag_descr2 <- merge(mag_descr, mag_ann, by.x = "genome.y", by.y = "bin_oid")
mag_descr2$group <- ifelse(mag_descr2$genome.y %in% more2, "WT", "OIL")
mag_descr2$plant <- ifelse(mag_descr2$Orig_name.y == "WT (CP88-1762)", "WT", "OIL")
write.csv(mag_descr2, file="mag_info.txt", quote=FALSE)

####this is after you get the list of unique cogs
targets_cog <- merged_abund %>% select(Sample.y, genome.y, OTU, MEAN, Orig_name.y, T1.y, COG)
#This is the averaged abundance of genoem overlaid on each function and then used to make the functional plots
cog1 <- read.csv(file = "cog-20.def.tab.txt", header=FALSE, sep="\t")
cog2 <- read.csv(file = "fun-20.tab.txt", header=FALSE, sep="\t")
cog_merge <- merge(cog1, cog2, by.x = "V2", by.y = "V1")
ps_cog <- merge(targets_cog, cog_merge, by.x = "COG", by.y = "V1")
oil <- less2
wt <- more2

#mirror genome chart
annot2
annot3 <- subset(annot2, genome %in% c(more2, less2))
annot4 <- subset(annot3, COG != "")
annot4$ORF <- rownames(annot4)
df <- ddply(annot4, .(genome, COG), summarize, count = length(ORF))  #count is the number of ORFs in each COG
df$group <- ifelse(df$genome %in% more2, "WT", "OIL") 
df2 <- ddply(df, .(group, COG), summarize, MEAN = mean(count), SE=sd(count)/sqrt(length(count)))
df3 <- merge(df2, cog_merge, by.x = "COG", by.y = "V1")


df <- df3 %>% mutate(MEAN_mirrored = ifelse(group == "WT", -MEAN, MEAN)) #makes it mirrored for teh graphic
ggplot(df, aes(x = COG, y = MEAN_mirrored, fill = group)) +
  geom_bar(stat = "identity", position = "identity") +
  coord_flip() +
  scale_y_continuous(labels = abs) +
  labs(title = "Mirrored Bar Plot of COG Enrichment",
       y = "Mean Abundance",
       x = "COG") + facet_wrap(~V3.y, scale="free_y")+
  theme_minimal() +
  scale_fill_manual(values = c("OIL" = "red", "WT" = "steelblue"))



df2 <- df %>% #checks to see if the COG have the same number of ORFs, 25%
  group_by(`COG`) %>%
  mutate(
    flagged = if (n() == 2 &&
                  all(abs(MEAN_mirrored) == abs(MEAN_mirrored[1])) &&
                  any(MEAN_mirrored < 0) &&
                  any(MEAN_mirrored > 0)) TRUE else FALSE
  ) %>%
  ungroup()

df3 <- df2 %>% #COG is only in one group, 1036, 18%
  group_by(`COG`) %>%
  mutate(flag_missing_other = if (n() == 1 && MEAN_mirrored != 0) TRUE else FALSE) %>%
  ungroup()

df4 <- df3[df3$flag_missing_other,]
ggplot(df4, aes(x = COG, y = MEAN_mirrored, fill = group)) +
  geom_bar(stat = "identity", position = "identity") +
  coord_flip() +
  scale_y_continuous(labels = abs) +
  labs(y = "Mean Abundance",
       x = "COG") + facet_wrap(~V3.y, scale="free_y")+
  theme_minimal() +
  scale_fill_manual(values = c("OIL" = "red", "WT" = "steelblue"))+
  theme(axis.title.x = element_text(size=15),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 15),
        axis.title = element_text(size = 15),
        strip.text.x = element_text(size=15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))
  


df_summary <- df4 %>%
  group_by(V3.y, group) %>%
  summarise(COG_count = n(), .groups = "drop")
ggplot(df_summary, aes(x = V3.y, y = COG_count, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(y = "Mean Abundance",
       x = "COG") +
  theme_minimal() +
  scale_fill_manual(values = c("OIL" = "red", "WT" = "steelblue"))+
  theme(axis.title.x = element_text(size=15),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 15),
        axis.title = element_text(size = 15),
        strip.text.x = element_text(size=15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))+
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5))

list_of_unique_cogs <- unique(df4$COG)



f_cog <- subset(ps_cog, COG %in% list_of_unique_cogs)
f2 <- ddply(f_cog, .(V3.y, genome.y, Orig_name.y, T1.y, Sample.y), summarize, SUM = sum(MEAN)) 
f2b <- subset(f2, genome.y != "3300056388_325")
f3 <- ddply(f2b, .(V3.y, Orig_name.y, genome.y, T1.y), summarize, MEAN = mean(SUM) , SE = sd(SUM)/sqrt(length(SUM)))
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)

t1_colors <- c(
  
  "#E69F00", # orange
  "#56B4E9", # sky blue
  "#009E73", # bluish green
  "#F0E442", # yellow
  "#0072B2", # blue
  "#D55E00", # vermillion
  "#CC79A7", # reddish purple
  "#999999", # gray
  "#000000", # black
  "#8B4513", # saddle brown (extra)
  "#800080"  # purple (extra)
  
)
list_t1 <- unique(f3$T1.y)
taxon_col <- setNames(t1_colors, list_t1)


p = ggplot(f3, aes(x=Orig_name.y, y=MEAN, fill=T1.y))        
p = p+geom_bar(stat="identity", position="stack") +theme_bw()+theme(axis.text.x = element_text(angle = 90))
p+facet_grid(col=vars(V3.y), scales="free_x")+
  theme(axis.title.x = element_text(size=15),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.title = element_text(size = 15),
        strip.text.x = element_text(size=15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15), strip.text = element_text(angle = 90) )+ylab('Abundance')+ xlab('Functional Genes')+
        scale_fill_manual(values=taxon_col)

f_cog <- subset(ps_cog, COG %in% list_of_unique_cogs)
f_cog$group <- ifelse(f_cog$genome %in% more2, "WT", "OIL")
f_cog$plant <- ifelse(f_cog$Orig_name.y == "WT (CP88-1762)", "WT", "OIL")
f_cogb <- subset(f_cog, group == "WT")
f2 <- ddply(f_cogb, .(V3.y, COG, group, plant, genome.y, Sample.y), summarize, SUM = sum(MEAN)) 
f3 <- ddply(f2, .(V3.y, group, plant, Sample.y), summarize, SUM2 = sum(SUM)) 
f3 <- ddply(f2, .(V3.y, group, plant, Sample.y), summarize, SUM2 = sum(SUM)) 
f4 <- ddply(f3, .(V3.y, group, plant), summarize, MEAN = mean(SUM2), SE = sd(SUM2)/sqrt(length(SUM2)))

f4_wt <- f4

p = ggplot(f5, aes(x=V3.y, y=MEAN, color = plant, shape=plant)  )      
p = p+geom_point(size=5) +theme_bw()+theme(axis.text.x = element_text(angle = 90))+geom_errorbar(limits, width=0)
p+ theme(axis.title.x = element_text(size=15),
         axis.text.y = element_text(size = 15),
         axis.text.x = element_text(size = 15, vjust=0.5, hjust=1),
         axis.title = element_text(size = 15),
         strip.text.x = element_text(size=15),
         legend.text = element_text(size = 15),
         legend.title = element_text(size = 15))+ylab('Abundance')+ xlab('Functional Genes')+facet_wrap(~group)+scale_color_manual(values=c("red", "black"))


f_cog <- subset(ps_cog, COG %in% list_of_unique_cogs)
f_cog$group <- ifelse(f_cog$genome %in% more2, "WT", "OIL")
f_cog$plant <- ifelse(f_cog$Orig_name.y == "WT (CP88-1762)", "WT", "OIL")
f_cogb <- subset(f_cog, group == "OIL")
f2 <- ddply(f_cogb, .(V3.y, COG, group, plant, genome.y, Sample.y), summarize, SUM = sum(MEAN)) 
f3 <- ddply(f2, .(V3.y, group, plant, Sample.y), summarize, SUM2 = sum(SUM)) 
f4 <- ddply(f3, .(V3.y, group, plant), summarize, MEAN = mean(SUM2), SE = sd(SUM2)/sqrt(length(SUM2)))

f4_oil <- f4

f5 <- rbind(f4_wt, f4_oil)
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)

p = ggplot(f5, aes(x=V3.y, y=MEAN, color = plant, shape=plant)  )      
p = p+geom_point(size=5) +theme_bw()+theme(axis.text.x = element_text(angle = 90))+geom_errorbar(limits, width=0)
p+ theme(axis.title.x = element_text(size=15),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15, vjust=0.5, hjust=1),
        axis.title = element_text(size = 15),
        strip.text.x = element_text(size=15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))+ylab('Abundance')+ xlab('Functional Genes')+facet_wrap(~group)+scale_color_manual(values=c("red", "black"))



t1_colors <- c(
  
  "#E69F00", # orange
  "#56B4E9", # sky blue
  "#009E73", # bluish green
  "#F0E442", # yellow
  "#0072B2", # blue
  "#D55E00", # vermillion
  "#CC79A7", # reddish purple
  "#999999", # gray
  "#000000", # black
  "#8B4513", # saddle brown (extra)
  "#800080"  # purple (extra)
  
)
list_t1 <- unique(f3$T1.y)
taxon_col <- setNames(t1_colors, list_t1)


head(ps_cog)
f_cog <- subset(ps_cog, COG %in% list_of_unique_cogs)
f_cog$group <- ifelse(f_cog$genome %in% more2, "WT", "OIL")
f_cog$plant <- ifelse(f_cog$Orig_name.y == "WT (CP88-1762)", "WT", "OIL")
f2 <- ddply(f_cog, .(V3.y, COG, group, plant, genome.y, T1.y, Sample.y), summarize, SUM = sum(MEAN)) 
f3 <- ddply(f2, .(COG, group, plant, genome.y, T1.y), summarize, MEAN = mean(SUM)) 
f4 <- merge(f3, cog_merge, by.x = "COG", by.y = "V1")
f4b <- subset(f4, genome.y == "3300056388_133")
list_t1 <- unique(f4$T1.y)
taxon_col <- setNames(t1_colors, list_t1)
l <- unique(f4$V3.y)
i = 24
f4b <- subset(f4, V3.y == l[i])
head(f4b)

unique(subset(f4b, group == "WT")$V3.x)
unique(subset(f4b, group == "OIL")$V3.x)
# Calculate mean MEAN per COG
mean_vals <- tapply(f4b$MEAN, f4b$V3.x, mean)

# Order COG levels by mean
ordered_levels <- names(sort(mean_vals))

# Set factor levels
f4b$V3.x <- factor(f4b$V3.x, levels = ordered_levels)

p = ggplot(f4b, aes(x=V3.x, y=MEAN, color = T1.y, shape=plant))   
p = p+geom_point(size=5) +theme_bw()+theme(axis.text.x = element_text(angle = 90))+facet_wrap(~group, scale="free_x")
p+ theme(axis.title.x = element_text(size=15),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15, vjust=0.5, hjust=1),
        axis.title = element_text(size = 15),
        strip.text.x = element_text(size=15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))+ylab('Abundance')+ xlab(l[i])+scale_color_manual(values=taxon_col)
ggsave(filename='foo.png', width=30, height=20)
####ps_cog#####################################################

merged_abund <- merge(targets, f, by = "name2") #so now there is avg mean assigned to each abundance, MEAN
targets_cog <-merged_abund %>% select(Sample.y, genome.y, OTU, MEAN, Orig_name.y, T1.y, COG)
#This is the averaged abundance of genoem overlaid on each function and then used to make the functional plots
targets_cog2 <- subset(targets_cog, genome.y %in% less2)
cog1 <- read.csv(file = "cog-20.def.tab.txt", header=FALSE, sep="\t")
cog2 <- read.csv(file = "fun-20.tab.txt", header=FALSE, sep="\t")
cog_merge <- merge(cog1, cog2, by.x = "V2", by.y = "V1")
ps_cog <- merge(targets_cog2, cog_merge, by.x = "COG", by.y = "V1")
ps_cog2 <- subset(ps_cog, V3.y != "")
f <- ddply(ps_cog2, .(COG, V3.y, genome.y, T1.y, Sample.y,Orig_name.y), summarise, SUM=sum(MEAN)) #sums COGs over ORFs
f2 <- ddply(f, .(V3.y, genome.y, Orig_name.y, T1.y, Sample.y), summarize, SUM2 = sum(SUM)) 
f2b <- subset(f2, genome.y != "3300056388_325")

t1_colors <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
  "#D55E00", "#CC79A7", "#999999", "#000000", "#FFFFFF"
)


list_t1 <- unique(f3$T1.y)
taxon_col <- setNames(t1_colors, list_t1)

f3 <- ddply(f2b, .(V3.y, Orig_name.y, genome.y, T1.y), summarize, MEAN = mean(SUM2) , SE = sd(SUM2)/sqrt(length(SUM2)))
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
p = ggplot(f3, aes(x=Orig_name.y, y=MEAN, color=T1.y))        
p = p+geom_point(size=3) +theme_bw()+theme(axis.text.x = element_text(angle = 90))+geom_errorbar(limits, width=0)
p+facet_grid(col=vars(V3.y), scales="free_x")+
  theme(axis.title.x = element_text(size=15),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        axis.title = element_text(size = 15),
        strip.text.x = element_text(size=15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15), strip.text = element_text(angle = 90) )+ylab('Abundance')+ xlab('Functional Genes')+
        scale_colour_manual(values=taxon_col)

oil <- less2
wt <- more2
targets <- subset(ps_WT_all, genome %in% c(oil, wt))
targets$group <- ifelse(targets$genome %in% wt, "WT", "oilcane")
f <- ddply(targets, .(genome,  T1, Sample, Orig_name, group), summarise, MEAN=mean(norm)) #Avg coverage of a genome
f_wt <- subset(f, group == "WT")
df <- f_wt
f_oil <- subset(f, group == "oilcane")
df2 <- f_oil
df <- f
####Data for Sankey
head(targets_cog)
targets_cog$group <- ifelse(targets_cog$genome %in% more2, "WT", "OIL")
targets_cog$plant <- ifelse(targets_cog$Orig_name.y == "WT (CP88-1762)", "WT", "OIL")
targets2_gen <- unique(targets_cog %>% select(genome.y, COG, group, OTU))
targets3 <- merge(targets2_gen, cog_merge, by.x = "COG", by.y = "V1")
write.csv(targets3, file="targets3.csv", quote=FALSE)

targets4 <- ddply(targets3, .(genome.y,  COG, group, V3.x, V3.y, V5), summarise, count_orfs=length(unique(OTU))) #Avg coverage of a genome

write.csv(targets4, file="targets4.csv", quote=FALSE)
targets5 <- unique(targets_cog %>% select(genome.y, Sample.y, MEAN, group, plant))
write.csv(targets4, file="targets4.csv", quote=FALSE)


targets6 <- targets5 %>%  group_by(plant, genome.y) %>% mutate(Group_MEAN = mean(MEAN, na.rm = TRUE)) %>% ungroup()
targets7 <- unique(targets6 %>% select(genome.y, group, Group_MEAN, plant))
targets8 <- merge(targets7, mag_ann, by.x = 'genome.y', by.y = "bin_oid")
write.csv(targets8, file="targets8.csv", quote=FALSE)



library(data.table)
dt <- as.data.table(df)
sample_sums_wt <- dt[, .(sum_MEAN_wt = sum(MEAN, na.rm = TRUE)), by = Sample]
wt <- as.data.frame(sample_sums_wt)
dt <- as.data.table(df2)
sample_sums_oil <- dt[, .(sum_MEAN_oil = sum(MEAN, na.rm = TRUE)), by = Sample]
oil <- as.data.frame(sample_sums_oil)

combo <- merge(wt, oil, by.x ="Sample", by.y = "Sample")
ggplot(combo, aes(x= sum_MEAN_wt, y = sum_MEAN_oil)) + geom_point()
meta1 <- read.csv(file="~/Box Sync/Oilcane/metadata_oilcane2.csv", header=TRUE, row.names = 1)
meta2 <- read.csv(file="~/Box Sync/Oilcane/jgi-oilcane.csv", header=TRUE)
meta3 <- read.csv(file="~/Box Sync/Oilcane/sample_key.txt", header=FALSE, sep="\t")
meta1$X
meta2$Seq.Project.Name
meta2$Seq.Project.Name <- gsub(" metagenome", "", meta2$Seq.Project.Name)
meta2$Seq.Project.Name
meta_merge1 <- merge(meta2, meta1, by.x = c("Seq.Project.Name"), by.y = c("sample_id"))
meta <- merge(meta_merge1, meta3, by.x=c("Seq.Project.Name"), by.y =  c("V2"))
head(meta)
combo2 <- merge(combo, meta, by.x = "Sample", by.y = "V3")
combo3 <- subset(combo2, sum_MEAN_wt < 0.0050)
combo_s <- subset(combo2, Orig_name == "WT (CP88-1762)")
ggplot(combo2, aes(x= sum_MEAN_wt, y = sum_MEAN_oil, color = Orig_name)) + geom_point(size=6)+
  geom_smooth(method = "lm", se = TRUE, color = "black") +xlab('Average abundance (norm) of WT')
model <- lm(sum_MEAN_oil ~ sum_MEAN_wt, data = combo2)
summary(model)
ggplot(combo3, aes(x= sum_MEAN_wt, y = sum_MEAN_oil, color=Orig_name)) + geom_point(size=6)+
  geom_smooth(method = "lm", se = TRUE, color = "black") +xlab('Average abundance (norm) of WT-enriched Taxa')+ylab('Average abundance (norm) of oilcane-enriched MAGs')
model <- lm(sum_MEAN_oil ~ sum_MEAN_wt, data = combo3)
summary(model)
ggplot(combo2, aes(x= sum_MEAN_wt, y = sum_MEAN_oil, color=Orig_name, shape=Sampling.time)) + geom_point(size=6)+
  geom_smooth(method = "lm", se=FALSE, linetype='dotted') +xlab('Average abundance (norm) of WT-enriched Taxa')+ylab('Average abundance (norm) of oilcane-enriched taxa')+theme_bw()+
  theme(axis.text = element_text(size = 14), axis.title = element_text(size=14))
model <- lm(sum_MEAN_oil ~ sum_MEAN_wt, data = combo2)
summary(model)

######special#####
f #this is the mean abundance of each genome per sample
f2 <- ddply(f, .(genome, Orig_name, T1, group), summarise, MEAN2=mean(MEAN), SE2=sd(MEAN)/sqrt(length(MEAN))) #Average of the observation of genomes 
#f2 is a good dataset and is the mean per cultivar of each genome.
#this will dicate the size of the dots.
genome_abund <- f2

ps_WT_subset <- subset(ps_WT_all, genome %in%   c(more2, less2))
ps_WT_sig <- ps_WT_subset
cog1 <- read.csv(file = "cog-20.def.tab.txt", header=FALSE, sep="\t")
cog2 <- read.csv(file = "fun-20.tab.txt", header=FALSE, sep="\t")
cog_merge <- merge(cog1, cog2, by.x = "V2", by.y = "V1")
ps_cog <- merge(ps_WT_sig, cog_merge, by.x = "COG", by.y = "V1")
ps_cog2 <- subset(ps_cog, V5 != "")
f <- ddply(ps_cog2, .(V3.y, genome, T1, Sample,Orig_name), summarise, SUM=sum(norm)) #sums COGs over ORFs
f2 <- ddply(f, .(V3.y, genome, T1, Orig_name), summarise, MEAN = mean(SUM))
f2$pa <- ifelse(f2$MEAN > 0, 1, 0)
f3 <- f2 %>% select(V3.y, genome, T1, Orig_name, pa)
genome_abund$uniq_id <- paste0(genome_abund$genome, "_", genome_abund$Orig_name)
f3$uniq_id <- paste0(f3$genome, "_", f3$Orig_name)
f4 <- merge(f3, genome_abund, by = "uniq_id")
ggplot(f4, aes(x = Orig_name.x, y = genome.x)) + geom_point(aes(size = MEAN2, fill = T1.y), shape = 21, color = "black") +
  scale_size(range = c(1, 10)) +
  theme_minimal() +
  labs(x = "Cultivar",y = "MAG", size = "MAG Coverage", fill = "Phyla") +
  theme(panel.grid.major = element_line(color = "grey90"))+facet_grid(group~V3.y, scales="free_y")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0))
# 
# as <- read.csv(file = "~/Box Sync/Papers/Oilcane/antismash.txt", header=FALSE, sep="\t")
# as2 <- as %>% select(V1, V3)
# as3 <- unique(as2)
# as4 <- as.data.frame(as3 %>% mutate(value=1) %>% distinct(V1,V3, .keep_all = TRUE) %>% pivot_wider(names_from= V1, values_from = value, values_fill=0))
# as5 <- reshape2::melt(as4)
# as5$genome <- as5$variable
# as6 <- merge(genome_abund, as5, by = "genome")
# ggplot(as6, aes(x = Orig_name, y = genome)) + geom_point(aes(size = MEAN2, fill = T1), shape = 21, color = "black") +
#   scale_size(range = c(1, 10)) +
#   theme_minimal() +
#   labs(x = "Cultivar",y = "MAG", size = "MAG Coverage", fill = "Phyla") +
#   theme(panel.grid.major = element_line(color = "grey90"))+facet_grid(group~V3, scales="free_y")+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0))


#Project Make functional diagram
ps_WT_subset <- subset(ps_WT_all, genome %in%   c(more2, less2))
ps_WT_sig <- ps_WT_subset
cog1 <- read.csv(file = "cog-20.def.tab.txt", header=FALSE, sep="\t")
cog2 <- read.csv(file = "fun-20.tab.txt", header=FALSE, sep="\t")
cog_merge <- merge(cog1, cog2, by.x = "V2", by.y = "V1")
ps_cog <- merge(ps_WT_sig, cog_merge, by.x = "COG", by.y = "V1")
ps_cog2 <- subset(ps_cog, V5 != "")
f <- ddply(ps_cog2, .(V3.y, COG, genome, T1, Sample,Orig_name), summarise, SUM=sum(norm)) #sums COGs over ORFs
f2 <- ddply(f, .(V3.y, COG, genome, T1, Orig_name), summarise, MEAN = mean(SUM))
f2$group <- ifelse(f2$genome %in% less2, "OIL", "WT")
f3 <- ddply(f2, .(V3.y, COG, genome, T1, group), summarise, MEAN2 = mean(MEAN))
f3$pa <- ifelse(f3$MEAN > 0, 1, 0)
f3$genome <- factor(f3$genome, levels = c("3300056388_1321", "3300056388_1831", "3300056388_398", "3300056388_1866", "3300056388_1875", "3300056388_1894", "3300056388_1898", 
                                          "3300056388_2006", "3300056388_2063", "3300056388_2120", "3300056388_229", 
                                          "3300056388_2568", "3300056388_2701", "3300056388_638"))
ggplot(f3, aes(x = COG, y = genome)) + geom_point(aes(fill = T1), shape = 21, size=1.5, color = "black") +
  facet_grid(~V3.y, scale = "free_y")+
  theme_minimal() +
  labs(x = "COG",y = "MAG", size = "MAG Coverage", fill = "Phyla") +
  theme(axis.text.x =element_blank(), strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0))


ggplot(f3, aes(x = V3.y, y = MEAN2, fill=T1)) + geom_bar(stat="identity", position="stack")+
  theme_minimal() + facet_wrap(~group)+
  labs(x = "COG",y = "MAG", size = "MAG Coverage", fill = "Phyla") +
  theme(axis.text.x =element_text(angle=90), strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0))


#f2 here has the per orginam
ggplot(f2, aes(x = COG, y = MEAN, )) + geom_point(aes(fill = T1), shape = 21, size=1.5, color = "black") +
  facet_grid(Orig_name~V3.y, scale = "free_y")+
  theme_minimal() + 
  labs(x = "COG",y = "MAG", size = "MAG Coverage", fill = "Phyla") +
  theme(axis.text.x =element_blank(), strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0))



f4 <- f3 %>%
  group_by(COG) %>%
  mutate(in_both_groups = all(c("WT", "OIL") %in% unique(group))) %>%
  ungroup()
f4$genome <- factor(f4$genome, levels = rev(c("3300056388_1321", "3300056388_1831", "3300056388_398", "3300056388_1866", "3300056388_1875", "3300056388_1894", "3300056388_1898", 
                                          "3300056388_2006", "3300056388_2063", "3300056388_2120", "3300056388_229", 
                                          "3300056388_2568", "3300056388_2701", "3300056388_638")))

# Map group to fill color: WT -> white, OIL -> black

ggplot(f4, aes(x = COG, y = genome)) + geom_point(aes(fill = in_both_groups), shape = 21, size=2, color = "black") +
  scale_fill_manual(values = c("TRUE" = "white", "FALSE" = "red"))+
  facet_grid(~V3.y, scale = "free_y")+
  theme_minimal() +
  labs(x = "COG",y = "MAG", size = "MAG Coverage", fill = "Phyla") +
  theme(axis.text.x =element_blank(), strip.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0))


head(f4)
#this subsets the in both groups to only the unique ones
f5 <- subset(f4, in_both_groups == FALSE)
f5 <- subset(f5, genome != "NA")
f5 <- merge(f5, cog_merge, by.x= "COG", by.y = "V1")
#f6 <- subset(f5, V3.y.y)
ggplot(f5, aes(x = COG, y = genome)) + geom_point(aes(fill = T1), shape = 21, size=4, color = "black") +
  theme_minimal() + facet_grid(~V3.y.x)+
  labs(x = "COG",y = "MAG", size = "MAG Coverage", fill = "Phyla") +
  theme(axis.text.x = element_text(angle = 90, size=12), strip.text.x = element_text(angle = 90, size=12, vjust = 0.5, hjust = 0))

list_of_unique_cogs <- unique(f5$COG)
#select specific panels 
f5b <- subset(f5, V3.y.x %in% c("Amino acid transport and metabolism", "Carbohydrate transport and metabolism", "Coenzyme transport and metabolism",
                                "Energy production and conversion", "Lipid transport and metabolism", "Nucleotide transport and metabolism", "Translation, ribosomal structure and biogenesis" ))

ggplot(f5b, aes(x = V3.x, y = genome)) + geom_point(aes(fill = T1), shape = 21, size=4, color = "black") +
  theme_minimal() + facet_grid(~V3.y.x, scales="free_x")+
  labs(x = "V3.x",y = "MAG", size = "MAG Coverage", fill = "Phyla") +
  theme(axis.text.x = element_text(angle = 90, size=10, vjust = 0.5, hjust = 1), strip.text.x = element_text(angle = 90, size=12, vjust = 1, hjust = 0))

f5c <- subset(f5b, T1 == " Thermoproteota")

h <- subset(f5, !genome %in% c("3300056388_1831", "3300056388_1321", "3300056388_398"))
l <- unique(f5$V3.y.x)
i = 2
subset_df <- subset(f5, V3.y.x == l[i])
all_genomes <- rev(c("3300056388_1321", "3300056388_1831", "3300056388_398", "3300056388_1866", "3300056388_1875", "3300056388_1894", "3300056388_1898", 
  "3300056388_2006", "3300056388_2063", "3300056388_2120", "3300056388_229", 
  "3300056388_2568", "3300056388_2701", "3300056388_638"))




# Create a full grid of all genomes and Sample_IDs
full_grid <- expand.grid(genome = all_genomes, V3.x = subset_df$V3.x)

# Join with your subset and fill missing MEAN2 with 0
filled_df <- full_grid %>%
  left_join(subset_df, by = c("genome", "V3.x")) %>%
  mutate(MEAN2 = ifelse(is.na(MEAN2), 0, MEAN2))

f2 <- ddply(filled_df, .(genome, V5))
ggplot(filled_df, aes(x = V5, y = genome)) + geom_point(aes(fill = T1, size=MEAN2), shape = 21, color = "black") +
  theme_minimal() +
  labs(x = l[i],y = "MAG", size = "MAG Coverage", fill = "Phyla") +
  theme(axis.text.x = element_text(angle = 90, size = 12, vjust=0.5, hjust=1))


# Making a heatmap of genome abundance
f <- ddply(targets, .(genome,  T1, Sample, Orig_name, Sampling.time, group), summarise, MEAN=mean(norm)) #Avg coverage of a genome
f2 <- ddply(f, .(genome, Orig_name, T1, Sampling.time, group), summarise, MEAN2=mean(MEAN), SE2=sd(MEAN)/sqrt(length(MEAN))) #Average of the observation of genomes 
f2a <- subset(f2, Sampling.time == "T1")
f2b <- subset(f2, Sampling.time == "T2")


desired_order <- c("3300056388_1321", "3300056388_1831", "3300056388_398", "3300056388_1866", 
                   "3300056388_1875", "3300056388_1894", "3300056388_1898", "3300056388_2006", 
                   "3300056388_2063", "3300056388_2120", "3300056388_229", "3300056388_2568", 
                   "3300056388_2701", "3300056388_638")

heatmap_data <- f2a %>%
  select(genome, Orig_name, MEAN2) %>%
  pivot_wider(names_from = Orig_name, values_from = MEAN2)

# Reorder rows based on desired order
heatmap_data <- heatmap_data %>%
  filter(genome %in% desired_order) %>%
  mutate(genome = factor(genome, levels = desired_order)) %>%
  arrange(genome)

heatmap_matrix <- as.data.frame(heatmap_data)
rownames(heatmap_matrix) <- heatmap_matrix$genome
heatmap_matrix$genome <- NULL

heatmap_matrix <- as.matrix(heatmap_matrix)

# Step 4: Create heatmap
p <- pheatmap(heatmap_matrix,  color = brewer.pal(10, 'Spectral'),
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         scale = "none", clustering_distance_rows = 'euclidean', clustering_method = 'average',
         main = "Mean Values per Genome")


# Include the sampling.time
df <- f2 %>%
  mutate(Sample_ID = paste(Orig_name, Sampling.time, sep = "_"))

# Reshape the data to wide format
desired_order <- c("3300056388_1321", "3300056388_1831", "3300056388_398", "3300056388_1866", 
                   "3300056388_1875", "3300056388_1894", "3300056388_1898", "3300056388_2006", 
                   "3300056388_2063", "3300056388_2120", "3300056388_229", "3300056388_2568", 
                   "3300056388_2701", "3300056388_638")


# Reshape the data to wide format
heatmap_data <- df %>%
  select(genome, Sample_ID, MEAN2) %>%
  pivot_wider(names_from = Sample_ID, values_from = MEAN2)


# Reorder rows based on desired order
heatmap_data <- heatmap_data %>%
  filter(genome %in% desired_order) %>%
  mutate(genome = factor(genome, levels = desired_order)) %>%
  arrange(genome)

rownames(heatmap_data) <- heatmap_data$genome
heatmap_data$genome <- NULL
heatmap_matrix <- as.matrix(heatmap_data)
heatmap_matrix2 <- log(heatmap_matrix)

# Plot the heatmap
pheatmap(heatmap_matrix2,color = brewer.pal(10, 'Spectral'),
         cluster_rows = FALSE,
         cluster_cols = FALSE, scale="none", 
         fontsize_row = 10,
         fontsize_col = 10)

#f5 here are the functions that are only in one or the other
#oil only
f5b <- subset(f5, group == "OIL")
ggplot(f5b, aes(x = genome, y = MEAN2, fill=)) + geom_point() +
  theme_minimal() + facet_wrap(~V3.y.y)
  labs(x = l[i],y = "MAG", size = "MAG Coverage", fill = "Phyla") +
  theme(axis.text.x = element_text(angle = 90, size = 12, vjust=0.5, hjust=1))

