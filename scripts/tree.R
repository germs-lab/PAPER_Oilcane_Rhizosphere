## ============================================================
## Libraries
## ============================================================


library(ggtree)
library(ggtreeExtra)
library(treeio)
library(ggplot2)


## ============================================================
## Loading in Data
## ============================================================
phy_filter4
ps <- psmelt(phy_filter4)
bin_info <- read.csv('exported_img_data.tsv', sep="\t")
ps_merged <- merge(ps, bin_info, by.x = "genome", by.y = "bin_oid")
x <- scan("./list-of-all-hkgs-in-mags.uniq.list2", what="", sep="\n")
hkg_phy <- prune_taxa(x, phy_filter2)
phy_melt_hkg = psmelt(hkg_phy)
f_hkg <- ddply(phy_melt_hkg, .(Sample), summarise, SUM=sum(Abundance))
ps_merged2 <- merge(ps_merged, f_hkg, by.x = "Sample", by.y = "Sample")
ps_merged2$norm <- ps_merged2$Abundance/ps_merged2$SUM
f <- ddply(ps_merged2, .(Sample, genome), summarise, MEAN=mean(Abundance), SE = sd(Abundance)/sqrt(length(Abundance)))
filtered_list <- unique(subset(f, MEAN > 5)$genome) 
ps_merged3 <- subset(ps_merged2, genome %in% filtered_list)
ps_WT_all <- ps_merged3 %>% mutate(is_WT = ifelse(Orig_name == "WT (CP88-1762)", "WT (CP88-1762)", "Oilcane"))
f <- ddply(ps_WT_all, .(genome, Sampling.time, T1, Sample, Orig_name), summarise, MEAN=mean(norm)) #Avg coverage of a genome
library(dplyr)
library(tidyr)
df <- f %>% select(genome, Sample, MEAN)

otu_table <- as.data.frame(df %>%
  pivot_wider(names_from = Sample, values_from = MEAN))
rownames(otu_table) <- otu_table$genome
otu_table$genome <- NULL
OTU <- otu_table(as.matrix(otu_table), taxa_are_rows = TRUE)
df <- f %>% select(genome, T1, )
mag_ann <- read.csv(file="exported_img_data.tsv", sep="\t")
rownames(mag_ann) <- mag_ann$bin_oid
mag_ann$bin_oid <- NULL
TAX <- tax_table(as.matrix(mag_ann))
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
META <- sample_data(meta)

ps <- phyloseq(TAX, META, OTU)
saveRDS(ps, 'ps_57.RDS')
tree <- read.tree('~/Box Sync/Oilcane/MAG-based/mags-25.tre')
tree$tip.label <- sub("\\.fna$", "", tree$tip.label)

ps <- merge_phyloseq(ps, tree)


ps_melt <- psmelt(ps)
melt_simple <- psmelt(ps)  %>% 
  select(OTU, val=Abundance)

ps_WT_all <- ps_melt %>% mutate(is_WT = ifelse(Orig_name == "WT (CP88-1762)", "WT (CP88-1762)", "Oilcane"))


## ============================================================
## Building Tree - WT vs Oilcane Greenhouse and Field - Figure 3
## ============================================================

subset_select_phyloseq <- function(data, group_value, time_value, group_column = "is_WT", time_column = "Sampling.time") {
  # Subset based on input values
  subset_data <- subset(data, 
                        get(group_column) == group_value & 
                          get(time_column) == time_value)
  
  # Select OTU and Abundance, rename Abundance to val
  selected_data <- subset_data %>% select(OTU, T1, val = Abundance)
  
  return(selected_data)
}

wt_t1 <-  subset_select_phyloseq(ps_WT_all, group_value = "WT (CP88-1762)", time_value = "T1")
wt_t2 <-  subset_select_phyloseq(ps_WT_all, group_value = "WT (CP88-1762)", time_value = "T2")
oil_t1 <-  subset_select_phyloseq(ps_WT_all, group_value = "Oilcane", time_value = "T1")
oil_t2 <-  subset_select_phyloseq(ps_WT_all, group_value = "Oilcane", time_value = "T2")


all_vals <- c(wt_t1$val, wt_t2$val, oil_t1$val, oil_t2$val)
global_min <- min(all_vals)
global_max <- max(all_vals)


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

list_t1 <- unique(ps_melt$T1)
taxon_col <- setNames(t1_colors, list_t1)

p_tree = ggtree(ps)
p = ggtree(ps)+ geom_tippoint(aes(color = T1), size = 3)+  theme(legend.position = "left")+
  scale_color_manual(values =taxon_col)



otu_order <- unique(p$data %>%
  filter(isTip) %>%
  arrange(y) %>%
  pull(label))

# 3. Set OTU factor levels in boxplot data to match tree order
wt_t1$OTU <- factor(wt_t1$OTU, levels = otu_order)
wt_t2$OTU <- factor(wt_t2$OTU, levels = otu_order)
oil_t1$OTU <- factor(oil_t1$OTU, levels = otu_order)
oil_t2$OTU <- factor(oil_t2$OTU, levels = otu_order)

p1 = ggplot(wt_t1, aes(x = val, y = OTU, fill=T1)) +
  geom_boxplot() +
  coord_cartesian(xlim = c(global_min, global_max)) +  # Controls x-axis range
  theme_minimal() +
  scale_fill_manual(values =taxon_col)+
  theme(legend.position = "none")+
  labs(title = "WT Greenhouse") +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )
p2 = ggplot(wt_t2, aes(x = val, y = OTU, fill=T1)) +
  geom_boxplot() +
  coord_cartesian(xlim = c(global_min, global_max)) +  # Controls x-axis range
  theme_minimal() +
  scale_fill_manual(values =taxon_col)+
  theme(legend.position = "none")+
  labs(title = "WT Field") +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )
p3 = ggplot(oil_t1, aes(x = val, y = OTU, , fill=T1)) +
  geom_boxplot() +
  coord_cartesian(xlim = c(global_min, global_max)) +  # Controls x-axis range
  theme_minimal() +
  scale_fill_manual(values =taxon_col)+
  theme(legend.position = "none")+
  labs(title = "Oilcane Greenhouse") +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )
p4 = ggplot(oil_t2, aes(x = val, y = OTU, fill=T1)) +
  geom_boxplot() +
  coord_cartesian(xlim = c(global_min, global_max)) +  # Controls x-axis range
  theme_minimal() +
  scale_fill_manual(values =taxon_col)+
  theme(legend.position = "none")+
  labs(title = "Oilcane Field") +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )
library(patchwork)
p + p1+p2+p3+p4+plot_layout(nrow=1)


## ============================================================
## Building Tree - By Accession
## ============================================================

subset_select_phyloseq <- function(data, group_value, time_value, group_column = "Orig_name") {
  # Subset based on input values
  subset_data <- subset(data, 
                        get(group_column) == group_value)
  
  # Select OTU and Abundance, rename Abundance to val
  selected_data <- subset_data %>% select(OTU, T1, val = Abundance)
  
  return(selected_data)
}

subset_not_select_phyloseq <- function(data, group_value, time_value, group_column = "Orig_name") {
  # Subset based on input values
  subset_data <- subset(data, 
                        get(group_column) != group_value)
  
  # Select OTU and Abundance, rename Abundance to val
  selected_data <- subset_data %>% select(OTU, T1, val = Abundance)
  
  return(selected_data)
}


foo1 <-  subset_select_phyloseq(ps_WT_all, group_value = "WT (CP88-1762)")
foo2 <-  subset_not_select_phyloseq(ps_WT_all, group_value = "WT (CP88-1762)")
p1 = ggplot(foo1, aes(x = val, y = OTU, fill=T1)) +
  geom_boxplot() +
  coord_cartesian(xlim = c(global_min, global_max)) +  # Controls x-axis range
  theme_minimal() +
  scale_fill_manual(values =taxon_col)+
  theme(legend.position = "none")+
  labs(title = "WT Greenhouse") +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )
p2 = ggplot(foo2, aes(x = val, y = OTU, fill=T1)) +
  geom_boxplot() +
  coord_cartesian(xlim = c(global_min, global_max)) +  # Controls x-axis range
  theme_minimal() +
  scale_fill_manual(values =taxon_col)+
  theme(legend.position = "none")+
  labs(title = "WT Field") +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

p + p1+p2+plot_layout(nrow=1)


foo_wt <-  subset_select_phyloseq(ps_WT_all, group_value = "WT (CP88-1762)")
foo1 <-  subset_select_phyloseq(ps_WT_all, group_value = "1565")
foo2 <-  subset_select_phyloseq(ps_WT_all, group_value = "1566")
foo3 <-  subset_select_phyloseq(ps_WT_all, group_value = "1569")
foo4 <-  subset_select_phyloseq(ps_WT_all, group_value = "1580")

foo_wt$OTU <- factor(foo_wt$OTU, levels = otu_order)
foo1$OTU <- factor(foo1$OTU, levels = otu_order)
foo2$OTU <- factor(foo2$OTU, levels = otu_order)
foo3$OTU <- factor(foo3$OTU, levels = otu_order)
foo4$OTU <- factor(foo4$OTU, levels = otu_order)

p1 = ggplot(foo_wt, aes(x = val, y = OTU, fill=T1)) +
  geom_boxplot() +
  coord_cartesian(xlim = c(global_min, global_max)) +  # Controls x-axis range
  theme_minimal() +
  scale_fill_manual(values =taxon_col)+
  theme(legend.position = "none")+
  labs(title = "WT Greenhouse") +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )
p2 = ggplot(foo1, aes(x = val, y = OTU, fill=T1)) +
  geom_boxplot() +
  coord_cartesian(xlim = c(global_min, global_max)) +  # Controls x-axis range
  theme_minimal() +
  scale_fill_manual(values =taxon_col)+
  theme(legend.position = "none")+
  labs(title = "1565") +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )
p3 = ggplot(foo2, aes(x = val, y = OTU, fill=T1)) +
  geom_boxplot() +
  coord_cartesian(xlim = c(global_min, global_max)) +  # Controls x-axis range
  theme_minimal() +
  scale_fill_manual(values =taxon_col)+
  theme(legend.position = "none")+
  labs(title = "1566") +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )
p4 = ggplot(foo3, aes(x = val, y = OTU, fill=T1)) +
  geom_boxplot() +
  coord_cartesian(xlim = c(global_min, global_max)) +  # Controls x-axis range
  theme_minimal() +
  scale_fill_manual(values =taxon_col)+
  theme(legend.position = "none")+
  labs(title = "1569") +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )
p5 = ggplot(foo4, aes(x = val, y = OTU, fill=T1)) +
  geom_boxplot() +
  coord_cartesian(xlim = c(global_min, global_max)) +  # Controls x-axis range
  theme_minimal() +
  scale_fill_manual(values =taxon_col)+
  theme(legend.position = "none")+
  labs(title = "1580") +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

p + p1+p2+p3+p4+p5+plot_layout(nrow=1)
