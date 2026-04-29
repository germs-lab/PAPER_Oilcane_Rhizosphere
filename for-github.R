## ============================================================
## Libraries
## ============================================================

library(plyr)
library(ggplot2)
library(phyloseq)
library(ggthemes)
library(vegan)

## ============================================================
## Load and build phyloseq object
## ============================================================
annot <- read.csv('./data/all-mags-derep-renamed.gff.annot2', sep='\t', header=FALSE)
cog_list <- read.csv2(file="./data/cog-20.def.tab.txt", sep="\t", header=FALSE)
annot2 <- merge(annot, cog_list, by.x = "V3", by.y = "V1", all.x = TRUE)
rownames(annot2) <- annot2$V1
annot2$V1 <- NULL
annot2 <- annot2[, colSums(is.na(annot2)) < nrow(annot2)]
colnames(annot2) <- c("COG", "genome", "COG_letters", "description", "gene", "descr2", 
                      "cog_num", "alphnum")
annot3 <- tax_table(as.matrix(annot2))
foo <- read.csv('./data/summary-count.txt')
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
meta1 <- read.csv(file="./data/metadata_oilcane2.csv", header=TRUE, row.names = 1)
meta2 <- read.csv(file="./data/jgi-oilcane.csv", header=TRUE)
meta3 <- read.csv(file="./data/sample_key.txt", header=FALSE, sep="\t")
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
phy_filter = prune_samples(sample_sums(phy)>=10000, phy)
phy_filter2 = phy
phy_filter2
phy_filter3 <- subset_samples(phy_filter2, Original.material..Soil..Roots..Stalks..Leaves. == "Root-associated soils")
phy_filter5 <- subset_samples(phy_filter3, Orig_name != "17T")
phy_filter3 <- phy_filter5

ps <- psmelt(phy_filter3)
bin_info <- read.csv('./data/exported_img_data.tsv', sep="\t")
ps_merged <- merge(ps, bin_info, by.x = "genome", by.y = "bin_oid")
save(ps_merged, file='./results/rds/ps_merged.save.rds')

x <- scan("./data/list-of-all-hkgs-in-mags.uniq.list2", what="", sep="\n")
hkg_phy <- prune_taxa(x, phy_filter2)
phy_melt_hkg = psmelt(hkg_phy)
f_hkg <- ddply(phy_melt_hkg, .(Sample), summarise, SUM=sum(Abundance))

ps_merged2 <- merge(ps_merged, f_hkg, by.x = "Sample", by.y = "Sample")
save(ps_merged2, file='./results/rds/ps_merged2.save.rds')


## ============================================================
## Genome filtering (RAW abundance)
## ============================================================
f <- ddply(ps_merged2, .(Sample, genome), summarise, MEAN=mean(Abundance), SE = sd(Abundance)/sqrt(length(Abundance)))
filtered_list <- unique(subset(f, MEAN > 5)$genome)  #this reduces the 377 to 179 genomes based on average coverage
ps_merged3 <- subset(ps_merged2, genome %in% filtered_list)
ps_merged3$norm <- ps_merged3$Abundance/ps_merged3$SUM
phy_filter4 <- subset_taxa(phy_filter3, genome %in% filtered_list)
saveRDS(phy_filter4, file='phy_filter4.rds')

# f_tab <- ddply(ps_WT_subset, .(genome, Sampling.time, T1, Sample, Orig_name), summarise, MEAN=mean(norm)) #Avg coverage of a genome
# df_wide <- f_tab2 %>%
#   pivot_wider(
#     id_cols = c(genome, T1),
#     names_from = Orig_name,
#     values_from = c(MEAN2, SE2),
#     names_glue = "{Orig_name}_{.value}"
#   )

## ============================================================
## PCoA
## ============================================================
plant_col <- c(
  "#E69F00", # Orange
  "#56B4E9", # Sky Blue
  "#009E73", # Bluish Green
  "#0072B2", # Blue
  "#D55E00")  # Vermillion
p_list <- c("WT (CP88-1762)", "1565", "1566", "1569", "1580")
plant_col2 <- setNames(plant_col, p_list)

ord <- ordinate(
  phy_filter4,
  method = "PCoA",
  distance = "bray"
)

p1 <- plot_ordination(phy_filter4, ord) +
  geom_point(
    aes(
      color = Orig_name,
      shape = Sampling.time
    ),
    size = 8
  ) +
  theme_bw() +
  scale_color_manual(
    values = plant_col2,
    name = "Plant genotype"
  ) +
  scale_shape_discrete(
    name = "Sampling time"
  )


ggsave("./results/figures/pcoa-test.png", plot = p1, width = 6, height = 5, dpi = 300)


bray <- phyloseq::distance(phy_filter4, method = "bray")
adonis_result <- adonis2(bray ~ Genotype*Sampling.time, data =  data.frame(sample_data(phy_filter4)), by = "terms")
print(adonis_result)
adonis2(formula = bray ~ Genotype * Sampling.time, data = data.frame(sample_data(phy_filter4)), by = "terms")

#                        Df SumOfSqs      R2      F Pr(>F)
# Genotype                4   1.7345 0.22285 2.3120  0.001 ***
# Sampling.time           1   1.2947 0.16634 6.9029  0.001 ***
# Genotype:Sampling.time  4   1.0029 0.12885 1.3368  0.064 .
# Residual               20   3.7511 0.48195
# Total                  29   7.7832 1.00000
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#load(file="heatmap_ps_merged2.rds") #This was previoulsy ps_merged2



## ============================================================
## Heatmap
## ============================================================
library(dplyr)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)


ps_merged4 <- psmelt(phy_filter4)
f <- ddply(ps_merged3, .(Sample, genome, Genotype), summarise, MEAN=mean(norm), SE = sd(norm)/sqrt(length(norm)))
for_hp <- f %>% select(Sample, genome, MEAN)
for_hp2 <- spread(for_hp, key = Sample, value = MEAN)
mags <- for_hp2$genome
matrix_data <- for_hp2[, -1]
matrix_data2 <- as.matrix(matrix_data)
rownames(matrix_data2) <- mags
matrix_data3 <- log10(matrix_data2 + 1e-5)
mat_col <- colnames(matrix_data3)
cn <- colnames(matrix_data3)
#cn_updated <- ifelse(cn %in% names(mapping_vec), mapping_vec[cn], cn)
#colnames(matrix_data3) <- cn_updated


tol21 <- c(
  "#4477AA", "#117733", "#DDCC77", "#CC6677", "#88CCEE",
  "#AA4499", "#44AA99", "#999933", "#882255", "#661100",
  "#6699CC", "#888888", "#AA4466", "#DDDDDD", "#77AADD",
  "#EE8866", "#EEDD88", "#FFAABB", "#99DDFF", "#44BB99", "#AAAA00"
)


mapping <- meta %>% select(V3, Orig_name)
ann <- meta %>% select(V3, Growing.condition)
cn_updated <- ifelse(cn %in% names(ann), ann[cn], cn)
ann_colors <- list(
  Group = c("Field" = "darkgreen",
            "Greenhouse" = "blueviolet"))

order_of_col <- c("Ga0501033", "Ga0501032", "Ga0501031", "Ga0501037", "Ga0501039", "Ga0501038", "Ga0501040", "Ga0501041", "Ga0501042", "Ga0501045", "Ga0501041", "Ga0501044", "Ga0501048", "Ga0501046", "Ga0501047", 
                  "Ga0501067", "Ga0501068", "Ga0501069", "Ga0501073", "Ga0501074", "Ga0501075", "Ga0501076", "Ga0501077", "Ga0501078", "Ga0501079", "Ga0501080", "Ga0501081", "Ga0501082", "Ga0501083", "Ga0501084")


mapping_vec <- setNames(mapping$Orig_name, mapping$V3)
cn_updated <- ifelse(cn %in% names(mapping_vec), mapping_vec[cn], cn)
cn_updated_man <- c("WT (CP88-1762) GH", "WT (CP88-1762) GH", "WT (CP88-1762) GH", "1565 GH", "1565 GH", "1565 GH", "1566 GH", "1566 GH", "1566 GH", 
                    "1569 GH", "1569 GH", "1569 GH", "1580 GH", "1580 GH", "1580 GH", "WT (CP88-1762) F", 
                    "WT (CP88-1762) F", "WT (CP88-1762) F", "1565 F", 
                    "1565 F", "1565 F", "1566 F", "1566 F", "1566 F", "1569 F", "1569 F", "1569 F", 
                    "1580 F", "1580 F", "1580 F")
colnames(matrix_data3) <- cn_updated_man
#rn <- rownames(matrix_data3[out$tree_row[["order"]],])
rn <- rownames(matrix_data3)
df1 <- as.data.frame(rn)
mag_ann <- read.csv(file="exported_img_data.tsv", sep="\t")
df2 <- merge(df1, mag_ann, by.x = "rn", by.y = "bin_oid")
annotation_r <-as.data.frame(df2$T1)
rownames(annotation_r) <- rn
colnames(annotation_r) <- "T1"
annotation_colors <- list(T1 = setNames(tol21, unique(annotation_r$T1)))

phylum_map <- c(
  # Major bacterial phyla
  " Proteobacteria"   = "Pseudomonadota",
  " Firmicutes"       = "Bacillota",
  " Bacteroidetes"    = "Bacteroidota",
  " Actinobacteria"   = "Actinobacteriota",
  " Chloroflexi"      = "Chloroflexota",
  " Acidobacteria"    = "Acidobacteriota",
  " Verrucomicrobia"  = "Verrucomicrobiota",
  " Cyanobacteria"    = "Cyanobacteriota",
  " Nitrospirae"      = "Nitrospirota",
  " Spirochaetes"     = "Spirochaetota",
  " Fusobacteria"     = "Fusobacteriota",
  " Tenericutes"      = "Mycoplasmatota",
  # Archaeal examples (optional, include if relevant)
  " Euryarchaeota"    = "Euryarchaeota",
  " Thaumarchaeota"  = "Thermoproteota",
  "Crenarchaeota"   = "Thermoproteota"
)

annotation_r$T1_new <- ifelse(
  annotation_r$T1 %in% names(phylum_map),
  phylum_map[annotation_r$T1],
  annotation_r$T1
)

annotation_r$T1 <- annotation_r$T1_new
annotation_r$T1_new <- NULL
annotation_colors <- list(T1 = setNames(tol21, unique(annotation_r$T1)))

out <- pheatmap(matrix_data3, color = brewer.pal(10, 'RdYlGn'), cluster_cols = FALSE,  colnames_col = order_of_col, clustering_distance_rows = 'euclidean', clustering_method = 'average', fontsize_row = 3, annotation_row = annotation_r, annotation_colors = annotation_colors, show_rownames=FALSE)
ggsave(out, file='./results/figures/heatmap.png', width=10, height=15, dpi=300)
save(df2, file='./results/rds/df2.rds')
