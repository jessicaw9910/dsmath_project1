library(magrittr)
library(ggplot2)
library(ggrepel)
library(reshape)
library(dplyr)
library(scales)
library(stringr)
library(svglite)

df_drug_info <- read.csv('../assets/drug_data.csv', stringsAsFactors = FALSE)
df_cell_info <- read.csv('../assets/cell_list.csv', stringsAsFactors = FALSE)

df_drug <- read.csv('../assets/drug.csv', stringsAsFactors = FALSE)
df_cell <- read.csv('../assets/cell.csv', stringsAsFactors = FALSE)

df_drug_lrm <- read.csv('../assets/drug_lrm.csv', stringsAsFactors = FALSE)
df_cell_lrm <- read.csv('../assets/cell_lrm.csv', stringsAsFactors = FALSE)

#########################
#### Clean Drug Data ####
#########################

## create index for later use
df_drug$idx <- seq(1, length(df_drug$drug))

## two entries require manual fixes
idx <- which(df_drug$drug %in% "KRAS\nG12C Inhibitor-12\n1855")
df_drug$drug[idx] <- "KRAS (G12C) Inhibitor-12\n1855"
idx <- which(df_drug$drug %in% "Nutlin-3a\n-\n1047")
df_drug$drug[idx] <- "Nutlin-3a (-)\n1047"

## extract drug name alone
df_drug$drug_name <- sapply(df_drug$drug,
                            function(x) gsub("[\r\n].*", "", x))

## create column with log10(1/IC50)
df_drug$inv_IC50 <- log10(1/exp(df_drug$lnIC50_mean))

## start cluster at 1
df_drug$knn_all <- df_drug$knn_all + 1

## merge with drug pathway info
df_drug <- merge(df_drug, df_drug_info[, c('drug_name', 'pathway_name')], 
                 by = 'drug_name', all.x = TRUE)

## remove duplicates
df_drug <- df_drug[-which(duplicated(df_drug$drug)), ]

#########################
#### Drug Data Plots ####
#########################

df_drug %>% ggplot(aes(x=tsne1, y=tsne2, 
                       color=pathway_name, shape=factor(knn_all))) +
  geom_point(aes(size = inv_IC50)) +
  geom_text_repel(size = 3, aes(label = toupper(drug_name)), 
                  point.padding = unit(0.3, "lines"), show.legend = FALSE) +
  scale_shape_manual(values = 0:10)
  labs(title = "Relative size indicates log10(1/Average IC50)") +
  xlab("PC1") +
  ylab("PC2") +
  scale_size(guide="none") +
  guides(fill=guide_legend(title="New Legend Title")) +
  guides(color=guide_legend(title="Drug Pathway", ncol=1, 
                            override.aes = list(size=6)))
  
#############################
#### Low-Rank Drug Model ####
#############################
  
dim(df_drug_lrm)

colnames(df_drug_lrm) <- paste(rep('V', 10), 
                               seq(1, dim(df_drug_lrm)[2]), sep = "")
df_drug_lrm$idx <- seq(1, dim(df_drug_lrm)[1])

df_drug_top10 <- as.data.frame(apply(df_drug_lrm[, 2:11], 2, 
                                     function(x) order(x, decreasing = TRUE)))
df_drug_top10$rank <- seq(1, dim(df_drug_top10)[1])
df_drug_top10 <- melt(df_drug_top10, id.vars = 'rank')
colnames(df_drug_top10) <- c('rank', 'V', 'idx')
df_drug_top10 <- merge(df_drug_top10, df_drug[, c('idx', 'inv_IC50', 'drug_name', 'pathway_name')],
                       by = 'idx', all.x = TRUE)

idx_top <- which(df_drug_top10$rank <= 10)
idx_bottom <- which(df_drug_top10$rank >= 189)

ggplot(df_drug_top10[idx_top, ], aes(x=rank, y=V, size=inv_IC50, color=pathway_name)) +
  geom_point()

ggplot(df_drug_top10[idx_bottom, ], aes(x=rank, y=V, size=inv_IC50, color=pathway_name)) +
  geom_point()

ggplot(df_drug_top10[idx_top, ], aes(x=V, y=inv_IC50)) + 
  geom_boxplot() +
  geom_jitter(aes(colour=pathway_name, size=inv_IC50), alpha=0.5)

ggplot(df_drug_top10[idx_bottom, ], aes(x=V, y=inv_IC50)) + 
  geom_boxplot() +
  geom_jitter(aes(colour=pathway_name, size=inv_IC50), alpha=0.5)

# df_stacked_top <- df_drug_top10[idx_top, ] %>% 
#   count(V, pathway_name, sort = TRUE)
# df_stacked_bottom <- df_drug_top10[idx_bottom, ] %>% 
#   count(V, pathway_name, sort = TRUE)

df_stacked_top <- df_drug_top10[idx_top, ] %>%
  group_by(V, pathway_name) %>%
  summarise(all_names = paste(drug_name, collapse = ", "))
df_stacked_top$n <- str_count(df_stacked_top$all_names, ",") + 1

df_stacked_bottom <- df_drug_top10[idx_bottom, ] %>%
  group_by(V, pathway_name) %>%
  summarise(all_names = paste(drug_name, collapse = ", "))
df_stacked_bottom$n <- str_count(df_stacked_bottom$all_names, ",") + 1

ggplot(df_stacked_top, aes(x = V, y = n, fill = pathway_name)) + 
  geom_bar(stat = "identity")

ggplot(df_stacked_bottom, aes(x = V, y = n, fill = pathway_name)) + 
  geom_bar(stat = "identity")

df_stacked_top$order <- "Top 10"
df_stacked_bottom$order <- "Bottom 10"

df_stacked <- rbind(df_stacked_top, df_stacked_bottom)
df_stacked$percent <- df_stacked$n / 10
df_stacked <- df_stacked[order(df_stacked$n), ]

df_path <- data.frame(pathway_name = unique(df_drug_info$pathway_name))
df_path$pathway_short <- c('PI3K/MTOR', 'Mitosis', 'DNA rep', 'Apoptosis',
                           'EGFR', 'Misc kinase', 'Unclassified', 'Other',
                           'Genome', 'CHA', 'IGF1R', 'MAPK', 'RTK',
                           'Cell cycle', 'Chromatin other', 'Cytoskeleton',
                           'JNK/p38', 'ABL', 'p53')

df_stacked <- merge(df_stacked, df_path, by = 'pathway_name', all.x = TRUE)
df_stacked$order <- relevel(as.factor(df_stacked$order), 'Top 10')

ggplot(df_stacked, aes(x = V, y = percent, fill = pathway_name)) + 
  geom_bar(stat = "identity") +
  facet_wrap(~order, dir = 'v') +
  geom_text(aes(label=stringr::str_wrap(all_names, 5)), size = 3, position = position_fill(vjust = 0.5))
  # geom_text(aes(label=paste(pathway_short, percent(percent, accuracy = 1), sep=": ")),
  #           size = 3, position = position_fill(vjust = 0.5)) +
  # theme(legend.position = "none") 

ggsave('../src/plots/Top10_Stacked Bar_Drug.svg', width = 12, height = 6)

ggplot(df_stacked, aes(x = V, y = percent, fill = pathway_name)) + 
  geom_bar(stat = "identity") +
  facet_wrap(~order, dir = 'v') +
  geom_text(aes(label=paste(pathway_short, percent(percent, accuracy = 1), sep=": ")),
            size = 3, position = position_fill(vjust = 0.5)) +
  theme(legend.position = "none") 