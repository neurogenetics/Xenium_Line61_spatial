library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(broom)
library(ggrepel)
library(extrafont)
library(svglite)


xenium_hip <- readRDS("/xenium_Line61_final/hip/xenium_HIP_metadata_pSyn.rds")

################################################################################

meta.df <- as.data.frame(xenium_hip@meta.data)

################################################################################

# number of cells per sample (Supplementary Fig. 5a)

num_cells.df <- meta.df %>%
  group_by(orig.ident) %>%
  summarise(cell_count = n())

num_cells.df$orig.ident <- factor(num_cells.df$orig.ident, levels = c("nTg_1", "nTg_2", "nTg_3", "nTg_4", 
                                                                      "Tg_1", "Tg_2", "Tg_3", "Tg_4")) 

cols <- c("blue", "blue", "blue", "blue", "red", "red", "red", "red")

p1 <- ggplot(num_cells.df, aes(x = orig.ident, y = cell_count, fill = orig.ident)) +
  geom_bar(stat = 'identity') + 
  labs(x = "",
       y = "# cells", 
       title = "# cells per sample") +
  theme_bw() + 
  scale_fill_manual(values = cols) +
  guides(fill = "none") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12, face = "bold"),
        text = element_text(family = "Arial")) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 21000)) 

svglite("/xenium_Line61_final/hip/outs/xenium_hip_numcellspersample.svg", 
        width = 6, height = 4.5)
plot(p1)
dev.off()

################################################################################

# number of transcripts per cell by sample (Supplementary Fig. 5b)

numts.df <- meta.df[c("nCount_Xenium", "orig.ident")]

cols <- c("blue", "blue", "blue", "blue", "red", "red", "red", "red")

p2 <- ggplot(numts.df, aes(x = nCount_Xenium, fill = orig.ident)) + 
  geom_histogram(binwidth = 10) +
  facet_wrap(~orig.ident, nrow = 2) +
  labs(title = "Number of transcripts per cell",
       y = "Number of cells",
       x = "Number of transcripts") +
  scale_fill_manual(values = cols) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        text = element_text(family = "Arial"),
        legend.position = "none",
        strip.text = element_text(face = "bold"))

svglite("/xenium_Line61_final/hip/outs/xenium_hip_numtranscripts.svg", 
        width = 10.5, height = 5.5)
plot(p2)
dev.off()

################################################################################

# number of genes per cell per sample (Supplementary Fig. 5c)

numgenes.df <- meta.df[c("nFeature_Xenium", "orig.ident")]

cols <- c("blue", "blue", "blue", "blue", "red", "red", "red", "red")

p3 <- ggplot(numgenes.df, aes(x = nFeature_Xenium, fill = orig.ident)) + 
  geom_histogram(binwidth = 5) +
  facet_wrap(~orig.ident, nrow = 2) +
  labs(title = "Number of unique genes per cell",
       y = "Number of cells",
       x = "Number of unique genes") +
  scale_fill_manual(values = cols) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        text = element_text(family = "Arial"),
        legend.position = "none",
        strip.text = element_text(face = "bold"))

svglite("/xenium_Line61_final/hip/outs/xenium_hip_numgenes.svg", 
        width = 10.5, height = 5.5)
plot(p3)
dev.off()

################################################################################

# negative probe rate per sample (Supplementary Fig. 5d)

negprobe.df <- meta.df[, c("orig.ident", "nCount_Xenium", "nCount_BlankCodeword")]

negprobe.df <- negprobe.df %>%
  mutate(negproberate = nCount_BlankCodeword / nCount_Xenium * 100)

cols <- c("blue", "blue", "blue", "blue", "red", "red", "red", "red")

p4 <- ggplot(negprobe.df, aes(x = negproberate, fill = orig.ident)) + 
  geom_histogram(binwidth = 0.2) +
  facet_wrap(~orig.ident, nrow = 2) +
  labs(title = "Negative probe rate per cell",
       y = "Number of cells",
       x = "% negative control probes") +
  scale_fill_manual(values = cols) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        text = element_text(family = "Arial"),
        legend.position = "none",
        strip.text = element_text(face = "bold")) 

svglite("/xenium_Line61_final/hip/outs/xenium_hip_negproberate.svg", 
        width = 10.5, height = 5.5)
plot(p4)
dev.off()
