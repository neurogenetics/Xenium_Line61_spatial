library(ggplot2)
library(tidyverse)
library(dplyr)
library(circlize)
library(extrafont)
library(svglite)

################################################################################

# figure to show DE for different cell types for percentiles analysis (Fig. 5c)

genes_to_keep <- read.csv2("/xenium_Line61_final/xenium_upload1/general/genes_to_keep_DE.csv", header = F)
genes_to_keep <- genes_to_keep$V1
res_Tg_pos_vs_neg.df <- readRDS("/xenium_Line61_final/ctx/outs/DE_ctx_Tg_pSyn_pos_neg_percentiles.rds")
res_Tg_pos_vs_neg.df <- res_Tg_pos_vs_neg.df[res_Tg_pos_vs_neg.df$gene %in% genes_to_keep, ]

sigs.df <- res_Tg_pos_vs_neg.df[res_Tg_pos_vs_neg.df$padj <= 0.05, ]
sigs.df <- na.omit(sigs.df)

sigs.df$logFDR <- -log10(sigs.df$padj)

sigs.df$shape <- ifelse(sigs.df$gene == "hSNCA", "hSNCA", 
                        ifelse(sigs.df$gene == "Plk2", "Plk2", "Other"))

sigs.df$shape <- factor(sigs.df$shape, levels = c("hSNCA", "Plk2", "Other"))

sigs.df$celltype <- factor(sigs.df$celltype, 
                           levels = rev(c("L2/3 IT ExN", "L4/5 IT ExN", "L5 IT ExN", "L5 ET ExN", 
                                          "L5 NP ExN", "L6 IT ExN", "L6 CT ExN", "L6b ExN")))


p2 <- ggplot(sigs.df, aes(x = log2FoldChange, y = celltype)) +
  geom_point(aes(size = logFDR, fill = log2FoldChange, shape = shape)) +
  scale_fill_gradient2(low = "#05409e", mid = "white", high = "#e64a02") +
  scale_shape_manual(values = c(23, 22, 21)) + 
  scale_size_continuous(range = c(1, 8)) + 
  theme_bw() +
  labs(x = bquote(~Log[2]~ 'fold change'),
       y = "Cell type",
       title = "pSyn+ / hSNCA-low vs. pSyn- / hSNCA-high") +
  theme(text = element_text(family ="Arial"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 10, face = "bold")) +
  guides(size = guide_legend(title = bquote(~-Log[10]~'FDR')),
         fill = guide_colourbar(title = bquote(~Log[2]~ 'fold change')),
         shape = guide_legend(title = "Gene"))


svglite("/xenium_Line61_final/ctx/outs/percentiles_DE_dotplot.svg", 
        width = 8, height = 4.5)
plot(p2)
dev.off()
