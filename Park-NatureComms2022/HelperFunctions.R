set.seed(625)
theme_pub = theme_linedraw() + 
  theme(strip.background=element_rect(fill="grey80", colour="grey50", size=0.2), 
        strip.text.x=element_text(colour="black"), 
        strip.text.y=element_text(colour="black")) + 
  theme(
    legend.position="right", 
    legend.title=element_text(size=15), legend.text=element_text(size=14), 
    axis.text.x = element_text(size=12, angle=-90, hjust=0, vjust=0.5), 
    axis.text.y=element_text(size=12), axis.title=element_text(size=15), 
    axis.title.y=element_text(vjust=1), plot.title = element_text(size=18, vjust=1.5), 
    strip.background = element_rect(fill="#EEEEEE"), strip.text = element_text(size = 11), 
    panel.grid.major = element_line(colour = "grey98"), panel.grid.minor = element_blank())

# functions
volcano_plot <- function(df0 = voldf, subdf = voldf_subset, plotname = "example plot name") {
  df <- df0[complete.cases(df0$padj),]
  ggplot(df, aes(x = log2FoldChange, y = -log10(padj))) + 
    geom_point(data = df, color = "gray50", alpha = 0.2) +
    geom_point(data = subset(df, log2FoldChange > FC_cutoff & padj < adjp_cutoff), color = "darkred", alpha = 0.3) + 
    geom_point(data = subset(df, log2FoldChange < -(FC_cutoff) & padj < adjp_cutoff), color = "navy", alpha = 0.3) + 
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = -log10(adjp_cutoff), linetype = "dashed") +
    geom_vline(xintercept = c(-(FC_cutoff), FC_cutoff), linetype = "dashed") + 
    theme_linedraw() + 
    theme(panel.grid = element_blank(), legend.position = "none", plot.title = element_text(hjust = 0.5)) +
    xlab("log2 Fold Change") + ylab("-log10 adjusted P-Value") + 
    scale_fill_manual(values = c("down" = "navy", "up" = "darkred")) + 
    ggrepel::geom_label_repel(data = subdf, 
                              aes(label = GENENAME, x = log2FoldChange, y = -log10(padj), fill = as.factor(direction)),
                              max.overlaps = 25,
                              box.padding = unit(.3, "lines"), color = "white", segment.color = "gray30", size = 3, alpha = 0.9) + 
    ggtitle(plotname)
}
process_voldf <- function(df = res) {
  res1 <- df
  res1$GENENAME <- rownames(res1)
  vol_subdf <- rbind(res1[((res1$log2FoldChange > FC_cutoff) & (res1$padj < adjp_cutoff)),],
                     res1[((res1$log2FoldChange < (-(FC_cutoff))) & (res1$padj < adjp_cutoff)),])
  vol_subdf <- vol_subdf %>% mutate(direction = case_when(log2FoldChange > 0 ~ "up", log2FoldChange < 0 ~ "down"))
  return(list(res1, vol_subdf))
}

run_volcano <- function(genedf = DEG_OE, dfname = "OE"){
  processed_df <- genedf %>% process_voldf()
  subdf <- processed_df[[2]] %>% as.data.frame()
  subdf <- subdf[complete.cases(subdf$pvalue),]
  subdf_degs <- c(nrow(subdf[subdf$direction == "up",]), nrow(subdf[subdf$direction == "down",]))
  volcano_plot(df = processed_df[[1]], subdf = processed_df[[2]], 
               plotname = paste("Volcano plot: ", "overall Post vs. Pre flight", ", |FC| >", FC_cutoff, " and adj.p-value <", adjp_cutoff, sep="")) + 
    labs(subtitle = paste(dfname, " (", subdf_degs[1], " up-regulated and ", subdf_degs[2], " down-regulated genes.)", sep="")) + 
    theme(plot.subtitle = element_text(hjust = 1))}

pathway_analysis <- function(df = dfname, pathway = pathways.hallmark) {
  stats <- data.frame(name = df$GENENAME, stats = df$stat)
  stats1 <- tibble::deframe(stats)
  fgseaRes <- fgsea(pathways = pathway, stats = stats1, maxSize = 500, eps = 0)
  fgseaRes <- fgseaRes %>% dplyr::select(-leadingEdge, -ES) %>% arrange(pval) %>% as.data.frame()
  #fgseaRes <- fgseaRes %>% top_n(50, pval)
  return(fgseaRes)
}
plot_pathways <- function(fgseaRes = fgseaRes, p_cut = 0.1) {
  plotdf <- fgseaRes %>% as_tibble() %>% arrange(pval) %>% dplyr::filter(pval <= p_cut)
  ggplot(plotdf, aes(reorder(pathway, NES), NES)) + 
    geom_col(aes(fill=-log10(padj))) +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title="Pathway NES from GSEA") + 
    theme_bw() + 
    scale_fill_gradient2(midpoint=-log10(p_cut), low="navy", mid="gray", high="darkred") +
    coord_flip()
}