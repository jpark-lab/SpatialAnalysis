library(ggplot2)
library(ggcorrplot)
library(ggpubr)
library(ggrepel)
library(ComplexHeatmap)

load("Fig3_data1.RData")
# contains dataframes below: 
# geomx_data
# geomx_avgprop
# geomx_avgprop_zoom
# statistics_normal_vs_highcov
# statistics_normal_vs_lowcov
# statistics_highcov_vs_lowcov
# entropy_ct1
# entropy_comb1
# entropy_df
# ssgseares_T
# ssgseares_neutrophil
# ssgseares_macrophage
# volcano_plot function (shown below)
# plot_ssgsea function (shown below)

# fig 3A: correlation by cell types
b <- as.character(unique(geomx_data$diseases))
df_forcorplot <- data.frame(geomx_data[,c(19,17,18,15,13,16,11,20,14,8,9,7,12,6)]) # cell type order
corplot_bydiseases <- vector('list', length(b))
for(i in 1:5) {
  corplot_bydiseases[[i]] <- ggcorrplot(cor((df_forcorplot[geomx_data$diseases == b[i],])), 
                                        hc.order = FALSE, 
                                        colors = c("blue", "grey", "red"),
                                        insig = "blank",
                                        p.mat = cor_pmat((df_forcorplot[geomx_data$diseases == b[i],]))) + 
    theme_pubclean() + 
    ggtitle(paste(b[i])) +
    theme(axis.text.x = element_text(angle=45, hjust=1), 
          axis.title.x = element_blank(), 
          axis.title.y = element_blank())
}
ggarrange(plotlist = corplot_bydiseases, nrow = 2, ncol = 3)

# fig 3 supplementary: quantification of correlations 
volcano_plot <- function(df = statistics_normal_vs_highcov, plotname = "Normal_vs_COVID-19_High") {
  ggplot(df, aes(x = zstat, y = -log10(pvalue))) + 
    theme_pubr() + 
    scale_fill_gradient(low = "lightgray", high = "navy") +
    scale_color_gradient(low = "lightgray", high = "navy") +
    geom_point(data = df, color = "grey", alpha = 0.5, size = 2) +
    geom_point(data = subset(df, zstat > 2), color = "red", alpha = 0.5) + 
    geom_point(data = subset(df, zstat < -2), color = "navy", alpha = 0.5) + 
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    ggtitle(paste(plotname)) +
    theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5)) +
    xlab("Fisher transformed z") +
    ylab("-log10 P-Value") +
    geom_label_repel(data = subset(df, (zstat>2&-log10(pvalue)>2.5)), aes(label = name),
                     box.padding   = 0.35, point.padding = 0.5, segment.color = 'darkred') +
    geom_label_repel(data = subset(df, (zstat<(-2)&-log10(pvalue)>2.5)), aes(label = name),
                     box.padding   = 0.35, point.padding = 0.5, segment.color = 'navy')
}
volcano_plot(df = statistics_normal_vs_highcov, plotname = "Normal vs. Covid-19 High")
volcano_plot(df = statistics_normal_vs_lowcov, plotname = "Normal vs. Covid-19 Low")
volcano_plot(df = statistics_highcov_vs_lowcov, plotname = "Covid-19 Low vs. Covid-19 High")

# fig 3B: average proportions
ggplot(geomx_avgprop, aes(x = celltypes, y = value, group = disease)) +
  scale_fill_brewer(palette = "Set2") +
  ylab("average proportions (%)") + 
  geom_bar(stat = "identity", position = "dodge", aes(fill = disease)) + 
  geom_errorbar(stat = "identity", color = "darkgrey", position = "dodge", width = 0.9, aes(ymin = value, ymax = value + sd)) +
  theme_pubr() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        panel.grid.major.x = element_line(size = 10, linetype = "solid", colour = "grey90"), 
        panel.grid.major.y = element_blank())
# zoom in version of small proportions
geomx_avgprop_zoom <- geomx_avgprop[(geomx_avgprop$celltypes == "mast" | geomx_avgprop$celltypes == "plasma" | geomx_avgprop$celltypes == "B" | geomx_avgprop$celltypes == "pDC"),]
ggplot(geomx_avgprop_zoom, aes(x = celltypes, y = value, group = disease)) +
  scale_fill_brewer(palette = "Set2") +
  ylab("average proportions (%)") + 
  geom_bar(stat = "identity", position = "dodge", aes(fill = disease)) + 
  geom_errorbar(stat = "identity", color = "darkgrey", position = "dodge", width = 0.9, aes(ymin = value, ymax = value + sd)) +
  theme_pubr() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# fig 3C: Entropy calculations
# within the tissue type
ComplexHeatmap::Heatmap(entropy_ct1, cluster_columns = T, cluster_rows = F,
                        rect_gp = gpar(col = "white", lwd = 2), 
                        column_title = "Entropy estimates")
# within the cell types
# I removed plasma cells because the normal was insignificant, but can keep it as well
ComplexHeatmap::Heatmap(entropy_comb1, cluster_columns = F, 
                        rect_gp = gpar(col = "white", lwd = 2), 
                        column_title = "Entropy estimates, total")

# fig 3 supplementary: tissue-type specific entropy visualizations
ht1 <- ComplexHeatmap::Heatmap(as.matrix(entropy_df[[1]]), cluster_columns = F, 
                               rect_gp = gpar(col = "white", lwd = 2), 
                               column_title = "Large Airway")
ht2 <- ComplexHeatmap::Heatmap(as.matrix(entropy_df[[2]]), cluster_columns = F, 
                               rect_gp = gpar(col = "white", lwd = 2), 
                               show_heatmap_legend = F,
                               column_title = "Alveolar") 
ht3 <- ComplexHeatmap::Heatmap(as.matrix(entropy_df[[3]]), cluster_columns = F, 
                               rect_gp = gpar(col = "white", lwd = 2),
                               show_heatmap_legend = F,
                               column_title = "Vascular")
htlist <- ht1 + ht2 + ht3
draw(htlist, column_title = "Entropy Calculations (ML) for cell-type counts, relative to normal")

# fig 3 supplementary: ssgsea results
col_fun = colorRamp2(c(-1.5, 0, 1.5), c("#6D9EC1", "white", "#E46726"))
plot_ssgsea <- function(df = ssgseares_macrophage) {
  ComplexHeatmap::Heatmap(df %>% scale() %>% t(),
                          show_column_names = TRUE, 
                          show_row_dend = T, 
                          cluster_rows = T,
                          cluster_columns = T, 
                          na_col = "grey", col = col_fun, 
                          row_names_gp = gpar(fontsize = 10), 
                          row_names_max_width = unit(15, "cm")) %>% 
    draw(heatmap_legend_side = "left")
}
plot_ssgsea(ssgseares_T)
plot_ssgsea(ssgseares_neutrophil)
plot_ssgsea(ssgseares_macrophage)
