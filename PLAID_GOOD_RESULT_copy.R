require("biclust")
require("RColorBrewer")
layout(matrix(seq(1:6),1,6))

data_with_title = read.table("/Users/jgu1/WorkSpace/lab/curr/clustering/cluster_input_q0.5_log10p.txt", fill=TRUE,header=TRUE)
GWAS_names = colnames(data_with_title)[-1]# -1 means remove the first element
#gene_names = data_with_title$gwas[-1]
gene_names = data_with_title$gwas
data_with_title_no_gwas = subset(data_with_title,select = -gwas)

full_matrix = data.matrix(data_with_title_no_gwas)
#replace missing value with very small value
full_matrix[is.na(full_matrix)] = 1e-08
rownames(full_matrix) = gene_names


bics <- biclust(full_matrix, method=BCPlaid(), cluster="b", fit.model = y ~ m + a + b,
         background = FALSE, row.release = 0.7, 
         col.release = 0.7, shuffle = 3, back.fit = 0, max.layers = 20, iter.startup = 5,
         iter.layer = 10, verbose = TRUE)

#my_palette <- colorRampPalette(c("#ffeda0", "#feb24c", "#f03b20"),bias=0.8)(n = 299)
my_palette <- colorRampPalette(c("#000000","#800000","#ff0000"),bias=100)(n = 299)
#my_palette <- colorRampPalette( colors = rev(brewer.pal(9,"Reds")) )
#my_palette = heat.colors(100, alpha = 1)

#full_matrix = ifelse(full_matrix > 4, 4, full_matrix)

for(i in seq(1:bics@Number)){
  drawHeatmap(full_matrix,bics,i, beamercolor=TRUE, paleta = my_palette,local=TRUE)
}
