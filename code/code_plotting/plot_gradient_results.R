
library(tidyverse)
library(ggrepel)
library(plotly)


net_cols <- c("#77AADD", "#CC6677", "#44BB99", "#EEDD88")

output <- '/home/mgell/Work/FC/plots/'

g_vals <- read_csv('/home/mgell/Work/FC/text_files/res/gradient_values_pca_80_zscore.csv', col_names = c('G1','G2','G3'))
tab <- read_delim('/home/mgell/Work/FC/nets/All_localmax.node', delim = '\t', col_names = c('x','y','z','net','netn','roi'))
tab <- tab[c(-2,-17,-24),]
tab$net[2] <- 'I+SV'

net_communities <- read_csv('/home/mgell/Work/FC/hcp/text_files/res/replic_5mm_net_communities.csv', col_names = 'netn')
# net_communities$netn[net_communities$netn == 1] = 'Delayed con sensitivity'
# net_communities$netn[net_communities$netn == 2] = 'Salience'
# net_communities$netn[net_communities$netn == 3] = 'Temporoparietal'
# net_communities$netn[net_communities$netn == 4] = 'Frontoparietal'
tab$netn <- net_communities$netn


g_vals$roi <- tab$roi
g_vals$netn <- as.factor(tab$netn)


# g_vals$G1 <- g_vals$G1 + rnorm(length(g_vals$G1), mean=0, sd=0.0005)
# g_vals$G2 <- g_vals$G2 + rnorm(length(g_vals$G2), mean=0, sd=0.005)
# g_vals$G3 <- g_vals$G3 + rnorm(length(g_vals$G3), mean=0, sd=0.005)

fig <- plot_ly(g_vals, x = ~G1, y = ~G2, z = ~G3, color = ~netn, colors = net_cols)#, width = 500, height = 500)

fig <- fig %>% add_markers()

#fig <- fig %>% layout(showlegend = FALSE)

#fig




fig1 <- ggplot(g_vals, aes(G1, G2)) +
  geom_point(aes(col=netn), size = 4) + #position = position_jitter(width = .002, height = 0.002), size = 4) +
  scale_colour_manual(values = net_cols, name = "Community",
                      labels = c("Fronto-medial", "Cingulo-insular", "Frontoparietal", "Temporoparietal")) +
  theme_classic()

#ggsave(paste0(output, 'pca_gradients_g12_no_labels.png'),fig1,width = 6,height = 4)

fig1 <- fig1 + geom_text_repel(aes(label=roi), max.overlaps = 15)

#ggsave(paste0(output, 'dm_gradients_g12.png'),fig1,width = 6,height = 4)
ggsave(paste0(output, 'pca_gradients_g12.png'),fig1,width = 6,height = 4)



fig2 <- ggplot(g_vals, aes(G1, G3)) +
  geom_point(aes(col=netn), size = 4) + #position = position_jitter(width = .0035), size = 4) +
  scale_colour_manual(values = net_cols, name = "Community",
                      labels = c("Delayed con. sensitivity", "Salience", "Frontoparietal", "Temporoparietal")) +
  theme_classic()

#ggsave(paste0(output, 'pca_gradients_g13_no_labels.png'),fig2,width = 6,height = 4)

fig2 <- fig2 + geom_text_repel(aes(label=roi), max.overlaps = 20)

ggsave(paste0(output, 'pca_gradients_g13.png'),fig2,width = 6,height = 4)



fig3 <- ggplot(g_vals, aes(G3, G2)) +
  geom_point(aes(col=netn), size = 4) +
  scale_colour_manual(values = net_cols, name = "Community",
                      labels = c("Delayed con. sensitivity", "Salience", "Frontoparietal", "Temporoparietal")) +
  theme_classic()

#ggsave(paste0(output, 'pca_gradients_g23_no_labels.png'),fig3,width = 6,height = 4)

fig3 <- fig3 + geom_text_repel(aes(label=roi), max.overlaps = 20)

ggsave(paste0(output, 'pca_gradients_g23.png'),fig3,width = 6,height = 4)

