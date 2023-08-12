
library(tidyverse)
library(ggrepel)


#-----------------------SET UP--------------------------
net_cols <- c("#CC6677", "#77AADD", "#44BB99", "#EEDD88")

# load gradient values
g_vals <- read_csv('/home/mgell/Work/FC/hcp/text_files/res/replic_5mm_mean_all_FC_community_prop_21.csv')
g_vals <- g_vals[,7:9]

tab <- read_delim('/home/mgell/Work/FC/nets/All_localmax.node', delim = '\t', col_names = c('x','y','z','net','netn','roi'))
tab <- tab[c(-2,-17,-24),]

net_cummunities <- read_csv('/home/mgell/Work/FC/hcp/text_files/res/replic_5mm_net_communities.csv', col_names = 'netn')

g_vals$roi <- tab$roi
g_vals$net <- as.factor(net_cummunities$netn)

# WD
wd <- '/home/mgell/Work/FC/PET/'

#NAT_file = "NAT_iso2_tri.txt"
sert_file = "5HT1a_WAY_HC36_2x2x2_iso2_tri_FSL.txt"





### 5HT1a
# load
#receptor <- read_delim(paste0(wd, sert_file), delim = '  ', col_names = FALSE)
#receptor <- receptor[,c(-2,-24,-25)]
receptor <- read_delim(paste0(wd, sert_file), delim = '\t',  trim_ws = TRUE)
receptor <- receptor[,c(-22)]
g_vals$receptor <- as.numeric(receptor[])

# Spearman
print(sert_file)
c <- cor.test(g_vals$Z, g_vals$receptor, method = 'spearman')
print(c)

# plot
plt = ggplot(g_vals, aes(x = rank(Z), y = rank(receptor))) +
  geom_point(colour = 'skyblue2', size = 3) +
  geom_smooth(method = lm, alpha = 0.2, color = 'slategray') +
  theme_classic() +
  ylab(paste0('rank transformed 5HT1a density')) +
  xlab('rank transformed within-module degree z-score') +
  ylim(c(0,22)) + xlim(c(0,22))

ggsave(paste0('/home/mgell/Work/FC/PET/Association_plots/xhcp_', sert_file, '_by_Z_spearman_',
              as.character(round(c$estimate, digits = 2)), '.png'), plt, width = 4.5, height = 3.5)

plt <- plt + geom_text_repel(aes(label=roi), max.overlaps = 20)

ggsave(paste0('/home/mgell/Work/FC/PET/Association_plots/hcp_labelled_', sert_file, '_by_Z_spearman_',
              as.character(round(c$estimate, digits = 2)), '.png'), plt, width = 4.5, height = 3.5)


# Pearson 
c <- cor.test(g_vals$Z, g_vals$receptor, method = 'pearson')
print(c)

# plot
plt = ggplot(g_vals, aes(x = Z, y = receptor)) +
  geom_point(colour = 'skyblue2', size = 3) +
  geom_smooth(method = lm, alpha = 0.2, color = 'slategray') + 
  theme_classic() +
  ylab('5HT1a') +
  xlab('rank transformed within-module degree z-score') +
  #scale_y_continuous(breaks = c(5,8,11,14,17), limits = (c(4,18.5))) + xlim(c(0.35,0.73)) + 
  geom_text_repel(aes(label=roi), max.overlaps = 20)

ggsave(paste0('/home/mgell/Work/FC/PET/Association_plots/xhcp_labelled_', sert_file,'_by_Z_pearson_',
              as.character(round(c$estimate, digits = 2)), '.png'), plt, width = 4.5, height = 3.5)


### SPATIAL PERMUTATION
permpath <- '/home/mgell/Work/FC/PET/permutation/Sert5HT/'
files <- dir(permpath, pattern = 'Sert5HT*')
res <- rep(0,length(files))
i = 1

for (file_i in files) {
  nodes_i <- read_delim(paste0(permpath, file_i), delim = '\t', trim_ws = TRUE)
  nodes_i <- nodes_i[-22]
  
  # save
  res[i] <- cor(as.numeric(nodes_i[1,]), g_vals$ppos, method = 'spearman')
  i <- i+1
}

d <- data.frame('perm' = res)

# Get empirical correlation
empirical <- cor(g_vals$Z, g_vals$receptor, method = 'spearman')

sum(d$perm >= empirical) / 1000

sort_perm = sort(d$perm)

plt1 <- ggplot(d, aes(perm)) + geom_histogram(fill="skyblue3", alpha=0.5, bins = 40) +
  #scale_y_continuous(expand=c(0,0), limits = c(0,610), breaks = seq(0,600,100))+
  scale_y_continuous(limits = c(0,80), breaks = seq(0,80,20))+
  theme_classic()+
  theme(legend.position="none") + 
  geom_segment(linewidth= 1,
               aes(x = empirical, y = 0, 
                   xend = empirical, 
                   yend = 30, colour= "red")) +
  geom_segment(linewidth= 0.5, aes(x = sort_perm[950], y = 0, xend = sort_perm[950], yend = 60), colour = "gray60") +
  xlab('rho') + 
  ylab('count')

ggsave(paste0('/home/mgell/Work/FC/plots/hcp_5HT_perm_cor_Z.png'), plt1, width = 2, height = 2)






# ### NAT
# # load
# receptor <- read_delim(paste0(wd, NAT_file), delim = '  ', col_names = FALSE)
# receptor <- receptor[,c(-2,-24,-25)]
# g_vals$receptor <- as.numeric(receptor[])
# 
# # Spearman correlation
# print(NAT_file) #name
# c <- cor.test(g_vals$ppos, g_vals$receptor, method = 'spearman')
# print(c)
# 
# # Plot
# plt = ggplot(g_vals, aes(x = rank(ppos), y = rank(receptor))) +
#   geom_point(colour = 'skyblue2', size = 3) +
#   geom_smooth(method = lm, alpha = 0.2, color = 'slategray') + 
#   theme_classic() +
#   ylab(paste0('rank transformed NAT')) +
#   xlab('rank transformed participation coefficient') +
#   ylim(c(0,22)) + xlim(c(0,22))
# 
# ggsave(paste0('/home/mgell/Work/FC/PET/Association_plots/xhcp_', NAT_file, '_by_ppos_spearman_',
#               as.character(round(c$estimate, digits = 2)), '.png'), plt, width = 4.5, height = 3.5)
# 
# plt <- plt + geom_text_repel(aes(label=roi), max.overlaps = 20)
# ggsave(paste0('/home/mgell/Work/FC/PET/Association_plots/xhcp_labelled_', NAT_file, '_by_ppos_spearman_',
#               as.character(round(c$estimate, digits = 2)), '.png'), plt, width = 4.5, height = 3.5)
# 
# 
# # Pearson 
# c <- cor.test(g_vals$ppos, g_vals$receptor, method = 'pearson')
# print(c)
# 
# # plot
# plt = ggplot(g_vals, aes(x = ppos, y = receptor)) +
#   geom_point(colour = 'skyblue2', size = 3) +
#   geom_smooth(method = lm, alpha = 0.2, color = 'slategray') + 
#   theme_classic() +
#   ylab('NAT') +
#   xlab('participation coefficient') +
#   scale_y_continuous(breaks = c(5,8,11,14,17), limits = (c(5,18.5))) + xlim(c(0.05,0.65)) + 
#   geom_text_repel(aes(label=roi), max.overlaps = 20)
# 
# ggsave(paste0('/home/mgell/Work/FC/PET/Association_plots/xhcp_labelled_', NAT_file,'_by_ppos_pearson_',
#               as.character(round(c$estimate, digits = 2)), '.png'), plt, width = 4.5, height = 3.5)
# 
# 
# # Permutation
# #perm <- read_csv('/home/mgell/Work/FC/brainsmash/NAT_surrogate_maps.csv', col_names = FALSE)
# #surrogate_cors <- apply(perm, 1, function(x) cor(g_vals$ppos, x, method = 'spearman'))
# #d <- data.frame('perm' = surrogate_cors)
# 
# permpath <- '/home/mgell/Work/FC/PET/permutation/NAT/'
# files <- dir(permpath, pattern = 'NAT*')
# res <- rep(0,length(files))
# i = 1
# 
# for (file_i in files) {
#   nodes_i <- read_delim(paste0(permpath, file_i), delim = '\tab', trim_ws = TRUE)
#   nodes_i <- nodes_i[-23]
#   
#   # save
#   res[i] <- cor(as.numeric(nodes_i[1,]), g_vals$ppos, method = 'spearman')
#   i <- i+1
# }
# 
# d <- data.frame('perm' = res)
# 
# # Get empirical correlation
# empirical <- cor(g_vals$ppos, g_vals$receptor, method = 'spearman')
# 
# sum(d$perm >= empirical) / 1000
# 
# remove(receptor,c)
# #########
