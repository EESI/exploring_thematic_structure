meta$FEMALE <- ifelse(meta$SEX == 'female',1,0)
meta$AGE <- as.numeric(meta$AGE_YEARS)
colnames(taxa) <- c('Kingdom','Phylum','Class','Order','Family','Genus','Species')
PS <- phyloseq(otu_table(otu,taxa_are_rows=FALSE),
               sample_data(meta),
               tax_table(taxa))
PSLOG <- transform_sample_counts(PS, function(x) log10(x+1))
CCA <- ordinate(PSLOG,'CCA',formula=PSLOG ~ FEMALE + AGE)


tax <- tax_table(PS)@.Data %>% data.frame(stringsAsFactors=FALSE)
tax$otu_id <- seq_len(ncol(otu_table(PS)))

ps_scores <- vegan::scores(CCA)
sites <- data.frame(ps_scores$sites)
sites$SampleID <- rownames(sites)
sites <- sites %>% left_join(sample_data(PS))

species <- data.frame(ps_scores$species)
species$otu_id <- seq_along(colnames(otu_table(PS)))
species <- species %>% left_join(tax)



evals_prop <- 100 * CCA$CCA$eig[1:2] / sum(CCA$CA$eig)

colourCount <- length(unique(species$Phylum)) + 2
getPalette <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))

ggplot() +
  geom_point(data = sites, aes(x = CCA1, y = CCA2, color=SEX), shape = 2, alpha = 0.8) +
  geom_point(data = species, aes(x = CCA1, y = CCA2, color = Phylum), size = 0.5) +
  geom_text_repel(data = species %>% filter(CCA2 < -2),
                  aes(x = CCA1, y = CCA2, label = otu_id),
                  size = 1.5, segment.size = 0.1) +
  facet_wrap( ~ AGE_DECADE_RANK) +
  guides(col = guide_legend(override.aes = list(size = 3))) +
  labs(x = sprintf("Axis1 [%s%% variance]", round(evals_prop[1], 2)),
       y = sprintf("Axis2 [%s%% variance]", round(evals_prop[2], 2))) +
  scale_color_manual(values=getPalette(colourCount)) + 
  coord_fixed(sqrt(CCA$CCA$eig[2] / CCA$CCA$eig[1])) +
  theme_bw() +
  theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)))















meta$FEMALE <- ifelse(meta$SEX == 'female',1,0)
meta$AGE <- as.numeric(meta$AGE_YEARS)
meta$MILESTONE <- ifelse(meta$AGE <= 10, 'Prepubescent',
                         ifelse(meta$AGE > 10 & meta$AGE <= 17, 'Pubescent',
                                ifelse(meta$AGE > 17 & meta$AGE <= 54, 'Adult',
                                       ifelse(meta$AGE > 54 & meta$AGE <= 69, 'Postmenopausal',
                                              'Elderly'))))
meta$MILESTONE <- factor(meta$MILESTONE,levels=c('Prepubescent','Pubescent','Adult','Postmenopausal','Elderly'))
colnames(taxa) <- c('Kingdom','Phylum','Class','Order','Family','Genus','Species')
PS <- phyloseq(otu_table(otu,taxa_are_rows=FALSE),
               sample_data(meta),
               tax_table(taxa))
PSLOG <- transform_sample_counts(PS, function(x) log10(x+1))
CCA <- ordinate(PSLOG,'CCA',formula=PSLOG ~ FEMALE + AGE)


tax <- tax_table(PS)@.Data %>% data.frame(stringsAsFactors=FALSE)
tax$otu_id <- seq_len(ncol(otu_table(PS)))

ps_scores <- vegan::scores(CCA)
sites <- data.frame(ps_scores$sites)
sites$SampleID <- rownames(sites)
sites <- sites %>% left_join(sample_data(PS))

species <- data.frame(ps_scores$species)
species$otu_id <- seq_along(colnames(otu_table(PS)))
species <- species %>% left_join(tax)



evals_prop <- 100 * CCA$CCA$eig[1:2] / sum(CCA$CA$eig)

colourCount <- length(unique(species$Phylum)) + 2
getPalette <- colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))

ggplot() +
  geom_point(data = sites, aes(x = CCA1, y = CCA2), shape = 2, alpha = 0.8) +
  geom_point(data = species, aes(x = CCA1, y = CCA2, color = Phylum), size = 0.5) +
#   geom_text_repel(data = species %>% filter(CCA1 > 1),
#                   aes(x = CCA1, y = CCA2, label = otu_id),
#                   size = 4.5, segment.size = 2.0) +
  facet_grid(SEX ~ MILESTONE) +
  guides(col = guide_legend(override.aes = list(size = 3))) +
  labs(x = sprintf("Axis1 [%s%% variance]", round(evals_prop[1], 2)),
       y = sprintf("Axis2 [%s%% variance]", round(evals_prop[2], 2))) +
  scale_color_manual(values=getPalette(colourCount)) + 
  coord_fixed(sqrt(CCA$CCA$eig[2] / CCA$CCA$eig[1])) +
  theme_bw() +
  theme(panel.border = element_rect(color = "#787878", fill = alpha("white", 0)))


