TAXA <- dat$taxa
colnames(TAXA) <- c('Kingdom','Phylum','Class','Order','Family','Genus','Species')

colnames(beta2) <- paste0('T',1:K)
tt <- order(coef_k,decreasing=TRUE)[c(1:5,31:35)]
B1 <- beta1[,tt]
B2 <- beta2[,tt]
B <- as.data.frame(cbind(B1,B2))



TAXON <- apply(TAXA[rownames(B),1:2], 1, function(x) paste0(x,collapse='.'))
TAXONP <- pretty_taxa_names(TAXON,'\\.',rownames(B))
TAXONPU <- unique(data.frame(TAXONP=TAXONP))
TAXONPU$OTU <- rownames(TAXONPU)

TAX <- 6

BB <- as.data.frame(t(B)) %>%
  mutate(Topic=rownames(.),
         Sex=rep(c('female','male'),each=NCOL(B)/2)) %>%
  gather(OTU,Abundance,-Topic,-Sex) %>%
  left_join(data.frame(TAXONP,OTU=names(TAXONP)),by='OTU') %>%
  group_by(Topic,Sex,TAXONP) %>%
  summarise(Abundance=sum(Abundance)) %>%
  ungroup() %>%
  filter(Abundance > 0) %>%
  group_by(Topic,Sex) %>%
  mutate(Frequency = Abundance/sum(Abundance)) %>%
  ungroup() %>%
  left_join(TAXONPU,by='TAXONP') %>%
  left_join(data.frame(TAXA[,1:TAX],OTU=rownames(TAXA)),by='OTU')

BB %>%
  mutate(Topic=factor(Topic,levels=paste0('T',tt))) %>%
  arrange(Phylum) %>%
  ggplot(aes(x=Sex,y=Frequency,fill=Phylum)) + 
  geom_bar(colour='black',stat='identity') +
  facet_grid(.~Topic)


BB %>%
  mutate(Topic=factor(Topic,levels=paste0('T',tt))) %>%
  group_by(Phylum) %>%
  mutate(Group=Abundance > sort(Abundance,decreasing=TRUE)[10]) %>%
  ungroup() %>%
  mutate(Phylum = ifelse(Group,as.character(Phylum),'Other')) %>%
  arrange(Phylum) %>%
  ggplot(aes(x=Sex,y=Frequency,fill=Phylum)) + 
  geom_bar(colour='black',stat='identity') +
  facet_grid(.~Topic)







TAXON <- apply(TAXA[rownames(B),1:3], 1, function(x) paste0(x,collapse='.'))
TAXONP <- pretty_taxa_names(TAXON,'\\.',rownames(B))
TAXONPU <- unique(data.frame(TAXONP=TAXONP))
TAXONPU$OTU <- rownames(TAXONPU)

BB <- as.data.frame(t(B)) %>%
  mutate(Topic=rownames(.),
         Sex=rep(c('female','male'),each=NCOL(B)/2)) %>%
  gather(OTU,Abundance,-Topic,-Sex) %>%
  left_join(data.frame(TAXONP,OTU=names(TAXONP)),by='OTU') %>%
  group_by(Topic,Sex,TAXONP) %>%
  summarise(Abundance=sum(Abundance)) %>%
  ungroup() %>%
  filter(Abundance > 0) %>%
  group_by(Topic,Sex) %>%
  mutate(Frequency = Abundance/sum(Abundance)) %>%
  ungroup() %>%
  left_join(TAXONPU,by='TAXONP') %>%
  left_join(data.frame(TAXA[,1:TAX],OTU=rownames(TAXA)),by='OTU')

BB %>%
  mutate(Topic=factor(Topic,levels=paste0('T',tt))) %>%
  group_by(Phylum) %>%
  mutate(GroupSum=sum(Abundance)) %>%
  ungroup() %>%
  filter(dense_rank(desc(GroupSum)) <= 5) %>%
  ggplot(aes(x=Sex,y=Frequency,fill=Class)) + 
  geom_bar(colour='black',stat='identity') +
  facet_grid(Phylum ~ Topic)




TAX <- 6
TAXA <- dat$taxa
colnames(TAXA) <- c('Kingdom','Phylum','Class','Order','Family','Genus','Species')

tt <- order(coef_k,decreasing=TRUE)[c(1:3)]
B1 <- beta1[,tt]
B2 <- beta2[,tt]
B <- as.data.frame(cbind(B1,B2))

TAXON <- apply(TAXA[rownames(B),1:5], 1, function(x) paste0(x,collapse='.'))
TAXONP <- pretty_taxa_names(TAXON,'\\.',rownames(B))
TAXONPU <- unique(data.frame(TAXONP=TAXONP))
TAXONPU$OTU <- rownames(TAXONPU)

BB <- as.data.frame(t(B)) %>%
  mutate(Topic=rownames(.),
         Sex=rep(c('female','male'),each=NCOL(B)/2)) %>%
  gather(OTU,Abundance,-Topic,-Sex) %>%
  left_join(data.frame(TAXONP,OTU=names(TAXONP)),by='OTU') %>%
  group_by(Topic,Sex,TAXONP) %>%
  summarise(Abundance=sum(Abundance)) %>%
  ungroup() %>%
  filter(Abundance > 0) %>%
  group_by(Topic,Sex) %>%
  mutate(Frequency = Abundance/sum(Abundance)) %>%
  ungroup() %>%
  left_join(TAXONPU,by='TAXONP') %>%
  left_join(data.frame(TAXA[,1:TAX],OTU=rownames(TAXA)),by='OTU')

BB %>%
  mutate(Topic=factor(Topic,levels=paste0('T',tt))) %>%
  group_by(Order) %>%
  mutate(GroupSum=sum(Abundance)) %>%
  ungroup() %>%
  filter(dense_rank(desc(GroupSum)) <= 5) %>%
  ggplot(aes(x=Sex,y=Frequency,fill=Family)) + 
  geom_bar(colour='black',stat='identity') +
  facet_grid(Order ~ Topic)


BBF <- BB %>%
  mutate(Topic=factor(Topic,levels=paste0('T',tt))) %>%
  group_by(Order) %>%
  mutate(GroupSum=sum(Abundance)) %>%
  ungroup() %>%
  filter(dense_rank(desc(GroupSum)) <= 5) %>%
  group_by(Family) %>%
  mutate(GroupSum=sum(Abundance)) %>%
  ungroup() %>%
  filter(dense_rank(desc(GroupSum)) <= 10) %>%
  ggplot(aes(x=Sex,y=Frequency,fill=Family)) + 
  geom_bar(colour='black',stat='identity') +
  facet_grid(Order ~ Topic)

