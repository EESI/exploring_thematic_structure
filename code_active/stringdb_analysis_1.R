library(stm)
library(biom)
library(readr)
library(tidyr)
library(dplyr)
library(fastICA)
library(randomForest)
library(stringr)
library(kernlab)
library(Rcpp)
library(parallel)
library(foreach)
library(ape)
library(phyloseq)
library(doParallel)
library(stm)
library(LDAvis)
library(glmnet)
library(ggplot2)
library(vegan)


# source('~/Dropbox/stm_microbiome/code_active/nav_froz_fxns_3.R') 
# R <- create_ref_set() # no omission
# saveRDS(R,'~/Dropbox/stm_microbiome/data_active/string_ref_set.rds')
# 
# R <- readRDS('~/Dropbox/stm_microbiome/data_active/string_ref_set.rds')
# dict <- create_dict(R)
# saveRDS(dict,'~/Dropbox/stm_microbiome/data_active/string_dict.rds')




params <- expand.grid(K=c(25,50,75,125,200),
            content=c(FALSE,TRUE),
            variable=c('DIAGNOSIS','ISOLATION_SOURCE'))
params <- params %>% 
   filter(!(content == FALSE & variable == 'ISOLATION_SOURCE')) %>%
   arrange(K,variable)

for (param in 1:NROW(params)){
#START LOOP FOR PARAMETERS

   

   
   
source('~/Dropbox/stm_microbiome/code_active/stm_functions.R')
source('~/Dropbox/stm_microbiome/code_active/nav_froz_fxns_3.R') 
source('~/Dropbox/stm_microbiome/code_active/performance_1.R')
source('~/Dropbox/stm_microbiome/code_active/framework.R')

# R <- readRDS('~/Dropbox/stm_microbiome/data_active/string_ref_set.rds')
dict <- readRDS('~/Dropbox/stm_microbiome/data_active/string_dict.rds')



random_seed <- FALSE
load_fit <- TRUE
save_fit <- FALSE
K <- params[param,]$K
cn_normalize <- TRUE
content <- params[param,]$content
seq_sim <- 's97'
variable <- as.character(params[param,]$variable)



prepare_framework(random_seed,K,cn_normalize,content,variable,seq_sim)
output_dir <- file.path('~/Dropbox/stm_microbiome/output_active/',save_fit_foldername)

biom_file <- read_biom(file.path(save_dir,save_fit_foldername,'beta_predicted_cogs.biom'))
beta_cogs <- prep_string_data(biom_file,ref_set=R)
biom_file <- read_biom(file.path(save_dir,save_fit_foldername,'beta_permuted_predicted_cogs.biom'))
beta_cogs_permuted <- prep_string_data(biom_file)

matches_mat <- matrix(NA,K,3,dimnames=list(1:K,c('Total','Direct','Indirect')))
matches_mat_permuted <- matrix(NA,K,3,dimnames=list(1:K,c('Total','Direct','Indirect')))
cutoff <- .001
cat('\n')
for (k in 1:K){
   matches <- find_matching_pairs(dict,beta_cogs$cogs_ra,k,min_cutoff=cutoff)
   matches_permuted <- find_matching_pairs(dict,beta_cogs_permuted$cogs_ra,k,min_cutoff=cutoff)
   
   matches_mat[k,] <- matches
   matches_mat_permuted[k,] <- matches_permuted
}

sink(file.path(output_dir,'string_direct_sig.txt'))
cat('Direct Connections\n\n')
matches_hits <- data.frame(True=matches_mat[,2],Permuted=matches_mat_permuted[,2])
if (shapiro.test(matches_hits$True)$p.value < .05 & shapiro.test(matches_hits$Permuted)$p.value < .05){
   cat('Log transformted data.\n\n')
   print(t.test(log(matches_hits$True),log(matches_hits$Permuted),paired=TRUE,alternative='greater'))
}else{
   print(t.test(matches_hits$True,matches_hits$Permuted,paired=TRUE,alternative='greater'))
}
print(wilcox.test(matches_hits$True,matches_hits$Permuted,paired=TRUE,exact=FALSE,alternative='greater'))
sink()

pdf(file.path(output_dir,'string_direct_sig.pdf'), height=8, width=12)
matches_hits <- data.frame(True=matches_mat[,2],Permuted=matches_mat_permuted[,2])
if (shapiro.test(matches_hits$True)$p.value < .05 & shapiro.test(matches_hits$Permuted)$p.value < .05){
   print(log(matches_hits) %>% gather(set,count) %>% ggplot(aes(set,count)) + geom_boxplot() + ggtitle('Direct Connections (Log transformed)'))
}else{
   print(matches_hits %>% gather(set,count) %>% ggplot(aes(set,count)) + geom_boxplot() + ggtitle('Direct Connections'))
}
dev.off()


sink(file.path(output_dir,'string_coverage_sig.txt'))
cat('Direct & Indirect Connections\n\n')
matches_hits <- data.frame(True=matches_mat[,1],Permuted=matches_mat_permuted[,1])
if (shapiro.test(matches_hits$True)$p.value < .05 & shapiro.test(matches_hits$Permuted)$p.value < .05){
   cat('Log transformted data.\n\n')
   print(t.test(log(matches_hits$True),log(matches_hits$Permuted),paired=TRUE,alternative='greater'))
}else{
   print(t.test(matches_hits$True,matches_hits$Permuted,paired=TRUE,alternative='greater'))
}
print(wilcox.test(matches_hits$True,matches_hits$Permuted,paired=TRUE,exact=FALSE,alternative='greater'))
sink()

pdf(file.path(output_dir,'string_coverage_sig.pdf'), height=8, width=12)
cat('Direct & Indirect Connections')
matches_hits <- data.frame(True=matches_mat[,1],Permuted=matches_mat_permuted[,1])
if (shapiro.test(matches_hits$True)$p.value < .05 & shapiro.test(matches_hits$Permuted)$p.value < .05){
   print(log(matches_hits) %>% gather(set,count) %>% ggplot(aes(set,count)) + geom_boxplot() + ggtitle('Direct & Indirect Connections (Log transformed)'))
}else{
   print(matches_hits %>% gather(set,count) %>% ggplot(aes(set,count)) + geom_boxplot() + ggtitle('Direct & Indirect Connections'))
}
dev.off()



load_fits(file.path(save_dir,save_fit_foldername,save_fit_filename))
train_fit <- fit_frozen[[1]]
beta_frozen <- t(exp(train_fit$beta$logbeta[[1]]))
beta_frozen_ra <- beta_frozen


beta <- as.matrix(biom_data(read_biom(file.path(save_dir,save_fit_foldername,'beta_table.biom'))))
beta_permuted <- as.matrix(biom_data(read_biom(file.path(save_dir,save_fit_foldername,'beta_table_permuted.biom'))))

beta <- t(t(beta)/colSums(beta))
beta_permuted <- t(t(beta_permuted)/colSums(beta_permuted))

df <- data.frame(True=matches_mat[,2],Permuted=matches_mat_permuted[,2]) %>% gather(Set,Direct)
df$Coverage <- data.frame(True=matches_mat[,1],Permuted=matches_mat_permuted[,1]) %>% gather(Set,Coverage) %>% dplyr::select(Coverage)
df$UniqueTaxa <- c(apply(beta,2,function(x) sum(x > 10e-5)),apply(beta_permuted,2,function(x) sum(x > 10e-5)))
df$ShannonDiversity <- c(diversity(beta,index='shannon',MARGIN=2),diversity(beta_permuted,index='shannon',MARGIN=2))

df <- df %>% mutate(Coverage=scale(Coverage),
                    Direct=scale(Direct),
                    UniqueTaxa=scale(UniqueTaxa),
                    Set=ifelse(Set == 'True',1,0))


sink(file.path(output_dir,'string_glm_sig.txt'))
print(summary(glm1 <- glm(Set ~ Direct + UniqueTaxa, family=binomial(link='logit'),data=df)))
print(summary(glm2 <- glm(Set ~ Coverage + UniqueTaxa, family=binomial(link='logit'),data=df)))

print(summary(glm1 <- glm(Set ~ Direct + ShannonDiversity, family=binomial(link='logit'),data=df)))
print(summary(glm2 <- glm(Set ~ Coverage + ShannonDiversity, family=binomial(link='logit'),data=df)))


print(summary(lm1 <- lm(Direct ~ UniqueTaxa,filter(df,Set==0))))
print(summary(lm2 <- lm(Coverage ~ UniqueTaxa,filter(df,Set==1))))

print(summary(lm1 <- lm(Direct ~ ShannonDiversity,filter(df,Set==0))))
print(summary(lm2 <- lm(Coverage ~ ShannonDiversity,filter(df,Set==1))))
sink()

pdf(file.path(output_dir,'string_glm_sig.pdf'), height=8, width=12)
print(
df %>% 
   gather(Network,Count,-Set,-UniqueTaxa,-ShannonDiversity) %>%
   mutate(Set=ifelse(Set == 1,'True','Permuted')) %>%
   ggplot(aes(x=UniqueTaxa,y=Count,colour=Network)) + geom_point(alpha=.5) + facet_wrap(~Set) + geom_smooth(method=lm,se=FALSE)
)
print(
   df %>% 
      gather(Network,Count,-Set,-UniqueTaxa,-ShannonDiversity) %>%
      mutate(Set=ifelse(Set == 1,'True','Permuted')) %>%
      ggplot(aes(x=ShannonDiversity,y=Count,colour=Network)) + geom_point(alpha=.5) + facet_wrap(~Set) + geom_smooth(method=lm,se=FALSE)
)
dev.off()





rm(list = ls()[!(ls() %in% c('param','params'))]) ####
} ## end save loop ####

