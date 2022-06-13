library(readr)
library(readxl)
library(dplyr)
library(survival)
library(ggplot2)
cohort <- read_excel("COHORT_February2022_latest.xlsx")


cohort$IID <- paste(cohort$accession, cohort$SUBJECT_ID, 'normal',sep='_')
#cohort$IID <- paste(cohort$accession, cohort$SUBJECT_ID,sep='_')
#cohort = subset(cohort, accession != "phs000452")

PGS_1e4_RAPIDO_panIO <- read_table2("/Users/skleeman/Dropbox (CSHL Dropbox Team)/PhD/Data/ukbiobank/UKB380_PGS_LDPRED2_panIO_exome.sscore")
PGS_1e4_RAPIDO_panIO = inner_join(PGS_1e4_RAPIDO_panIO, cohort, by="IID")

PGS_1e4_RAPIDO_panIO$SCORE = PGS_1e4_RAPIDO_panIO$SCORE1_AVG

exome_gwas_PC_EUR <- read_table2("/Users/skleeman/Dropbox (CSHL Dropbox Team)/Scratch/Downloads_old/exome_gwas_PC_EUR.tsv")
names(exome_gwas_PC_EUR)[1]="IID"
PGS_1e4_RAPIDO_panIO = inner_join(PGS_1e4_RAPIDO_panIO, exome_gwas_PC_EUR, by="IID")

scale2 <- function(x, na.rm = TRUE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)

PGS_1e4_RAPIDO_panIO$score_scale = scale2(PGS_1e4_RAPIDO_panIO$SCORE)

PGS_1e4_RAPIDO_panIO$score_decile <- as.factor(ntile(PGS_1e4_RAPIDO_panIO$SCORE, 10))
PGS_1e4_RAPIDO_panIO$score_quartile <- as.factor(ntile(PGS_1e4_RAPIDO_panIO$SCORE, 4))


#Cancer sensitivity analysis

comb3 = data.frame()
for(i in 1:length(unique(PGS_1e4_RAPIDO_panIO$cancer))) {
  
  cancerx = unique(PGS_1e4_RAPIDO_panIO$cancer)[i]
  sub = subset(PGS_1e4_RAPIDO_panIO, cancer==cancerx)
  print(cancerx)
  res.cox <- coxph(Surv(os_days, os_stat) ~ score_scale + sex + PC1 + PC2 + PC3 + PC4, data = sub)
  print(summary(res.cox))
  mylogit <- glm(dcb ~ score_scale + sex + PC1 + PC2 + PC3 + PC4, data = sub, family = "binomial")
  print(summary(mylogit))
  
  add = data.frame(type="logistic", comparison=cancerx, HR = exp(mylogit$coefficients[2]), lower = exp(confint(mylogit, level=0.95)[2,1]),
                   upper = exp(confint(mylogit, level=0.95)[2,2]))
  comb3 = rbind(comb3, add)
  
}



comb3 %>% filter(comparison != "hnscc") %>% ggplot(
  aes(x = comparison,y = HR, ymin = lower, ymax = upper ))+
  geom_pointrange(aes(col=comparison))+
  geom_hline(aes(fill=comparison),yintercept =1, linetype=2)+
  xlab('Group')+ ylab("Hazard Ratio (95% Confidence Interval)")+
  geom_errorbar(aes(ymin=lower, ymax=upper,col=comparison),width=0.5,cex=1)+ 
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"))+
  coord_flip()+theme_bw()+xlab("")+ theme(legend.position = "none") 

ggsave('/Users/skleeman/Google Drive/PhD/Janowitz/cystatinc_share/figures/sub/figure5c.pdf', width=4, height=6)





#Overall survival
comb2= data.frame()
res.cox <- coxph(Surv(os_days, os_stat) ~ score_scale + sex + PC1 + PC2 + PC3 + PC4 + cancer, data = PGS_1e4_RAPIDO_panIO)

test = summary(res.cox)

add = data.frame(type="cox", comparison="overall_survival", HR = test$conf.int[1,1], lower = test$conf.int[1,3],
                 upper = test$conf.int[1,4])
comb2 = rbind(comb2, add)

#PFS
res.cox <- coxph(Surv(pfs_days, pfs_stat) ~ score_scale + sex + PC1 + PC2 + PC3 + PC4 + cancer, data = PGS_1e4_RAPIDO_panIO)

test = summary(res.cox)

add = data.frame(type="cox", comparison="pf_survival", HR = test$conf.int[1,1], lower = test$conf.int[1,3],
                 upper = test$conf.int[1,4])
comb2 = rbind(comb2, add)

#DCB
mylogit <- glm(dcb ~ score_scale + sex + PC1 + PC2 + PC3 + PC4 + cancer, data = PGS_1e4_RAPIDO_panIO, family = "binomial")

add = data.frame(type="logistic", comparison="dcb", HR = exp(mylogit$coefficients[2]), lower = exp(confint(mylogit, level=0.95)[2,1]),
                 upper = exp(confint(mylogit, level=0.95)[2,2]))
comb2 = rbind(comb2, add)

#Plot

comb2 %>% ggplot(
  aes(x = comparison,y = HR, ymin = lower, ymax = upper ))+
  geom_pointrange(aes(col=comparison))+
  geom_hline(aes(fill=comparison),yintercept =1, linetype=2)+
  xlab('Group')+ ylab("Hazard Ratio (95% Confidence Interval)")+
  geom_errorbar(aes(ymin=lower, ymax=upper,col=comparison),width=0.5,cex=1)+ 
  facet_wrap(~type,strip.position="left",nrow=9,scales = "free_y") +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"))+
  coord_flip()+theme_bw()+xlab("")+ theme(legend.position = "none") 

ggsave('/Users/skleeman/Google Drive/PhD/Janowitz/cystatinc_share/figures/sub/figure5b.pdf', width=8, height=5)



