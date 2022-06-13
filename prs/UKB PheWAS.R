library(readxl)
library(reshape2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(stringr)
library(survival)
library(fst)
library(tidyr)
library(ggplot2)
library(anytime)
library(lubridate)

##Phenome-wide analysis (Figure 2c)

#Import files
df <- fst("data/biobank_update_new.fst")
map_tte <- read_excel("~/Downloads/map_tte.xlsx")
icd_coding <- read_delim("data/icd_coding.tsv",
                         "\t", escape_double = FALSE, trim_ws = TRUE)
covariates <- read_table2("~/Google Drive/PhD/Data/ukbiobank/covariates_eur.tsv")
Phecode_map_v1_2_icd10_beta <- read_csv("data/Phecode_map_v1_2_icd10_beta.csv")
phecode_definitions1_2 <- read_csv("data/phecode_definitions1.2.csv")
phecode_definitions1_3 <- read_excel("phecode_cancer.xlsx") #My filtered set

#Load CST3 polygenic risk score
scale2 <- function(x, na.rm = TRUE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)
PGS_1e4_RAPIDO <- read_table2("~/Google Drive/PhD/Data/ukbiobank/UKB380_PGS_LDPRED2_UKB_exome.sscore")
names(PGS_1e4_RAPIDO)[2]<-"id"
PGS_1e4_RAPIDO$SCORE = scale2(PGS_1e4_RAPIDO$SCORE1_AVG)

populations <- read_csv("~/Google Drive/PhD/Data/ukbiobank/populations.csv")
PGS_1e4_RAPIDO = subset(PGS_1e4_RAPIDO, id %in% subset(populations, pop=="EUR")$s)

validation_cohort <- read_delim("~/Google Drive/PhD/Data/ukbiobank/validation_cohort.tsv", 
                                "\t", escape_double = FALSE, trim_ws = TRUE)

PGS_1e4_RAPIDO = subset(PGS_1e4_RAPIDO, id %in% validation_cohort$s)
names(PGS_1e4_RAPIDO) = c("fid","s","NMISS_ALLELE_CT","NAMED_ALLELE_DOSAGE_SUM","SCORE1_AVG","call")

#Load Biobank codes
map_tte=subset(map_tte, type=="Date")
map_tte$id_code = paste('f.',map_tte$id_code,sep='')

list = str_split(map_tte$description, " ")
map_tte$coding = sapply(list, "[[", 2)

#Load phecodes

icd_coding$coding<-gsub( " .*$", "", icd_coding$meaning)
names(phecode_definitions1_2)[1]<-"PHECODE"
Phecode_map_v1_2_icd10_beta<-left_join(Phecode_map_v1_2_icd10_beta, phecode_definitions1_2, by="PHECODE")

names(Phecode_map_v1_2_icd10_beta)[1]<-"coding"
phecode<-left_join(icd_coding, Phecode_map_v1_2_icd10_beta, by="coding")

map_tte1 = inner_join(map_tte, phecode, by="coding")
map_tte1 = subset(map_tte1, is.na(PHECODE)==FALSE)
map_tte$coding = paste(map_tte$coding,'.0',sep='')
map_tte2 = inner_join(map_tte, phecode, by="coding")
map_tte2 = subset(map_tte2, is.na(PHECODE)==FALSE)

map_tte = rbind(map_tte1,map_tte2)
map_tte = map_tte[!duplicated(map_tte$id_code),]

#Load phenotypes


id_codes=c('f.34\\.',map_tte$id_code)

sub<-df[,grepl(paste(c("f.eid",id_codes), collapse="|"), names(df))]

#End date per patient
id_codes=c('f.21003\\.','f.40007\\.')
subx<-df[,grepl(paste(c("f.eid",id_codes), collapse="|"), names(df))]

max_age = data.frame(max_age = do.call(pmax, c(subx[,-1], na.rm=TRUE)))

#Add cancer codes
id_codes=c('f.40005\\.','f.40006\\.')
my_data <- list()
for(i in 1:length(id_codes)) {
  subx<-df[,grepl(paste(c("f.eid",id_codes[i]), collapse="|"), names(df))]
  long <- melt(subx, id.vars="f.eid")
  list = str_split(long$variable,"[.]")
  long$variable = sapply(list, "[[", 3)
  names(long)<-c("eid","visit",id_codes[i])
  my_data[[i]] <- long
}
collated<-my_data %>% purrr:reduce(full_join, by = c("eid","visit"))
collated<-collated[rowSums(!is.na(collated)) > 3,]
names(collated) = c("eid","visit","date_of_diagnosis","icd_code")

##Load phecodes again

phecode$coding = gsub('\\.', '', phecode$coding)

names(phecode_definitions1_3)[1]<-"PHECODE"
phecode_definitions1_3<-subset(phecode_definitions1_3, include=="yes")
phecode_definitions1_3 = phecode_definitions1_3 %>% select('PHECODE','type')
phecode_definitions1_3$PHECODE<-as.character(phecode_definitions1_3$PHECODE)

phecode = inner_join(phecode, phecode_definitions1_3, by="PHECODE")

m<-match(collated$icd_code, phecode$coding)
collated$phecode = phecode$type[m]
collated = subset(collated, is.na(phecode)==FALSE)
collated$phecode = paste(collated$phecode, 'cancer',sep=' ')

collated = collated %>%
  group_by(eid) %>%
  arrange(date_of_diagnosis) %>%
  slice(1L) %>% select('eid','date_of_diagnosis','phecode')

cancer = pivot_wider(collated, names_from=phecode, values_from=date_of_diagnosis)
names(cancer)[1]='f.eid'
cancer$anycancer = do.call(pmax, c(cancer[,c(2:(ncol(cancer)))], na.rm=TRUE))

#Merge to phecodes
unique_phes = unique(map_tte$PHECODE)
sub2 = data.frame(f.eid=sub[,1],f.34.0.0=sub[,2])
for(i in 1:length(unique_phes)) {
  selectx = subset(map_tte, PHECODE == unique_phes[i])$id_code
  extractx = sub %>% select(selectx)
  addx = data.frame(column = do.call(pmin, c(extractx, na.rm=TRUE) ))
  names(addx) = subset(map_tte, PHECODE == unique_phes[i])$phenotype[1]
  sub2 = cbind(sub2, addx)
}

sub = sub2
sub = cbind(max_age, sub)
sub = left_join(sub, cancer, by="f.eid") #Add cancer fields


#Import covariates

names(covariates)[2]<-"s"

input = inner_join(PGS_1e4_RAPIDO, covariates, by="s")


#Prepare input

sub <- subset(sub, f.eid %in% input$s)
m<-match(sub$f.eid, input$s)
sub$call <- input$call[m]
sub$sex <- as.factor(input$sex[m])
sub$PC1 <- input$PC1[m]
sub$PC2 <- input$PC2[m]
sub$PC3 <- input$PC3[m]
sub$PC4 <- input$PC4[m]
sub$PC5 <- input$PC5[m]
sub$PC6 <- input$PC6[m]
sub$PC7 <- input$PC7[m]
sub$PC8 <- input$PC8[m]
sub$PC9 <- input$PC9[m]
sub$PC10 <- input$PC10[m]
sub$array <- as.factor(input$array[m])
sub$centre <- as.factor(input$centre[m])

#sub$call = as.factor(ifelse(sub$call>0,1,0))

#Max date across ICD codes we are considering
sub[ , c(4:(ncol(sub)-14))] <- lapply(sub[ ,c(4:(ncol(sub)-14))], anydate) #Sort out date format
max_age2=do.call(pmax, c(sub[,c(4:(ncol(sub)-14))], na.rm=TRUE))
sub$max_age = apply(cbind((decimal_date(max_age2) - sub$f.34.0.0), sub$max_age), 1, max, na.rm=TRUE) #Either (most recent ICD-DOB) or (most recent attendance/death-DOB)
sub$max_age = ifelse(sub$max_age > (2021-sub$f.34.0.0), (2021-sub$f.34.0.0), sub$max_age) #Set maximum possible age to 2021-DOB

framex = data.frame()

sex_map = data.frame(id=names(sub[,4:(ncol(sub)-14)]))
m<-match(sex_map$id, Phecode_map_v1_2_icd10_beta$phenotype)
sex_map$gender = Phecode_map_v1_2_icd10_beta$sex[m]
list_female = c("endometrial cancer","cervical cancer","ovary cancer")
list_male = c("prostate cancer")
sex_map$gender[sex_map$id %in% list_female]="Female"
sex_map$gender[sex_map$id %in% list_male]="Male"

pb = txtProgressBar(min = 4, (ncol(sub)-14), initial = 0) 

#for(i in 4:(ncol(sub)-14)) {
for(i in list) {
  
  col=sub[,i]
  name = names(sub)[i]
  stat = ifelse(is.na(col),0,1)
  max_age = sub$max_age
  time = sub$max_age
  time[stat==1] = decimal_date(col[stat==1]) - sub$f.34.0.0[stat==1]
  time = ifelse(time > max_age, max_age, time)
  time = ifelse(time<0, 0, time)
  
  use = data.frame(sex = sub$sex, PC1 = sub$PC1, PC2 = sub$PC2, PC3=sub$PC3, PC4=sub$PC4, PC5 = sub$PC5,
                   PC6 = sub$PC6, PC7 = sub$PC7, PC8=sub$PC8, PC9 = sub$PC9, PC10 = sub$PC10, array=sub$array,
                   centre=sub$centre, time=time, stat=stat, call=sub$call, year_of_birth=sub$f.34.0.0, max_age = max_age)
  
  sex_match = subset(sex_map, id==name)$gender
  sex_match = ifelse(is.na(sex_match)==TRUE,"Both",sex_match)
  
  if(sex_match == "Female" | sex_match=="Male") {
    res.cox <- coxph(Surv(time, stat) ~ call + year_of_birth + PC1 + PC2 + PC3 + PC4,data=use)
    test = summary(res.cox)
    model <- glm(stat ~ call + year_of_birth + PC1 + PC2 + PC3 + PC4,family=binomial(link='logit'),data=use)
    test2 = summary(model)
  } else {
    res.cox <- coxph(Surv(time, stat) ~ call + year_of_birth + sex + PC1 + PC2 + PC3 + PC4,data=use)
    test = summary(res.cox)
    model <- glm(stat ~ call + year_of_birth + sex + PC1 + PC2 + PC3 + PC4,family=binomial(link='logit'),data=use)
    test2 = summary(model)
  }
  

  
  p_value = as.numeric(test$coefficients[1,5])
  p_value2 = as.numeric(test2$coefficients[2,4])
  add = data.frame(name=name, events=sum(stat), p_value = p_value, hr=test$coefficients[1,2], p_value_logit = p_value2,
                   beta_logit = test2$coefficients[2,1])
  #add = data.frame(name=name, events=sum(stat), p_value = p_value, hr=test$coefficients[1,2])
  framex=rbind(framex,add)
  setTxtProgressBar(pb,i)
}

m<-match(framex$name, Phecode_map_v1_2_icd10_beta$phenotype)
framex$phecode = Phecode_map_v1_2_icd10_beta$PHECODE[m]
matched = subset(framex, is.na(phecode)==FALSE)
unmatched = subset(framex, is.na(phecode)==TRUE)


phecode_definitions1_3$compare = paste(phecode_definitions1_3$type, 'cancer',sep=' ')
m<-match(unmatched$name, phecode_definitions1_3$compare)
unmatched$phecode = phecode_definitions1_3$PHECODE[m]

together = rbind(matched, unmatched)

library(PheWAS)
input = together %>% select('phecode', 'hr', 'p_value')
names(input) = c("phenotype","OR",'p')
#input$OR = exp(input$OR)

phewasManhattan(input, OR.direction	= TRUE)


##Cancer specific analysis (Figure 2d)
merge_sub_2000_survivalx = merge_sub_2000_survival

merge_sub_2000_survivalx = subset(merge_sub_2000_survivalx, primary != "cervical")
merge_sub_2000_survivalx = subset(merge_sub_2000_survivalx, primary != "prostate")
names(merge_sub_2000_survivalx)[1]='s'
input = inner_join(PGS_1e4_RAPIDO, covariates, by="s")

merge_sub_2000_survivalx = inner_join(merge_sub_2000_survivalx, input, by='s')

ages = data.frame(id = df$eid, age = df$age_at_recruitment_f21022_0_0)

m<-match(merge_sub_2000_survivalx$s, ages$id)
merge_sub_2000_survivalx$age_recruitment = ages$age[m]

res.cox <- coxph(Surv(days_to_death, status_os) ~ call + sex.x + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + operation + primary + age_of_diagnosis + array + centre,data=merge_sub_2000_survivalx)
summary(res.cox)

res.cox <- coxph(Surv(days_to_death, status_os) ~ call + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + primary + age_of_diagnosis + array + centre,data=subset(merge_sub_2000_survivalx, operation=="no"))
summary(res.cox)



framex = data.frame()

prs_survival = merge_sub_2000_survivalx

for(i in 1:length(unique(prs_survival$primary))) {
  tryCatch({
    cancer = unique(prs_survival$primary)[i]
    use = subset(prs_survival, primary==cancer)
    if(length(unique(use$sex.x)) ==2 ) {
      res.cox = coxph(Surv(days_to_death, status_dss) ~ call + sex.x + PC1 + PC2 + PC3 + PC4 + operation + age_recruitment + array + centre, data=use)
    } else {
      res.cox = coxph(Surv(days_to_death, status_dss) ~ call + PC1 + PC2 + PC3 + PC4 +  operation + age_recruitment + array + centre, data=use)
    }
    test = summary(res.cox)
    p_value = as.numeric(test$coefficients[1,5])
    add = data.frame(name=cancer, p_value = p_value, hr=test$coefficients[1,2], se=test$coefficients[1,3], lower = test$conf.int[1,3],
                     upper = test$conf.int[1,4])
    test = summary(res.cox)
    framex=rbind(framex,add)
  }, error=function(e){})
}


library(meta)
meta = metagen(log(framex$hr), framex$se, sm = "HR", studlab=framex$name)

pdf(file = "/Users/skleeman/Google Drive/PhD/Janowitz/cystatinc_share/figures/sub/figure4cx.pdf", width = 8, height = 6)
forest.meta(meta,leftlabs = c("Cancer primary"), leftcols = c("studlab"),rightcols=c("effect", "ci"),sortvar = TE)

dev.off()

library(ggplot2)


add = data.frame(name='total', p_value = 0.0586, hr = 1.0950, se=0, lower =0.9967,upper= 1.2030)

#framex=rbind(framex,add)

framex  %>% ggplot(
  aes(x = name,y = hr, ymin = lower, ymax = upper))+
  geom_pointrange(aes(col=name))+
  geom_hline(aes(fill=name),yintercept =1, linetype=2)+
  xlab('Group')+ ylab("Hazard Ratio (95% Confidence Interval)")+
  geom_errorbar(aes(ymin=lower, ymax=upper,col=name),width=0.5,cex=1)+ 
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"))+
  coord_flip()+theme_bw()+xlab("")+ theme(legend.position = "none") 