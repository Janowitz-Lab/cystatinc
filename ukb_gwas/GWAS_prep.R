library(fst)
library(dplyr)
library(readr)
library(tidyr)
library(reshape2)
library(tidyverse)
library(nephro)

df <- fst("data/biobank_whole.fst")

#Import static data

qc = c("f22000_","f22007_","f22008_","f22001_","f22021_","f22006_","f22019_","f22027_","f22028_", "f31_0_0", "f50_0_0")
add_values <- grep(paste(qc,collapse="|"), names(df), value=TRUE)

fixed<-df[,names(df) %in% c("eid",add_values)]
fixed <- subset(fixed, is.na(use_in_phasing_chromosomes_122_f22028_0_0)==FALSE) #Only include samples where genomic data was used for phasing
fixed <- subset(fixed, is.na(sex_chromosome_aneuploidy_f22019_0_0) == TRUE) #Remove sex chromosome aneuploidy
fixed <- subset(fixed, is.na(outliers_for_heterozygosity_or_missing_rate_f22027_0_0) == TRUE) #Remove outliers for het/missing rate
fixed <- subset(fixed, sex_f31_0_0 == genetic_sex_f22001_0_0) #Remove discordant sex

names(fixed)[1]<-"id"

#Import relevant covariates

id_codes=c("date_of_attending_assessment_centre_f53", "weight_f21002","standing_height_f50","body_mass_index_bmi_f21001","creatinine_f30700","cystatin_c_f30720","age_when_attended_assessment_centre_f21003","uk_biobank_assessment_centre_f54_0_0")


my_data <- list()

for(i in 1:length(id_codes)) {
  sub<-df[,grepl(paste(c("eid",id_codes[i]), collapse="|"), names(df))]
  long <- melt(sub, id.vars="eid")
  long$variable <- sub('.*_', '', gsub('_0(.*)','\\1',long$variable))
  
  names(long)<-c("eid","visit",id_codes[i])
  my_data[[i]] <- long
}


collated<-my_data %>% reduce(full_join, by = c("eid","visit"))
collated<-collated[rowSums(!is.na(collated)) > 2,]
names(collated)[1]<-"id"

#Remove any NA values, take earliest available complete datapoint

collated<-collated[complete.cases(collated), ]
collated<-collated[!duplicated(collated$id),]

#Merge

collated<-inner_join(fixed, collated, by="id")

#Derive relevant phenotypes

collated <- collated %>% mutate(sex_f31_0_0 = ifelse(sex_f31_0_0=="Male",1,0)) %>% mutate(egfr_creatinine=CKDEpi.creat((creatinine_f30700*0.0113),sex_f31_0_0, age_when_attended_assessment_centre_f21003, replicate(nrow(collated), 0)), egfr_cystatin = CKDEpi.cys((cystatin_c_f30720),sex_f31_0_0, age_when_attended_assessment_centre_f21003))

#Save file
collated$genotype_measurement_batch_f22000_0_0<-ifelse(collated$genotype_measurement_batch_f22000_0_0 >= 0, "A","B")
write.csv(collated, "phenotype_input_cystatinc_gwas.csv")
