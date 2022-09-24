library(reshape2)
library(data.table)
library(ggplot2)
library(lme4)
# library("ggpubr")
library(emmeans)
library(lmerTest)
library(stringi)
library(stringr)
library(dplyr)
library(purrr)
library(tidyverse)

autist_path<- "/Volumes/My Passport for Mac/theta_4_7/trained/autism/dataframe_for_LMEM"
norm_path<- "/Volumes/My Passport for Mac/theta_4_7/trained/normotypical/df_for_lmem"


read_plus <- function(flnm) {
  read_csv(flnm) %>% 
    mutate(filename = flnm)
}

aut_df <-list.files(autist_path,pattern = "*.csv", full.names = T) %>% 
  map_df(~read_plus(.))
setDT(aut_df)
aut_df$"group"<- "autists"
aut_df$criterion_learning<- NULL
aut_df$...1<-NULL
norm_df <-list.files(norm_path,pattern = "*.csv", full.names = T) %>% 
  map_df(~read_plus(.))
setDT(norm_df)
norm_df$...1<-NULL
norm_df$`beta power [-1.7 -1.5]`<-NULL
norm_df$`beta power [-1.5 -1.3]`<-NULL
norm_df$`beta power [-1.3 -1.1]`<-NULL
norm_df$`beta power [-1.1 -0.9]`<- NULL
#norm_df$scheme<-NULL
norm_df$"group"<- "normotypical"

df<- rbind(norm_df, aut_df)

sensor_info <- fread('/Users/kristina/Documents/stc/theta_sensors/sensors.csv', header = TRUE)

names(sensor_info)[1]<- "sensor"
files <- data.table(full_filename=list.files(norm_path, pattern = '*.csv', full.names = T))
files_2<-data.table(full_filename=list.files(autist_path, pattern = '*.csv', full.names = T))
files$short_filename <- list.files(norm_path, pattern = '*.csv', full.names = F)
files_2$short_filename <- list.files(autist_path, pattern = '*.csv', full.names = F)


files[, sensor:=stri_extract_first_regex(short_filename,"[0-9]+")]
files[, interval:=str_extract(short_filename,'[0-9]+_[0-9]+.csv')]
# files[,interval:=gsub('.csv','',interval)]
files$sensor <- as.integer(files$sensor)
files <- merge.data.table(files,sensor_info,by = c('sensor'))
files$interval<-NULL
files$Name<-NULL

files_2[, sensor:=stri_extract_first_regex(short_filename,"[0-9]+")]
files_2[, interval:=str_extract(short_filename,'[0-9]+_[0-9]+.csv')]
# files[,interval:=gsub('.csv','',interval)]
files_2$sensor <- as.integer(files_2$sensor)
files_2 <- merge.data.table(files_2,sensor_info,by = c('sensor'))
files_2$Name<-NULL
files_2$interval<-NULL

all_sensors<- rbind(files,files_2)
names(df)[25]<- "full_filename"
df<- merge.data.table(all_sensors,df,by = c('full_filename'))

df$trial_type<- as.factor(df$trial_type)
df$feedback_cur<- as.factor(df$feedback_cur)
df$group<- as.factor(df$group)
df_lp_hp<- subset(df, trial_type =="risk" | trial_type =="norisk")

################# Subset the participants just with all trial_types and feedbacks ###############

####### All pairs ###########
new_df<- subset(df, subject== 'P301'|subject=='P304'| subject=='P307'|subject=='P332'|subject=='P312'| subject=='P313'|subject== 'P314'|subject=='P316'| subject=='P318'|subject=='P320'| subject=='P321'|subject== 'P322'|subject=='P323'|
                    subject=='P325'|subject=='P326'|subject=='P327'|subject=='P328'|subject=='P329'|subject=='P333'| subject=='P334'|subject=='P335'|subject=='P341'|subject=='P342'|
                    subject=='P324'|subject=='P023'|subject== 'P053'|subject=='P022'|subject=='P016'|subject=='P040'|subject=='P065'|subject=='P001'|subject=='P064'|subject=='P055'|subject=='P060'|subject=='P019'|subject== 'P034'| 
                    subject=='P004'|subject=='P039'|subject=='P035'|subject== 'P008'|subject=='P047'|subject== 'P059'| subject=='P021'|subject== 'P063'|subject=='P032'|subject== 'P044'| subject==  'P061'|subject== 'P067')
                    
#### Pairs with all condition ############
new_df<- subset(df_lp_hp,subject=='P301'|subject== 'P304'|subject== 'P307'|subject== 'P332'|subject== 'P318'| 
                  subject=='P321'|subject== 'P323'|subject=='P325'|subject=='P327'|subject=='P328'|
                  subject=='P329'|subject=='P334'| subject=='P335'| subject=='P308'|subject=='P341'|
                  subject=='P023'|subject=='P053'|subject=='P022'|subject=='P016'|subject=='P055'|subject=='P019'|subject=='P004'|
                  subject=='P039'|subject=='P008'|subject=='P047'|subject=='P059'|subject=='P063'|
                  subject=='P032'|subject== 'P017'|subject=='P044')
                    
                   

cols <- colnames(new_df)[grep('[0-9]',colnames(df))]
sensors<- unique(new_df$sensor)
p_vals <- data.table()

############### for green heads (main_effects) ##############
for (i in 1:length(sensors)) {
  
  temp <-subset(new_df, sensor == sensors[i])
  #print(temp)
  for (j in cols){
    m <- lmer(get(j) ~ group*trial_type*feedback_cur + (1|subject), data = temp)
    an <- anova(m)
    an <- data.table(an,keep.rownames = TRUE)
    an_cols <- c('rn','Pr(>F)') 
    an <- an[, ..an_cols]
    an$`Pr(>F)` <- format(an$`Pr(>F)`, digits = 3)
    an$interval <- j
    an$interval <- gsub('theta power','',an$interval)
    an <- dcast(an,formula = interval~rn,value.var = 'Pr(>F)')
    an$sensor <- sensors[i] 
    #an$sensor_name <- files[sensor==i]$Name
    p_vals <- rbind(p_vals,an)
  }
}

setwd('/Users/kristina/Documents/stc/theta_sensors/theta_4_7/compare')
write.csv(p_vals, "trial_type_lp_hp.csv")

####################### Tukey tests ##################
p_vals <- data.table()
for (i in 1:length(sensors)) {
  temp <-subset(new_df, sensor == sensors[i])
  
  
  for (j in cols) {
    m <- lmer(get(j) ~ group*feedback_cur*trial_type+ (1|subject), data = temp)
    #m_sc <- update(m,data=temp_resc)
    Tuk <- data.table(summary(emmeans(m, pairwise ~ group|trial_type|feedback_cur, adjust = 'tukey',lmer.df = "satterthwaite",lmerTest.limit=8000))$contrasts)
    #print(Tuk)
    Tuk[,contrast:=gsub(' - ','_',contrast)] 
    Tuk[,p.value:=format(p.value, digits = 3)]
    #Tuk[,contrast:=paste0(scheme,'-',contrast)]
    Tuk[,contrast:=paste0(trial_type,'-',feedback_cur,'-',contrast)] ####### If you need triple interaction !!!!!!!!!
    #Tuk[,contrast:=paste0(trial_type,'-',contrast)]  ####### If you need double interaction (for example trial type and feedback) !!!!!!!!!
    columns <- c('contrast','p.value')
    Tuk <- Tuk[,..columns]
    
    Tuk$interval <- j
    Tuk$interval <- gsub('theta power','',Tuk$interval)
    
    Tuk <- dcast(Tuk,formula = interval~contrast,value.var = 'p.value',fun.aggregate = NULL)
    Tuk$sensor <- sensors[i] 
    p_vals <- rbind(p_vals,Tuk)
  }
}

setwd('/Users/kristina/Documents/stc/theta_sensors/theta_4_7/compare')
write.csv(p_vals, "lp_hp_tukey.csv")








