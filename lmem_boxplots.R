library(reshape2)
library(data.table)
library(ggplot2)
library(ggpubr)
library(lme4)
library(ggpubr)
library(emmeans)
library(lmerTest)
library(plotrix)
library(stringi)
library(gridExtra)
library(ggprism)
library(dplyr)
options(scipen = 999)

path_p <- '/Users/kristina/Documents/stc/lmem_sensors' # path to files with p-vals
path <- '/Users/kristina/Documents/stc/dataframe_for_LMEM/df_lmem'
setwd('/Users/kristina/Documents/stc/lmem_sensors')
out_path <- '/Users/kristina/Documents/stc/lmem_sensors' ## path to save tables and pictures

## load subj_list ##
subj_list <- fread('/Users/kristina/Documents/stc/subj_list.csv')


df <- fread ('/Users/kristina/Documents/stc/LMEM_trial_cur_fb.csv')
df$V1 <- NULL

# fdr 
df[, interaction_fdr:=p.adjust(`trial_type:feedback_cur`, method = 'fdr')]

# choose sensors signigicant in both intervals
sensors_1 <- unique(df[interaction_fdr <=0.05 | interval== '[1.5 1.7]']$sensor_name)
sensors_2 <- unique(df[df$interaction_fdr <=0.05 | df$interval== '[1.7 1.9]']$sensor_name)

sensors_all <- intersect(sensors_1,sensors_2) 

sensors_all<- fread('/Users/kristina/Documents/stc/dataframe_for_LMEM/sensor_ttype.csv', header=TRUE)
sensors_all$V6<- NULL

####### make dataframe with files info ########

sensor_info <- fread('/Users/kristina/Documents/stc/sensors.csv', header = TRUE)
names(sensor_info)[1] <- 'sensor'

files <- data.table(full_filename=list.files(path, pattern = '*.csv', full.names = T))
files$short_filename <- list.files(path, pattern = '*.csv', full.names = F)

# files$short_filename <- gsub('planar2','',files$short_filename)

files[, sensor:=stri_extract_first_regex(short_filename,"[0-9]+")]
# files[, interval:=str_extract(short_filename,'[0-9]+_[0-9]+.csv')]
# files[,interval:=gsub('.csv','',interval)]
files$sensor <- as.integer(files$sensor)
files <- merge.data.table(files,sensor_info,by = c('sensor'))
names(files)[4] <- 'sensor_name'
files$effect <- NULL


##### filter files and leave needed sensors only ######

files <- files[sensor_name %in% sensors_all$sensor_name] 


######## collect data and average #############

temp <- fread(files[sensor==5]$full_filename) #donor of colnames for empty datatable "beta"
temp$V1 <- NULL
beta <- setNames(data.table(matrix(nrow = 0, ncol = length(colnames(temp))+2)), c(colnames(temp),'sensor','sensor_name'))

for (i in files$sensor){
  temp <- fread(files[sensor==i]$full_filename)
  temp$V1 <- NULL
  temp <- as.data.table(temp)
  temp <- temp[subject %in% subj_list$subj_list]
  
  temp$sensor <- i
  temp$sensor_name <- files[sensor==i]$sensor_name
  
  beta <- rbind(beta,temp)
}

beta[,`mean beta power [1.5 1.7]`:=rowMeans(beta[,.SD,.SDcol=c("beta power [1.5 1.7]","beta power [1.7 1.9]")])] # mean of intervals 


beta[, index := 1:.N, by=c('subject','sensor')] #set indexes of all trials of each subject, same for all sensors 


means <- beta[, mean(`mean beta power [1.5 1.7]`),by=c('subject','index')] # compute means of sensors

cols <- c("subject","round","trial_type","feedback_cur","feedback_prev","scheme",'index')
means<- merge.data.table(means,beta[sensor==6, ..cols], by = c('subject','index'), all.x = TRUE) # take trial classification from "beta", sensor is random you can take any

means$subject <- as.factor(means$subject)
means$round <- as.factor(means$round)
means$feedback_cur <-as.factor(means$feedback_cur)
means$feedback_prev <-as.factor(means$feedback_prev)
means$scheme <-as.factor(means$scheme)
means$trial_type <- as.factor(means$trial_type)

setnames(means,'V1','mean_beta')

###########СОЗДАЕМ МОДЕЛЬ#############
f <- seq(from=-2.5, to=-0.5, by=0.1)
interval <- 'mean_beta' #name of dependent variable

m <- lmer(get(interval) ~ trial_type*feedback_cur + (1|subject), data = means) # main part, fit model!
m2 <- m # it is because we missed step, so final model=initial model


####Усредняем, способ номер 1 (внутри испытуемого и трайлов усредняем фидбэк, затем между фидбэками)
fb<- means %>% 
  group_by(subject, trial_type, feedback_cur) %>% 
  summarize(mean_beta=mean(mean_beta), sd = sd(mean_beta), se = sd / sqrt(n())) #Усредняем внутри каждого испытумого, внутри трайлов , внутри фидбэков

fb_1<- fb%>% 
  group_by(trial_type, subject)%>% 
  summarize(mean_beta = mean(mean_beta), sd = mean(sd), se = mean(se)) #Усредняем между фидбэками

fb_1<- fb_1%>% 
  group_by(trial_type)%>% summarize(mean_beta= mean(mean_beta), sd = mean(sd), se = mean(se)) #Усредняем между испытуемыми

plot_fb1<- ggplot(data = fb_1, aes(x = factor(trial_type,level = c("norisk","prerisk","risk","postrisk")), y = mean_beta, label=mean_beta, group = 1))+geom_point() +
  geom_line()+geom_text(hjust=0, vjust=0)+ ylim(-2.5, -0.5)+ labs(y = "Усреднение (вариант 1)")

plot_fb1

####Усредняем, способ номер 2 (усредняем фидбэки насквозь по траелам)

fb_2<- means %>% 
  group_by(trial_type, feedback_cur) %>% 
  summarize(mean_beta = mean(mean_beta),sd = sd(mean_beta), se = sd / sqrt(n()))# усредняем фидбэки внутри трайлов

fb_2<- means %>% 
  group_by(trial_type) %>% 
  summarize(mean_beta = mean(mean_beta),sd = sd(mean_beta), se = sd / sqrt(n())) #усредняем между фидбэками

plot_fb2<- ggplot(data = fb_2, aes(x = factor(trial_type,level = c("norisk","prerisk","risk","postrisk")), y = mean_beta, label=mean_beta, group = 1))+geom_point() +
  geom_line()+geom_text(hjust=0, vjust=0)+ ylim(-2.5, -0.5)+ labs(y = "Усреднение альтернативное (вариант 2)")

plot_fb2 
          
######Усредняем насквозь 3 ##### 
fun_mean <- function(x){
  return(data.frame(y=mean(x),label=mean(x,na.rm=T)))}

p<- ggplot(data = means, aes(x = factor(trial_type,level = c("norisk","prerisk","risk","postrisk")), y = mean_beta))+geom_line()+ 
  stat_summary(fun =mean, geom="point",colour="darkred", size=3) +
  stat_summary(fun.data = fun_mean, geom="text", vjust=-0.7) 
p
p+ ylim(-2.5, -0.5)#тут лучше не приводить к общей шкале, тк реальный разброс выходит от -5 до 5

####### Маргинальные значения из emmeans(усреднение 4)
emm_options(pbkrtest.limit = 5000)
marginal_em <- emmeans(m2, ~ as.factor(trial_type), level = 0.99)
marginal_em<- as.data.frame(marginal_em)

plot_emmean<-ggplot(data = marginal_em, aes(x = factor(trial_type,level = c("norisk","prerisk","risk","postrisk")), y = emmean, label=emmean, group = 1))+
  geom_point() +geom_pointrange(aes(ymax = upper.CL, ymin = lower.CL))+
  geom_line()+geom_text(hjust=0, vjust=0)+ ylim(-2.5, -0.5)+ labs(y = "EMMEN")

plot_emmean






#marginal_split <- emmeans(m2, ~ as.factor(trial_type|feedback_cur), level = 0.99)






