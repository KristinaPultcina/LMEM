############# Create the plot of beta power and response time interaction based on lmem stat #########

library(reshape2)
library(data.table)
library(ggplot2)
library(ggpubr)
library(lme4)
library(emmeans)
library(lmerTest)
library(plotrix)
library(stringi)
library(dplyr)
library(tibble)
library(sjPlot)
library(glmmTMB)
library(MuMIn)
require(rgl)
library(ggiraphExtra)
library(predict3d)
library(ggeffects)
options(scipen = 999)
if(!require(devtools)) install.packages("devtools")
devtools::install_github("cardiomoon/predict3d")


rt<- read.csv("/Users/kristina/Documents/stc/dataframe_for_LMEM/rt/events_classification_clean_final.csv", sep= ";")
rt<- na.omit(rt)
rt<-rt[!( rt$label.type=="Stimulus (right side)"|rt$label.type == "Button (left)"| rt$label.type == "Stimulus (left side)"|rt$label.type == "	
Feedback LOSE INV"| rt$label.type == "Button (right)"| rt$label.type == "Feedback REW INV" | rt$label.type == "Feedback LOSE INV"|
            rt$response.class == "risk and prerisk"| rt$response.class == "risk and postrisk"| rt$response.class == "postrisk and prerisk"| rt$response.class=="last" | rt$label.type == "Show Bank"|
            rt$response.class=="risk, postrisk and prerisk"| rt$response.class=="first"),]



my_df <- rt %>%
  filter(response.time > 300 & response.time < 4000)

my_df<- filter(my_df, subj %in% c("P001", "P002", "P004", "P006", "P007", "P008", "P011", "P014", "P015", "P016", "P017",
                                  "P019", "P021","P022", "P023","P024", "P025", "P028", "P029", "P030", "P031", "P032",
                                  "P033", "P034","P035","P039", "P040", "P042", "P043", "P044","P045", "P047", "P048", "P052",
                                  "P053", "P055", "P057", "P059", "P060", "P062"))


#Create dataset with trained trials###########

trained_df<-my_df[!(my_df$learning.criteria == "learning"),]
trained_df$learning.criteria<- NULL
setwd('/Users/kristina/Documents/stc')
write_csv(trained_df, "response_time_df.csv")
############ Correlation beta&RT ###########
path <- "/Users/kristina/Documents/stc/dataframe_for_LMEM/rt/dataframe_for_LMEM_beta_16_30_trf_early_log_with_times_of_resp"
#path<- '/Users/kristina/Documents/stc/dataframe_for_LMEM/not_trained/dataframe_for_LMEM_beta_16_30_trf_early_log_not_trained'
out_path <- '/Users/kristina/Documents/stc/lmem_sensors' ## path to save tables and pictures
## load subj_list ##
subj_list <- fread('/Users/kristina/Documents/stc/subj_list.csv')


df <- fread ('/Users/kristina/Documents/stc/LMEM_trial_cur_fb.csv')
df$V1 <- NULL

# fdr 
df[, interaction_fdr:=p.adjust(`trial_type:feedback_cur`, method = 'fdr')]

######################## CHOOSE sensors signigicant in both intervals #############

#sensors_all<- fread('/Users/kristina/Documents/stc/dataframe_for_LMEM/sensors/occipital_cluster_1.5_1.9.csv', header=TRUE)
#sensors_all<- fread('/Users/kristina/Documents/stc/dataframe_for_LMEM/sensors/anterior_cluster(fig.6).csv', header=TRUE)
#sensors_all<- fread('/Users/kristina/Documents/stc/dataframe_for_LMEM/sensors/1.100_1.500.csv', header=TRUE)
#sensors_all<- fread('/Users/kristina/Documents/stc/dataframe_for_LMEM/sensors/cluster_feedback_anticipation(fig.4).csv', header=TRUE)
sensors_all<- fread('/Users/kristina/Documents/stc/dataframe_for_LMEM/sensors/cluster_decision_making(Fig.3).csv', header=TRUE)
#sensors_all<- fread('/Users/kristina/Documents/stc/dataframe_for_LMEM/sensors/0.1_0.5(before_feed).csv', header=TRUE)
#sensors_all<- fread('/Users/kristina/Documents/stc/dataframe_for_LMEM/sensors/1500_1900_all.csv', header=TRUE)




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

beta[,`mean beta power [-0.9 -0.7]`:=rowMeans(beta[,.SD,.SDcol=c("beta power [-0.9 -0.7]","beta power [-0.7 -0.5]", "beta power [-0.5 -0.3]")])] # mean of intervals 


beta[, index := 1:.N, by=c('subject','sensor')] #set indexes of all trials of each subject, same for all sensors 


means <- beta[, mean(`mean beta power [-0.9 -0.7]`),by=c('subject','index')] # compute means of sensors


cols <- c("subject","round","trial_type","trial_number","feedback_cur","feedback_prev","scheme",'index', "times_of_response")
means<- merge.data.table(means,beta[sensor==5, ..cols], by = c('subject','index'), all.x = TRUE) # take trial classification from "beta", sensor is random you can take any

means$subject <- as.factor(means$subject)
means$round <- as.factor(means$round)
means$feedback_cur <-as.factor(means$feedback_cur)
means$feedback_prev <-as.factor(means$feedback_prev)
means$scheme <-as.factor(means$scheme)
means$trial_type <- as.factor(means$trial_type)
means$times_of_response <- as.factor(means$times_of_response)
setnames(means,'V1','mean_beta')
setnames(means,'real.time','times_of_response')


setnames(trained_df, "run", "round")
setnames(trained_df,"response.class", "trial_type")
#setnames(trained_df,"trial_number", "trial_type")
setnames(trained_df, "cur_fb", "feedback_cur")
setnames(trained_df, "prev_fb", "feedback_prev")
setnames(trained_df, "subj", "subject")
setnames(trained_df, "response.time", "rt")
setnames(trained_df, "event.0...time.", "times_of_response")

trained_df$rt<- log10(trained_df$rt)

trained_df$times_of_response <- as.factor(trained_df$times_of_response)
new_dataset <- trained_df %>%inner_join(means, by=c("subject", "round", "trial_type","feedback_cur","times_of_response"))
new_df<- new_dataset %>%
  filter(mean_beta > 5 & mean_beta < -10)
#write.csv(new_dataset, "anterior_1500_1900.csv")

######Create lmem model (подставляем тип трайла, который интересен)
lmem<- lmer(data=new_dataset, mean_beta ~ rt*trial_type+(1|subject)) # модель только для RT и mean_beta

gg<- ggpredict(lmem,terms=c("rt"),rawdata = TRUE,se=TRUE,xpos=0.5)
gg<-gg[-5,]
MuMIn::r.squaredGLMM(lmem)
summary(lmem)
ggeffect(lmem)
mydf <- ggpredict(lmem,terms="mean_beta")
plot(my_df, rawdata = TRUE)
ggplot(mydf, aes(x, predicted)) + geom_line()


p<-gg %>% 
  plot(add.data = TRUE, dot.size = 0.6)+labs(y = "Beta power change,dB", x = "Response time,ms", title = NULL)+
  theme_classic()+
  theme(text = element_text(size = 27))+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"))+
  scale_y_continuous(round(seq(min(new_dataset$mean_beta), max(new_dataset$mean_beta, by=1.5))))+
  scale_x_continuous(breaks = seq(min(new_dataset$rt), max(new_dataset$rt), by = 0.3), labels = c("350","500","1200","2300"))+
  geom_hline(yintercept= 0, linetype='dashed', col = 'black', size = 2)+scale_y_continuous("Beta power change,dB")  

p<- ggpar(p,
            font.ytickslab = 30,
            font.xtickslab = 25,
            font.main = 35,
            font.submain = 27,
            font.x = 30,
            font.y = 27)
p
setwd('/Users/kristina/Documents/stc/plots_article')
ggsave('anterior_corr.png',p,width =  6, height = 5)

round(seq(min(new_dataset$mean_beta), max(new_dataset$mean_beta), by = 1.5))
seq(min(new_dataset$rt), max(new_dataset$rt), by = 0.2)



###### SPEARMAN CORR ##########


cor.test(new_dataset$rt,new_dataset$mean_beta) ### просто распечатываем корреляции

g <- ggscatter(new_dataset, x = "rt", y = "mean_beta", 
               add = "reg.line", conf.int = TRUE, 
               cor.coef = TRUE, cor.method = "spearman",
               xlab = "RT,ms", ylab = "Beta power, dB") 
g


