
########## Create the boxplots based on lmem statistics (marginal values) ################

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
library(ggsignif)
options(scipen = 999)
library(ggpattern)
library(cowplot)

sterr <- function(x) sd(x)/sqrt(length(x))

path_p <- '/Users/kristina/Documents/stc/lmem_sensors' # path to files with p-vals
path <- '/Users/kristina/Documents/stc/dataframe_for_LMEM/dataframe_for_LMEM_beta_16_30_trf_early_log'
setwd('/Users/kristina/Documents/stc/lmem_sensors')
out_path <- '/Users/kristina/Documents/stc/lmem_sensors' ## path to save tables and pictures

## load subj_list ##
subj_list <- fread('/Users/kristina/Documents/stc/subj_list.csv')


#df[, interaction_fdr:=p.adjust(`trial_type:feedback_cur`, method = 'fdr')]

# choose sensors signigicant in both intervals

############## CHOOSE SENSOR GROUP #############
#sensors_all <- intersect(sensors_1,sensors_2) 
#post response : 0.700, 0.900
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


######## collect data and average (don't forget to choose needed sensor number and time interval) #############

temp <- fread(files[sensor==9]$full_filename) #donor of colnames for empty datatable "beta"
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

beta[,`mean beta power [-0.9 -0.7]`:=rowMeans(beta[,.SD,.SDcol=c("beta power [-0.9 -0.7]","beta power [-0.7 -0.5]","beta power [-0.5 -0.3]")])] # mean of intervals 


beta[, index := 1:.N, by=c('subject','sensor')] #set indexes of all trials of each subject, same for all sensors 


means <- beta[, mean(`mean beta power [-0.9 -0.7]`),by=c('subject','index')] # compute means of sensors

cols <- c("subject","round","trial_type","feedback_cur","feedback_prev","scheme",'index')
means<- merge.data.table(means,beta[sensor==9, ..cols], by = c('subject','index'), all.x = TRUE) # take trial classification from "beta", sensor is random you can take any

means$subject <- as.factor(means$subject)
means$round <- as.factor(means$round)
means$feedback_cur <-as.factor(means$feedback_cur)
means$feedback_prev <-as.factor(means$feedback_prev)
means$scheme <-as.factor(means$scheme)
means$trial_type <- as.factor(means$trial_type)

setnames(means,'V1','mean_beta')


####### Set the model ########
emm_options(lmerTest.limit = 6000)

m <- lmer(mean_beta ~ trial_type*feedback_cur+ (1|subject), data = means) # main part, fit model!
summary(m)

### Take the marginal means from the model ########
emm_options(pbkrtest.limit = 5000)
marginal_em <- emmeans(m2, ~ as.factor(trial_type|feedback_cur), level = 0.99)
marginal_em<- as.data.frame(marginal_em)

######## Take statistics from model ######
an <- NULL
an <- anova(m)
an <- data.table(an,keep.rownames = TRUE)
#an[, eta2:=F_to_eta2(`F value`, NumDF, DenDF)$Eta2_partial]
an[`Pr(>F)`<0.001, stars:='***']
an[`Pr(>F)`<0.01 & `Pr(>F)`>0.001 , stars:='**']
an[`Pr(>F)`<0.05 & `Pr(>F)`>0.01 , stars:='*']
an[`Pr(>F)`>0.05 & `Pr(>F)`<0.1 , stars:='#']

####### POST HOC tests #########
Tuk<- NULL
thr1 <- max(means[, mean(mean_beta) + sterr(mean_beta), by=c('trial_type')]$V1) 
thr1 <- thr1+0.02 #for RT

thr1_min <- min(means[!is.na(mean_beta), mean(mean_beta) - sterr(mean_beta), by=c('trial_type')]$V1) 

Tuk<-data.table(summary(emmeans(m, pairwise ~ feedback_cur|trial_type, adjust = 'tukey',lmer.df = "satterthwaite"))$contrasts)
Tuk <- Tuk[, group1:=gsub(' -.*', '', contrast)][, group2:=gsub('.*- ', '', contrast)]
Tuk <- Tuk[p.value<0.1, p_significant:=format(p.value, digits = 3)]

n <- Tuk[!is.na(p_significant), .N]

Tuk[p.value<0.001, stars:='***']
Tuk[p.value<0.01 & p.value>0.001 , stars:='**']
Tuk[p.value<0.05 & p.value>0.01 , stars:='*']
Tuk[p.value>0.05 & p.value<0.1 , stars:='#']

if (n>1){
  Tuk <- Tuk[!is.na(p_significant), y.position := seq((thr1+0.01), (thr1+0.3), 0.29/(n-1))]
} else {
  Tuk <- Tuk[!is.na(p_significant), y.position := thr1+0.1]
}
y.position<-Tuk$y.position

Tuk$emmean<-y.position


###### Create the plot just for trial types ########

plot_emmean<-ggplot(data = marginal_em, aes(x = factor(trial_type,level = c("norisk","prerisk","risk","postrisk")), 
                                            y = emmean,  ymin=emmean-SE, ymax = emmean+SE, group = 1))+
  scale_x_discrete(labels = c('HP','Pre-LP','LP', 'post-LP'))+
  geom_point() + geom_errorbar(width = 0.1, size=1.5)+geom_line(size=1.5)+labs(y = "Beta power change,dB", x = "Choice type")+
  theme_classic()+theme(text = element_text(size=20))+ ylim (-0.5,2.5) +
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"))+
  geom_hline(yintercept= 0.0, linetype='dashed', col = 'black', size = 1.0)+
  stat_pvalue_manual(Tuk, label = 'stars', size = 12, bracket.size = 1.5, tip.length = 0.01,y.position =c(1.3,1.5,1.7,1.9,2.1,2.3),inherit.aes = FALSE)

plot_emmean <- ggpar(plot_emmean,
                     ylim = c(-0.5,2.5),
                     font.ytickslab = 30,
                     font.xtickslab = 27,
                     font.main = 25,
                     font.submain = 25,
                     font.x = 27,
                     font.y = 30)

plot_emmean

setwd('/Users/kristina/Documents/stc/plots_article')
ggsave('mean_beta_anterior.png',plot_emmean,width =  6, height = 5)

########## Split current feedbacks ##########
signif <- Tuk[!is.na(stars)]

sequence <-data.table(trial_type=c("norisk","prerisk","risk", "postrisk"),number=c(1,2,3,4))
sterr <- function(x) sd(x)/sqrt(length(x))

y_values_rew <- means[feedback_cur == 'positive',
                      mean(mean_beta)+sterr(mean_beta)+0.7, by='trial_type']
setnames(y_values_rew,'V1','y_values_rew')

y_values_lose <-  means[feedback_cur == 'negative',
                        mean(mean_beta)+sterr(mean_beta)+0.7, by='trial_type']

setnames(y_values_lose,'V1','y_values_lose')

y_values <- merge(y_values_lose,y_values_rew,by='trial_type')
y_values <- merge(y_values,sequence,by='trial_type')
y_values[,y_max:=max(y_values_lose,y_values_rew),by=trial_type]
y_values[,y_min:=min(y_values_lose,y_values_rew),by=trial_type]

# ylim1 <- min(y_values$y_min)
# ylim2 <- max(y_values$y_max)

y_values <- merge(y_values,signif,by='trial_type')

setnames(marginal_em, 'y_values_lose', "emmean")

p1 <- ggplot(marginal_em, aes(x = factor(trial_type,level = c("norisk","prerisk","risk","postrisk")),
                              y = emmean,  ymin=emmean-SE, ymax = emmean+SE, color = feedback_cur,group = feedback_cur))+
  scale_x_discrete(labels = c('HP','Pre-LP','LP', 'post-LP'))+ geom_line(size=1.5)+
  geom_point(position=position_dodge(0.1)) + geom_errorbar(width = 0.1,  position=position_dodge(0.1), size=1.5)+labs(y = "Beta power change,dB", x = "Choice type")+
  theme_classic()+ theme(text = element_text(size=20))+scale_color_discrete(name = "Current feedback", labels = c("Loss", "Gain"))+theme(legend.position="bottom")+
  ylim(-1.0,2.6) +
  geom_hline(yintercept=-0.0, linetype='dashed', col = 'black', size = 1.0)+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"))

p1 <- p1+geom_signif(y_position=c(y_values$y_max +0.05),
                     xmin=c(y_values$number-0.075), xmax=c(y_values$number+0.075),
                     annotation=c(y_values$stars),col='black',
                     tip_length=0.003,textsize = 12,vjust = 0.4,size = 1.2)

p1<- ggpar(p1,
           ylim = c(-1.0,2.6),
           font.ytickslab = 30,
           font.xtickslab = 27,
           font.main = 25,
           font.submain = 25,
           font.x = 27,
           font.y = 27)

p1

setwd('/Users/kristina/Documents/stc/plots_article')
ggsave('mean_beta_anterior_1500_1900_fb.png',p1,width =  6, height = 5)


################### Compare feedback prev & feedback cur ######################
emm_options(lmerTest.limit = 6000)
lp<- filter(means, trial_type %in% "risk")
hp<- filter(means, trial_type %in% "norisk")
prelp<- filter(means, trial_type %in% "prerisk")
postlp<- filter(means, trial_type %in% "postrisk")

lp_hp<- filter(means, trial_type %in% "risk" | trial_type %in% "norisk")

lp_hp_lmem<- lmer(mean_beta ~ trial_type*feedback_prev*feedback_cur+ (1|subject), data = lp_hp)
lp_lmem <- lmer(mean_beta ~ feedback_prev*feedback_cur+ (1|subject), data = lp)
hp_lmem <- lmer(mean_beta ~ feedback_prev*feedback_cur+(1|subject), data = hp)
prelp_lmem <- lmer(mean_beta ~ feedback_prev*feedback_cur+ (1|subject), data = prelp)
postlp_lmem <- lmer(mean_beta ~ feedback_prev*feedback_cur+ (1|subject), data = postlp)
# main part, fit model!

an <- NULL
an <- anova(lp_hp_lmem)
an <- data.table(an,keep.rownames = TRUE)
#an[, eta2:=F_to_eta2(`F value`, NumDF, DenDF)$Eta2_partial]
an[`Pr(>F)`<0.001, stars:='***']
an[`Pr(>F)`<0.01 & `Pr(>F)`>0.001 , stars:='**']
an[`Pr(>F)`<0.05 & `Pr(>F)`>0.01 , stars:='*']
an[`Pr(>F)`>0.05 & `Pr(>F)`<0.1 , stars:='#']

marginal_lp_hp <- emmeans(lp_hp_lmem, ~ as.factor(trial_type|feedback_prev|feedback_cur), level = 0.99)
marginal_lp_hp<- as.data.frame(marginal_lp)
marginal_lp<- filter(marginal_lp_hp, trial_type %in% 'risk')
marginal_hp<- filter(marginal_lp_hp, trial_type %in% 'norisk')

marginal_hp <- emmeans(hp_lmem, ~ as.factor(feedback_cur|feedback_prev), level = 0.99)
marginal_hp<- as.data.frame(marginal_hp)
marginal_postlp <- emmeans(postlp_lmem, ~ as.factor(feedback_prev|feedback_cur), level = 0.99)
marginal_postlp<- as.data.frame(marginal_postlp)
marginal_prelp <- emmeans(prelp_lmem, ~ as.factor(feedback_prev|feedback_cur), level = 0.99)
marginal_prelp<- as.data.frame(marginal_prelp)


hp_lp <- filter(means, trial_type %in% 'risk' | trial_type %in% 'norisk')
lmem<- lmer(mean_beta ~ trial_type*feedback_prev*feedback_cur+ (1|subject), data = hp_lp)

marginal <- emmeans(lmem, ~ as.factor(trial_type|feedback_prev|feedback_cur), level = 0.99)
marginal<- as.data.frame(marginal)

marginal_lp<- filter(marginal, trial_type %in% "risk")
marginal_hp<- filter(marginal, trial_type %in% "norisk")

Tuk<- NULL
thr1 <- max(hp_lp[, mean(mean_beta) + sterr(mean_beta), by=c('trial_type')]$V1) 
thr1 <- thr1+0.02 #for RT

thr1_min <- min(hp_lp[!is.na(mean_beta), mean(mean_beta) - sterr(mean_beta), by=c('trial_type')]$V1) 

Tuk<-data.table(summary(emmeans(lmem, pairwise ~ feedback_prev|feedback_cur|trial_type, adjust = 'tukey',lmer.df = "satterthwaite"))$contrasts)

Tuk <- Tuk[, group1:=gsub(' -.*', '', contrast)][, group2:=gsub('.*- ', '', contrast)]
Tuk <- Tuk[p.value<0.1, p_significant:=format(p.value, digits = 3)]

n <- Tuk[!is.na(p_significant), .N]

Tuk[p.value<0.001, stars:='***']
Tuk[p.value<0.01 & p.value>0.001 , stars:='**']
Tuk[p.value<0.05 & p.value>0.01 , stars:='*']
Tuk[p.value>0.05 & p.value<0.1 , stars:='#']

tuk_hp<-filter(Tuk, trial_type %in% "norisk")
tuk_lp<-filter(Tuk, trial_type %in% "risk")


if (n>1){
  Tuk <- Tuk[!is.na(p_significant), y.position := seq((thr1+0.01), (thr1+0.3), 0.29/(n-1))]
} else {
  Tuk <- Tuk[!is.na(p_significant), y.position := thr1+0.1]
}
y.position<-Tuk$y.position

Tuk$emmean<-y.position

emm_options(pbkrtest.limit = 4363)
marginal_lp <- emmeans(lmem, ~ feedback_cur|feedback_prev, level = 0.95)
marginal_lp<- as.data.frame(marginal_lp)
#marginal_fc<- marginal_feedback_choice
marginal_lp_neg<- filter(marginal_hp, feedback_prev %in% "negative")
marginal_lp_pos<- filter(marginal_hp, feedback_prev %in% "positive")
marginal_lp_prev_neg<- filter(marginal_hp,  feedback_prev %in% "negative")
marginal_lp_prev_pos<- filter(marginal_hp, feedback_prev %in% "positive")


signif <- tuk_hp[!is.na(stars)]
signif_neg<- signif[1,]
signif_pos<- signif[2,]



y_values_rew <- lp_hp[feedback_cur == 'positive',
                      mean(mean_beta)+sterr(mean_beta)+0.7, by='feedback_cur']
setnames(y_values_rew,'V1','y_values_rew')

y_values_lose <-  lp_hp[feedback_cur == 'negative',
                        mean(mean_beta)+sterr(mean_beta)+0.7, by='feedback_cur']

setnames(y_values_lose,'V1','y_values_lose')
setnames(signif_neg,'feedback_prev','feedback_cur')
setnames(signif_pos,'feedback_prev','feedback_cur')
y_values_rew[,y_max:=max(y_values_rew)]
y_values_rew[,y_min:=min(y_values_rew)]

y_values_lose[,y_max:=max(y_values_lose)]
y_values_lose[,y_min:=min(y_values_lose)]
# ylim1 <- min(y_values$y_min)
# ylim2 <- max(y_values$y_max)

y_values_rew <- merge(y_values_rew, signif_pos)
y_values_lose <- merge(y_values_lose, signif_neg)

neg <- ggplot(marginal_lp_neg, aes(x=as.factor(feedback_cur), y=emmean,  ymin=emmean-SE, ymax = emmean+SE,fill=feedback_cur))+ 
  geom_col_pattern(aes(pattern = feedback_prev, fill=feedback_cur),width= 1, 
                   fill= c('coral2','darkturquoise'), 
                   colour =  c('black'),pattern="stripe",
                   pattern_aspect_ratio = 1, pattern_fill= "black",
                   pattern_density= 0.1)+labs(y = "Beta power change,dB",col="black")+geom_line(size=1.5)+
  geom_point(position=position_dodge(0.07)) + geom_errorbar(width = 0.1,  position=position_dodge(0.1), size=1.5)+
  theme_classic()+ scale_x_discrete(labels=c("Loss", "Gain"))+
  geom_hline(yintercept= 0.0, linetype='dashed', col = 'black', size = 2.0)+ylim(-1, 2.0)+
  theme(axis.text.x = element_text(colour="black"), axis.text.y = element_text(colour="black"))+
  
  theme(legend.position = "none")

neg<- ggpar(neg,
            ylim = c(-1, 2.0),
            font.ytickslab = 30,
            font.xtickslab = 27,
            font.main = 25,
            font.submain = 25,
            font.x = 27,
            font.y = 27)

neg<- neg+geom_signif(y_position=c(1.40),
                      xmin=c(1.15), xmax=c(1 + 0.75),
                      annotation=c(y_values_lose$stars),col='black',
                      tip_length=0.003,textsize = 8,vjust = 0.7,size = 1.2,label.size = 10)

neg


pos <-ggplot(marginal_hp_pos, aes(x=as.factor(feedback_cur), y=emmean,  ymin=emmean-SE, ymax = emmean+SE,fill=feedback_cur))+ 
  geom_col_pattern(aes(pattern = feedback_prev, fill=feedback_cur),width= 1, 
                   fill= c('coral2','darkturquoise'), 
                   colour =  c('black'),pattern="crosshatch",
                   pattern_aspect_ratio = 1, pattern_fill= "black",
                   pattern_density= 0.1)+labs(x="Current feedback",col="black")+geom_line(size=1.5)+
  geom_point(position=position_dodge(0.07)) + geom_errorbar(width = 0.1,  position=position_dodge(0.1), size=1.5)+
  theme_classic()+ scale_x_discrete(labels=c("Loss", "Gain"))+
  geom_hline(yintercept= 0.0, linetype='dashed', col = 'black', size = 2.0)+ylim(-1, 2.0)+
  theme(axis.text.x = element_text(colour="black"), axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  theme(legend.position="none")+labs(y=NULL)+guides(y = "none")

pos<- ggpar(pos,
            ylim = c(-1, 2.0),
            font.ytickslab = 30,
            font.xtickslab = 27,
            font.main = 25,
            font.submain = 25,
            font.x = 27,
            font.y = 27)

pos<- pos+geom_signif(y_position=c(1.68),
                      xmin=c(1.15), xmax=c(1 + 0.99),
                      annotation=c(y_values_rew$stars),col='black',
                      tip_length=0.003,textsize = 8,vjust = 0.7,size = 1.2,label.size = 10)

pos



pn<- plot_grid(neg,pos+ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                                      axis.ticks.y = ggplot2::element_blank(),
                                      axis.title.y = ggplot2::element_blank()), nrow = 1, align = "hv")


pn<- pn+geom_bracket(xmin = 0.20, xmax = 0.70, y.position = 0.96, label = "*", col='black',
                     size =1.2,label.size = 10,vjust = 0.7)

pn<- pn+geom_bracket(xmin = 0.40, xmax = 0.85, y.position = 0.92, label = "*",col='black',
                     size =1.2,label.size = 10,vjust = 0.5)                   
pn


setwd('/Users/kristina/Documents/stc/plots_article')
ggsave('HP_prev_cur_anterior.png',pn,width =  6, height = 5)
plotl <- function(...) {
  x <- seq(0, 30, 0.01)
  plot(besselJ(x, 0), col = 2, type = "l",
       lwd = 2, ylab = "Jn(x)", xlab = "", ...)
  lines(besselJ(x, 2), col = 3, type = "l", lwd = 2, lty = 2)
}

plotl()
legend("topright",
       legend = c(3, 4, 5),
       fill = 2:4,       # Color of the squares
       border = "black")

legend("topright",title="Current feedback", legend = c("Loss", "Gain"),
       fill = c("coral2","darkturquoise"))

