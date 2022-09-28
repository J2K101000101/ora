rm(list = ls())
library(pacman)
p_load(rjags,ggplot2,bayesplot,ggmcmc,rstan,rstatix)
load.module('wiener','/usr/lib/JAGS/modules-4/')

workspace = getwd()
outputFolder = paste(workspace, "output", sep = "/") 

# data ####
workspace <- getwd()
read <- paste(workspace,'cleanedTraining','cleanedTrainingData.csv',sep = '/')
df <- read.csv(read)[,-1]
# Convert RTs to seconds: rt of first responses; RT of second responses (in dual)
# a,b,c,d is the training groups: simpleRT,choiceRT, switch, dual
df$rt[df$group %in% c('a','b')] <- df$reactionTime[df$group %in% c('a','b')]/1000
df$rt[df$group == 'c'& df$shiftingType == 'repetition'] <- df$reactionTime[df$group == 'c' & df$shiftingType == 'repetition']/1000
df$RT[df$group == 'c'& df$shiftingType == 'switch'] <- df$reactionTime[df$group == 'c' & df$shiftingType == 'switch']/1000
df$rt[df$group == 'd'] <- as.numeric(df$reactionTime1[df$group == 'd'])/1000
df$RT[df$group == 'd'] <- as.numeric(df$reactionTime2[df$group == 'd'])/1000

# Convert score: score of first responses; SCORE second responses (in dual)
df$score[df$group == 'd'] <- df$score1[df$group == 'd']
df$SCORE[df$group == 'd'] <- df$score2[df$group == 'd']

# Code RTs negatively if error response
df$rt[df$score == 0] = -df$rt[df$score == 0]
df$RT[df$SCORE == 0 & df$group == 'd'] = -df$RT[df$SCORE == 0 & df$group == 'd']

# Fit in Data ####

# MCMC settings
nBurnSamples = 500     # How many samples to discard as burn-in
nChains = 4            # How many samplers
nThin = 10             # How often to thin
nSamples = 10000       # How many samples
nSim = 10000          # How many trials to simulate for posterior predictive check

# Groups = unique(df$group)

Fit <- function(df,Groups,type){
  allDf <- c()
  for (g in Groups){# for each group
    data_by_group <- subset(df,df$group == g)
    for (id in unique(data_by_group$extId)){# for each participant, as everyone might have different nSessions
      data_by_id <- subset(data_by_group,data_by_group$extId == id)
      
      # re-code sessionId
      nSessions <- length(unique(data_by_id$sessionId))
      data_by_id <- data_by_id %>% mutate(session = dense_rank(sessionId))
      
      
      # set up data
      if (type == 1){
        jagsData <- list(session = data_by_id$session, RT = data_by_id$RT, nSessions = nSessions, nTrials = length(data_by_id$RT))
      }else{
        jagsData <- list(session = data_by_id$session, RT = data_by_id$rt, nSessions = nSessions, nTrials = length(data_by_id$rt))
      }
      
      # Define initial parameters
      inits = list(threshold = rep(0.5, nSessions),
                   ndt = rep(0.00001, nSessions),
                   rate = rep(0.5, nSessions))
      
      # Initialise the model
      model = jags.model("DM_wiener", data = jagsData, inits = inits, n.chains = nChains)
      
      
      # Burn-in
      update(model, n.iter = nBurnSamples)
      
      # Get posterior samples
      mcmcList = coda.samples(model, c("threshold", "ndt", "rate"), nSamples, nThin)
      
      sample <- as.matrix(mcmcList)
      
      if (type == 1){
        sample_file_name <- paste('sample1','_group_',g,'_id_',id,'.csv',sep = '')
      }else{
        sample_file_name <- paste('sample','_group_',g,'_id_',id,'.csv',sep = '')
      }
      sample_file_path <- paste(outputFolder,sample_file_name,sep = '/')
      write.csv(sample, sample_file_path)
      
      # for each participant 
      ndt <-  c()
      threshold <-  c()
      rate <-  c()
      for (n in  1: nSessions){
        ndtColName = paste('ndt','[',n,']',sep = '')
        thresColName = paste('threshold','[',n,']',sep = '')
        rateColName = paste('rate','[',n,']',sep = '')
        ndt[n]=mean(sample[,ndtColName])
        threshold[n]=mean(sample[,thresColName])
        rate[n]=mean(sample[,rateColName])
      }
      parameters <- data_frame(ndt,threshold,rate)
      parmDf <- cbind(extId = id,group = g,session = c(1:nSessions),parameters)
      allDf <- rbind(allDf,parmDf)
      
    }
    if (type == 1){
      group_file_name <- paste('meanPars1','_group_',g,'.csv',sep = '')  
    }else{
      group_file_name <- paste('meanPars','_group_',g,'.csv',sep = '')  
    }
    group_file_path <- paste(outputFolder,group_file_name, sep = '/')
    write.csv(allDf, group_file_path)
  }
  return(allDf)
  
}

allDF_a_rt = Fit(df, c('a'),0)
allDF_a_rt = cbind(allDF_a_rt, type = c('rt'))
allDF_b_rt = Fit(df, c('b'),0)
allDF_b_rt = cbind(allDF_b_rt, type = c('rt'))
allDF_c_rt = Fit(df,c('c'),0)
allDF_c_rt = cbind(allDF_c_rt, type = c('rt'))
allDF_c_RT = Fit(df, c('c'),1)
allDF_c_RT = cbind(allDF_c_RT, type = c('RT'))
allDF_d_rt = Fit(df, c('d'),0)
allDF_d_rt = cbind(allDF_d_rt, type = c('rt'))
allDF_d_RT = Fit(df, c('d'),1)
allDF_d_RT = cbind(allDF_d_RT, type = c('RT'))
training_mean_Pars_group_all <- rbind(allDF_a_rt,allDF_b_rt,allDF_c_rt,allDF_c_RT,allDF_d_rt,allDF_d_RT)
write.csv(training_mean_Pars_group_all,'all_training.csv')
