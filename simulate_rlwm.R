# Load fit object and data
#library(nloptr)
library(tidyverse)
`%notin%` =Negate(`%in%`)

fitData <- read.csv('fitData.csv')
pes = read.csv('pes.csv')
pes = filter(pes, participant %in% fitData$participant)

#include = mss$participant
#pes = filter(pes, participant %in% include)

subData = select(pes, participant, rl_resp.keys, key, load, block, rl_trials.thisIndex, stimRep, trialN)
subData$key_press = ifelse(subData$rl_resp.keys=='j', 1, ifelse(subData$rl_resp.keys=='k', 2, 3))
subData$corkey = ifelse(subData$key=='j', 1, ifelse(subData$key=='k', 2, 3))
subData=select(subData, -rl_resp.keys, -key)

colnames(subData)[2]='nSet'
colnames(subData)[3]='qNum'
colnames(subData)[4]='stimNum'
colnames(subData)[5]='stimIter'

subData$stimNum = subData$stimNum+1
subData$trialN = subData$trialN+1
subData$stimIter = subData$stimIter+1

fit = group_by(fitData, participant) %>% summarize(negLL=min(negLL))
fit = merge(fitData, fit, by=c('participant','negLL'))
fit = group_by(fit, participant) %>% 
  summarize(negLL=mean(negLL), alpha=mean(alpha), forget=mean(forget), rho=mean(rho), noise=mean(noise),
            neg_alpha=mean(neg_alpha), K=mean(K))
fit = filter(fit, participant %in% pes$participant)
fit = arrange(fit, participant)

subData = filter(subData, participant %in% fit$participant)

# Subject loop
subs <- unique(fit$participant)
nsubs <- length(subs) # number of subs
nblocks <- 12
ntrials_block <- matrix(NA, nrow = nsubs, ncol = nblocks)

# Define the softmax function
mcdougle_softmax_func <- function(Q, tau) {
  exp(Q * tau) / sum(exp(Q * tau))
}

sim_correct <- array(NA, dim = c(nsubs, nblocks, 1000))
sim_accu <- array(NA, dim = c(nsubs, nblocks, 1000))

subData = drop_na(subData)
# Subject loop
for (si in 1:nsubs) {
  idx_s <- subData$participant == subs[si]
  data <- subData[idx_s, ]
  
  # Set simulation parameters for each subject via the fit structure
  alpha <- fit$alpha[si]
  forget <- fit$forget[si]
  rho <- fit$rho[si]
  noise <- fit$noise[si]
  neg_alpha <- fit$neg_alpha[si]
  K <- fit$K[si]
  beta <- 100
  na <- 3 # number of possible button responses
  
  # Loop over simulation iterations
  iterations <- 1000
  blocks <- max(data$qNum)
  
  sim_correct_block <- array(NA, dim = c(iterations, blocks, 100))
  sim_accu_block <-array(NA, dim = c(iterations, blocks, 100))
  
  for (k in 1:iterations) {
    # Loop over blocks
    for (b in 1:blocks) {
      # Set parameters from stim that participant saw
      idx <- data$qNum == b
      ns <- unique(data$nSet[idx])
      stim <- data$stimNum[idx]
      sub_action <- data$key_press[idx]
      corresp <- data$corkey[idx]
      
      ntrials <- sum(idx)
      ntrials_block[si, b] <- ntrials
      
      # Model initialization
      q_rl <- matrix(1 / na, nrow = ns, ncol = na)
      q_wm <- matrix(1 / na, nrow = ns, ncol = na)
      p <- rep(NA, ntrials)
      weight <- rho * min(1, K / ns)
      
      # Trial loop
      for (t in 1:ntrials) {
        # Modular policies
        p_rl <- mcdougle_softmax_func(q_rl[stim[t], ], beta)
        p_wm <- mcdougle_softmax_func(q_wm[stim[t], ], beta)
        
        # Weighted policy
        pi <- (1 - weight) * p_rl + weight * p_wm
        pi <- noise * (1 / na) + (1 - noise) * pi
        #pi[is.na(pi)]<- 0
        
        # Now select an action
        x <- runif(1)
        sim_action <- which.max(cumsum(pi) >= x)
        
        # Determine if simulation correct
        reward <- ifelse(sim_action == corresp[t], 1, 0)
        
        # Q-learning
        if (reward == 1) { # pos rpe
          q_rl[stim[t], sim_action] <- q_rl[stim[t], sim_action] + alpha * (reward - q_rl[stim[t], sim_action])
          q_wm[stim[t], sim_action] <- q_wm[stim[t], sim_action] + 1 * (reward - q_wm[stim[t], sim_action])
        } else { # neg rpe
          q_rl[stim[t], sim_action] <- q_rl[stim[t], sim_action] + neg_alpha * alpha * (reward - q_rl[stim[t], sim_action])
          q_wm[stim[t], sim_action] <- q_wm[stim[t], sim_action] + neg_alpha * 1 * (reward - q_wm[stim[t], sim_action])
        }
        
        if (!is.na(sub_action[t])) {
          if (sim_action==sub_action[t]) {
            sim_corr = 1
          } else {
            sim_corr = 0
          }
        }
        
        
        # Forgetting
        q_wm <- q_wm + forget * ((1 / na) - q_wm)
        
        # Store simulation learning curves
        sim_correct_block[k, b, t] <- reward
        sim_accu_block[k, b, t] <- sim_corr
      }
    }
  }
  
  # Plotting
  par(mfrow = c(4, 3))
  par(mar=c(1,1,1,1))

  for (i in 1:blocks) {
    sim_correct[si, i, 1:ntrials_block[si, i]] <- colMeans(sim_correct_block[, i, 1:ntrials_block[si, i]], na.rm = TRUE)
    sim_accu[si, i, 1:ntrials_block[si, i]] <- colMeans(sim_accu_block[, i, 1:ntrials_block[si, i]], na.rm = TRUE)
    plot(1:ntrials_block[si, i], sim_correct[si, i, 1:ntrials_block[si, i]], type = "l",
         main = sprintf("subject: %.3d; block %.2d", si, i))
  }
}
#save(sim_correct, file = "sim_correct.RData")

# Tack simulation results onto the other data
#fitData_wSim <- read.csv('scarcityRL_v1_cleanedData_n32.csv')
#load('sim_correct.RData')

subData_wSim = subData
# Get some numbers for the loops
subs <- unique(subData_wSim$participant)
nsubs <- length(subs) # number of subs
nRL <- max(subData_wSim$qNum)

subData_wSim$simResults <- rep(NA, nrow(subData_wSim))
subData_wSim$simAccu<- rep(NA, nrow(subData_wSim))

# Loop over subjects and blocks
for (si in 1:nsubs) {
  for (bi in 1:nRL) {
    idx <- which(subData_wSim$participant == subs[si] & subData_wSim$qNum == bi)
    ntrials <- length(idx)
    
    simResults_block <- sim_correct[si, bi, 1:ntrials]
    simAccu <- sim_accu[si, bi, 1:ntrials]
    subData_wSim$simResults[idx] <- simResults_block
    subData_wSim$simAccu[idx]<- simAccu
  }
}
#save(fitData_wSim, file = "fitData_wSim.RData")

pes_rl = read.csv('pes_rl.csv')
pes_rl = filter(pes_rl, checks==3, N_skip<3, P_skip<3, D_skip<3)
acc = group_by(pes_rl, participant) %>% summarize(acc=mean(accuracy))
e = filter(acc, acc<=0.33 | acc>1)
e = e$participant
pes_rl = filter(pes_rl, participant %notin% e)

simAccuracy=group_by(subData_wSim, participant) %>% summarize(pred_acc=mean(simAccu, na.rm=TRUE))
mean(simAccuracy$pred_acc, na.rm = TRUE)
ggplot(simAccuracy, aes(x=pred_acc)) +geom_histogram()

subData_wSim$nSet = as.factor(subData_wSim$nSet)
ggplot(subData_wSim, aes(x=trialN, y=simResults, color=nSet)) + geom_smooth()
mean(subData_wSim$simResults, na.rm=TRUE)

subData_wSim$resp_corr = ifelse(subData_wSim$key_press==subData_wSim$corkey, 1, 0)

# Plotting
stimIter <- 15
nSetTwo <- matrix(NA, nrow = nsubs, ncol = stimIter)
nSetFour <- matrix(NA, nrow = nsubs, ncol = stimIter)

# Loop over subjects
for (si in 1:nsubs) {
  subjData <- subData_wSim[subData_wSim$participant == subs[si], ]
  nTwo_subj <- subjData[subjData$nSet == 2, ]
  nFour_subj <- subjData[subjData$nSet == 4, ]
  # Loop over iterations
  for (ii in 1:stimIter) {
    nSetTwo[si, ii] <- mean(nTwo_subj$simResults[nTwo_subj$stimIter == ii], na.rm = TRUE)
    nSetFour[si, ii] <- mean(nFour_subj$simResults[nFour_subj$stimIter == ii], na.rm = TRUE)
  }
}

nSet2 <- matrix(NA, nrow = nsubs, ncol = stimIter)
nSet4 <- matrix(NA, nrow = nsubs, ncol = stimIter)

for (si in 1:nsubs) {
  subjData <- subData_wSim[subData_wSim$participant == subs[si], ]
  nTwo_subj <- subjData[subjData$nSet == 2, ]
  nFour_subj <- subjData[subjData$nSet == 4, ]
  
  # Loop over iterations
  for (ii in 1:stimIter) {
    nSet2[si, ii] <- mean(nTwo_subj$resp_corr[nTwo_subj$stimIter == ii], na.rm = TRUE)
    nSet4[si, ii] <- mean(nFour_subj$resp_corr[nFour_subj$stimIter == ii], na.rm = TRUE)
  }
}

# Plot results
plot(1:15, colMeans(nSet2, na.rm = TRUE), type = "l", ylim = c(0, 1), xlab = "Iteration", ylab = "Accuracy")
lines(1:15, colMeans(nSet4, na.rm = TRUE), col = "purple")

lines(1:15, colMeans(nSetTwo, na.rm = TRUE), col='black', lty=5)
lines(1:15, colMeans(nSetFour, na.rm = TRUE), col = "purple", lty=5)
legend("bottomright", legend = c("Set Size 2", "Set Size 4"), col = c("black", "purple"), lty = 1)


sim_sum = ungroup(subData_wSim) %>% select(participant, nSet, stimIter, simResults, resp_corr) %>%
  gather(key=simResults, resp_corr, 4:5)
colnames(sim_sum)[4:5]=c('type', 'correct')
sim_sum$type = ifelse(sim_sum$type=='resp_corr', 'response', 'simulated')

ggplot(sim_sum, aes(x=stimIter, y=correct, color=nSet, linetype=type)) + geom_smooth() + 
  labs(color='Load', linetype='Type') + xlab('Stimulus iteration') + 
  theme(text=element_text(size = 19))
ggsave('sim_results.png')
