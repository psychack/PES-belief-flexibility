library(nloptr)
library(tidyverse)

pes = read.csv('pes.csv')
pes_rl = read.csv('pes_rl.csv')

#`%notin%` = Negate(`%in%`)
include = pes_rl$participant
pes = filter(pes, participant %in% include)

# preparing data
subs = unique(pes$participant)
nsubs = length(subs)
fitData <- matrix(NA, nrow=nsubs*250, ncol=8)
fitData = data.frame(fitData)
colnames(fitData) = c('participant','alpha','forget', 'rho','noise','neg_alpha', 'K','negLL')
pes$qNum = pes$block
#fitData = read.csv('fitData.csv')

lower_bounds <- c(0, 0, 0, 0, 0, 1)
upper_bounds <- c(1, 1, 1, 1, 1, 5)

# data to fit RLWM
pes_fit = select(pes, participant, rl_resp.keys, rl_resp.corr, load, qNum, rl_trials.thisIndex)
pes_fit$key_press = ifelse(pes_fit$rl_resp.keys=='j', 1, ifelse(pes_fit$rl_resp.keys=='k', 2, 3))

colnames(pes_fit)[3] = 'correct'
colnames(pes_fit)[4]='nSet'
colnames(pes_fit)[6]='stimNum'

pes_fit$stimNum = pes_fit$stimNum+1
pes_fit=select(pes_fit, -rl_resp.keys)

mcdougle_softmax_func <- function(Q, tau) {
  # Input: Q is estimated payoff values
  #        tau = inverse "temperature." High values make
  #        probs farther from each other (more greedy), low values less greedy
  # Output: The probabilities of choosing each action
  
  # Compute the exponentiated and scaled Q values
  exp_Q_tau <- exp(Q * tau)
  
  # Normalize to get probabilities
  Probs <- exp_Q_tau / sum(exp_Q_tau)
  
  return(Probs)
}

func_rlwm <- function(params, df, beta) {
  
  data <- df # data
  blocks <- max(data$qNum) # get number of blocks
  # inputted param values (K is input up top)
  alpha <- params[1] # keep track of ordering based on the "fit" script
  forget <- params[2]
  rho <- params[3]
  noise <- params[4]
  neg_alpha <- params[5]
  K<- params[6]
  
  na <- 3 # number of possible button responses
  
  negLL <- numeric(blocks)
  
  # loop over each block 
  for (b in 1:blocks) {
    
    # extract data for this block
    idx <- data$qNum == b
    ns <- unique(data$nSet[idx]) # set size of this block
    reward <- data$correct[idx] # rewarded? (1:ntrials)
    stim <- data$stimNum[idx] # stimuli identifiers (1:ntrials)
    sub_action <- data$key_press[idx] # subject's action (1:ntrials), subtract 73 to get 1,2,or3
    ntrials <- sum(idx) # sum of idx instead of max(blocktrial) incase trials were excluded
    
    # model initialization
    q_rl <- matrix(1 / na, nrow = ns, ncol = na) # q values of RL system (init at 0.33)
    q_wm <- matrix(1 / na, nrow = ns, ncol = na) # q values of WM system (init at 0.33)
    p <- numeric(ntrials) # policy (prob of each action each trial)
    weight <- rho * min(1, K / ns) # initial weighting of WM (implements rho and k)
    
    # trial loop
    for (t in 1:ntrials) {
      
      if (!is.na(sub_action[t])) { # valid trial?
        
        # modular policies
        p_rl <- mcdougle_softmax_func(q_rl[stim[t], ], beta) # put q_rl through softmax
        p_wm <- mcdougle_softmax_func(q_wm[stim[t], ], beta) # put q_wm through softmax
        
        # weighted policy
        pi <- (1 - weight) * p_rl + weight * p_wm # policy vector weighting RL and WM
        pi <- noise * (1 / na) + (1 - noise) * pi # add undirected noise
        p[t] <- pi[sub_action[t]] # get model's probability of action taken by subject
        
        # Q-learning
        if (reward[t] == 1) { # pos rpe
          q_rl[stim[t], sub_action[t]] <- q_rl[stim[t], sub_action[t]] + 
            alpha * (reward[t] - q_rl[stim[t], sub_action[t]])
          q_wm[stim[t], sub_action[t]] <- q_wm[stim[t], sub_action[t]] + 
            1 * (reward[t] - q_wm[stim[t], sub_action[t]]) # perfect learning rate for wm
        } else { # neg rpe
          q_rl[stim[t], sub_action[t]] <- q_rl[stim[t], sub_action[t]] + 
            neg_alpha * alpha * (reward[t] - q_rl[stim[t], sub_action[t]]) # less learning on non-rewarded trials
          q_wm[stim[t], sub_action[t]] <- q_wm[stim[t], sub_action[t]] + 
            neg_alpha * 1 * (reward[t] - q_wm[stim[t], sub_action[t]])
        }
        
        # forgetting
        q_wm <- q_wm + forget * ((1 / na) - q_wm)
        
      }
    }
    
    epsilon <- 0.000001 # prevent overflow in likelihood function
    p <- epsilon / 2 + (1 - epsilon) * p
    
    negLL[b] <- -sum(log(p), na.rm = TRUE) # build up block log-likelihood
  }
  
  ll <- sum(negLL) # this is the output to be minimized (i.e., negative log likelihood)
  return(ll)
}

for (s in 1:nsubs) {
  id <- pes_fit$participant == subs[s]
  fit_data <- pes_fit[id, ]
    
  for (i in 1:300) {
    obj_func <- function(params) {
      df <- fit_data
      beta <- 100
      ll <- func_rlwm(params, df, beta)
      return(ll)
    }
    
    initial_params <- c(alpha = runif(1,0,1), forget = runif(1,0,1), rho = runif(1,0,1), 
                        noise = runif(1,0,1), neg_alpha = runif(1,0,1), K=runif(1, 1, 5))
    result <- nloptr(
      x0 = initial_params,
      eval_f = obj_func,
      lb = lower_bounds,
      ub = upper_bounds,
      opts = list(
        algorithm = "NLOPT_LN_BOBYQA",
        xtol_rel = 1e-10,
        maxeval = 1e+6
      )
    )
    fitData[(s-1)*300+i, 1]=subs[s]
    fitData[(s-1)*300+i, 2:7]=result$solution
    fitData[(s-1)*300+i, 8]=result$objective
    cat("Participant:", subs[s], "Iteration:", i, "NegLL:", result$objective, '\n')
  }
}
