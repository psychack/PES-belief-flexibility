library(tidyverse)
library(robustbase)
library(lmerTest)
library(lme4)
library(effectsize)

pes_rl = read.csv('pes_rlwm.csv')
pes_bade = read.csv('pes_belief.csv')
fitData = read.csv('rlwm_params.csv')
pes = read.csv('pes_data.csv')

# effect of positive schizotypy on accuracy
lmrob(accuracy ~ load + z_P + z_N + z_D + z_P:load + z_N:load + z_D:load, pes_rl) %>% summary()
lmrob(accuracy ~ load + z_P + z_N + z_D + z_P:load + z_N:load + z_D:load, pes_rl) %>% effectsize()

# effect of positive schizotypy on PES
lmer(PES ~ load + z_P + z_N + z_D + z_reward + z_delay + z_reward:z_delay + z_reward:load + 
       z_delay:load + z_P:z_reward + z_N:z_reward + z_D:z_reward + z_P:load + z_N:load + 
       z_D:load + z_P:z_delay + z_N:z_delay + z_D:z_delay + (1|participant), pes) %>% summary()

# positive schizotypy and RLWM parameters
lmrob(alpha ~ z_P + z_N + z_D, fitData) %>% summary() # no sig.

lmrob(forget ~ z_P + z_N + z_D, fitData) %>% summary() # higher WM decay in Pos SZ
lmrob(forget ~ z_P + z_N + z_D, fitData) %>% effectsize()

lmrob(persever ~ z_P + z_N + z_D, fitData) %>% summary() # less receptive to negative feedback in Pos SZ
lmrob(persever ~ z_P + z_N + z_D, fitData) %>% effectsize()

lmrob(rho ~ z_P + z_N + z_D, fitData) %>% summary()

lmrob(noise ~ z_P + z_N + z_D, fitData) %>% summary() # higher noise in Dis SZ
lmrob(noise ~ z_P + z_N + z_D, fitData) %>% effectsize()

lmrob(K ~ z_P + z_N + z_D, fitData) %>% summary() # lower WM capacity in Pos SZ
lmrob(K ~ z_P + z_N + z_D, fitData) %>% effectsize()

# PES and belief flexibility
lmrob(flexi ~ PES_trad + z_P + z_N + z_D + PES_trad:z_P + PES_trad:z_N + PES_trad:z_D, pes_bade) %>% 
  summary()
lmrob(flexi ~ PES_trad + z_P + z_N + z_D + PES_trad:z_P + PES_trad:z_N + PES_trad:z_D, pes_bade) %>% 
  effectsize()

# RLWM parameters and PES
cor.test(pes_rl$forget, pes_rl$PES_trad, method='spearman') #*** #p<0.001 with robust PES
cor.test(pes_rl$alpha, pes_rl$PES_trad, method='spearman') #***
cor.test(pes_rl$persever, pes_rl$PES_trad, method='spearman') #*** # p=0.01 with robust PES
cor.test(pes_rl$rho, pes_rl$PES_trad, method='spearman') #*
cor.test(pes_rl$noise, pes_rl$PES_trad, method='spearman')
cor.test(pes_rl$K, pes_rl$PES_trad, method='spearman') #**
