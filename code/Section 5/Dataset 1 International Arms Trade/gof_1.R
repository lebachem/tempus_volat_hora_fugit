# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ------                         Survey Paper                           ------ #
# ------                      (Fitting a (S)TERGM)                      ------ #
# ------         Data Set 1: International Arms Trade Network           ------ #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

# Preliminaries ----
library(statnet)
library(texreg)
library(PRROC)

rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
load("load_data.RData")
load("stergm_tergm_1.RData")
load("rem_1.Rdata")
source("../help_functions.R")

#TERGM
pdf("gof_tergm.pdf")
layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE))
my_gof(tergm_gof, plotlogodds = T, new_names = c("Rep","Edges","Reci","GWID","GWOD",
                                                 "GWESP","Polity Score","log(GDP) S","log(GDP) R"))
dev.off()

# STERGM
pdf("gof_stergm.pdf")
layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE))
my_gof(stergm_gof$formation, plotlogodds = T, new_names = c("Edges","Reci","GWID","GWOD",
                                                 "GWESP","Polity Score","log(GDP) S","log(GDP) R"))
my_gof(stergm_gof$dissolution, plotlogodds = T, new_names = c("Edges","Reci","GWID","GWOD",
                                                            "GWESP","Polity Score","log(GDP) S","log(GDP) R"))
dev.off()


# PR and ROC ----

# TERGM

# First get the change statistics 
change_info_tergm = ergmMPLE(net~edgecov(lnet)+
           edges+
           edgecov(lrecip)+
           gwidegree(log(2),fixed=T)+
           gwodegree(log(2),fixed=T)+
           gwesp(log(2),fixed=T)+
           edgecov(poldiff)+
           nodeocov("lgdp")+
           nodeicov("lgdp"),output = "array")
# Then construct the linear predictor for calculating the probability of 
# ties conditional on the rest of the network 
eta = lapply(1:length(tergm$coef), function(x){tergm$coef[x] * change_info_tergm$predictor[,,x]})
eta = Reduce("+",eta)
eta = as.numeric(eta)
eta = eta[-is.na(eta)]
P_tergm = exp(eta)/(1 + exp(eta))


# STERGM
# Save a matrix of the network at t where the diagonals are NA not 0 
tmp_mat = (as.matrix(nets[[1]]))
diag(tmp_mat) = NA
# Use this matrix to indicate which ties are governed 
# by the formation or dissolution model
formation_network = which(tmp_mat == 0, arr.ind = T)
dissolution_network = which(tmp_mat == 1, arr.ind = T)

# The covaraites slightly differ -> bew change statistics 
change_info_stergm = ergmMPLE(net~ edges+
                                mutual+
                                gwidegree(log(2),fixed=T)+
                                gwodegree(log(2),fixed=T)+
                                gwesp(log(2),fixed=T)+
                                edgecov(poldiff)+
                                nodeocov("lgdp")+
                                nodeicov("lgdp"),output = "array")


# This part is unchanged
eta_formation = lapply(1:length(stergm$formation.fit$coef), 
                       function(x){stergm$formation.fit$coef[x] * 
                           change_info_stergm$predictor[cbind(formation_network,x)]})
eta_formation = Reduce("+",eta_formation)
eta_dissolution = lapply(1:length(stergm$dissolution.fit$coef), 
                         function(x){stergm$dissolution.fit$coef[x] * 
                             change_info_stergm$predictor[cbind(dissolution_network,x)]})
eta_dissolution = Reduce("+",eta_dissolution)
eta_stergm = c(eta_formation, eta_dissolution)
P_stergm = exp(eta_stergm)/(1 + exp(eta_stergm))
real_stergm = c(stergm$formation.fit$network[formation_network],
                stergm$dissolution.fit$network[dissolution_network])


# REM

Lambda_rem = as.numeric(exp(coefs[1]*rem_arms_trade_data$dep[[1]]$change.contributions[,,1]  + 
                           coefs[2]*rem_arms_trade_data$dep[[1]]$change.contributions[,,2] +
                           coefs[3]*rem_arms_trade_data$dep[[1]]$change.contributions[,,3] +
                           coefs[4]*rem_arms_trade_data$dep[[1]]$change.contributions[,,4] +
                           coefs[5]*rem_arms_trade_data$dep[[1]]$change.contributions[,,5] + 
                           coefs[6]*rem_arms_trade_data$dep[[1]]$change.contributions[,,6] +
                           coefs[7] *rem_arms_trade_data$dep[[1]]$change.contributions[,,7] +
                             coefs[8] *rem_arms_trade_data$dep[[1]]$change.contributions[,,8]))


# Calculate the ROC values
roc_tergm = roc_my(P_tergm,networks[[2]])
roc_stergm = roc_my(P_stergm,real_stergm)
roc_rem = roc.curve(scores.class0 = Lambda_rem,weights.class0 = networks[[2]],curve=T)

df_roc_tergm = data.frame("specificity" =1- roc_tergm$curve[,1],
                          "sensitivity" = roc_tergm$curve[,2], "Model" = "TERGM")
df_roc_stergm = data.frame("specificity" =1- roc_stergm$curve[,1],
                           "sensitivity" = roc_stergm$curve[,2], "Model" = "STERGM")
df_roc_rem = data.frame("specificity" =1- roc_rem$curve[,1],
                        "sensitivity" = roc_rem$curve[,2], "Model" = "REM")
df_roc= rbind(df_roc_tergm, df_roc_stergm, df_roc_rem)

pdf("roc_comparison.pdf",width = 9,height =9)

ggplot(data = df_roc ,mapping = aes(x = specificity,y = sensitivity, lty = Model)) +
  geom_step(alpha = 0.8, size=1) +
  scale_x_reverse() +
  geom_abline(slope = 1, intercept = 1, lty = 2, color = "grey40") +
  ylab("Sensitivity")+
  xlab("Specificity") +
  theme_pubr() +
  scale_linetype_manual(name = "", 
                     labels = c(paste0("TERGM (",round(roc_tergm$auc,digits = 3) ,")"),
                                paste0("STERGM (",round(roc_stergm$auc,digits = 3) ,")"),
                                paste0("REM (",round(roc_rem$auc,digits = 3) ,")")), 
                     values = c(1,4,2),guide = guide_legend(keywidth = 4)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.x = element_text(size = 30)) +
  theme(axis.title.y = element_text(size = 30)) +
  theme(axis.text.x = element_text(size = 20))  +
  theme(axis.text.y = element_text(size = 20)) +
  theme(plot.title = element_text(size = 30), 
        legend.text = element_text(size = 15) )

dev.off()

# Calculate the pr values
pr_tergm = pr_my(P_tergm,networks[[2]])
pr_stergm = pr_my(P_stergm,real_stergm)
pr_rem = pr.curve(scores.class0 = Lambda_rem,weights.class0 = networks[[2]],curve=T)

df_pr_tergm = data.frame("Precision" =pr_tergm$curve[,1],
                         "Recall" = pr_tergm$curve[,2], "Model" = "TERGM")
df_pr_stergm = data.frame("Precision" = pr_stergm$curve[,1],
                          "Recall" = pr_stergm$curve[,2], "Model" = "STERGM")
df_pr_rem = data.frame("Precision" =pr_rem$curve[,1],
                       "Recall" = pr_rem$curve[,2], "Model" = "REM")
df_pr= rbind(df_pr_tergm, df_pr_stergm, df_pr_rem)

pdf("pr_comparison.pdf",width = 9,height =9)
ggplot(data = df_pr ,mapping = aes(x = Precision,y = Recall, lty = Model)) +
  geom_step(alpha = 0.8, size=1) +
  ylab("Precision")+
  xlab("Recall") +
  theme_pubr() +
  scale_linetype_manual(name = "", 
                        labels = c(paste0("TERGM (",round(pr_tergm$auc.integral,digits = 3) ,")"),
                                   paste0("STERGM (",round(pr_stergm$auc.integral,digits = 3) ,")"),
                                   paste0("REM (",round(pr_rem$auc.integral,digits = 3) ,")")), 
                        values = c(1,4,2),guide = guide_legend(keywidth = 3)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.title.x = element_text(size = 30)) +
  theme(axis.title.y = element_text(size = 30)) +
  theme(axis.text.x = element_text(size = 20))  +
  theme(axis.text.y = element_text(size = 20)) +
  theme(plot.title = element_text(size = 30), 
      legend.text = element_text(size = 15) )

dev.off()
