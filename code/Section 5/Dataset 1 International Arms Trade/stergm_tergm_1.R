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

# mat to vec
mat_to_vec=function(mat){
  resp=as.matrix(mat)
  diag(resp)=NA
  resp=c(resp)
  return(resp[is.na(resp)==F])
}
# set

# ROC for matrix predictions
roc_pr<-function(Pred,Response,ret_val=T){

  p_hat=mat_to_vec(Pred)
  response=mat_to_vec(Pred)
  
  roc<-roc.curve(scores.class0 = p_hat,weights.class0 = response,curve=T)
  pr<-pr.curve(scores.class0 = p_hat,weights.class0 = response,curve=T)

  plot(roc,color=F)
  plot(pr,color=F)
  
  if (ret_val==T){
    return(list(roc<-roc$auc,pr<-pr$auc.integral))
  }
}


set.seed(123)


# Start and End
beginning = 2016
ending = 2017

# Save the subgraph of those networks 
networks = list()
countries_ind =which(rowSums(EX[,beginning:ending -1949]) == 2)

# Discretizing
networks[[1]] = amk[[beginning - 1949]][countries_ind, countries_ind]>0
networks[[2]] = amk[[ending - 1949]][countries_ind, countries_ind]>0

# Response network ----
net<-network(networks[[2]])


# Covariates----

# log GDP
lgdp<-log(real_gdp[countries_ind,67])

# Absolute difference polity score
poldiff<-outer(polity[countries_ind,67],polity[countries_ind,67],FUN=function(x,y) abs(x-y))

# Lagged network
lnet<-network(networks[[1]])
lrecip<-network(t(networks[[1]]))
# set attributes
net <- set.vertex.attribute(net, "lgdp", lgdp)
lgdp<-log(real_gdp[countries_ind,66])
lnet <- set.vertex.attribute(lnet, "lgdp", lgdp)

# Combine them in a list
nets<-network.list(list(lnet,net))

# Estimation ----
nsim=5000
# TERGM

tergm<-ergm(net~edgecov(lnet)+
              edges+
              edgecov(lrecip)+
              gwidegree(log(2),fixed=T)+
              gwodegree(log(2),fixed=T)+
              gwesp(log(2),fixed=T)+
              edgecov(poldiff)+
              nodeocov("lgdp")+
              nodeicov("lgdp"),control = control.ergm(MCMC.samplesize = nsim))

# STERGM ----

stergm<-stergm(nets,formation = ~edges+
                 mutual+
                 gwidegree(log(2),fixed=T)+
                 gwodegree(log(2),fixed=T)+
                 gwesp(log(2),fixed=T)+
                 edgecov(poldiff)+
                 nodeocov("lgdp")+
                 nodeicov("lgdp"),
               dissolution  = ~edges+
                 mutual+
                 gwidegree(log(2),fixed=T)+
                 gwodegree(log(2),fixed=T)+
                 gwesp(log(2),fixed=T)+
                 edgecov(poldiff)+
                 nodeocov("lgdp")+
                 nodeicov("lgdp"),times=1:2,estimate = "CMLE",control = control.stergm(CMLE.MCMC.interval = nsim))



# PR and ROC ----

# TERGM

change_we_believe_in<-ergmMPLE(net~edgecov(lnet)+
  edges+
  edgecov(lrecip)+
  gwidegree(log(2),fixed=T)+
  gwodegree(log(2),fixed=T)+
  gwesp(log(2),fixed=T)+
  edgecov(poldiff)+
  nodeocov("lgdp")+
  nodeicov("lgdp"))
logit=change_we_believe_in$predictor%*%tergm$coef
p_hat=1/(1+exp(-logit))

roc<-roc.curve(scores.class0 = p_hat,weights.class0 = change_we_believe_in$response,curve=T)
pr<-pr.curve(scores.class0 = p_hat,weights.class0 = change_we_believe_in$response,curve=T)
plot(roc)
plot(pr)

sim_tergm<-simulate(tergm,nsim=nsim,nw.start=lnet)
P_tergm=Reduce("+",lapply(sim_tergm,as.matrix))/nsim

pdf("roc_pr_tergm.pdf")
par(mfrow=c(1,2))
roc_pr(P_tergm,networks[[2]])
dev.off()

# STERGM
sim_stergm<-simulate(stergm,nw.start=lnet,nsim=nsim)
P_stergm=Reduce("+",lapply(sim_stergm,as.matrix))/nsim

pdf("roc_pr_stergm.pdf")
par(mfrow=c(1,2))
roc_pr(P_stergm,networks[[2]])
dev.off()

# GOF statistics ----

#TERGM
tergm_gof<-gof(tergm,verbose=T)
summary(tergm_gof)

pdf("gof_tergm.pdf")
layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE))
plot(tergm_gof, plotlogodds=T)
dev.off()

# STERGM
stergm_gof<-gof(stergm)
summary(stergm_gof)

pdf("gof_stergm.pdf")
layout(matrix(c(1,2,3,4,5,5), 3, 2, byrow = TRUE))
plot(stergm_gof, plotlogodds=T)
dev.off()

# Save the results
save.image("stergm_tergm.RData")


# Create a Table ----
theta_stergm<-coef(stergm)
theta_tergm<-coef(tergm)
Mat<-cbind(theta_tergm,c(0,theta_stergm$formation),c(0,theta_stergm$dissolution))

# Save the P-values 
p_val<-list()
p_val[[2]]<-c(0,summary(stergm)$formation$coefs[,4])
p_val[[3]]<-c(0,summary(stergm)$dissolution$coefs[,4])
p_val[[1]]<-summary(tergm)$coefs[,4]

# Save GOF-Stats----
lgof<-list()
lgof[[2]] <- c(summary(stergm)$formation$aic,summary(stergm)$formation$mle.lik)
lgof[[3]] <- c(summary(stergm)$dissolution$aic,summary(stergm)$dissolution$mle.lik)
lgof[[1]] <- c(AIC(tergm),logLik(tergm))
gof.names <- c("AIC", "Log Likelihood") #names of GOFs

est <- Mat
se <- cbind(summary(tergm)$coefs[,2],c(0,summary(stergm)$formation$coefs[,2]),c(0,summary(stergm)$dissolution$coefs[,2]))
colnames(est) <- c("TERGM","STERGM, Formation","STERGM, Dissolution")

# add row labels:
rownames(est) <- c("lagged Network",rownames(Mat)[-1])

# create a texreg object:
tr <- list()
for (j in 1:ncol(est)) {
  tr[[j]] <- createTexreg(
    coef.names = rownames(est), 
    coef = est[, j], 
    se = se[, j],
    pvalues = p_val[[j]],
    ci.low = numeric(0),
    ci.up = numeric(0),
    gof.names <- c("AIC",  "Log Likelihood"),
    gof =lgof[[j]]
    
  )
}

# for text output:
screenreg(tr, custom.model.names = colnames(est), 
          custom.note = "",stars = c(0.001,  0.01, 0.05),gof=gof,digits = 3)

# for LaTeX output:
texreg(tr, custom.model.names = colnames(est),
       custom.note = "",digits=3,stars = c(0.001,  0.01, 0.05))

