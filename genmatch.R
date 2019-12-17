library(Matching)
library(fastDummies)
rm(list=ls())
load("prematch_data_all.RData")
data <- data.prematch

# Set treatment groups
data$TREATV <- NA
data[data$VOLATQ==5,"TREATV"] <- 1
data[data$VOLATQ==1,"TREATV"] <- 0
data$TREATB <- NA
data[data$BETAQ==5,"TREATB"] <- 1
data[data$BETAQ==1,"TREATB"] <- 0

# Take log of cap and volume
data$LOGCAP <- log(data$CAP)
data$LOGVOL <- log(data$VOL)

# Remove rows with missing data. Make volatility related data
data.vol <- data[,c("TREATV","PERIOD","RETURN","CAP","LOGCAP","VOL","LOGVOL","GROUP","INDUSTRY")]
data.vol <- data.vol[complete.cases(data.vol),]
data.vol <- data.vol[data.vol$CAP>0,]
data.vol <- data.vol[data.vol$VOL>0,]

group_coarse <- rep(0, nrow(data.vol))
group_col <- data.vol$GROUP
for (i in (1:length(group_col))) {
  if (group_col[i] >= 1 & group_col[i] <= 9) {
    group_coarse[i] <- 1
  }
  if (group_col[i] >= 10 & group_col[i] <= 14) {
    group_coarse[i] <- 2
  }
  if (group_col[i] >= 15 & group_col[i] <= 17) {
    group_coarse[i] <- 3
  }
  if (group_col[i] >= 20 & group_col[i] <= 39) {
    group_coarse[i] <- 4
  }
  if (group_col[i] >= 40 & group_col[i] <= 49) {
    group_coarse[i] <- 5
  }
  if (group_col[i] == 50 | group_col[i] == 51) {
    group_coarse[i] <- 6
  }
  if (group_col[i] >= 52 & group_col[i] <= 59) {
    group_coarse[i] <- 7
  }
  if (group_col[i] >= 60 & group_col[i] <= 67) {
    group_coarse[i] <- 8
  }
  if (group_col[i] >= 70 & group_col[i] <= 89) {
    group_coarse[i] <- 9
  }
  if (group_col[i] >= 90 & group_col[i] <= 99) {
    group_coarse[i] <- 10
  }
}

data.vol$GROUP_COARSE <- group_coarse
data.vol <- data.vol[data.vol$GROUP_COARSE != 0,]
data.vol <- dummy_cols(data.vol, select_columns = 'GROUP_COARSE')

# Make beta related data
data.beta <- data[,c("TREATB","PERIOD","RETURN","CAP","LOGCAP","VOL","LOGVOL","GROUP","INDUSTRY")]
data.beta <- data.beta[complete.cases(data.beta),]
data.beta <- data.beta[data.beta$CAP>0,]
data.beta <- data.beta[data.beta$VOL>0,]

start <- 61
end <- 552
data.prunedv <- data.frame()
data.prunedb <- data.frame()


data.curv <- data.vol[data.vol$PERIOD==62,]
Xv <- cbind(data.curv$LOGCAP, data.curv$LOGVOL, data.curv$GROUP_COARSE_5,
            data.curv$GROUP_COARSE_8, data.curv$GROUP_COARSE_2, data.curv$GROUP_COARSE_4,
            data.curv$GROUP_COARSE_9, data.curv$GROUP_COARSE_6, data.curv$GROUP_COARSE_7,
            data.curv$GROUP_COARSE_3, data.curv$GROUP_COARSE_10)
genoutv <- GenMatch(Tr=data.curv$TREATV, X=Xv, estimand="ATT", M=1, pop.size=500,
                    max.generations=20, wait.generations=10, caliper = 0.1, print.level = 0)
moutv <- Match(Tr=data.curv$TREATV, X=Xv, estimand="ATT", Weight.matrix=genoutv, caliper = 0.1)
mb <- MatchBalance(TREATV ~ LOGCAP + LOGVOL + GROUP_COARSE_5 + GROUP_COARSE_8+ GROUP_COARSE_2 + GROUP_COARSE_4+ GROUP_COARSE_9+ GROUP_COARSE_6+ GROUP_COARSE_7+ GROUP_COARSE_3+ GROUP_COARSE_10,
             data = data.curv, match.out = moutv, nboots = 10000, print.level = 0)
mb$AMsmallest.p.value
summary(moutv)
moutv



data.curb <- data.beta[data.beta$PERIOD==61,]
Xb <- cbind(data.curb$LOGCAP, data.curb$LOGVOL, data.curb$GROUP)
genoutb <- GenMatch(Tr=data.curb$TREATB, X=Xb, estimand="ATT", M=1, pop.size=100,
                    max.generations=10, wait.generations=5, exact = c(F, F, T))
moutb <- Match(Tr=data.curb$TREATB, X=Xb, estimand="ATT", Weight.matrix=genoutb)
##### works up to here (well not really because the match balance is horrible) #####


outputv <- match.data(moutv)
outputb <- match.data(genoutb)
data.prunedv <- rbind(data.prunedv,moutv)
data.prunedb <- rbind(data.prunedb,moutb)
