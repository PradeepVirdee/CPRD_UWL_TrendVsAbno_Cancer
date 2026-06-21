
install.packages("survival")
install.packages("joineRML")
install.packages("pROC")
library("joineRML")
library("pROC")
library("survival")


# Read in haemoglobin data for 10-year trend
data = read.csv("R/Data/DataForAnalysis_CohortUWLTrends_haemoglobin10y.csv")


# Keep co-occurring/latest test for abnormality model
datastatic = data[data$cooccurtest == 1,]


# Keep those with 2+ tests for JM
data = data[data$numtests > 1,]


# Linear splines (mixed model)
data$time0 = (data$time - 9.5)*I(data$time > 9.5)
data$time1 = (data$time - 9)*I(data$time > 9)
data$time2 = (data$time - 8)*I(data$time > 8)
data$time3 = (data$time - 7)*I(data$time > 7)
data$age1 = (data$age - 40)*I(data$age > 40)
data$age2 = (data$age - 80)*I(data$age > 80)


# Fit joint model
set.seed(123456)
jointmodel = mjoint(
  formLongFixed = list("hb" = haemoglobin ~ time*age + time*age1 + time*age2 + 
                         time0*age + time0*age1 + time0*age2 + 
                         time1*age + time1*age1 + time1*age2 + 
                         time2*age + time2*age1 + time2*age2 + 
                         time3*age + time3*age1 + time3*age2 + 
                         sex*age + sex*age1 + sex*age2),
  formLongRandom = list("hb" = ~ time | e_patid),
  formSurv = Surv(futime, cancer) ~ age + sex,
  data = list(data),
  timeVar = "time",
  control = list(covariance = "unstructured"),
)
summary(jointmodel)


# Run a loop that predicts the survival for each patient
datasurvpred = subset(data, select = e_patid)
datasurvpred = unique(datasurvpred)
datasurvpred$surv = NA
list = c(datasurvpred$e_patid)
for(i in list) {
  tryCatch({
    datapt = droplevels(data[data$e_patid == i, ])
    sfit = dynSurv(jointmodel, datapt, type = "first-order", u = 10.5)
    datasurvpred$surv[datasurvpred$e_patid == i] = sfit[["pred"]][["surv"]]
    datasurvpred$blup1int[datasurvpred$e_patid == i] = sfit[["b.hat"]][1]
    datasurvpred$blup1slope[datasurvpred$e_patid == i] = sfit[["b.hat"]][2]
  },
  error = function(e) {
    cat(i, "ERROR :", conditionMessage(e), "\n")
  })
}


# AUC - JM for trend
temp = data[!duplicated(data$e_patid), ]
dataforroc = merge(datasurvpred, temp, by = "e_patid")
dataforroc$risk = 1 - dataforroc$surv
roc(data = dataforroc, response = cancer, predictor = risk, auc = TRUE, ci = TRUE)


# AUC - Cox model for abnormality
datacox = datastatic
datacox$abnorm[datacox$haemoglobin < 130 & datacox$sex == "Male"] = 1
datacox$abnorm[is.na(datacox$abnorm) & datacox$haemoglobin < 115 & datacox$sex == "Female"] = 1
datacox$abnorm[is.na(datacox$abnorm)] = 0
coxfit = coxph(Surv(futime, cancer) ~ abnorm + age + sex, data = datacox)
datacox$lp = predict(coxfit, type = "lp", se.fit = FALSE)
bh = basehaz(coxfit)
bh$time2 = as.character(bh$time)
bh$time2 = substr(bh$time2, start=1, stop=4)
rows = which(grepl(10.5, bh$time2))
bh$time2 = NULL
num = rows[1] - 1
baseS = exp(-(bh$hazard[num]))
datacox$surv = baseS^(exp(datacox$lp))
datacox$risk = 1 - datacox$surv
roc(data = datacox, response = cancer, predictor = risk, auc = TRUE, ci = TRUE)


# Repeat entire code for each other blood test and trend window, changing code where relevant



