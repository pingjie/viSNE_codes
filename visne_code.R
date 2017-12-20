#setwd("C:\\Users\\pingj2\\Dropbox\\Work\\Vanderbilt\\viSNE\\data")

library(tsne)
library(Hmisc)
library(anytime)
library(randomForest)
library(lubridate)
library(ggplot2)
library(rpart)

timeinternal <- function(dob, age.day = today(), units = "years", floor = TRUE) {
  calc.age = new_interval(dob, age.day) / duration(num = 1, units = units)
  if (floor)
    return(as.integer(floor(calc.age)))
  return(calc.age)
}

PCA_3 <- function(data, pch_point = 1, mar = c(4, 4, 4, 2), cex_point = 2, cex_name = 2, main = "PCA", reverse = T, 
                  scale = T, add_name = T, col_name = "red", col_point = "black" , barplot = F, 
                  pos = NULL, cex.axis = 1.2, cex.main = 1.2, cex.lab = 1.2, ...) {
  if (barplot == T) {
    op <- par(mfrow = c(2, 1), mar = .1 + c(2, 2, 2, 2))
  }
  if (reverse == T) {
    data <- t(data)
  }
  result <- prcomp(data, scale = scale)
  
  par(mar = mar)
  plot(result$x[, 1], result$x[, 2], type = "n",
       xlab = paste("PC1 (", round(((result$sdev)^2/sum((result$sdev)^2))[1], 2)*100, "%)", sep = ""),
       ylab = paste("PC2 (", round(((result$sdev)^2/sum((result$sdev)^2))[2], 2)*100, "%)", sep = ""),
       main  =main, cex.axis = cex.axis, cex.main = cex.main, cex.lab = cex.lab, ...)
  points(result$x[, 1], result$x[, 2], pch = pch_point, cex = cex_point, col = col_point)
  
  if (add_name == T) {
    text(result$x[, 1], result$x[, 2], labels = rownames(result$x), cex = cex_name, col = col_name, pos = pos)
  }    
  
  #abline(h=0,lty=3)
  #abline(v=0,lty=3)
  
  if (barplot == T) {
    barplot((result$sdev)^2/sum((result$sdev)^2), names = colnames(result$x), main = "Proportion of Variance")
  }
  
  return (result)
}


setwd("Z:/JiePing/viSNE/data") ### Research Drive
setwd("/Users/pingjie/Dropbox/Work/Vanderbilt/viSNE/data/") ### macOS setwd

all.input.data <- read.csv("2017-09-18baseline_covariates.csv", header = T, row.names = 1)
all.outcome.data <- read.csv("2017-09-18index_procedure_and_recurrence.csv", header = T, row.names = 1)

all.data <- merge(all.input.data, all.outcome.data, by = "row.names", all.x = T, all.y = T)

case_with_procedure <- all.data[which(!is.na(all.data$procedureType)), ]

case_with_procedure <- case_with_procedure[, -(14:16)] ## remove BIPQ.illnessCause1, BIPQ.illnessCause2, BIPQ.illnessCause3
case_with_procedure <- case_with_procedure[, -224]  ### remove MOS.closeFriends

age_at_procedure <- timeinternal(anydate(case_with_procedure$dob), anydate(case_with_procedure$indexProcedureDate))
procedure_to_lastfollowup_days <- timeinternal(anydate(case_with_procedure$indexProcedureDate), anydate(case_with_procedure$lastFollowUpDate), units = "day")

case_with_procedure <- cbind(case_with_procedure, age_at_procedure, procedure_to_lastfollowup_days)

q1data <- all.data[, grep(pattern = "BIPQ|CDQ|EAT10|FoPQ|MOS|SDMQ9|SF12v2|VHI10|procedureType|recurrence.*", x = colnames(all.data))]
q1data <- q1data[, -128] ### remove procedureTypeOther
q1data <- q1data[, -(9:11)] ## remove BIPQ.illnessCause1, BIPQ.illnessCause2, BIPQ.illnessCause3
q1data <- q1data[, -72]  ### remove MOS.closeFriends

## Summarized original data, correlation between features and procedure type or recurrence
summaryTable <- matrix(nrow = 122, ncol = 5)
rownames(summaryTable) <- colnames(q1data)[1:122]

for (n in 2:123) {
  summaryTable[n - 1, 1] <- summary(q1data[, n])["NA's"]
  
  lm_procedure <- lm(q1data[, n] ~ q1data$procedureType)
  
  summaryTable[n - 1, 2] <- summary(lm_procedure)$r.squared
  summaryTable[n - 1, 3] <- summary(lm_procedure)$adj.r.squared
  
  lm_recurrence <- lm(q1data[, n] ~ q1data$recurrence)
  
  summaryTable[n - 1, 4] <- summary(lm_recurrence)$r.squared
  summaryTable[n - 1, 5] <- summary(lm_recurrence)$adj.r.squared
}

colnames(summaryTable) <- c("No of NAs", "R squared Procedure", "adjusted R squared Procedure", "R squared recurrence", "adjusted R squared recurrence")
write.csv(summaryTable, file = "variable_summary_cor.csv")

#### tSNE with only original data

q1data.complete <- q1data[complete.cases(q1data), ]

tsne_q1_complete <- tsne(q1data.complete[, 1:122], perplexity = 50)

jpeg(filename = "tsne_q1_complete_procedure.jpg", width = 2400, height = 1600, units = "px", pointsize = 8, quality = 100, bg = "white", res = 300)
ggplot(as.data.frame(tsne_q1_complete), aes(x = V1, y = V2, color = q1data.complete$procedureType)) + 
  geom_point(cex = 2) + scale_colour_manual(name = "procedureType", values = rainbow(length(unique(q1data.complete$procedureType)))) +
  ggtitle("tSNE of procedureType from only original data") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

jpeg(filename = "tsne_q1_complete_recurrence.jpg", width = 2400, height = 1600, units = "px", pointsize = 8, quality = 100, bg = "white", res = 300)
ggplot(as.data.frame(tsne_q1_complete), aes(x = V1, y = V2, color = q1data.complete$recurrence)) + 
  geom_point(cex = 2) + scale_colour_manual(name = "recurrence", values = rainbow(length(unique(q1data.complete$recurrence)))) +
  ggtitle("tSNE of recurrence from only original data") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

### Two procedure types
q1_endo_trac_procedure_data <- subset(q1data, procedureType %in% c("Endoscopic procedure", "Tracheal resection"))
q1_endo_trac_procedure_data$procedureType <- droplevels(q1_endo_trac_procedure_data$procedureType)

q1endo.complete <- q1_endo_trac_procedure_data[complete.cases(q1_endo_trac_procedure_data), ]

tsne_q1 <- tsne(q1endo.complete[, 1:122], perplexity = 50)

jpeg(filename = "tsne_q1_complete_two_procedure.jpg", width = 2400, height = 1600, units = "px", pointsize = 8, quality = 100, bg = "white", res = 300)
ggplot(as.data.frame(tsne_q1), aes(x = V1, y = V2, color = q1endo.complete$procedureType)) + 
  geom_point(cex = 2) + scale_colour_manual(name = "procedureType", values = rainbow(length(unique(q1endo.complete$procedureType)))) +
  ggtitle("tSNE of two procedureType with only original data") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

####impute by Hmisc###

imputeNames <- colnames(case_with_procedure)[grep(pattern = "BIPQ|CDQ|EAT10|FoPQ|MOS|SDMQ9|SF12v2|VHI10.*", x = colnames(case_with_procedure))]

imputef_withoutage <- as.formula(paste("~", paste(imputeNames, collapse = "+")))
imputef_withage <- as.formula(paste("~", paste(imputeNames, collapse = "+"), "+ age_at_procedure + procedure_to_lastfollowup_days" ))

set.seed(1)

## without age and follow-up days
q1impute_withoutage <- aregImpute(imputef_withoutage, data = case_with_procedure, n.impute = 10, nk = 0)

q1impute_withoutage_n1 <- impute.transcan(q1impute_withoutage, imputation = 1, list.out = T, data = case_with_procedure)

q1impute_withoutage_n1 <- as.data.frame(q1impute_withoutage_n1)

rownames(q1impute_withoutage_n1) <- rownames(case_with_procedure)

q1impute_withoutage_n1 <- cbind(q1impute_withoutage_n1, case_with_procedure$procedureType)
colnames(q1impute_withoutage_n1)[123] <- "procedureType"

q1impute_withoutage_n1_endo_trac_procedure_data <- subset(q1impute_withoutage_n1, procedureType %in% c("Endoscopic procedure", "Tracheal resection"))
q1impute_withoutage_n1_endo_trac_procedure_data$procedureType <- droplevels(q1impute_withoutage_n1_endo_trac_procedure_data$procedureType)

tsne_q1impute_withoutage_n1 <- tsne(q1impute_withoutage_n1[, 1:122], perplexity = 50)
tsne_q1impute_twotype_withoutage_n1 <- tsne(q1impute_withoutage_n1_endo_trac_procedure_data[, 1:122], perplexity = 50)

jpeg(filename = "tsne_q1_imputed_all_procedure.jpg", width = 2400, height = 1600, units = "px", pointsize = 8, quality = 100, bg = "white", res = 300)
ggplot(as.data.frame(tsne_q1impute_withoutage_n1), aes(x = V1, y = V2, color = q1impute_withoutage_n1$procedureType)) + 
  geom_point(cex = 2) + scale_colour_manual(name = "procedureType", values = rainbow(length(unique(q1impute_withoutage_n1$procedureType)))) +
  ggtitle("tSNE of all procedure type with imputed data from 8 groups of features") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

jpeg(filename = "tsne_q1_imputed_two_procedure.jpg", width = 2400, height = 1600, units = "px", pointsize = 8, quality = 100, bg = "white", res = 300)
ggplot(as.data.frame(tsne_q1impute_twotype_withoutage_n1), aes(x = V1, y = V2, color = q1impute_withoutage_n1_endo_trac_procedure_data$procedureType)) + 
  geom_point(cex = 2) + scale_colour_manual(name = "procedureType", values = rainbow(length(unique(q1impute_withoutage_n1_endo_trac_procedure_data$procedureType)))) +
  ggtitle("tSNE of two procedure type with imputed data from 8 groups of features") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

q1impute_withoutage_n1_recur <- cbind(q1impute_withoutage_n1, case_with_procedure$recurrence)
colnames(q1impute_withoutage_n1_recur)[124] <- "recurrence"
tsne_q1impute_withoutage_n1_recur <- tsne(q1impute_withoutage_n1_recur[, 1:122], perplexity = 50)

jpeg(filename = "tsne_q1_imputed_recurrence.jpg", width = 2400, height = 1600, units = "px", pointsize = 8, quality = 100, bg = "white", res = 300)
ggplot(as.data.frame(tsne_q1impute_withoutage_n1_recur), aes(x = V1, y = V2, color = q1impute_withoutage_n1_recur$recurrence)) + 
  geom_point(cex = 2) + scale_colour_manual(name = "recurrence", values = rainbow(length(unique(q1impute_withoutage_n1_recur$recurrence)))) +
  ggtitle("tSNE of recurrence status with imputed data from 8 groups of features") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

## with age and follow-up days

q1impute_withage <- aregImpute(imputef_withage, data = case_with_procedure, n.impute = 10, nk = 0)

q1impute_withage_n1 <- impute.transcan(q1impute_withage, imputation = 1, list.out = T, data = case_with_procedure)

q1impute_withage_n1 <- as.data.frame(q1impute_withage_n1)

rownames(q1impute_withage_n1) <- rownames(case_with_procedure)

q1impute_withage_n1 <- cbind(q1impute_withage_n1, case_with_procedure$procedureType)
colnames(q1impute_withage_n1)[125] <- "procedureType"

q1impute_withage_n1_endo_trac_procedure_data <- subset(q1impute_withage_n1, procedureType %in% c("Endoscopic procedure", "Tracheal resection"))
q1impute_withage_n1_endo_trac_procedure_data$procedureType <- droplevels(q1impute_withage_n1_endo_trac_procedure_data$procedureType)

tsne_q1impute_withage_n1 <- tsne(q1impute_withage_n1_endo_trac_procedure_data[,1:122], perplexity = 50)

ggplot(as.data.frame(tsne_q1impute_withage_n1), aes(x = V1, y = V2, color = q1impute_withage_n1_endo_trac_procedure_data$procedureType)) + 
  geom_point(cex = 2) + scale_colour_manual(name = "procedureType", values = rainbow(length(unique(q1impute_withage_n1_endo_trac_procedure_data$procedureType)))) +
  ggtitle("tSNE of procedureType with imputed data from ages at procedure, fellow-up days and 8 groups of features") +
  theme(plot.title = element_text(hjust = 0.5))

#### Random Forest Feature Selection ####

#rfimputef <- as.formula(paste("procedureType ~", paste(imputeNames, collapse = "+"), "+ age_at_procedure + procedure_to_lastfollowup_days" ))
#d1rfimpute <- rfImpute(rfimputef, data = q1_endo_trac_procedure_data)

rfimputef <- as.formula(paste("procedureType ~", paste(imputeNames, collapse = "+")))

q1rf_impute_withoutage <- aregImpute(rfimputef, data = q1_endo_trac_procedure_data, n.impute = 10, nk = 0)

q1rf_impute_withoutage_n1 <- impute.transcan(q1rf_impute_withoutage, imputation = 1, list.out = T, data = q1_endo_trac_procedure_data)

q1rf_impute_withoutage_n1 <- as.data.frame(q1rf_impute_withoutage_n1)

rownames(q1rf_impute_withoutage_n1) <- rownames(q1_endo_trac_procedure_data)

d1rf <- randomForest(rfimputef, data = q1rf_impute_withoutage_n1, importance = TRUE, keep.forest = TRUE)

d1rf.imp.table <- importance(d1rf)
write.csv(d1rf.imp.table, "random_forest_q1_procedure_ImportanceTable.csv")

jpeg(filename = "random_forest_q1_procedure_variable_importance.jpg", width = 1600, height = 1600, units = "px", pointsize = 8, quality = 100, bg = "white", res = 300)
varImpPlot(d1rf, sort = T, main="Variable Importance", n.var = 30, type = 1)
dev.off()

d1rf.cv <- rfcv(q1rf_impute_withoutage_n1[, 2:ncol(q1rf_impute_withoutage_n1)], q1rf_impute_withoutage_n1$procedureType, cv.fold = 10, recursive = T)

jpeg(filename = "random_forest_q1_procedure_error_rate.jpg", width = 1600, height = 1600, units = "px", pointsize = 8, quality = 100, bg = "white", res = 300)
with(d1rf.cv, plot(n.var, error.cv, log = "x", type = "o", lwd = 2, 
                   xlab = "Number of Variables", ylab = "Error Rate"))
title(main="Estimated Error Rate")
dev.off()

d1cv.result <- replicate(5, rfcv(q1rf_impute_withoutage_n1[, 2:ncol(q1rf_impute_withoutage_n1)], q1rf_impute_withoutage_n1$procedureType), simplify=FALSE)
error.cv <- sapply(d1cv.result, "[[", "error.cv")

jpeg(filename = "random_forest_q1_procedure_cross_validation_error_rate.jpg", width = 1600, height = 1600, units = "px", pointsize = 8, quality = 100, bg = "white", res = 300)
matplot(d1cv.result[[1]]$n.var, cbind(rowMeans(error.cv), error.cv), type="l",
        lwd=c(2, rep(1, ncol(error.cv))), col=1, lty=1, log="x",
        xlab="Number of variables", ylab="CV Error")
title(main="Cross-Validation Error Rate")
dev.off()

q1rf_tsne_data <- q1rf_impute_withoutage_n1[, colnames(q1rf_impute_withoutage_n1) %in% rownames(d1rf.imp.table[order(d1rf.imp.table[, 3], decreasing = T), ])[1:30]]
tsne_q1rf <- tsne(q1rf_tsne_data, perplexity = 50)

jpeg(filename = "random_forest_q1_procedure_tsne.jpg", width = 2400, height = 1600, units = "px", pointsize = 8, quality = 100, bg = "white", res = 300)
ggplot(as.data.frame(tsne_q1rf), aes(x = V1, y = V2, color = q1rf_impute_withoutage_n1$procedureType)) + 
  geom_point(cex = 2) + scale_colour_manual(name = "procedureType", values = rainbow(length(unique(q1rf_impute_withoutage_n1$procedureType)))) +
  ggtitle("tSNE of procedureType with top 30 features selected by random forest") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


######### decision tree

q1_dt <- rpart(rfimputef, data = q1_endo_trac_procedure_data)

jpeg(filename = "decision_tree_q1_procedure_tsne.jpg", width = 2400, height = 1600, units = "px", pointsize = 8, quality = 100, bg = "white", res = 300)
plot(q1_dt)
text(q1_dt)
dev.off()


#q1data.complete <- q1data[complete.cases(q1data), ]
#q1c_outcome <- merge(q1data.complete, outcome.data, by = "row.names")

#rownames(q1c_outcome) <- q1c_outcome[, 1]
#q1c_outcome <- q1c_outcome[, -1]

#colors <- rainbow(length(unique(q1data.complete$procedureType)))
#ecb <- function(x, y) {
#  plot(x, col = colors[q1data.complete$procedureType])
#}
#tsne_q1 <- tsne(q1data.complete[,1:122], epoch_callback = ecb, perplexity=50)


############################## Recurrence ################

case_with_recurrence <- all.data[which(!is.na(all.data$recurrence)), ]

case_with_recurrence <- case_with_recurrence[, -(14:16)] ## remove BIPQ.illnessCause1, BIPQ.illnessCause2, BIPQ.illnessCause3
case_with_recurrence <- case_with_recurrence[, -224]  ### remove MOS.closeFriends

q1data_recur <- case_with_recurrence[, grep(pattern = "BIPQ|CDQ|EAT10|FoPQ|MOS|SDMQ9|SF12v2|VHI10|recurrence.*", x = colnames(case_with_recurrence))]

q1recur.complete <- q1data_recur[complete.cases(q1data_recur), ]

tsne_q1_recur <- tsne(q1recur.complete[, 1:122], perplexity = 50)

ggplot(as.data.frame(tsne_q1_recur), aes(x = V1, y = V2, color = q1recur.complete$recurrence)) + 
  geom_point(cex = 2) + scale_colour_manual(name = "Recurrence", values = rainbow(length(unique(q1recur.complete$recurrence)))) +
  ggtitle("tSNE of recurrence with only original data") +
  theme(plot.title = element_text(hjust = 0.5))








