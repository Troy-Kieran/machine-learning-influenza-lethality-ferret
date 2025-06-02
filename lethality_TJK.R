########################################################################################
########################################################################################
### Lethality spinoff - lethality model improvements
### Influenza-Ferret Predictive Analytics
### 31 October 2023 - 27 March 2024
### Troy J. Kieran
########################################################################################
########################################################################################
### Load packages
library(tidyverse)
library(tidylog)
library(caret)
library(gbm)

########################################################################################
### Import Data

## replace file.csv with name of file

## download data
# Pathogenesis Laboratory Team, Influenza Division, CDC. 
# An aggregated dataset of serially collected influenza A virus morbidity and titer measurements from virus-infected ferrets.  
# https://data.cdc.gov/National-Center-for-Immunization-and-Respiratory-D/An-aggregated-dataset-of-serially-collected-influe/cr56-k9wj/about_data
# and
# An aggregated dataset of day 3 post-inoculation viral titer measurments from influenza A virus-infected ferret tissues.
# https://data.cdc.gov/National-Center-for-Immunization-and-Respiratory-D/An-aggregated-dataset-of-day-3-post-inoculation-vi/d9u6-mdu6/about_data
#
# fullData <- read.csv("file1.csv", header = TRUE, check.names = FALSE)

########################################################################################

fullData %>% 
  select(-c(contains('DC'), contains('RD'))) %>%
  filter(expt == 'path') %>%
  skimr::skim()

## Convert the NA tissue values to respective LOD value
fullData <- fullData %>%
  mutate(Lg_avg = ifelse(is.na(Lg_avg) & Lg_rep == "no" & units == "EID", 1.5,
                     ifelse(is.na(Lg_avg) & Lg_rep == "no" & units == "London", 1,
                            Lg_avg)),
         BnOB_avg = ifelse(is.na(BnOB_avg) & BnOB_rep == "no" & units == "EID", 1.5,
                     ifelse(is.na(BnOB_avg) & BnOB_rep == "no" & units == "London", 1,
                            BnOB_avg)),
         Bn_avg = ifelse(is.na(Bn_avg) & Bn_rep == "no" & units == "EID", 1.5,
                     ifelse(is.na(Bn_avg) & Bn_rep == "no" & units == "London", 1,
                            Bn_avg)),
         Int_avg = ifelse(is.na(Int_avg) & Int_rep == "no" & units == "EID", 1.5,
                     ifelse(is.na(Int_avg) & Int_rep == "no" & units == "London", 1,
                            Int_avg)))

########################################################################################

## Predictive Power Scores
## lethality
fullData %>%
  dplyr::select(lethal, HPAI_MBAA, RBS, PBS, HA, wt_loss, temp_5, AUC_6_f, 
                NT_avg, Lg_avg, BnOB_avg, Bn_avg) %>%
  ppsr::score_predictors(df = ., y = 'lethal')

fullData %>%
  dplyr::select(lethal, HPAI_MBAA, RBS, PBS, HA, 
                NT_avg, Lg_avg, BnOB_avg, Bn_avg) %>%
  ppsr::score_predictors(df = ., y = 'lethal')

########################################################################################

fitControl <- trainControl(method = "repeatedcv",   
                           number = 10,  # number of folds
                           repeats = 2,  # repeated two times = 20 folds
                           savePredictions = 'final',
                           classProbs = TRUE,
                           summaryFunction = multiClassSummary)

###

leth_test_base <- fullData %>% 
  drop_na(lethal) %>%
  #filter(expt == 'path') %>%
  dplyr::select(c(lethal, HPAI_MBAA, RBS, PBS, HA,
                  wt_loss, temp_5, AUC_6_f))

leth_test1a <- fullData %>% 
  drop_na(lethal) %>%
  #filter(expt == 'path') %>%
  dplyr::select(c(lethal, HPAI_MBAA, RBS, PBS, HA,
                  wt_loss, temp_5, AUC_6_f, NT_avg))

leth_test1b <- fullData %>% 
  drop_na(lethal) %>%
  #filter(expt == 'path') %>%
  dplyr::select(c(lethal, HPAI_MBAA, RBS, PBS, HA,
                  wt_loss, temp_5, AUC_6_f, NT_rep))

leth_test2a <- fullData %>% 
  drop_na(lethal) %>%
  #filter(expt == 'path') %>%
  dplyr::select(c(lethal, HPAI_MBAA, RBS, PBS, HA,
                  wt_loss, temp_5, AUC_6_f, Lg_avg))

leth_test2b <- fullData %>% 
  drop_na(lethal) %>%
  #filter(expt == 'path') %>%
  dplyr::select(c(lethal, HPAI_MBAA, RBS, PBS, HA,
                  wt_loss, temp_5, AUC_6_f, Lg_rep))

leth_test3a <- fullData %>% 
  drop_na(lethal) %>%
  #filter(expt == 'path') %>%
  dplyr::select(c(lethal, HPAI_MBAA, RBS, PBS, HA,
                  wt_loss, temp_5, AUC_6_f, BnOB_avg))

leth_test3b <- fullData %>% 
  drop_na(lethal) %>%
  #filter(expt == 'path') %>%
  dplyr::select(c(lethal, HPAI_MBAA, RBS, PBS, HA,
                  wt_loss, temp_5, AUC_6_f, BnOB_rep))

leth_test4a <- fullData %>% 
  drop_na(lethal) %>%
  #filter(expt == 'path') %>%
  dplyr::select(c(lethal, HPAI_MBAA, RBS, PBS, HA,
                  wt_loss, temp_5, AUC_6_f, Bn_avg))

leth_test4b <- fullData %>% 
  drop_na(lethal) %>%
  #filter(expt == 'path') %>%
  dplyr::select(c(lethal, HPAI_MBAA, RBS, PBS, HA,
                  wt_loss, temp_5, AUC_6_f, Bn_rep))

leth_test5 <- fullData %>% 
  drop_na(lethal) %>%
  #filter(expt == 'path') %>%
  dplyr::select(c(lethal, HPAI_MBAA, RBS, PBS, HA,
                  wt_loss, temp_5, AUC_6_f, 
                  NT_avg, Lg_avg, BnOB_avg, Bn_avg))

leth_test6 <- fullData %>% 
  drop_na(lethal) %>%
  #filter(expt == 'path') %>%
  dplyr::select(c(lethal, HPAI_MBAA, RBS, PBS, HA,
                  NT_avg, Lg_avg, BnOB_avg, Bn_avg))

leth_test7a <- fullData %>% 
  drop_na(lethal) %>%
  #filter(expt == 'path') %>%
  dplyr::select(c(lethal, HPAI_MBAA, RBS, PBS, HA,
                  NT_avg))

leth_test7b <- fullData %>% 
  drop_na(lethal) %>%
  #filter(expt == 'path') %>%
  dplyr::select(c(lethal, HPAI_MBAA, RBS, PBS, HA,
                  Lg_avg))

leth_test7c <- fullData %>% 
  drop_na(lethal) %>%
  #filter(expt == 'path') %>%
  dplyr::select(c(lethal, HPAI_MBAA, RBS, PBS, HA,
                  BnOB_avg))

leth_test7d <- fullData %>% 
  drop_na(lethal) %>%
  #filter(expt == 'path') %>%
  dplyr::select(c(lethal, HPAI_MBAA, RBS, PBS, HA,
                  Bn_avg))

leth_test_base_dummy <- fastDummies::dummy_cols(
  leth_test_base, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'HA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

leth_test1a_dummy <- fastDummies::dummy_cols(
  leth_test1a, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'HA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

leth_test2a_dummy <- fastDummies::dummy_cols(
  leth_test2a, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'HA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

leth_test3a_dummy <- fastDummies::dummy_cols(
  leth_test3a, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'HA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

leth_test4a_dummy <- fastDummies::dummy_cols(
  leth_test4a, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'HA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

leth_test1b_dummy <- fastDummies::dummy_cols(
  leth_test1b, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'HA', 'NT_rep'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

leth_test2b_dummy <- fastDummies::dummy_cols(
  leth_test2b, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'HA', 'Lg_rep'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

leth_test3b_dummy <- fastDummies::dummy_cols(
  leth_test3b, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'HA', 'BnOB_rep'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

leth_test4b_dummy <- fastDummies::dummy_cols(
  leth_test4b, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'HA', 'Bn_rep'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

leth_test5_dummy <- fastDummies::dummy_cols(
  leth_test5, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'HA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

leth_test6_dummy <- fastDummies::dummy_cols(
  leth_test6, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'HA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

leth_test7a_dummy <- fastDummies::dummy_cols(
  leth_test7a, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'HA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

leth_test7b_dummy <- fastDummies::dummy_cols(
  leth_test7b, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'HA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

leth_test7c_dummy <- fastDummies::dummy_cols(
  leth_test7c, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'HA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

leth_test7d_dummy <- fastDummies::dummy_cols(
  leth_test7d, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'HA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

leth_test_base_dummy$lethal <- as.factor(leth_test_base_dummy$lethal)
leth_test1a_dummy$lethal <- as.factor(leth_test1a_dummy$lethal)
leth_test1b_dummy$lethal <- as.factor(leth_test1b_dummy$lethal)
leth_test2a_dummy$lethal <- as.factor(leth_test2a_dummy$lethal)
leth_test2b_dummy$lethal <- as.factor(leth_test2b_dummy$lethal)
leth_test3a_dummy$lethal <- as.factor(leth_test3a_dummy$lethal)
leth_test3b_dummy$lethal <- as.factor(leth_test3b_dummy$lethal)
leth_test4a_dummy$lethal <- as.factor(leth_test4a_dummy$lethal)
leth_test4b_dummy$lethal <- as.factor(leth_test4b_dummy$lethal)
leth_test5_dummy$lethal <- as.factor(leth_test5_dummy$lethal)
leth_test6_dummy$lethal <- as.factor(leth_test6_dummy$lethal)
leth_test7a_dummy$lethal <- as.factor(leth_test7a_dummy$lethal)
leth_test7b_dummy$lethal <- as.factor(leth_test7b_dummy$lethal)
leth_test7c_dummy$lethal <- as.factor(leth_test7c_dummy$lethal)
leth_test7d_dummy$lethal <- as.factor(leth_test7d_dummy$lethal)

set.seed(9595)
leth_base_dumSplit <- rsample::initial_split(leth_test_base_dummy, prop = 0.70)
leth_1a_dumSplit <- rsample::initial_split(leth_test1a_dummy, prop = 0.70)
leth_1b_dumSplit <- rsample::initial_split(leth_test1b_dummy, prop = 0.70)
leth_2a_dumSplit <- rsample::initial_split(leth_test2a_dummy, prop = 0.70)
leth_2b_dumSplit <- rsample::initial_split(leth_test2b_dummy, prop = 0.70)
leth_3a_dumSplit <- rsample::initial_split(leth_test3a_dummy, prop = 0.70)
leth_3b_dumSplit <- rsample::initial_split(leth_test3b_dummy, prop = 0.70)
leth_4a_dumSplit <- rsample::initial_split(leth_test4a_dummy, prop = 0.70)
leth_4b_dumSplit <- rsample::initial_split(leth_test4b_dummy, prop = 0.70)
leth_5_dumSplit <- rsample::initial_split(leth_test5_dummy, prop = 0.70)
leth_6_dumSplit <- rsample::initial_split(leth_test6_dummy, prop = 0.70)
leth_7a_dumSplit <- rsample::initial_split(leth_test7a_dummy, prop = 0.70)
leth_7b_dumSplit <- rsample::initial_split(leth_test7b_dummy, prop = 0.70)
leth_7c_dumSplit <- rsample::initial_split(leth_test7c_dummy, prop = 0.70)
leth_7d_dumSplit <- rsample::initial_split(leth_test7d_dummy, prop = 0.70)

trainData_base <- rsample::training(leth_base_dumSplit)
testData_base <- rsample::testing(leth_base_dumSplit)
trainData_1a <- rsample::training(leth_1a_dumSplit)
testData_1a <- rsample::testing(leth_1a_dumSplit)
trainData_1b <- rsample::training(leth_1b_dumSplit)
testData_1b <- rsample::testing(leth_1b_dumSplit)
trainData_2a <- rsample::training(leth_2a_dumSplit)
testData_2a <- rsample::testing(leth_2a_dumSplit)
trainData_2b <- rsample::training(leth_2b_dumSplit)
testData_2b <- rsample::testing(leth_2b_dumSplit)
trainData_3a <- rsample::training(leth_3a_dumSplit)
testData_3a <- rsample::testing(leth_3a_dumSplit)
trainData_3b <- rsample::training(leth_3b_dumSplit)
testData_3b <- rsample::testing(leth_3b_dumSplit)
trainData_4a <- rsample::training(leth_4a_dumSplit)
testData_4a <- rsample::testing(leth_4a_dumSplit)
trainData_4b <- rsample::training(leth_4b_dumSplit)
testData_4b <- rsample::testing(leth_4b_dumSplit)
trainData_5 <- rsample::training(leth_5_dumSplit)
testData_5 <- rsample::testing(leth_5_dumSplit)
trainData_6 <- rsample::training(leth_6_dumSplit)
testData_6 <- rsample::testing(leth_6_dumSplit)
trainData_7a <- rsample::training(leth_7a_dumSplit)
testData_7a <- rsample::testing(leth_7a_dumSplit)
trainData_7b <- rsample::training(leth_7b_dumSplit)
testData_7b <- rsample::testing(leth_7b_dumSplit)
trainData_7c <- rsample::training(leth_7c_dumSplit)
testData_7c <- rsample::testing(leth_7c_dumSplit)
trainData_7d <- rsample::training(leth_7d_dumSplit)
testData_7d <- rsample::testing(leth_7d_dumSplit)

set.seed(2626)
mbase <- train(lethal ~ ., data = trainData_base,
             method = "gbm",
             na.action = na.exclude,
             preProcess = c("nzv", "scale", "center"),
             trControl = fitControl)

set.seed(2626)
m1a <- train(lethal ~ ., data = trainData_1a,
               method = "gbm",
               na.action = na.exclude,
               preProcess = c("nzv", "scale", "center"),
               trControl = fitControl)

set.seed(2626)
m1b <- train(lethal ~ ., data = trainData_1b,
             method = "gbm",
             na.action = na.exclude,
             preProcess = c("nzv", "scale", "center"),
             trControl = fitControl)

set.seed(2626)
m2a <- train(lethal ~ ., data = trainData_2a,
             method = "gbm",
             na.action = na.exclude,
             preProcess = c("nzv", "scale", "center"),
             trControl = fitControl)

set.seed(2626)
m2b <- train(lethal ~ ., data = trainData_2b,
             method = "gbm",
             na.action = na.exclude,
             preProcess = c("nzv", "scale", "center"),
             trControl = fitControl)

set.seed(2626)
m3a <- train(lethal ~ ., data = trainData_3a,
             method = "gbm",
             na.action = na.exclude,
             preProcess = c("nzv", "scale", "center"),
             trControl = fitControl)

set.seed(2626)
m3b <- train(lethal ~ ., data = trainData_3b,
             method = "gbm",
             na.action = na.exclude,
             preProcess = c("nzv", "scale", "center"),
             trControl = fitControl)

set.seed(2626)
m4a <- train(lethal ~ ., data = trainData_4a,
             method = "gbm",
             na.action = na.exclude,
             preProcess = c("nzv", "scale", "center"),
             trControl = fitControl)

set.seed(2626)
m4b <- train(lethal ~ ., data = trainData_4b,
             method = "gbm",
             na.action = na.exclude,
             preProcess = c("nzv", "scale", "center"),
             trControl = fitControl)

set.seed(2626)
m5 <- train(lethal ~ ., data = trainData_5,
             method = "gbm",
             na.action = na.exclude,
             preProcess = c("nzv", "scale", "center"),
             trControl = fitControl)

set.seed(2626)
m6 <- train(lethal ~ ., data = trainData_6,
            method = "gbm",
            na.action = na.exclude,
            preProcess = c("nzv", "scale", "center"),
            trControl = fitControl)

set.seed(2626)
m7a <- train(lethal ~ ., data = trainData_7a,
            method = "gbm",
            na.action = na.exclude,
            preProcess = c("nzv", "scale", "center"),
            trControl = fitControl)

set.seed(2626)
m7b <- train(lethal ~ ., data = trainData_7b,
             method = "gbm",
             na.action = na.exclude,
             preProcess = c("nzv", "scale", "center"),
             trControl = fitControl)

set.seed(2626)
m7c <- train(lethal ~ ., data = trainData_7c,
             method = "gbm",
             na.action = na.exclude,
             preProcess = c("nzv", "scale", "center"),
             trControl = fitControl)

set.seed(2626)
m7d <- train(lethal ~ ., data = trainData_7d,
             method = "gbm",
             na.action = na.exclude,
             preProcess = c("nzv", "scale", "center"),
             trControl = fitControl)

## resample metrics
resamps <- resamples(list(Base = mbase,
                          NT_avg = m1a,
                          NT_rep = m1b,
                          Lg_avg = m2a,
                          Lg_rep = m2b,
                          BnOB_avg = m3a,
                          BnOB_rep = m3b,
                          Bn_avg = m4a,
                          Bn_rep = m4b,
                          All_avg = m5,
                          All_avg2 = m6,
                          NT_avg2 = m7a,
                          Lg_avg2 = m7b,
                          BnOB_avg2 = m7c,
                          Bn_avg2 = m7d))
summary(resamps)
bwplot(resamps)
dotplot(resamps)
modelCor(resamps)
splom(resamps)

## check for variable importance
#ggplot(varImp(mbase))
varImp(mbase)
varImp(m1a)
varImp(m1b)
varImp(m2a)
varImp(m2b)
varImp(m3a)
varImp(m3b)
varImp(m4a)
varImp(m4b)
varImp(m5)
varImp(m6)
varImp(m7a)
varImp(m7b)
varImp(m7c)
varImp(m7d)

## rotate through
#testData <- testData_base
#model <- mbase
#testData <- testData_1a
#model <- m1a
#testData <- testData_1b
#model <- m1b
#testData <- testData_2a
#model <- m2a
#testData <- testData_2b
#model <- m2b
#testData <- testData_3a
#model <- m3a
#testData <- testData_3b
#model <- m3b
#testData <- testData_4a
#model <- m4a
#testData <- testData_4b
#model <- m4b
#testData <- testData_5
#model <- m5
#testData <- testData_6
#model <- m6
#testData <- testData_7a
#model <- m7a
#testData <- testData_7b
#model <- m7b
#testData <- testData_7c
#model <- m7c
testData <- testData_7d
model <- m7d

probTest <- predict(model, testData, type = 'prob')
probTest <- factor(ifelse(probTest$no >= 0.94, 'no', 'yes')) %>%
  as.data.frame()
probTruth <- testData %>%
  na.omit() %>%
  dplyr::select(lethal) %>%
  cbind(., probTest)
confusionMatrix(probTruth$., probTruth$lethal, mode = "everything")

## check a range of probability values
thresholds <- seq(0.5, 0.95, by = 0.05)

## empty list to populate
prob_results <- list()
## function
for (threshold in thresholds) {
  probTest <- predict(model, testData, type = 'prob')
  probTest <- factor(ifelse(probTest$no >= threshold, 'no', 'yes')) %>%
    as.data.frame()
  
  probTruth <- testData %>%
    na.omit() %>%
    dplyr::select(lethal) %>%
    cbind(., probTest)
  
  confusion_matrix <- confusionMatrix(probTruth$., probTruth$lethal, mode = "everything")
  
  prob_results[[as.character(threshold)]] <- list(
    threshold = threshold,
    confusion_matrix = confusion_matrix)}

## Access the results for each threshold, e.g., results[['0.5']]
prob_results[['0.5']]
prob_results[['0.55']]
prob_results[['0.6']]
prob_results[['0.65']]
prob_results[['0.7']]
prob_results[['0.75']]
prob_results[['0.8']]
prob_results[['0.85']]
prob_results[['0.9']]
prob_results[['0.95']]

########################################################################################

## Matthew's Correlation Coefficient

## manual calculation/function
## T = truth, F = false, P = positive, N = negative
## use data from confusion matrix
mcc_func <- function(TP, FP, FN, TN){
  ((TP*TN) - (FP*FN))/
    sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))}

## path expt only
## base
mcc_func(47, 4, 2, 17)
## NT
mcc_func(32, 2, 0, 17)
mcc_func(33, 4, 2, 14)
## Lg
mcc_func(30, 5, 2, 14)
mcc_func(46, 3, 0, 9)
## BnOB
mcc_func(37, 2, 0, 8)
mcc_func(33, 4, 0, 13)
## Bn
mcc_func(32, 2, 2, 14)
mcc_func(41, 4, 3, 12)
## All_avg
mcc_func(23, 2, 2, 17)
## All_avg2
mcc_func(29, 4, 1, 18)
##
mcc_func(36, 4, 3, 11)
mcc_func(32, 6, 2, 16)
mcc_func(29, 3, 0, 16)
mcc_func(39, 2, 5, 13)

## All
mcc_func(171, 13, 0, 22)
mcc_func(120, 8, 1, 23)
mcc_func(126, 12, 0, 30)
mcc_func(133, 0, 0, 26)
mcc_func(136, 7, 1, 24)
mcc_func(127, 7, 0, 26)
mcc_func(118, 9, 8, 34)
mcc_func(129, 19, 5, 26)
mcc_func(132, 12, 4, 23)
mcc_func(133, 4, 8, 21)
mcc_func(125, 22, 2, 24)


########################################################################################
########################################################################################

### Section below combines data from lethality & morbidity for plot making

########################################################################################
########################################################################################

### Sample Sizes, Model Metrics, and Features Heatmaps/Figures

features <- read.csv("ResultsSummary_Rinputs/Features4Heatmap.csv", header = TRUE)
modelMetrics <- read.csv("ResultsSummary_Rinputs/ModelMetrics4Heatmap.csv", header = TRUE)
sample_size_data <- read.csv("ResultsSummary_Rinputs/SampleSizeData.csv", header = TRUE)
SelectedImportance <- read.csv("ResultsSummary_Rinputs/SelectedImportance.csv", header = TRUE)

###

## Select relevant columns for creating the presence-absence matrix
presence_matrix <- features[, 5:11] != ""
## Convert to binary matrix format (True for presence, False for absence)
binary_matrix <- as.matrix(presence_matrix)
## drop text Feature columns
features <- features[, 1:4]
## combine binary to features
features_binary <- cbind(features, binary_matrix)   
## convert TRUE/FALSE to 1/0
features_binary <- features_binary %>%
  mutate(across(where(is.logical), ~ if_else(.x, 1, 0)))

gathered_data <- features_binary %>%
  gather(key = "Feature", value = "Presence", 5:11)

## should separate out the Types and then patchwork together rather than facet.
gathered_data$Test <- as.factor(gathered_data$Test)

gathered_data$Type <- factor(gathered_data$Type, 
                             levels=c('Base', 'Bn', 'BnOB', 'Lg', 
                                      'NT', 'All', 'Tissue'))

gathered_data$Feature <- factor(gathered_data$Feature, 
                             levels=c('Bn', 'BnOB', 'Lg', 'NT', 
                                      'AUC_6', 'temp_5', 'wt_loss'))
## subset data
Lethality_data <- gathered_data %>% filter(Classification == 'Lethality')
Morbidity_data <- gathered_data %>% filter(Classification == 'Morbidity')

## make subplots
Lethality_plot <- 
  ggplot(data = Lethality_data) +
  geom_tile(aes(x = Feature, y = Type, fill = Presence)) +
  scale_fill_viridis_c(direction = -1, end = 0.7) +
  theme_minimal() +
  facet_grid(~ Classification, scales = 'free', space = 'free') +
  ggtitle('Features') +
  theme(legend.position = 'none', axis.title.y = element_blank(),
        strip.text.x = element_blank(), axis.text.y = element_blank(),
        plot.title = element_text(hjust = 0.5), axis.title.x = element_blank()) +
  scale_x_discrete(position = "top")

Morbidity_plot <- 
  ggplot(data = Morbidity_data) +
  geom_tile(aes(x = Feature, y = Type, fill = Presence)) +
  scale_fill_viridis_c(direction = -1, end = 0.7) +
  theme_minimal() +
  facet_grid(~ Classification, scales = 'free', space = 'free') +
  theme(legend.position = 'none', axis.title.y = element_blank(),
        strip.text.x = element_blank(), axis.text.y = element_blank(),
        axis.title.x = element_blank()) +
  scale_x_discrete(position = "top")

###

## final models metric heatmap
modelMetrics <- modelMetrics %>% 
  gather(L.base:M.Tissue, key = Model, value = Value, factor_key = TRUE) %>%
  dplyr::filter(Metric != 'test') %>%
  drop_na()

modelMetrics$Metric <- factor(modelMetrics$Metric, 
                              levels=c('BA', 'F1', 'MCC'))

modelMetrics_L <- modelMetrics %>% filter(str_detect(Model, "^L"))
modelMetrics_M <- modelMetrics %>% filter(str_detect(Model, "^M"))

modelMetrics_L$Model <- factor(modelMetrics_L$Model, 
                               levels=c('L.base', 'L.Bn', 'L.BnOB', 'L.Lg', 
                                        'L.NT', 'L.All', 'L.Tissue'))

modelMetrics_M$Model <- factor(modelMetrics_M$Model, 
                               levels=c('M.base', 'M.Bn', 'M.BnOB', 'M.Lg', 
                                        'M.NT', 'M.All', 'M.Tissue'))

modelMetrics_L$Value <- as.numeric(modelMetrics_L$Value)
modelMetrics_M$Value <- as.numeric(modelMetrics_M$Value)

viridis_scale <- scale_fill_viridis_c(limits = c(0.3, 1), direction = -1, end = 0.9)

modelMetrics_Lplot <- 
  ggplot(modelMetrics_L, aes(Metric, Model, fill = Value)) + 
  geom_tile() +
  geom_text(aes(label = Value), color = "white", size = 3) +
  viridis_scale + 
  theme_minimal() +
  ggtitle('Metrics') +
  theme(legend.position = 'none', plot.title = element_text(hjust = 0.5), 
        axis.title.x = element_blank(), axis.title.y = element_blank()) +
  scale_x_discrete(position = "top")

modelMetrics_Mplot <- 
  ggplot(modelMetrics_M, aes(Metric, Model, fill = Value)) + 
  geom_tile() +
  geom_text(aes(label = Value), color = "white", size = 3) +
  viridis_scale + 
  #scale_fill_viridis_d(direction = -1, begin = 0.1, end = 0.85) +
  theme_minimal() +
  theme(legend.position = 'none', axis.title.x = element_blank(), 
        axis.title.y = element_blank()) +
  scale_x_discrete(position = "top")

###

sample_size_data_L <- sample_size_data %>% filter(str_detect(Model, "^L"))
sample_size_data_M <- sample_size_data %>% filter(str_detect(Model, "^M"))

sample_size_data_L$Model <- factor(sample_size_data_L$Model, 
                               levels=c('L-base', 'L-Bn', 'L-BnOB', 'L-Lg', 
                                        'L-NT', 'L-All', 'L-Tissue'))

sample_size_data_M$Model <- factor(sample_size_data_M$Model, 
                                   levels=c('M-base', 'M-Bn', 'M-BnOB', 'M-Lg', 
                                            'M-NT', 'M-All', 'M-Tissue'))

sample_size_Lplot <- ggplot(sample_size_data_L, aes(x = x, y = Model)) +
  geom_tile(fill = 'dodgerblue4') +
  geom_text(aes(label = Virus), color = 'white', size = 3) +
  theme_minimal() +
  ggtitle('Virus') +
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        plot.title = element_text(size = 9, hjust = 0.5))

sample_size_Lplot2 <- ggplot(sample_size_data_L, aes(x = x, y = Model)) +
  geom_tile(fill = 'dodgerblue3') +
  geom_text(aes(label = Obs_yes), color = 'white', size = 2.8) +
  theme_minimal() +
  ggtitle('Obs(yes)') +
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        plot.title = element_text(size = 9, hjust = 0.5))

sample_size_Mplot <- ggplot(sample_size_data_M, aes(x = x, y = Model)) +
  geom_tile(fill = 'dodgerblue4') +
  geom_text(aes(label = Virus), color = 'white', size = 3) +
  theme_minimal() +
  ggtitle('Virus') +
  theme(axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size = 9, hjust = 0.5))

sample_size_Mplot2 <- ggplot(sample_size_data_M, aes(x = x, y = Model)) +
  geom_tile(fill = 'dodgerblue3') +
  geom_text(aes(label = Obs_yes), color = 'white', size = 2.8) +
  theme_minimal() +
  ggtitle('Obs(yes)') +
  theme(axis.text = element_blank(), axis.title = element_blank(),
        plot.title = element_text(size = 9, hjust = 0.5))

layout <- ((sample_size_Lplot | sample_size_Lplot2) / 
             (sample_size_Mplot | sample_size_Mplot2) |
             (modelMetrics_Lplot / modelMetrics_Mplot) | 
             (Lethality_plot / Morbidity_plot))

layout + patchwork::plot_layout(widths = c(0.55, 0.8, 1.4)) 

########################################################################################

## Feature Importance Figure

## barplots with selected importance variables

SelectedImportance_L <- read.csv("ResultsSummary_Rinputs/SelectedImportance_L.csv", header = TRUE)
SelectedImportance_M <- read.csv("ResultsSummary_Rinputs/SelectedImportance_M.csv", header = TRUE)

## Lethality
SelectedImportance_L$Model <- factor(SelectedImportance_L$Model, 
                                     levels=c('L-base', 'L-All', 'L-Tissue'))

SelectedImportance_L_num <- SelectedImportance_L %>%
  filter(Catergory == 'Numeric') %>%
  mutate(Feature = fct_reorder(factor(Feature), Importance, .desc = FALSE))

SelectedImportance_L_fac <- SelectedImportance_L %>%
  filter(Catergory != 'Numeric')

SelectedImportance_L_num_plot <- 
  ggplot(SelectedImportance_L_num, 
         aes(x = Feature, y = Importance, fill = Feature)) +
  geom_col() +
  facet_wrap(~ Model, scales = "free_y") +
  coord_flip() +
  ylab("Relative Ranked Importance") +
  scale_fill_viridis_d(option = 'B', begin = 0.1, end = 0.8) +
  theme_minimal() +
  theme(legend.position = 'none')

SelectedImportance_L_fac_plot <- 
  ggplot(SelectedImportance_L_fac, 
         aes(x = Feature, y = Importance, fill = Feature)) +
  geom_col() +
  facet_wrap(~ Model) +
  coord_flip() +
  ylab("Relative Ranked Importance") +
  scale_fill_viridis_d(option = 'G', end = 0.8) +
  theme_minimal() +
  theme(legend.position = 'none') +
  scale_y_continuous(expand = expansion(add = c(0, 25)))

# SelectedImportance_L_num_plot / 
#   SelectedImportance_L_fac_plot + 
#   patchwork::plot_annotation(tag_levels = 'A')

###

## Morbidity
SelectedImportance_M$Model <- factor(SelectedImportance_M$Model, 
                                     levels=c('M-base', 'M-All', 'M-Tissue'))

SelectedImportance_M_num <- SelectedImportance_M %>%
  filter(Catergory == 'Numeric') %>%
  mutate(Feature = fct_reorder(factor(Feature), Importance, .desc = FALSE))

SelectedImportance_M_fac <- SelectedImportance_M %>%
  filter(Catergory != 'Numeric')

SelectedImportance_M_num_plot <- 
  ggplot(SelectedImportance_M_num, 
         aes(x = Feature, y = Importance, fill = Feature)) +
  geom_col() +
  facet_wrap(~ Model, scales = "free_y") +
  coord_flip() +
  ylab("Relative Ranked Importance") +
  scale_fill_viridis_d(option = 'B', begin = 0.1, end = 0.8) +
  theme_minimal() +
  theme(legend.position = 'none')

SelectedImportance_M_fac_plot <- 
  ggplot(SelectedImportance_M_fac, 
         aes(x = Feature, y = Importance, fill = Feature)) +
  geom_col() +
  facet_wrap(~ Model) +
  coord_flip() +
  ylab("Relative Ranked Importance") +
  scale_fill_viridis_d(option = 'G', end = 0.8) +
  theme_minimal() +
  theme(legend.position = 'none')

# SelectedImportance_M_num_plot / 
#   SelectedImportance_M_fac_plot + 
#   patchwork::plot_annotation(tag_levels = 'A')


## combine lethality & morbidity plots
(SelectedImportance_L_num_plot | SelectedImportance_L_fac_plot) / 
  (SelectedImportance_M_num_plot | SelectedImportance_M_fac_plot) + 
  patchwork::plot_annotation(tag_levels = 'A')

########################################################################################
########################################################################################
### End of Code
########################################################################################
########################################################################################
