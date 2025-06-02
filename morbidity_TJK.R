########################################################################################
########################################################################################
### Lethality spinoff - morbidity model improvements
### Influenza-Ferret Predictive Analytics
### 31 October 2023 - 30 January 2024
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
## morbidity
fullData %>%
  dplyr::select(wt_loss_high, HPAI_MBAA, RBS, PBS, HA, wt_loss, temp_5, AUC_6_f, 
                NT_avg, Lg_avg, BnOB_avg, Bn_avg) %>%
  ppsr::score_predictors(df = ., y = 'wt_loss_high')

fullData %>%
  dplyr::select(wt_loss_high, HPAI_MBAA, RBS, PBS, HA, 
                NT_avg, Lg_avg, BnOB_avg, Bn_avg) %>%
  ppsr::score_predictors(df = ., y = 'wt_loss_high')

########################################################################################

fitControl <- trainControl(method = "repeatedcv",   
                           number = 10,  # number of folds
                           repeats = 2,  # repeated two times = 20 folds
                           savePredictions = 'final',
                           classProbs = TRUE,
                           summaryFunction = multiClassSummary)

###

wt_test_base <- fullData %>% 
  drop_na(wt_loss_high) %>%
  #filter(expt == 'path') %>%
  dplyr::select(c(wt_loss_high, HPAI_MBAA, RBS, PBS, HA,
                  temp_5, AUC_6_f))

wt_test1a <- fullData %>% 
  drop_na(wt_loss_high) %>%
  #filter(expt == 'path') %>%
  dplyr::select(c(wt_loss_high, HPAI_MBAA, RBS, PBS, HA,
                  temp_5, AUC_6_f, NT_avg))

wt_test1b <- fullData %>% 
  drop_na(wt_loss_high) %>%
  #filter(expt == 'path') %>%
  dplyr::select(c(wt_loss_high, HPAI_MBAA, RBS, PBS, HA,
                  temp_5, AUC_6_f, NT_rep))

wt_test2a <- fullData %>% 
  drop_na(wt_loss_high) %>%
  #filter(expt == 'path') %>%
  dplyr::select(c(wt_loss_high, HPAI_MBAA, RBS, PBS, HA,
                  temp_5, AUC_6_f, Lg_avg))

wt_test2b <- fullData %>% 
  drop_na(wt_loss_high) %>%
  #filter(expt == 'path') %>%
  dplyr::select(c(wt_loss_high, HPAI_MBAA, RBS, PBS, HA,
                  temp_5, AUC_6_f, Lg_rep))

wt_test3a <- fullData %>% 
  drop_na(wt_loss_high) %>%
  #filter(expt == 'path') %>%
  dplyr::select(c(wt_loss_high, HPAI_MBAA, RBS, PBS, HA,
                  temp_5, AUC_6_f, BnOB_avg))

wt_test3b <- fullData %>% 
  drop_na(wt_loss_high) %>%
  #filter(expt == 'path') %>%
  dplyr::select(c(wt_loss_high, HPAI_MBAA, RBS, PBS, HA,
                  temp_5, AUC_6_f, BnOB_rep))

wt_test4a <- fullData %>% 
  drop_na(wt_loss_high) %>%
  #filter(expt == 'path') %>%
  dplyr::select(c(wt_loss_high, HPAI_MBAA, RBS, PBS, HA,
                  temp_5, AUC_6_f, Bn_avg))

wt_test4b <- fullData %>% 
  drop_na(wt_loss_high) %>%
  #filter(expt == 'path') %>%
  dplyr::select(c(wt_loss_high, HPAI_MBAA, RBS, PBS, HA,
                  temp_5, AUC_6_f, Bn_rep))

wt_test5 <- fullData %>% 
  drop_na(wt_loss_high) %>%
  #filter(expt == 'path') %>%
  dplyr::select(c(wt_loss_high, HPAI_MBAA, RBS, PBS, HA,
                  temp_5, AUC_6_f, 
                  NT_avg, Lg_avg, BnOB_avg, Bn_avg))

wt_test6 <- fullData %>% 
  drop_na(wt_loss_high) %>%
  #filter(expt == 'path') %>%
  dplyr::select(c(wt_loss_high, HPAI_MBAA, RBS, PBS, HA,
                  NT_avg, Lg_avg, BnOB_avg, Bn_avg))

wt_test7a <- fullData %>% 
  drop_na(wt_loss_high) %>%
  #filter(expt == 'path') %>%
  dplyr::select(c(wt_loss_high, HPAI_MBAA, RBS, PBS, HA,
                  NT_avg))

wt_test7b <- fullData %>% 
  drop_na(wt_loss_high) %>%
  #filter(expt == 'path') %>%
  dplyr::select(c(wt_loss_high, HPAI_MBAA, RBS, PBS, HA,
                  Lg_avg))

wt_test7c <- fullData %>% 
  drop_na(wt_loss_high) %>%
  #filter(expt == 'path') %>%
  dplyr::select(c(wt_loss_high, HPAI_MBAA, RBS, PBS, HA,
                  BnOB_avg))

wt_test7d <- fullData %>% 
  drop_na(wt_loss_high) %>%
  #filter(expt == 'path') %>%
  dplyr::select(c(wt_loss_high, HPAI_MBAA, RBS, PBS, HA,
                  Bn_avg))

## alternative one hot coding using recipes
# wt_test_base %>% 
#   recipes::recipe(wt_loss_high ~ .) %>% 
#   recipes::step_dummy(HPAI_MBAA, RBS, PBS, HA, one_hot = TRUE) %>% 
#   recipes::prep() %>% 
#   recipes::bake(wt_test_base)

wt_test_base_dummy <- fastDummies::dummy_cols(
  wt_test_base, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'HA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

wt_test1a_dummy <- fastDummies::dummy_cols(
  wt_test1a, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'HA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

wt_test2a_dummy <- fastDummies::dummy_cols(
  wt_test2a, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'HA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

wt_test3a_dummy <- fastDummies::dummy_cols(
  wt_test3a, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'HA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

wt_test4a_dummy <- fastDummies::dummy_cols(
  wt_test4a, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'HA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

wt_test1b_dummy <- fastDummies::dummy_cols(
  wt_test1b, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'HA', 'NT_rep'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

wt_test2b_dummy <- fastDummies::dummy_cols(
  wt_test2b, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'HA', 'Lg_rep'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

wt_test3b_dummy <- fastDummies::dummy_cols(
  wt_test3b, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'HA', 'BnOB_rep'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

wt_test4b_dummy <- fastDummies::dummy_cols(
  wt_test4b, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'HA', 'Bn_rep'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

wt_test5_dummy <- fastDummies::dummy_cols(
  wt_test5, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'HA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

wt_test6_dummy <- fastDummies::dummy_cols(
  wt_test6, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'HA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

wt_test7a_dummy <- fastDummies::dummy_cols(
  wt_test7a, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'HA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

wt_test7b_dummy <- fastDummies::dummy_cols(
  wt_test7b, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'HA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

wt_test7c_dummy <- fastDummies::dummy_cols(
  wt_test7c, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'HA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

wt_test7d_dummy <- fastDummies::dummy_cols(
  wt_test7d, select_columns =
    c('HPAI_MBAA', 'RBS', 'PBS', 'HA'),
  remove_selected_columns = TRUE,
  remove_first_dummy = FALSE,
  ignore_na = TRUE)

wt_test_base_dummy$wt_loss_high <- as.factor(wt_test_base_dummy$wt_loss_high)
wt_test1a_dummy$wt_loss_high <- as.factor(wt_test1a_dummy$wt_loss_high)
wt_test1b_dummy$wt_loss_high <- as.factor(wt_test1b_dummy$wt_loss_high)
wt_test2a_dummy$wt_loss_high <- as.factor(wt_test2a_dummy$wt_loss_high)
wt_test2b_dummy$wt_loss_high <- as.factor(wt_test2b_dummy$wt_loss_high)
wt_test3a_dummy$wt_loss_high <- as.factor(wt_test3a_dummy$wt_loss_high)
wt_test3b_dummy$wt_loss_high <- as.factor(wt_test3b_dummy$wt_loss_high)
wt_test4a_dummy$wt_loss_high <- as.factor(wt_test4a_dummy$wt_loss_high)
wt_test4b_dummy$wt_loss_high <- as.factor(wt_test4b_dummy$wt_loss_high)
wt_test5_dummy$wt_loss_high <- as.factor(wt_test5_dummy$wt_loss_high)
wt_test6_dummy$wt_loss_high <- as.factor(wt_test6_dummy$wt_loss_high)
wt_test7a_dummy$wt_loss_high <- as.factor(wt_test7a_dummy$wt_loss_high)
wt_test7b_dummy$wt_loss_high <- as.factor(wt_test7b_dummy$wt_loss_high)
wt_test7c_dummy$wt_loss_high <- as.factor(wt_test7c_dummy$wt_loss_high)
wt_test7d_dummy$wt_loss_high <- as.factor(wt_test7d_dummy$wt_loss_high)

set.seed(9595)
wt_base_dumSplit <- rsample::initial_split(wt_test_base_dummy, prop = 0.70)
wt_1a_dumSplit <- rsample::initial_split(wt_test1a_dummy, prop = 0.70)
wt_1b_dumSplit <- rsample::initial_split(wt_test1b_dummy, prop = 0.70)
wt_2a_dumSplit <- rsample::initial_split(wt_test2a_dummy, prop = 0.70)
wt_2b_dumSplit <- rsample::initial_split(wt_test2b_dummy, prop = 0.70)
wt_3a_dumSplit <- rsample::initial_split(wt_test3a_dummy, prop = 0.70)
wt_3b_dumSplit <- rsample::initial_split(wt_test3b_dummy, prop = 0.70)
wt_4a_dumSplit <- rsample::initial_split(wt_test4a_dummy, prop = 0.70)
wt_4b_dumSplit <- rsample::initial_split(wt_test4b_dummy, prop = 0.70)
wt_5_dumSplit <- rsample::initial_split(wt_test5_dummy, prop = 0.70)
wt_6_dumSplit <- rsample::initial_split(wt_test6_dummy, prop = 0.70)
wt_7a_dumSplit <- rsample::initial_split(wt_test7a_dummy, prop = 0.70)
wt_7b_dumSplit <- rsample::initial_split(wt_test7b_dummy, prop = 0.70)
wt_7c_dumSplit <- rsample::initial_split(wt_test7c_dummy, prop = 0.70)
wt_7d_dumSplit <- rsample::initial_split(wt_test7d_dummy, prop = 0.70)

trainData_base <- rsample::training(wt_base_dumSplit)
testData_base <- rsample::testing(wt_base_dumSplit)
trainData_1a <- rsample::training(wt_1a_dumSplit)
testData_1a <- rsample::testing(wt_1a_dumSplit)
trainData_1b <- rsample::training(wt_1b_dumSplit)
testData_1b <- rsample::testing(wt_1b_dumSplit)
trainData_2a <- rsample::training(wt_2a_dumSplit)
testData_2a <- rsample::testing(wt_2a_dumSplit)
trainData_2b <- rsample::training(wt_2b_dumSplit)
testData_2b <- rsample::testing(wt_2b_dumSplit)
trainData_3a <- rsample::training(wt_3a_dumSplit)
testData_3a <- rsample::testing(wt_3a_dumSplit)
trainData_3b <- rsample::training(wt_3b_dumSplit)
testData_3b <- rsample::testing(wt_3b_dumSplit)
trainData_4a <- rsample::training(wt_4a_dumSplit)
testData_4a <- rsample::testing(wt_4a_dumSplit)
trainData_4b <- rsample::training(wt_4b_dumSplit)
testData_4b <- rsample::testing(wt_4b_dumSplit)
trainData_5 <- rsample::training(wt_5_dumSplit)
testData_5 <- rsample::testing(wt_5_dumSplit)
trainData_6 <- rsample::training(wt_6_dumSplit)
testData_6 <- rsample::testing(wt_6_dumSplit)
trainData_7a <- rsample::training(wt_7a_dumSplit)
testData_7a <- rsample::testing(wt_7a_dumSplit)
trainData_7b <- rsample::training(wt_7b_dumSplit)
testData_7b <- rsample::testing(wt_7b_dumSplit)
trainData_7c <- rsample::training(wt_7c_dumSplit)
testData_7c <- rsample::testing(wt_7c_dumSplit)
trainData_7d <- rsample::training(wt_7d_dumSplit)
testData_7d <- rsample::testing(wt_7d_dumSplit)

set.seed(2626)
mbase <- train(wt_loss_high ~ ., data = trainData_base,
             method = "gbm",
             na.action = na.exclude,
             preProcess = c("nzv", "scale", "center"),
             trControl = fitControl)

set.seed(2626)
m1a <- train(wt_loss_high ~ ., data = trainData_1a,
               method = "gbm",
               na.action = na.exclude,
               preProcess = c("nzv", "scale", "center"),
               trControl = fitControl)

set.seed(2626)
m1b <- train(wt_loss_high ~ ., data = trainData_1b,
             method = "gbm",
             na.action = na.exclude,
             preProcess = c("nzv", "scale", "center"),
             trControl = fitControl)

set.seed(2626)
m2a <- train(wt_loss_high ~ ., data = trainData_2a,
             method = "gbm",
             na.action = na.exclude,
             preProcess = c("nzv", "scale", "center"),
             trControl = fitControl)

set.seed(2626)
m2b <- train(wt_loss_high ~ ., data = trainData_2b,
             method = "gbm",
             na.action = na.exclude,
             preProcess = c("nzv", "scale", "center"),
             trControl = fitControl)

set.seed(2626)
m3a <- train(wt_loss_high ~ ., data = trainData_3a,
             method = "gbm",
             na.action = na.exclude,
             preProcess = c("nzv", "scale", "center"),
             trControl = fitControl)

set.seed(2626)
m3b <- train(wt_loss_high ~ ., data = trainData_3b,
             method = "gbm",
             na.action = na.exclude,
             preProcess = c("nzv", "scale", "center"),
             trControl = fitControl)

set.seed(2626)
m4a <- train(wt_loss_high ~ ., data = trainData_4a,
             method = "gbm",
             na.action = na.exclude,
             preProcess = c("nzv", "scale", "center"),
             trControl = fitControl)

set.seed(2626)
m4b <- train(wt_loss_high ~ ., data = trainData_4b,
             method = "gbm",
             na.action = na.exclude,
             preProcess = c("nzv", "scale", "center"),
             trControl = fitControl)

set.seed(2626)
m5 <- train(wt_loss_high ~ ., data = trainData_5,
             method = "gbm",
             na.action = na.exclude,
             preProcess = c("nzv", "scale", "center"),
             trControl = fitControl)

set.seed(2626)
m6 <- train(wt_loss_high ~ ., data = trainData_6,
            method = "gbm",
            na.action = na.exclude,
            preProcess = c("nzv", "scale", "center"),
            trControl = fitControl)

set.seed(2626)
m7a <- train(wt_loss_high ~ ., data = trainData_7a,
             method = "gbm",
             na.action = na.exclude,
             preProcess = c("nzv", "scale", "center"),
             trControl = fitControl)

set.seed(2626)
m7b <- train(wt_loss_high ~ ., data = trainData_7b,
             method = "gbm",
             na.action = na.exclude,
             preProcess = c("nzv", "scale", "center"),
             trControl = fitControl)

set.seed(2626)
m7c <- train(wt_loss_high ~ ., data = trainData_7c,
             method = "gbm",
             na.action = na.exclude,
             preProcess = c("nzv", "scale", "center"),
             trControl = fitControl)

set.seed(2626)
m7d <- train(wt_loss_high ~ ., data = trainData_7d,
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
probTest <- factor(ifelse(probTest$no >= 0.56, 'no', 'yes')) %>%
  as.data.frame()
probTruth <- testData %>%
  na.omit() %>%
  dplyr::select(wt_loss_high) %>%
  cbind(., probTest)
confusionMatrix(probTruth$., probTruth$wt_loss_high, mode = "everything")

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
    dplyr::select(wt_loss_high) %>%
    cbind(., probTest)
  
  confusion_matrix <- confusionMatrix(probTruth$., probTruth$wt_loss_high, mode = "everything")
  
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

## path only
## base
mcc_func(38, 12, 7, 13)
## NT
mcc_func(25, 8, 2, 16)
mcc_func(20, 9, 6, 18)
## Lg
mcc_func(26, 7, 1, 17)
mcc_func(28, 12, 6, 12)
## BnOB
mcc_func(25, 10, 0, 12)
mcc_func(24, 6, 5, 15)
## Bn
mcc_func(19, 11, 4, 16)
mcc_func(26, 13, 7, 14)
## All_avg
mcc_func(23, 8, 2, 14)
## All_avg2
mcc_func(26, 7, 2, 17)
##
mcc_func(28, 9, 5, 12)
mcc_func(29, 12, 2, 13)
mcc_func(22, 12, 3, 11)
mcc_func(31, 12, 7, 9)

## All
mcc_func(122, 37, 16, 29)
mcc_func(105, 8, 18, 16)
mcc_func(94, 33, 13, 29)
mcc_func(86, 30, 16, 30)
mcc_func(102, 21, 15, 25)
mcc_func(77, 38, 14, 31)
mcc_func(89, 23, 20, 37)
mcc_func(98, 26, 19, 33)
mcc_func(79, 37, 15, 40)
mcc_func(99, 25, 22, 27)
mcc_func(99, 29, 18, 25)

########################################################################################
########################################################################################
### End of Code
########################################################################################
########################################################################################
