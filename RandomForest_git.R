##        Script Description      ##
# This script conducts 
        # 1. Random Forest prediction;
        # 2. Merging all models prediction all together;

library(tidyverse)
library(ranger)


## define command
'%ni%' <- Negate('%in%')
select <- dplyr::select
## define file path
figure_dir <- "./plot/"
data_dir <- "./data/"

## Cross-Validation 10 folds 
# Outputs a list with test (or train) indices
kfoldSplit <- function(x, k=10, train=TRUE){
    x <- sample(x, size = length(x), replace = FALSE)
    out <- suppressWarnings(split(x, factor(1:k)))
    if(train) out <- lapply(out, FUN = function(x, len) (1:len)[-x], len=length(unlist(out)))
    return(out)
}



rsds_data <- read_csv("./data/extended_data_sample.csv", 
                      na = c(""))
head(rsds_data)

auxilary_v <- c('CLD','DTR','FRS', 'PRE', 'TMN', 'TMP', 'TMX', 'VAP', 'WET','URB','LAT', 'LON','MON','ALT')
neighbor_v <- c('z_k1', 'zw_k1', 'zsig_k1', 'z_k2', 'zw_k2', 'zsig_k2','z_k3', 'zw_k3', 'zsig_k3', 'z_k1k2', 'z_k2k3', 'zw_k1k2', 'zw_k2k3','RAD_tm1', 'RAD_t1', 'RAD_t2', 'RAD_t3', 'RAD_mam2', 'RAD_mam3',  'RAD_may2', 'RAD_may3')
regressors <- c(auxilary_v, neighbor_v)

###    Random Forest Prediction ### --------------------------------
con_name <- 'EU'
con_folder <- 'Europe/'

## Continent Time Period Coverage TC
TC <- switch(con_name, 
             'EU' = c(1967, 2013),
             'AF' = c(1964, 1998),
             'NrA' = c(1964, 2006),
             'SA' = c(1967, 1998),
             'AS' = c(1964, 2010),
             'OC' = c(1971, 1991),
             c(1964, 2013) # default case
)

EU <- rsds_data %>% 
    filter(CON==con_name & between(YR, TC[1], TC[2])) 
    

# For North America
con_name <- 'NrA'
con_folder <- 'NorthAmerica/'
cast <- is.na(rsds_data['CON']) %>% which()
EU <- rsds_data[cast,] 
EU <- EU %>% filter(LON<130 & between(YR, 1964, 2006)) 
nrow(EU)
distinct(EU, ID, .keep_all = T)


## Random Forest 10-fold prediction ----------------------------------
RF_validation <- matrix(NA, ncol=8)
colnames(RF_validation) <- c('X1', 'ID', 'YR', 'MON', 
                             'LAT', 'LON', "RAD", "pred_RF")

set.seed(123)
kfolds_EU <- kfoldSplit(1:nrow(EU), k = 10, train = TRUE)


for (i in 1:10){
    print (i)
    idx <- kfolds_EU[[i]]
    
    RF_ranger <- ranger(as.formula(
        paste("RAD", 
              paste(c(auxilary_v, neighbor_v), collapse = " + "), 
              sep = " ~ ")), data=EU[idx, ])
    rf.pred.test <- predict(RF_ranger, data=EU[-idx, c(auxilary_v, neighbor_v)], num.trees=700, mtry=6)$predictions
    
    obs.test <- EU[-idx,'RAD'] %>% pull()
    
    new_validation <- cbind(EU[-idx,c('X1', 'ID', 'YR', 'MON', 'LAT', 'LON')], 
                            obs.test, rf.pred.test)
    new_validation <- as.matrix(new_validation)
    RF_validation <- rbind(RF_validation, new_validation)
}



RF_validation <- RF_validation[-1,]
RF_validation <- as.tibble(RF_validation) %>% arrange(ID,YR,MON)
head(RF_validation)
dim(RF_validation)

# f_name <- paste0(data_dir, con_folder, 
#                  sprintf('RF_validation_%s.csv', con_name) )
# f_name
# write_csv(RF_validation, f_name)



## read in Regression Kriging predictions      ##
## and merge with RF prediction                ##
f_name <- paste0(data_dir, con_folder, 
                 sprintf('%s_prediction.csv', con_name) )
f_name
continent_prediction <- read_csv(f_name)
head(continent_prediction)
nrow(continent_prediction)
unique(continent_prediction$X1) %>% length()

continent_prediction_RF <- left_join(RF_validation, continent_prediction[,c(1, 8:14)], by='X1')


f_name <- paste0(data_dir, con_folder, 
                 sprintf('%s_prediction_RF.csv', con_name))

f_name
# write_csv(continent_prediction_RF, f_name)

