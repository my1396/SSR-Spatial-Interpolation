## Script Description
# This file performs Regression Kriging.
# Logistics: 1. Regression models to predict trends;
#            2. Kriging to predict spatial correlated residuals;
#            3. Combinding trends with residuals to obtain final
#               predictions


library(tidyverse)
library(latex2exp)
library(sp)
library(raster)
library(gstat)
library(randomForest)
library(mgcv)
library(ggplot2)

# define custom commands
select <- dplyr::select
'%ni%' <- Negate('%in%')
# define file path
figure_dir <- "./plot/"
data_dir <- "./data/"

## 10-fold Cross-Validation data partition function
kfoldSplit <- function(x, k=10, train=TRUE){
    # Outputs a list with test (or train) indices
    x <- sample(x, size = length(x), replace = FALSE)
    out <- suppressWarnings(split(x, factor(1:k)))
    if(train) out <- lapply(out, FUN = function(x, len) (1:len)[-x], len=length(unlist(out)))
    return(out)
}


##  program starts       ----------------------------------
rsds_data <- read_csv("./data/extended_data.csv", 
                      na = c(""))
head(rsds_data)


# unique stations in GEBA
stations_all <- distinct(rsds_data, ID, .keep_all = T)
nrow(stations_all) # 1227

auxilary_v <- c('CLD','DTR','FRS', 'PRE', 'TMN', 'TMP', 'TMX', 'VAP', 'WET','URB','LAT', 'LON','MON','ALT')
neighbor_v <- c('z_k3', 'RAD_may3', 'RAD_tm1', 'zw_k3', 'RAD_may2', 'RAD_t1', 'z_k2')


##  subset continent data ----------------------------------
# Continent other than NrA
con_name <- 'EU'
label <- 'Europe'
con_folder <- 'Europe/'

# Subset period 1964-2013 (50)
EU <- rsds_data %>% 
  filter(CON==con_name) %>% filter(YR>=1964)

# for North America
con_name <- 'NrA'
con_folder <- 'NorthAmerica/'
EU <- rsds_data %>% 
    filter(is.na(CON)) %>% filter(YR>=1964)
EU <- EU %>% filter(LON<130) 


##  check No. stations per year   ##  -----------------------------------
sprintf("No. of Years: %s", EU$YR %>% unique() %>% length())  # No. Years
EU$YR %>% unique() %>% sort()    # YR sequence
EU$ID %>% unique() %>% length()  # Total No. Stations across time

EU_yr <- EU %>% group_by(YR) %>% group_split()

station_count <- EU %>% group_by(YR) %>% 
  summarise( ID_count = unique(ID)%>%length() )

head(station_count)


cast <- seq(1964,2013) %ni% station_count$YR
if (sum(cast)!=0){
  # fill no observation year with count zero
  station_count <- station_count %>% 
    add_row(YR=seq(1964,2013)[cast], 
            ID_count=rep(0,sum(cast)))
}



dev.new()
p <- ggplot(data=station_count, aes(x=YR, y=ID_count)) +
  geom_bar(stat="identity") +
  labs(y="# of stations") + 
  annotate(geom = 'label', label = label, x = Inf, y = Inf, hjust = 2, vjust = 2, parse=TRUE) +
  theme_bw() +
  theme(axis.title.x=element_blank())

p

f_name <- paste0(data_dir, con_folder, 
                 sprintf("No_observation_%s.png", con_name))
f_name
png(f_name, width=7*ppi, height=5*ppi, res=ppi)
p
dev.off()


## remove years with too little stations
head(station_count)
tail(station_count, 20)
## Continent Time Period Coverage TC
TC <- switch(con_name, 
             'EU' = c(1967, 2013),
             'AF' = c(1964, 1998),
             'NrA' = c(1964, 2006),
             'SA' = c(1965, 1998),
             'AS' = c(1964, 2010),
             'OC' = c(1971, 1991),
             c(1964, 2013) # default case
             )
EU <- EU %>% filter(between(YR, TC[1], TC[2]))
nrow(EU)


## Exploratory data analysis -----------------------------
radiation <- EU[,c('YR','MON','RAD')]
radiation <- radiation %>% group_by(YR) %>% summarise(rad=mean(RAD))

f_name <- paste0(figure_dir, "MON_boxplot.png")
f_name
png(f_name, width=7*ppi, height=5*ppi, res=ppi)
par(oma=c(0.5,0.5,0,2), mar=c(3, 4, 1, 0), mgp=c(2,0.8,0))
dev.new()
boxplot(RAD~MON, data= EU_subset, xlab="MON", ylab=TeX("RAD Wm^{-2}"))
dev.off()



###   Linear regression to predict trends [ten-fold] ### -------------------
set.seed(123)
kfolds_EU <- kfoldSplit(1:nrow(EU), k = 10, train = TRUE)

# define regression formula
f_ols <- as.formula(
  paste("RAD", 
        paste(c('CLD','DTR','FRS', 'PRE', 'TMP', 'WET', 'LAT', 'LON', 'ALT','URB','MON','YR'), collapse = " + "), 
        sep = " ~ "))

f_gam <- as.formula(
  paste("RAD", 
        paste(c('CLD','DTR','FRS', 'PRE', 'TMP', 'WET', 'LAT', 'LON', 'ALT','URB','s(MON)','YR'), collapse = " + "), 
        sep = " ~ "))

f_lsdv <- as.formula(
  paste("RAD", 
         paste(c('CLD','DTR','FRS', 'PRE', 'TMP', 'WET', 'LAT', 'LON', 'ALT','URB','factor(MON)','YR'), collapse = " + "), 
         sep = " ~ "))

Regression_validation <- matrix(NA, ncol=10)
colnames(Regression_validation) <- c('X1','ID', 'YR', 'MON', "LAT", "LON", "RAD", "pred_ols", "pred_gam", "pred_lsdv")

for (i in 1:10){
    print (i)
    idx <- kfolds_EU[[i]]
    obs.test <- EU[-idx,'RAD'] %>% pull()
    
    OLS <- lm(f_ols, data=EU[idx, ])
    ols.pred.test <- predict(OLS, newdata=EU[-idx,])

    GAM <- gam(f_gam,  data=EU[idx,])
    gam.pred.test <- predict(GAM, newdata=EU[-idx,])
    
    LSDV <- lm(f_lsdv, data=EU[idx, ])
    lsdv.pred.test <- predict(LSDV, newdata=EU[-idx,])
    
    new_validation <- cbind(EU[-idx,c('X1','ID', 'YR', 'MON', 'LAT', 'LON')],  obs.test, ols.pred.test, gam.pred.test, lsdv.pred.test)
    new_validation <- as.matrix(new_validation)
    Regression_validation <- rbind(Regression_validation, new_validation)
    
}


Regression_validation <- Regression_validation[-1,]
head(Regression_validation,5)
dim(Regression_validation)


EU_mon <- as.tibble(Regression_validation) %>% 
                mutate(resid_ols=RAD-pred_ols,
                       resid_gam=RAD-pred_gam,
                       resid_lsdv=RAD-pred_lsdv)
EU_mon <- EU_mon %>% arrange(YR, MON) 
months <- EU_mon %>% group_by(YR, MON) %>% group_keys()
months["idx"] <- seq(1,nrow(months))
months

EU_mon <- EU_mon %>% group_by(YR, MON) %>% group_split

##  Kriging to predict spatial correlation [ten-fold] ## ---------------------------

## define Variogram formula
formMod.ok <- RAD ~ 1
formMod.ols.ok <- resid_ols ~ 1
formMod.gam.ok <- resid_gam ~ 1
formMod.lsdv.ok <- resid_lsdv ~ 1
proj4Str <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
mod <- vgm(model=c("Exp", "Mat", "Sph"))

Regression_validation2 <- matrix(NA,ncol=14) # 4 more columns to include regression trends + OK residuals
set.seed(123)

for (i in 1:nrow(months)){
  print (i)
  
  rk <- EU_mon[[i]]
  
  if (nrow(rk)<=10) {
    next
  }
  
  # longitude convertion from 0~360 to -180~180 
  rk$LON <- rk$LON-180 
  # remove duplication locations
  rk <- distinct(rk, LAT, LON, .keep_all = TRUE)
  statPoints <- SpatialPointsDataFrame(coords=rk[,c("LON","LAT")],
                                       data=rk,
                                       proj4string = CRS(proj4Str))
  
  # Kriging for each month
  kfolds_mon <- kfoldSplit(1:nrow(rk), k=10, train=TRUE)
  
  for (j in 1:10){
  to_train <- kfolds_mon[[j]]
  to_train_bool <-  1:nrow(rk) %in% to_train                                        
  
  # Ordinary kriging predictions
  variog <- variogram(formMod.ok, statPoints)
  variogFitOLS<-fit.variogram(variog, model = mod)
  
  if (variogFitOLS$range>0){
    OK <- krige(formula = formMod.ok ,
                locations = statPoints[to_train_bool, ], 
                model = variogFitOLS,
                newdata = statPoints[!to_train_bool, ],
                debug.level = 0)
    ok.pred.test <- OK@data$var1.pred
  } else {ok.pred.test <- NA}
  
  
  
  # OLS-kriging predictions
  variog <- variogram(formMod.ols.ok, statPoints)
  variogFitOLS.ok<-fit.variogram(variog, model = mod)
  
  if (variogFitOLS.ok$range>0){
    OK <- krige(formula = formMod.ols.ok ,
                locations = statPoints[to_train_bool, ], 
                model = variogFitOLS.ok,
                newdata = statPoints[!to_train_bool, ],
                debug.level = 0)
    ok.ols.pred.test <- OK@data$var1.pred+pull(rk[!to_train_bool, 'pred_ols'])
  } else {ok.ols.pred.test <- pull(rk[!to_train_bool, 'pred_ols'])}
  
  
  # GAM-kriging predictions
  variog <- variogram(formMod.gam.ok, statPoints)
  variogFitgam.ok<-fit.variogram(variog, model = mod)
  
  if (variogFitgam.ok$range>0){
    OK <- krige(formula = formMod.gam.ok ,
                locations = statPoints[to_train_bool, ], 
                model = variogFitgam.ok,
                newdata = statPoints[!to_train_bool, ],
                debug.level = 0)
    ok.gam.pred.test <- OK@data$var1.pred + pull(rk[!to_train_bool, 'pred_gam'])
  } else {ok.gam.pred.test <- pull(rk[!to_train_bool, 'pred_gam'])}
  
  
  # LSDV-kriging predictions
  variog <- variogram(formMod.lsdv.ok, statPoints)
  variogFitlsdv.ok<-fit.variogram(variog, model = mod)
  
  if (variogFitlsdv.ok$range>0){
    OK <- krige(formula = formMod.lsdv.ok ,
                locations = statPoints[to_train_bool, ], 
                model = variogFitlsdv.ok,
                newdata = statPoints[!to_train_bool, ],
                debug.level = 0)
    ok.lsdv.pred.test <- OK@data$var1.pred + pull(rk[!to_train_bool, 'pred_lsdv'])
  } else {ok.lsdv.pred.test <- pull(rk[!to_train_bool, 'pred_lsdv'])}
  
  
  
  
  prediction <- rk[!to_train_bool, 1:10]
  prediction <- prediction %>% mutate(pred_ok = ok.pred.test,
                                      pred_ols_ok = ok.ols.pred.test,
                                      pred_gam_ok = ok.gam.pred.test,
                                      pred_lsdv_ok = ok.lsdv.pred.test)
  
  Regression_validation2 <- rbind(Regression_validation2, as.matrix(prediction))
  } # end of 10 folds prediction for each month
  
}


Regression_validation2 <- Regression_validation2[-1,]
head(Regression_validation2, 10)
tail(Regression_validation2, 20)
dim(Regression_validation2)

f_name <- paste0(data_dir, con_folder, 
                 sprintf('%s_prediction.csv', con_name))
f_name
write_csv(as.data.frame(Regression_validation2), f_name)

