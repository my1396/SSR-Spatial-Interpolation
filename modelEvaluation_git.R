##        Script Description      ##
# This script conducts 
        # 1. Models evaluation (various error metrics);

library(tidyverse)
library(ggplot2)
library(latex2exp)
library(gridExtra)
library(lubridate)
library(zoo)


## define command
'%ni%' <- Negate('%in%')
select <- dplyr::select
## define file path
figure_dir <- "./plot/"
data_dir <- "./data/"


f_name <- paste0(data_dir, con_folder, 
                 sprintf('%s_prediction_RF.csv', con_name))
f_name
continent_prediction_RF <- read_csv(f_name)

###   Error inspection   ### ----------------------------------
head(continent_prediction_RF)
continent_error <- continent_prediction_RF %>% select(7:15)
continent_error <- continent_error - matrix(c(rep(0, nrow(continent_error)), rep(continent_error$RAD, 8)), ncol=9)

head(continent_error)

error_inspection <- cbind(continent_prediction_RF[,1:6], continent_error)

# drop station c(2108,1314) in EU (suspicious records)
if (con_name=='EU'){
    error_inspection <- error_inspection %>% filter(ID %ni% c(2108,1314))    
}


cast <- apply(error_inspection, 1, function(x) sum(x %in% c(0, NA, NaN, Inf, -Inf)))
cast <- (cast == 0)
error_inspection <- error_inspection[cast,]


##    inspect MAE per month     ##   ---------------------------
MAE_df <- error_inspection %>% group_by(MON) %>% 
    summarise(RF=mean(abs(pred_RF)),
              OK=mean(abs(pred_ok)),
              GAM_OK=mean(abs(pred_gam_ok)),
              LSDV_OK=mean(abs(pred_lsdv_ok)))

error_inspection %>% group_by(MON) %>% 
    summarise(RF=mean(abs(pred_RF)),
              OK=mean(abs(pred_ok)),
              GAM_OK=mean(abs(pred_gam_ok)),
              LSDV_OK=mean(abs(pred_lsdv_ok))) %>% 
    colMeans()


MAE_df

long_format <- gather(MAE_df, key="Model", value="MAE", -MON)
long_format$MON <- factor(long_format$MON, as.character(seq(1, 12)))

head(long_format)
str(long_format$MON)
dev.new()
f_name <- paste0(figure_dir, "MAE_EU.png")
f_name
png(f_name, width=7*ppi, height=5*ppi, res=ppi)
ggplot(long_format, aes(x=MON, y=MAE, fill=Model)) +
    geom_bar(position="dodge", stat="identity") +
    ylab(TeX("MAE \\[$Wm^{-2}$\\]")) +
    theme_minimal()
dev.off()

##      MAPE     ##   ---------------------------
MAPE_df <- error_inspection %>% group_by(MON) %>% 
    summarise(RF=mean(abs(pred_RF/RAD)),
              OK=mean(abs(pred_ok/RAD)),
              GAM_OK=mean(abs(pred_gam_ok/RAD)),
              LSDV_OK=mean(abs(pred_lsdv_ok/RAD)))

error_inspection %>% group_by(MON) %>% 
    summarise(RF=mean(abs(pred_RF/RAD)),
              OK=mean(abs(pred_ok/RAD)),
              GAM_OK=mean(abs(pred_gam_ok/RAD)),
              LSDV_OK=mean(abs(pred_lsdv_ok/RAD))) %>%
    colMeans()

abs(error_inspection$pred_ok/error_inspection$RAD) %>% mean()
abs(error_inspection$pred_gam_ok/error_inspection$RAD) %>% mean()

long_format <- gather(MAPE_df, key="Model", value="MAPE", -MON)
long_format$MON <- factor(long_format$MON, as.character(seq(1, 12)))
head(long_format)

dev.new()
f_name <- paste0(figure_dir, "MAPE_EU.png")
f_name
png(f_name, width=7*ppi, height=5*ppi, res=ppi)
ggplot(long_format, aes(x=MON, y=MAPE*100, fill=Model)) +
    geom_bar(position="dodge", stat="identity") +
    ylab("MAPE %") +
    theme_minimal()
dev.off()

##    RMSE      ##   ---------------------------
rmse <- function(y_resid){
    (y_resid)^2 %>% mean() %>% sqrt()
}


RMSE_df <- error_inspection %>% group_by(MON) %>% 
    summarise(RF=rmse(pred_RF),
              OK=rmse(pred_ok),
              GAM_OK=rmse(pred_gam_ok),
              LSDV_OK=rmse(pred_lsdv_ok) )

long_format <- gather(RMSE_df, key="Model", value="RMSE", -MON)
long_format$MON <- factor(long_format$MON, as.character(seq(1, 12)))
head(long_format)

dev.new()
f_name <- paste0(figure_dir, "RMSE_EU.png")
f_name
png(f_name, width=7*ppi, height=5*ppi, res=ppi)
ggplot(long_format, aes(x=MON, y=RMSE, fill=Model)) +
    geom_bar(position="dodge", stat="identity") +
    ylab("RMSE") +
    theme_minimal()
dev.off()


##    RMSPE     ##   ---------------------------
rmspe <- function(y_resid, obs){
    (y_resid/obs)^2 %>% mean() %>% sqrt()
}


RMSPE_df <- error_inspection %>% group_by(MON) %>% 
    summarise(RF=rmspe(pred_RF, RAD),
              OK=rmspe(pred_ok, RAD),
              GAM_OK=rmspe(pred_gam_ok, RAD),
              LSDV_OK=rmspe(pred_lsdv_ok, RAD) )

long_format <- gather(RMSPE_df, key="Model", value="RMSPE", -MON)
long_format$MON <- factor(long_format$MON, as.character(seq(1, 12)))
head(long_format)

dev.new()
f_name <- paste0(figure_dir, "RMSPE_EU.png")
f_name
png(f_name, width=7*ppi, height=5*ppi, res=ppi)
ggplot(long_format, aes(x=MON, y=RMSPE*100, fill=Model)) +
    geom_bar(position="dodge", stat="identity") +
    ylab("RMSPE %") +
    theme_minimal()
dev.off()

##   R-quared     ##   ---------------------------
r_squared <- function(y_resid, obs){
    1- (y_resid^2 %>% sum()) / ((y_resid-mean(obs))^2 %>% sum())
}

R_squared_df <- error_inspection %>% group_by(MON) %>% 
    summarise(RF=r_squared(pred_RF, RAD),
              OK=r_squared(pred_ok, RAD),
              GAM_OK=r_squared(pred_gam_ok, RAD),
              LSDV_OK=r_squared(pred_lsdv_ok, RAD) )


long_format <- gather(R_squared_df, key="Model", value="R_squared", -MON)
long_format$MON <- factor(long_format$MON, as.character(seq(1, 12)))
head(long_format)

dev.new()
f_name <- paste0(figure_dir, "R_squared_EU.png")
f_name
png(f_name, width=7*ppi, height=5*ppi, res=ppi)
ggplot(long_format, aes(x=MON, y=R_squared, fill=Model)) +
    geom_bar(position="dodge", stat="identity") +
    coord_cartesian(ylim=c(0.9,1)) +
    ylab(TeX("$R^2$")) +
    theme_minimal()
dev.off()



NaRV.omit <- function(x){
    # not a Regular Value
    cast <- x %in% c(NA, NaN, Inf, -Inf)
    x[!cast]
}

Performance_stats <- function(x){
    obs <- x[,1]
    pred <- x[,-1]
    apply(pred, 2, function(col) {
        R_squared <- 1- (col^2 %>% sum()) / ((obs-mean(obs))^2 %>% sum())
        MAE <- abs(col) %>% mean()
        MAPE <- abs(col/obs) %>% NaRV.omit() %>% mean()
        RMSE <- (col)^2 %>% NaRV.omit() %>% mean() %>% sqrt()
        RMSPE <- (col/obs)^2 %>% NaRV.omit() %>% mean() %>% sqrt()
        c("R_squared"=R_squared, "MAE"=MAE, "MAPE"=MAPE, "RMSE"=RMSE, "RMSPE"=RMSPE)
        })
    }


##   Evaluation metrics for prediction error  ## --------------------
evaluation_mat <- error_inspection[,7:15] %>% Performance_stats()
evaluation_mat %>% round(4) 


f_name <- paste0(data_dir, con_folder, 
                 sprintf('eval_%s.csv', con_name))
f_name
# evaluation_mat %>% round(2) %>% as.data.frame() %>% write.csv(f_name)













