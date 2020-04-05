##        Script Description      ##
# This script conducts 
        # 1. Trend inspection and visualization;
        # 2. Simulation v.s. Observation scatter plot with R^2 reported

library(tidyverse)
library(sp)
library(raster)
library(gstat)
library(ranger)
library(mgcv)
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


monthly_anomaliy <- function(the_month, all=1){
    # Calculate monthly anomalies by 
    # subtracting original value by corresponding monthly mean
    if (all==1){
        anomaly <- the_month %>% 
            select("RAD","pred_RF","pred_ok","pred_gam_ok","pred_lsdv_ok")
        anomaly_2 <- anomaly-matrix(rep(colMeans(anomaly,na.rm=T), 
                                        nrow(anomaly)), ncol=5,  byrow=TRUE)
    } else if (all==0) {
        anomaly <- the_month %>% select("RAD","pred_RF") 
        anomaly_2 <- anomaly-matrix(rep(colMeans(anomaly,na.rm=T), 
                                        nrow(anomaly)), ncol=2,  byrow=TRUE)
    }
    
    cbind(the_month[,'date'], anomaly_2)
}

con_name <- 'EU'
con_folder <- 'Europe/'

f_name <- paste0(data_dir, con_folder, 
                 sprintf('%s_prediction_RF.csv', con_name))
f_name

continent_prediction_RF <- read_csv(f_name)



# observations for each station
ID_count <- continent_prediction_RF %>% group_by(ID) %>% summarise(No=n())
ID_count %>% arrange(desc(No))

head(continent_prediction_RF)
nrow(continent_prediction_RF)
continent_prediction_RF$ID %>% unique() %>% length()



##     Visulization            ##  -----------------------------------------
##      1. Monthly Anomalieis       ##  -------------------------
##        1.1 Station anomaliy plot   ## -----------------
station <- continent_prediction_RF %>% filter(ID==1188)
station <- station %>% arrange(YR, MON)
station$day <- 1
station$date <- with(station, ymd(sprintf('%04d%02d%02d', YR, MON, day)))

head(station)

station_mon <- station %>% group_by(MON) %>% group_split()
station_ano <- station_mon %>% map_dfr(monthly_anomaliy)
station_ano <- station_ano[order(station_ano$date),]
cast <- is.na(station_ano$pred_ok)
station_ano[cast, 'pred_ok'] <- station_ano[cast, 'RAD']
n <- nrow(station_ano)

station_ano <- station_ano %>% mutate(
    ## Gaussian kernel smoother: ksmooth with 12 months BW ##
    smooth_RAD=ksmooth(1:n, RAD, kernel="normal", bandwidth=6)$y,
    smooth_RF=ksmooth(1:n, pred_RF, kernel="normal", bandwidth=6)$y,
    smooth_ok=ksmooth(1:n, pred_ok, kernel="normal", bandwidth=6)$y,
    smooth_gam_ok=ksmooth(1:n, pred_gam_ok, kernel="normal", bandwidth=6)$y,
    smooth_lsdv_ok=ksmooth(1:n, pred_lsdv_ok, kernel="normal", bandwidth=6)$y
    ) # end mutate

head(station_ano)


long_format <- gather(station_ano, key="Model", value="prediction", -date)
nrow(long_format)
head(long_format)

p1 <- ggplot(long_format %>% filter(Model %in% c('RAD','pred_RF','smooth_RAD','smooth_RF')) %>% 
                 mutate(Model=factor(Model, c('RAD','pred_RF','smooth_RAD','smooth_RF'))), aes(x=date, y=prediction, color = Model, linetype = Model, size = Model)) +
    geom_line() +
    ylab(TeX("SSR \\[$Wm^{-2}$\\]")) +
    # xlab('') +
    labs(subtitle="Panel A") +
    guides(col=guide_legend(nrow = 1),
           linetype=guide_legend(nrow = 1)) +
    theme_bw() + # theme_bw() will overwirte all theme() settings
    theme(
        axis.title.x=element_blank(),
        plot.subtitle = element_text(size=12, vjust=-.05),
        axis.title.y = element_text(size = 12), 
        legend.justification=c(1,0),
        legend.position = c(1, 0),
        legend.background = element_blank(),
        legend.title=element_blank(),
        plot.margin = unit(c(5.5, 5.5, 0, 5.5), "pt")
    ) +
    scale_color_manual(values=c("RAD" = "red", "pred_RF" = "blue", "smooth_RAD" = "red", "smooth_RF" = "blue")) +
    scale_linetype_manual(values=c("RAD" = "dashed", "pred_RF" = "dashed", "smooth_RAD" = "solid", "smooth_RF" = "solid")) +
    scale_size_manual(values=c("RAD" = 0.3, "pred_RF" = 0.3, "smooth_RAD" = 0.7, "smooth_RF" = 0.7))

p2 <- ggplot(long_format %>% filter(Model %in% c('RAD','pred_ok','smooth_RAD','smooth_ok')) %>% 
                 mutate(Model=factor(Model, c('RAD','pred_ok','smooth_RAD','smooth_ok'))), aes(x=date, y=prediction, group=Model, size = Model)) +
    geom_line(aes(color = Model, linetype = Model)) +
    ylab(TeX("SSR \\[$Wm^{-2}$\\]")) +
    # xlab('') +
    labs(subtitle="Panel B") +
    guides(col=guide_legend(nrow = 1),
           linetype=guide_legend(nrow = 1)) +
    theme_bw() + # theme_bw() will overwirte all theme() settings
    theme(
        axis.title.x=element_blank(),
        plot.subtitle = element_text(size=12, vjust=-.05),
        axis.title.y = element_text(size = 12), 
        legend.justification=c(1,0),
        legend.position = c(1, 0),
        legend.background = element_blank(),
        legend.title=element_blank(),
        plot.margin = unit(c(3.5, 5.5, 3.5, 5.5), "pt")
    ) +
    scale_color_manual(values=c("RAD" = "red", "pred_ok" = "blue", "smooth_RAD" = "red", "smooth_ok" = "blue")) +
    scale_linetype_manual(values=c("RAD" = "dashed", "pred_ok" = "dashed", "smooth_RAD" = "solid", "smooth_ok" = "solid")) +
    scale_size_manual(values=c("RAD" = 0.3, "pred_ok" = 0.3, "smooth_RAD" = 0.7, "smooth_ok" = 0.7))


p3 <- ggplot(long_format %>% filter(Model %in% c('RAD','pred_gam_ok', 'smooth_RAD', 'smooth_gam_ok')) %>% 
                 mutate(Model=factor(Model, c('RAD','pred_gam_ok', 'smooth_RAD', 'smooth_gam_ok'))), aes(x=date, y=prediction, color = Model, linetype = Model, size = Model)) +
    geom_line() +
    ylab(TeX("SSR \\[$Wm^{-2}$\\]")) +
    xlab('Date') +
    labs(subtitle="Panel C") +
    guides(col=guide_legend(nrow = 1),
           linetype=guide_legend(nrow = 1)) +
    theme_bw() + # theme_bw() will overwirte all theme() settings
    theme(
        plot.subtitle = element_text(size=12, vjust=-.05),
        axis.title.y = element_text(size = 12), 
        legend.justification=c(1,0),
        legend.position = c(1, 0),
        legend.background = element_blank(),
        legend.title=element_blank(),
        plot.margin = unit(c(0, 5.5, 5.5, 5.5), "pt")
    ) +
    scale_color_manual(values=c("RAD" = "red", "pred_gam_ok" = "blue", "smooth_RAD" = "red", "smooth_gam_ok" = "blue")) +
    scale_linetype_manual(values=c("RAD" = "dashed", "pred_gam_ok" = "dashed", "smooth_RAD" = "solid", "smooth_gam_ok" = "solid")) +
    scale_size_manual(values=c("RAD" = 0.3, "pred_gam_ok" = 0.3, "smooth_RAD" = 0.7, "smooth_gam_ok" = 0.7))


p4 <- ggplot(long_format %>% 
                 filter(Model %in% c('RAD','pred_lsdv_ok', 'smooth_RAD', 'smooth_lsdv_ok')) %>% 
                 mutate(Model=factor(Model, c('RAD','pred_lsdv_ok', 'smooth_RAD', 'smooth_lsdv_ok'))), 
             aes(x=date, y=prediction, color = Model, linetype = Model, size = Model)) +
    geom_line() +
    ylab(TeX("SSR \\[$Wm^{-2}$\\]")) +
    xlab('Date') +
    guides(col=guide_legend(nrow = 1),
           linetype=guide_legend(nrow = 1)) +
    theme_bw() + # theme_bw() will overwirte all theme() settings
    theme(
        plot.subtitle = element_text(size=12, vjust=-.05),
        axis.title.y = element_text(size = 12), 
        legend.justification=c(1,0),
        legend.position = c(1, 0),
        legend.background = element_blank(),
        legend.title=element_blank(),
        plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt")
    ) +
    scale_color_manual(values=c("RAD" = "red", "pred_lsdv_ok" = "blue", "smooth_RAD" = "red", "smooth_lsdv_ok" = "blue")) +
    scale_linetype_manual(values=c("RAD" = "dashed", "pred_lsdv_ok" = "dashed", "smooth_RAD" = "solid", "smooth_lsdv_ok" = "solid")) +
    scale_size_manual(values=c("RAD" = 0.3, "pred_lsdv_ok" = 0.3, "smooth_RAD" = 0.7, "smooth_lsdv_ok" = 0.7))



f_name <- paste0(figure_dir, "station_1188_lsdv.png")
f_name
png(f_name, width=10*ppi, height=8*ppi, res=ppi)
dev.new()
p4
dev.off()

f_name <- paste0(figure_dir, "station_1188.png")
f_name
png(f_name, width=10*ppi, height=8*ppi, res=ppi)

dev.new()
grid.arrange(p1, p2, p3, ncol=1)
dev.off()



##     1.2 Continent level anomalies plot       ##  ----------------
station <- continent_prediction_RF

station$day <- 1
station$date <- with(station, 
                     ymd(sprintf('%04d%02d%02d', YR, MON, day)))

# aggregate to continent level
station <- station %>% group_by(date) %>% 
        summarise(RAD=mean(RAD),
              pred_RF=mean(pred_RF),
              pred_ok=mean(pred_ok, na.rm=T),
              pred_gam_ok=mean(pred_gam_ok, na.rm=T),
              pred_lsdv_ok=mean(pred_lsdv_ok, na.rm=T)
              )


##     Monthly Anomalies   ##  -----------
station$MON <- station$date %>% month()

# group by MON and subtract respective monthly mean --> Monthly Anomalies
station_mon <- station %>% group_by(MON) %>% group_split()

# station_ano <- station_mon %>% map_dfr(monthly_anomaliy, all=0)
station_ano <- station_mon %>% map_dfr(monthly_anomaliy, all=1)
station_ano <- station_ano[order(station_ano$date),]
n <- nrow(station_ano)


cast <- is.na(station_ano)
cast <- is.na(station_ano$pred_ok)
station_ano[cast, 'pred_ok'] <- station_ano[cast, 'RAD']

station_ano <- station_ano %>% mutate(
    ## Gaussian kernel smoother: ksmooth with 12 months BW ##
    smooth_RAD=ksmooth(1:n, RAD, kernel="normal", bandwidth=6)$y,
    smooth_RF=ksmooth(1:n, pred_RF, kernel="normal", bandwidth=6)$y,
    smooth_ok=ksmooth(1:n, pred_ok, kernel="normal", bandwidth=6)$y,
    smooth_gam_ok=ksmooth(1:n, pred_gam_ok, kernel="normal", bandwidth=6)$y,
    smooth_lsdv_ok=ksmooth(1:n, pred_lsdv_ok, kernel="normal", bandwidth=6)$y
) # end mutate



long_format <- gather(station_ano, key="Model", value="prediction", -date)
nrow(long_format)
head(long_format)

p1 <- ggplot(long_format %>% filter(Model %in% c('RAD','pred_RF','smooth_RAD','smooth_RF')) %>% 
                 mutate(Model=factor(Model, c('RAD','pred_RF','smooth_RAD','smooth_RF'))), aes(x=date, y=prediction, color = Model, linetype = Model, size = Model)) +
    geom_line() +
    ylab(TeX("SSR \\[$Wm^{-2}$\\]")) +
    # xlab('') +
    labs(subtitle="Panel A") +
    guides(col=guide_legend(nrow = 1),
           linetype=guide_legend(nrow = 1)) +
    theme_bw() + # theme_bw() will overwirte all theme() settings
    theme(
        axis.title.x=element_blank(),
        plot.subtitle = element_text(size=12, vjust=-.05),
        axis.title.y = element_text(size = 12), 
        legend.justification=c(1,0),
        legend.position = c(1, 0),
        legend.background = element_blank(),
        legend.title=element_blank(),
        plot.margin = unit(c(5.5, 5.5, 0, 5.5), "pt")
    ) +
    scale_color_manual(values=c("RAD" = "red", "pred_RF" = "blue", "smooth_RAD" = "red", "smooth_RF" = "blue")) +
    scale_linetype_manual(values=c("RAD" = "dashed", "pred_RF" = "dashed", "smooth_RAD" = "solid", "smooth_RF" = "solid")) +
    scale_size_manual(values=c("RAD" = 0.3, "pred_RF" = 0.3, "smooth_RAD" = 0.7, "smooth_RF" = 0.7))

p2 <- ggplot(long_format %>% filter(Model %in% c('RAD','pred_ok','smooth_RAD','smooth_ok')) %>% 
                 mutate(Model=factor(Model, c('RAD','pred_ok','smooth_RAD','smooth_ok'))), aes(x=date, y=prediction, group=Model, size = Model)) +
    geom_line(aes(color = Model, linetype = Model)) +
    ylab(TeX("SSR \\[$Wm^{-2}$\\]")) +
    # xlab('') +
    labs(subtitle="Panel B") +
    guides(col=guide_legend(nrow = 1),
           linetype=guide_legend(nrow = 1)) +
    theme_bw() + # theme_bw() will overwirte all theme() settings
    theme(
        axis.title.x=element_blank(),
        plot.subtitle = element_text(size=12, vjust=-.05),
        axis.title.y = element_text(size = 12), 
        legend.justification=c(1,0),
        legend.position = c(1, 0),
        legend.background = element_blank(),
        legend.title=element_blank(),
        plot.margin = unit(c(3.5, 5.5, 3.5, 5.5), "pt")
    ) +
    scale_color_manual(values=c("RAD" = "red", "pred_ok" = "blue", "smooth_RAD" = "red", "smooth_ok" = "blue")) +
    scale_linetype_manual(values=c("RAD" = "dashed", "pred_ok" = "dashed", "smooth_RAD" = "solid", "smooth_ok" = "solid")) +
    scale_size_manual(values=c("RAD" = 0.3, "pred_ok" = 0.3, "smooth_RAD" = 0.7, "smooth_ok" = 0.7))


p3 <- ggplot(long_format %>% filter(Model %in% c('RAD','pred_gam_ok', 'smooth_RAD', 'smooth_gam_ok')) %>% 
                 mutate(Model=factor(Model, c('RAD','pred_gam_ok', 'smooth_RAD', 'smooth_gam_ok'))), aes(x=date, y=prediction, color = Model, linetype = Model, size = Model)) +
    geom_line() +
    ylab(TeX("SSR \\[$Wm^{-2}$\\]")) +
    xlab('Date') +
    labs(subtitle="Panel C") +
    guides(col=guide_legend(nrow = 1),
           linetype=guide_legend(nrow = 1)) +
    theme_bw() + # theme_bw() will overwirte all theme() settings
    theme(
        plot.subtitle = element_text(size=12, vjust=-.05),
        axis.title.y = element_text(size = 12), 
        legend.justification=c(1,0),
        legend.position = c(1, 0),
        legend.background = element_blank(),
        legend.title=element_blank(),
        plot.margin = unit(c(0, 5.5, 5.5, 5.5), "pt")
    ) +
    scale_color_manual(values=c("RAD" = "red", "pred_gam_ok" = "blue", "smooth_RAD" = "red", "smooth_gam_ok" = "blue")) +
    scale_linetype_manual(values=c("RAD" = "dashed", "pred_gam_ok" = "dashed", "smooth_RAD" = "solid", "smooth_gam_ok" = "solid")) +
    scale_size_manual(values=c("RAD" = 0.3, "pred_gam_ok" = 0.3, "smooth_RAD" = 0.7, "smooth_gam_ok" = 0.7))


p4 <- ggplot(long_format %>% 
                 filter(Model %in% c('RAD','pred_lsdv_ok', 'smooth_RAD', 'smooth_lsdv_ok')) %>% 
                 mutate(Model=factor(Model, c('RAD','pred_lsdv_ok', 'smooth_RAD', 'smooth_lsdv_ok'))), 
             aes(x=date, y=prediction, color = Model, linetype = Model, size = Model)) +
    geom_line() +
    ylab(TeX("SSR \\[$Wm^{-2}$\\]")) +
    xlab('Date') +
    guides(col=guide_legend(nrow = 1),
           linetype=guide_legend(nrow = 1)) +
    theme_bw() + # theme_bw() will overwirte all theme() settings
    theme(
        plot.subtitle = element_text(size=12, vjust=-.05),
        axis.title.y = element_text(size = 12), 
        legend.justification=c(1,0),
        legend.position = c(1, 0),
        legend.background = element_blank(),
        legend.title=element_blank(),
        plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt")
    ) +
    scale_color_manual(values=c("RAD" = "red", "pred_lsdv_ok" = "blue", "smooth_RAD" = "red", "smooth_lsdv_ok" = "blue")) +
    scale_linetype_manual(values=c("RAD" = "dashed", "pred_lsdv_ok" = "dashed", "smooth_RAD" = "solid", "smooth_lsdv_ok" = "solid")) +
    scale_size_manual(values=c("RAD" = 0.3, "pred_lsdv_ok" = 0.3, "smooth_RAD" = 0.7, "smooth_lsdv_ok" = 0.7))




dev.new()
f_name <- paste0(data_dir, con_folder, "monthly_anomalies_EU_lsdv.png")
f_name
png(f_name, width=10*ppi, height=6*ppi, res=ppi)
p4
dev.off()

dev.new()
f_name <- paste0(data_dir, con_folder, "monthly_anomalies_EU.png")
f_name
png(f_name, width=10*ppi, height=8*ppi, res=ppi)
grid.arrange(p1, p2, p3, ncol=1)
dev.off()


p <- ggplot(long_format %>% filter(Model %in% c('RAD','pred_RF','smooth_RAD','smooth_RF')) %>% 
                mutate(Model=factor(Model, c('RAD','pred_RF','smooth_RAD','smooth_RF'))), aes(x=date, y=prediction, color = Model, linetype = Model, size = Model)) +
    geom_line() +
    ylab(TeX("SSR \\[$Wm^{-2}$\\]")) +
    xlab('Year') +
    guides(col=guide_legend(nrow = 1),
           linetype=guide_legend(nrow = 1)) +
    theme_bw() + # theme_bw() will overwirte all theme() settings
    theme(
        plot.subtitle = element_text(size=12, vjust=-.05),
        axis.title.y = element_text(size = 12), 
        legend.justification=c(1,0),
        legend.position = c(1, 0),
        legend.background = element_blank(),
        legend.title=element_blank(),
        plot.margin = unit(c(5.5, 5.5, 0, 5.5), "pt")
    ) +
    scale_color_manual(values=c("RAD" = "red", "pred_RF" = "blue", "smooth_RAD" = "red", "smooth_RF" = "blue")) +
    scale_linetype_manual(values=c("RAD" = "dashed", "pred_RF" = "dashed", "smooth_RAD" = "solid", "smooth_RF" = "solid")) +
    scale_size_manual(values=c("RAD" = 0.3, "pred_RF" = 0.3, "smooth_RAD" = 0.7, "smooth_RF" = 0.7))

dev.new()
f_name <- paste0(data_dir, con_folder, 
                 sprintf("monthly_anomalies_%s.png", con_name))
f_name
png(f_name, width=10*ppi, height=6*ppi, res=ppi)
p
dev.off()




##    2. Yearly Anomalies    ##  -----------
station$YR <- station$date %>% year()
station$MON <- station$date %>% month()
head(station)
tail(station)
# View(station)

yearly_mean <- station %>% group_by(YR) %>% 
    summarise(RAD=mean(RAD),
              pred_RF=mean(pred_RF),
              pred_ok=mean(pred_ok, na.rm=T),
              pred_gam_ok=mean(pred_gam_ok, na.rm=T),
              pred_lsdv_ok=mean(pred_lsdv_ok, na.rm=T)
    )


yearly_mean <- yearly_mean %>% mutate(
    RAD_ma = rollapply(RAD, 5, mean, align='right', partial=TRUE),
    pred_RF_ma = rollapply(pred_RF, 5, mean, align='right', partial=TRUE)
    )


yearly_trend <- yearly_mean[,-1] - matrix(rep(colMeans(yearly_mean[,-1]) %>% as.numeric(), nrow(yearly_mean)),
       nrow=nrow(yearly_mean), byrow=TRUE)
yearly_trend$YR <- yearly_mean$YR


long_format <- gather(yearly_trend, key="Model", value="prediction", -YR)

p1 <- ggplot(long_format %>% filter(Model %in% c('RAD','pred_RF','RAD_ma','pred_RF_ma')) %>% 
                 mutate(Model=factor(Model, c('RAD','pred_RF','RAD_ma','pred_RF_ma'))), aes(x=YR, y=prediction, color = Model, linetype = Model, size = Model)) +
    geom_line() +
    ylab(TeX("SSR \\[$Wm^{-2}$\\]")) +
    # xlab('') +
    labs(subtitle="Panel A") +
    guides(col=guide_legend(nrow = 1),
           linetype=guide_legend(nrow = 1)) +
    theme_bw() + # theme_bw() will overwirte all theme() settings
    theme(
        axis.title.x=element_blank(),
        plot.subtitle = element_text(size=12, vjust=-.05),
        axis.title.y = element_text(size = 12), 
        legend.justification = c(0, 1),
        legend.position = c(0, 1),
        legend.background = element_blank(),
        legend.title=element_blank(),
        plot.margin = unit(c(5.5, 5.5, 0, 5.5), "pt")
    ) +
    scale_color_manual(values=c("RAD" = "red", "pred_RF" = "blue", "RAD_ma" = "red", "pred_RF_ma" = "blue")) +
    scale_linetype_manual(values=c("RAD" = "dashed", "pred_RF" = "dashed", "RAD_ma" = "solid", "pred_RF_ma" = "solid")) +
    scale_size_manual(values=c("RAD" = 0.3, "pred_RF" = 0.3, "RAD_ma" = 0.7, "pred_RF_ma" = 0.7))

dev.new()
p1




##     R sqaured simulation vs observation    ## --------------------
lm_pred_RF <- lm(pred_RF~RAD, continent_prediction_RF)
lm_pred_RF$coef[2]
summary(lm_pred_RF)$adj.r.squared

lm_pred_ok <- lm(pred_ok~RAD, continent_prediction_RF)
summary(lm_pred_ok)$adj.r.squared

lm_pred_ols <- lm(pred_ols~RAD, continent_prediction_RF)
summary(lm_pred_ols)$r.squared

lm_pred_ols_ok <- lm(pred_ols_ok~RAD, continent_prediction_RF)
summary(lm_pred_ols_ok)$adj.r.squared

lm_pred_gam_ok <- lm(pred_gam_ok~RAD, continent_prediction_RF)
summary(lm_pred_gam_ok)$adj.r.squared

lm_pred_lsdv_ok <- lm(pred_lsdv_ok~RAD, continent_prediction_RF)
summary(lm_pred_lsdv_ok)$adj.r.squared


p1 <- ggplot(continent_prediction_RF, aes(x=RAD, y=pred_RF)) + 
    geom_point() + 
    geom_smooth(method = "lm", se = FALSE, color="red") +
    xlab(TeX("Observed SSR \\[$Wm^{-2}$\\]")) + 
    ylab(TeX("Simulated SSR \\[$Wm^{-2}$\\]")) +  
    ggtitle("Model: RF") + 
    annotate("text", x=0, y=380, 
             label=TeX(sprintf("y=%s*x+%s", 
                       round(lm_pred_RF$coef[2],2), 
                       round(lm_pred_RF$coef[1],2) )),
                       parse=TRUE, hjust = 0) + 
    annotate("text", x=0, y=350, 
             label=TeX(sprintf("$R^2$=%s",
                               round(summary(lm_pred_RF)$adj.r.squared,2)) ),
             parse=TRUE, hjust = 0) + 
    theme_bw() + 
    theme(plot.title = element_text(hjust=0.5)) # center title


p2 <- ggplot(continent_prediction_RF, aes(x=RAD, y=pred_ok)) + 
    geom_point() + 
    geom_smooth(method = "lm", se = FALSE, color="red") +
    xlab(TeX("Observed SSR \\[$Wm^{-2}$\\]")) + 
    ylab(TeX("Simulated SSR \\[$Wm^{-2}$\\]")) +  
    ggtitle("Model: OK") + 
    annotate("text", x=0, y=380, 
             label=TeX(sprintf("y=%s*x+%s", 
                               round(lm_pred_ok$coef[2],2), 
                               round(lm_pred_ok$coef[1],2) )),
             parse=TRUE, hjust = 0) + 
    annotate("text", x=0, y=350, 
             label=TeX(sprintf("$R^2$=%s",
                               round(summary(lm_pred_ok)$adj.r.squared,2)) ),
             parse=TRUE, hjust = 0) + 
    theme_bw() + 
    theme(plot.title = element_text(hjust=0.5)) # center title


p3 <- ggplot(continent_prediction_RF, aes(x=RAD, y=pred_gam_ok)) + 
    geom_point() + 
    geom_smooth(method = "lm", se = FALSE, color="red") +
    xlab(TeX("Observed SSR \\[$Wm^{-2}$\\]")) + 
    ylab(TeX("Simulated SSR \\[$Wm^{-2}$\\]")) +  
    ggtitle("Model: GAM+OK") + 
    annotate("text", x=0, y=380, 
             label=TeX(sprintf("y=%s*x+%s", 
                               round(lm_pred_gam_ok$coef[2],2), 
                               round(lm_pred_gam_ok$coef[1],2) )),
             parse=TRUE, hjust = 0) + 
    annotate("text", x=0, y=350, 
             label=TeX(sprintf("$R^2$=%s",
                               round(summary(lm_pred_gam_ok)$adj.r.squared,2)) ),
             parse=TRUE, hjust = 0) + 
    theme_bw() + 
    theme(plot.title = element_text(hjust=0.5)) # center title


p4 <- ggplot(continent_prediction_RF, aes(x=RAD, y=pred_lsdv_ok)) + 
    geom_point() + 
    geom_smooth(method = "lm", se = FALSE, color="red") +
    xlab(TeX("Observed SSR \\[$Wm^{-2}$\\]")) + 
    ylab(TeX("Simulated SSR \\[$Wm^{-2}$\\]")) +  
    ggtitle("Model: LSDV+OK") + 
    annotate("text", x=0, y=380, 
             label=TeX(sprintf("y=%s*x+%s", 
                               round(lm_pred_lsdv_ok$coef[2],2), 
                               round(lm_pred_lsdv_ok$coef[1],2) )),
             parse=TRUE, hjust = 0) + 
    annotate("text", x=0, y=350, 
             label=TeX(sprintf("$R^2$=%s",
                               round(summary(lm_pred_lsdv_ok)$adj.r.squared,2)) ),
             parse=TRUE, hjust = 0) + 
    theme_bw() + 
    theme(plot.title = element_text(hjust=0.5)) # center title


dev.new()
f_name <- paste0(data_dir, con_folder, 
                 sprintf("simu_obs_%s.png", con_name))
f_name
ppi <- 300
png(f_name, width=10*ppi, height=8*ppi, res=ppi)
grid.arrange(p1, p2, p3, p4, ncol=2, nrow = 2)
dev.off()













