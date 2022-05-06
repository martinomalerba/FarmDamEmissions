##############################################################################################################
#  Analyses for "Fencing farm dams to exclude livestock halves methane emissions and improves water quality" #
##############################################################################################################


rm(list = ls(all = TRUE))

# General libraries
library(PerformanceAnalytics) # Pairplots
library(nlme)
library(ggplot2)
library(cowplot)
library(effects)
library(visreg)
library(car)
library(MuMIn)
library(dplyr)
library(lsmeans)
library(ggsignif)
library(sf)
library(reshape2)
library(readxl)

# For log10 axis annotation
library(scales)

# gglocator
library(ggmap)

# Check non-linearity using gamm
library(mgcv)

# Partial residual plots
library(remef)

# Citations to include in the ms
citation("plyr")
citation("tidyverse")
citation("nlme")
citation("effects")

##############
# CHANGE REQUIRED
# Set where the R code is located
dir.script = ""
setwd(dir.script)
##############

# Load the data
load("../Data/AllData.RData")



# DESCRIPTION OF THE DATA (AllData)
# -------------
# dim(AllData) = 699x32

# Each row is a rate of methane emissions recorded at a farm dam 


# COLUMNS:

# names(AllData) =
# "SiteCode" (unique site ID)
# "Location" (14 locations, replicate readings from 2 to 8 at each site)                      
# "Study_region" (3 locations, rep 7 to 34)
# "Dam_type" (Unfenced Fenced - rep 7 to 32)   
# "AirTemp" (measured onsite)
# "WaterTemp"  (measured onsite)
# "Oxygen_mg"  (measured onsite)
# "Oxygen_perc" (do not use - almost perfect corection with mg)
# "Conductivity" (5 NA)                        

# New nutrient measurements:
# "Total Nitrogen as N"
# "Total Phosphorus as P"

# "Date_UGGA" (10 levels, rep 5 to 8)

#"FID" (another ID for farm dams)
#"Lat" (Lat farm dams)
#"Lon" (Lon farm dams)
#"Area_M2_gmap" (Area measured with gmap)
#"Stock_tCha" (Carbon stock size from CN analysis)
#"Flux_mgm2day" (Flux of GHG emitted)
#"Type" (diffusive emissions)
#"ppm_dif" (Difference in ppm between start and finish)
#"Gas" (Is it CH4 or CO2)



# Create the dataset for the statistical models
Data.lm = data.frame(
			Dam_type = AllData$Dam_type,
			Dam_type2 = AllData$Dam_type2,
			id = AllData$id,
			Owner = AllData$Owner,
			AirTemp = scale(AllData$AirTemp),
			WaterTemp = scale((AllData$WaterTemp)),
			WaterTemp_raw = ((AllData$WaterTemp)),
			Oxygen_mg = scale((AllData$Oxygen_mg)),
			Oxygen_mg_raw = AllData$Oxygen_mg,
			pH_raw = AllData$pHMean,
			Area_M2_gmap = scale(log10(AllData$Area_M2_gmap)),
			Area_M2_gmap_raw = (log10(AllData$Area_M2_gmap)),
			Conductivity = scale(log10(AllData$Conductivity)),
			Conductivity_raw = (log10(AllData$Conductivity)),
			TotalN = scale(log10(AllData$TotalN)),
			TotalP = scale(log10(AllData$TotalP)),
			Stock_tCha = scale(log10(AllData$Stock_tCha)),
			Flux_mgm2day_tr = AllData$Flux_mgm2day,
			Type = AllData$Type,
			SiteCode = AllData$SiteCode,
			Gas = AllData$Gas,
			Lon = AllData$Lon,
			Lat = AllData$Lat,
			TotalN_raw = (log10(AllData$TotalN)),
			TotalP_raw = (log10(AllData$TotalP)),
			Stock_tCha_raw = log10(AllData$Stock_tCha)
	)





# Create four datasets for CO2 and CH4
# CH4 Diffusion
Data.lm_CH4_D = Data.lm %>% subset(Gas == "CH4_ppm" & Type == "D") %>%
				dplyr::mutate(Flux_mgm2day = log10(Flux_mgm2day_tr+2)) %>%
				# Create an ID for the replicates within each site for each Gas
				dplyr::mutate(Rep.group = as.numeric(as.factor(interaction(Gas, SiteCode)))) %>% 
				group_by(Rep.group) %>% 
				mutate(id.rep = 1:length(Rep.group))



# CO2 Diffusion
Data.lm_CO2_D = Data.lm %>% subset(Gas == "CO2_ppm" & Type == "D") %>%
				dplyr::mutate(Flux_mgm2day = log10(Flux_mgm2day_tr+1800)) %>%
				# Create an ID for the replicates within each site for each Gas
				dplyr::mutate(Rep.group = as.numeric(as.factor(interaction(Gas, SiteCode)))) %>% 
				group_by(Rep.group) %>% 
				mutate(id.rep = 1:length(Rep.group))



# Calculate the mean at each site for each dataset
MeanData_CH4_D = Data.lm_CH4_D %>% 
		   dplyr::select(Flux_mgm2day, Oxygen_mg,Oxygen_mg_raw,Conductivity_raw,WaterTemp,WaterTemp_raw,
										 Stock_tCha,Stock_tCha_raw,
		                 TotalN_raw,TotalN,
		                 TotalP_raw,TotalP,
		                 Area_M2_gmap,SiteCode,Dam_type,Dam_type2,Owner,
		                 Lon, Lat, Area_M2_gmap_raw) %>%
		   dplyr::group_by(SiteCode, Dam_type, Dam_type2, Owner) %>%
		   dplyr::summarise(across(everything(), list(mean))) %>%
			 # Backtransform
			 dplyr::mutate(Flux_mgm2day_tr_1 = (10^Flux_mgm2day_1)-2)



MeanData_CO2_D = Data.lm_CO2_D %>% 
		   dplyr::select(Flux_mgm2day, Oxygen_mg,Oxygen_mg_raw,Conductivity_raw,WaterTemp,WaterTemp_raw,
										 Stock_tCha,Stock_tCha_raw,
		                 TotalN_raw,TotalN,
		                 TotalP_raw,TotalP,
		                 Area_M2_gmap,SiteCode,Dam_type,Dam_type2,Owner,
		                 Lon, Lat, Area_M2_gmap_raw) %>%
		   dplyr::group_by(SiteCode, Dam_type, Dam_type2, Owner) %>%
		   dplyr::summarise(across(everything(), list(mean))) %>%
			 # Backtransform
			 dplyr::mutate(Flux_mgm2day_tr_1 = (10^(Flux_mgm2day_1))-1800)


MeanData_CO2eq_D = 
		MeanData_CH4_D %>% 
			 dplyr::rename(Flux_mgm2day_CH4_D = Flux_mgm2day_1) %>% 
			 dplyr::select(Flux_mgm2day_CH4_D, Oxygen_mg_1,Oxygen_mg_raw_1,Conductivity_raw_1,WaterTemp_1,WaterTemp_raw_1,
										 Stock_tCha_1,Stock_tCha_raw_1,
		                 TotalN_raw_1,TotalN_1,
		                 TotalP_raw_1,TotalP_1,
		                 Area_M2_gmap_1,SiteCode, Dam_type, Dam_type2, Owner,
		                 Lon_1, Lat_1, Area_M2_gmap_raw_1) %>%
		merge(
			 MeanData_CO2_D %>% 
			 dplyr::rename(Flux_mgm2day_CO2_D = Flux_mgm2day_1) %>% 
			 dplyr::select(SiteCode, Flux_mgm2day_CO2_D,Dam_type,Dam_type2,Owner),
			 by = c("SiteCode","Dam_type","Dam_type2","Owner")
			 ) %>%

		# Backtransform
		mutate(Flux_mgm2day_CH4_D_tr = (10^(Flux_mgm2day_CH4_D))-2,
					 Flux_mgm2day_CO2_D_tr = (10^(Flux_mgm2day_CO2_D))-1800) %>%

		# Transform into CO2eq
		# using the 20‐year global warming potential of CH4 from Neubauer and Megonigal (2015)
		mutate(Flux_mgm2day_CO2eq_D_tr = Flux_mgm2day_CO2_D_tr + 96*(Flux_mgm2day_CH4_D_tr)) %>%

		# Get ratio of CH4 to CO2
		mutate(Flux_mgm2day_CH4toCO2_D_tr = (96*(Flux_mgm2day_CH4_D_tr))/Flux_mgm2day_CO2_D_tr) %>%

		# Transform to log10
		mutate(Flux_mgm2day_CO2eq_D = log10(Flux_mgm2day_CO2eq_D_tr + 1800))





# ------------------------
# MANAGEMENT ON CH4_D
# ------------------------

# Aim: test the effects of dam management on CH4 emissions
# Dep.Var: Flux_mgm2day
# Exp.Var: Dam_type
# Mod.name: CH4_D
# Data: Data.lm_CH4_D



CH4_D.full1 = lme(Flux_mgm2day ~  Dam_type2,
				random = ~1|Owner/SiteCode,
				weights = varIdent(form = ~1|Dam_type2),
				method = "ML",
				data = Data.lm_CH4_D)

# Not working
#CH4_D.full2 = lme(Flux_mgm2day ~  Dam_type2,
#				random = ~1|Owner/SiteCode,
#				weights = varComb(varExp(),varIdent(form = ~1|Dam_type2)),
#				method = "ML",
#				data = Data.lm_CH4_D)

# Not working
#CH4_D.full3 = lme(Flux_mgm2day ~  Dam_type2,
#				random = ~1|Owner/SiteCode,
#				weights = varComb(varPower()),
#				method = "ML",
#				data = Data.lm_CH4_D)

CH4_D.full4 = lme(Flux_mgm2day ~  Dam_type2,
				random = ~1|Owner/SiteCode,
				method = "ML",
				data = Data.lm_CH4_D)

AICc(CH4_D.full1, CH4_D.full4)

CH4_D.best = CH4_D.full1

summary(CH4_D.best)
Anova(CH4_D.best, type=3)
plot(allEffects(CH4_D.best))
plot(CH4_D.best)
capture.output(Anova(CH4_D.best, type = "III"), file = "../Results/CH4_D.best.txt")
save(CH4_D.best, file = "../Results/CH4_D.best.RData")
lsmeans(object=CH4_D.best, pairwise ~ Dam_type2, adjust= "tukey")



CH4_D.full1b = lme(Flux_mgm2day ~  -1 + Dam_type2,
				random = ~1|Owner/SiteCode,
				weights = varIdent(form = ~1|Dam_type2),
				method = "REML",
				data = Data.lm_CH4_D)


Data.CH4_D.full1b <- as.data.frame(cbind("Dam_type2" = c("Unfenced", "Fenced"), intervals(CH4_D.full1b)$fixed)) # Data to plot
Data.CH4_D.full1b.trs <- as.data.frame(cbind("Dam_type2" = c("Unfenced", "Fenced"), ((10^intervals(CH4_D.full1b)$fixed)-2))) # Data to plot



CH4_D.full1b.plot <- ggplot(data = Data.CH4_D.full1b.trs, aes(x=Dam_type2, y=as.numeric(est.))) +
  geom_bar(stat='identity', width=0.8, fill=c("darkred","forestgreen"), alpha=0.7) +
  geom_errorbar(aes(x=Dam_type2, ymin=as.numeric(lower), ymax=as.numeric(upper)), width=0.1, colour="black", alpha=0.7, size=.5) +
  labs(x = "Farm dam management",
       y = expression(paste("CH"[4]," diffusion (mg m"^-2, " day"^-1, ")"))) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  geom_signif(comparisons=list(c("Unfenced", "Fenced")), annotations="*", textsize=10,
              y_position = 14, tip_length = 0, vjust=0.6) +
  scale_y_continuous(expand = expansion(mult = c(0, .1)))



# WINNER PLOT
CH4_D.full1b.plot_RawData_log10 <- 
	ggplot(data = Data.lm_CH4_D, aes(x = Dam_type2)) +
  geom_jitter(aes(y = 10^Flux_mgm2day), width = 0.2, alpha = 0.5) +
	#ggplot(data = Data.CH4_D.full1b.trs, aes(x=Dam_type2, y=as.numeric(est.))) +
  #geom_bar(stat='identity', width=0.8, fill=c("darkred","forestgreen"), alpha=0.7) +
  geom_errorbar(data = Data.CH4_D.full1b, aes(ymin=10^as.numeric(lower), ymax=10^as.numeric(upper)), width=0.1, colour="red", alpha= 1, size= 2) +
	geom_point(data = Data.CH4_D.full1b, aes(y=10^as.numeric(est.)), colour="red", size= 8) +
	#stat_summary(fun = "mean", colour = "red", size = 2, geom = "point") +
	#stat_summary(fun.data = "mean_cl_boot", colour = "red", size = 2) +
  labs(x = "Farm dam management",
       y = expression(paste("CH"[4]," diffusion (mg m"^-2, " day"^-1, ")"))) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  geom_signif(aes(y = Flux_mgm2day_tr), comparisons=list(c("Unfenced", "Fenced")), annotations="*", textsize=10,
              y_position = 2.5, tip_length = 0, vjust=0.6) +#
  #scale_y_log10(expand = expansion(mult = c(0, .1)), limits = c(0.05, 1000))
  scale_y_log10(expand = expansion(mult = c(0, .1)), limits = c(0.05, 1000),
  						  breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "l")

CH4_D.full1b.plot_RawData_log10b <- 
	ggplot(data = Data.lm_CH4_D, aes(x = Dam_type2)) +
  geom_jitter(aes(y = 10^(Flux_mgm2day)-2), width = 0.2, alpha = 0.2, size = 3) +
	#ggplot(data = Data.CH4_D.full1b.trs, aes(x=Dam_type2, y=as.numeric(est.))) +
  #geom_bar(stat='identity', width=0.8, fill=c("darkred","forestgreen"), alpha=0.7) +
  #geom_errorbar(data = Data.CH4_D.full1b.trs, aes(ymin=as.numeric(lower), ymax=as.numeric(upper)), width=0.1, colour="red", alpha= 1, size= 2) +
	#geom_point(data = Data.CH4_D.full1b.trs, aes(y=as.numeric(est.)), colour="red", size= 8) +
	#stat_summary(fun = "mean", colour = "red", size = 2, geom = "point") +
	#stat_summary(fun.data = "mean_cl_boot", colour = "red", size = 2) +
  geom_pointrange(data = Data.CH4_D.full1b.trs, aes(y=as.numeric(est.), ymin=as.numeric(lower), ymax=as.numeric(upper)), size = 2, fatten = 2) +
	labs(x = "Farm dam management",
       y = expression(paste("CH"[4]," diffusion (mg m"^-2, " day"^-1, " + 2)"))) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  geom_signif(aes(y = Flux_mgm2day_tr), comparisons=list(c("Unfenced", "Fenced")), annotations="*", textsize=10,
              y_position = 2.5, tip_length = 0, vjust=0.6) +#
  #scale_y_log10(expand = expansion(mult = c(0, .1)), limits = c(0.05, 1000))
  scale_y_log10(expand = expansion(mult = c(0, .1)), limits = c(0.05, 1000),
  						  breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "l")


CH4_D.full1b.plot_RawData <- 
	ggplot(data = Data.lm_CH4_D, aes(x = Dam_type2)) +
  geom_jitter(aes(y = Flux_mgm2day_tr), width = 0.1, alpha = 0.5) +
	#ggplot(data = Data.CH4_D.full1b.trs, aes(x=Dam_type2, y=as.numeric(est.))) +
  #geom_bar(stat='identity', width=0.8, fill=c("darkred","forestgreen"), alpha=0.7) +
  geom_errorbar(data = Data.CH4_D.full1b.trs, aes(ymin=as.numeric(lower), ymax=as.numeric(upper)), width=0.1, colour="red", alpha= 1, size= 2) +
	geom_point(data = Data.CH4_D.full1b.trs, aes(y=as.numeric(est.)), colour="red", size= 4) +
	#stat_summary(fun = "mean", colour = "red", size = 2, geom = "point") +
	#stat_summary(fun.data = "mean_cl_boot", colour = "red", size = 2) +
  labs(x = "Farm dam management",
       y = expression(paste("CH"[4]," diffusion (mg m"^-2, " day"^-1, ")"))) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  geom_signif(aes(y = Flux_mgm2day_tr), comparisons=list(c("Unfenced", "Fenced")), annotations="*", textsize=10,
              y_position = 250, tip_length = 0, vjust=0.6)


AllCH4Plots = plot_grid(CH4_D.full1b.plot,CH4_D.full1b.plot_RawData_log10,CH4_D.full1b.plot_RawData, nrow = 1)
save_plot(AllCH4Plots, base_height = 7.33, base_width = 11.05, file = "../Results/AllCH4Plots.pdf")

AllCH4Plots2 = plot_grid(CH4_D.full1b.plot,CH4_D.full1b.plot_RawData_log10, nrow = 1)
save_plot(AllCH4Plots2, base_height = 7.33, base_width = 11.05, file = "../Results/AllCH4Plots2.pdf")

 
CH4_D.full1b.plot
ggsave("../Results/CH4_D.best.pdf", width = 5.24, height = 4.06)



# Model validation
# ================
if(FALSE){
	Data.lm_CH4_D$resid_CH4_D = resid(CH4_D.best, type = "pearson")
	Data.lm_CH4_D$pred_CH4_D = predict(CH4_D.best)

	scatterplot(data = Data.lm_CH4_D, resid_CH4_D~pred_CH4_D)
	boxplot(data = Data.lm_CH4_D, resid_CH4_D~Dam_type)

	scatterplot(data = Data.lm_CH4_D, Flux_mgm2day~pred_CH4_D)
	scatterplot(data = Data.lm_CH4_D, resid_CH4_D~WaterTemp)
	scatterplot(data = Data.lm_CH4_D, resid_CH4_D~Oxygen_mg)
	scatterplot(data = Data.lm_CH4_D, resid_CH4_D~Conductivity)
	scatterplot(data = Data.lm_CH4_D, resid_CH4_D~Area_M2_gmap)
	scatterplot(data = Data.lm_CH4_D, resid_CH4_D~Nitrate)
	scatterplot(data = Data.lm_CH4_D, resid_CH4_D~Phosphate)
	scatterplot(data = Data.lm_CH4_D, resid_CH4_D~NDVI.veg)

	boxplot(data = Data.lm_CH4_D, Flux_mgm2day~Dam_type)
	boxplot(data = Data.lm_CH4_D, resid_CH4_D~Owner)
	boxplot(data = Data.lm_CH4_D, resid_CH4_D~SiteCode)
}









# ------------------------
# MANAGEMENT ON CO2_D
# ------------------------

# Aim: test the effects of dam management on CO2 emissions
# Dep.Var: Flux_mgm2day
# Exp.Var: Dam_type
# Mod.name: CO2_D
# Data: Data.lm_CO2_D

CO2_D.full1 = lme(Flux_mgm2day ~  Dam_type2,
				random = ~1|Owner/SiteCode,
				weights = varIdent(form = ~1|Dam_type2),
				method = "ML",
				data = Data.lm_CO2_D)


CO2_D.full2 = lme(Flux_mgm2day ~  Dam_type2,
				random = ~1|Owner/SiteCode,
				weights = varComb(varPower(),varIdent(form = ~1|Dam_type2)),
				method = "REML",
				data = Data.lm_CO2_D)

CO2_D.full3 = lme(Flux_mgm2day ~  Dam_type2,
				random = ~1|Owner/SiteCode,
				weights = varPower(),
				method = "ML",
				data = Data.lm_CO2_D)

CO2_D.full4 = lme(Flux_mgm2day ~  Dam_type2,
				random = ~1|Owner/SiteCode,
				method = "ML",
				data = Data.lm_CO2_D)

AICc(CO2_D.full1, CO2_D.full2, CO2_D.full3, CO2_D.full4)

CO2_D.best = CO2_D.full3

summary(CO2_D.best)
Anova(CO2_D.best, type=3) # Chi2 = 223.49, P < 2.2e-16 | Martino: Chi2 = same
plot(allEffects(CO2_D.best))
plot(CO2_D.best)
capture.output(Anova(CO2_D.best, type = "III"), file = "../Results/CO2_D.best.txt")
save(CO2_D.best, file = "../Results/CO2_D.best.RData")
lsmeans(object=CO2_D.best, pairwise ~ Dam_type2, adjust= "tukey")



CO2_D.full3b = lme(Flux_mgm2day ~  -1 + Dam_type2,
				random = ~1|Owner/SiteCode,
				weights = varComb(varPower()),
				method = "ML",
				data = Data.lm_CO2_D)

Data.CO2_D.full3b <- as.data.frame(cbind("Dam_type2" = c("Unfenced", "Fenced"), intervals(CO2_D.full3b)$fixed)) # Data to plot
Data.CO2_D.full3b.trs <- as.data.frame(cbind("Dam_type2" = c("Unfenced", "Fenced"), ((10^intervals(CO2_D.full3b)$fixed)-1800))) # Data to plot


CO2_D.full3b.plot <- ggplot(data = Data.CO2_D.full3b.trs, aes(x=Dam_type2, y=as.numeric(est.))) +
  geom_bar(stat='identity', width=0.8, fill=c("darkred","forestgreen"), alpha=0.7) +
  geom_errorbar(aes(x=Dam_type2, ymin=as.numeric(lower), ymax=as.numeric(upper)), width=0.1, colour="black", alpha=0.7, size=.5) +
  labs(x = "Farm dam management",
       y = expression(paste("CO"[2]," diffusion (mg m"^-2, " day"^-1, ")"))) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))# +
#  geom_signif(comparisons=list(c("Unfenced", "Fenced")), annotations="*", textsize=10,
#              y_position = 4000, tip_length = 0, vjust=0.6)
CO2_D.full3b.plot
ggsave("../Results/CO2_D.full3b.pdf", width = 5.24, height = 4.06)








# WINNER PLOT
CO2_D.full1b.plot_RawData_log10 <- 
	ggplot(data = Data.lm_CO2_D, aes(x = Dam_type2)) +
  geom_jitter(aes(y = Flux_mgm2day_tr), width = 0.2, alpha = 0.5) +
	#ggplot(data = Data.CO2_D.full3b.trs, aes(x=Dam_type2, y=as.numeric(est.))) +
  #geom_bar(stat='identity', width=0.8, fill=c("darkred","forestgreen"), alpha=0.7) +
  geom_errorbar(data = Data.CO2_D.full3b.trs, aes(ymin=as.numeric(lower), ymax=as.numeric(upper)), width=0.1, colour="red", alpha= 1, size= 2) +
	geom_point(data = Data.CO2_D.full3b.trs, aes(y=as.numeric(est.)), colour="red", size= 8) +
	#stat_summary(fun = "mean", colour = "red", size = 2, geom = "point") +
	#stat_summary(fun.data = "mean_cl_boot", colour = "red", size = 2) +
  labs(x = "Farm dam management",
       y = expression(paste("CO"[2]," diffusion (mg m"^-2, " day"^-1, ")"))) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  #geom_signif(aes(y = Flux_mgm2day_tr), comparisons=list(c("Unfenced", "Fenced")), annotations="*", textsize=10,
  #            y_position = 2.5, tip_length = 0, vjust=0.6) #+#
  scale_y_log10()#expand = expansion(mult = c(0, .1)), limits = c(0.05, 1000)



CO2_D.full1b.plot_RawData_log10b <- 
	ggplot(data = Data.lm_CO2_D, aes(x = Dam_type2)) +
  geom_jitter(aes(y = Flux_mgm2day_tr), width = 0.2, alpha = 0.2, size = 3) +
	#ggplot(data = Data.CO2_D.full3b.trs, aes(x=Dam_type2, y=as.numeric(est.))) +
  #geom_bar(stat='identity', width=0.8, fill=c("darkred","forestgreen"), alpha=0.7) +
  #geom_errorbar(data = Data.CO2_D.full3b.trs, aes(ymin=as.numeric(lower), ymax=as.numeric(upper)), width=0.1, colour="red", alpha= 1, size= 2) +
	#geom_point(data = Data.CO2_D.full3b.trs, aes(y=as.numeric(est.)), colour="red", size= 8) +
	#stat_summary(fun = "mean", colour = "red", size = 2, geom = "point") +
	#stat_summary(fun.data = "mean_cl_boot", colour = "red", size = 2) +
	geom_pointrange(data = Data.CO2_D.full3b.trs, aes(y=as.numeric(est.), ymin=as.numeric(lower), ymax=as.numeric(upper)), size = 2, fatten = 2) +
  labs(x = "Farm dam management",
       y = expression(paste("CO"[2]," diffusion (mg m"^-2, " day"^-1, ")"))) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_y_log10()  
  #geom_signif(aes(y = Flux_mgm2day_tr), comparisons=list(c("Unfenced", "Fenced")), annotations="*", textsize=10,
  #            y_position = 2.5, tip_length = 0, vjust=0.6) #+#
  #scale_y_log10()#expand = expansion(mult = c(0, .1)), limits = c(0.05, 1000)





CO2_D.full1b.plot_RawData_log10c <- 
	ggplot(data = Data.lm_CO2_D, aes(x = Dam_type2)) +
  geom_jitter(aes(y = (10^Flux_mgm2day)/1000), width = 0.2, alpha = 0.2, size = 3) +
	#ggplot(data = Data.CO2_D.full3b.trs, aes(x=Dam_type2, y=as.numeric(est.))) +
  #geom_bar(stat='identity', width=0.8, fill=c("darkred","forestgreen"), alpha=0.7) +
  #geom_errorbar(data = Data.CO2_D.full3b.trs, aes(ymin=as.numeric(lower), ymax=as.numeric(upper)), width=0.1, colour="red", alpha= 1, size= 2) +
	#geom_point(data = Data.CO2_D.full3b.trs, aes(y=as.numeric(est.)), colour="red", size= 8) +
	#stat_summary(fun = "mean", colour = "red", size = 2, geom = "point") +
	#stat_summary(fun.data = "mean_cl_boot", colour = "red", size = 2) +
	geom_pointrange(data = Data.CO2_D.full3b, aes(y=(10^as.numeric(est.))/1000, ymin=(10^as.numeric(lower))/1000, ymax=(10^as.numeric(upper))/1000), size = 2, fatten = 2) +
  labs(x = "Farm dam management",
       y = expression(paste("CO"[2]," diffusion (g m"^-2, " day"^-1, " + 1.8)"))) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
	scale_y_log10(
	        breaks = scales::trans_breaks("log10", function(x) 10^x),
	        labels = scales::trans_format("log10", scales::math_format(10^.x))
	      ) +
  annotation_logticks(sides = "l")







# Model validation
# ================
if(FALSE){
	Data.lm_CO2_D$resid_CO2_D = resid(CO2_D.best, type = "pearson")
	Data.lm_CO2_D$pred_CO2_D = predict(CO2_D.best)

	scatterplot(data = Data.lm_CO2_D, resid_CO2_D~pred_CO2_D)
	boxplot(data = Data.lm_CO2_D, resid_CO2_D~Dam_type)

	scatterplot(data = Data.lm_CO2_D, Flux_mgm2day~pred_CO2_D)
	scatterplot(data = Data.lm_CO2_D, resid_CO2_D~WaterTemp)
	scatterplot(data = Data.lm_CO2_D, resid_CO2_D~Oxygen_mg)
	scatterplot(data = Data.lm_CO2_D, resid_CO2_D~Conductivity)
	scatterplot(data = Data.lm_CO2_D, resid_CO2_D~Area_M2_gmap)
	scatterplot(data = Data.lm_CO2_D, resid_CO2_D~Nitrate)
	scatterplot(data = Data.lm_CO2_D, resid_CO2_D~Phosphate)
	scatterplot(data = Data.lm_CO2_D, resid_CO2_D~NDVI.veg)

	boxplot(data = Data.lm_CO2_D, Flux_mgm2day~Dam_type)
	boxplot(data = Data.lm_CO2_D, resid_CO2_D~Owner)
	boxplot(data = Data.lm_CO2_D, resid_CO2_D~SiteCode)
}










# ------------------------
# MANAGEMENT ON CO2eq_D (based on averages)
# ------------------------


# Aim: test the effects of dam management on CO2+CH4 emissions
# Dep.Var: Flux_mgm2day_CO2eq_D
# Exp.Var: Dam_type2
# Mod.name: CO2eq_D
# Data: MeanData_CO2eq_D

MeanData_CO2eq_D2 = MeanData_CO2eq_D


CO2eq_D.full1 = lme(Flux_mgm2day_CO2eq_D ~  Dam_type2,
				random = ~1|Owner,
				method = "ML",
				data = MeanData_CO2eq_D2)

CO2eq_D.full2 = lme(Flux_mgm2day_CO2eq_D ~  Dam_type2,
				random = ~1|Owner,
				weights = varIdent(form = ~1|Dam_type2),
				method = "ML",
				data = MeanData_CO2eq_D2)

CO2eq_D.full3 = lme(Flux_mgm2day_CO2eq_D ~  Dam_type2,
				random = ~1|Owner,
				weights = varPower(),
				method = "ML",
				data = MeanData_CO2eq_D2)

CO2eq_D.full4 = lme(Flux_mgm2day_CO2eq_D ~  Dam_type2,
				random = ~1|Owner,
				weights = varComb(varPower(),varIdent(form = ~1|Dam_type2)),
				method = "ML",
				data = MeanData_CO2eq_D2)

CO2eq_D.full5 = lme(Flux_mgm2day_CO2eq_D ~  Dam_type2,
				random = ~1|Owner,
				weights = varIdent(form = ~1|Owner),
				method = "ML",
				data = MeanData_CO2eq_D2)

AICc(CO2eq_D.full1,CO2eq_D.full2,CO2eq_D.full3,CO2eq_D.full4,CO2eq_D.full5)


CO2eq_D.best = CO2eq_D.full1

plot(CO2eq_D.best)
Anova(CO2eq_D.best, type=3)
plot(allEffects(CO2eq_D.best))

capture.output(Anova(CO2eq_D.best, type = "III"), file = "../Results/CO2eq_D.best.txt")
save(CO2eq_D.best, file = "../Results/CO2eq_D.best.RData")



CO2eq_D.best_b = lme(Flux_mgm2day_CO2eq_D ~ -1 + Dam_type2,
				random = ~1|Owner,
				method = "ML",
				data = MeanData_CO2eq_D2)

Anova(CO2eq_D.best_b, type=3)


Data.CO2eq_D.best <- as.data.frame(cbind("Dam_type2" = c("Unfenced", "Fenced"), intervals(CO2eq_D.best_b, which = "fixed")$fixed)) # Data to plot
Data.CO2eq_D.best_tr <- as.data.frame(cbind("Dam_type2" = c("Unfenced", "Fenced"), (10^intervals(CO2eq_D.best_b, which = "fixed")$fixed)-1800)) # Data to plot


# Data.CO2eq_D.full3
lsmeans(object=CO2eq_D.best_b, pairwise ~ Dam_type2, adjust= "tukey")

CO2eq_D.best.plot <- ggplot(data = Data.CO2eq_D.best_tr, aes(x=Dam_type2, y=as.numeric(est.))) +
  geom_bar(stat='identity', width=0.8, fill=c("darkred", "forestgreen"), alpha=0.7) +
  geom_errorbar(aes(x=Dam_type2, ymin=as.numeric(lower), ymax=as.numeric(upper)), width=0.1, colour="black", alpha=0.7, size=.5) +
  labs(x = "Farm dam management",
       y = expression(paste("CO"[2],"-eq. diffusion (mg m"^-2, " day"^-1, ")"))) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  geom_signif(comparisons=list(c("Unfenced", "Fenced")), annotations="p = 0.07", textsize=3,
              y_position = 4000, tip_length = 0, vjust=0.1)
CO2eq_D.best.plot
ggsave("../Results/CO2eq_D.best.plot.pdf", width = 5.24, height = 4.06)






# WINNER PLOT
CO2eq_D.full1b.plot_RawData_log10 <- 
	ggplot(data = MeanData_CO2eq_D2, aes(x = Dam_type2)) +
  geom_jitter(aes(y = 10^(Flux_mgm2day_CO2eq_D)), width = 0.2, alpha = 0.5) +
	#ggplot(data = Data.CO2eq_D.best_tr, aes(x=Dam_type2, y=as.numeric(est.))) +
  #geom_bar(stat='identity', width=0.8, fill=c("darkred","forestgreen"), alpha=0.7) +
  geom_errorbar(data = Data.CO2eq_D.best, aes(ymin=10^as.numeric(lower), ymax=10^as.numeric(upper)), width=0.1, colour="red", alpha= 1, size= 2) +
	geom_point(data = Data.CO2eq_D.best, aes(y=10^as.numeric(est.)), colour="red", size= 8) +
	#stat_summary(fun = "mean", colour = "red", size = 2, geom = "point") +
	#stat_summary(fun.data = "mean_cl_boot", colour = "red", size = 2) +
  labs(x = "Farm dam management",
       y = expression(paste("CO"[2],"-eq. diffusion (mg m"^-2, " day"^-1, ")"))) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_y_log10()



CO2eq_D.full1b.plot_RawData_log10b <- 
	ggplot(data = MeanData_CO2eq_D2, aes(x = Dam_type2)) +
  geom_jitter(aes(y = 10^(Flux_mgm2day_CO2eq_D)/1000), width = 0.2, alpha = 0.2, size = 3) +
	#ggplot(data = Data.CO2eq_D.best, aes(x=Dam_type2, y=as.numeric(est.))) +
  #geom_bar(stat='identity', width=0.8, fill=c("darkred","forestgreen"), alpha=0.7) +
  #geom_errorbar(data = Data.CO2eq_D.best, aes(ymin=as.numeric(lower), ymax=as.numeric(upper)), width=0.1, colour="red", alpha= 1, size= 2) +
	#geom_point(data = Data.CO2eq_D.best, aes(y=as.numeric(est.)), colour="red", size= 8) +
	geom_pointrange(data = Data.CO2eq_D.best, aes(y=10^as.numeric(est.)/1000, ymin=10^as.numeric(lower)/1000, ymax=10^as.numeric(upper)/1000), size = 2, fatten = 2) +
	#stat_summary(fun = "mean", colour = "red", size = 2, geom = "point") +
	#stat_summary(fun.data = "mean_cl_boot", colour = "red", size = 2) +
  labs(x = "Farm dam management",
       y = expression(paste("CO"[2],"-eq. diffusion (g m"^-2, " day"^-1, " + 1.8)"))) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_y_log10(#expand = expansion(mult = c(0, .1)), limits = c(0.05, 1000),
  						  breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "l")




# ----------------
# Model validation
# ----------------

if(FALSE){

MeanData_CO2eq_D2$resid_CO2eq_D = resid(CO2eq_D.best, type = "pearson")
MeanData_CO2eq_D2$pred_CO2eq_D = predict(CO2eq_D.best)


boxplot(data = MeanData_CO2eq_D2, resid_CO2eq_D~Dam_type2)
boxplot(data = MeanData_CO2eq_D2, resid_CO2eq_D~Owner)
boxplot(data = MeanData_CO2eq_D2, Flux_mgm2day_CO2eq_D~Dam_type2)
}




















# ------------------------
# MANAGEMENT ON OXYGEN
# ------------------------
# Aim: test the effects of dam management on CH4 emissions
# Dep.Var: Oxygen_mg_1
# Exp.Var: Dam_type
# Mod.name: CH4_Oxy
# Data: MeanData_CH4_D


CH4_Oxy.full1 = lme(Oxygen_mg_1 ~  -1 + Dam_type2,
				random = ~1|Owner,
				method = "ML",
				data = MeanData_CH4_D)


CH4_Oxy.full2 = lme(Oxygen_mg_1 ~ Dam_type2,
				random = ~1|Owner,
				weights = varIdent(form = ~1|Dam_type2),
				method = "ML",
				data = MeanData_CH4_D)

AICc(CH4_Oxy.full1, CH4_Oxy.full2)

CH4_Oxy.best = CH4_Oxy.full2

plot(CH4_Oxy.best)
Anova(CH4_Oxy.best, type=3)
plot(allEffects(CH4_Oxy.best))

capture.output(Anova(CH4_Oxy.best, type = "III"), file = "../Results/CH4_Oxy.best.txt")
save(CH4_Oxy.best, file = "../Results/CH4_Oxy.best.RData")



CH4_Oxy.best_plot = lme(Oxygen_mg_raw_1 ~ -1 + Dam_type2,
				random = ~1|Owner,
				weights = varIdent(form = ~1|Dam_type2),
				method = "ML",
				data = MeanData_CH4_D)

Anova(CH4_Oxy.best_plot, type=3)

#mean(MeanData_CH4_D$Oxygen_mg_raw_1[MeanData_CH4_D$Dam_type2 == "Fenced"])
#mean(MeanData_CH4_D$Oxygen_mg_raw_1[MeanData_CH4_D$Dam_type2 == "Unfenced"])

intervals(CH4_Oxy.best_plot)$fixed

Data.CH4_Oxy.best <- as.data.frame(cbind("Dam_type2" = c("Unfenced", "Fenced"), intervals(CH4_Oxy.best_plot)$fixed)) 
Data.CH4_Oxy.best_tr <- as.data.frame(cbind("Dam_type2" = c("Unfenced", "Fenced"), intervals(CH4_Oxy.best_plot)$fixed)) 



# Data.CH4_Oxy.best
lsmeans(object=CH4_Oxy.best, pairwise ~ Dam_type2, adjust= "tukey")

CH4_Oxy.best.plot <- ggplot(data = Data.CH4_Oxy.best_tr, aes(x=Dam_type2, y=as.numeric(est.))) +
  geom_bar(stat='identity', width=0.8, fill=c("darkred", "forestgreen"), alpha=0.7) +
  geom_errorbar(aes(x=Dam_type2, ymin=as.numeric(lower), ymax=as.numeric(upper)), width=0.1, colour="black", alpha=0.7, size=.5) +
  labs(x = "Farm dam management",
       y = expression(paste("Oxygen (mg L"^-1, ")"))) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  geom_signif(comparisons=list(c("Unfenced", "Fenced")), annotations="*", textsize=10,
              y_position = 9, tip_length = 0, vjust=0.6)
CH4_Oxy.best.plot
ggsave("../Results/CH4_Oxy.best.plot.pdf", width = 5.24, height = 4.06)








# WINNER PLOT
CH4_Oxy.plot_RawData_log10 <- 
	ggplot(data = MeanData_CH4_D, aes(x = Dam_type2)) +
  geom_jitter(aes(y = Oxygen_mg_raw_1), width = 0.2, alpha = 0.5) +
	#ggplot(data = Data.CH4_Oxy.best_tr, aes(x=Dam_type2, y=as.numeric(est.))) +
  #geom_bar(stat='identity', width=0.8, fill=c("darkred","forestgreen"), alpha=0.7) +
  geom_errorbar(data = Data.CH4_Oxy.best_tr, aes(ymin=as.numeric(lower), ymax=as.numeric(upper)), width=0.1, colour="red", alpha= 1, size= 2) +
	geom_point(data = Data.CH4_Oxy.best_tr, aes(y=as.numeric(est.)), colour="red", size= 8) +
	#stat_summary(fun = "mean", colour = "red", size = 2, geom = "point") +
	#stat_summary(fun.data = "mean_cl_boot", colour = "red", size = 2) +
  labs(x = "Farm dam management",
       y = expression(paste("Oxygen (mg L"^-1, ")"))) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  geom_signif(aes(y = Oxygen_mg_raw_1), comparisons=list(c("Unfenced", "Fenced")), annotations="*", textsize=10,
              y_position = 18, tip_length = 0, vjust=0.6) #+
  #scale_y_continuous(expand = expansion(mult = c(0, .1)), limits = c(-3000, 18000))


CH4_Oxy.plot_RawData_log10b <- 
	ggplot(data = MeanData_CH4_D, aes(x = Dam_type2)) +
  geom_jitter(aes(y = Oxygen_mg_raw_1), width = 0.2, alpha = 0.2, size = 3) +
	#ggplot(data = Data.CH4_Oxy.best_tr, aes(x=Dam_type2, y=as.numeric(est.))) +
  #geom_bar(stat='identity', width=0.8, fill=c("darkred","forestgreen"), alpha=0.7) +
  #geom_errorbar(data = Data.CH4_Oxy.best_tr, aes(ymin=as.numeric(lower), ymax=as.numeric(upper)), width=0.1, colour="red", alpha= 1, size= 2) +
	#geom_point(data = Data.CH4_Oxy.best_tr, aes(y=as.numeric(est.)), colour="red", size= 8) +
	geom_pointrange(data = Data.CH4_Oxy.best_tr, aes(y=as.numeric(est.), ymin=as.numeric(lower), ymax=as.numeric(upper)), size = 2, fatten = 2) +
	#stat_summary(fun = "mean", colour = "red", size = 2, geom = "point") +
	#stat_summary(fun.data = "mean_cl_boot", colour = "red", size = 2) +
  labs(x = "Farm dam management",
       y = expression(paste("Oxygen (mg L"^-1, ")"))) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  geom_signif(aes(y = Oxygen_mg_raw_1), comparisons=list(c("Unfenced", "Fenced")), annotations="*", textsize=10,
              y_position = 18, tip_length = 0, vjust=0.6) #+
  #scale_y_continuous(expand = expansion(mult = c(0, .1)), limits = c(-3000, 18000))





# ----------------
# Model validation
# ----------------

if(FALSE){
	scatterplot(residuals(CH4_Oxy.best, type = "pearson")~predict(CH4_Oxy.best))
	scatterplot(residuals(CH4_Oxy.best, type = "pearson")~MeanData_CH4_D$Oxygen_mg_1)
	boxplot(MeanData_CH4_D$Oxygen_mg_1~MeanData_CH4_D$Dam_type)
}



















# ------------------------
# MANAGEMENT ON NITROGEN
# ------------------------
# Aim: test the effects of dam management on N and P
# Dep.Var: TotalN_1
# Exp.Var: Dam_type
# Mod.name: CH4_N
# Data: MeanData_CH4_D


CH4_Ntot.full1 = lme(TotalN_1 ~  Dam_type2,
				random = ~1|Owner,
				method = "ML",
				data = MeanData_CH4_D)

CH4_Ntot.full2 = lme(TotalN_1 ~  Dam_type2,
				weights = varIdent(form = ~1|Dam_type2),
				random = ~1|Owner,
				method = "ML",
				data = MeanData_CH4_D)
AICc(CH4_Ntot.full1,CH4_Ntot.full2)

CH4_Ntot.best = CH4_Ntot.full1

plot(CH4_Ntot.best)
Anova(CH4_Ntot.best, type=3)
summary(CH4_Ntot.best)
plot(allEffects(CH4_Ntot.best))

capture.output(Anova(CH4_Ntot.best, type = "III"), file = "../Results/CH4_Ntot.best.txt")
save(CH4_Ntot.best, file = "../Results/CH4_Ntot.best.RData")


CH4_Ntot.best_plot = lme(TotalN_raw_1 ~ -1 + Dam_type2,
				random = ~1|Owner,
				#method = "ML",
				data = MeanData_CH4_D)

AICc(CH4_Ntot.best_plot)


CH4_Ntot.best.df <- as.data.frame(cbind("Dam_type2" = c("Unfenced", "Fenced"), 10^(intervals(CH4_Ntot.best_plot)$fixed))) # Data to plot



# Sort unfenced and fenced
CH4_Ntot.best.df$Dam_type <- factor(CH4_Ntot.best.df$Dam_type, levels = c("Unfenced", "Transition", "Fenced"))

# Data.CH4_Ntot.full2
lsmeans(object=CH4_Ntot.best, pairwise ~ Dam_type2, adjust= "tukey") # sig *

CH4_Ntot.best.plot <- ggplot(data = CH4_Ntot.best.df, aes(x=Dam_type2, y=as.numeric(est.))) +
  geom_bar(stat='identity', width=0.8, fill=c("darkred", "forestgreen"), alpha=0.7) +
  geom_errorbar(aes(x=Dam_type2, ymin=as.numeric(lower), ymax=as.numeric(upper)), width=0.1, colour="black", alpha=0.7, size=.5) +
  labs(x = "Farm dam management",
       y = expression(paste("Nitrogen (mg L"^-1, ")"))) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  geom_signif(comparisons=list(c("Unfenced", "Fenced")), annotations="*", textsize=10,
              y_position = 4, tip_length = 0, vjust=0.6)
CH4_Ntot.best.plot
ggsave("../Results/CH4_Ntot.best.plot.pdf", width = 5.24, height = 4.06)





# WINNER PLOT
CH4_Ntot.plot_RawData_log10 <- 
	ggplot(data = MeanData_CH4_D, aes(x = Dam_type2)) +
  geom_jitter(aes(y = 10^TotalN_raw_1), width = 0.2, alpha = 0.5) +
	#ggplot(data = CH4_Ntot.best.df, aes(x=Dam_type2, y=as.numeric(est.))) +
  #geom_bar(stat='identity', width=0.8, fill=c("darkred","forestgreen"), alpha=0.7) +
  geom_errorbar(data = CH4_Ntot.best.df, aes(ymin=as.numeric(lower), ymax=as.numeric(upper)), width=0.1, colour="red", alpha= 1, size= 2) +
	geom_point(data = CH4_Ntot.best.df, aes(y=as.numeric(est.)), colour="red", size= 8) +
	#stat_summary(fun = "mean", colour = "red", size = 2, geom = "point") +
	#stat_summary(fun.data = "mean_cl_boot", colour = "red", size = 2) +
  labs(x = "Farm dam management",
       y = expression(paste("Nitrogen (mg L"^-1, ")"))) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  geom_signif(aes(y = 10^TotalN_raw_1), comparisons=list(c("Unfenced", "Fenced")), annotations="*", textsize=10,
              y_position = 1, tip_length = 0, vjust=0.6) +
  scale_y_log10(expand = expansion(mult = c(0, .1)),
  						  breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
                #limits = c(0.005, 2)) +
  annotation_logticks(sides = "l")


CH4_Ntot.plot_RawData_log10b <- 
	ggplot(data = MeanData_CH4_D, aes(x = Dam_type2)) +
  geom_jitter(aes(y = 10^TotalN_raw_1), width = 0.2, alpha = 0.2, size = 3) +
	#ggplot(data = CH4_Ntot.best.df, aes(x=Dam_type2, y=as.numeric(est.))) +
  #geom_bar(stat='identity', width=0.8, fill=c("darkred","forestgreen"), alpha=0.7) +
  #geom_errorbar(data = CH4_Ntot.best.df, aes(ymin=as.numeric(lower), ymax=as.numeric(upper)), width=0.1, colour="red", alpha= 1, size= 2) +
	#geom_point(data = CH4_Ntot.best.df, aes(y=as.numeric(est.)), colour="red", size= 8) +
	geom_pointrange(data = CH4_Ntot.best.df, aes(y=as.numeric(est.), ymin=as.numeric(lower), ymax=as.numeric(upper)), size = 2, fatten = 2) +
	#stat_summary(fun = "mean", colour = "red", size = 2, geom = "point") +
	#stat_summary(fun.data = "mean_cl_boot", colour = "red", size = 2) +
  labs(x = "Farm dam management",
       y = expression(paste("Nitrogen (mg L"^-1, ")"))) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  geom_signif(aes(y = 10^TotalN_raw_1), comparisons=list(c("Unfenced", "Fenced")), annotations="*", textsize=10,
              y_position = 1, tip_length = 0, vjust=0.6) +
  scale_y_log10(expand = expansion(mult = c(0, .1)),
  						  breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
                #limits = c(0.005, 2)) +
  annotation_logticks(sides = "l")




# ----------------
# Model validation
# ----------------

if(FALSE){
scatterplot(residuals(CH4_Ntot.best, type = "pearson")~predict(CH4_Ntot.best))
scatterplot(residuals(CH4_Ntot.best, type = "pearson")~MeanData_CH4_D$TotalN_1)
boxplot(MeanData_CH4_D$TotalN_1~MeanData_CH4_D$Dam_type)
}










# ------------------------
# MANAGEMENT ON PHOSPHORUS
# ------------------------
# Aim: test the effects of dam management on TotP
# Dep.Var: TotalP_1
# Exp.Var: Dam_type
# Mod.name: CH4_Ptot
# Data: MeanData_CH4_D

# Remove odd value
MeanData_CH4_D3 = subset(MeanData_CH4_D, !(Dam_type == "Fenced" & TotalP_1 > 2))

CH4_Ptot.full1 = lme(TotalP_1 ~  Dam_type2,
				random = ~1|Owner,
				#method = "ML",
				data = MeanData_CH4_D3)

CH4_Ptot.full2 = lme(TotalP_1 ~  Dam_type2,
				weights = varIdent(form = ~1|Dam_type2),
				random = ~1|Owner,
				#method = "ML",
				data = MeanData_CH4_D3)

AICc(CH4_Ptot.full1,CH4_Ptot.full2)

CH4_Ptot.best = CH4_Ptot.full1

plot(CH4_Ptot.best)
Anova(CH4_Ptot.best, type=3)
summary(CH4_Ptot.best)
plot(allEffects(CH4_Ptot.best))

capture.output(Anova(CH4_Ptot.best, type = "III"), file = "../Results/CH4_Ptot.best.txt")
save(CH4_Ptot.best, file = "../Results/CH4_Ptot.best.RData")



CH4_Ptot.best_plot = lme(TotalP_raw_1 ~ -1 + Dam_type2,
				random = ~1|Owner,
				#method = "ML",
				data = MeanData_CH4_D3)

CH4_Ptot.best.df <- as.data.frame(cbind("Dam_type2" = c("Unfenced", "Fenced"), 10^(intervals(CH4_Ptot.best_plot)$fixed))) # Data to plot


CH4_Ptot.best.plot <- ggplot(data = CH4_Ptot.best.df, aes(x=Dam_type2, y=as.numeric(est.))) +
  geom_bar(stat='identity', width=0.8, fill=c("darkred", "forestgreen"), alpha=0.7) +
  geom_errorbar(aes(x=Dam_type2, ymin=as.numeric(lower), ymax=as.numeric(upper)), width=0.1, colour="black", alpha=0.7, size=.5) +
  labs(x = "Farm dam management",
       y = expression(paste("Phosphorous (mg L"^-1, ")"))) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  geom_signif(comparisons=list(c("Unfenced", "Fenced")), annotations="*", textsize=10,
              y_position = 0.15, tip_length = 0, vjust=0.6)
CH4_Ptot.best.plot
ggsave("../Results/CH4_Ptot.best.plot.pdf", width = 5.24, height = 4.06)






# WINNER PLOT
	CH4_Ptot.plot_RawData_log10 <- 
		ggplot(data = MeanData_CH4_D, aes(x = Dam_type2)) +
	  geom_jitter(aes(y = 10^TotalP_raw_1), width = 0.2, alpha = 0.5) +
		#ggplot(data = CH4_Ptot.best.df, aes(x=Dam_type2, y=as.numeric(est.))) +
	  #geom_bar(stat='identity', width=0.8, fill=c("darkred","forestgreen"), alpha=0.7) +
	  geom_errorbar(data = CH4_Ptot.best.df, aes(ymin=as.numeric(lower), ymax=as.numeric(upper)), width=0.1, colour="red", alpha= 1, size= 2) +
		geom_point(data = CH4_Ptot.best.df, aes(y=as.numeric(est.)), colour="red", size= 8) +
		#stat_summary(fun = "mean", colour = "red", size = 2, geom = "point") +
		#stat_summary(fun.data = "mean_cl_boot", colour = "red", size = 2) +
	  labs(x = "Farm dam management",
       y = expression(paste("Phosphorous (mg L"^-1, ")"))) +
	  theme_classic() +
	  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
	  geom_signif(aes(y = 10^TotalP_raw_1), comparisons=list(c("Unfenced", "Fenced")), annotations="*", textsize=10,
	              y_position = 0.1, tip_length = 0, vjust=0.6) +
  scale_y_log10(expand = expansion(mult = c(0, .1)),
  						  breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(0.005, 2)) +
  annotation_logticks(sides = "l")



	CH4_Ptot.plot_RawData_log10b <- 
		ggplot(data = MeanData_CH4_D, aes(x = Dam_type2)) +
	  geom_jitter(aes(y = 10^TotalP_raw_1), width = 0.2, alpha = 0.2, size = 3) +
		#ggplot(data = CH4_Ptot.best.df, aes(x=Dam_type2, y=as.numeric(est.))) +
	  #geom_bar(stat='identity', width=0.8, fill=c("darkred","forestgreen"), alpha=0.7) +
	  #geom_errorbar(data = CH4_Ptot.best.df, aes(ymin=as.numeric(lower), ymax=as.numeric(upper)), width=0.1, colour="red", alpha= 1, size= 2) +
		#geom_point(data = CH4_Ptot.best.df, aes(y=as.numeric(est.)), colour="red", size= 8) +
		geom_pointrange(data = CH4_Ptot.best.df, aes(y=as.numeric(est.), ymin=as.numeric(lower), ymax=as.numeric(upper)), size = 2, fatten = 2) +
		#stat_summary(fun = "mean", colour = "red", size = 2, geom = "point") +
		#stat_summary(fun.data = "mean_cl_boot", colour = "red", size = 2) +
	  labs(x = "Farm dam management",
       y = expression(paste("Phosphorous (mg L"^-1, ")"))) +
	  theme_classic() +
	  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
	  geom_signif(aes(y = 10^TotalP_raw_1), comparisons=list(c("Unfenced", "Fenced")), annotations="*", textsize=10,
	              y_position = 0.1, tip_length = 0, vjust=0.6) +
  scale_y_log10(expand = expansion(mult = c(0, .1)),
  						  breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(0.005, 2)) +
  annotation_logticks(sides = "l")








# ----------------
# Model validation
# ----------------

if(FALSE){
scatterplot(residuals(CH4_Ptot.best, type = "pearson")~predict(CH4_Ptot.best))
scatterplot(residuals(CH4_Ptot.best, type = "pearson")~MeanData_CH4_D3$TotalP_1)
boxplot(10^MeanData_CH4_D3$TotalP_1 ~ MeanData_CH4_D3$Dam_type)
}










# ------------------------
# MANAGEMENT ON TEMP
# ------------------------
# Aim: test the effects of dam management on Temp
# Dep.Var: WaterTemp_1
# Exp.Var: Dam_type
# Mod.name: CH4_Temp
# Data: MeanData_CH4_D

CH4_Temp.full1 = lme(WaterTemp_1 ~  Dam_type,
				random = ~1|Owner,
				#method = "ML",
				data = MeanData_CH4_D)

CH4_Temp.full2 = lme(WaterTemp_1 ~ Dam_type2,
				random = ~1|Owner,
				weights = varIdent(form = ~1|Dam_type2),
				method = "ML",
				data = MeanData_CH4_D)


AICc(CH4_Temp.full1, CH4_Temp.full2)

CH4_Temp.best = CH4_Temp.full2

plot(CH4_Temp.best)
Anova(CH4_Temp.best)
plot(allEffects(CH4_Temp.best))

capture.output(Anova(CH4_Temp.best, type = "III"), file = "../Results/CH4_Temp.best.txt")
save(CH4_Temp.best, file = "../Results/CH4_Temp.best.RData")


CH4_Temp.best_plot = lme(WaterTemp_raw_1 ~ -1 + Dam_type2,
				random = ~1|Owner,
				weights = varIdent(form = ~1|Dam_type2),
				data = MeanData_CH4_D3)

CH4_Temp.best.df <- as.data.frame(cbind("Dam_type2" = c("Unfenced", "Fenced"), (intervals(CH4_Temp.best_plot)$fixed))) # Data to plot


CH4_Temp.best.plot <- ggplot(data = CH4_Temp.best.df, aes(x=Dam_type2, y=as.numeric(est.))) +
  geom_bar(stat='identity', width=0.8, fill=c("darkred", "forestgreen"), alpha=0.7) +
  geom_errorbar(aes(x=Dam_type2, ymin=as.numeric(lower), ymax=as.numeric(upper)), width=0.1, colour="black", alpha=0.7, size=.5) +
  labs(x = "Farm dam management",
       y = expression("Temperature ("*~degree*C*")")) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) #+
  #geom_signif(comparisons=list(c("Unfenced", "Fenced")), annotations="*", textsize=10,
  #            y_position = 0.15, tip_length = 0, vjust=0.6)
CH4_Temp.best.plot
ggsave("../Results/CH4_Temp.best.plot.pdf", width = 5.24, height = 4.06)


# WINNER PLOT
CH4_Temp.plot_RawData_log10 <- 
		ggplot(data = MeanData_CH4_D, aes(x = Dam_type2)) +
	  geom_jitter(aes(y = WaterTemp_raw_1), width = 0.2, alpha = 0.5) +
		#ggplot(data = CH4_Temp.best.df, aes(x=Dam_type2, y=as.numeric(est.))) +
	  #geom_bar(stat='identity', width=0.8, fill=c("darkred","forestgreen"), alpha=0.7) +
	  geom_errorbar(data = CH4_Temp.best.df, aes(ymin=as.numeric(lower), ymax=as.numeric(upper)), width=0.1, colour="red", alpha= 1, size= 2) +
		geom_point(data = CH4_Temp.best.df, aes(y=as.numeric(est.)), colour="red", size= 8) +
		#stat_summary(fun = "mean", colour = "red", size = 2, geom = "point") +
		#stat_summary(fun.data = "mean_cl_boot", colour = "red", size = 2) +
	  labs(x = "Farm dam management",
       y = expression("Temperature ("*~degree*C*")")) +
	  theme_classic() +
	  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
	  

CH4_Temp.plot_RawData_log10b <- 
		ggplot(data = MeanData_CH4_D, aes(x = Dam_type2)) +
	  geom_jitter(aes(y = WaterTemp_raw_1), width = 0.2, alpha = 0.2, size = 3) +
		#ggplot(data = CH4_Temp.best.df, aes(x=Dam_type2, y=as.numeric(est.))) +
	  #geom_bar(stat='identity', width=0.8, fill=c("darkred","forestgreen"), alpha=0.7) +
	  #geom_errorbar(data = CH4_Temp.best.df, aes(ymin=as.numeric(lower), ymax=as.numeric(upper)), width=0.1, colour="red", alpha= 1, size= 2) +
		#geom_point(data = CH4_Temp.best.df, aes(y=as.numeric(est.)), colour="red", size= 8) +
		geom_pointrange(data = CH4_Temp.best.df, aes(y=as.numeric(est.), ymin=as.numeric(lower), ymax=as.numeric(upper)), size = 2, fatten = 2) +
		#stat_summary(fun = "mean", colour = "red", size = 2, geom = "point") +
		#stat_summary(fun.data = "mean_cl_boot", colour = "red", size = 2) +
	  labs(x = "Farm dam management",
       y = expression("Temperature ("*~degree*C*")")) +
	  theme_classic() +
	  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
	  







# ----------------
# Model validation
# ----------------

if(FALSE){
scatterplot(residuals(CH4_Temp.best, type = "pearson")~predict(CH4_Temp.best))
scatterplot(residuals(CH4_Temp.best, type = "pearson")~MeanData_CH4_D$WaterTemp_1)
boxplot(MeanData_CH4_D$WaterTemp_1~MeanData_CH4_D$Dam_type)
boxplot(MeanData_CH4_D$WaterTemp_1~MeanData_CH4_D$Dam_type2)
}

















# ------------------------
# MANAGEMENT ON CARBON STOCK
# ------------------------
# Aim: test the effects of dam management on carbon stock from CN analysis of cores
# Dep.Var: Stock_tCha_1
# Exp.Var: Dam_type
# Mod.name: CStock
# Data: MeanData_CH4_D


CStock.full1 = lme(Stock_tCha_1 ~  Dam_type2,
				random = ~1|Owner,
				weights = varIdent(form = ~1|Dam_type2),
				method = "ML",
				data = MeanData_CH4_D, na.action = na.exclude)

CStock.full2 = lme(Stock_tCha_1 ~  Dam_type2,
				random = ~1|Owner,
				method = "ML",
				data = MeanData_CH4_D, na.action = na.exclude)


# See if there are some covariances to correct for
CStock.full3 = lme(Stock_tCha_1 ~ Dam_type2,
				 random = ~ 1|Owner,
				 method = "ML",
                 data=MeanData_CH4_D, na.action = na.exclude,
                 correlation = corSpatial(form = ~ Lon_1 + Lat_1, type='gaussian'))

# Check covariance
CStock.full4 = lme(Stock_tCha_1 ~  Dam_type2 + NDVI.veg_1,
				random = ~1|Owner,
				method = "ML",
				data = MeanData_CH4_D, na.action = na.exclude)



AICc(CStock.full1, CStock.full2,CStock.full3,CStock.full4)
Anova(CStock.full2)

CStock.best = CStock.full2


# No effects
plot(allEffects(CStock.best))
summary(CStock.best)

# Have a look at the data
boxplot(data = MeanData_CH4_D, Stock_tCha_1 ~ Dam_type2)
boxplot(data = MeanData_CH4_D, Stock_tCha_1 ~ interaction(Dam_type2,Owner))
MeanData_CH4_D %>% dplyr::group_by(Dam_type2, Owner) %>%
dplyr::summarise(Stock_tCha_1 = mean(Stock_tCha_1)) %>%
ggplot(aes(x = Dam_type2, y = Stock_tCha_1, col = Owner, group = Owner)) +
geom_line()


plot(CStock.best)
Anova(CStock.best)

capture.output(Anova(CStock.best, type = "III"), file = "../Results/CStock.best.txt")
save(CStock.best, file = "../Results/CStock.best.RData")



CStock.full1_plot = lme(Stock_tCha_raw_1 ~  -1 + Dam_type2,
				random = ~1|Owner,
				method = "REML",
				data = MeanData_CH4_D, na.action = na.exclude)

Anova(CStock.full1_plot, type=3)


Data.CStock.full <-  as.data.frame(cbind("Dam_type2" = c("Unfenced", "Fenced"), 10^intervals(CStock.full1_plot)$fixed)) # Data to plot



CStock.full.plot <- ggplot(data = Data.CStock.full, aes(x=Dam_type2, y=as.numeric(est.))) +
  geom_bar(stat='identity', width=0.8, fill=c("darkred", "forestgreen"), alpha=0.7) +
  geom_errorbar(aes(x=Dam_type2, ymin=as.numeric(lower), ymax=as.numeric(upper)), width=0.1, colour="black", alpha=0.7, size=.5) +
  labs(x = "Farm dam management",
       y = expression(paste("Organic carbon in soil (tC ha"^-1, ")"))) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))# +
  #geom_signif(comparisons=list(c("Unfenced", "Fenced")), annotations="*", textsize=10,
  #            y_position = 0.56, tip_length = 0, vjust=0.6)

ggsave("../Results/CStock.full.plot.pdf", width = 5.24, height = 4.06)




# WINNER PLOT
	CH4_CStock.plot_RawData_log10 <- 
		ggplot(data = MeanData_CH4_D, aes(x = Dam_type2)) +
	  geom_jitter(aes(y = 10^Stock_tCha_raw_1), width = 0.2, alpha = 0.5) +
		#ggplot(data = Data.CStock.full, aes(x=Dam_type2, y=as.numeric(est.))) +
	  #geom_bar(stat='identity', width=0.8, fill=c("darkred","forestgreen"), alpha=0.7) +
	  geom_errorbar(data = Data.CStock.full, aes(ymin=as.numeric(lower), ymax=as.numeric(upper)), width=0.1, colour="red", alpha= 1, size= 2) +
		geom_point(data = Data.CStock.full, aes(y=as.numeric(est.)), colour="red", size= 8) +
		#stat_summary(fun = "mean", colour = "red", size = 2, geom = "point") +
		#stat_summary(fun.data = "mean_cl_boot", colour = "red", size = 2) +
	  labs(x = "Farm dam management",
       y = expression(paste("Organic carbon in soil (tC ha"^-1, ")"))) +
	  theme_classic() +
	  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_y_log10(expand = expansion(mult = c(0, .1)),
  						  breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+#,
  #              limits = c(0.005, 2)) +
  annotation_logticks(sides = "l")



	CH4_CStock.plot_RawData_log10b <- 
		ggplot(data = MeanData_CH4_D, aes(x = Dam_type2)) +
	  geom_jitter(aes(y = 10^Stock_tCha_raw_1), width = 0.2, alpha = 0.2, size = 3) +
		#ggplot(data = Data.CStock.full, aes(x=Dam_type2, y=as.numeric(est.))) +
	  #geom_bar(stat='identity', width=0.8, fill=c("darkred","forestgreen"), alpha=0.7) +
	  #geom_errorbar(data = Data.CStock.full, aes(ymin=as.numeric(lower), ymax=as.numeric(upper)), width=0.1, colour="red", alpha= 1, size= 2) +
		#geom_point(data = Data.CStock.full, aes(y=as.numeric(est.)), colour="red", size= 8) +
		geom_pointrange(data = Data.CStock.full, aes(y=as.numeric(est.), ymin=as.numeric(lower), ymax=as.numeric(upper)), size = 2, fatten = 2) +
		#stat_summary(fun = "mean", colour = "red", size = 2, geom = "point") +
		#stat_summary(fun.data = "mean_cl_boot", colour = "red", size = 2) +
	  labs(x = "Farm dam management",
       y = expression(paste("Organic carbon in soil (tC ha"^-1, ")"))) +
	  theme_classic() +
	  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_y_log10(expand = expansion(mult = c(0, .1)),
  						  breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+#,
  #              limits = c(0.005, 2)) +
  annotation_logticks(sides = "l")












# Summary of our results:
CH4_Ntot.best.df
CH4_Ptot.best.df
Data.CH4_Oxy.best_tr
Data.CH4_D.full1b.trs
Data.CO2_D.full3b.trs
Data.CStock.full

# Calculate percentage-change
PerChange = function(x){
	#x = Data.CH4_D.full1b.trs
	x[,"est."] = as.numeric(x[,"est."])
	round((subset(x, Dam_type2 == "Fenced")$est. - subset(x, Dam_type2 == "Unfenced")$est.)/subset(x, Dam_type2 == "Unfenced")$est.,2)
}

PerChange(CH4_Ntot.best.df)
PerChange(CH4_Ptot.best.df)
PerChange(Data.CH4_Oxy.best_tr)
PerChange(Data.CH4_D.full1b.trs)
PerChange(Data.CO2_D.full3b.trs)
PerChange(Data.CStock.full)




# Plot all together

AllEffectPlots = plot_grid(CH4_D.full1b.plot,CO2_D.full3b.plot,CO2eq_D.best.plot,CStock.full.plot,CH4_Oxy.best.plot,CH4_Ntot.best.plot,CH4_Ptot.best.plot,CH4_Temp.best.plot, ncol = 4)
save_plot(AllEffectPlots)




AllEffectPlots_RawData = plot_grid(
							# Carbon
							CH4_D.full1b.plot_RawData_log10,
							CO2_D.full1b.plot_RawData_log10,
							CO2eq_D.full1b.plot_RawData_log10,
							CH4_CStock.plot_RawData_log10,
							# Water quality
							CH4_Ntot.plot_RawData_log10,
							CH4_Ptot.plot_RawData_log10,
							CH4_Oxy.plot_RawData_log10,
							CH4_Temp.plot_RawData_log10,
							ncol = 2)

# Divided in two plots
AllEffectPlots_RawData_carbon = plot_grid(
							CH4_D.full1b.plot_RawData_log10,
							CO2_D.full1b.plot_RawData_log10,
							CO2eq_D.full1b.plot_RawData_log10,
							CH4_CStock.plot_RawData_log10,
							ncol = 2, 
							labels = c('A', 'B', 'C', 'D'),
							align="hv")
save_plot(filename = "../Results/AllEffectPlots_RawData_carbon.pdf", plot = AllEffectPlots_RawData_carbon, base_width = 9.9, base_height = 8.5)

AllEffectPlots_RawData_wq = plot_grid(
							CH4_Ntot.plot_RawData_log10,
							CH4_Ptot.plot_RawData_log10,
							CH4_Oxy.plot_RawData_log10,
							CH4_Temp.plot_RawData_log10,
							ncol = 2, 
							labels = c('A', 'B', 'C', 'D'),
							align="hv")
save_plot(filename = "../Results/AllEffectPlots_RawData_wq.pdf", plot = AllEffectPlots_RawData_wq, base_width = 9.9, base_height = 8.5)




# Use wisker plot
AllEffectPlots_RawData_carbon_b = plot_grid(
							CH4_D.full1b.plot_RawData_log10b,
							CO2_D.full1b.plot_RawData_log10c,
							CO2eq_D.full1b.plot_RawData_log10b,
							CH4_CStock.plot_RawData_log10b,
							ncol = 2, 
							labels = c('A', 'B', 'C', 'D'),
							align="hv")
save_plot(filename = "../Results/AllEffectPlots_RawData_carbon_wisker.pdf", plot = AllEffectPlots_RawData_carbon_b, base_width = 9.9, base_height = 8.5)
save_plot(filename = "../Results/AllEffectPlots_RawData_carbon_wisker.png", plot = AllEffectPlots_RawData_carbon_b, base_width = 9.9, base_height = 8.5)




AllEffectPlots_RawData_wq_b = plot_grid(
							CH4_Ntot.plot_RawData_log10b,
							CH4_Ptot.plot_RawData_log10b,
							CH4_Oxy.plot_RawData_log10b,
							CH4_Temp.plot_RawData_log10b,
							ncol = 2,
							labels = c('A', 'B', 'C', 'D'),
							align="hv")
save_plot(filename = "../Results/AllEffectPlots_RawData_wq_wisker.pdf", plot = AllEffectPlots_RawData_wq_b, base_width = 9.9, base_height = 8.5)
save_plot(filename = "../Results/AllEffectPlots_RawData_wq_wisker.png", plot = AllEffectPlots_RawData_wq_b, base_width = 9.9, base_height = 8.5)







# ------------------------
# ENVIRONMENTAL VARIABLES
# ------------------------

# Aim: test the effects of environmental parameters on CH4
# Dep.Var: Flux_mgm2day
# Exp.Var: WaterTemp, Oxygen_mg, Conductivity, Area_M2_gmap, Rain, FA, NDVI.veg
# Mod.name: CH4_D.env
# Data: Data.lm_CH4_D




CO2eq_D.env.full_n10.averages1 = lme(Flux_mgm2day_CO2eq_D ~  			WaterTemp_1+
												Oxygen_mg_1+
												Stock_tCha_1+
												Area_M2_gmap_1+
												#FA_1+
												#NDVI.veg_1+
												TotalP_1+
												TotalN_1,#+
				random = ~1|Owner,
				#weights = varExp(),
				method = "ML",
				data = MeanData_CO2eq_D,
				na.action = na.exclude)

CO2eq_D.env.full_n10.averages2 = lme(Flux_mgm2day_CO2eq_D ~  			WaterTemp_1+
												Oxygen_mg_1+
												Stock_tCha_1+
												Area_M2_gmap_1+
												#FA_1+
												#NDVI.veg_1+
												TotalP_1+
												TotalN_1,#+
				random = ~1|Owner,
				weights = varExp(),
				method = "ML",
				data = MeanData_CO2eq_D,
				na.action = na.exclude)



CO2eq_D.env.full_n10.averages3 = lme(Flux_mgm2day_CO2eq_D ~  			
											  WaterTemp_1*Oxygen_mg_1+
												Stock_tCha_1+
												Area_M2_gmap_1+
												#FA_1+
												#NDVI.veg_1+
												TotalP_1+
												TotalN_1,#+
				random = ~1|Owner,
				weights = varExp(),
				method = "ML",
				data = MeanData_CO2eq_D,
				na.action = na.exclude)



CO2eq_D.env.full_n10.averages4 = lme(Flux_mgm2day_CO2eq_D ~  			
											  WaterTemp_1*Oxygen_mg_1+
												Stock_tCha_1+
												Area_M2_gmap_1+
												#FA_1+
												#NDVI.veg_1+
												TotalP_1+
												TotalN_1,#+
				random = ~1|Owner,
				weights = varPower(),
				method = "ML",
				data = MeanData_CO2eq_D,
				na.action = na.exclude)


CO2eq_D.env.full_n10.averages5 = lme(Flux_mgm2day_CO2eq_D ~  			WaterTemp_1+
												Oxygen_mg_1+
												Stock_tCha_1+
												#Area_M2_gmap_1+
												#FA_1+
												#NDVI.veg_1+
												TotalP_1+
												TotalN_1,#+
				random = ~1|Owner,
				weights = varExp(),
				method = "ML",
				data = MeanData_CO2eq_D,
				na.action = na.exclude)


# What about looking at CO2eq instead?
CO2eq_D.env.full_n10.averages6 = lme(Flux_mgm2day_CO2eq_D ~  			WaterTemp_1+
												Oxygen_mg_1+
												Stock_tCha_1+
												#Area_M2_gmap_1+
												#FA_1+
												#NDVI.veg_1+
												TotalP_1+
												TotalN_1,#+
				random = ~1|Owner,
				#weights = varExp(),
				method = "ML",
				data = MeanData_CO2eq_D,
				na.action = na.exclude)


# Remote stock
CO2eq_D.env.full_n10.averages7 = lme(Flux_mgm2day_CO2eq_D ~  			WaterTemp_1+
												Oxygen_mg_1+
												#Stock_tCha_1+
												#Area_M2_gmap_1+
												#FA_1+
												#NDVI.veg_1+
												TotalP_1+
												TotalN_1,#+
				random = ~1|Owner,
				weights = varExp(),
				method = "ML",
				data = MeanData_CO2eq_D,
				na.action = na.exclude)


AIC(CO2eq_D.env.full_n10.averages1, CO2eq_D.env.full_n10.averages2, CO2eq_D.env.full_n10.averages3, CO2eq_D.env.full_n10.averages4, CO2eq_D.env.full_n10.averages5, CO2eq_D.env.full_n10.averages6, CO2eq_D.env.full_n10.averages7)


# Choose overall best
CO2eq_D.env.best = CO2eq_D.env.full_n10.averages5


summary(lm(residuals(CO2eq_D.env.best, type = "pearson")~predict(CO2eq_D.env.best)))

vif(CO2eq_D.env.best)

Anova(CO2eq_D.env.best)
plot(allEffects(CO2eq_D.env.best))
capture.output(Anova(CO2eq_D.env.best, type = "III"), file = "../Results/CO2eq_D.env.best.txt")
save(CO2eq_D.env.best, file = "../Results/CO2eq_D.env.best.RData")

plot(CO2eq_D.env.best)
mean(resid(CO2eq_D.env.best))





# MODEL VALIDATION
if(FALSE){


MeanData_CO2eq_D$Resid = residuals(CO2eq_D.env.best, type = "pearson")
MeanData_CO2eq_D$Fitted = predict(CO2eq_D.env.best)

scatterplot(data = MeanData_CO2eq_D, Resid~Fitted)

# Verify the final model
scatterplot(data = MeanData_CO2eq_D, Resid~WaterTemp_1)
scatterplot(data = MeanData_CO2eq_D, Resid~TotalN_1)
scatterplot(data = MeanData_CO2eq_D, Resid~TotalP_1)
scatterplot(data = MeanData_CO2eq_D, Resid~Oxygen_mg_1)
scatterplot(data = MeanData_CO2eq_D, Resid~Flux_mgm2day_CO2eq_D)

boxplot(data = MeanData_CO2eq_D, Resid~Owner/SiteCode)
boxplot(data = MeanData_CO2eq_D, Resid~Owner) # This is quite ok
boxplot(data = MeanData_CO2eq_D, Resid~SiteCode) # This not so much...

}








# ------------
# 2) Plot all effects of the model
# ------------


# First, remove NA in Stock from the dataset because makes issues below
MeanData_CO2eq_D_noNA = MeanData_CO2eq_D[!is.na(MeanData_CO2eq_D$Stock_tCha_raw_1),]


# Fit an equivalent model with all raw values
CO2eq_D.env.best.plot = lme(Flux_mgm2day_CO2eq_D ~ 
															WaterTemp_raw_1+
															Oxygen_mg_raw_1+
															Stock_tCha_raw_1+
															TotalP_raw_1+
															TotalN_raw_1,
				random = ~1|Owner,
				weights = varExp(),
				method = "ML",
				data = MeanData_CO2eq_D_noNA
				)




# 1) Plot WaterTemp
CO2eq_D.WaterTemp.res = Effect("WaterTemp_raw_1", partial.residuals=T, CO2eq_D.env.best.plot)
WaterTemp.plot = plot(CO2eq_D.WaterTemp.res, smooth.residuals=F, colors = c('black', 'black'),
     xlab=expression("Temperature ("*~degree*C*")"), ylab=expression(CO[2]-eq ~ diffusion ~ (mg ~ m^-2 ~ day^-1)),
     main='')
pdf("../Results/CO2eq_D.WaterTemp.plot.pdf", width = 6.24, height = 5.06)
WaterTemp.plot
dev.off()




# 2) Plot of Oxygen_mg
CO2eq_D.Oxygen_mg.res = Effect("Oxygen_mg_raw_1", partial.residuals=T, CO2eq_D.env.best.plot)
Oxygen_mg.plot = plot(CO2eq_D.Oxygen_mg.res, smooth.residuals=F, colors = c('black', 'black'),
                      xlab=expression(Oxygen ~ (mg ~ L^-1)), ylab=expression(CO[2]-eq ~ diffusion ~ (mg ~ m^2 ~ day^-1)),
                      main='')
pdf("../Results/CO2eq_D.Oxygen_mg.plot.pdf", width = 6.24, height = 5.06)
Oxygen_mg.plot
dev.off()



# 3) Plot of Total N
CO2eq_D.TotalN.res = Effect("TotalN_raw_1", partial.residuals=T, CO2eq_D.env.best.plot)
TotalN.plot = plot(CO2eq_D.TotalN.res, smooth.residuals=F, colors = c('black', 'black'),
                      xlab=expression(Total ~ Nitrogen ~ (mg ~ L^-1)), ylab=expression(CO[2]-eq ~ diffusion ~ (mg ~ m^2 ~ day^-1)),
                      main='')
pdf("../Results/CO2eq_D.TotalN.plot.pdf", width = 6.24, height = 5.06)
TotalN.plot
dev.off()




# 4) Plot of Total P
CO2eq_D.TotalP.res = Effect("TotalP_raw_1", partial.residuals=T, CO2eq_D.env.best.plot)
TotalP.plot = plot(CO2eq_D.TotalP.res, smooth.residuals=F, colors = c('black', 'black'),
                      xlab=expression(Total ~ Phosphorous ~ (mg ~ L^-1)), ylab=expression(CO[2]-eq ~ diffusion ~ (mg ~ m^2 ~ day^-1)),
                      main='')
pdf("../Results/CO2eq_D.TotalP.plot.pdf", width = 6.24, height = 5.06)
TotalP.plot
dev.off()





# 5) Organic carbon
CO2eq_D.Stock_tCha_1.res = Effect("Stock_tCha_raw_1", partial.residuals=T, CO2eq_D.env.best.plot)
Stock_tCha_1.plot = plot(CO2eq_D.Stock_tCha_1.res, smooth.residuals=F, colors = c('black', 'black'),
                      xlab= expression(paste("Organic carbon in soil (tC ha"^-1, ")")),
                      ylab=expression(CO[2]-eq ~ diffusion ~ (mg ~ m^2 ~ day^-1)),
                      main='')
pdf("../Results/CO2eq_D.Stock_tCha_1.plot.pdf", width = 6.24, height = 5.06)
Stock_tCha_1.plot
dev.off()














# ------------------------
# ENVIRONMENTAL VARIABLES ON CH4
# ------------------------

# Aim: test the effects of environmental parameters on average CH4
# Dep.Var: Flux_mgm2day_1
# Exp.Var: Oxygen_mg_1,Conductivity,Area_M2_gmap,FA,NDVI.veg,TotalN_1,TotalP_1
# Mod.name: CH4_D.env.full_n10.averages
# Data: MeanData_CH4_D



CH4_D.env.full_n10.averages1 = lme(Flux_mgm2day_1 ~  			WaterTemp_1+
												Oxygen_mg_1+
												#Conductivity+
												#Area_M2_gmap+
												#FA+
												#NDVI.veg+
												TotalN_1+
												TotalP_1,
				random = ~1|Owner,
				#weights = varExp(),
				method = "ML",
				data = MeanData_CH4_D,
				na.action = na.exclude)

CH4_D.env.full_n10.averages2 = lme(Flux_mgm2day_1 ~  			WaterTemp_1+
												Oxygen_mg_1+
												#Conductivity+
												#Area_M2_gmap+
												#FA+
												#NDVI.veg+
												TotalN_1+
												TotalP_1,
				random = ~1|Owner,
				weights = varExp(),
				method = "ML",
				data = MeanData_CH4_D,
				na.action = na.exclude)

CH4_D.env.full_n10.averages3 = lme(Flux_mgm2day_1 ~  			WaterTemp_1+
												Oxygen_mg_1+
												#Conductivity+
												#Area_M2_gmap+
												#FA+
												#NDVI.veg+
												TotalN_1+
												TotalP_1,
				random = ~1|Owner,
				weights = varPower(),
				method = "ML",
				data = MeanData_CH4_D,
				na.action = na.exclude,
				control = lmeControl(msTol = 1e-5, msVerbose = TRUE, maxIter = 5000, msMaxIter = 5000, opt = "nlminb"))


CH4_D.env.full_n10.averages4 = lme(Flux_mgm2day_1 ~  			WaterTemp_1+
												Oxygen_mg_1+
												Stock_tCha_1+
												#Conductivity+
												#Area_M2_gmap+
												#FA+
												#NDVI.veg+
												TotalN_1+
												TotalP_1,
				random = ~1|Owner,
				weights = varExp(),
				method = "ML",
				data = MeanData_CH4_D,
				na.action = na.exclude,
				control = lmeControl(msTol = 1e-5, msVerbose = TRUE, maxIter = 5000, msMaxIter = 5000, opt = "nlminb"))


AIC(CH4_D.env.full_n10.averages1,CH4_D.env.full_n10.averages2,CH4_D.env.full_n10.averages3,CH4_D.env.full_n10.averages4)




# Choose overall best
CH4_D.env.best = CH4_D.env.full_n10.averages4
Anova(CH4_D.env.best)



Anova(CH4_D.env.best)
plot(allEffects(CH4_D.env.best))
capture.output(Anova(CH4_D.env.best, type = "III"), file = "../Results/CH4_D.env.best.txt")
save(CH4_D.env.best, file = "../Results/CH4_D.env.best.RData")

plot(CH4_D.env.best)
mean(resid(CH4_D.env.best))




# MODEL VALIDATION
if(FALSE){


MeanData_CH4_D$Resid = residuals(CH4_D.env.best, type = "pearson")
MeanData_CH4_D$Fitted = predict(CH4_D.env.best)

scatterplot(data = MeanData_CH4_D, Resid~Fitted)

# Verify the final model
scatterplot(data = MeanData_CH4_D, Resid~WaterTemp_1)
scatterplot(data = MeanData_CH4_D, Resid~TotalN_1)
scatterplot(data = MeanData_CH4_D, Resid~TotalP_1)
scatterplot(data = MeanData_CH4_D, Resid~Oxygen_mg_1)
scatterplot(data = MeanData_CH4_D, Resid~Flux_mgm2day_CO2eq_D)

boxplot(data = MeanData_CH4_D, Resid~Owner/SiteCode)
boxplot(data = MeanData_CH4_D, Resid~Owner)
boxplot(data = MeanData_CH4_D, Resid~SiteCode)

}





# Fit an equivalent model with all raw values
CH4_D.env.best.plot = lme(Flux_mgm2day_1 ~  			WaterTemp_raw_1+
												Oxygen_mg_raw_1+
												Stock_tCha_raw_1+
												#Conductivity+
												#Area_M2_gmap+
												#FA+
												#NDVI.veg+
												TotalN_raw_1+
												TotalP_raw_1,
				random = ~1|Owner,
				weights = varExp(),
				method = "ML",
				data = na.omit(MeanData_CH4_D),
				#na.action = na.exclude,
				control = lmeControl(msTol = 1e-5, msVerbose = TRUE, maxIter = 5000, msMaxIter = 5000, opt = "nlminb"))


# Check it is equivalent to model with scaled&centered variables
round(AIC(CH4_D.env.best.plot),3) == round(AIC(CH4_D.env.best),3)
# Almost there... this missmatch is due to the tollerance











# ------------------------
# ENVIRONMENTAL VARIABLES ON CO2
# ------------------------

# Aim: test the effects of environmental parameters on average CO2
# Dep.Var: Flux_mgm2day_1
# Exp.Var: Oxygen_mg_1,Conductivity,Area_M2_gmap,FA,NDVI.veg,TotalN_1,TotalP_1
# Mod.name: CO2_D.env.full_n10.averages
# Data: MeanData_CO2_D


# There is a very low value that looks weird
MeanData_CO2_D2 = MeanData_CO2_D

# Check what happens if you analyse averages
CO2_D.env.full_n10.averages1 = lme(Flux_mgm2day_1 ~  			WaterTemp_1+
												Oxygen_mg_1+
												Stock_tCha_1+
												#Conductivity+
												#Area_M2_gmap+
												#FA+
												#NDVI.veg+
												TotalN_1+
												TotalP_1,
				random = ~1|Owner,
				#weights = varExp(),
				method = "ML",
				data = MeanData_CO2_D2,
				na.action = na.exclude)

CO2_D.env.full_n10.averages2 = lme(Flux_mgm2day_1 ~  			WaterTemp_1+
												Oxygen_mg_1+
												Stock_tCha_1+
												#Conductivity+
												#Area_M2_gmap+
												#FA+
												#NDVI.veg+
												TotalN_1+
												TotalP_1,
				random = ~1|Owner,
				weights = varExp(),
				method = "ML",
				data = MeanData_CO2_D2,
				na.action = na.exclude)

CO2_D.env.full_n10.averages3 = lme(Flux_mgm2day_1 ~  			WaterTemp_1+
												Oxygen_mg_1+
												Stock_tCha_1+
												Area_M2_gmap_1+
												#FA_1+
												#NDVI.veg_1+
												TotalPO_1+
												TotalN_1,#+
				random = ~1|Owner,
				weights = varPower(),
				method = "ML",
				data = MeanData_CO2_D2,
				na.action = na.exclude)

# Add stock
CO2_D.env.full_n10.averages4 = lme(Flux_mgm2day_1 ~  			WaterTemp_1+
												Oxygen_mg_1+
												Stock_tCha_1+
												Area_M2_gmap_1+
												#FA_1+
												#NDVI.veg_1+
												TotalPO_1+
												TotalN_1,#+
				random = ~1|Owner,
				#weights = varExp(),
				method = "ML",
				data = MeanData_CO2_D2,
				na.action = na.exclude)


AIC(CO2_D.env.full_n10.averages1,CO2_D.env.full_n10.averages2,CO2_D.env.full_n10.averages3,CO2_D.env.full_n10.averages4)




# Choose overall best
CO2_D.env.best = CO2_D.env.full_n10.averages1
Anova(CO2_D.env.best)

plot(allEffects(CO2_D.env.best))


Anova(CO2_D.env.best)
plot(allEffects(CO2_D.env.best))
capture.output(Anova(CO2_D.env.best, type = "III"), file = "../Results/CO2_D.env.best.txt")
save(CO2_D.env.best, file = "../Results/CO2_D.env.best.RData")

plot(CO2_D.env.best)
mean(na.omit(resid(CO2_D.env.best)))




# MODEL VALIDATION
if(FALSE){


MeanData_CO2_D2$Resid = residuals(CO2_D.env.best, type = "pearson")
MeanData_CO2_D2$Fitted = predict(CO2_D.env.best)

scatterplot(data = MeanData_CO2_D2, Resid~Fitted)

# Verify the final model
scatterplot(data = MeanData_CO2_D2, Resid~WaterTemp_1)
scatterplot(data = MeanData_CO2_D2, Resid~TotalN_1)
scatterplot(data = MeanData_CO2_D2, Resid~TotalP_1)
scatterplot(data = MeanData_CO2_D2, Resid~Oxygen_mg_1)
scatterplot(data = MeanData_CO2_D2, Resid~Flux_mgm2day_CO2eq_D)

boxplot(data = MeanData_CO2_D2, Resid~Owner/SiteCode)
boxplot(data = MeanData_CO2_D2, Resid~Owner)
boxplot(data = MeanData_CO2_D2, Resid~SiteCode)

}




# Fit an equivalent model with all raw values
CO2_D.env.best.plot = lme(Flux_mgm2day_1 ~ 
															WaterTemp_raw_1+
															Oxygen_mg_raw_1+
															Stock_tCha_raw_1+
															TotalP_raw_1+
															TotalN_raw_1,
				random = ~1|Owner,
				#weights = varExp(),
				method = "ML",
				data = na.omit(MeanData_CO2_D2)
				#na.action = na.exclude
				)

# Check it is equivalent to model with scaled&centered variables
round(AIC(CO2_D.env.best.plot),4) == round(AIC(CO2_D.env.best),4)



CO2_D.env.best.plot2 = lmer(Flux_mgm2day_1 ~ 
															WaterTemp_raw_1+
															Oxygen_mg_raw_1+
															Stock_tCha_raw_1+
															TotalP_raw_1+
															TotalN_raw_1 + (1|Owner),
				data = na.omit(MeanData_CO2_D2)
				)

visreg(CO2_D.env.best.plot2)

















# ===================
# Calculate relative importance in the two models
# ===================

# The three models of interest
#CO2eq_D.env.best.plot
#CO2_D.env.best.plot
#CH4_D.env.best.plot

# Each variable has 100 replicates
Rep = 100



# MODEL: CO2eq_D.env.best.plot
# --------------------

# Get the real predictors
RealPredict_CO2eq = predict(CO2eq_D.env.best.plot, level = 0)


# Vector of variables to shuffle
Shuffle_CO2eq = c("Oxygen_mg_raw_1",
				  "Stock_tCha_raw_1",
				  "TotalN_raw_1")


cor_CO2eq = matrix(nrow = length(Shuffle_CO2eq), ncol = Rep)
row.names(cor_CO2eq) = Shuffle_CO2eq


if(FALSE){

for(v in Shuffle_CO2eq){
	#v = Shuffle_CO2eq[1]

	# Fit the model
	for(r in 1:Rep){
	# r = 1

	# Create a copy of the data and shuffle the variable of interest	
	LoopData = CO2eq_D.env.best.plot$data
	LoopData[,as.character(v)] = LoopData[sample(dim(LoopData)[1]),v]

	# Fit model with shuffled data
	Shuffle.lm = try(lme(formula(CO2eq_D.env.best.plot),
					random = ~1|Owner,
					#weights = varExp(),
					method = "ML",
					data = LoopData
					))

	if(class(Shuffle.lm) != "try-error") {

		Shuffle.predict = predict(Shuffle.lm, level = 0)
		cor_CO2eq[v,r] = 1-cor(RealPredict_CO2eq, Shuffle.predict)

	}

	rm(Shuffle.lm, Shuffle.predict, LoopData)

	print(paste0(round(sum(!is.na(as.numeric(cor_CO2eq)))/length(as.numeric(cor_CO2eq)),2)*100, "% completed"))

	}

}

save(cor_CO2eq, file = "../Results/Correlation_coef_CO2eq.RData")
write.csv(cor_CO2eq, file = "../Results/Correlation_coef_CO2eq.csv")

}











# MODEL: CH4_D.env.best.plot
# --------------------

# Get the real predictors
RealPredict_CH4 = predict(CH4_D.env.best.plot, level = 0)


# Vector of variables to shuffle
Shuffle_CH4 = c("WaterTemp_raw_1",
				  "Oxygen_mg_raw_1",
				  "Stock_tCha_raw_1",
				  "TotalN_raw_1",
				  "TotalP_raw_1")


cor_CH4 = matrix(nrow = length(Shuffle_CH4), ncol = Rep)
row.names(cor_CH4) = Shuffle_CH4


if(FALSE){

for(v in Shuffle_CH4){
	#v = Shuffle_CH4[1]

	# Fit the model
	for(r in 1:Rep){
	# r = 1

	# Create a copy of the data and shuffle the variable of interest	
	LoopData = CH4_D.env.best.plot$data
	LoopData[,as.character(v)] = LoopData[sample(dim(LoopData)[1]),v]

	# Fit model with shuffled data
	Shuffle.lm = try(lme(formula(CH4_D.env.best.plot),
					random = ~1|Owner,
					#weights = varExp(),
					method = "ML",
					data = LoopData
					))

	if(class(Shuffle.lm) != "try-error") {

		Shuffle.predict = predict(Shuffle.lm, level = 0)
		cor_CH4[v,r] = 1-cor(RealPredict_CH4, Shuffle.predict)

	}

	rm(Shuffle.lm, Shuffle.predict, LoopData)

	print(paste0(round(sum(!is.na(as.numeric(cor_CH4)))/length(as.numeric(cor_CH4)),2)*100, "% completed"))

	}

}


save(cor_CH4, file = "../Results/Correlation_coef_CH4.RData")
write.csv(cor_CH4, file = "../Results/Correlation_coef_CH4.csv")
}


load(file = "../Results/Correlation_coef_CH4.RData")
load(file = "../Results/Correlation_coef_CO2eq.RData")







# Melt the dataframe
cor.melt_CO2eq = melt(cor_CO2eq)
cor.melt_CH4 = melt(cor_CH4)

# Calculate average
cor.lm_CO2eq = lm(data = cor.melt_CO2eq, value~-1+Var1)
cor.lm_CH4 = lm(data = cor.melt_CH4, value~-1+Var1)


# DF for sensitivity analysis
# ----
Sensitivity_df_CO2eq = data.frame(Vars = names(coef(cor.lm_CO2eq)), Mean = percent(unname(coef(cor.lm_CO2eq)),accuracy = 0.1), LCI = percent(unname(confint(cor.lm_CO2eq))[,1], accuracy = 0.1), UCI = percent(unname(confint(cor.lm_CO2eq))[,2], accuracy = 0.1))
Sensitivity_df_CH4 = data.frame(Vars = names(coef(cor.lm_CH4)), Mean = percent(unname(coef(cor.lm_CH4)),accuracy = 0.1), LCI = percent(unname(confint(cor.lm_CH4))[,1], accuracy = 0.1), UCI = percent(unname(confint(cor.lm_CH4))[,2], accuracy = 0.1))
Sensitivity_df_CO2 = data.frame(Vars = "Var1Oxygen_mg_raw_1", Mean = percent(1), LCI = NA, UCI = NA)


# Add the model
Sensitivity_df_CO2eq$Model = "CO2eq"
Sensitivity_df_CH4$Model = "CH4"
Sensitivity_df_CO2$Model = "CO2"

# Merge together
Sensitivity_df = rbind(Sensitivity_df_CO2eq, Sensitivity_df_CH4, Sensitivity_df_CO2)


# Parse the labels
Sensitivity_df$Val_parsed_log = factor(Sensitivity_df$Vars, 
						levels = c("Var1Oxygen_mg_raw_1", "Var1TotalN_raw_1", "Var1Stock_tCha_raw_1", "Var1WaterTemp_raw_1", "Var1TotalP_raw_1"),
						labels = 
								# Create the parsed lables in the right order
								  c(expression(paste("Oxygen (mg L"^-1,")")),
										expression(paste("Nitrogen (log"[10], " N mg L"^-1,")")),
										expression(paste("C Stock (log"[10], " tons C ha"^-1,")")),
										expression("Water temp ("*~degree*C*")"),
										expression(paste("Phosphorous (log"[10], " P mg L"^-1,")"))
						))



# Parse the models
Sensitivity_df$Model_parsed_log = factor(Sensitivity_df$Model, 
	levels = c("CO2eq","CH4","CO2"),
	labels = c(
						expression(paste("CO"[2],"-eq flux (log"[10], " g m"^2," day"^-1," + 1.8)")),
				expression(paste("CH"[4]," flux (log"[10], " g m"^2," day"^-1," + 0.002)")),
				expression(paste("CO"[2]," flux (log"[10], " g m"^2," day"^-1," + 1.8)"))
				)
	)


Vars_all = c("Var1Oxygen_mg_raw_1", "Var1TotalN_raw_1", "Var1Stock_tCha_raw_1", "Var1WaterTemp_raw_1", "Var1TotalP_raw_1")
Models_all = c("CO2eq","CH4","CO2")

All_ModelVars = expand.grid(Model = Models_all,Vars = Vars_all)

Sensitivity_df_all = merge(All_ModelVars, Sensitivity_df, all = T)

write.csv(Sensitivity_df_all, file = "../Results/Correlation_coef_summary_all.csv")









# ======================= #
# PLOT ALL RESULTS IN ONE #
# ======================= #


Anova(CO2eq_D.env.best.plot)
Anova(CO2_D.env.best.plot)
Anova(CH4_D.env.best.plot)



# Create a function to get the dataset to plot the effects in arithmetic scale
# -------------------
PlotArithmetic = function(mod){
	# mod = CH4_D.env.best.plot

	Terms = names(fixef(mod))[-1]

	PredictedValues.all = c()
	Residuals.all = c()

	for(t in 1:length(Terms))	{
		#t = 1
		
		# Associate a term
		term = Terms[t]

		# Get the effects
		options("na.action")
		Effects = Effect(term, partial.residuals=T, mod, na.action = na.exclude)

		# Get the predicted values
		PredictedValues = as.data.frame(Effects)
		names(PredictedValues)[1] = "x"

		# Add the name of the character
		PredictedValues$Val = as.factor(term)

		# Create a new dataset with the partial residuals
		Residuals = data.frame("y" = c(predict(mod, level = 0) + resid(Effects)))
		Residuals$y_raw = Effects$data[,Effects$response]
		Residuals$x = as.numeric(unlist(c(mod$data[,term])))
		Residuals$Val = as.factor(term)

		# Residuals from this post: https://stackoverflow.com/questions/43950459/use-ggplot-to-plot-partial-effects-obtained-with-effects-library
		closest <- function(x, x0) apply(outer(x, x0, FUN=function(x, x0) abs(x - x0)), 1, which.min)
		x.fit <- unlist(Effects$x.all)
		x <- data.frame(lower = Effects$lower, upper = Effects$upper, fit = Effects$fit, term = Effects$x[,term])
		xy <- data.frame(x = x.fit, y = x$fit[closest(x.fit, x$term)] + Effects$residuals)
		Residuals$y_NewRes = xy$y

		# Merge all values together	
		PredictedValues.all = rbind(PredictedValues.all, PredictedValues)
		Residuals.all = rbind(Residuals.all, Residuals)

	}

	return(list(Predicted = PredictedValues.all, Residuals = Residuals.all))

}


# Run the function
PlotData_CO2eq = PlotArithmetic(CO2eq_D.env.best.plot)
PlotData_CO2 = PlotArithmetic(CO2_D.env.best.plot)
PlotData_CH4 = PlotArithmetic(CH4_D.env.best.plot)


# Add the info on which model
PlotData_CO2eq$Predicted$Model = "CO2eq"
PlotData_CO2eq$Residuals$Model = "CO2eq"

PlotData_CO2$Predicted$Model = "CO2"
PlotData_CO2$Residuals$Model = "CO2"

PlotData_CH4$Predicted$Model = "CH4"
PlotData_CH4$Residuals$Model = "CH4"


# Combine together
AllPlotData_Residuals = rbind(PlotData_CO2eq$Residuals,PlotData_CO2$Residuals,PlotData_CH4$Residuals)
AllPlotData_Predicted = rbind(PlotData_CO2eq$Predicted,PlotData_CO2$Predicted,PlotData_CH4$Predicted)


# Save p values
Pvalues_CO2eq = as.data.frame(Anova(CO2eq_D.env.best.plot))
Pvalues_CO2eq$levels = rownames(Pvalues_CO2eq)

Pvalues_CO2 = as.data.frame(Anova(CO2_D.env.best.plot))
Pvalues_CO2$levels = rownames(Pvalues_CO2)

Pvalues_CH4 = as.data.frame(Anova(CH4_D.env.best.plot))
Pvalues_CH4$levels = rownames(Pvalues_CH4)

# Add info on models
Pvalues_CO2eq$Model = "CO2eq"
Pvalues_CO2$Model = "CO2"
Pvalues_CH4$Model = "CH4"

# Combine
Pvalues = rbind(Pvalues_CO2eq,Pvalues_CO2,Pvalues_CH4)
rownames(Pvalues) = NULL

names(Pvalues)[3] = "p"
Pvalues$p_round = round(Pvalues$p,3)
Pvalues$p_round2 = ifelse(Pvalues$p_round < 0.01, "p < 0.01", paste0("p = ", Pvalues$p_round))



# Finally we identify and remove the non-sign. fitted lines
AllPlotData_Predicted$ModelVal = interaction(AllPlotData_Predicted$Model, AllPlotData_Predicted$Val)
ModelVal_NonSign = c("CO2eq.WaterTemp_raw_1", "CO2eq.TotalP_raw_1", "CO2.WaterTemp_raw_1", "CO2.Stock_tCha_raw_1", "CO2.TotalP_raw_1", "CO2.TotalN_raw_1", "CH4.WaterTemp_raw_1")
AllPlotData_Predicted$fit_2 = ifelse(AllPlotData_Predicted$ModelVal %in% ModelVal_NonSign, NA, AllPlotData_Predicted$fit)
AllPlotData_Predicted$lower_2 = ifelse(AllPlotData_Predicted$ModelVal %in% ModelVal_NonSign, NA, AllPlotData_Predicted$lower)
AllPlotData_Predicted$upper_2 = ifelse(AllPlotData_Predicted$ModelVal %in% ModelVal_NonSign, NA, AllPlotData_Predicted$upper)





# ----------------
# Labelling - Log scale
# ----------------

# Create the parsed lables in the right order
ParseLabeled_log = c(expression(paste("Oxygen (mg L"^-1,")")),
				expression(paste("Nitrogen (log"[10], " N mg L"^-1,")")),
				expression(paste("C Stock (log"[10], " tons C ha"^-1,")")),
				expression("Water temp ("*~degree*C*")"),
				expression(paste("Phosphorous (log"[10], " P mg L"^-1,")"))
		)

ParseModel_log = c(
				expression(paste("CO"[2],"-eq flux (log"[10], " g m"^2," day"^-1," + 1.8)")),
				expression(paste("CH"[4]," flux (log"[10], " g m"^2," day"^-1," + 0.002)")),
				expression(paste("CO"[2]," flux (log"[10], " g m"^2," day"^-1," + 1.8)"))
		)


# Label the residuals
AllPlotData_Residuals$Val_parsed_log = factor(AllPlotData_Residuals$Val, 
	levels = c("Oxygen_mg_raw_1","TotalN_raw_1","Stock_tCha_raw_1","WaterTemp_raw_1","TotalP_raw_1"),
	labels = ParseLabeled_log
	)

# Label the predicted values
AllPlotData_Predicted$Val_parsed_log = factor(AllPlotData_Predicted$Val, 
	levels = c("Oxygen_mg_raw_1","TotalN_raw_1","Stock_tCha_raw_1","WaterTemp_raw_1","TotalP_raw_1"),
	labels = ParseLabeled_log
	)

# Label the models
AllPlotData_Predicted$Model_parsed_log = factor(AllPlotData_Predicted$Model, 
	levels = c("CO2eq","CH4","CO2"),
	labels = ParseModel_log
	)

AllPlotData_Residuals$Model_parsed_log = factor(AllPlotData_Residuals$Model, 
	levels = c("CO2eq","CH4","CO2"),
	labels = ParseModel_log
	)

# Label the p-values
Pvalues$Val_parsed_log = factor(Pvalues$levels, 
	levels = c("Oxygen_mg_raw_1","TotalN_raw_1","Stock_tCha_raw_1","WaterTemp_raw_1","TotalP_raw_1"),
	labels = ParseLabeled_log
	)

Pvalues$Model_parsed_log = factor(Pvalues$Model, 
	levels = c("CO2eq","CH4","CO2"),
	labels = ParseModel_log
	)



# Turn from mg to g
AllEffectPlots_g_log10 = 	
	ggplot() +
	geom_point(data = AllPlotData_Residuals, aes(x = x, y = y_NewRes-3)) +
	geom_line(data = AllPlotData_Predicted, aes(x = x, y = fit_2-3)) +
	geom_ribbon(data = AllPlotData_Predicted, aes(x = x, ymin = lower_2-3, ymax = upper_2-3), alpha = 0.2) +
	facet_grid(Model_parsed_log~Val_parsed_log, scale = "free", labeller = labeller(Model_parsed_log = label_parsed, Val_parsed_log = label_parsed), switch="y") +
	theme_bw() +
	xlab("")+
	ylab("") +
	theme(axis.text = element_text(size = 12),
			axis.title = element_text(size = 14),
		  panel.grid.major = element_blank(),
		  panel.grid.minor = element_blank(),
			strip.text.x = element_text(size = 14),
			strip.text.y = element_text(size = 14)) +
	geom_text(data=Pvalues, aes(x = -Inf, y = -Inf, label=p_round2), parse = F, inherit.aes=FALSE,hjust = -0.25, vjust = -1, size = 6, colour = " dark grey") +
	scale_y_continuous(position = "right") +
	geom_text(data=subset(Sensitivity_df_all, !is.na(Mean)), aes(x = -Inf, y = -Inf, label=Mean), parse = F, inherit.aes=FALSE,hjust = -0.3, vjust = -2.5, size = 6, colour = " dark grey")

ggsave(AllEffectPlots_g_log10, file = "../Results/AllEffectPlots_g_log10_3mods.pdf", width = 15.6, height  = 11)
ggsave(AllEffectPlots_g_log10, file = "../Results/AllEffectPlots_g_log10_3mods.png", width = 15.6, height  = 11)






# ----------------
# Arithmetic scale
# ----------------


# Create the parsed lables in the right order
ParseLabeled = c(expression(paste("Oxygen (mg L"^-1,")")),
				expression(paste("Nitrogen (N mg L"^-1,")")),
				expression(paste("C Stock (tons C ha"^-1,")")),
				expression("Water temp ("*~degree*C*")"),
				expression(paste("Phosphorous (P mg L"^-1,")"))
		)

ParseModel = c(
				expression(paste("CO"[2],"-eq flux (g m"^2," day"^-1,")")),
				expression(paste("CH"[4]," flux (g m"^2," day"^-1,")")),
				expression(paste("CO"[2]," flux (g m"^2," day"^-1,")"))
		)


# Label the residuals
AllPlotData_Residuals$Val_parsed = factor(AllPlotData_Residuals$Val, 
	levels = c("Oxygen_mg_raw_1","TotalN_raw_1","Stock_tCha_raw_1","WaterTemp_raw_1","TotalP_raw_1"),
	labels = ParseLabeled
	)

# Label the predicted values
AllPlotData_Predicted$Val_parsed = factor(AllPlotData_Predicted$Val, 
	levels = c("Oxygen_mg_raw_1","TotalN_raw_1","Stock_tCha_raw_1","WaterTemp_raw_1","TotalP_raw_1"),
	labels = ParseLabeled
	)

# Label the models
AllPlotData_Predicted$Model_parsed = factor(AllPlotData_Predicted$Model, 
	levels = c("CO2eq","CH4","CO2"),
	labels = ParseModel
	)

AllPlotData_Residuals$Model_parsed = factor(AllPlotData_Residuals$Model, 
	levels = c("CO2eq","CH4","CO2"),
	labels = ParseModel
	)

# Label the p-values
Pvalues$Val_parsed = factor(Pvalues$levels, 
	levels = c("Oxygen_mg_raw_1","TotalN_raw_1","Stock_tCha_raw_1","WaterTemp_raw_1","TotalP_raw_1"),
	labels = ParseLabeled
	)

Pvalues$Model_parsed = factor(Pvalues$Model, 
	levels = c("CO2eq","CH4","CO2"),
	labels = ParseModel
	)



# Calculate the arithmetic scale for the dependent variable
AllPlotData_Predicted$fit_2_aritm = ifelse(AllPlotData_Predicted$Model == "CH4", ((10^AllPlotData_Predicted$fit_2) - 2) / 1000, ((10^AllPlotData_Predicted$fit_2) - 1800) / 1000)
AllPlotData_Predicted$lower_2_aritm = ifelse(AllPlotData_Predicted$Model == "CH4", ((10^AllPlotData_Predicted$lower_2) - 2) / 1000, ((10^AllPlotData_Predicted$lower_2) - 1800) / 1000)
AllPlotData_Predicted$upper_2_aritm = ifelse(AllPlotData_Predicted$Model == "CH4", ((10^AllPlotData_Predicted$upper_2) - 2) / 1000, ((10^AllPlotData_Predicted$upper_2) - 1800) / 1000)
AllPlotData_Residuals$y_raw_aritm = ifelse(AllPlotData_Residuals$Model == "CH4", ((10^AllPlotData_Residuals$y_raw) - 2) / 1000, ((10^AllPlotData_Residuals$y_raw) - 1800) / 1000)
AllPlotData_Residuals$y_NewRes_aritm = ifelse(AllPlotData_Residuals$Model == "CH4", ((10^AllPlotData_Residuals$y_NewRes) - 2) / 1000, ((10^AllPlotData_Residuals$y_NewRes) - 1800) / 1000)


# Now also convert the X-axis
AllPlotData_Residuals$x_aritm = with(AllPlotData_Residuals, ifelse(Val %in% c("Stock_tCha_raw_1","TotalP_raw_1","TotalN_raw_1"), 10^x, x))
AllPlotData_Predicted$x_aritm = with(AllPlotData_Predicted, ifelse(Val %in% c("Stock_tCha_raw_1","TotalP_raw_1","TotalN_raw_1"), 10^x, x))




# Plot Y on arithmetic scale
AllEffectPlots_g_aritm = 	
	ggplot() +
	geom_point(data = AllPlotData_Residuals, aes(x = x_aritm, y = y_NewRes_aritm)) +
	geom_line(data = AllPlotData_Predicted, aes(x = x_aritm, y = fit_2_aritm)) +
	geom_ribbon(data = AllPlotData_Predicted, aes(x = x_aritm, ymin = lower_2_aritm, ymax = upper_2_aritm), alpha = 0.2) +
	facet_grid(Model_parsed~Val_parsed, scale = "free", labeller = labeller(Model_parsed = label_parsed, Val_parsed = label_parsed), switch="y") +
	theme_bw() +
	scale_y_continuous(lim = c(-2, 18)) +
	xlab("")+
	ylab("") +
	#ylab(expression(paste("CO"[2],"-eq flux (mg m"^-2," day"^-1,")"))) +
	theme(axis.text = element_text(size = 12),
			axis.title = element_text(size = 14),
			strip.text.x = element_text(size = 14),
		  panel.grid.major = element_blank(),
		  panel.grid.minor = element_blank()) +
	geom_text(data=Pvalues, aes(x = -Inf, y = -Inf, label=p_round2), parse = F, inherit.aes=FALSE,hjust = -0.2, vjust = -1, size = 6, colour = " dark grey") +
	geom_hline(yintercept=0, linetype="dashed")

ggsave(AllEffectPlots_g_aritm, file = "../Results/AllEffectPlots_g_aritm.pdf", width = 11.9, height  = 6.5)
ggsave(AllEffectPlots_g_aritm, file = "../Results/AllEffectPlots_g_aritm.png", width = 11.9, height  = 6.5)





# Plot only the effects of Oxygen on CO2-eq

AllEffectPlots_g_aritm_O2 =
ggplot() +
geom_point(data = subset(AllPlotData_Residuals, Model == "CO2eq" & Val == "Oxygen_mg_raw_1"), aes(x = x_aritm, y = y_NewRes_aritm)) +
geom_line(data = subset(AllPlotData_Predicted, Model == "CO2eq" & Val == "Oxygen_mg_raw_1"), aes(x = x_aritm, y = fit_2_aritm)) +
geom_ribbon(data = subset(AllPlotData_Predicted, Model == "CO2eq" & Val == "Oxygen_mg_raw_1"), aes(x = x_aritm, ymin = lower_2_aritm, ymax = upper_2_aritm), alpha = 0.2) +
#facet_grid(Model_parsed~Val_parsed, scale = "free", labeller = labeller(Model_parsed = label_parsed, Val_parsed = label_parsed), switch="y") +
theme_bw() +
scale_y_continuous(lim = c(-2, 11)) +
xlab(expression(paste("Oxygen (mg L"^-1, ")")))+
ylab(expression(paste("CO"[2], "-eq flux (g m"^2, " day"^-1, ")"))) +
theme(axis.text = element_text(size = 12),
			axis.title = element_text(size = 14),
			strip.text.x = element_text(size = 14),
		  panel.grid.major = element_blank(),
		  panel.grid.minor = element_blank()) +
#geom_text(data=Pvalues, aes(x = -Inf, y = -Inf, label=p_round2), parse = F, inherit.aes=FALSE,hjust = -0.2, vjust = -1, size = 6, colour = " dark grey") +
geom_hline(yintercept=0, linetype="dashed")

ggsave(AllEffectPlots_g_aritm_O2, file = "../Results/AllEffectPlots_g_aritm_O2.pdf", width = 6.146341, height  = 6.609756)
ggsave(AllEffectPlots_g_aritm_O2, file = "../Results/AllEffectPlots_g_aritm_O2.png", width = 6.146341, height  = 6.609756)





# Results reported in the MS


# -----
# CH4
# -----

# Doubling oxygen from 5 to 10 
Oxygen_5_CH4 = with(subset(AllPlotData_Predicted, Val == "Oxygen_mg_raw_1" & Model == "CH4"), fit_2_aritm[which.min(abs(x_aritm - 5))])
Oxygen_10_CH4 = with(subset(AllPlotData_Predicted, Val == "Oxygen_mg_raw_1" & Model == "CH4"), fit_2_aritm[which.min(abs(x_aritm - 10))])
(Oxygen_5_CH4-Oxygen_10_CH4)/Oxygen_5_CH4


# -----
# CO2
# -----

# Doubling oxygen from 5 to 10 
Oxygen_5_CO2 = with(subset(AllPlotData_Predicted, Val == "Oxygen_mg_raw_1" & Model == "CO2"), fit_2_aritm[which.min(abs(x_aritm - 5))])
Oxygen_10_CO2 = with(subset(AllPlotData_Predicted, Val == "Oxygen_mg_raw_1" & Model == "CO2"), fit_2_aritm[which.min(abs(x_aritm - 10))])
(Oxygen_5_CO2-Oxygen_10_CO2)/Oxygen_5_CO2


# -----
# CO2eq
# -----

# Doubling oxygen from 5 to 10 
Oxygen_5_CO2eq = with(subset(AllPlotData_Predicted, Val == "Oxygen_mg_raw_1" & Model == "CO2eq"), fit_2_aritm[which.min(abs(x_aritm - 5))])
Oxygen_10_CO2eq = with(subset(AllPlotData_Predicted, Val == "Oxygen_mg_raw_1" & Model == "CO2eq"), fit_2_aritm[which.min(abs(x_aritm - 10))])
(Oxygen_5_CO2eq-Oxygen_10_CO2eq)/Oxygen_5_CO2eq












  


# ===================
# Relationship between oxygen and pH
# ===================

# Get data on pH and Oxygem
MeanData_pH_a = Data.lm_CO2_D %>% 
		   dplyr::select(WaterTemp, Oxygen_mg_raw, Oxygen_mg, Stock_tCha, TotalN,TotalP, pH_raw, SiteCode, Dam_type, Dam_type2, Owner, TotalN_raw,TotalP_raw) %>%
		   dplyr::group_by(SiteCode, Dam_type, Dam_type2, Owner) %>%
		   dplyr::summarise(across(everything(), list(mean)))

# Add fluxes of CO2 and CH4
MeanData_pH = merge(MeanData_pH_a, 
	MeanData_CO2eq_D %>% select(SiteCode,Flux_mgm2day_CO2_D,Flux_mgm2day_CH4_D_tr,Flux_mgm2day_CO2_D_tr,Flux_mgm2day_CO2eq_D_tr,Flux_mgm2day_CH4toCO2_D_tr,Flux_mgm2day_CO2eq_D))
names(MeanData_pH)


VIFtest = lme(Flux_mgm2day_CO2eq_D ~  			WaterTemp_1+
												Oxygen_mg_1+
												Stock_tCha_1+
												#Area_M2_gmap_1+
												#FA_1+
												#NDVI.veg_1+
												TotalP_1+
												TotalN_1 +
												pH_raw_1,#+
				random = ~1|Owner,
				weights = varExp(),
				method = "ML",
				data = MeanData_pH,
				na.action = na.exclude)

vif(VIFtest)

# No effect of management nor interaction
Oxygen.lm1 = lm(data = MeanData_pH, pH_raw_1~Oxygen_mg_raw_1*Dam_type2)
summary(Oxygen.lm1)
#plot(Oxygen.lm1)
plot(allEffects(Oxygen.lm1))

Oxygen.lm2 = lm(data = MeanData_pH, pH_raw_1~Dam_type2)
summary(Oxygen.lm2)
#plot(Oxygen.lm2)
plot(allEffects(Oxygen.lm2))



# Plot Oxygen vs pH
cor.test(MeanData_pH$Oxygen_mg_raw_1, MeanData_pH$pH_raw_1)

pH_Oxygen =
ggplot(data = MeanData_pH, aes(y = Oxygen_mg_raw_1, x = pH_raw_1)) +
geom_point() +
geom_smooth(method = "lm", se = F, col = "black") +
theme_bw() +
ylab(expression(paste("Oxygen (mg L"^-1, ")")))+
xlab("pH") +
theme(axis.text = element_text(size = 12),
			axis.title = element_text(size = 14),
			strip.text.x = element_text(size = 14),
		  panel.grid.major = element_blank(),
		  panel.grid.minor = element_blank()) +
annotate("text", x = 6.8, y = 18, label = "r = 0.72", colour = " dark grey", size = 6) +
annotate("text", x = 6.85, y = 17, label = "p < 0.001", colour = " dark grey", size = 6)


ggsave(pH_Oxygen, file = "../Results/pH_Oxygen.pdf", width = 6.146341, height  = 6.609756)
ggsave(pH_Oxygen, file = "../Results/pH_Oxygen.png", width = 6.146341, height  = 6.609756)



# Plot CO2 vs pH
cor.test(log10(MeanData_pH$Flux_mgm2day_CO2_D_tr+1800)-3, MeanData_pH$pH_raw_1)

pH_CO2 =
ggplot(data = MeanData_pH, aes(y = log10(Flux_mgm2day_CO2_D_tr+1800)-3, x = pH_raw_1)) +
geom_point() +
geom_smooth(method = "lm", se = F, col = "black") +
theme_bw() +
ylab(expression(paste("CO"[2]," flux (log"[10], " g m"^2," day"^-1," + 1.8)")))+
xlab("pH") +
theme(axis.text = element_text(size = 12),
			axis.title = element_text(size = 14),
			strip.text.x = element_text(size = 14),
		  panel.grid.major = element_blank(),
		  panel.grid.minor = element_blank()) +
annotate("text", x = 6.8, y = -0.7, label = "r = - 0.76", colour = " dark grey", size = 6) +
annotate("text", x = 6.82, y = -0.85, label = "p < 0.001", colour = " dark grey", size = 6)

ggsave(pH_CO2, file = "../Results/pH_CO2.pdf", width = 6.146341, height  = 6.609756)
ggsave(pH_CO2, file = "../Results/pH_CO2.png", width = 6.146341, height  = 6.609756)








# Plot CO2eq vs pH
cor.test(log10(MeanData_pH$Flux_mgm2day_CO2eq_D_tr+1800)-3, MeanData_pH$pH_raw_1)

pH_CO2eq =
ggplot(data = MeanData_pH, aes(y = log10(Flux_mgm2day_CO2eq_D_tr+1800)-3, x = pH_raw_1)) +
geom_point() +
geom_smooth(method = "lm", se = F, col = "black") +
theme_bw() +
ylab(expression(paste("CO"[2],"-eq flux (log"[10], " g m"^2," day"^-1," + 1.8)")))+
xlab("pH") +
theme(axis.text = element_text(size = 12),
			axis.title = element_text(size = 14),
			strip.text.x = element_text(size = 14),
		  panel.grid.major = element_blank(),
		  panel.grid.minor = element_blank()) +
annotate("text", x = 6.8, y = -0.7, label = "r = - 0.60", colour = " dark grey", size = 6) +
annotate("text", x = 6.82, y = -0.85, label = "p < 0.001", colour = " dark grey", size = 6)



pH_CO2eq2 =
ggplot(data = MeanData_pH, aes(y = Flux_mgm2day_CO2eq_D_tr, x = pH_raw_1)) +
geom_point() +
geom_smooth(method = "loess", se = T, col = "black") +
theme_bw() +
ylab(expression(paste("CO"[2],"-eq flux (log"[10], " g m"^2," day"^-1," + 1.8)")))+
xlab("pH") +
theme(axis.text = element_text(size = 12),
			axis.title = element_text(size = 14),
			strip.text.x = element_text(size = 14),
		  panel.grid.major = element_blank(),
		  panel.grid.minor = element_blank())


ggsave(pH_CO2eq, file = "../Results/pH_CO2eq.pdf", width = 6.146341, height  = 6.609756)
ggsave(pH_CO2eq, file = "../Results/pH_CO2eq.png", width = 6.146341, height  = 6.609756)







# Plot CH4 vs pH
cor.test(log10(MeanData_pH$Flux_mgm2day_CH4_D_tr+2)-3, MeanData_pH$pH_raw_1)

pH_CH4 =
ggplot(data = MeanData_pH, aes(y = log10(Flux_mgm2day_CH4_D_tr+2)-3, x = pH_raw_1)) +
geom_point() +
theme_bw() +
ylab(expression(paste("CH"[4]," flux (log"[10], " g m"^2," day"^-1," + 0.002)")))+
xlab("pH") +
theme(axis.text = element_text(size = 12),
			axis.title = element_text(size = 14),
			strip.text.x = element_text(size = 14),
		  panel.grid.major = element_blank(),
		  panel.grid.minor = element_blank())

ggsave(pH_CH4, file = "../Results/pH_CH4.pdf", width = 6.146341, height  = 6.609756)
ggsave(pH_CH4, file = "../Results/pH_CH4.png", width = 6.146341, height  = 6.609756)




# All three plots together
All_pH_plots = plot_grid(pH_Oxygen,pH_CO2,pH_CH4, ncol = 1, align = "hv", labels = c("A", "B", "C"))
save_plot(filename = "../Results/All_pH_plots.pdf", plot = All_pH_plots, base_width = 5, base_height = 11)




# Plot Oxygen vs CO2
cor.test(MeanData_pH$Oxygen_mg_raw_1, log10(MeanData_pH$Flux_mgm2day_CO2_D_tr+1800)-3)

CO2_Oxygen =
ggplot(data = MeanData_pH, aes(y = Oxygen_mg_raw_1, x = log10(Flux_mgm2day_CO2_D_tr+1800)-3)) +
geom_point() +
geom_smooth(method = "lm", se = F, col = "black") +
theme_bw() +
ylab(expression(paste("Oxygen (mg L"^-1, ")")))+
xlab(expression(paste("CO"[2]," flux (log"[10], " g m"^2," day"^-1," + 1.8)")))+
theme(axis.text = element_text(size = 12),
			axis.title = element_text(size = 14),
			strip.text.x = element_text(size = 14),
		  panel.grid.major = element_blank(),
		  panel.grid.minor = element_blank()) +
annotate("text", x = 1, y = 18, label = "r = -0.82", colour = " dark grey", size = 6) +
annotate("text", x = 1, y = 17, label = "p < 0.001", colour = " dark grey", size = 6)


ggsave(CO2_Oxygen, file = "../Results/CO2_Oxygen.pdf", width = 6.146341, height  = 6.609756)
ggsave(CO2_Oxygen, file = "../Results/CO2_Oxygen.png", width = 6.146341, height  = 6.609756)



# Another version
All_pH_plots2 = plot_grid(pH_Oxygen,CO2_Oxygen, pH_CO2, ncol = 2, align = "hv", labels = c("A", "B", "C"))
save_plot(filename = "../Results/All_pH_plots2.pdf", plot = All_pH_plots2, base_width = 10.231707, base_height = 9.646341)







# Check for effect of farm dam treatment
# No effect of management on pH
Management_pH.lm1 = lme(pH_raw_1 ~ Dam_type2,
				random = ~1|Owner,
				method = "ML",
				weights = varIdent(form = ~1|Dam_type2),
				data = MeanData_pH,
				na.action = na.omit)
Anova(Management_pH.lm1)






# Pair plot
PairPlot = with(MeanData_pH, 
	data.frame(CO2eq = log10(Flux_mgm2day_CO2eq_D_tr+1800)-3,
						 CO2 = log10(Flux_mgm2day_CO2_D_tr+1800)-3,
						 CH4 = log10(Flux_mgm2day_CH4_D_tr + 2) -3,
						 pH = pH_raw_1,
						 Oxygen = Oxygen_mg_raw_1
	)) 
chart.Correlation(PairPlot, histogram=TRUE, method = c("pearson"))









# Correction of the p-value
# ------------------------- 
# Use the function "p.adjust"

pVals = c(
	0.021,
	0.2,
	0.21,
	0.42,
	0.007,
	0.010,
	0.018,
	0.518,
	0.001,
	0.001,
	0.001,
	0.001,
	0.012,
	0.001,
	0.472,
	0.222,
	0.617,
	0.689,
	0.001,
	0.001,
	0.001,
	0.692,
	0.253
	)

p.adjust(pVals, "fdr")





# ===================
# Relationship between oxygen and phosphoros and nitrogen
# ===================

Oxygen_Nutrients.lm1 = lme(pH_raw_1 ~ Dam_type2,
				random = ~1|Owner,
				method = "ML",
				weights = varIdent(form = ~1|Dam_type2),
				data = MeanData_pH,
				na.action = na.omit)


plot(MeanData_pH$TotalP_raw_1,MeanData_pH$Oxygen_mg_raw_1)
cor.test(MeanData_pH$TotalN_raw_1,MeanData_pH$Oxygen_mg_raw_1)











# Create a summary table with all statistics
# ------------------------- 

# Add more variables (e.g., pH, conductivity, )

# Keep only numeric vars
MeanData_CO2eq_D_num = MeanData_CO2_D %>%
	#dplyr::select_if(is.numeric) %>%
	dplyr::select(SiteCode,
								Oxygen_mg_raw_1,
								TotalN_raw_1,
								TotalP_raw_1,
								Conductivity_raw_1,
								WaterTemp_raw_1,
								Stock_tCha_raw_1,
								Area_M2_gmap_raw_1,
								Lon_1,
								Lat_1,
								Flux_mgm2day_tr_1) %>%
	dplyr::rename(Flux_mgm2day_CO2_D_tr = Flux_mgm2day_tr_1) %>%
	merge(MeanData_CO2eq_D[,c("SiteCode", "Flux_mgm2day_CO2eq_D_tr", "Flux_mgm2day_CH4_D_tr")], all = T) %>%
	merge(MeanData_pH[,c("SiteCode", "pH_raw_1")], all = T) %>%
		dplyr::mutate(
								Oxygen = Oxygen_mg_raw_1,
								WaterTemp = WaterTemp_raw_1,
								Stock = 10^Stock_tCha_raw_1,
								Area = 10^Area_M2_gmap_raw_1,
								Lon = Lon_1,
								Lat = Lat_1,
								CH4 = Flux_mgm2day_CH4_D_tr/1000,
								CO2 = Flux_mgm2day_CO2_D_tr/1000,
								CO2_eq = Flux_mgm2day_CO2eq_D_tr/1000,
								TotalN = 10^TotalN_raw_1,
								TotalP = 10^TotalP_raw_1,
								pH = pH_raw_1,
								Conductivity = 10^Conductivity_raw_1
								) %>%
	dplyr::select(
								Area,
								Lon,
								Lat,
								TotalN,
								TotalP,
								Oxygen,
								WaterTemp,
								pH,
								Conductivity,
								CH4,
								CO2,
								CO2_eq,
								Stock)


SummaryTable = NULL

for(i in 1:dim(MeanData_CO2eq_D_num)[2]) {
#i = 1
		loop.var = na.omit(MeanData_CO2eq_D_num[,i])

		loop = data.frame(
			Name = names(MeanData_CO2eq_D_num)[i], 
			Rep = length(loop.var), 
			Min = min(loop.var),
			Mean = mean(loop.var),
			Mean.log10 = mean(na.omit(log10(loop.var))),
			Median = median(loop.var),
			Max = max(loop.var),
			Var = var(loop.var),
			Var.log10 = var(na.omit(log10(loop.var)))
						)

	SummaryTable = rbind(SummaryTable, loop)

}



write.csv(SummaryTable, file = "../Results/SummaryTable.csv")








# ===================
# Lit review
# ===================

LitReview = read_xlsx("/Volumes/GoogleDrive/My Drive/Fellowship BCL/Modelling/14.Farm dam fieldwork with ANU/Analysis/12th try_raw values and new narrative_updated fd classification/Lit review/Lit Review Data 25Mar2022_MM.xlsx", sheet = "Dataframe2")

# Parse the models
LitReview$Gas_parsed = factor(LitReview$Gas, 
	levels = c("CH4","CO2"),
	labels = c(
				expression(paste("CH"[4]," flux (g m"^2," day"^-1,")")),
				expression(paste("CO"[2]," flux (g m"^2," day"^-1,")"))
				)
	)


# Parse the X axis
LitReview$Author_Country = paste0(LitReview$Country,"\n","(",LitReview$Author,")")



LitReview_plot =
	ggplot(data = LitReview, aes(x = Author_Country, col = LU)) + 
	#geom_pointrange(aes(ymin = Min/1000, ymax = Max/1000, size = N), position=position_dodge(width=0.2), fatten = 2) +	
	geom_linerange(aes(ymin = Min/1000, ymax = Max/1000), position=position_dodge(width=0.2), size = 1.5, show.legend = F) +
	geom_point(aes(y = Mean/1000, size = N), position=position_dodge(width=0.2)) +
	facet_grid(Gas_parsed~., scale = "free_y", labeller = labeller(Gas_parsed = label_parsed), switch="y") +
	labs(x = "",
       y = "",
       col = "") +
	theme_bw() +
	theme(axis.text = element_text(size = 12),
			axis.title = element_text(size = 14),
		  panel.grid.major = element_blank(),
		  panel.grid.minor = element_blank(),
			strip.text.x = element_text(size = 14),
			strip.text.y = element_text(size = 14),
			legend.position = c(0.9,0.8),
			legend.text = element_text(size = 14),
			legend.key.size = unit(1, 'cm')) +
	scale_y_continuous(position = "right") +
	scale_size(range = c(2, 10), trans="log10") +
	guides(col = guide_legend(override.aes = list(size = 5)))


ggsave(LitReview_plot, file = "../Results/LitReview_plot.pdf", width = 13, height  = 7.5)




# Try log10 axis
# --------------
LitReview_log10 = LitReview
LitReview_log10_CH4 = subset(LitReview_log10, Gas == "CH4")
LitReview_log10_CH4$Min = ifelse(LitReview_log10_CH4$Min < 0.1, 0.1, LitReview_log10_CH4$Min)

LitReview_log10_CO2 = subset(LitReview_log10, Gas == "CO2")
LitReview_log10_CO2$Mean = LitReview_log10_CO2$Mean + 2300
LitReview_log10_CO2$Min = LitReview_log10_CO2$Min + 2300
LitReview_log10_CO2$Max = LitReview_log10_CO2$Max + 2300

LitReview_log10_all = rbind(LitReview_log10_CH4, LitReview_log10_CO2)

# Parse the models
LitReview_log10_all$Gas_parsed = factor(LitReview$Gas, 
	levels = c("CH4","CO2"),
	labels = c(
				expression(paste("CH"[4]," flux (g m"^2," day"^-1,")")),
				expression(paste("CO"[2]," flux (g m"^2," day"^-1," + 2.3)"))
				)
	)


Overall.range.CH4 = round(with(LitReview_log10_CH4,
					c(Mean = mean(Mean), 
					  Min = mean(Min),
					  Max = mean(Max)
					  )),1
				)/1000

Overall.range.CO2 = round(with(na.omit(LitReview_log10_CO2),
					c(Mean = mean(Mean), 
					  Min = mean(Min),
					  Max = mean(Max)
					  )),1
				)/1000



LitReview_plot_log10 =
	ggplot(data = LitReview_log10_all, aes(x = Author_Country, col = LU)) + 
	#geom_pointrange(aes(ymin = Min/1000, ymax = Max/1000, size = N), position=position_dodge(width=0.2), fatten = 2) +	
	geom_linerange(aes(ymin = Min/1000, ymax = Max/1000), position=position_dodge(width=0.2), size = 1.5, show.legend = F) +
	geom_point(aes(y = Mean/1000), position=position_dodge(width=0.2), size = 5) +
	facet_grid(Gas_parsed~., scale = "free_y", labeller = labeller(Gas_parsed = label_parsed), switch="y") +
	labs(x = "",
       y = "",
       col = "") +
	theme_bw() +
	theme(axis.text = element_text(size = 12),
			axis.title = element_text(size = 14),
		  panel.grid.major = element_blank(),
		  panel.grid.minor = element_blank(),
			strip.text.x = element_text(size = 14),
			strip.text.y = element_text(size = 14),
			legend.position = c(0.9,0.2),
			legend.text = element_text(size = 14),
			legend.key.size = unit(1, 'cm')) +
	scale_y_continuous(position = "right") +
	#scale_size(range = c(2, 10), trans="log10") +
	guides(col = guide_legend(override.aes = list(size = 5))) +
	scale_y_log10(#expand = expansion(#mult = c(0, .1)), #limits = c(0.05, 1000),
  						  breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "l") +
  annotate("text", x = 1.179, y = 1.29,
       label = paste0("CH4 = ", Overall.range.CH4["Mean"], " (", Overall.range.CH4["Min"], "; ", Overall.range.CH4["Max"], ")")
       ) +
  annotate("text", x = 1, y = 1.29,
       label = paste0("CO2 = ", Overall.range.CO2["Mean"], " (", Overall.range.CO2["Min"], "; ", Overall.range.CO2["Max"], ")")
       )

ggsave(LitReview_plot_log10, file = "../Results/LitReview_plot_log10.pdf", width = 13, height  = 7.5)




