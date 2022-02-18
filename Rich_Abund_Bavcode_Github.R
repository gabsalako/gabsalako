library(raster) #Work with raster data
library(rgdal) #Export GeoTIFFs and other core GIS functions
library (randomForest)
library(rpart)
library(sp)	
library(caret)#for confusion matrix
library(corrplot)	
library(cluster)
library(caTools)#for data split
library(biomod2)
library(usdm)
library(dismo)#for GBM dismo
library(mapview) # for view data on openstreet and other mapping platforms
library(virtualspecies)	
library(data.table)
library(ggpubr)
library(kriging)
library(gbm)#for variable contribution, gbm plot
library(gstat)
library(ecospat) # to plot corr matrix with graph
library(mgcv) #to run GAM
library(moments)
library(dplyr)
library(GGally)# for ggpairs
library(rms) # for C AUC
library(ggplot2)
library(factoextra)
library(FactoMineR)
library(pROC)#for evaluating model e.g AUC(ROC)'
library(SSDM)
library(performance)#comparing model performance
library(dominanceanalysis)
setwd("~/GBM_code_Abundance_2")
list.files()
#environmental variables
Alt <- raster("Alt_res.asc")
SoilMoisture <- raster("ResMsk_SM.asc")
LULC <- raster("re_LULC2_output.asc")
Bio1 <- raster("Bio1_output.asc")
Bio12 <- raster( "Bio12_output.asc")
Bio13 <- raster("Bio13_output.asc")
Bio14 <- raster("Bio14_output.asc")
GDDO <- raster("GDDO_output.asc")
NGDO <- raster("NGDO_output.asc")
Soil_Comp <- raster("Soil_Compgrfres.asc")
Soil_M <- raster("ResMsk_SM.asc")
SOM_npp <- raster("Npp_output.asc")
Soil_AirC <- raster( "Soil_Air.asc")
Soil_Compact <- raster("Soil_Compact.asc")
soil_Depth <- raster("Soil_Depth.asc")
Soil_AirC2 <- raster("Soil_AirC.asc")
Soil_sand <- raster("Soil_sand.tif")
Soil_clay <- raster("Soil_Clay.tif")
Soil_silt <- raster( "Soil_silt_Clip1.tif")
SoilMatter <- raster("SOM.asc")
Habitat <- raster("EUNIS_Hab.asc")
habitat_F <- as.factor(Habitat)
Habitat_factor <- raster("EUNIS_Factor.asc")
HaB_N <- raster("EUNIS_Fac_N.asc")
Mod_Sand <- raster("Mod_Sand.asc")
Mod_Silt <- raster("Mod_Silt.asc")
Mod_Clay <- raster("Mod_Clay.asc")
Abund_T <- raster("Abundance_RF_som.asc")
Pred_Soil_Ph <- raster("New_Soil_Ph.asc")
ModSoilPh <- raster("Mod_Ph.asc")
Soil_cla_sil_stack=stack(Soil_clay, Soil_sand, Soil_silt)
Pro_Soil_Text <- projectRaster(Soil_cla_sil_sa, crs="+proj=longlat +datum=WGS84")
res_Soil_Text <- resample(Pro_Soil_Text, Bio12)
res_Soil_Text
setwd("~/GBM_code_Abundance_2")
list.files()

#Abund_BAVaria
#predictions of abundance and richness
Abund_BavRfT_pred <- raster("Abundance_BavT_Pred.asc")
Abund_BavRf_pred <- raster("Abundance_Bav_Pred.asc")#use
Abund_BavGLM <- raster("Abundance_BavGAM.asc")
Abund_Bav_GBM <- raster("Abundance_BavGBM.asc")
Rich_GBM <- raster( "Richness_GBM.asc")
Rich_RF <- raster("Richness_RF_som.asc")
Rich_GAM <- raster("Richness_GAM.asc")
Rich_GLM <- raster("Rich_GLM.asc")

proj4Str <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
Abund_Bav_csv <- read.csv("Abund_Bav.csv")#load and read excel abundance data
Abund_Bav_csv[is.na(Abund_Bav_csv[])] <- 0 #remove NA
Abund_Bav_Points <- SpatialPointsDataFrame(coords = Abund_Bav_csv[,c("Longitude","Latitude")], data = Abund_Bav_csv, proj4string = CRS(proj4Str))
plot(Abund_Bav_Points)
EnvStackFac_Som1=stack(HaB_N, SoilMatter, Bio1, Bio12, Soil_Compact, Soil_AirC2, soil_Depth, ModSoilPh, Mod_Clay, Mod_Silt, SoilMoisture)#stack env. variables
RasAbund_BavExt=extract(EnvStackFac_Som1, Abund_Bav_Points)#extract values from rasters to points
Abund_Bav_Comb=cbind(Abund_Bav_Points, RasAbund_BavExt)#combine Abundance points and extract values
DF_Abund_Bav_comb <- data.frame(Abund_Bav_Comb)#creaate a dataframe 
DF_Abund_Bav_comb[is.na(DF_Abund_Bav_comb[])] <- 0 #some algorithms do not tolerate NA, so remove posibele NA
DF_Abund_Bav_comb[,'EUNIS_Fac_N'] = as.factor(DF_Abund_Bav_comb[,'EUNIS_Fac_N'])# add EUNIS habitat as categorical data
summary(DF_Abund_Bav_comb)
head(DF_Abund_Bav_comb)
#Latitude       Longitude      Ave..Abundance_Bav EUNIS.habitat.type.classification  EUNIS_Fac_N       SOM           Bio1_output      Bio12_output     Soil_Compact  
#Min.   : 0.00   Min.   : 0.000   Min.   :   0.0     Length:760                        10     : 88   Min.   :-0.5133   Min.   :  0.00   Min.   :   0.0   Min.   :0.000  
#1st Qu.:49.09   1st Qu.: 9.121   1st Qu.:  40.0     Class :character                  3      : 85   1st Qu.: 5.0000   1st Qu.: 88.00   1st Qu.: 604.8   1st Qu.:1.559  
#Median :51.25   Median :10.484   Median :  97.7     Mode  :character                  5      : 81   Median : 6.2169   Median : 92.00   Median : 728.5   Median :1.622  
#Mean   :50.85   Mean   :10.588   Mean   : 159.5                                       9      : 80   Mean   : 6.3872   Mean   : 89.42   Mean   : 747.0   Mean   :1.626  
#3rd Qu.:52.17   3rd Qu.:11.916   3rd Qu.: 215.8                                       6      : 76   3rd Qu.: 7.5635   3rd Qu.: 95.00   3rd Qu.: 832.2   3rd Qu.:1.729  
#Max.   :55.05   Max.   :14.995   Max.   :1520.0                                       4      : 73   Max.   :13.0000   Max.   :112.00   Max.   :1923.0   Max.   :2.420  
#(Other):277                                                                      
#Soil_AirC        Soil_Depth        Mod_Ph         Mod_Clay        Mod_Silt       ResMsk_SM      Longitude.1       Latitude.1    optional      
#Min.   : 0.000   Min.   : 0.00   Min.   :0.000   Min.   : 0.00   Min.   : 0.00   Min.   : 0.00   Min.   : 0.000   Min.   : 0.00   Mode:logical  
#1st Qu.: 4.635   1st Qu.: 6.50   1st Qu.:4.673   1st Qu.:13.29   1st Qu.:35.21   1st Qu.:13.00   1st Qu.: 9.121   1st Qu.:49.09   TRUE:760      
#Median : 6.373   Median :14.42   Median :5.055   Median :15.63   Median :48.13   Median :14.61   Median :10.484   Median :51.25                 
#Mean   : 7.946   Mean   :13.32   Mean   :5.077   Mean   :18.16   Mean   :45.08   Mean   :14.86   Mean   :10.588   Mean   :50.85                 
#3rd Qu.:10.131   3rd Qu.:20.00   3rd Qu.:5.473   3rd Qu.:23.25   3rd Qu.:54.47   3rd Qu.:16.40   3rd Qu.:11.916   3rd Qu.:52.17                 
#Max.   :22.581   Max.   :20.00   Max.   :6.847   Max.   :38.66   Max.   :65.55   Max.   :24.82   Max.   :14.995   Max.   :55.05  

knitr::kable(head(DF_Abund_Bav_comb, n=10))
#Latitude| Longitude| Ave..Abundance_Bav|EUNIS.habitat.type.classification                                                         |EUNIS_Fac_N |       SOM| Bio1_output| Bio12_output| Soil_Compact| Soil_AirC| Soil_Depth|   Mod_Ph|  Mod_Clay| Mod_Silt| ResMsk_SM| Longitude.1| Latitude.1|optional |
#|--------:|---------:|------------------:|:-----------------------------------------------------------------------------------------|:-----------|---------:|-----------:|------------:|------------:|---------:|----------:|--------:|---------:|--------:|---------:|-----------:|----------:|:--------|
#  |  47.4131|   11.2566|             341.20|Grasslands and lands dominated by forbs, mosses or lichens (Code: E)                      |5           | 12.000000|          70|          898|     1.508144|  0.937500|   5.851029| 4.917789| 10.311453| 38.90531|  4.380593|     11.2566|    47.4131|TRUE     |
#  |  47.5797|    9.6192|             192.80|Shrub plantations for ornamental purposes or for fruit, other than vineyards (Code: FB.3) |2           | 12.000000|         101|         1408|     1.321762|  5.036703|  17.805536| 4.495896|  9.502277| 40.68776| 11.668816|      9.6192|    47.5797|TRUE     |
#  |  47.5832|   12.9210|              51.40|Fir and spruce woodland (Code: G3.1)                                                      |1           | 12.000000|          33|         1764|     1.806078|  5.070937|   3.461475| 4.432995| 11.669026| 42.95113| 14.635960|     12.9210|    47.5832|TRUE     |
##  |  47.5825|   12.9230|              60.53|Fir and spruce woodland (Code: G3.1)                                                      |1           | 12.000000|          33|         1764|     1.806078|  5.070937|   3.461475| 4.432995| 11.669026| 42.95113| 14.635960|     12.9230|    47.5825|TRUE     |
#  |  47.6302|   10.0521|             174.00|Grasslands and lands dominated by forbs, mosses or lichens (Code: E)                      |5           | 12.000000|          73|         1923|     1.903132|  5.039444|   9.282712| 4.763350| 10.311453| 31.45209| 10.195517|     10.0521|    47.6302|TRUE     |
#  |  47.6325|   10.8214|              95.00|Fir and spruce woodland (Code: G3.1)                                                      |3           |  7.000000|          77|         1448|     1.970000|  4.991380|   2.336837| 4.464911| 12.493094| 29.22180|  7.093822|     10.8214|    47.6325|TRUE     |
#  |  47.6391|    9.8960|             386.80|Grasslands and lands dominated by forbs, mosses or lichens (Code: E)                      |2           |  8.000000|          87|         1573|     1.805333|  6.225870|  12.877278| 5.150857|  8.913555| 30.32485| 12.736215|      9.8960|    47.6391|TRUE     |
#  |  47.6452|   11.0270|             375.60|Grasslands and lands dominated by forbs, mosses or lichens (Code: E)                      |6           |  8.368593|          73|         1362|     1.964913|  8.267462|   2.709611| 4.107277| 15.630113| 33.86195|  6.912134|     11.0270|    47.6452|TRUE     |
#  |  47.6754|    9.5886|              59.81|Grasslands and lands dominated by forbs, mosses or lichens (Code: E)                      |5           |  8.000000|          98|         1281|     1.314000|  5.006034|  19.341015| 4.742937| 10.708864| 39.36699| 12.890745|      9.5886|    47.6754|TRUE     |
#  |  47.6816|   12.9309|             263.60|Grasslands and lands dominated by forbs, mosses or lichens (Code: E)                      |6           | 11.465095|          81|         1208|     1.804304|  4.780263|   2.684650| 5.052328| 10.998776| 42.95113| 14.392983|     12.9309|    47.6816|TRUE     |

#fitting the model
#Randomforest
RF_Abund_Bav <- randomForest(Ave..Abundance_Bav ~ EUNIS_Fac_N + Mod_Ph + Soil_AirC + ResMsk_SM
                             + Bio12_output + SOM + Bio1_output + Mod_Silt + Mod_Clay + Soil_Depth + Soil_Compact, data = DF_Abund_Bav_comb)
Abund_Bav_RF2 <- randomForest(x= DF_Abund_Bav_comb[,5:15], y= DF_Abund_Bav_comb[,3], ntree = 1000,nodesize = 10, importance = T)
pred_Abund_Bav <- predict(EnvStackFac_Som1, RF_Abund_Bav)
plot(pred_Abund_Bav)
pred_abund_BavRF2 <- predict(EnvStackFac_Som1, Abund_Bav_RF2)
plot(pred_abund_BavRF2)
#GLM
GLM_Abund_Bav <- glm(Ave..Abundance_Bav ~ EUNIS_Fac_N + Mod_Ph + Soil_AirC + ResMsk_SM + Bio12_output 
                     + SOM + Bio1_output + Mod_Silt + Mod_Clay + Soil_Depth, data = DF_Abund_Bav_comb)#Use GLM
pred_GLM_AbundBav <- predict(EnvStackFac_Som1, GLM_Abund_Bav)
plot(pred_GLM_AbundBav)
summary(GAM_Abund_Bav)
GLMgau_Abund_Bav <- glm(Ave..Abundance_Bav ~ EUNIS_Fac_N + Mod_Ph + Soil_AirC + ResMsk_SM + Bio12_output 
                        + SOM + Bio1_output + Mod_Silt + Mod_Clay + Soil_Depth,family = "gaussian", data = DF_Abund_Bav_comb)

summary(GLMgau_Abund_Bav)
Call:
  glm(formula = Ave..Abundance_Bav ~ EUNIS_Fac_N + Mod_Ph + Soil_AirC + 
        ResMsk_SM + Bio12_output + SOM + Bio1_output + Mod_Silt + 
        Mod_Clay + Soil_Depth, family = "gaussian", data = DF_Abund_Bav_comb)

#Deviance Residuals: 
#  Min       1Q   Median       3Q      Max  
#-271.75  -105.73   -48.16    51.23  1415.22  

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)   222.47136   67.69994   3.286  0.00106 ** 
#EUNIS_Fac_N1  151.79041  117.89029   1.288  0.19830    
#EUNIS_Fac_N2  115.17327  118.11316   0.975  0.32982    
#EUNIS_Fac_N3  127.64243  118.60521   1.076  0.28219    
#EUNIS_Fac_N4  105.95011  116.24412   0.911  0.36236    
#EUNIS_Fac_N5  168.53601  118.00601   1.428  0.15366    
#EUNIS_Fac_N6  144.74586  118.97723   1.217  0.22415    
#EUNIS_Fac_N7  116.38867  118.30293   0.984  0.32553    
#EUNIS_Fac_N8  123.21688  118.18025   1.043  0.29747    
#EUNIS_Fac_N9  134.24085  118.03931   1.137  0.25580    
#EUNIS_Fac_N10 114.96175  118.33260   0.972  0.33161    
#Mod_Ph         -2.43567   10.93032  -0.223  0.82372    
#Soil_AirC       1.96135    1.72714   1.136  0.25649    
#ResMsk_SM      -7.56095    1.85626  -4.073 5.14e-05 ***
#Bio12_output   -0.09633    0.04129  -2.333  0.01991 *  
#SOM             7.69364    2.58946   2.971  0.00306 ** 
#Bio1_output     0.73845    0.98486   0.750  0.45361    
#Mod_Silt       -2.21764    0.79037  -2.806  0.00515 ** 
#Mod_Clay       -0.14149    1.40106  -0.101  0.91959    
#Soil_Depth     -1.73363    1.13361  -1.529  0.12662    
#---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#(Dispersion parameter for gaussian family taken to be 31639.93)

par(mfrow = c(1, 1))
par(mfrow = c(4, 3))
partialPlot(RF_Abund_Bav, DF_Abund_Bav_comb, "Mod_Ph", plot = TRUE)#the model , the train data, the variable to plot
partialPlot(RF_Abund_Bav, DF_Abund_Bav_comb, "Mod_Clay", plot = TRUE)
partialPlot(RF_Abund_Bav, DF_Abund_Bav_comb, "SOM", plot = TRUE)
partialPlot(RF_Abund_Bav, DF_Abund_Bav_comb, "Bio1_output", plot = TRUE)
partialPlot(RF_Abund_Bav, DF_Abund_Bav_comb, "Bio12_output", plot = TRUE)
partialPlot(RF_Abund_Bav, DF_Abund_Bav_comb, "Mod_Silt", plot = TRUE)
partialPlot(RF_Abund_Bav, DF_Abund_Bav_comb, "ResMsk_SM", plot = TRUE)
partialPlot(RF_Abund_Bav, DF_Abund_Bav_comb, "Soil_AirC", plot = TRUE)
partialPlot(RF_Abund_Bav, DF_Abund_Bav_comb, "Soil_Depth", plot = TRUE)
partialPlot(RF_Abund_Bav, DF_Abund_Bav_comb, "Soil_Compact", plot = TRUE)
partialPlot(RF_Abund_Bav, DF_Abund_Bav_comb, "EUNIS_Fac_N", plot = TRUE)
varImpPlot(RF_Abund_Bav)

##Train data set Abund_Bav
sample_Abund_Bav = floor(0.7*nrow(DF_Abund_Bav_comb))#spit data into train (70) and test (30)
set.seed(123)
picked_Abund_Bav = sample(seq_len(nrow(DF_Abund_Bav_comb)),size = sample_Abund_Bav)
Development_Abund_Bav =DF_Abund_Bav_comb[picked_Abund_Bav,]
holdout =DF_Abund_Bav_comb[-picked_Abund_Bav,]
Abund_Bav_Trtest =DF_Abund_Bav_comb[picked_Abund_Bav,]
Abund_BavTest2 =DF_Abund_Bav_comb[-picked_Abund_Bav,]
Abund_BavTrtest_DF <- data.frame(Abund_Bav_Trtest)
Abund_Bav_Test2_DF <- data.frame(Abund_BavTest2)
Abund_BavTrtest_DF[is.na(Abund_BavTrtest_DF[])] <- 0 #very essential for the algorithm does not tolerate NA
Abund_BavTrtest_DF[,'EUNIS_Fac_N'] = as.factor(Abund_BavTrtest_DF[,'EUNIS_Fac_N'])#categorical
Abund_Bav_Test2_DF[is.na(Abund_Bav_Test2_DF[])] <- 0
Abund_Bav_Test2_DF[,'EUNIS_Fac_N'] = as.factor(Abund_Bav_Test2_DF[,'EUNIS_Fac_N'])#categorical
summary(Abund_BavTrtest_DF)
summary(Abund_Bav_Test2_DF)
## fiting with train dataset Abund_Bav
#Randomforest
Abund_Bav_RFTr <- randomForest(x= Abund_BavTrtest_DF[,5:15], y= Abund_BavTrtest_DF[,3], ntree = 1000,nodesize = 10, importance = T)
pred_Abund_BavTr <- predict(EnvStackFac_Som1, Abund_Bav_RFTr)
Abundance_RF <- pred_Abund_BavTr
plot(pred_Abund_BavTr)
crs (pred_Abund_BavTr) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 "
pred_abund_BavRF2
#GBM
GBM_Dismo_Bav <- gbm.step(data=Abund_Bav_Trtest, gbm.x = 5:15, gbm.y = 3,
                          family = "gaussian", tree.complexity = 5,
                          learning.rate = 0.005, bag.fraction = 0.5)
gbm.plot(GBM_Dismo_Bav)
gbm.plot.fits(GBM_Dismo_Bav)#use for habitat plot
pred_GBM_Bav <- predict(EnvStackFac_Som1, GBM_Dismo_Bav)
Abundance_GBM <- pred_GBM_Bav

writeRaster(pred_Abund_BavTr, filename = "Abundance_BavT_Pred.asc", format='ascii', overwrite=TRUE)#convert to ASCII and save line 1153
writeRaster(pred_Abund_Bav, filename = "Abundance_Bav_Pred.asc", format='ascii', overwrite=TRUE)#convert to ASCII and save line 1110
writeRaster(pred_GLM_AbundBav, filename = "Abundance_BavGAM.asc", format='ascii', overwrite=TRUE)#convert to ASCII and save line 1116
writeRaster(pred_GBM_Bav, filename = "Abundance_BavGBM.asc", format='ascii', overwrite=TRUE)#convert to ASCII and save line 1116

#rename and reload
Abund_BavRfT_pred <- raster("Abundance_BavT_Pred.asc")
Abund_BavRf_pred <- raster("Abundance_Bav_Pred.asc")#use best R2
Abund_BavGLM <- raster("Abundance_BavGAM.asc")
Abund_Bav_GBM <- raster("Abundance_BavGBM.asc")

#prediction dataframe for evaluation
PredAbunds_Bavtack=stack(Abund_BavRfT_pred, Abund_BavRf_pred, Abund_BavGLM, Abund_Bav_GBM)#stack the predicted raster saved into ascii and load back (see line)
RasPred_Abund_Bav_ext=extract(PredAbunds_Bavtack, Abund_Bav_Points)#extract point to raster
Pred_Abund_BavComb=cbind(Abund_Bav_Points, RasPred_Abund_Bav_ext)
DF_Prd_Abund_Bav <- data.frame(Pred_Abund_BavComb)
DF_Prd_Abund_Bav[is.na(DF_Prd_Abund_Bav[])] <- 0
summary(DF_Prd_Abund_Bav)
#Latitude       Longitude      Ave..Abundance_Bav EUNIS.habitat.type.classification Abundance_BavT_Pred Abundance_Bav_Pred Abundance_BavGAM Abundance_BavGBM
#Min.   : 0.00   Min.   : 0.000   Min.   :   0.0     Length:760                        Min.   :  0.00      Min.   :  0.00     Min.   :  0.0    Min.   :  0.0   
#1st Qu.:49.09   1st Qu.: 9.121   1st Qu.:  40.0     Class :character                  1st Qu.: 97.43      1st Qu.: 80.78     1st Qu.:119.6    1st Qu.:154.2   
#Median :51.25   Median :10.484   Median :  97.7     Mode  :character                  Median :147.62      Median :127.40     Median :154.3    Median :158.6   
#Mean   :50.85   Mean   :10.588   Mean   : 159.5                                       Mean   :164.65      Mean   :158.05     Mean   :155.6    Mean   :160.3   
#3rd Qu.:52.17   3rd Qu.:11.916   3rd Qu.: 215.8                                       3rd Qu.:211.09      3rd Qu.:201.06     3rd Qu.:193.3    3rd Qu.:167.4   
#Max.   :55.05   Max.   :14.995   Max.   :1520.0                                       Max.   :730.99      Max.   :839.42     Max.   :347.4    Max.   :213.8   
#Longitude.1       Latitude.1    optional      
#Min.   : 0.000   Min.   : 0.00   Mode:logical  
#1st Qu.: 9.121   1st Qu.:49.09   TRUE:760      
#Median :10.484   Median :51.25                 
#Mean   :10.588   Mean   :50.85                 
#3rd Qu.:11.916   3rd Qu.:52.17                 
#Max.   :14.995   Max.   :55.05                 

##AUC and R2 from lrm
AUCLRM_Abund_BavRFT <- lrm(Ave..Abundance_Bav  ~ Abundance_BavT_Pred  , data = DF_Prd_Abund_Bav, x= TRUE, y = TRUE)#R2 between the observed and predicted values in lrm
AUCLRM_Abund_BavRFT
#R2= 0.577 AUC = 0.768
AUCLRM_Abund_BavRF <- lrm(Ave..Abundance_Bav  ~ Abundance_Bav_Pred, data = DF_Prd_Abund_Bav, x= TRUE, y = TRUE)#R2 between the observed and predicted values in lrm=0.840
#R2= 0.840 AUC = 0.861
AUCLRM_Abund_BavGLM <- lrm(Ave..Abundance_Bav  ~ Abundance_BavGAM, data = DF_Prd_Abund_Bav, x= TRUE, y = TRUE)#R2 between the observed and predicted values in lrm
#R2= 0.089 AUC = 0.603
AUCLRM_Abund_BavGBM <- lrm(Ave..Abundance_Bav  ~ Abundance_BavGBM, data = DF_Prd_Abund_Bav, x= TRUE, y = TRUE)#R2 between the observed and predicted values in lrm
#R2= 0.050, AUC = 0.647

AUCLRM_Abund_BavRF <- lrm(Ave..Abundance_Bav  ~ Abundance_Bav, data = DF_Prd_Abund_Bav, x= TRUE, y = TRUE)#R2 between the observed and predicted values in lrm, high R2

#Model Likelihood    Discrimination    Rank Discrim.    
#Ratio Test           Indexes          Indexes    
#Obs           760    LR chi2    1395.06    R2       0.840    C       0.961*    
#max |deriv| 4e-08    d.f.             1    g        5.326    Dxy     0.722    
#Pr(> chi2) <0.0001    gr     205.661    gamma   0.727    
#gp       0.447    tau-a   0.721    
#Brier    0.071                     

#plot the predicted and observed graph with ggplot (Dataframe, the prediction on x axis and observed on y axis)
ggplot(DF_Prd_Abund_Bav,aes(x = Abundance_BavT_Pred, y = Ave..Abundance_Bav)) + geom_point() +geom_smooth(method = "lm")
ggplot(DF_Prd_Abund_Bav,aes(x = Abundance_Bav_Pred , y = Ave..Abundance_Bav)) + geom_point() +geom_smooth(method = "lm")#very good r2 and graph
ggplot(DF_Prd_Abund_Bav,aes(x = Abundance_BavGAM, y = Ave..Abundance_Bav)) + geom_point() +geom_smooth(method = "lm")
ggplot(DF_Prd_Abund_Bav,aes(x = Abundance_BavGBM, y = Ave..Abundance_Bav)) + geom_point() +geom_smooth(method = "lm")

#AUC overfit with this
roc.list_Bav <- roc(Ave..Abundance_Bav  ~ Abundance_BavT_Pred + Abundance_Bav_Pred + Abundance_BavGAM + Abundance_BavGBM, data = DF_Prd_Abund_Bav) #list the fitted prediction of all the models
roc.list_Bav %>% 
  map(~tibble(AUC = .x$auc)) %>% 
  bind_rows(.id = "name") -> data.auc

data.auc %>% 
  mutate(label_long=paste0(name," , AUC = ",paste(round(AUC,2))),
         label_AUC=paste0("AUC = ",paste(round(AUC,2)))) -> data.labels
ggroc(roc.list_Bav) +
  scale_color_discrete(labels=data.labels$label_long)


AUC_GmBAv <- gam(Ave..Abundance_Bav  ~ Abundance_BavT_Pred , data = DF_Prd_Abund_Bav)
AUC_Bav <- roc(DF_Prd_Abund_Bav$Ave..Abundance_Bav ~ AUC_GmBAv$fitted)
plot(AUC_Bav, legacy.axes = TRUE)

#for levelplot
#stack the new bavarian predictions
PredAbunds_Bavtack=stack(Abund_BavRfT_pred, Abund_BavRf_pred, Abund_BavGLM, Abund_Bav_GBM)#remove one of Rf (Rft)
levelplot(PredAbunds_Bavtack)#defaullt
levelplot(PredAbunds_Bavtack, at= unique(c(seq(600, 150, length=600), seq(0,150,length=600))), col.regions = colorRampPalette(c("white", "green", "yellow","red"))(1e4))
levelplot(PredAbunds_Bavtack, at= unique(c(seq(600, 100, length=600), seq(0,100,length=600))), col.regions = colorRampPalette(c("white", "white", "yellow","green","red"))(1e5))
levelplot(PredAbunds_Bavtack, at= unique(c(seq(600, 150, length=600), seq(0,150,length=600))), col.regions = colorRampPalette(c("white", "yellow", "red","green"))(1e4))

plot(PredAbunds_Bavtack)

GLMp_Abund_Bav <- glm(Ave..Abundance_Bav ~ EUNIS_Fac_N + Mod_Ph + Soil_AirC + ResMsk_SM + Bio12_output 
                      + SOM + Bio1_output + Mod_Silt + Mod_Clay + Soil_Depth,family = "poisson", data = DF_Abund_Bav_comb)
summary(GLMp_Abund_Bav)

GLMgau_Abund_Bav <- glm(Ave..Abundance_Bav ~ EUNIS_Fac_N + Mod_Ph + Soil_AirC + ResMsk_SM + Bio12_output 
                        + SOM + Bio1_output + Mod_Silt + Mod_Clay + Soil_Depth,family = "gaussian", data = DF_Abund_Bav_comb)

summary(GLMgau_Abund_Bav)

#Richness
######################################################################################################################################################################
#Richness with new som, full datasets and factor variable
Rich_csv <- read.csv("Gab_David_Rich.csv")
Av_Rich_Points <- SpatialPointsDataFrame(coords = Rich_csv[,c("Longitude","Latitude")], data = Rich_csv, proj4string = CRS(proj4Str))
head(Av_Rich_Points, 3)
RasModFac_Richsom_Ext=extract(EnvStackFac_Som1, Av_Rich_Points)
RichMod_CombFacsom=cbind(Av_Rich_Points, RasModFac_Richsom_Ext)
DFMod_Richsom_combFac <- data.frame(RichMod_CombFacsom)
DFMod_Richsom_combFac[is.na(DFMod_Richsom_combFac[])] <- 0 #some algorithms do not tolerate NA
DFMod_Richsom_combFac[,'EUNIS_Fac_N'] = as.factor(DFMod_Richsom_combFac[,'EUNIS_Fac_N'])#categorical
summary(DFMod_Richsom_combFac)
GBM_DISMO_RichFacsom <- gbm.step(data=DFMod_Richsom_combFac, gbm.x = 5:15, gbm.y = 3,
                                 family = "gaussian", tree.complexity = 5,
                                 learning.rate = 0.001, bag.fraction = 0.5) #trees adding;look good
gbm.plot(GBM_DISMO_RichFacsom, write.title = TRUE)#group plot
gbm.plot.fits(GBM_DISMO_RichFacsom)

GLM_Rich <- glm(Ave..Spp..Richness ~ EUNIS_Fac_N + Mod_Ph + Soil_AirC + ResMsk_SM + Bio12_output + SOM
                + Bio1_output + Mod_Silt + Mod_Clay + Soil_Depth, data = DFMod_Richsom_combFac)
summary(GLM_Rich)#GLM equation for richness extrapolation, full data som
Family: gaussian 
#Link function: identity 

Formula:
  Ave..Spp..Richness ~ EUNIS_Fac_N + Mod_Ph + Soil_AirC + ResMsk_SM + 
  Bio12_output + SOM + Bio1_output + Mod_Silt + Mod_Clay + 
  Soil_Depth

#Parametric coefficients:
#                  Estimate Std. Error t value Pr(>|t|)    
#  (Intercept)    3.673132   0.727102   5.052 5.24e-07 ***
#  EUNIS_Fac_MAR  -2.826857   1.267370  -2.230  0.02595 *  
#  EUNIS_Fac_COA -2.849673   1.278944  -2.228  0.02610 *  
#  EUNIS_Fac_ILW -3.149067   1.283581  -2.453  0.01433 *  
#  EUNIS_Fac_MBF -3.311925   1.265725  -2.617  0.00902 ** 
#  EUNIS_Fac_GRS -2.979800   1.276123  -2.335  0.01975 *  
#  EUNIS_Fac_HST -3.337402   1.281658  -2.604  0.00936 ** 
#  EUNIS_Fac_FOR -3.441007   1.280515  -2.687  0.00733 ** 
#  EUNIS_Fac_SpVeg -3.130856   1.288856  -2.429  0.01532 *  
#  EUNIS_Fac_CRO-3.071311   1.277282  -2.405  0.01638 *  
#  EUNIS_Fac_URB -2.968900   1.282145  -2.316  0.02079 *  
#  Mod_Ph         0.236466   0.109802   2.154  0.03152 *  
#  Soil_AirC     -0.045805   0.019664  -2.329  0.02005 *  
#  ResMsk_SM      0.008268   0.019568   0.423  0.67272    
#  Bio12_output   0.001505   0.000457   3.293  0.00103 ** 
#  SOM            0.054683   0.028090   1.947  0.05186 .  
#  Bio1_output    0.010449   0.010023   1.042  0.29745    
#  Mod_Silt      -0.006174   0.009237  -0.668  0.50406    
#  Mod_Clay       0.035451   0.015757   2.250  0.02469 *  
#  Soil_Depth     0.005720   0.012682   0.451  0.65207    
#---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
pred_Rich_GLM <- predict(EnvStackFac_Som1, GLM_Rich)
writeRaster(pred_Rich_GLM, filename = "Rich_GLM.asc", format='ascii', overwrite=TRUE)#convert to ASCII and save

#Richness with new som, train datasets and factor variable
##########################################################################################################################################################################
sample_Rich_ModfacTRsom = floor(0.7*nrow(DFMod_Richsom_combFac))#spit data into train (70) and test (30)
set.seed(123)
picked_RichModTRFsom = sample(seq_len(nrow(DFMod_Richsom_combFac)),size = sample_Rich_ModfacTRsom)
Development_RichModTRFsom =DFMod_Richsom_combFac[picked_RichModTRFsom,]
holdout =DFMod_Richsom_combFac[-picked_RichModTRFsom,]
Rich_Mod_TrainFacsom =DFMod_Richsom_combFac[picked_RichModTRFsom,]
Rich_Test_ModFacsom =DFMod_Richsom_combFac[-picked_RichModTRFsom,]
Rich_Modfac_TrainDFsom <- data.frame(Rich_Mod_TrainFacsom)
Rich_Test_ModFacsom <- data.frame(Rich_Test_ModFacsom)
Rich_Modfac_TrainDFsom[is.na(Rich_Modfac_TrainDFsom[])] <- 0 #very essential for some algorithm does not tolerate NA
Rich_Modfac_TrainDFsom[,'EUNIS_Fac_N'] = as.factor(Rich_Modfac_TrainDFsom[,'EUNIS_Fac_N'])#categorical
summary(Rich_Modfac_TrainDFsom)
#GBM
GBMfac_DISMO_RichTRsom <- gbm.step(data=Rich_Modfac_TrainDFsom, gbm.x = 5:15, gbm.y = 3,
                                   family = "gaussian", tree.complexity = 5,
                                   learning.rate = 0.001, bag.fraction = 0.5) #trees adding, the selected GBM
gbm.plot(GBMfac_DISMO_RichTRsom)
gbm.plot.fits(GBMfac_DISMO_RichTRsom)
summary(GBMfac_DISMO_RichTRsom)
#var   rel.inf
#Mod_Clay         Mod_Clay 18.964067
#EUNIS_Fac_N   EUNIS_Fac_N 17.149826
#Bio12_output Bio12_output 11.386023
#Bio1_output   Bio1_output 10.654035
#Mod_Ph             Mod_Ph  8.738080
#ResMsk_SM       ResMsk_SM  8.305669
#Soil_AirC       Soil_AirC  6.352979
#Soil_Compact Soil_Compact  6.144690
#Mod_Silt         Mod_Silt  5.013346
#SOM                   SOM  4.802054
#Soil_Depth     Soil_Depth  2.489230

predGBMRichT <- predict(EnvStackFac_Som1, GBMfac_DISMO_RichTRsom)
plot(predGBMRichT)
#RandomForest
RF_RichFTR <- randomForest(Ave..Spp..Richness ~ EUNIS_Fac_N + Mod_Ph + Soil_AirC + ResMsk_SM + Bio12_output + SOM
                           + Bio1_output + Mod_Silt + Mod_Clay + Soil_Depth + Soil_Compact, data = Rich_Modfac_TrainDFsom)
PredRichRF <- predict(EnvStackFac_Som1, RF_RichFTR)
plot(PredRichRF)
varImpPlot(RF_RichFTR)
par(mfrow = c(1, 1))
par(mfrow = c(4, 3))
#the model , the train data, the variable to plot
partialPlot(RF_RichFTR,  Rich_Modfac_TrainDFsom, "Mod_Ph", plot = TRUE)
partialPlot(RF_RichFTR,  Rich_Modfac_TrainDFsom, "Mod_Clay", plot = TRUE)
partialPlot(RF_RichFTR,  Rich_Modfac_TrainDFsom, "SOM", plot = TRUE)
partialPlot(RF_RichFTR,  Rich_Modfac_TrainDFsom,"Bio1_output", plot = TRUE)
partialPlot(RF_RichFTR,  Rich_Modfac_TrainDFsom, "Bio12_output", plot = TRUE)
partialPlot(RF_RichFTR,  Rich_Modfac_TrainDFsom, "Mod_Silt", plot = TRUE)
partialPlot(RF_RichFTR,  Rich_Modfac_TrainDFsom, "ResMsk_SM", plot = TRUE)
partialPlot(RF_RichFTR,  Rich_Modfac_TrainDFsom, "Soil_AirC", plot = TRUE)
partialPlot(RF_RichFTR,  Rich_Modfac_TrainDFsom, "Soil_Depth", plot = TRUE)
partialPlot(RF_RichFTR,  Rich_Modfac_TrainDFsom, "Soil_Compact", plot = TRUE)
partialPlot(RF_RichFTR,  Rich_Modfac_TrainDFsom,"EUNIS_Fac_N", plot = TRUE)
varImpPlot(RF_RichFTR)

Rich_RF2 <- randomForest(x= Rich_Modfac_TrainDFsom[,5:15], y= Rich_Modfac_TrainDFsom[,3], ntree = 1000,nodesize = 10, importance = T)
Pred_RichRF2 <- predict(EnvStackFac_Som1, Rich_RF2)
plot(Pred_RichRF2)
#GLM
GLM_RichFTR <- glm(Ave..Spp..Richness ~ EUNIS_Fac_N + Mod_Ph + Soil_AirC + ResMsk_SM + Bio12_output + SOM  + Bio1_output + Mod_Silt + Mod_Clay + Soil_Depth, data = Rich_Modfac_TrainDFsom)
summary(GLM_RichFTR)
pred_Rich_GLM <- predict(EnvStackFac_Som1, GLM_RichFTR)
plot(pred_Rich_GLM)
#Stack richness predictions
Rich_stackPred=stack(pred_Rich_GLM, predGBMRichT, PredRichRF)
names(Rich_stackPred) <- c("Rich_GLM", "Rich_GBM", "Rich_RF")
levelplot(Rich_stackPred, at= unique(c(seq(12, 1, length=12), seq(1, 12,length=12))), col.regions = colorRampPalette(c("white", "green", "red", "yellow"))(1e4))#good
levelplot(Rich_stackPred, at= unique(c(seq(12, 0, length=12), seq(0, 15,length=15))), col.regions = colorRampPalette(c("white", "red", "green", "yellow"))(1e4))#good

#Richness predictions of model algorithms previously saved in ascii and reloaded
Rich_GBM <- raster( "Richness_GBM.asc")
Rich_RF <- raster("Richness_RF_som.asc")
Rich_GLM <- raster("Rich_GLM.asc")

#prediction dataframe for evaluation
Pred_Rich_stack=stack(Rich_GBM, Rich_GLM, Rich_RF)#stack the predicted raster saved into ascii and load back (see line)
Ras_Rich=extract(Pred_Rich_stack, Av_Rich_Points)#extract point to raster
Pred_RichComb=cbind(Av_Rich_Points, Ras_Rich)
DF_Pred_Rich<- data.frame(Pred_RichComb)
DF_Pred_Rich[is.na(DF_Pred_Rich[])] <- 0
summary(DF_Pred_Rich)

##AUC and R2 from lrm
AUCLRM_Rich_RF <- lrm(Ave..Spp..Richness ~ Richness_RF_som  , data = DF_Pred_Rich, x= TRUE, y = TRUE)#R2 between the observed and predicted values in lrm
AUCLRM_Rich_RF
#Model Likelihood    Discrimination    Rank Discrim.    
#Ratio Test           Indexes          Indexes    
#Obs           979    LR chi2     832.55    R2       0.574    C       0.819    
#max |deriv| 0.005    d.f.             1    g        2.698    Dxy     0.638    
#Pr(> chi2) <0.0001    gr      14.855    gamma   0.639    
#gp       0.357    tau-a   0.590    
#Brier    0.134                     

AUCLRM_Rich_GBM <- lrm(Ave..Spp..Richness ~ Richness_GBM, data = DF_Pred_Rich, x= TRUE, y = TRUE)#R2 between the observed and predicted values in lrm
AUCLRM_Rich_GBM
#Model Likelihood    Discrimination    Rank Discrim.    
#Ratio Test           Indexes          Indexes    
#Obs           979    LR chi2     119.47    R2       0.115    C       0.658    
#max |deriv| 6e-10    d.f.             1    g        0.668    Dxy     0.316    
#Pr(> chi2) <0.0001    gr       1.951    gamma   0.317    
#gp       0.139    tau-a   0.292    
#Brier    0.216                     

AUCLRM_Rich_GLM <- lrm(Ave..Spp..Richness ~ Rich_GLM, data = DF_Pred_Rich, x= TRUE, y = TRUE)#R2 between the observed and predicted values in lrm
AUCLRM_Rich_GLM
#Model Likelihood    Discrimination    Rank Discrim.    
#Ratio Test           Indexes          Indexes    
#Obs           979    LR chi2      28.61    R2       0.029    C       0.583    
#max |deriv| 3e-11    d.f.             1    g        0.291    Dxy     0.165    
#Pr(> chi2) <0.0001    gr       1.338    gamma   0.166    
#gp       0.066    tau-a   0.153    
#Brier    0.235          

#ggplot
ggplot(DF_Pred_Rich,aes(x =Richness_RF_som , y = Ave..Spp..Richness)) + geom_point() +geom_smooth(method = "lm")
ggplot(DF_Pred_Rich,aes(x =Richness_GBM, y = Ave..Spp..Richness)) + geom_point() +geom_smooth(method = "lm")
ggplot(DF_Pred_Rich,aes(x =Rich_GLM , y = Ave..Spp..Richness)) + geom_point() +geom_smooth(method = "lm")

Abund_HaB_stack=stack(Abund_BavRfT_pred, Abund_BavRf_pred, Abund_BavGLM, Abund_Bav_GBM, HaB_N)
RasAbund_Hab=extract(Abund_HaB_stack, Abund_Bav_Points)
AbundHab_Comb=cbind(Abund_Bav_Points, RasAbund_Hab)
DF_AbundHab <- data.frame(AbundHab_Comb)
DF_AbundHab[is.na(DF_AbundHab[])] <- 0 
summary(DF_AbundHab)

