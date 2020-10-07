#### Lithothamnion crispatum model ----


#### setup environment ----
setwd("C:/Users/gabri/Google Drive/Mestrado/crispatum/biomod2")

## load the required packages
library(biomod2)
library(ggplot2)
library(gridExtra)
library(raster)
library(rasterVis)
library(dplyr)
library(reshape2)

#### read data ----
# Occurrence data
occ <- read.csv("envfilter_lcrispatum_pa.csv")
summary(occ)

# Environmental data
setwd("C:/Users/gabri/Google Drive/Scripts/biooracle")
envdata = Sys.glob("*.tif") #Or whatever identifies your files
stck = raster::stack() #empty stack for raster
e <- raster::extent(-80, -20, -40, 10)  #for world usually crop (-180, 180, -90, 90)
for(i in 1:NROW(envdata)){
  tempraster = raster::raster(envdata[i])
  cp_tempraster= raster::crop(tempraster, e)
  stck = raster::stack(stck, cp_tempraster)
}
names(stck)
envpred <- dropLayer(stck, c(1,3,5,7))
#envpred <- stck
names(envpred)

#------------------------------------------------------------
#### Formating the data ----
set.seed(123)
data <- 
  BIOMOD_FormatingData(
    resp.var = occ['lcrispatum'],
    resp.xy = occ[, c('lon', 'lat')],
    expl.var = envpred,
    resp.name = "lcrispatum",
    PA.nb.rep = 10,
    PA.nb.absences = 60,
    PA.strategy = 'disk',
    PA.dist.min = 200000, #2 degrees of distance
    PA.dist.max = 250000)

## Check data
data

#### Defining individual models options ---- 
# GBM according to Elith (2017), RF to Breiman (2003) and ANN default options
opt <- 
  BIOMOD_ModelingOptions(
    GBM = list(distribution = "bernoulli",
               n.trees = 10000,
               interaction.depth = 2, #tree complexity
               n.minobsinnode = 5,
               shrinkage = 0.001, #learning rate
               bag.fraction = 0.5,
               cv.folds = 10),
    RF = list(do.classif = TRUE, #if True classification random.forest computed, else regression
               #random.forest will be done
               ntree = 3500, #Number of trees to grow. If you want auxiliary information like variable importance or proximities grow
               #a lot of trees--say a 1000 or more.Run out to 5000 trees if there are many variables and I want the variables
               #importances to be stable.
               mtry = 'default', #Number of variables randomly sampled as candidates at each split
               #default = sqrt(predictors); setting mtry equal to the square root of mdim gives generally near optimum results
               nodesize = 1, #Minimum size of terminal nodes. The default that always gives good performances is 1 
               #Note that the default values are different for classification (1) and regression (5).
               maxnodes = NULL #Maximum number of terminal nodes trees in the forest can have. If not
               #given, trees are grown to the maximum possible (subject to limits by nodesize).
               ),
    ANN = list(NbCV = 10, #number of cross validation to find best size and decay parameters
               size = NULL, #number of units in the hidden layer, if NULL then size parameter
               # will be optimised by cross validation based on model AUC
               decay = NULL, #parameter for weight decay. If NULL then decay parameter will be
               #optimised by corss validation on model AUC
               rang = 0.1, #Initial random weights on [-rang,rang]
               maxit = 200 #maximum number of iterations
                 )
  )

#### Running individual models ----
models <- 
  BIOMOD_Modeling(
    data = data,
    models = c("GBM", "RF", "ANN"),
    models.options = opt,
    NbRunEval = 10, #number of repeated partitions of testing and training datasets,
    DataSplit = 70,
    VarImport = 3,
    SaveObj = T,
    Prevalence = 0.5, #absences will be weighted equally to the presences (i.e. 
    #the weighted sum of presence equals the weighted sum of absences). 
    modeling.id = "lcrispatum_bat"
    )

#### Assessing individual models quality ----

## Evaluation scores
models_scores <- get_evaluations(models)

## Variable importance
models_var_import <- get_variables_importance(models)
models_var_import

## Mean of variable importance by algorithm
apply(models_var_import, c(1,2), mean)

## Running the ensemble models ----
ensemble_models <- 
  BIOMOD_EnsembleModeling(
    modeling.output = models,
    em.by = 'all',
    eval.metric = 'TSS',
    eval.metric.quality.threshold = 0.8,
    models.eval.meth = c('TSS','ROC'),
    prob.mean = FALSE,
    prob.cv = TRUE, 
    committee.averaging = TRUE,
    prob.mean.weight = TRUE,
    prob.mean.weight.decay = 'proportional',
    VarImport = 3 #three permutations to check variable importance
  )

## Ensemble models evaluations ----
ensemble_models_scores <- get_evaluations(ensemble_models)
ensemble_models_scores

## Getting the model I want to create response plots from ----
ens <- BIOMOD_LoadModels(ensemble_models, full.name = "lcrispatum_EMcaByTSS_mergedAlgo_mergedRun_mergedData")

ens_eval_strip <- 
  biomod2::response.plot2(
    models  = ens,
    Data = get_formal_data(models,'expl.var'), 
    show.variables= get_formal_data(models,'expl.var.names'),
    do.bivariate = FALSE, #predicted response curves are plotter for every single variable independently (2 dimension)
    fixed.var.metric = 'mean', #statistic used to fix as constant the remaining variables
    legend = FALSE,
    display_title = FALSE,
    data_species = get_formal_data(models,'resp.var'),
    plot=F
  )

# Individual response plots
nit_rp <- ens_eval_strip %>%
  filter(expl.name == "BO2_nitrateltmax_bdmin_lonlat")
nit_plot <- ggplot(data=nit_rp, aes(x=expl.val, y=pred.val)) +
  geom_line(color = 'gray40', size=2) +
  labs(x="Maximum nitrate", y = "Probability of occurrence") +
  ylim(0.5, 1) +
  theme_minimal()
tiff("resp_nit.tiff", width=250, height=160, bg="white", res=500,unit="mm")
nit_plot + theme(text = element_text(size=20)) +
  scale_x_continuous(breaks = seq(0,35,5))
dev.off()

sal_rp <- ens_eval_strip %>%
  filter(expl.name == "BO2_salinityltmax_bdmin_lonlat")
sal_plot <- ggplot(data=sal_rp, aes(x=expl.val, y=pred.val)) +
  geom_line(color = 'gray40',size=2) +
  labs(x="Maximum salinity", y = "Probability of occurrence") +
  ylim(0.5, 1) +
  theme_minimal()
tiff("resp_sal.tiff", width=250, height=160, bg="white", res=500,unit="mm")
sal_plot + theme(text = element_text(size=20)) +
  scale_x_continuous(breaks = seq(30,40,1))
dev.off()

temp_rp <- ens_eval_strip %>%
  filter(expl.name == "BO2_templtmax_bdmin_lonlat")
temp_plot <- ggplot(data=temp_rp, aes(x=expl.val, y=pred.val)) +
  geom_line(color = 'gray40',size=2) +
  labs(x="Maximum temperature", y = "Probability of occurrence") +
  ylim(0.5, 1) +
  theme_minimal()
tiff("resp_temp.tiff", width=250, height=160, bg="white", res=500,unit="mm")
temp_plot + theme(text = element_text(size=20)) +
  scale_x_continuous(breaks = seq(0,31,5))
dev.off()

cal_rp <- ens_eval_strip %>%
  filter(expl.name == "calcsat_av_1951_to_2000")
cal_plot <- ggplot(data=cal_rp, aes(x=expl.val, y=pred.val)) +
  geom_line(color = 'gray40',size=2) +
  labs(x="Mean calcite saturation", y = "Probability of occurrence") +
  ylim(0.5, 1) +
  theme_minimal()
tiff("resp_cal.tiff", width=250, height=160, bg="white", res=500,unit="mm")
cal_plot + theme(text = element_text(size=20)) +
  scale_x_continuous(breaks = seq(1,6,1))
dev.off()

#### Projecting models ----

## current projections
models_proj_current <- 
  BIOMOD_Projection(
    modeling.output = models,
    new.env = envpred,
    proj.name = "current_bat",
    binary.meth = "TSS",
    output.format = ".img",
    do.stack = FALSE
  )

ensemble_models_proj_current <- 
  BIOMOD_EnsembleForecasting(
    EM.output = ensemble_models,
    projection.output = models_proj_current,
    binary.meth = "TSS",
    output.format = ".img",
    do.stack = FALSE
  )

#### future projections ----

## load 2100 variables - RCP 2.6
# Has to be in the same name as variables used for model construction
# For the variables not available in the future, using the current conditions
setwd("C:/Users/gabri/Google Drive/Scripts/2100_RCP26_bdmin")
envdata = Sys.glob("*.tif") #Or whatever identifies your files
rcp26 = raster::stack() #empty stack for raster
e <- raster::extent(-80, -20, -40, 10)  #for world usually crop (-180, 180, -90, 90)
for(i in 1:NROW(envdata)){
  tempraster = raster::raster(envdata[i])
  cp_tempraster= raster::crop(tempraster, e)
  rcp26 = raster::stack(rcp26, cp_tempraster)
}
names(rcp26)

setwd("C:/Users/gabri/Google Drive/Mestrado/crispatum/biomod2")
models_proj_2100_rcp26 <- 
  BIOMOD_Projection(
    modeling.output = models,
    new.env = rcp26,
    proj.name = "2100_rcp26_bat",
    binary.meth = "TSS",
    output.format = ".img",
    do.stack = FALSE
  )

ensemble_models_proj_2100_rcp26 <- 
  BIOMOD_EnsembleForecasting(
    EM.output = ensemble_models,
    projection.output = models_proj_2100_rcp26,
    binary.meth = "TSS",
    output.format = ".img",
    do.stack = FALSE
  )

## load 2100 variables - RCP 8.5
setwd("C:/Users/gabri/Google Drive/Scripts/2100_RCP85_bdmin")
envdata = Sys.glob("*.tif") #Or whatever identifies your files
rcp85 = raster::stack() #empty stack for raster
e <- raster::extent(-80, -20, -40, 10)  #for world usually crop (-180, 180, -90, 90)
for(i in 1:NROW(envdata)){
  tempraster = raster::raster(envdata[i])
  cp_tempraster= raster::crop(tempraster, e)
  rcp85 = raster::stack(rcp85, cp_tempraster)
}
names(rcp85)

setwd("C:/Users/gabri/Google Drive/Mestrado/crispatum/biomod2")

models_proj_2100_rcp85 <- 
  BIOMOD_Projection(
    modeling.output = models,
    new.env = rcp85,
    proj.name = "2100_rcp85_bat",
    binary.meth = "TSS",
    output.format = ".img",
    do.stack = FALSE
  )

ensemble_models_proj_2100_rcp85 <- 
  BIOMOD_EnsembleForecasting(
    EM.output = ensemble_models,
    projection.output = models_proj_2100_rcp85,
    binary.meth = "TSS",
    output.format = ".img",
    do.stack = FALSE
  )

#### Species Range Change (SRC) ----
## loading binary projections
bin_proj_current <- 
  stack( 
    c(
      ca = "lcrispatum/proj_current_bat/individual_projections/lcrispatum_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img",
      wm = "lcrispatum/proj_current_bat/individual_projections/lcrispatum_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img"
    )
  )

bin_proj_2100_rcp26 <- 
  stack( 
    c(
      ca = "lcrispatum/proj_2100_rcp26_bat/individual_projections/lcrispatum_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img",
      wm = "lcrispatum/proj_2100_rcp26_bat/individual_projections/lcrispatum_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img"
    )
  )

bin_proj_2100_rcp85 <- 
  stack( 
    c(
      ca = "lcrispatum/proj_2100_rcp85_bat/individual_projections/lcrispatum_EMcaByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img",
      wm = "lcrispatum/proj_2100_rcp85_bat/individual_projections/lcrispatum_EMwmeanByTSS_mergedAlgo_mergedRun_mergedData_TSSbin.img"
    )
  )

## SRC current -> 2100 RCP 26
SRC_current_2100_rcp26 <- 
  BIOMOD_RangeSize(
    bin_proj_current,
    bin_proj_2100_rcp26
  )

SRC_current_2100_rcp26$Compt.By.Models

## SRC current -> 2100 RCP85
SRC_current_2100_rcp85 <- 
  BIOMOD_RangeSize(
    bin_proj_current,
    bin_proj_2100_rcp85
  )

#raster of committed average
R85 <- SRC_current_2100_rcp85[["Diff.By.Pixel"]]@layers[[1]]
R26 <- SRC_current_2100_rcp26[["Diff.By.Pixel"]]@layers[[1]]

writeRaster(R85, filename="SRC_current_2100_rcp85.tif", overwrite = T)
writeRaster(R26, filename="SRC_current_2100_rcp26.tif", overwrite = T)
#diff.by.pixel = -2 given pixel is predicted to be lost by the species
# -1 given pixel is predicted to be stable for the species
# 0 is the given pixel was not occupied, and will not be in the future (absence)
# 1 given pixel was not occupied, and is predicted to be into the future


###################################################################################
#### Model evaluation graph
## get models
models <- get(load("lcrispatum/lcrispatum.lcrispatum_bat.models.out"))

## get models evaluation scores
models_scores <- get_evaluations(models)

## rearrange data
ModEval = melt(models_scores)
ModEval = ModEval[which(ModEval$Var2 == "Testing.data"), ]
head(ModEval)

plasma_pal <- c("#0D0887FF", "#6A00A8FF", "#B12A90FF", "#E16462FF", "#FCA636FF")
fun_col <- colorRampPalette(plasma_pal)
tiff("modeval_ind_run.tiff", width=250, height=160, bg="white", res=500,unit="mm")
ggplot(ModEval, aes(x = Var1, y = value, color = interaction(Var5,Var3))) +
  scale_color_manual("", values = rep(fun_col(3), each = 10)) +
  geom_boxplot(size = 0.5) +
  labs(x = "", y = "") +
  theme_minimal(base_size = 10)
dev.off()

###################################################################################
#### Ensemble model evaluation
## get ensemble model
ensemble_models <- get(load("lcrispatum/lcrispatum.lcrispatum_batensemble.models.out"))

## get ensemble model evaluation scores
ensemble_models_scores <- get_evaluations(ensemble_models)

## rearrange data
EnsModEval = melt(ensemble_models_scores)
EnsModEval = EnsModEval[which(EnsModEval$Var2 == "Testing.data"), ]
head(EnsModEval)

###################################################################################
#### Variable importance - Ensemble model
## get ensemble model variable importance
ensemble_models_var_import <- get_variables_importance(ensemble_models)

## rearrange data
EnsVarImport = melt(ensemble_models_var_import)
EnsVarImport$Var3 = sub(".*EM", "", EnsVarImport$Var3)
EnsVarImport$Var3 = sub("_merged.*", "", EnsVarImport$Var3)

varimp_ca <- EnsVarImport %>%
  filter(Var3 == "caByTSS")

head(varimp_ca)
varimp_ca <- varimp_ca[,c(1,2,4)]

#I want the importance in percentage
rand1 <- varimp_ca %>%
  filter(Var2 == "rand1") 
sum_rand1 <- sum(rand1$value)
rand1$perc <- rand1$value*100/sum_rand1

rand2 <- varimp_ca %>%
  filter(Var2 == "rand2") 
sum_rand2 <- sum(rand2$value)
rand2$perc <- rand2$value*100/sum_rand2

rand3 <- varimp_ca %>%
  filter(Var2 == "rand3") 
sum_rand3 <- sum(rand3$value)
rand3$perc <- rand3$value*100/sum_rand3

varimp <- rbind(rand1, rand2, rand3)
varimp <- varimp[,c(1:2,4)]

varimp <- varimp %>%
          mutate(Var1 = case_when(
            Var1 == "BO2_nitrateltmax_bdmin_lonlat" ~ "Maximum Nitrate",
            Var1 == "BO2_templtmax_bdmin_lonlat" ~ "Maximum Temperature",
            Var1 == "BO2_salinityltmax_bdmin_lonlat" ~ "Maximum Salinity",
            Var1 == "BO_bathymin_lonlat" ~ "Bathymetry",
            Var1 == "calcsat_av_1951_to_2000" ~ "Calcite Saturation",
            Var1 == "ph_av_1951_to_2000" ~ "Mean pH",
            Var1 == "heatwave_median_intensity" ~ "Heatwave Median Intensity")) %>%
            arrange(-perc)

varimp$Var1 <- factor(varimp$Var1, levels = c("Heatwave Median Intensity",
                                              "Mean pH",
                                              "Calcite Saturation",
                                              "Bathymetry",
                                              "Maximum Salinity",
                                              "Maximum Temperature",
                                              "Maximum Nitrate"))

varplot <- ggplot(data=varimp, aes(x=perc, y=Var1, fill=Var2)) +
  geom_bar(stat="identity", position = position_dodge2(reverse = T))+
  scale_fill_brewer(palette = "YlGnBu", labels = c("1", "2", "3")) +
  labs(x = "Relative Variable Importance (%)", y = "Environmental Predictor", fill = "Permutation run") +
  theme_minimal()

tiff("varimp.tiff", width=250, height=160, bg="white", res=500,unit="mm")
varplot + theme(legend.position="bottom",
                text = element_text(size=16)) 
dev.off()
###################################################################################

## get coefficient of variation
ras.cv = raster("lcrispatum/proj_current_bat/individual_projections/lcrispatum_EMcvByTSS_mergedAlgo_mergedRun_mergedData.img")

ras.cv.pts = rasterToPoints(ras.cv)

writeRaster(ras.cv, filename="cv_lcrispatum.tif", overwrite = T)



