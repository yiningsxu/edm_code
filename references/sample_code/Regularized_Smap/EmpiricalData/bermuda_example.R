rm(list=ls(all=TRUE)) 
suppressMessages(library(Matrix))
suppressMessages(library(quantreg))
suppressMessages(library(parallel))
suppressMessages(library(compiler))
suppressMessages(library(lars))
suppressMessages(library(elasticnet))
suppressMessages(library(caret))
options(warn=-1)
#####################################################################################
source('Auxiliar.r')
source('elastic_net_fit.r')
source('ridge.r')
source('LOOCV.r')
source('KernelFunctions.r')
source('OutOfSample.r')
source('TrainingError.r')
###########
ShowPlot = FALSE
lags = TRUE
ModelName = 'bacteria_ts_bermuda2011'
FileName = paste(ModelName, '.txt', sep = '')
###################################
logspace <- function(d1, d2, n) exp(log(10)*seq(d1, d2, length.out=n)) 
std_err <- function(x) sd(x)/sqrt(length(x))
############# Choose the kernel
Kernel.Options = c('Exponential.Kernel', 'Epanechnikov.Kernel', 'TriCubic.Kernel', 'Matern.Kernel')
Regression.Kernel = Kernel.Options[1]
############# Parameters for cross validation
lambda = logspace(-3,0,15) 
if(Regression.Kernel == 'Exponential.Kernel'){
  tht = seq(from = 0., to = 10, length = 30)         
}else{
  tht = seq(from = 0.1, to = 3, length = 20)     
}
parameters_on_grid = expand.grid(tht, lambda)     
### Read Time series
d = as.matrix(read.table(FileName, header= T))
######################
original.Embedding = c('Pro', 'Syn','Piceu')
original.TargetList = original.Embedding
d = d[, original.Embedding]
#### Here you take combinations of lags (best lag are 1 - 2)
x.lag = 1; y.lag = 2; z.lag = 1
sp.lag.selection = c(x.lag, y.lag, z.lag)
lagged.time.series = make.lagged.ts(d, sp.lag.selection)
d = lagged.time.series$time.series
original.col = lagged.time.series$original.variables
if(lags == TRUE){ var.sel = original.col; }else{ var.sel = colnames(d)}
##### Names and embedding in the laged dataset
if(lags == TRUE){ colnames(d) = Embedding =  TargetList = LETTERS[1:ncol(d)]}else{
  Embedding =  TargetList = original.Embedding
}
##### length of training and test set
length.testing = 2
length.training = nrow(d) - length.testing
#### Preserve training for the interactions
ts.train.preserved = d[1:length.training, var.sel]
std.ts.train = Standardizza(ts.train.preserved)
#### Preserve testing for the test (you want your algorithm to learn the real structure of the model)
ts.test.preserved = d[(length.training + 1):nrow(d), var.sel]
#### Training set:
d.training = Standardizza(d[1:length.training, ])
#### You now need to standardize the test set using mean and sd of the training set
d.testing = Standardizza.test(ts.test.preserved,ts.train.preserved)
############## Prepare for parallel computing
Lavoratori = detectCores() - 2
cl <- makeCluster(Lavoratori, type = "FORK")
####
RegressionType = 'ELNET_fit'
#RegressionType = 'ridge_fit'
alpha = 0.85

### should you compute all the variables or not?
BestModel = BestModelLOOCV(cl, d.training, TargetList, Embedding, parameters_on_grid, RegressionType,alpha)
BestCoefficients = BestModel$BestCoefficients
BestParameters = BestModel$BestParameters
### Forecast
out.of.samp.ELNET = out_of_sample_sequence(cl, BestCoefficients, 
                                           BestParameters$BestTH,
      	                                   BestParameters$BestLM, 
              	                           d.training, length.testing)
prd = out.of.samp.ELNET$out_of_samp
stopCluster(cl)
###### Now check the quality of the out-of-sample forecast of only the unlaged variables
prd =  prd[,var.sel]
###############################
jacobiano.inference = take.coeff(BestCoefficients, var.sel, original.Embedding)
###### Now check the in-sample error
Reconstruction = ReconstructionOfTrainingSet(d.training, BestCoefficients)
Reconstruction = Reconstruction[,var.sel]
##### RMSE
Training.RMSE = ComputeTrainingError(d.training, BestCoefficients, var.sel)$rmse
rmse.test = compute.rmse.test(d.testing, prd)
rmse.test.naive = naive.forecast(std.ts.train[nrow(d.training),],d.testing)
##### R^2
R2.training = as.numeric(postResample(pred = Reconstruction, 
                                      obs = d.training[1:nrow(d.training)-1,var.sel])['Rsquared'])
R2.testing = as.numeric(postResample(pred = prd, obs = d.testing)['Rsquared'])
###############################
if(ShowPlot == TRUE){
  source('PlotFunctions.r')
  Plot_out(d.testing, prd)
  Plot_AllSpeciesTraining(std.ts.train[1:nrow(d.training)-1,], Reconstruction)
  Plot_SingleSpeciesTraining(std.ts.train[1:nrow(d.training)-1,], 
                             Reconstruction, 1)
}
for(sp in 1:ncol(prd)){
  cat('Error on Species', sp, ':', sqrt(mean((d.testing[,sp] - prd[,sp])^2)), '\n')
}
