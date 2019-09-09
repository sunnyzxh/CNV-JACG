library(randomForest)
args <- commandArgs(T)

trained    <- args[1]
inputFile  <- args[2]
outputFile <- args[3]

set.seed(12)

topredict <- read.table(inputFile, header=T)
trained_model <- readRDS(trained)

Prediction <- predict(trained_model,topredict)
Output = cbind(Prediction, topredict)

write.table(Output,file=outputFile,row.names = FALSE,sep = "\t",quote = FALSE)
