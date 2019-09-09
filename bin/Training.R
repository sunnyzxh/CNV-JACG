library(randomForest)
pdf(file="../lib/training.pdf")
del <- read.table("../lib/DEL.training.txt", header=T)
dim(del)
set.seed(12)
del_classifier = randomForest(Decision ~ ., data=del, ntree=200, mtry=10, importance=TRUE)
del_classifier
saveRDS(del_classifier, "../lib/trained_model_del.rds")
plot(del_classifier)
varImpPlot(del_classifier)

dup <- read.table("../lib/DUP.training.txt", header=T)
dim(dup)
set.seed(12)
dup_classifier = randomForest(Decision ~ ., data=dup, ntree=450, mtry=10, importance=TRUE)
dup_classifier
saveRDS(dup_classifier, "../lib/trained_model_dup.rds")
plot(dup_classifier)
varImpPlot(del_classifier)
