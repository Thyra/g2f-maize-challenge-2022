library(rjson)
library(data.table)
library(randomForest)

d = setDT(fread("../../processed_data/combined_mat.csv"))
x = data.frame(d[, 12:ncol(d)])
x = x[,colSums(is.na(x))==0]
# Adjust yield to 15.5% moisture
y = d$Yield_Mg_ha-(d$Yield_Mg_ha*(d$Grain_Moisture-15.5)/100)

split = fromJSON(file = "../../processed_data/train_test_split_v2.json")

rmses = data.frame(t.set = character(), rmse=numeric())
for(t.set in names(split)) {
  print(t.set)
  train.inds = split[[t.set]]$train
  test.inds = split[[t.set]]$test
  x.train = x[train.inds,]
  y.train = y[train.inds]
    
  rf = randomForest(x.train, y.train)
  predValues <- predict(rf, x[test.inds,])
  rmse = sqrt(mean((y[test.inds] - predValues)^2))
  print(rmse)
  rmses = rbind(rmses, list(t.set = t.set, rmse = sqrt(mean((y[test.inds] - predValues)^2))))
}

write.csv(rmses, "rmses.csv", row.names=F)
