set.seed(123)

library("randomForest")

model <- readRDS("input/model.RDS")
test_dat_pca <- readRDS("input/pca.RDS")
output <- read.csv("/input/template.csv", row.names = 1)

predict <- predict(model, newdata = as.data.frame(test_dat_pca))

for(cell_line in rownames(output)){
  for(drug in colnames(output)){
    output[cell_line, drug] <- predict[grepl(paste0(cell_line, "_", drug), rownames(test_dat_pca))][[1]]
  }
}
print("ajsdlkfjj")
write.csv(output, "DreamChallenge/predictions.csv")
