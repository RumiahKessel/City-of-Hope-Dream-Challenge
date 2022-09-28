library('stringr')
library("randomForest")

set.seed(123)

data <- read.csv("DreamChallenge/data_achilles_shRNA_DREAMv2.csv", header = TRUE, row.names = 1)
cellLine_revealed <- vector(mode = "list", length = 11)
cellLine_concealed <- vector(mode = "list", length = 11)
my_files <- list.files(path = "DreamChallenge/drug_activity_challenge_rnaseq/rnaseq/", full.names = T)
my_files_concealed <- list.files(path = "DreamChallenge/rnaseq_concealed/", full.names = T)

linMap <- function(x, from, to){
  (x - min(x)) / max(x - min(x)) * (to - from) + from
}

subset_test_dat_by_excluded_drug <- function(dat, drug){
  subset <- dat[!grepl(drug, row.names(dat)),]
  return(subset)
}

subset_cell_line_by_excluded_drug <- function(dat, drug){
  subset <- dat[, !grepl(drug, names(dat))]
  return(subset)
}


for (i in seq_along(my_files)) {
  cellLine_revealed[[i]] <- read.csv(file = my_files[i], header = TRUE, row.names = 1)
}

for(i in seq(1:length(cellLine_revealed))){
  print(i)
  cellLine_revealed[[i]] <- cellLine_revealed[[i]][order(rownames(cellLine_revealed[[i]])), order(colnames(cellLine_revealed[[i]]))]
}

for(i in seq(1:length(cellLine_revealed))){
  print(i)
  cellLine_revealed[[i]] <- subset_cell_line_by_excluded_drug(cellLine_revealed[[i]], "UNTREATED")
  cellLine_revealed[[i]] <- subset_cell_line_by_excluded_drug(cellLine_revealed[[i]], "DMSO")
}

#concealed processing
for (i in seq_along(my_files_concealed)) {
  cellLine_concealed[[i]] <- read.csv(file = my_files_concealed[i], header = TRUE, row.names = 1)
}

for(i in seq(1:length(cellLine_concealed))){
  print(i)
  cellLine_concealed[[i]] <- cellLine_concealed[[i]][order(rownames(cellLine_concealed[[i]])), order(colnames(cellLine_concealed[[i]]))]
}

for(i in seq(1:length(cellLine_concealed))){
  print(i)
  cellLine_concealed[[i]] <- subset_cell_line_by_excluded_drug(cellLine_concealed[[i]], "untreated")
  cellLine_concealed[[i]] <- subset_cell_line_by_excluded_drug(cellLine_concealed[[i]], "dmso")
}

data <- data[order(rownames(data)), ]

ctrp_data <- read.table("DreamChallenge/CTRPV2.2.csv", sep = "\t", header = FALSE, quote="", row.names = 2)

ctrp_data <- ctrp_data[,order(ctrp_data[2,])]

ctrp_data <- ctrp_data[-c(1,3,4,5),-c(2,3)]

ctrp_data <- ctrp_data[order(rownames(ctrp_data)),order(ctrp_data[1,])]

intersect <- rownames(cellLine_concealed[[1]])[unlist(rownames(cellLine_concealed[[1]])) %in% unlist(rownames(data))]
intersect <- intersect[!str_detect(intersect, "&")]

drugs <- unique(str_extract(colnames(cellLine_revealed[[1]]), regex("[:alnum:]+")))

superfinal <- data.frame(matrix(ncol = 23339, nrow = 0))
colnames(superfinal) <- rownames(cellLine_revealed[[1]])

test_dat <- data.frame(matrix(ncol = 23339, nrow = 0))
colnames(test_dat) <- rownames(cellLine_concealed[[1]])

testable_values <- which(colnames(data) %in% ctrp_data[1,])

i <- 1
for(cell_line in colnames(data)[testable_values]){
  print(i)
  tested <- !is.na(data[intersect,cell_line]) & (data_ccle[intersect, cell_line] != 0)
  genes <- intersect[tested]
  genes <- genes[order(genes)]
  length_genes <- length(genes)
  for(drug in drugs){
    avg <- data.frame(matrix(ncol= 22, nrow = length_genes))
    for(k in seq(1:11)){
      avg[,2*k-1] <- cellLine_revealed[[k]][genes,grep(colnames(cellLine_revealed[[k]]), pattern = drug)][1]
      avg[,2*k] <- cellLine_revealed[[k]][genes,grep(colnames(cellLine_revealed[[k]]), pattern = drug)][2]
    }
    
    #final_training_data[[i]][drug, ] <- apply(avg, 1, mean) * data[genes, cell_line]
    superfinal[paste0(cell_line, "_", drug),genes] <- apply(avg, 1, mean) * data[genes, cell_line]
     
  }
  i <- i + 1
}

drugs <- unique(str_extract(colnames(cellLine_concealed[[1]]), regex(	"cmpd_[[:alpha:]]+")))

i <- 1
for(cell_line in colnames(data)){
  print(i)
  tested <- !is.na(data[intersect,cell_line])
  genes <- intersect[tested]
  genes <- genes[order(genes)]
  length_genes <- length(genes)
  for(drug in drugs){
    avg <- data.frame(matrix(ncol= 22, nrow = length_genes))
    for(k in seq(1:11)){
      avg[,2*k-1] <- cellLine_concealed[[k]][genes,grep(colnames(cellLine_concealed[[k]]), pattern = drug)][1]
      avg[,2*k] <- cellLine_concealed[[k]][genes,grep(colnames(cellLine_concealed[[k]]), pattern = drug)][2]
    }
    
    #final_training_data[[i]][drug, ] <- apply(avg, 1, mean) * data[genes, cell_line]
    test_dat[paste0(cell_line, "_", drug),genes] <- apply(avg, 1, mean) * data[genes, cell_line]
    
  }
  i <- i + 1
}

superfinal[is.na(superfinal)] <- 0

superfinal.pca <- prcomp(superfinal, Center = TRUE, scale = TRUE)
superfinal <- superfinal.pca$x
superfinal <- superfinal[,1:1000]

#creat solution set

testing <- vector(mode = "numeric", length = 14144)

i <- 1
for(cell_line in colnames(data)[testable_values]){
  for(drug in drugs){
    testing[i] <- ctrp_data[tolower(drug), match(tolower(cell_line), tolower(ctrp_data[1,]))]
    i <- i + 1
  }
}

testing <- as.numeric(testing)

superfinal <- superfinal[!is.na(testing),]

testing <- testing[!is.na(testing)]

testing <- testing[!is.nan(testing)]

superfinal <- superfinal[!is.nan(testing),]

testing <- linMap(testing,1,0)

#setting up the training

train_rows <- sample(1:7673, .9*7673)

superfinal.train <- superfinal[train_rows,]
superfinal.test <- superfinal[-train_rows,]

testing.train <- testing[train_rows]
testing.test <- testing[-train_rows]

model <- randomForest(formula = testing.train ~ ., data = as.data.frame(superfinal.train))

test_dat[is.na(test_dat)] <- 0

test_dat_pca <- predict(superfinal.pca, newdata = test_dat)

test_dat_pca <- test_dat_pca[,1:1000]

predict <- predict(model, newdata = as.data.frame(test_dat_pca))

output <- read.csv("/DreamChallenge/template_final.csv")

for(cell_line in rownames(output)){
  for(drug in colnames(output)){
    output[cell_line, drug] <- predict[grepl(past0(cell_line, "_", drug), rownames(test_dat))]
  }
}

write.csv(output, "/DreamChallenge/predictions.csv")

saveRDS(superfinal, "DreamChallenge/Training_Data/Final_Trainingv3.RDS")
saveRDS(testing, "DreamChallenge/Training_Data/Final_Testingv1.RDS")
saveRDS(superfinal.pca, "DreamChallenge/superfinal_pca.RDS")
saveRDS(test_dat, "DreamChallenge/test_dat.RDS")
saveRDS(pca.train, "DreamChallenge/pca_train.RDS")
saveRDS(model, "DreamChallenge/model.RDS")