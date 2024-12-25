library(devtools)
library(usethis)
library(SPARK)

# get coord and count
sp_count <- read.csv("D:\\CRC_meta1_counts_sp.csv", header = TRUE, sep = ",", stringsAsFactors = TRUE)
location <- read.csv("D:\\CRC_meta1_coord_sp.csv", header = TRUE, sep = ",", stringsAsFactors = TRUE)
sp_count <- as.matrix(sp_count)
location <- as.matrix(location)
dim(sp_count)
dim(location)
# Removing mitochondrial genes
mt_idx      <- grep("mt-",rownames(sp_count))
if(length(mt_idx)!=0){
  sp_count    <- sp_count[-mt_idx,]
}
# sparkX
sparkX <- sparkx(sp_count,location,numCores=1,option="mixture")
head(sparkX$res_mtest)
write.csv(sparkX$res_mtest,"D:\\CRC_meta1_Spark.csv")