library(hierfstat)
library(gaston)

### Open VCF file
vcf <- read.VCF("/Users/perrinekergoat/Master_Project/Code/Data/Simu.vcf")

print(vcf)
dim(vcf)

### Convert VCF file to matrix
vcf_mat <- as.matrix(vcf)

### Get dosage matrix 
dosage <- beta.dosage(vcf_mat)

### Matrix output
png(file="/Users/perrinekergoat/Master_Project/Code/Data/HM.png",
    width=20000, height=20000)
heatmap(dosage, scale = "none", keep.dendro = FALSE) # enlever le tri des individus 
dev.off()

