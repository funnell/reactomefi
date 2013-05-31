library(RCurl)
library(XML)

source("AllClasses.R")
source("AllGenerics.R")
source("ReactomeFIService.R")

genes <- scan("genes.txt", what = "character")

service <- new("ReactomeFIService", version="2012")
fis <- queryFIs(service, genes)
cl <- cluster(service, fis)
annot <- annotateGeneSet(service, genes, "Pathway")
