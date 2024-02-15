# Download and extract the data (from command line in the newly created directory)
# `wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE122nnn/GSE122126/suppl/GSE122126_RAW.tar`
# `tar -xvf GSE122126_RAW.tar`

library(minfi)
library(here)

setwd(here("data"))
dir.create("Moss_et_al-NatComm-2018", showWarnings=FALSE)
setwd("Moss_et_al-NatComm-2018")
if(length(list.files(".")) == 0) {
    print("Downloading GEO files")
    cmd <- "curl ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE122nnn/GSE122126/suppl/GSE122126_RAW.tar -o GSE122126_RAW.tar"
    system(cmd)
    print("Done")
} else {
    print("Files are already downloaded.")
}

if(length(list.files(".")) == 1) {
    cmd <- "tar -xvf GSE122126_RAW.tar"
    system(cmd)
}

dir.create("results", showWarnings=FALSE)

if (!file.exists("results/healthy.rds")) {
  files <- list.files(".")
  foi <- c("GSM3455776", "GSM3455777", "GSM3455778", "GSM3455779", "GSM3455780", "GSM3455781", "GSM3455782", "GSM3455783")
  foi_legend <- c("young_male_1", "young_male_2", "young_female_1", "young_female_2", "old_male_1", "old_male_2", "old_female_1", "old_female_2")
  array <- c()
  f <- c()
  slide <- c()
  for (i in foi) {
    array <- c(array, unlist(strsplit(files[grep("Grn", files)][grep(i, files[grep("Grn", files)])], "_"))[3])
    f <- c(f, files[grep("Grn", files)][grep(i, files[grep("Grn", files)])],files[grep("Red", files)][grep(i, files[grep("Red", files)])])
    slide <- c(slide, as.numeric(unlist(strsplit(files[grep("Grn", files)][grep(i, files[grep("Grn", files)])], "_"))[2]))
  }

  targets <- data.frame("Sample_Name" = paste0(slide, "_", array),
                        "Sample_Well" = NA,
                        "Sample_Plate" = NA,
                        "Sample_Group" = NA,
                        "Pool_ID" = foi_legend,
                        "Array" = array,
                        "Slide" = slide,
                        "Basename" = file.path(paste0(foi, "_", slide, "_", array)))
  RGset <- read.metharray.exp(targets = targets)
  GRset.funnorm <- preprocessFunnorm(RGset, bgCorr = TRUE, dyeCorr = TRUE, ratioConvert=T)
  Mdata.mapped <- mapToGenome(GRset.funnorm)
  beta <- getBeta(Mdata.mapped)
  annotation <- getAnnotation(Mdata.mapped)
  position <- annotation[1:2]
  total <- merge(position, beta, by=0)
  saveRDS(total, file.path("results", "healthy.rds"))
}
