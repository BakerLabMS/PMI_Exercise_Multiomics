library("httr")
library("jsonlite")
library(fingerprint)
library(rcdk)
library(rprojroot)

setwd(find_root("Github_PMI_Exercise.Rproj"))
setwd("PMI_EX_analyses/PMI_EX_data")


df <- read.csv("PMI_multiomic_SMILES.csv", stringsAsFactors = F)

head(df)

sum(df$smi == "")

df.c <- df[df$smi != "", ]
df.na <- df[df$smi == "", ]

head(df.c)

fing <- parse.smiles(df.c$smi)
fing <- lapply(fing, get.fingerprint, type = "circular")
fp.sim <- fp.sim.matrix(fing)
row.names(fp.sim) <- df.c$abbrev
fp.dist <- as.dist(1 - fp.sim)
cls <- hclust(fp.dist)
plot(cls)

save.image("fp_PMI_multiomics.rdata")

