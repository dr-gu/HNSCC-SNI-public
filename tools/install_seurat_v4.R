dir.create("~/R_envs/seurat_v4", recursive = TRUE)
setwd("~/R_envs/seurat_v4")

install.packages("renv")
renv::init(bare = TRUE)

renv::install("Seurat@4.3.0")

setwd("~/R_envs/seurat_v4")
renv::activate()

install.packages("leiden")
renv::install("Seurat@4.3.0")

renv::snapshot()

file.edit("~/.Rprofile")
