args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    paste(
      "Usage:",
      "Rscript run_immunophenotype.R <input_file> <output_dir> [perm_times]",
      "\nExample:",
      "Rscript run_immunophenotype.R",
      "'tables/08_Tumor Purity/TIP_input.txt'",
      "'tables/08_Tumor Purity/TIP_immunophenotype' 100"
    )
  )
}

input_file <- normalizePath(args[1], mustWork = TRUE)
output_dir <- args[2]
perm_times <- if (length(args) >= 3) as.integer(args[3]) else 100L

if (is.na(perm_times) || perm_times < 1) {
  stop("perm_times must be a positive integer.")
}

repo_dir <- getwd()
input_dir <- dirname(input_file)
input_name <- basename(input_file)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

source(file.path(repo_dir, "1.MainFunction", "ErrorProcess.R"))
source(file.path(repo_dir, "2.ImmuneActivityScore", "immunityScore.server.R"))

check_result <- checkError(
  filePath = input_dir,
  fileName = input_name,
  codePath = repo_dir,
  format.of.file = "TPM",
  type.of.data = "RNA-seq",
  saveDir = output_dir
)

if (!is.null(check_result)) {
  stop(paste("Input check failed:", check_result))
}

expr <- get(load(file.path(output_dir, "expression.afterIDConvert.RData")))
sample_number <- ncol(expr)

immunityScore.server(
  codePath = file.path(repo_dir, "2.ImmuneActivityScore"),
  filePath = input_dir,
  fileName = input_name,
  saveDir = output_dir,
  sampleNumber = sample_number,
  permTimes = perm_times,
  type.of.data = "RNA-seq",
  format.of.file = "TPM"
)

cat("Analysis finished.\n")
cat("Input:", input_file, "\n")
cat("Output:", normalizePath(output_dir), "\n")
cat("Samples:", sample_number, "\n")
cat("Key result:", file.path(output_dir, "ssGSEA.normalized.score.txt"), "\n")
