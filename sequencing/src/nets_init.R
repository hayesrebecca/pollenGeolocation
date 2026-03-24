mainDir <- "../pnw_survey/data"
subDir <- "networks"

if (file.exists(subDir)){
} else {
    dir.create(file.path(mainDir, subDir),
               showWarnings = FALSE)
}
save.dir <- file.path(mainDir, subDir)

source('../pnw_survey/data_prep/sequence_prep/src/illumSplit.R')
source('../pnw_survey/data_prep/sequence_prep/src/calcFuncUniqOrig.R')
source('../pnw_survey/data_prep/sequence_prep/src/misc.R')
load('../pnw_survey/data/spec_RBCL_16s.Rdata')
