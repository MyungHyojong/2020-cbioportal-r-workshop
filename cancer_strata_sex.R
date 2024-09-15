## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(cache = TRUE)

## ----include=TRUE,results="hide",message=FALSE,warning=FALSE------------------
library(cBioPortalData)
library(AnVIL)

## ----eval=FALSE---------------------------------------------------------------
#  citation("MultiAssayExperiment")
#  citation("cBioPortalData")

## -----------------------------------------------------------------------------
cbio <- cBioPortal()
studies <- getStudies(cbio, buildReport = TRUE)
head(studies)

## ----message=FALSE,warning=FALSE----------------------------------------------
## Use ask=FALSE for non-interactive use
laml <- cBioDataPack("laml_tcga", ask = FALSE)
laml

## ----warning=FALSE------------------------------------------------------------
acc <- cBioPortalData(api = cbio, by = "hugoGeneSymbol", studyId = "acc_tcga",
                      genePanelId = "IMPACT341",
                      molecularProfileIds = c("acc_tcga_rppa", "acc_tcga_linear_CNA")
)
acc

## -----------------------------------------------------------------------------
metadata(acc)
acc[[1]]

## ----echo=FALSE---------------------------------------------------------------
cat(
  "Our testing shows that '%s' is not currently building.\n",
  "  Use 'downloadStudy()' to manually obtain the data.\n",
  "  Proceed anyway? [y/n]: y"
)

## ----eval=FALSE---------------------------------------------------------------
#  removeCache("laml_tcga")

## ----eval=FALSE---------------------------------------------------------------
#  unlink("~/.cache/cBioPortalData/")

## ----message=FALSE,warning=FALSE----------------------------------------------
library(survival)
library(survminer)

## -----------------------------------------------------------------------------
table(colData(laml)$OS_STATUS)
class(colData(laml)$OS_MONTHS)

## -----------------------------------------------------------------------------
collaml <- colData(laml)
collaml[collaml$OS_MONTHS == "[Not Available]", "OS_MONTHS"] <- NA
collaml$OS_MONTHS <- as.numeric(collaml$OS_MONTHS)
colData(laml) <- collaml

## -----------------------------------------------------------------------------
fit <- survfit(
  Surv(OS_MONTHS, as.numeric(substr(OS_STATUS, 1, 1))) ~ SEX,
  data = colData(laml)
)
ggsurvplot(fit, data = colData(laml), risk.table = TRUE)

## -----------------------------------------------------------------------------
sessionInfo()
