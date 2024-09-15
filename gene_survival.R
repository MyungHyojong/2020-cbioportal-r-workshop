# 특정 연구의 데이터를 불러오기 위한 API 연결
cbio <- cBioPortal()

# 특정 연구 ID와 분자 프로필 ID
studyId <- "acc_tcga"
profileId <- "acc_tcga_mutations"

# 연구에서 사용 가능한 샘플 목록을 불러옵니다
sample_list <- allSamples(api = cbio, studyId = studyId)

# TP53 유전자의 돌연변이 데이터를 불러옵니다 (샘플 ID를 사용)
entrez_gene_id <- 7157  # TP53의 Entrez ID
mutation_data <- mutationData(
  api = cbio,
  molecularProfileIds = profileId,
  entrezGeneIds = entrez_gene_id,
  sampleIds = sample_list$sampleId  # 샘플 ID 목록 추가
)

# 임상 데이터 가져오기 (생존 정보 포함)
clinical_data <- clinicalData(
  api = cbio,
  studyId = studyId
)

# 유전자 돌연변이 여부 추가 (TP53 돌연변이 여부)
clinical_data$TP53_status <- ifelse(clinical_data$sampleId %in% mutation_data[[1]]$sampleId, "Mutant", "Wild-Type")
table(clinical_data$TP53_status)

# 생존 분석을 위한 Surv 객체 생성
surv_object <- Surv(time = as.numeric(clinical_data$OS_MONTHS), event = clinical_data$OS_STATUS == "1:DECEASED")

# 생존 분석 모델 적합 (TP53 돌연변이 여부에 따라)
fit <- survfit(surv_object ~ TP53_status, data = clinical_data)

# 생존 곡선 시각화
ggsurvplot(
  fit,
  data = clinical_data,
  risk.table = TRUE,        # 위험 테이블 추가
  pval = TRUE,              # p-value 표시
  conf.int = TRUE,          # 신뢰 구간 표시
  legend.labs = c("Wild-Type", "Mutant"), # 범례 라벨
  xlab = "Time (Months)",   # x축 라벨
  ylab = "Survival Probability"  # y축 라벨
)
