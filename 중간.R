getwd()
# 1. raw 폴더의 파일 목록 읽기
files <- list.files("raw", full.names = TRUE)

# 2. sorted_R 폴더로 파일 자동 정렬 및 복사
for (f in files) {
  fname <- basename(f)
  
  cond <- ifelse(grepl("Control", fname), "Control", "Treated")
  r_dir <- ifelse(grepl("R1", fname), "R1", "R2")
  
  target_dir <- file.path("sorted_R", cond, r_dir)
  if (!dir.exists(target_dir)) {
    dir.create(target_dir, recursive = TRUE)
  }
  
  file.copy(f, file.path(target_dir, fname), overwrite = TRUE)
}

# 3. Bash 결과와 R 결과의 파일 목록/개수 일치 여부 비교
bash_files <- list.files("sorted_bash", recursive = TRUE)
r_files <- list.files("sorted_R", recursive = TRUE)

cat("파일 개수 일치 여부:", length(bash_files) == length(r_files), "\n")
cat("파일 경로 및 이름 완전 일치 여부:", all(bash_files == r_files), "\n")




sequencing <- data.frame(
  platform = c("Illumina", "Oxford Nanopore", "PacBio HiFi"),
  yield_Gb = c(18, 30, 12),
  host_fraction = c(0.20, 0.70, 0.10),
  cost_USD = c(150, 200, 300)
)
budget <- 600
# 1. 예산 내 정수 샘플 수 계산 (floor 함수 사용)
sequencing$samples <- floor(budget / sequencing$cost_USD)

# 2. 샘플당 순수 미생물 데이터량 계산 (1 - host_fraction)
sequencing$microbial_fraction <- 1 - sequencing$host_fraction
sequencing$microbial_yield_per_sample <- sequencing$yield_Gb * sequencing$microbial_fraction

# 3. 플랫폼별 총 미생물 데이터량(Gb) 계산
sequencing$total_microbial_Gb <- sequencing$samples * sequencing$microbial_yield_per_sample

# 결과 출력
print(sequencing[, c("platform", "samples", "total_microbial_Gb")])

library(dada2)
rawpath <- system.file("extdata", package = "dada2")
fnFs <- file.path(rawpath,
                  c("sam1F.fastq.gz", "sam2F.fastq.gz"))
fnRs <- file.path(rawpath,
                  c("sam1R.fastq.gz", "sam2R.fastq.gz"))
sample.names <- c("sam1", "sam2")
plotQualityProfile(fnFs)
plotQualityProfile(fnRs)

filtFs <- file.path(tempdir(), paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(tempdir(), paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen = c(230, 160), maxEE = c(2, 2),
                     maxN = 0, truncQ = 2, rm.phix = TRUE,
                     compress = TRUE, multithread = FALSE)
# 에러 모델 학습
errF <- learnErrors(filtFs, multithread = FALSE)
errR <- learnErrors(filtRs, multithread = FALSE)
# Denoising (ASV 추정)
dadaFs <- dada(filtFs, errF, multithread = FALSE)
dadaRs <- dada(filtRs, errR, multithread = FALSE)
# Paired-end 병합 (Merge)
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)
# ASV 테이블 생성 및 키메라 제거
seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = FALSE, verbose = TRUE)

# 단계별 Read 수 추적 표 작성
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

cat("=== Read Tracking Table ===\n")
print(track)

# 키메라 제거 개수 계산
num_asvs_raw <- ncol(seqtab)
num_asvs_nochim <- ncol(seqtab.nochim)
removed_asvs <- num_asvs_raw - num_asvs_nochim

reads_raw <- sum(seqtab)
reads_nochim <- sum(seqtab.nochim)
removed_reads <- reads_raw - reads_nochim

cat("\n[키메라 제거 결과]\n")
cat("제거된 키메라 ASV 개수:", removed_asvs, "개\n")
cat("제거된 키메라 Read 개수:", removed_reads, "개\n")

BiocManager::install(c("DECIPHER", "phangorn"), update = FALSE, ask = FALSE)

library(Biostrings)
library(DECIPHER)
library(phangorn)

# 1. 리드 수 상위 5개 ASV 추출 및 이름 지정
top5_idx <- order(colSums(seqtab.nochim), decreasing = TRUE)[1:5]
top5_seqs <- colnames(seqtab.nochim)[top5_idx]
names(top5_seqs) <- paste0("ASV", 1:5)

# 2. DECIPHER 서열 정렬
dna <- Biostrings::DNAStringSet(top5_seqs)
alignment <- AlignSeqs(dna)

# 3. phangorn Neighbor-Joining 계통수 생성
phang <- phyDat(as(alignment, "matrix"), type = "DNA")
dm <- dist.ml(phang)
tree <- NJ(dm)

# 4. 트리 시각화 및 Newick(.nwk) 파일 저장
plot(tree, main = "Neighbor-Joining Tree of Top 5 ASVs")
write.tree(tree, file = "top5_ASVs.nwk")

#_____________________________________________

condition <- c("Untreated", "Treated")
total_reads <- c(20000000, 20000000)
host_fraction <- c(0.99, 0.90)
taxon_X_fraction <- c(0.10, 0.91)

# 1. 미생물 총 리드 수 계산 (Total reads * (1 - host_fraction))
microbial_reads <- total_reads * (1 - host_fraction)

# 2. Taxon_X 리드 수 계산 (Microbial reads * taxon_X_fraction)
taxon_X_reads <- microbial_reads * taxon_X_fraction

# 3. Taxon_X 제외 후 남은 유용한 미생물 리드 수
reads_excl_taxon_X <- microbial_reads - taxon_X_reads

# 4. Fold Change 계산 (Treated / Untreated)
fc_microbial <- microbial_reads[2] / microbial_reads[1]
fc_excl_taxon_X <- reads_excl_taxon_X[2] / reads_excl_taxon_X[1]

# 결과 출력
df_res <- data.frame(
  condition, total_reads, microbial_reads, taxon_X_reads, reads_excl_taxon_X
)
print(df_res)

cat("\nMicrobial Reads Fold Change (Treated / Untreated):", fc_microbial, "배\n")
cat("Excluding Taxon_X Fold Change (Treated / Untreated):", fc_excl_taxon_X, "배\n")

#________________________________________________


score <- c(0.20, 0.60)
delta_score <- score[2] - score[1] # 0.40

OR <- c(estimate = 0.47, lower = 0.22, upper = 0.97)
CASI <- c(estimate = -0.92, lower = -1.51, upper = -0.33)
p0 <- 0.40


# 1. Delta score = 0.40 일 때의 Odds Ratio (OR^0.40) 계산
OR_delta <- OR ^ delta_score

# 2. Delta score = 0.40 일 때의 CASI 차이 (-0.92 * 0.40) 계산
CASI_delta <- CASI * delta_score

# 3. 새로운 증상 확률 (p1) 계산
# initial odds = p0 / (1 - p0)
# new odds = initial odds * OR_delta$estimate
# p1 = new odds / (1 + new odds)
initial_odds <- p0 / (1 - p0)
new_odds <- initial_odds * OR_delta["estimate"]
p1 <- new_odds / (1 + new_odds)

# 결과 출력
cat("=== (b) 계산 결과 ===\n")
cat("1. OR for 0.40 unit increase:\n")
print(round(OR_delta, 3))

cat("\n2. CASI difference for 0.40 unit increase:\n")
print(round(CASI_delta, 3))

cat("\n3. New probability (p1):\n")
cat("Initial p0 =", p0, "--> New p1 =", round(p1, 4), "(or", round(p1 * 100, 2), "%)\n")

