setwd("/mnt/picea/projects/aspseq/nstreet/kmer_ml/sex/selected_kmers/")
library(tidyverse)
library(magrittr)
library(ranger)
library(caret)
library(GenomicAlignments)
library(see)
library(UpSetR)
counts <- scan("count_data.txt", "")


ss1 <- read_delim("../meta/sex_set_metadata.tsv", "\t")
ss2 <- read_delim("../meta/scotasp.txt", "\t")

ss1 %<>% mutate(Sex = na_if(Sex, "excluded"))
unique(ss1$Sex)
unique(ss2$PCR_Sex)
ss1$SampleSet <- "SwAsp_UmAsp"

ss2 %<>% select(NGI_ID, PCR_Sex)
colnames(ss2) <- c("Sample", "Sex")
ss2$SampleSet <- "ScotAsp"
meta <- bind_rows(ss1, ss2)
# 

meta$ShortName <- sapply(str_split(meta$Sample, "_"), function(x) {
  str_replace(paste(x[1], x[2], x[3], sep = "_"), "(_1.jf)|(_NA)", "")
})

assertthat::assert_that(! any(duplicated(meta$ShortName)))

meta$CountFile <- unlist(sapply(meta$ShortName, function(x) {
  res <- counts[which(grepl(x, counts, fixed = TRUE))]
  if (length(res) == 0) {
    return(NA)
  } else {
    return(res)
  }
}))

meta_f <- filter(meta, !is.na(Sex) & !is.na(CountFile))
meta_f$Sex <- as.factor(meta_f$Sex)
meta_f$Coverage <- str_replace(meta_f$CountFile, "\\.count", ".cov")

# Read all the data into a vector
count_matrix <- scan(
  pipe(
    paste("cat", paste(meta_f$CountFile, collapse = " ")
    )
  )
)
# Cast to matrix
count_matrix <- matrix(count_matrix, nrow = nrow(meta_f), byrow = TRUE)

coverage <- scan(
  pipe(
    paste("cat", paste(meta_f$Coverage, collapse = " ")
    )
  )
)
size_factor <- coverage / max(coverage)
count_matrix_norm <- sweep(count_matrix, 1, size_factor, "/")

# Need names for ranger()
kmers <- scan("top_k.txt", "")
colnames(count_matrix_norm) <- kmers
rownames(count_matrix_norm) <- meta_f$ShortName

train_x <- count_matrix_norm[which(meta_f$SampleSet == "SwAsp_UmAsp"), ]
test_x <- count_matrix_norm[which(meta_f$SampleSet != "SwAsp_UmAsp"), ]
train_y <- meta_f$Sex[which(meta_f$SampleSet == "SwAsp_UmAsp")]
test_y <- meta_f$Sex[which(meta_f$SampleSet != "SwAsp_UmAsp")]

head(train_x[, 1:6])

model <- ranger(x = train_x, y = train_y, num.trees = 50000,
                importance = "permutation", verbose = TRUE,
                oob.error = TRUE, classification = TRUE)
save(model, file = "model.RData")
preds <- predict(model, data = test_x, predict.all = TRUE)

weights <- rowSums(preds$predictions == 1) / 50000
outcomes <- ifelse(weights > 0.5, 0, 1)

confusionMatrix(data = preds$predictions, reference = test_y)

rox <- pROC::roc(response = as.numeric(factor(test_y)),
                 predictor = weights)
plot(rox)
tibble(Specificity = rox$specificities, Sensitivity = rox$sensitivities) %>%
ggplot(aes(x = Specificity, y = Sensitivity)) +
  geom_path()+
  scale_x_reverse() +
  geom_abline(slope = 1, intercept = 1, lty = 2, col = "red") +
  theme_bw()
ggsave("~/roc.pdf", height = 9, width = 9, dpi = 300)

model2 <- ranger(x = train_x, y = train_y, num.trees = 50000,
                importance = "impurity_corrected", verbose = TRUE,
                oob.error = TRUE, classification = TRUE)
save(model2, file = "model2.RData")
imp_pvals <- importance_pvalues(model2)

imp_pvals[imp_pvals[, 2] < 0.05, ]
hist(imp_pvals[, 2])


# alignments
potra_v2 <- readGAlignments("../alignments/Potra_v2_sorted.bam", use.names = TRUE)
mcols(potra_v2)$importance <- model$variable.importance[names(potra_v2)]
w52 <- readGAlignments("../alignments/W52_sorted.bam", use.names = TRUE)
mcols(w52)$importance <- model$variable.importance[names(w52)]

# Calculate various sets
miss_potra <- setdiff(rownames(imp_pvals), names(potra_v2))
miss_w52 <- setdiff(rownames(imp_pvals), names(w52))
aln_potra <- intersect(rownames(imp_pvals), names(potra_v2))
aln_w52 <- intersect(rownames(imp_pvals), names(w52))

miss_both <- intersect(miss_w52, miss_potra)
miss_w52_only <- setdiff(miss_w52, miss_potra)
miss_potra_only <- setdiff(miss_potra, miss_w52)
aln_both <- intersect(aln_potra, aln_w52)
aln_w52_only <- setdiff(aln_w52, aln_potra)
aln_potra_only <- setdiff(aln_potra, aln_w52)

upset(data = fromList(list("Aligned to Male" = names(w52), "Aligned to Female" = names(potra_v2), "Kmer superset" = rownames(imp_pvals))))

as_tibble(colSums(count_matrix_norm[, miss_both] > 0)) %>%
  ggplot(aes(x = value)) +
  geom_histogram(binwidth = 10)
ggsave("~/kmer_hist_missing.pdf", width = 16, height = 9)

list(
  tibble(Pvalue = imp_pvals[miss_potra, 2], Kmer = miss_potra, Genome = "asp201", Type = "Missing"),
  tibble(Pvalue = imp_pvals[miss_w52, 2], Kmer = miss_w52, Genome = "W52", Type = "Missing"),
  tibble(Pvalue = imp_pvals[aln_potra, 2], Kmer = aln_potra, Genome = "asp201", Type = "Aligned"),
  tibble(Pvalue = imp_pvals[aln_w52, 2], Kmer = aln_w52, Genome = "W52", Type = "Aligned")
) %>% bind_rows() -> imp_alns

filter(imp_alns, Pvalue < 0.1) %>%
  ggplot(aes(x = Genome, y = -log10(Pvalue + .Machine$double.eps), fill = Type)) +
  geom_violin() +
  scale_fill_material()

imp_alns %>% filter(Pvalue < 0.1) %>% group_by(Genome, Type) %>%
  summarise(N = length(Pvalue), Min = min(Pvalue), Max = max(Pvalue),
            Mean = mean(Pvalue), Median = median(Pvalue))

## Manhattan plot for Potra
potra_v2_chr <- potra_v2[str_detect(seqnames(potra_v2), "^chr")]
w52_chr <- w52[str_detect(seqnames(w52), "^chr")]


preds2 <- predict(model, data = test_x)
baseline <- confusionMatrix(data = preds2$predictions, reference = test_y)
baseline

#######
## Anything below here was experimental and not included
#######

model3 <- ranger(x = train_x[, ! colnames(train_x) %in% chr1kmers],
                 y = train_y, num.trees = 50000,
                verbose = TRUE, oob.error = TRUE, classification = TRUE)
save(model3, file = "model3.RData")
preds3 <- predict(model3, data = test_x)
no_chr_1 <- confusionMatrix(data = preds3$predictions, reference = test_y)


model4 <- ranger(x = train_x[, ! colnames(train_x) %in% chr19kmers],
                 y = train_y, num.trees = 50000,
                 verbose = TRUE, oob.error = TRUE, classification = TRUE)
save(model4, file = "model4.RData")
preds4 <- predict(model4, data = test_x)
no_chr_19 <- confusionMatrix(data = preds4$predictions, reference = test_y)


model5 <- ranger(x = train_x[, ! colnames(train_x) %in% miss_both],
                 y = train_y, num.trees = 50000,
                 verbose = TRUE, oob.error = TRUE, classification = TRUE)
save(model5, file = "model5.RData")
preds5 <- predict(model5, data = test_x)
no_mb <- confusionMatrix(data = preds4$predictions, reference = test_y)


model6 <- ranger(x = train_x[, ! colnames(train_x) %in% c(miss_both, chr19kmers, chr1kmers)],
                 y = train_y, num.trees = 50000,
                 verbose = TRUE, oob.error = TRUE, classification = TRUE)
save(model6, file = "model6.RData")
preds6 <- predict(model6, data = test_x)
no_mb <- confusionMatrix(data = preds4$predictions, reference = test_y)

pred_tbl <- tibble(
  Sample = meta_f$Sample[which(meta_f$SampleSet != "SwAsp_UmAsp")],
  Predicted = as.character(preds6$predictions),
  Actual = test_y,
  Match = as.character(preds6$predictions) == test_y
)
writexl::write_xlsx(pred_tbl, "~/test_set_predictions.xlsx")


as_tibble(imp_pvals, rownames = "Kmer") -> imps

as_tibble(count_matrix_norm[, head(arrange(imps, desc(importance))$Kmer, 16)]) -> tmp
tmp$Sex <- meta_f$Sex
tmp$Sample <- meta_f$ShortName
pivot_longer(tmp, -c(Sample, Sex)) -> plot_df

left_join(plot_df, pred_tbl) %>%
  ggplot(aes(x = Sex, y = value + .Machine$double.eps, label = Sample)) +
  geom_violin(outlier.alpha = 0) +
  geom_jitter(aes(col = Match), width = 0.2) +
  facet_wrap(~name, scales = "free_y") +
  theme_bw() +
  scale_y_log10()

plotly::ggplotly()
