test <- merge(as.data.frame(bobcat.df), as.data.frame(bbs.df) %>% mutate(Year = 1976:(1975+nrow(bbs.df))))
test <- merge(test, as.data.frame(awssi.df) %>%
                dplyr::rename(AWSSI.E = Eastern, AWSSI.W = Western) %>%mutate(Year = 1976:(1975+nrow(awssi.df))))
test <- merge(test, as.data.frame(econ_data) %>%
                mutate(Year = 1976:(1975+nrow(econ_data))))
test <- merge(test, as.data.frame(pdsi_df) %>%
                dplyr::rename(PDSI.E = Eastern, PDSI.W = Western) %>% mutate(Year = 1976:(1975+nrow(pdsi_df))))

cortest <- test %>% select(-Year) %>% na.omit()
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
# matrix of the p-value of the correlation
p.mat <- cor.mtest(cortest)
head(p.mat[, 1:5])


cor(cortest)
corrplot::corrplot(cor(cortest), method="color", type = "upper", order = "hclust",
                   p.mat = p.mat, sig.level = 0.01)
