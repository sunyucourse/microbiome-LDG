#!/usr/bin/env Rscript
# P0-3: per-habitat GAM with vs without study random effect (dominant study per grid).
suppressMessages(library(mgcv)); suppressMessages(library(MASS))

df <- read.csv("/workspace/data/processed/p03_grid_studies.csv")
gt <- read.csv("/mnt/results/tables/grid_table.csv")
d <- merge(gt, df, by = c("habitat", "grid_id"))
d <- d[!is.na(d$richness) & d$richness > 0, ]
d$log_rich <- log(d$richness)
d$log_effort <- log(d$n_samples)

set.seed(1)
out <- list()
for (hab in sort(unique(d$habitat))) {
  x <- d[d$habitat == hab, ]
  if (nrow(x) < 30 || length(unique(x$dom_study)) < 5) next
  k_lat <- max(4, min(10, floor(nrow(x) / 30)))
  f0 <- log_rich ~ s(lat_abs, bs = "tp", k = k_lat) + s(log_effort, bs = "tp", k = 4) + tech_frac_wgs
  f1 <- log_rich ~ s(lat_abs, bs = "tp", k = k_lat) + s(log_effort, bs = "tp", k = 4) + tech_frac_wgs +
         s(dom_study, bs = "re")
  g0 <- gam(f0, data = x, method = "REML")
  x$dom_study <- factor(x$dom_study)
  g1 <- tryCatch(gam(f1, data = x, method = "REML"), error = function(e) { cat(hab, "ERR:", conditionMessage(e), "\n"); NULL })
  # predicted ratio 10v50 from each model (median covariates)
  nd <- data.frame(lat_abs = c(10, 50), log_effort = median(x$log_effort),
                   tech_frac_wgs = median(x$tech_frac_wgs))
  r0 <- exp(diff(predict(g0, newdata = nd)))  # exp(fit10 - fit50) -> careful order
  r0 <- as.numeric(exp(predict(g0, newdata = nd)[1] - predict(g0, newdata = nd)[2]))
  r1 <- NA
  if (!is.null(g1)) {
    nd$dom_study <- names(which.max(table(x$dom_study)))
    r1 <- as.numeric(exp(predict(g1, newdata = nd)[1] - predict(g1, newdata = nd)[2]))
    edf_study <- sum(g1$edf[-(1:(length(g1$edf) - 1))])  # last smooth = study RE
    pval_study <- summary(g1)$s.table[nrow(summary(g1)$s.table), 4]
  } else { edf_study <- NA; pval_study <- NA }
  out[[hab]] <- data.frame(habitat = hab, n_grids = nrow(x), n_studies = length(unique(x$dom_study)),
    dev0 = round(summary(g0)$dev.expl, 3),
    ratio0 = round(r0, 2), ratio1 = ifelse(is.na(r1), NA, round(r1, 2)),
    dev1 = ifelse(is.null(g1), NA, round(summary(g1)$dev.expl, 3)),
    study_pval = signif(pval_study, 3))
}
res <- do.call(rbind, out)
write.csv(res, "/workspace/data/processed/p03_study_re.csv", row.names = FALSE)
print(res)
