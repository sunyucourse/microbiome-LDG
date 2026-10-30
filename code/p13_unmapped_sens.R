#!/usr/bin/env Rscript
# P1-3b: does adding grid-level unmapped fraction change the latitude gradient?
suppressMessages(library(mgcv))
d <- read.csv("/workspace/data/processed/p13_grid_unmapped.csv")
gt <- read.csv("/mnt/results/tables/grid_table.csv")
m <- merge(gt, d[, c("habitat_groups", "grid_id", "unmapped")],
           by.x = c("habitat", "grid_id"), by.y = c("habitat_groups", "grid_id"))
m <- m[!is.na(m$richness) & m$richness > 0 & !is.na(m$unmapped), ]
m$log_rich <- log(m$richness); m$log_effort <- log(m$n_samples)

out <- list()
for (hab in sort(unique(m$habitat))) {
  x <- m[m$habitat == hab, ]
  if (nrow(x) < 30) next
  k <- max(4, min(10, floor(nrow(x) / 30)))
  g0 <- gam(log_rich ~ s(lat_abs, bs = "tp", k = k) + s(log_effort, bs = "tp", k = 4) + tech_frac_wgs,
            data = x, method = "REML")
  g1 <- gam(log_rich ~ s(lat_abs, bs = "tp", k = k) + s(log_effort, bs = "tp", k = 4) + tech_frac_wgs +
              s(unmapped, bs = "tp", k = 4), data = x, method = "REML")
  nd <- data.frame(lat_abs = c(10, 50), log_effort = median(x$log_effort),
                   tech_frac_wgs = median(x$tech_frac_wgs), unmapped = median(x$unmapped))
  r0 <- as.numeric(exp(predict(g0, newdata = nd)[1] - predict(g0, newdata = nd)[2]))
  r1 <- as.numeric(exp(predict(g1, newdata = nd)[1] - predict(g1, newdata = nd)[2]))
  out[[hab]] <- data.frame(habitat = hab, n = nrow(x),
    ratio_no_unmapped = round(r0, 2), ratio_with_unmapped = round(r1, 2),
    dev0 = round(summary(g0)$dev.expl, 3), dev1 = round(summary(g1)$dev.expl, 3))
}
res <- do.call(rbind, out)
write.csv(res, "/workspace/data/processed/p13_unmapped_sensitivity.csv", row.names = FALSE)
print(res)
