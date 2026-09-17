ROOT <- "/Users/phil/Desktop/Tools/MicroFitGut"
SC   <- file.path(ROOT, "microfitgut-plugin/skills/microfitgut/scripts")
for (f in sort(list.files(SC, pattern = "[.]R$", full.names = TRUE))) suppressWarnings(source(f))
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
OUT <- file.path(ROOT, "microfitgut-plugin/validation/results/figures")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
sv <- function(p, n, w, h) { ggsave(file.path(OUT, paste0(n, ".png")), p,
                                    width = w, height = h, dpi = 300, bg = "white"); n }

# Status colors: green reproduced, red contradicted, amber not licensed,
# greys for everything outside the denominator. Okabe-Ito derived, CVD-safe.
# Scorable statuses carry hue. Unscorable statuses are grey, in five evenly
# separated lightness steps so they stay distinguishable from one another: they
# are different findings, and a legend that cannot be read pools them by accident.
SCOL <- c(reproduced = "#009E73", not_reproduced = "#D55E00",
          not_licensed = "#E69F00", not_tested = "#F0E442",
          data_absent = "#3F3F3F", contrast_mismatch = "#6B6B6B",
          schema_gap = "#979797", not_attempted = "#C3C3C3",
          out_of_scope = "#E2E2E2", pending_adjudication = "#F2F2F2")

cl <- read.csv("claims-rescored.csv", stringsAsFactors = FALSE)
cl$new_status <- factor(cl$new_status, levels = MFG_CLAIM_STATUSES)
cl$scorable   <- cl$new_status %in% MFG_CLAIM_SCORABLE
cl$label <- paste0(cl$study_id, " ", cl$study, " ", cl$year)

## F1. Claim status by study, ordered by adjudicability.
ord <- cl %>% group_by(label) %>% summarise(s = mean(scorable), n = n()) %>%
  arrange(s, n) %>% pull(label)
f1 <- cl %>% mutate(label = factor(label, levels = ord)) %>%
  ggplot(aes(label, fill = new_status)) +
  geom_bar(color = "grey30", linewidth = 0.25) + coord_flip() +
  scale_fill_manual(values = SCOL, drop = TRUE, name = NULL) +
  labs(title = "Claim outcomes by study",
       subtitle = "Coloured statuses enter the endpoint denominator; greys do not",
       x = NULL, y = "claims") +
  theme_canis(11) + theme(legend.position = "bottom")
sv(f1, "F1-claim-status-by-study", 8.5, 5.2)

## F2. The split the pilot actually found: adjudicability varies, reproduction does not.
f2d <- cl %>% group_by(label, technology) %>%
  summarise(n = n(), scorable = sum(scorable),
            repro = sum(new_status == "reproduced"), .groups = "drop") %>%
  mutate(adjud = scorable/n, rep_cond = ifelse(scorable > 0, repro/scorable, NA))
f2 <- f2d %>% select(label, technology, Adjudicability = adjud,
                     `Reproduction | adjudicable` = rep_cond) %>%
  pivot_longer(-c(label, technology)) %>%
  ggplot(aes(name, value)) +
  geom_line(aes(group = label), color = "grey75", linewidth = 0.4) +
  geom_point(aes(color = technology), size = 3.2, alpha = 0.85) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_color_manual(values = c("16S" = "#0072B2", "Shotgun" = "#CC79A7"), name = NULL) +
  labs(title = "Where the between-study variation lives",
       subtitle = paste("Each line is one study. Permutation test of between-study",
                        "homogeneity:\nadjudicability p = 0.06, reproduction p = 0.49.",
                        "Both spreads are wide at n = 10."),
       x = NULL, y = NULL) +
  theme_canis(12) + theme(legend.position = "bottom")
sv(f2, "F2-adjudicability-vs-reproduction", 6.2, 5.0)

## F3. Significance-first rule: which claims it moved, and why.
# Faceted by what the paper asserted. The shaded region means opposite things on
# the two panels, which is exactly why they cannot share one.
f3d <- cl %>% filter(!is.na(p_value), new_status %in% c("reproduced","not_reproduced")) %>%
  mutate(p = as.numeric(p_value),
         dir = ifelse(direction_agrees == "TRUE", "direction agrees", "direction reversed"),
         panel = factor(claim_polarity,
                        levels = c("difference", "null"),
                        labels = c("Paper asserted a difference\n(shaded = claim not supported)",
                                   "Paper asserted no difference\n(shaded = claim supported)")))
# Patchwork aligns axes across a stacked composite, so the widest label here
# sets the left edge of EVERY panel. Long claim ids cost the whole figure.
f3d$lab <- ifelse(nchar(f3d$claim_id) > 22,
                  paste0(substr(f3d$claim_id, 1, 21), "\u2026"), f3d$claim_id)
f3 <- ggplot(f3d, aes(p, reorder(paste0(study_id, ":", lab), -p))) +
  annotate("rect", xmin = 0.05, xmax = 2, ymin = -Inf, ymax = Inf,
           fill = "grey50", alpha = 0.10) +
  geom_vline(xintercept = 0.05, linetype = "22", color = "grey35") +
  geom_point(aes(color = new_status, shape = dir), size = 2.6) +
  facet_grid(panel ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_x_log10(breaks = c(1e-15, 1e-10, 1e-5, 0.05, 1),
                labels = c("1e-15", "1e-10", "1e-5", "0.05", "1"),
                limits = c(1e-17, 2)) +
  scale_color_manual(values = SCOL, name = NULL) +
  scale_shape_manual(values = c("direction agrees" = 16, "direction reversed" = 17), name = NULL) +
  labs(title = "Every scored claim, by the p-value of its reanalysis",
       subtitle = "Four claims previously scored as holding on direction alone now sit on the wrong side of 0.05.",
       x = "reanalysis p-value (log scale)", y = NULL) +
  theme_canis(10) + theme(legend.position = "bottom", legend.box = "vertical",
                          strip.placement = "outside",
                          strip.text.y.left = element_text(angle = 0, hjust = 0, size = 8))
sv(f3, "F3-significance-first", 9.0, 8.0)

## F4. Course material: agreement by problem set, with the causes named.
cw <- read.csv(file.path(ROOT,
  "output/regression-20260912-193250/tables/regression_results_full.csv"),
  stringsAsFactors = FALSE)
cw$set <- sub("^(ps[0-9]+|quiz[0-9]+).*$", "\\1", cw$source)
CCOL <- c(MATCH = "#009E73", DIFFER = "#D55E00", BLOCKED = "#E69F00",
          NODATA = "#999999", INFO = "#DDDDDD")
f4 <- cw %>% mutate(status = factor(status, levels = names(CCOL))) %>%
  ggplot(aes(reorder(set, ave(rep(1, nrow(cw)), set, FUN = length)), fill = status)) +
  geom_bar(color = "white", linewidth = 0.3) + coord_flip() +
  scale_fill_manual(values = CCOL, name = NULL) +
  labs(title = "Course material: 106 checks against the stated answers",
       subtitle = "All 8 discrepancies fall in ps10 and ps7, and trace to three documented causes",
       x = NULL, y = "checks") +
  theme_canis(11) + theme(legend.position = "bottom")
sv(f4, "F4-coursework-agreement", 7.5, 4.2)

## F5. The headline rates, side by side, WITHOUT pooling them.
# Body A is shown at both units because the row level is not the honest one:
# seven ps10 rows are the attribution for the six ps10 discrepancies, not
# independent successes, and one fitted model appears twice.
f5d <- data.frame(
  body = factor(c(rep("Course material\n(verification)", 2),
                  rep("Published papers\n(reproducibility)", 2)),
                levels = c("Course material\n(verification)",
                           "Published papers\n(reproducibility)")),
  what = c("per exercise\nquestion", "per row\n(not independent)",
           "primary claim\n(POST HOC)", "scorable claims\nreproduced"),
  x = c(28, 84, 10, 30), n = c(34, 92, 10, 40),
  stringsAsFactors = FALSE)
f5d$what <- factor(f5d$what, levels = f5d$what)
f5d$p <- f5d$x / f5d$n
# Clopper-Pearson: conservative, which is the right direction when the person
# reporting the validation also wrote the tool.
cp <- function(x, n) stats::binom.test(x, n)$conf.int[1:2]
f5d[c("lo","hi")] <- t(mapply(cp, f5d$x, f5d$n))
f5 <- ggplot(f5d, aes(what, p, color = body)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.12, linewidth = 0.7) +
  geom_point(size = 4) +
  geom_text(aes(label = sprintf("%d/%d", x, n)), nudge_x = 0.16, hjust = 0,
            size = 3.3, show.legend = FALSE) +
  scale_x_discrete(expand = expansion(add = c(0.6, 0.95))) +
  facet_grid(~ body, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = scales::percent, limits = c(0.4, 1.06)) +
  scale_color_manual(values = c("#0072B2", "#CC79A7"), guide = "none") +
  labs(title = "Two bodies of evidence, reported separately",
       subtitle = paste("95% Clopper-Pearson intervals. These measure different",
                        "things and are never\npooled into one accuracy figure."),
       x = NULL, y = NULL) +
  theme_canis(11)
sv(f5, "F5-two-bodies", 8.0, 4.8)

## F0. The stacked composite for the README: all five panels in one column.
# Stacking rather than gridding keeps every panel full width, which matters
# because F1 and F3 carry long study and claim labels that a two-column grid
# would compress to illegibility.
if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  # NOTE on operator precedence: `+` binds tighter than `&`, so writing
  #   x + plot_layout() & theme() + plot_annotation()
  # parses as (x + plot_layout()) & (theme() + plot_annotation()), which adds
  # the annotation to the THEME and silently drops the composite title. Build
  # it in separate statements instead.
  comp <- (f5 + labs(tag = "A")) /
          (f4 + labs(tag = "B")) /
          (f1 + labs(tag = "C")) /
          (f2 + labs(tag = "D")) /
          (f3 + labs(tag = "E"))
  comp <- comp + plot_layout(heights = c(1.0, 0.9, 1.15, 1.1, 1.9))
  # plot.tag.location defaults to "margin", which reserves a gutter outside the
  # plot for the tag. "plot" draws it over the plot area instead.
  comp <- comp & theme(plot.tag = element_text(face = "bold", size = 15),
                       plot.tag.location = "plot",
                       plot.tag.position = c(0.002, 0.985))
  comp <- comp + plot_annotation(
    title = "MicroFitGut validation: two bodies of evidence",
    subtitle = paste("A-B: 106 checks against the course material's stated answers.",
                     "C-E: 57 claims from ten published papers.",
                     "\nThe two are never pooled into one accuracy figure."),
    theme = theme_canis(13) +
      theme(plot.title = element_text(face = "bold", size = 17),
            plot.subtitle = element_text(size = 11, lineheight = 1.2)))
  ggsave(file.path(OUT, "F0-validation-stacked.png"), comp,
         width = 9.5, height = 30, dpi = 200, bg = "white", limitsize = FALSE)
  cat("stacked composite written\n")
} else {
  cat("patchwork not installed; stacked composite skipped\n")
}
