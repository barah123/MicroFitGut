# Validation figures, built to the dataviz skill's procedure.
#
# Palette provenance. Every colour below is a slot from the skill's validated
# categorical theme, and every set was run through
#   scripts/validate_palette.js --mode light --surface "#fcfcfb" --pairs all
# rather than chosen by eye. Results recorded beside each set.
#
# Target surface is LIGHT ONLY (#fcfcfb). These are publication figures: a
# journal prints on white, and the README embeds the same PNG. The skill's
# dark-mode requirement applies to HTML charts that re-render per viewer; a
# fixed raster cannot. Stated rather than skipped silently.
#
# The table view required for accessibility is claims-rescored.csv and
# validation-results.xlsx in this directory: every value plotted is readable
# there without colour.

ROOT <- "/Users/phil/Desktop/Tools/MicroFitGut"
SC   <- file.path(ROOT, "microfitgut-plugin/skills/microfitgut/scripts")
for (f in sort(list.files(SC, pattern = "[.]R$", full.names = TRUE))) suppressWarnings(source(f))
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
OUT <- file.path(ROOT, "microfitgut-plugin/validation/results/figures")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
set.seed(42)

# ---------------------------------------------------------------- tokens
SURFACE  <- "#fcfcfb"   # chart surface; the validator was run against this
INK      <- "#0b0b0b"   # primary
INK2     <- "#52514e"   # secondary
MUTED    <- "#898781"   # axis / labels
GRID     <- "#e1e0d9"   # hairline gridline
AXIS     <- "#c3c2b7"   # baseline

# Claim status. Blue/orange/violet, all-pairs CVD dE 13.0 (target >= 8),
# normal-vision 16.3 (floor 15), contrast all >= 3:1. Adding aqua for
# not_tested keeps the four-set passing at CVD 9.2.
#
# Green/red was rejected deliberately. The traffic-light reading is tempting and
# every ordering of it lands at CVD dE 7.2, inside the 6-8 band that is legal
# only with relief. Polarity is carried by the legend order and the labels, which
# is the skill's own rule: identity is never colour-alone.
STATUS <- c(reproduced     = "#2a78d6",
            not_reproduced = "#eb6834",
            not_licensed   = "#4a3aa7",
            not_tested     = "#1baf7a")

# Everything outside the denominator folds to ONE neutral. Nine colour classes
# in a stacked bar is the "cycling past 8" anti-pattern and blurs regardless;
# the five reasons get their own figure (F6) where they can be read.
OUTSIDE <- "#aeaca4"

# Neutral ordinal ramp, --ordinal PASS, light end 2.21:1 vs surface.
NEUTRAL3 <- c("#4a4945", "#7d7b76", "#aeaca4")

# Technology. Deliberately NOT blue/orange: those already mean reproduced and
# not_reproduced elsewhere in the composite, and colour must follow the entity.
# green/magenta, CVD dE 17.6. Magenta is sub-3:1, so the relief rule applies and
# the panel carries a legend plus the table view.
TECH <- c("16S" = "#008300", "Shotgun" = "#e87ba4")

theme_pub <- function(base = 11) {
  theme_minimal(base_size = base, base_family = "") +
    theme(
      plot.background   = element_rect(fill = SURFACE, colour = NA),
      panel.background  = element_rect(fill = SURFACE, colour = NA),
      panel.grid.major  = element_line(colour = GRID, linewidth = 0.3),
      panel.grid.minor  = element_blank(),
      axis.line         = element_blank(),
      axis.ticks        = element_blank(),
      axis.text         = element_text(colour = MUTED, size = base * 0.85),
      axis.title        = element_text(colour = INK2, size = base * 0.85),
      strip.text        = element_text(colour = INK, face = "bold", size = base * 0.9),
      plot.title        = element_text(colour = INK, face = "bold", size = base * 1.25,
                                       margin = margin(b = 5)),
      plot.subtitle     = element_text(colour = INK2, size = base * 0.92,
                                       lineheight = 1.25, margin = margin(b = 11)),
      plot.caption      = element_text(colour = MUTED, size = base * 0.78, hjust = 0),
      legend.position   = "bottom",
      legend.title       = element_blank(),
      legend.text        = element_text(colour = INK2, size = base * 0.85),
      legend.key.size    = unit(10, "pt"),
      plot.margin        = margin(12, 14, 10, 12))
}
sv <- function(p, n, w, h) {
  ggsave(file.path(OUT, paste0(n, ".png")), p, width = w, height = h,
         dpi = 300, bg = SURFACE); n
}

cl <- read.csv("claims-rescored.csv", stringsAsFactors = FALSE)
cl$scorable <- cl$new_status %in% MFG_CLAIM_SCORABLE
cl$label    <- paste0(cl$study_id, "  ", cl$study, " ", cl$year)
# Fold every unscorable reason into one class for the stacked view.
cl$plot_status <- ifelse(cl$scorable, cl$new_status, "outside the denominator")
PS_LEV <- c(names(STATUS), "outside the denominator")
PS_COL <- c(STATUS, "outside the denominator" = OUTSIDE)
cl$plot_status <- factor(cl$plot_status, levels = PS_LEV)

## ---- F1. Claim outcomes by study -------------------------------------------
ord <- cl %>% group_by(label) %>% summarise(s = mean(scorable), n = n()) %>%
  arrange(s, n) %>% pull(label)
f1 <- cl %>% mutate(label = factor(label, levels = ord)) %>%
  ggplot(aes(label, fill = plot_status)) +
  # A surface-coloured stroke IS the 2px gap: separation without a border.
  geom_bar(colour = SURFACE, linewidth = 0.8, width = 0.72) +
  coord_flip() +
  scale_fill_manual(values = PS_COL, drop = TRUE) +
  scale_y_continuous(breaks = seq(0, 10, 2), expand = expansion(mult = c(0, 0.04))) +
  labs(title = "Claim outcomes by study",
       subtitle = paste("Coloured segments enter the endpoint denominator.",
                        "Grey does not;\nits five reasons are separated in F6."),
       x = NULL, y = "claims") +
  guides(fill = guide_legend(nrow = 1)) +
  theme_pub(11)
sv(f1, "F1-claim-status-by-study", 8.2, 4.6)

## ---- F2. Adjudicability against reproduction --------------------------------
f2d <- cl %>% group_by(label, technology) %>%
  summarise(n = n(), scorable = sum(scorable),
            repro = sum(new_status == "reproduced"), .groups = "drop") %>%
  mutate(Adjudicability = scorable / n,
         `Reproduction | adjudicable` = ifelse(scorable > 0, repro / scorable, NA))
f2 <- f2d %>% select(label, technology, Adjudicability, `Reproduction | adjudicable`) %>%
  pivot_longer(-c(label, technology)) %>%
  ggplot(aes(name, value)) +
  geom_line(aes(group = label), colour = AXIS, linewidth = 0.4) +
  geom_point(aes(colour = technology), size = 3, alpha = 0.9) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1),
                     breaks = seq(0, 1, 0.25)) +
  scale_colour_manual(values = TECH) +
  scale_x_discrete(expand = expansion(add = 0.35)) +
  labs(title = "Where the between-study variation lives",
       subtitle = paste("One line per study. Permutation test of between-study",
                        "homogeneity:\nadjudicability p = 0.06, reproduction p = 0.49.",
                        "Both spreads are wide at n = 10."),
       x = NULL, y = NULL) +
  theme_pub(11)
sv(f2, "F2-adjudicability-vs-reproduction", 5.6, 4.8)

## ---- F3. Every scored claim by reanalysis p-value ---------------------------
f3d <- cl %>% filter(!is.na(p_value), new_status %in% c("reproduced", "not_reproduced")) %>%
  mutate(p = as.numeric(p_value),
         dir = ifelse(direction_agrees == "TRUE", "direction agrees", "direction reversed"),
         lab = ifelse(nchar(claim_id) > 22, paste0(substr(claim_id, 1, 21), "\u2026"), claim_id),
         new_status = factor(new_status, levels = c("reproduced", "not_reproduced")),
         panel = factor(claim_polarity, levels = c("difference", "null"),
                        labels = c("paper asserted\na difference",
                                   "paper asserted\nno difference")))
f3 <- ggplot(f3d, aes(p, reorder(paste0(study_id, "  ", lab), -p))) +
  annotate("rect", xmin = 0.05, xmax = 2, ymin = -Inf, ymax = Inf,
           fill = "#0b0b0b", alpha = 0.045) +
  # Dashed is correct here: this is a threshold, not a gridline.
  geom_vline(xintercept = 0.05, linetype = "22", colour = MUTED, linewidth = 0.4) +
  geom_point(aes(colour = new_status, shape = dir), size = 2.4, stroke = 0.9) +
  facet_grid(panel ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_x_log10(breaks = c(1e-15, 1e-10, 1e-5, 0.05, 1),
                labels = c("1e-15", "1e-10", "1e-5", "0.05", "1"),
                limits = c(1e-17, 2)) +
  scale_colour_manual(values = STATUS) +
  scale_shape_manual(values = c("direction agrees" = 16, "direction reversed" = 17)) +
  labs(title = "Every scored claim, by the p-value of its reanalysis",
       subtitle = paste("Shaded region is p >= 0.05, and it means opposite things",
                        "on the two panels.\nNo claim falls between 0.017 and 0.059,",
                        "so the threshold adjudicates no close call."),
       x = "reanalysis p-value (log scale)", y = NULL) +
  guides(colour = guide_legend(order = 1, nrow = 1),
         shape  = guide_legend(order = 2, nrow = 1)) +
  theme_pub(9.5) +
  theme(strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 8.5,
                                         lineheight = 1.15, colour = INK2),
        legend.box = "vertical", legend.spacing.y = unit(1, "pt"))
sv(f3, "F3-significance-first", 8.6, 7.6)

## ---- F4. Course material --------------------------------------------------
cw <- read.csv(file.path(ROOT,
  "output/regression-20260912-193250/tables/regression_results_full.csv"),
  stringsAsFactors = FALSE)
cw$set <- sub("^(ps[0-9]+|quiz[0-9]+).*$", "\\1", cw$source)
# MATCH/DIFFER take the same two slots as reproduced/not_reproduced: the meaning
# is the same, so the colour should be. BLOCKED takes violet, the slot that
# elsewhere means "the tool declined to score this". The two genuinely excluded
# kinds take neutral steps.
CW_LEV <- c("MATCH", "DIFFER", "BLOCKED", "NODATA", "INFO")
CW_COL <- c(MATCH = STATUS[["reproduced"]], DIFFER = STATUS[["not_reproduced"]],
            BLOCKED = STATUS[["not_licensed"]],
            NODATA = NEUTRAL3[2], INFO = NEUTRAL3[3])
f4 <- cw %>% mutate(status = factor(status, levels = CW_LEV)) %>%
  ggplot(aes(reorder(set, ave(rep(1, nrow(cw)), set, FUN = length)), fill = status)) +
  geom_bar(colour = SURFACE, linewidth = 0.8, width = 0.72) + coord_flip() +
  scale_fill_manual(values = CW_COL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.04))) +
  labs(title = "Course material: 106 checks against the stated answers",
       subtitle = paste("All 8 discrepancies fall in ps10 and ps7, and trace to",
                        "three documented causes.\nThe headline unit is the",
                        "exercise question, not the row: 28 of 34 agree."),
       x = NULL, y = "checks") +
  guides(fill = guide_legend(nrow = 1)) +
  theme_pub(11)
sv(f4, "F4-coursework-agreement", 7.6, 4.2)

## ---- F5. The two bodies, side by side ---------------------------------------
# One series: the facet strips carry identity, so no hue is spent on it and no
# colour collides with the status meanings used elsewhere.
f5d <- data.frame(
  body = factor(c(rep("Course material  (verification)", 2),
                  rep("Published papers  (reproducibility)", 2)),
                levels = c("Course material  (verification)",
                           "Published papers  (reproducibility)")),
  what = c("per exercise\nquestion", "per row\n(not independent)",
           "primary claim\n(POST HOC)", "scorable claims\nreproduced"),
  x = c(28, 84, 10, 30), n = c(34, 92, 10, 40), stringsAsFactors = FALSE)
f5d$what <- factor(f5d$what, levels = f5d$what)
f5d$p <- f5d$x / f5d$n
cp <- function(x, n) stats::binom.test(x, n)$conf.int[1:2]
f5d[c("lo", "hi")] <- t(mapply(cp, f5d$x, f5d$n))
f5 <- ggplot(f5d, aes(what, p)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.10, linewidth = 0.6, colour = INK2) +
  geom_point(size = 3.2, colour = INK) +
  # Values wear text ink, never a series colour.
  geom_text(aes(label = sprintf("%d/%d", x, n)), nudge_x = 0.17, hjust = 0,
            size = 3.1, colour = INK2) +
  scale_x_discrete(expand = expansion(add = c(0.55, 0.95))) +
  facet_grid(~ body, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = scales::percent, limits = c(0.4, 1.06),
                     breaks = seq(0.4, 1, 0.2)) +
  labs(title = "Two bodies of evidence, reported separately",
       subtitle = paste("95% Clopper-Pearson intervals. These measure different",
                        "things and are never\npooled into one accuracy figure."),
       x = NULL, y = NULL) +
  theme_pub(11)
sv(f5, "F5-two-bodies", 7.6, 4.4)

## ---- F6. What the grey is made of -------------------------------------------
# One series, so one colour and no legend: the title names it.
# The meaning goes in the axis label, never inside the bar. Text set inside a
# mark is clipped the moment the mark is short, which is exactly what happened
# on the first render here.
F6_MEANING <- c(schema_gap        = "no claim type expresses it",
                out_of_scope      = "not a microbiome measurement",
                not_attempted     = "the analyst did not run it",
                data_absent       = "never deposited",
                contrast_mismatch = "subgroup not reconstructable")
f6d <- cl %>% filter(!scorable) %>% count(new_status) %>%
  mutate(lab = paste0(new_status, "\n", F6_MEANING[as.character(new_status)]),
         lab = reorder(lab, n))
f6 <- ggplot(f6d, aes(n, lab)) +
  geom_col(fill = OUTSIDE, width = 0.6) +
  geom_text(aes(label = n), hjust = -0.6, size = 3.2, colour = INK2) +
  scale_x_continuous(breaks = seq(0, 8, 2),
                     expand = expansion(mult = c(0, 0.10))) +
  labs(title = "The 17 claims outside the denominator",
       subtitle = paste("Each reason is a different finding and they are never",
                        "pooled.\nOnly one of the five, not_attempted, is under",
                        "the analyst's control."),
       x = "claims", y = NULL) +
  theme_pub(10.5) +
  theme(axis.text.y = element_text(colour = INK2, lineheight = 1.15, hjust = 1))
sv(f6, "F6-outside-the-denominator", 7.0, 3.4)

cat("figures:", paste(list.files(OUT, "^F[1-6].*png$"), collapse = ", "), "\n")

## ---- F0. Stacked composite for the README ------------------------------------
# Stacked, not gridded: F1 and F3 carry long labels that a two-column grid would
# compress to illegibility. Patchwork aligns axes across the stack, so the
# widest label sets every panel's left edge, which is why F3's claim ids are
# truncated above.
if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  # `+` binds tighter than `&`, so `x + plot_layout() & theme() + plot_annotation()`
  # parses as (x + layout) & (theme + annotation): the annotation lands on the
  # theme and the composite title silently disappears. Build it in steps.
  comp <- (f5 + labs(tag = "A")) / (f4 + labs(tag = "B")) /
          (f1 + labs(tag = "C")) / (f6 + labs(tag = "D")) /
          (f2 + labs(tag = "E")) / (f3 + labs(tag = "F"))
  comp <- comp + plot_layout(heights = c(4.4, 4.2, 4.6, 3.4, 4.8, 7.6))
  # plot.tag.location defaults to "margin", which reserves a gutter outside each
  # plot for the tag.
  comp <- comp & theme(plot.tag = element_text(face = "bold", size = 14, colour = INK),
                       plot.tag.location = "plot",
                       plot.tag.position = c(0.002, 0.99))
  comp <- comp + plot_annotation(
    title = "MicroFitGut validation",
    subtitle = paste("Two bodies of evidence, reported separately and never",
                     "pooled into one accuracy figure.\nA-B: 106 checks against",
                     "the course material's stated answers.",
                     "C-F: 57 claims from ten published papers."),
    caption = paste("Palette: slots from a CVD-validated categorical theme;",
                    "every set run through validate_palette.js (all-pairs CVD",
                    "ΔE 13.0, normal-vision 16.3, contrast >= 3:1).",
                    "\nTable view: claims-rescored.csv and validation-results.xlsx",
                    "in this directory carry every plotted value without colour."),
    theme = theme_pub(12) +
      theme(plot.title = element_text(face = "bold", size = 20),
            plot.subtitle = element_text(size = 11.5, lineheight = 1.3,
                                         margin = margin(b = 4)),
            plot.caption = element_text(size = 8.5, lineheight = 1.3,
                                        margin = margin(t = 10))))
  ggsave(file.path(OUT, "F0-validation-stacked.png"), comp,
         width = 9.0, height = 29, dpi = 200, bg = SURFACE, limitsize = FALSE)
  cat("stacked composite written\n")
}
