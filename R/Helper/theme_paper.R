# ============================================================
# theme_paper.R
# Unified ggplot2 theme + color palette for the TNBC exosomal
# proteomics overview mega-figure (AHP-CDS paper)
#
# Usage:
#   source("theme_paper.R")
#   ggplot(df, aes(x, y)) + geom_point() + theme_paper()
#
# Designed at 2x scale for editing (base_size = 16). When you
# assemble the mega-figure in BioRender, each sub-plot should be
# exported large and then scaled down together — never mix panels
# drawn at different base_size, or the final figure will have
# mismatched font weight when everything is reduced to print size.
# ============================================================

library(ggplot2)

# ---- Base theme ---------------------------------------------
# OUP/Bioinformatics convention: panel labels A, B, C (uppercase),
# placed top-left of each panel in BioRender, same font/size
# throughout. Sans-serif (Arial/Helvetica) reads reliably down to
# ~6-7pt final print size — keep this consistent across every
# sub-plot regardless of which R package generated it.

theme_paper <- function(base_size = 16, base_family = "Arial") {
  theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # Panel
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
      panel.border     = element_blank(),
      axis.line        = element_line(color = "grey20", linewidth = 0.3),
      axis.ticks       = element_line(color = "grey20", linewidth = 0.3),
      
      # Text
      plot.title   = element_text(size = rel(1), face = "bold", hjust = 0, margin = margin(b = 4)),
      axis.title   = element_text(size = rel(0.9), color = "grey20"),
      axis.text    = element_text(size = rel(0.8), color = "grey30"),
      legend.title = element_text(size = rel(0.85)),
      legend.text  = element_text(size = rel(0.75)),
      
      # Legend — bottom, compact (mega-figures run out of horizontal
      # room fast once panels sit side by side)
      legend.position = "bottom",
      legend.key.size = unit(0.35, "cm"),
      legend.margin   = margin(t = 2),
      
      # Facet strips (GO BP/CC/MF panel, tier breakdowns, etc.)
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text       = element_text(size = rel(0.85), face = "bold", color = "grey20"),
      
      plot.margin = margin(4, 4, 4, 4)
    )
}

# ---- Color palettes (Okabe-Ito base — colorblind-safe, and the
# same 8-color set reproduces reasonably in grayscale, which
# matters if a reviewer prints the PDF) --------------------------

# Method comparators: reused across 5.1 rank plots, 5.3 portfolio
# plots — AHP-CDS should always be the same color everywhere it
# appears so a reader tracks it visually across panels.
pal_method <- c(
  "AHP-CDS"         = "#D55E00",  # your method — reserve this color for AHP-CDS only
  "FC-only"         = "#0072B2",
  "Centrality-only" = "#009E73",
  "Equal-weight"    = "#CC79A7",
  "RWR"             = "#56B4E9",
  "Null"            = "grey60"
)

# CTD tiers: 5.3 portfolio evaluation, 5.4 case study annotations. Tier 2
# and 3 share the same downstream meaning (scientific basis to keep
# investigating the protein further), so they stay one bucket here.
pal_tier <- c(
  "Tier 1"   = "#D55E00",
  "Tier 2/3" = "#E69F00",
  "Tier 4"   = "grey70"
)

# STRING network panel only: Tier 3's definition is itself network-based,
# so this panel needs Tier 2 and 3 split apart (unlike pal_tier above)
# to actually show how Tier 3 connects to Tier 1. Same ordinal
# red->orange->yellow->grey progression, just one more step.
pal_tier_string <- c(
  "Tier 1" = "#D55E00",
  "Tier 2" = "#E69F00",
  "Tier 3" = "#F0E442",
  "Tier 4" = "grey70"
)

# AHP criteria: 5.2 contribution-share and LOCO plots
pal_criterion <- c(
  "FC"  = "#D55E00",
  "FDR" = "#0072B2",
  "Deg" = "#009E73",
  "Bet" = "#CC79A7"
)

# GO ontology categories: supplementary GO enrichment panel
pal_go <- c(
  "Biological Process" = "#D55E00",
  "Cellular Component" = "#009E73",
  "Molecular Function"  = "#0072B2"
)

# Volcano / significance coloring: limma DE panel
pal_sig <- c("Up" = "#D55E00", "Down" = "#0072B2", "NS" = "grey70")

# Sample condition (missingness boxplot, PCA/QC panels): Okabe-Ito
# blue/orange — the standard two-condition pair from that set.
# Deliberately not vermillion (#D55E00), since that's reserved for
# AHP-CDS elsewhere in the figure.
pal_condition <- c(
  "MCF10A (Normal)"    = "#0072B2",
  "MDA-MB-231 (Tumor)" = "#E69F00"
)

# Volcano category (actual scheme used): NOT direction-of-change like
# pal_sig — this axis is "does this protein survive the pipeline's
# scoping filters", so vermillion here is the SAME referent as
# AHP-CDS's highlight color elsewhere (Candidate = the 143-protein
# pool that AHP-CDS ranks), not a reuse-collision.
pal_candidate <- c(
  "Non-significant" = "grey70",
  "DE only"          = "#0072B2",
  "Candidate"        = "#D55E00"
)

# ---- STRING network panel (igraph + ggraph) --------------------
# ggraph plots are ggplot objects but theme_paper() alone won't
# style node/edge geoms — set those explicitly so the network
# panel matches the rest of the figure.
#
# library(ggraph); library(igraph)
#
# ggraph(g, layout = "fr") +
#   geom_edge_link(color = "grey70", width = 0.3, alpha = 0.6) +
#   geom_node_point(aes(fill = tier), shape = 21, color = "white",
#                    stroke = 0.3, size = 3) +
#   scale_fill_manual(values = pal_tier) +
#   geom_node_text(aes(label = label), size = 2.2, repel = TRUE,
#                   family = "Arial", color = "grey20") +
#   theme_graph_paper()

theme_graph_paper <- function(base_size = 16, base_family = "Arial") {
  ggraph::theme_graph(base_size = base_size, base_family = base_family) +
    theme(
      legend.position = "bottom",
      legend.title     = element_text(size = rel(0.85)),
      legend.text      = element_text(size = rel(0.75)),
      legend.key.size  = unit(0.35, "cm"),
      plot.margin      = margin(4, 4, 4, 4)
    )
}
# Case-study nodes (COL18A1, GSTP1, VCL, CASP14, HLA-DRA) — don't
# give them a 7th ad hoc color. Keep the pal_tier fill and instead
# make them visually pop via geom_node_point(size =) larger, or a
# thicker white stroke, so the tier-color logic stays legible
# across the whole figure.

# ---- Panel labels (BioRender assembly) --------------------------
# Confirmed: uppercase A, B, C convention (OUP/Bioinformatics norm).
# Apply as a text layer in BioRender, not inside each R plot —
# keeps label font/size/position 100% consistent regardless of
# which panel came from ggplot2 vs ggraph vs a non-R source (e.g.
# a STRING-website export), since BioRender is the single place
# controlling that layer.

# ---- Notes on applying this across the pipeline figure --------
# - Missingness panel (sigmoid + boxplot): use grey60/grey20, no
#   categorical palette needed — it's not a comparator plot.
# - STRING network: edge stroke 0.3-0.5pt, grey70; highlight case
#   study nodes (COL18A1, GSTP1, VCL, CASP14, HLA-DRA) using
#   pal_tier by their CTD tier, not a separate ad hoc color.
# - Reactome GSEA bar/dot plot: use pal_sig-style two-tone
#   (enriched vs not) rather than a rainbow gradient.
# - Keep line width consistent at 0.3-0.5pt across every panel —
#   mixing thicker/thinner strokes between R-generated panels is
#   the most common tell that a mega-figure was assembled from
#   inconsistent sources.