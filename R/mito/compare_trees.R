#!/usr/bin/env Rscript

# ==============================================================================
# FALSE KILLER WHALE MITOGENOME COMPARISON SCRIPT (FIXED FOR LAYER Z-ORDER)
# ==============================================================================

library(treeio)
library(tidytree)
library(tidyverse)
library(ggforce)

# --- USER PARAMETERS ---
file_CDS_only  <- "false_killer_whale_CDS_only.tree"
file_w_Dloop   <- "false_killer_whale_w_Dloop.tree"
support_cutoff <- 0.85  


# ==============================================================================
# STEP 1: Load BEAST Trees
# ==============================================================================
message("Reading BEAST tree files...")
tree_CDS_only <- read.beast(file = 'BEAST/xml/aln.3.CDS-rRNA_constant_orc_logMRCA/aln.3.CDS-rRNA_constant_orc_logMRCA-MCC.tree')
tree_w_Dloop  <- read.beast(file = 'BEAST/xml/aln3.CDS-rRNA-Dloop_constant_orc_logMRCA/aln3.CDS-rRNA-Dloop_constant_orc_logMRCA_combined-MCC.tree')


# ==============================================================================
# STEP 2: Generate Clade-to-Tip Maps
# ==============================================================================
message("Mapping internal nodes to tip-clade compositions...")

get_clade_map <- function(beast_tree) {
  tree_df <- as_tibble(beast_tree)
  internal_nodes <- tree_df %>% filter(is.na(label)) %>% pull(node)
  
  map_dfr(internal_nodes, function(n) {
    tips <- offspring(tree_df, n) %>% 
      filter(!is.na(label)) %>% 
      pull(label) %>% 
      sort() %>% 
      paste(collapse = ",")
    
    tibble(node = n, clade_tips = tips)
  })
}

map_CDS_only <- get_clade_map(tree_CDS_only)
map_w_Dloop  <- get_clade_map(tree_w_Dloop)


# ==============================================================================
# STEP 3: Tidy and Filter Node Metadata
# ==============================================================================
message("Wrangling and filtering node data based on support thresholds...")

# A. Extract w_Dloop nodes that CLEAR the user support cutoff
nodes_w_Dloop <- as_tibble(tree_w_Dloop) %>%
  filter(tibble::is.tibble(.)) %>%
  filter(!is.na(posterior) & posterior >= support_cutoff) %>%
  mutate(Scheme = "w_Dloop") %>%
  select(node, posterior, height, height_0.95_HPD, Scheme) %>%
  left_join(map_w_Dloop, by = "node")

# B. Extract ALL internal nodes from CDS_only
nodes_CDS_only_all <- as_tibble(tree_CDS_only) %>%
  filter(tibble::is.tibble(.)) %>%
  filter(is.na(label)) %>% 
  mutate(Scheme = "CDS_only") %>%
  select(node, posterior, height, height_0.95_HPD, Scheme) %>%
  left_join(map_CDS_only, by = "node")


# ==============================================================================
# STEP 4: Merge Trees, Handle NAs, and Set Layer Z-Ordering
# ==============================================================================
message("Merging datasets and ordering layers to prevent occlusion...")

# Setup factor level labels dynamically
label_both  <- paste0("Well-supported in both (>= ", support_cutoff, ")")
label_drop  <- paste0("Supported in D-loop, but < ", support_cutoff, " in CDS")
label_novel <- "Novel Node (Resolved by D-loop only)"

plot_data_clean <- left_join(
  nodes_w_Dloop, 
  nodes_CDS_only_all, 
  by = "clade_tips", 
  suffix =  c("_w_Dloop", "_CDS_only")
) %>%
  mutate(
    # 1. Classify node status
    Node_Status_Raw = case_when(
      is.na(posterior_CDS_only) ~ label_novel,
      posterior_CDS_only < support_cutoff ~ label_drop,
      TRUE ~ label_both
    ),
    Node_Status = factor(Node_Status_Raw, levels = c(label_both, label_drop, label_novel)),
    
    # 2. Manually create size categories
    Support_Class_Raw = case_when(
      posterior_w_Dloop >= 1.0  ~ "1.00",
      posterior_w_Dloop >= 0.95 ~ "0.95 - 0.99",
      TRUE                      ~ "0.85 - 0.95"
    ),
    # Order the factor so the legend displays them from smallest to largest
    Support_Class = factor(Support_Class_Raw, levels = c("0.85 - 0.95", "0.95 - 0.99", "1.00")),
    
    # 3. Handle missing Y coordinates safely
    height_CDS_only_plot = ifelse(is.na(height_CDS_only), 0, height_CDS_only),
    
    # 4. Unpack HPD lists safely
    hpd_min_w_Dloop  = map_dbl(height_0.95_HPD_w_Dloop, 1),
    hpd_max_w_Dloop  = map_dbl(height_0.95_HPD_w_Dloop, 2),
    hpd_min_CDS_only = map_dbl(height_0.95_HPD_CDS_only, ~ ifelse(is.null(.x), NA, .x[1])),
    hpd_max_CDS_only = map_dbl(height_0.95_HPD_CDS_only, ~ ifelse(is.null(.x), NA, .x[2]))
  ) %>%
  # 5.
  # Sorting the dataframe by this factor ensures background dots are evaluated 
  # first by geom_point, leaving rare/novel nodes to be painted last (on top).
  arrange(Node_Status)


# ==============================================================================
# STEP 5: Generate and Display Plot
# ==============================================================================
message("Generating ggplot visualization...")

# Establish color assignments anchored explicitly to factor levels
color_palette <- c("#2b8cbe", "#fe9929", "#e31a1c")
names(color_palette) <- c(label_both, label_drop, label_novel)

# Define an explicit named vector for manual dot sizes
size_palette <- c("0.85 - 0.95" = 1.5, "0.95 - 0.99" = 3.5, "1.00" = 6.0)

divergence_plot <- ggplot(plot_data_clean, aes(x = height_w_Dloop, y = height_CDS_only_plot)) +
  
  # 1:1 line representing identical estimated node ages
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  
  # Error bars (Drawn underneath the actual points for clean aesthetics)
  geom_errorbar(aes(ymin = hpd_min_CDS_only, ymax = hpd_max_CDS_only, color = Node_Status), alpha = 0.3) +
  geom_errorbarh(aes(xmin = hpd_min_w_Dloop, xmax = hpd_max_w_Dloop, color = Node_Status), alpha = 0.3) +
  
  # Size is now mapped to our explicit discrete factor variable
  geom_point(aes(color = Node_Status, size = Support_Class), alpha = 0.8) +
  
  # Highlights the baseline axis where the unmapped novel nodes live
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray70") +
  
  # Formatting
  scale_color_manual(values = color_palette) +
  scale_size_manual(values = size_palette) + 
  
  labs(
    title = "False Killer Whale Mitogenome Divergence Time Comparison",
#    subtitle = paste0("Anchored on nodes with posterior >= ", support_cutoff, " in the w_Dloop tree.\nNodes on the Y=0 line are unique to the D-loop data."),
    x = "Divergence Time (with Dloop)",
    y = "Divergence Time (w/o Dloop)",
    size = "with Dloop Posterior Support",
    color = "Node Status"
  ) +
  theme_minimal() +
  theme(
    # --- LEGEND LAYOUT FIXES ---
    legend.position = "right",        # Moves the legend to the right side
    legend.direction = "vertical",    # Orients each internal legend vertically
    legend.box = "vertical",          # Stacks the Color and Size legends vertically
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.minor = element_blank()
  )

print(divergence_plot)

# Do fancy zooming where I end up with two plots
final_zoom_plot <- divergence_plot +
  facet_zoom(
    xlim = c(0, 0.05), 
    ylim = c(0, 0.05),
    horizontal = FALSE # Stacks the zoom panel vertically underneath the main plot
  )

print(final_zoom_plot)

# Standard zoom where the oldest node is excluded
zoomed_divergence_plot <- divergence_plot +
  coord_cartesian(xlim = c(0, 0.085), ylim = c(0, 0.085))

# Print the zoomed version
pdf(file = 'results-raw/compare_trees_zoomed_plot.pdf', height = 7, width = 9)
print(zoomed_divergence_plot)
dev.off()
