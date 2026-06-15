library(ape)
library(tidyverse)
library(viridis)
library(scales)
library(xml2)
library(rsvg)
load('../Pcra.database.data/data/ASPhyloCR.rda')
load('../Pcra.database.data/data/Pcra.strata.rda')
load('../Pcra.database.data/data/CR.haps.rda')
load('data/broad_stratum_colors.rda')

project <- 'Pcra_CR_Broad'

#######################################################################
# Generate PopArt input file
PopArt_data <- left_join(ASPhyloCR, CR.haps) |> 
  left_join(Pcra.strata) |> 
  mutate(Broad = ifelse(Broad %in% c('Atl-weird', 'Atl-normal', 'Atl-unknown'), 'N_Atlantic', Broad),
         Broad = ifelse(Broad == 'Atl-South', 'S_Atlantic', Broad),
         Broad = ifelse(Broad == 'Spac', 'S_Pacific', Broad),
         Broad = ifelse(Broad == 'Wpac', 'NPac_Western', Broad),
         Broad = ifelse(Broad == 'Epac', 'NPac_Eastern', Broad),
         Broad = ifelse(Broad == 'CNP', 'NPac_Central', Broad),
         Broad = ifelse(Broad == 'Indian Ocean', 'Indian_Ocean', Broad))

# add CR.haps for the 7 Aus residents that were sequenced in Australia
PopArt_data <- bind_rows(
  PopArt_data, 
  tibble(
    Animal.ID = paste0('AusRes', 1:7), 
    CR.hap = 45,
    Location = 'Australia',
    Fine = 'Australia',
    Broad = 'Aus_resident',
    Ocean_Basin = 'Indian_Ocean',
    mito.Broad = 'Aus_resident'
  )
)


# PopArt can't deal with ambiguous sights, so I have to collapse haplotypes with
# SNPs into their more common parent haplotype
PopArt_data$CR.hap[which(PopArt_data$CR.hap == 37)] <- 9 # hap 37 becomes 9
PopArt_data$CR.hap[which(PopArt_data$CR.hap %in% c(39, 40))] <- 38 # haps 39 and 40 become 38
PopArt_data$CR.hap[which(PopArt_data$CR.hap == 41)] <- 17 # hap 41 becomes 17

# Count haplotype frequencies by Broad stratification
haplotype_freq <- PopArt_data %>%
  group_by(CR.hap, Broad) %>%
  summarise(n = n(), .groups = 'drop') %>%
  pivot_wider(names_from = Broad, values_from = n, values_fill = 0) %>%
  arrange(CR.hap)

# Re-arrange columns so that strata are in the order I want them in the legend
haplotype_freq <- haplotype_freq[,c(1,11,8,4,9,6,5,7,2,3,12,10)]

# Get the stratum names (trait labels) in order
trait_labels <- setdiff(names(haplotype_freq), "CR.hap")
n_traits <- length(trait_labels)

# Create the output file
output_file <- paste0("PopArt/", project, "_popart_traits.txt")

# Write the trait block
cat("begin traits;\n", sep = "", file = output_file)
cat("Dimensions NTRAITS=", n_traits, ";\n", sep = "", file = output_file, append = TRUE)
cat("Format labels=yes missing=? separator=Comma;\n", file = output_file, append = TRUE)
cat("TraitLabels ", paste(trait_labels, collapse = " "), ";\n", sep = "", file = output_file, append = TRUE)
cat("Matrix\n", file = output_file, append = TRUE)

# Write the matrix data
haplotype_freq %>%
  rowwise() %>%
  mutate(row_text = paste(CR.hap, paste(c_across(all_of(trait_labels)), collapse = ","), sep = " ")) %>%
  pull(row_text) %>%
  walk(~cat(.x, "\n", file = output_file, append = TRUE))

cat(";\n", file = output_file, append = TRUE)

#####################################################
# Use PopArt to generate the MJN, export as a .svg
#
#  REMEMBER TO:
#  1) Change the color of the unknown stratum so that
#     it's different from N. Atlantic
#  2) Export as an .svg file
#
#####################################################

network_name <- paste0(project, '_MJN')

# Read the SVG
svg <- read_xml(paste0('PopArt/', network_name, '.svg'))

# First, let's examine what colors PopArt used
# Get all fill colors in the SVG
fills <- xml_find_all(svg, "//*[@fill]")
fill_colors <- unique(xml_attr(fills, "fill"))
strokes <- xml_find_all(svg, "//*[@stroke]")
stroke_colors <- unique(xml_attr(strokes, "stroke"))
names(fill_colors) <- fill_colors

# Plot them to see what they look like
show_col(fill_colors[!is.na(fill_colors) & fill_colors != "none"])

color_map <- fill_colors
color_map[3:7] <- viridis_pal()(5) #N_Alantic thru NPac_Western
color_map[3] <- "#481568"
color_map[8] <- "#caa611" # NPac_Central
color_map[9] <- "#E5C61BFF" # NPac_Eastern
color_map[10] <- "#BEBEBE" #MHI, will add cross-hatching
color_map[11] <- "#BEBEBE" #NWHI, will add cross-hatching
color_map[12] <- "#666666" #Aus resident, will add cross-hatching
color_map[13] <- "white" #Unknown
show_col(color_map[color_map!= "none"])

# Function to replace colors
replace_colors <- function(svg, color_map) {
  for (old_color in names(color_map)) {
    new_color <- color_map[old_color]
    
    # Replace in fill attributes
    fills <- xml_find_all(svg, sprintf("//*[@fill='%s']", old_color))
    xml_set_attr(fills, "fill", new_color)
    
    # Replace in stroke attributes
    strokes <- xml_find_all(svg, sprintf("//*[@stroke='%s']", old_color))
    xml_set_attr(strokes, "stroke", new_color)
    
    # Replace in style attributes (sometimes colors are in CSS style)
    styled <- xml_find_all(svg, "//*[@style]")
    for (node in styled) {
      style <- xml_attr(node, "style")
      if (!is.na(style) && grepl(old_color, style, fixed = TRUE)) {
        new_style <- gsub(old_color, new_color, style, fixed = TRUE)
        xml_set_attr(node, "style", new_style)
      }
    }
  }
  return(svg)
}

# Apply color replacements
svg_modified <- replace_colors(svg, color_map)

# Save the modified SVG
write_xml(svg_modified, paste0('PopArt/', network_name, '_recolored.svg'))
rsvg_pdf(paste0('PopArt/', network_name, '_recolored.svg'), 
         file = paste0('PopArt/', network_name, '_final.pdf'))

#######################################################
# Add cross-hatching to the resident strata
# Define a pattern (this needs to go in a <defs> section)
svg <- read_xml(paste0('PopArt/', network_name, '_recolored.svg'))

target_colors <- c("#C9A511", "#22948FFF")

# Remove any existing patterns first
defs <- xml_find_first(svg, "//defs")
if (length(defs) > 0) {
  for (i in 1:2) {
    existing_pattern <- xml_find_first(defs, sprintf("//pattern[@id='diagonalHatch%d']", i))
    if (length(existing_pattern) > 0) {
      xml_remove(existing_pattern)
    }
  }
}

# If no defs, create one
if (length(defs) == 0) {
  svg_root <- xml_root(svg)
  defs <- xml_add_child(svg_root, "defs", .where = 0)
}

# Create patterns for each color
for (i in seq_along(target_colors)) {
  pattern_str <- sprintf('
<pattern id="diagonalHatch%d" patternUnits="userSpaceOnUse" width="6" height="6">
  <rect width="6" height="6" style="fill:%s;stroke:none"/>
  <line x1="0" y1="6" x2="6" y2="0" style="stroke:black;stroke-width:1"/>
  <line x1="-1" y1="1" x2="1" y2="-1" style="stroke:black;stroke-width:1"/>
  <line x1="5" y1="7" x2="7" y2="5" style="stroke:black;stroke-width:1"/>
</pattern>', i, target_colors[i])
  
  pattern_node <- read_xml(pattern_str)
  xml_add_child(defs, pattern_node)
}

# Find and apply patterns for each color
all_elements <- xml_find_all(svg, "//*")

for (i in seq_along(target_colors)) {
  target_color <- target_colors[i]
  elements_to_pattern <- list()
  
  for (j in seq_along(all_elements)) {
    elem <- all_elements[[j]]
    attrs <- xml_attrs(elem)
    
    if (any(grepl(target_color, attrs, fixed = TRUE))) {
      elements_to_pattern[[length(elements_to_pattern) + 1]] <- elem
    }
  }
  
  print(paste("Found", length(elements_to_pattern), "elements for color", target_color))
  
  for (elem in elements_to_pattern) {
    xml_set_attr(elem, "fill", sprintf("url(#diagonalHatch%d)", i))
  }
}

write_xml(svg, paste0('PopArt/', network_name, '_with_patterns.svg'))
cat("Cross-hatching applied to both colors!\n")
