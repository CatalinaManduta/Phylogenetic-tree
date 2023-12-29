# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Project: Contraception
# Version: 01
# Author: Johanna Staerk
# Description: Phylogenetic tree
# Adaptation: Catalina Manduta 
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# ==== INSTALL PACKAGES: ====
list_of_packages <- c("tidyverse", "ape", "ggtree", "ggnewscale", "ggtreeExtra", "patchwork", "plotly", "ggimage")
install.packages(list_of_packages, dependencies = TRUE)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtreeExtra")

install.packages("ggstance")

library(tidyverse)
library(ape)
library(ggtree)
library(ggnewscale)
library(ggtreeExtra)
library(patchwork)
library(plotly)
library(ggimage)
library(dplyr)
library(ggstance)
library(readr)


# ==== DATA LOADING====
# life expectancies
dif <- read_csv("lifeExpBaSTA_new.csv")

# read consensus tree mammals
load("/work/CatalinaThesis/analysis/data/rdata/trees/maxCredTreeMammalia.RData")
# read consensus tree mammals
tree <- maxCred

# # order names from vertlife
tax <- MDD_v1_4_6533species

# ==== DATA PREPARATION: ====
# Subset columns
dif <- dif %>% select(species,Male_None_Mean,Male_Surgical_Mean)

# calculated the relative differences:
dif$dif_surg_male <-(dif$Male_None_Mean - dif$Male_Surgical_Mean) / dif$Male_None_Mean
dif$dif_surg_male_2 <- abs(dif$dif_surg_male)
dif$surgical_male <-
  ifelse(dif$dif_surg_male < 0, "surgical", "no contraception")

# convert NA's to no data
table(dif$surgical_male)
dif$surgical_male[which(is.na(dif$surgical_male))] <- "no data"
dif$dif_surg_male_2[which(is.na(dif$dif_surg_male_2))] <- 0

# subset df for which we have data
dif <- dif %>% filter(!is.na(dif_surg_male))

tree.label <- str_replace(tree$tip.label, "_", " ")
print(tree.label)
dif <- dif[which(dif$species %in% tree.label), ]

# Replace names on tree to match ZIMS names
dif[which(!(dif$species %in% tree.label)), ]

tree.label[which(!(tree.label %in% dif$species))]
tree$tip.label[which(tree$tip.label == "Bubalus_arnee")] <- "Bubalus_bubalis"
tree$tip.label[which(tree$tip.label == "Equus_africanus")] <- "Equus_asinus"
tree$tip.label[which(tree$tip.label == "Aonyx_cinerea")] <- "Aonyx_cinereus"

# NOTE Cervus canandensis is considered the same species as Cervus elaphus on the tree and cannot be matched to tree!
dif[which(!(dif$species %in% tree.label)), ]
dif <- dif %>% filter(!species == "Pseudocheirus peregrinus") # Outlier

# Convert the tibble to a data frame
dif <- as.data.frame(dif)

# Replace spaces with underscores in the species column
dif$species <- gsub(" ", "_", dif$species)

# taxonomy
tax <-
  tax %>% select(sciName,
                 majorType:majorSubtype,
                 vert.order = order,
                 mainCommonName)

# adjust taxonomy in the vertlife file so that it fits with our data
tax[which(tax$sciName == "Notamacropus_eugenii"), "sciName"] <-
  "Macropus eugenii"
tax[which(tax$sciName == "Notamacropus_parma"), "sciName"] <-
  "Macropus parma"
tax[which(tax$sciName == "Osphranter_robustus"), "sciName"] <-
  "Macropus robustus"
tax[which(tax$sciName == "Notamacropus_rufogriseus"), "sciName"] <-
  "Macropus rufogriseus"
tax[which(tax$sciName == "Osphranter_rufus"), "sciName"] <-
  "Macropus rufus"
tax[which(tax$sciName == "Dicotyles_tajacu"), "sciName"] <-
  "Pecari tajacu"
tax[which(tax$sciName == "Lama_pacos"), "sciName"] <-
  "Vicugna pacos"
tax[which(tax$sciName == "Lama_vicugna"), "sciName"] <-
  "Vicugna vicugna"

dif$id <- rownames(dif)
taxm <-
  dif %>% left_join(tax, by = c("species" = "sciName")) %>% select(tip = "species", order = vert.order)
taxm[which(taxm$tip == "Bison_bison"), "order"] <- "ARTIODACTYLA"
taxm[which(taxm$tip == "Vicugna_pacos"), "order"] <- "ARTIODACTYLA"
taxm[which(taxm$tip == "Vicugna_vicugna"), "order"] <- "ARTIODACTYLA"


# Subset the tree
keep <- dif$species

to_drop <-
  tree$tip.label[which(!(tree$tip.label %in% keep))]
tree_subset <- drop.tip(tree, to_drop)


# add order names to the correct positions on tree
# convert phylo object to tibble
tree.tibble <- as_tibble(tree_subset)

ord <- unique(taxm$order)

# find nodes of each order
dford <- data.frame(order = ord, node = NA)
dford <- dford[!is.na(dford$order), ]

for (i in ord) {
  tip = taxm[which(taxm$order == i), "tip"]
  if (length(tip) > 1) {
    dford[which(dford$order == i), "node"] <-
      getMRCA(tree_subset, tip = tip) # function to determine common ancestor node
  } else {
    # if order contains only one species we use the node of that species
    dford[which(dford$order == i), "node"] <-
      tree.tibble[which(tree.tibble$label == taxm[which(taxm$order == i), "tip"]), "node"]
  }
}
# convert to phylo obj
tree_subset <- as.phylo(tree_subset)


# ==== PLOT TREE: ====
d1 <- dif %>% select(species, dif_surg_male_2, dif_surg_male, surgical_male) %>% as.data.frame()
d1 <- d1[which(!(is.na(d1$dif_surg_male_2))), ]
d1 <- d1[which(!(is.na(d1$dif_surg_male))), ]
# This line removes rows from d1 where the dif_surg_male_2 column is NA

all_species <- unique(tree_subset$tip.label)

# Create a data frame of all species
all_species_df <- data.frame(species = all_species)

# Merge with d1, filling in NA for missing values
d1 <- merge(all_species_df, d1, by = "species", all = TRUE)

p <- ggtree(
  tree_subset,
  layout = "rectangular",
  ladderize = TRUE,
  size = 0.5,
  color = "#454747"
) + geom_tiplab(offset = 1, size = 2)
p
pal1 <- c("red", "#0099f9") 

# Plot without the legend for more space for the tree
p2 <- facet_plot(
  p,
  'Non-castrated vs castrated males',
  data = d1,
  geom = geom_segment,
  mapping = aes(
    x = 0,
    xend = dif_surg_male_2,
    y = y,
    yend = y,
    color = surgical_male),size = 2) +
  scale_color_manual(values = pal1) +
  xlim_expand(c(0, 200), panel = "Tree") +
  theme(legend.position = "none")


# Convert the tip.label to a data.frame
tree_subset_df <- data.frame(id = tree_subset$tip.label)

# Merge the dataframes for the tips of the tree
tree_species <- data.frame(species = tree_subset$tip.label, 
                           order = 1:length(tree_subset$tip.label))
merged_data <- merge(tree_species, d1, by = "species", all.x = TRUE)
merged_data <- merged_data[order(merged_data$order), ]
d2 <- data.frame(id = merged_data$species, location = merged_data$surgical_male)


d2 <- data.frame(id=tree_subset$tip.label, location = d1$surgical_male)
p2 <- p2 %<+% d2 + geom_tippoint(aes(color=location))


# Add the images from he pylopic website
phylopic_info <- data.frame(node = c(229,180,150,257,142,125,226,262,10,124),
                            phylopic = c("a1116e25-7b50-4666-bef5-de18b6e2778c",
                                         "a0afbb0a-2bb8-4d69-999d-3ed3a09e9966",
                                         "0174801d-15a6-4668-bfe0-4c421fbe51e8",
                                         "fb007db0-cc86-4215-afdd-18f954244e2f",
                                         "c4572239-3b7c-4259-8049-3f44e6d52e6f",
                                         "189b4f43-4caf-4c89-ac9c-bd4695a668de",
                                         "623da06a-1857-43a2-9d1b-8a29f5f378cf",
                                         "8a06d489-024f-4505-8ccb-f86e84e00e75",
                                         "4f228457-6c87-46d8-8748-5a33c98e0f16",
                                         "bae84fd8-937a-4390-b123-0a22bfcd7df7"))

p2 <- p2 %<+% phylopic_info + geom_nodelab(aes(image=phylopic), geom="phylopic", alpha=.9, color='black')

p2 <- p2 + theme(plot.background = element_rect(fill = "grey"),
                 panel.background = element_rect(fill = "grey"))

#Save as a png file
ggsave("plot5.png", p2, width = 10, height = 10, dpi = 300)
