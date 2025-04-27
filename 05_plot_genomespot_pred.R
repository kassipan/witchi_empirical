library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(patchwork)
library(forcats)

taxonomy.dt <- fread("data/taxonomy_map.tsv", sep = "\t", header = FALSE, 
                     col.names = c("genome_name", "ncbi_accesison", "classification"))

ind <- "genomespot_predictions/"
file.v <- list.files(ind, pattern = "*.predictions.tsv")

pred.dt <- data.table()
for (file in file.v){
  genome_name <- gsub(".predictions.tsv", "", file)
  tmp.dt <- fread(sprintf("%s/%s", ind, file), sep = "\t")
  tmp.dt[, genome_name := genome_name]
  pred.dt <- rbind(pred.dt, tmp.dt)
}

#Phyla per "Supergroups" 
eury <- c("p__JACRDV01", "p__Hydrothermarchaeota", "p__Hadarchaeota", "p__Halobacteriota",
         "p__Methanobacteriota_A", "p__Methanobacteriota_B", "p__Methanobacteriota",
         "p__Thermoplasmatota")
tack <- c("p__Thermoproteota")
dpann <- c("p__SpSt-1190", "p__EX4484-52", "p__B1Sed10-29", "p__Undinarchaeota",
           "p__Huberarchaeota", "p__Iainarchaeota", "p__Aenigmatarchaeota",
           "p__Altiarchaeota", "p__Micrarchaeota", "p__Nanoarchaeota", "p__Nanohaloarchaeota")
asgard <- c("p__Asgardarchaeota")

#Modify taxonomy  
pred.dt.ext <- pred.dt %>%
  left_join(taxonomy.dt, by = "genome_name") %>%
  mutate(Phylum = gsub(".*p__(.*?);.*", "p__\\1", classification)) %>%
  mutate(Class = gsub(".*c__(.*?);.*", "c__\\1", classification)) %>%
  mutate(Order = gsub(".*o__(.*?);.*", "o__\\1", classification)) %>%
  mutate(Family = gsub(".*f__(.*?);.*", "f__\\1", classification)) %>%
  filter(!grepl("genome missing features", warning)) %>%
  filter(!(grepl("temperature_min", target) & grepl("min_exceeded", warning))) %>%
  filter(!(grepl("oxygen", target) & error < 0.75)) %>%
  mutate(Group = case_when(
      grepl(paste(eury, collapse = "|"), Phylum) ~ "Eury",
      grepl(paste(tack, collapse = "|"), Phylum) ~ "TACK",
      grepl(paste(dpann, collapse = "|"), Phylum) ~ "DPANN",
      TRUE ~ "Asgard"  
    ))

write.table(pred.dt.ext, file = "data/genomespot_pred_arcGTDBr220.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

#Keep only optimum values
pred.optima <- pred.dt.ext %>%
  filter(!grepl("ph_max|ph_min|salinity_max|salinity_min|temperature_max|temperature_min", target))

split.dt <- split(pred.optima, pred.optima$target)
lapply(names(split.dt), function(x) assign(x, split.dt[[x]], envir = .GlobalEnv))

#Temperature
temperature_optimum$value <- as.numeric(temperature_optimum$value)
breaks_temp <- c(-Inf, 15, 30, 45, 60, 80, Inf)
labels_temp <- c("0", "15", "30", "45", "60", ">80")
temperature_optimum$category <- cut(temperature_optimum$value, breaks = breaks_temp, 
                                    labels = labels_temp, right = FALSE)
colors_temp <- c("0" = "#1A3360", "15" = "#4094C4", "30" = "#D1E7F3", "45" = "#FDDBC7", 
                 "60" = "#D9624E", ">80" = "#660821")
##Class
temperature_optimum_per_class <- temperature_optimum %>%
    mutate(category = factor(category, levels = rev(labels_temp))) %>%
    group_by(Class, category, Group) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(Class) %>%
    mutate(
      total_count = sum(count),
      percentage = (count / sum(count)) * 100,
      Class = paste0(Class, " (n=",total_count,")")) %>%
    ungroup() %>%
    mutate(Class = fct_rev(factor(Class))) %>%
  ggplot(aes(x = percentage, y = Class, fill = category)) +
    geom_col(position = "stack", width = 0.8) +
    scale_fill_manual(values = colors_temp) +
    theme_light() +
    ggtitle("Optimal growth temperature per Archaeal class") +
    theme(
      strip.text.x = element_text(size = 12, face = "bold"),
      panel.grid = element_blank(),
      legend.position = "right",
      plot.title = element_text(size = 14, face = "bold")) +
  labs(x = "no. of genomes (%)", y = NULL, 
       fill = "Temperature") +
  facet_grid(Group~., scales = "free_y", space = "free")

pdf(file = "figures/temperature_optimum_per_class.pdf", width = 8, height = 13)
temperature_optimum_per_class
dev.off() 

png(file = "figures/temperature_optimum_per_class.png", 
    width = 10000, height = 12000, res = 1200)
temperature_optimum_per_class
dev.off()


##Family
temperature_optimum_per_family <- temperature_optimum %>%
  mutate(category = factor(category, levels = rev(labels_temp))) %>%
  group_by(Family, category, Group) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Family) %>%
  mutate(
    total_count = sum(count),
    percentage = (count / sum(count)) * 100,
    Family = paste0(Family, " (n=",total_count,")")) %>%
  ungroup() %>%
  mutate(Family = fct_rev(factor(Family))) %>%
  ggplot(aes(x = percentage, y = Family, fill = category)) +
  geom_col(position = "stack", width = 0.8) +
  scale_fill_manual(values = colors_temp) +
  theme_light() +
  ggtitle("Optimal growth temperature per Archaeal family") +
  theme(
    strip.text.x = element_text(size = 12, face = "bold"),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.title = element_text(size = 14, face = "bold")) +
  labs(x = "no. of genomes (%)", y = NULL, 
       fill = "Temperature") +
  facet_grid(Group~., scales = "free_y", space = "free")

pdf(file = "figures/temperature_optimum_per_family.pdf", width = 15, height = 60)
temperature_optimum_per_family
dev.off()  

#Salinity
salinity_optimum$value <- as.numeric(salinity_optimum$value)
breaks_sal <- c(-Inf, 0.5, 1, 2.5, 5, 10, 15, Inf)
labels_sal <- c("0", "0.5", "1", "2.5", "5", "10", ">15")
salinity_optimum$category <- cut(salinity_optimum$value, breaks = breaks_sal, 
                                    labels = labels_sal, right = FALSE)
colors_sal <- c("0" = "#1A3360", "0.5" = "#4094C4", "1" = "#D1E7F3", "2.5" = "#F2F2F2", 
                "5" = "#FDDBC7", "10" = "#D9624E", ">15" = "#660821")

## Class
salinity_optimum_per_class <- salinity_optimum %>%
  mutate(category = factor(category, levels = rev(labels_sal))) %>%
  group_by(Class, category, Group) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Class) %>%
  mutate(
    total_count = sum(count),
    percentage = (count / sum(count)) * 100,
    Class = paste0(Class, " (n=",total_count,")")) %>%
  ungroup() %>%
  mutate(Class = fct_rev(factor(Class))) %>%
  ggplot(aes(x = percentage, y = Class, fill = category)) +
  geom_col(position = "stack", width = 0.8) +
  scale_fill_manual(values = colors_sal) +
  theme_light() +
  ggtitle("Optimum Salinity per archaeal class") +
  theme(
    strip.text.x = element_text(size = 12, face = "bold"),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.title = element_text(size = 14, face = "bold")) +
  labs(x = "no. of genomes (%)", y = NULL, 
       fill = "Salinity (%NaCl)") +
  facet_grid(Group~., scales = "free_y", space = "free")

#pdf(file = "figures/salinity_optimum_per_class.pdf", width = 15, height = 20)
#salinity_optimum_per_class
#dev.off()  

##Family
salinity_optimum_per_family <- salinity_optimum %>%
  mutate(category = factor(category, levels = rev(labels_sal))) %>%
  group_by(Family, category, Group) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Family) %>%
  mutate(
    total_count = sum(count),
    percentage = (count / sum(count)) * 100,
    Family = paste0(Family, " (n=",total_count,")")) %>%
  ungroup() %>%
  mutate(Family = fct_rev(factor(Family))) %>%
  ggplot(aes(x = percentage, y = Family, fill = category)) +
  geom_col(position = "stack", width = 0.8) +
  scale_fill_manual(values = colors_sal) +
  theme_light() +
  ggtitle("Optimum Salinity per archaeal family") +
  theme(
    strip.text.x = element_text(size = 12, face = "bold"),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.title = element_text(size = 14, face = "bold")) +
  labs(x = "no. of genomes (%)", y = NULL, 
       fill = "Salinity (%NaCl)") +
  facet_grid(Group~., scales = "free_y", space = "free")

#pdf(file = "figures/salinity_optimum_per_family.pdf", width = 15, height = 60)
#salinity_optimum_per_family
#dev.off() 

#pH
ph_optimum$value <- as.numeric(ph_optimum$value)
breaks_ph <- c(-Inf, 4, 5, 6, 7, 8, 9, 10, 11, Inf)
labels_ph <- c("<4", "4", "5", "6", "7", "8", "9", "10", ">11")
ph_optimum$category <- cut(ph_optimum$value, breaks = breaks_ph, labels = labels_ph,
                           right = FALSE)
colors_ph <- c("<4" = "#1A3360", "4" = "#2173B4", "5" = "#6BADD3", "6" = "#C3DFEF", 
               "7" = "#F9F9FA", "8" = "#FCCDB5", "9" = "#E68166", "10" = "#BB2933", 
               ">11" = "#660821")

##Class
ph_optimum_per_class <- 
  ph_optimum %>%
  mutate(category = factor(category, levels = rev(labels_ph))) %>%
  group_by(Class, category, Group) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Class) %>%
  mutate(
    total_count = sum(count),
    percentage = (count / sum(count)) * 100,
    Class = paste0(Class, " (n=",total_count,")")) %>%
  ungroup() %>%
  mutate(Class = fct_rev(factor(Class))) %>%
  ggplot(aes(x = percentage, y = Class, fill = category)) +
  geom_col(position = "stack", width = 0.8) +
  scale_fill_manual(values = colors_ph) +
  theme_light() +
  ggtitle("Optimum pH per archaeal class") +
  theme(
    strip.text.x = element_text(size = 12, face = "bold"),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.title = element_text(size = 14, face = "bold")) +
  labs(x = "no. of genomes (%)", y = NULL, 
       fill = "pH") +
  facet_grid(Group~., scales = "free_y", space = "free")

#pdf(file = "figures/ph_optimum_per_class.pdf", width = 15, height = 20)
#ph_optimum_per_class
#dev.off()  

##Family
ph_optimum_per_family <- 
  ph_optimum %>%
  mutate(category = factor(category, levels = rev(labels_ph))) %>%
  group_by(Family, category, Group) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Family) %>%
  mutate(
    total_count = sum(count),
    percentage = (count / sum(count)) * 100,
    Family = paste0(Family, " (n=",total_count,")")) %>%
  ungroup() %>%
  mutate(Family = fct_rev(factor(Family))) %>%
  ggplot(aes(x = percentage, y = Family, fill = category)) +
  geom_col(position = "stack", width = 0.8) +
  scale_fill_manual(values = colors_ph) +
  theme_light() +
  ggtitle("Optimum pH per archaeal family") +
  theme(
    strip.text.x = element_text(size = 12, face = "bold"),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.title = element_text(size = 14, face = "bold")) +
  labs(x = "no. of genomes (%)", y = NULL, 
       fill = "pH") +
  facet_grid(Group~., scales = "free_y", space = "free")

#pdf(file = "figures/ph_optimum_per_family.pdf", width = 15, height = 60)
#ph_optimum_per_family
#dev.off()  
