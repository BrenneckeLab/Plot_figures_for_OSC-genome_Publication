library(tidyverse)
library(cowplot)
library(foreach)
library(patchwork)
library(khroma)
library(paletteer)

theme_set(theme_bw())

setwd("~OSC-paper//Fig4/")

###############################################################################################
###############################################################################################
###############################################################################################
# ###########################   plot  IDR tething analysis - local examples   ##############################
###############################################################################################
###############################################################################################
###############################################################################################

# Read and pivot
RAW <- read_tsv('flam_tiles.sRNA.counts.txt', col_names = TRUE)%>%
  filter(uniqPERC>50, str_detect(LIBRARY,"days"), !LIBRARY=='1xIDR-3days',!str_detect(LIBRARY,"1xIDR"))%>%
  separate(ID, into = c("NAME", "POS"), sep = "_", remove = FALSE, convert = TRUE) %>% 
  mutate(POS=POS*100)%>%
  filter(POS<1000000)%>%
  mutate(
    CLASS = case_when(
      str_detect(LIBRARY, "1x") ~ "1x",
      str_detect(LIBRARY, "2x") ~ "2x",
      TRUE ~ "other"
    )
  )

TABLE = RAW %>% 
  select(LIBRARY,NAME, POS, sRNA.counts.perPOS) %>%
  pivot_wider(
    names_from = LIBRARY,
    values_from = sRNA.counts.perPOS,
    values_fill = 0
  )%>%
  filter(if_all(contains("days" ), ~ . > 1)) %>% # Filter out rows where all 'days' columns are zero
  mutate(
    #for all columns containing 1xIDR divide by 1xGAL4
    across(starts_with("1xIDR"), ~ . / `1xGAL4-3days`),
    #for all columns containing 2xIDR divide by 2xGAL4
    across(starts_with("2xIDR"), ~ . / `2xGal4-3days`)
  )%>%
  select(-c(`1xGAL4-3days`, `2xGal4-3days`))%>%
  pivot_longer(
    cols = -c(NAME, POS),
    names_to = "LIBRARY",
    values_to = "sRNA.counts.perPOS"
  ) %>%
  mutate(
    CLASS = case_when(
      str_detect(LIBRARY, "1x") ~ "1x",
      str_detect(LIBRARY, "2x") ~ "2x",
      TRUE ~ "other"
    )
  )

###############################################################################################
# flamenco
plotTABLE=TABLE %>%
  filter(NAME=="flam") 

DIST=100000
all_breaks <- seq(
  from = floor(min(plotTABLE$POS, na.rm = TRUE) / DIST) * DIST,
  to   = ceiling(max(plotTABLE$POS, na.rm = TRUE) / DIST) * DIST,
  by   = DIST
)

# Label only every 2nd break (odd indices: 1st, 3rd, 5th, etc.)
labeled_breaks <- all_breaks[seq(1, length(all_breaks), by = 2)]

p=plotTABLE %>%
  mutate(
    GROUP = case_when(
      str_detect(LIBRARY, "1x") ~ "1x",
      str_detect(LIBRARY, "2x") ~ "2x",
      TRUE ~ "other"
    )
  )%>%
  group_by(NAME, POS, GROUP) %>%
  summarise(
    sRNA.counts.perPOS = mean(sRNA.counts.perPOS, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  ggplot(aes(x = POS, y = sRNA.counts.perPOS)) +
  geom_point(aes(text=paste(NAME,POS,sep="_"))) +
  scale_x_continuous(
    breaks = all_breaks,           # All ticks every 100,000
    labels = function(x) {         # Custom labeling function
      ifelse(x %in% labeled_breaks, 
             label_number()(x),    # Label these breaks
             "")                   # Empty label for others
    }
  ) +
  scale_y_continuous(limits=c(0.1,1.2)) +
  theme_cowplot(13)+
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +
  labs(
    x = "Position in the flamenco locus",
    y = "fold change upon IDR tethering",
  ) +
  theme(
    legend.position = "none"
  )
p

ggsave("flamenco_silencing.flam.pdf", p, width = 4, height = 3)

# ggplotly(p, tooltip = "text") %>%
#   layout(hovermode = "closest") %>%
#   config(displayModeBar = FALSE)

###############################################################################################
# 20A
plotTABLE=TABLE %>%
  filter(NAME=="20A") 

DIST=5000
all_breaks <- seq(
  from = floor(min(plotTABLE$POS, na.rm = TRUE) / DIST) * DIST,
  to   = ceiling(max(plotTABLE$POS, na.rm = TRUE) / DIST) * DIST,
  by   = DIST
)

# Label only every 2nd break (odd indices: 1st, 3rd, 5th, etc.)
labeled_breaks <- all_breaks[seq(1, length(all_breaks), by = 2)]

p=plotTABLE %>%
  mutate(
    GROUP = case_when(
      str_detect(LIBRARY, "1x") ~ "1x",
      str_detect(LIBRARY, "2x") ~ "2x",
      TRUE ~ "other"
    )
  )%>%
  group_by(NAME, POS, GROUP) %>%
  summarise(
    sRNA.counts.perPOS = mean(sRNA.counts.perPOS, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  ggplot(aes(x = POS, y = sRNA.counts.perPOS)) +
  geom_point(aes(text=paste(NAME,POS,sep="_"))) +
  scale_x_continuous(
    breaks = all_breaks,           # All ticks every 100,000
    labels = function(x) {         # Custom labeling function
      ifelse(x %in% labeled_breaks, 
             label_number()(x),    # Label these breaks
             "")                   # Empty label for others
    }
  ) +
  scale_y_continuous(limits=c(0.1,1.2)) +
  theme_cowplot(13)+
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +
  labs(
    x = "Position in the 20A locus",
    y = "fold change upon IDR tethering",
  ) +
  theme(
    legend.position = "none"
  )
p

ggsave("flamenco_silencing.20A.pdf", p, width = 4, height = 3)


###############################################################################################
# 77B
plotTABLE=TABLE %>%
  filter(NAME=="77B") 

DIST=5000
all_breaks <- seq(
  from = floor(min(plotTABLE$POS, na.rm = TRUE) / DIST) * DIST,
  to   = ceiling(max(plotTABLE$POS, na.rm = TRUE) / DIST) * DIST,
  by   = DIST
)

# Label only every 2nd break (odd indices: 1st, 3rd, 5th, etc.)
labeled_breaks <- all_breaks[seq(1, length(all_breaks), by = 2)]

p=plotTABLE %>%
  mutate(
    GROUP = case_when(
      str_detect(LIBRARY, "1x") ~ "1x",
      str_detect(LIBRARY, "2x") ~ "2x",
      TRUE ~ "other"
    )
  )%>%
  group_by(NAME, POS, GROUP) %>%
  summarise(
    sRNA.counts.perPOS = mean(sRNA.counts.perPOS, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  ggplot(aes(x = POS, y = sRNA.counts.perPOS)) +
  geom_point(aes(text=paste(NAME,POS,sep="_"))) +
  scale_x_continuous(
    limits=c(0,25000),
    breaks = all_breaks,           # All ticks every 100,000
    labels = function(x) {         # Custom labeling function
      ifelse(x %in% labeled_breaks, 
             label_number()(x),    # Label these breaks
             "")                   # Empty label for others
    }
  ) +
  scale_y_continuous(limits=c(0.1,1.2)) +
  theme_cowplot(13)+
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +
  labs(
    x = "Position in the 77B locus",
    y = "fold change upon IDR tethering",
  ) +
  theme(
    legend.position = "none"
  )
p

ggsave("flamenco_silencing.77B.pdf", p, width = 4, height = 3)



###############################################################################################
# tj

plotTABLE=TABLE %>%
  filter(NAME=="tj") 

DIST=250
all_breaks <- seq(
  from = floor(min(plotTABLE$POS, na.rm = TRUE) / DIST) * DIST,
  to   = ceiling(max(plotTABLE$POS, na.rm = TRUE) / DIST) * DIST,
  by   = DIST
)

# Label only every 2nd break (odd indices: 1st, 3rd, 5th, etc.)
labeled_breaks <- all_breaks[seq(1, length(all_breaks), by = 2)]

p=plotTABLE %>%
  mutate(
    GROUP = case_when(
      str_detect(LIBRARY, "1x") ~ "1x",
      str_detect(LIBRARY, "2x") ~ "2x",
      TRUE ~ "other"
    )
  )%>%
  group_by(NAME, POS, GROUP) %>%
  summarise(
    sRNA.counts.perPOS = mean(sRNA.counts.perPOS, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  ggplot(aes(x = POS, y = sRNA.counts.perPOS)) +
  geom_point(aes(text=paste(NAME,POS,sep="_"))) +
  scale_x_continuous(
    limits=c(0,1500),
    breaks = all_breaks,           # All ticks every 100,000
    labels = function(x) {         # Custom labeling function
      ifelse(x %in% labeled_breaks, 
             label_number()(x),    # Label these breaks
             "")                   # Empty label for others
    }
  ) +
  scale_y_continuous(limits=c(0.1,1.2)) +
  theme_cowplot(13)+
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +
  labs(
    x = "Position in the tj 3' UTR",
    y = "fold change upon IDR tethering",
  ) +
  theme(
    legend.position = "none"
  )
p

ggsave("flamenco_silencing.tj.pdf", p, width = 4, height = 3)


###############################################################################################
# myc

plotTABLE=TABLE %>%
  filter(NAME=="myc") 

DIST=500
all_breaks <- seq(
  from = floor(min(plotTABLE$POS, na.rm = TRUE) / DIST) * DIST,
  to   = ceiling(max(plotTABLE$POS, na.rm = TRUE) / DIST) * DIST,
  by   = DIST
)

# Label only every 2nd break (odd indices: 1st, 3rd, 5th, etc.)
labeled_breaks <- all_breaks[seq(1, length(all_breaks), by = 2)]

p=plotTABLE %>%
  mutate(
    GROUP = case_when(
      str_detect(LIBRARY, "1x") ~ "1x",
      str_detect(LIBRARY, "2x") ~ "2x",
      TRUE ~ "other"
    )
  )%>%
  group_by(NAME, POS, GROUP) %>%
  summarise(
    sRNA.counts.perPOS = mean(sRNA.counts.perPOS, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  ggplot(aes(x = POS, y = sRNA.counts.perPOS)) +
  geom_point(aes(text=paste(NAME,POS,sep="_"))) +
  scale_x_continuous(
    breaks = all_breaks,           # All ticks every 100,000
    labels = function(x) {         # Custom labeling function
      ifelse(x %in% labeled_breaks, 
             label_number()(x),    # Label these breaks
             "")                   # Empty label for others
    }
  ) +
  scale_y_continuous(limits=c(0.1,1.2)) +
  theme_cowplot(13)+
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +
  labs(
    x = "Position in the myc 3' UTR",
    y = "fold change upon IDR tethering",
  ) +
  theme(
    legend.position = "none"
  )
p

ggsave("flamenco_silencing.myc.pdf", p, width = 4, height = 3)

###############################################################################################
###############################################################################################
###############################################################################################
# ###########################   plot  IDR tething analysis - full X chromosome   ##############################
###############################################################################################
###############################################################################################
###############################################################################################


# Read and pivot
RAW <- read_tsv('chrX.sRNA.counts.txt', col_names = TRUE)%>%
  filter(uniqPERC>50, str_detect(LIBRARY,"days"), !LIBRARY=='1xIDR-3days',!str_detect(LIBRARY,"1xIDR"))%>%
  separate(ID, into = c("NAME", "POS","STRAND"), sep = "_", remove = FALSE, convert = TRUE) %>% 
  # mutate(
  #   STRAND=case_when(
  #     str_detect(NAME,"_+") ~ "sense",
  #     str_detect(NAME,"_-") ~ "antisense",
  #     TRUE ~ "unknown"
  #   )
  # )%>%
  # mutate(POS=POS*100)%>%
  # filter(POS<1000000)%>%
  mutate(
    CLASS = case_when(
      str_detect(LIBRARY, "1x") ~ "1x",
      str_detect(LIBRARY, "2x") ~ "2x",
      TRUE ~ "other"
    )
  )

RAW %>% 
  select(LIBRARY,NAME, POS, STRAND,sRNA.counts.perPOS) %>%
  pivot_wider(
    id_cols = c(NAME, POS, STRAND),
    names_from = LIBRARY,
    values_from = sRNA.counts.perPOS,
    values_fill = 0
  )%>%
  filter(POS==9602750)

TABLE = RAW %>% 
  select(LIBRARY,NAME, POS, STRAND,sRNA.counts.perPOS) %>%
  pivot_wider(
    id_cols = c(NAME, POS, STRAND),
    names_from = LIBRARY,
    values_from = sRNA.counts.perPOS,
    values_fill = 0
  )%>%
  filter(if_all(contains("GAL4" ), ~ . >75)) %>% # Filter out rows where all 'days' columns are zero
  mutate(
    #for all columns containing 1xIDR divide by 1xGAL4
    across(starts_with("1xIDR"), ~ . / `1xGAL4-3days`),
    #for all columns containing 2xIDR divide by 2xGAL4
    across(starts_with("2xIDR"), ~ . / `2xGal4-3days`)
  )%>%
  select(-c(`1xGAL4-3days`, `2xGal4-3days`))%>%
  pivot_longer(
    cols = -c(NAME, POS, STRAND),
    names_to = "LIBRARY",
    values_to = "sRNA.counts.perPOS"
  ) %>%
  mutate(
    CLASS = case_when(
      str_detect(LIBRARY, "1x") ~ "1x",
      str_detect(LIBRARY, "2x") ~ "2x",
      TRUE ~ "other"
    )
  )

DIST=2000000
all_breaks <- seq(
  from = floor(min(TABLE$POS, na.rm = TRUE) / DIST) * DIST,
  to   = ceiling(max(TABLE$POS, na.rm = TRUE) / DIST) * DIST,
  by   = DIST
)

# Label only every 2nd break (odd indices: 1st, 3rd, 5th, etc.)
labeled_breaks <- all_breaks[seq(1, length(all_breaks), by = 2)]


plotTABLE = TABLE %>%
  mutate(
    GROUP = case_when(
      str_detect(LIBRARY, "1x") ~ "1x",
      str_detect(LIBRARY, "2x") ~ "2x",
      TRUE ~ "other"
    )
  )%>%
  group_by(NAME, POS, GROUP,STRAND) %>%
  summarise(
    sRNA.counts.perPOS = mean(sRNA.counts.perPOS, na.rm = TRUE),
    .groups = 'drop'
  )

p= plotTABLE %>%
  ggplot(aes(x = POS, y = sRNA.counts.perPOS)) +
  geom_rect(
    aes(xmin = 22673638, xmax = 24100000, ymin = 0.1, ymax = 3),
    fill = "orange",    # Color of the shaded region
    alpha = 0.005       # Transparency of the shaded region
  ) +
  geom_point(aes(text=paste(NAME,POS,sep="_")), size=0.3) +
  scale_x_continuous(
    breaks = all_breaks,           # All ticks every 100,000
    labels = function(x) {         # Custom labeling function
      ifelse(x %in% labeled_breaks, 
             label_number()(x),    # Label these breaks
             "")                   # Empty label for others
    }
  ) +
  # scale_y_continuous(limits=c(0.1,1.2)) +
  # facet_col(~LIBRARY)+
  scale_y_log10()+
  theme_cowplot(14)+
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +
  labs(
    x = "Position in the flamenco locus",
    y = "fold change upon IDR tethering",
  ) +
  theme(
    legend.position = "none"
  )
p
ggsave("chrX_sRNAcounts.pdf",p, width = 8, height = 4)


REGstart=22503638
REGend=24200000
DIST=100000
all_breaks <- seq(
  from = floor(REGstart / DIST) * DIST,
  to   = ceiling(REGend / DIST) * DIST,
  by   = DIST
)

# Label only every 2nd break (odd indices: 1st, 3rd, 5th, etc.)
labeled_breaks <- all_breaks[seq(1, length(all_breaks), by = 2)]

p = plotTABLE %>% 
  filter(POS> REGstart   & POS < REGend)%>%
  ggplot(aes(x = POS, y = sRNA.counts.perPOS)) +
  geom_rect(
    aes(xmin = 22673638, xmax = 24100000, ymin = 0.1, ymax = 3),
    fill = "orange",    # Color of the shaded region
    alpha = 0.005       # Transparency of the shaded region
  ) +
  geom_point(aes(text=paste(NAME,POS,sep="_")), size=0.6) +
  scale_x_continuous(
    breaks = all_breaks,           # All ticks every 100,000
    labels = function(x) {         # Custom labeling function
      ifelse(x %in% labeled_breaks, 
             label_number()(x),    # Label these breaks
             "")                   # Empty label for others
    }
  ) +
  # scale_y_continuous(limits=c(0.1,1.2)) +
  # facet_col(~LIBRARY)+
  scale_y_log10()+
  theme_cowplot(14)+
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 22923638, linetype = "dotted", color = "purple",alpha=0.5) +
  geom_vline(xintercept = 23173638, linetype = "dotted", color = "purple",alpha=0.5) +
  geom_vline(xintercept = 23423638, linetype = "dotted", color = "purple",alpha=0.5) +
  labs(
    x = "Position in the flamenco locus",
    y = "fold change upon IDR tethering",
  ) +
  theme_cowplot(14)+
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +
  labs(
    x = "Position in the flamenco locus",
    y = "fold change upon IDR tethering",
  ) +
  theme(
    legend.position = "none"
  )

p

ggsave("flamenco_silencing.flameenco_inclFlank.pdf", p, width = 5, height = 3)


###############################################################################################
###############################################################################################
###############################################################################################
# ###########################   plot  TE coverage in clusters   ##############################
###############################################################################################
###############################################################################################
###############################################################################################


#read in TE deregulation
DGE = read_tsv("RNAseq_DGE.txt")%>%
  filter(str_detect(gene_id, "TE:"), !str_detect(gene_id, "_AS"))%>%
  separate(gene_id, c("X","TE"), sep = ":", remove = FALSE)
  

###################################################################################################
#read in table
RAW = read_tsv("TEcoverage.merged.ALL.txt")%>%
  separate(`CHR:POS`, c("CHR", "POS"), sep = ":",convert = TRUE)

#select TE to plot
currTE="springer"

#determine max and min coverage for selected TE
maxVAL= RAW %>%
  filter(CHR==currTE)%>%
  pivot_longer(cols = -c(CHR, POS), names_to = "TE", values_to = "coverage")%>%
  summarize(maxVAL = max(coverage))%>%
  pull(maxVAL)
minVAL= RAW %>%
  filter(CHR==currTE)%>%
  pivot_longer(cols = -c(CHR, POS), names_to = "TE", values_to = "coverage")%>%
  summarize(minVAL = min(coverage))%>%
  pull(minVAL)

#determine max coverage for each TE and filter out those with low peak-coverage
max_coverage = RAW %>%
  pivot_longer(cols = -c(CHR, POS), names_to = "TYPE", values_to = "coverage")%>%
  mutate(coverage = abs(coverage))%>%
  group_by(CHR)%>%
  summarize(max_coverage = max(coverage))%>%
  filter(max_coverage > 1000)

med_coverage = RAW %>%
  pivot_longer(cols = -c(CHR, POS), names_to = "TYPE", values_to = "coverage")%>%
  filter(abs(coverage) > 5  )%>%
  mutate(coverage = abs(coverage))%>%
  group_by(CHR)%>%
  summarize(mean_coverage = mean(coverage))%>%
  filter(mean_coverage > 220)

#convert results to binary format
CUTOFF=10
TABLE = RAW %>%
  filter(!str_detect(CHR, "DNAREP"))%>%
  left_join(med_coverage, by = "CHR")%>%
  filter(CHR %in% med_coverage$CHR)%>%
  mutate(dm6.sense.bin = case_when(TEcoverage.dm6.sense.bedgraph > 0 ~ 1,
                                   TRUE ~ 0),
         OSC.sense.bin = case_when(TEcoverage.OSC_r1.01.sense.bedgraph > 0 ~ 1,
                                         TRUE ~ 0),
         dm6.antisense.bin = case_when(abs(TEcoverage.dm6.antisense.bedgraph) > 0 ~ 1,
                                       TRUE ~ 0),
         OSC.antisense.bin = case_when(abs(TEcoverage.OSC_r1.01.antisense.bedgraph) > 0 ~ 1,
                                             TRUE ~ 0),
         sRNA.sense.bin = case_when(TEcoverage.sRNA.sense.bedgraph > mean_coverage*0.1 ~ 1,
                                    TRUE ~ 0),
         sRNA.antisense.bin = case_when(abs(TEcoverage.sRNA.antisense.bedgraph) >  mean_coverage*0.1 ~ 1,
                                        TRUE ~ 0))
         

#########################################################################################
#analyze full TE regions

X=TABLE %>%
  group_by(CHR)%>%
  mutate(
    OSCsense = case_when(
      (sRNA.sense.bin == OSC.sense.bin ) ~ 1,
      TRUE ~ -1),
    OSCantisense = case_when(
      (sRNA.antisense.bin == OSC.antisense.bin ) ~ 1,
      TRUE ~ -1),
    DM6sense = case_when(
      (sRNA.sense.bin == dm6.sense.bin ) ~ 1,
      TRUE ~ -1),
    DM6antisense = case_when(
      (sRNA.antisense.bin == dm6.antisense.bin ) ~ 1,
      TRUE ~ -1)
  )

#deterrmine length of sRNA covered region
sRNAcovered = X %>%
  pivot_longer(cols = contains("TEcoverage.sRNA"), names_to = "TYPE", values_to = "coverage" )%>%
  filter(coverage >CUTOFF)%>%
  group_by(CHR,TYPE)%>%
  summarize(
    LENGTH = n()
  )

#summarize --> my own code
Y = X %>%
  group_by(CHR)%>%
  summarize(
    OSCsense = sum(OSCsense),
    OSCantisense = sum(OSCantisense),
    DM6sense = sum(DM6sense),
    DM6antisense = sum(DM6antisense),
    LENGTHsense=n(),
    LENGTHantisense=n(),
    sRNA_antisense = sum(TEcoverage.sRNA.antisense.bedgraph),
    sRNA_sense = sum(TEcoverage.sRNA.sense.bedgraph)
  )%>%
  mutate(
    OSCsense_AVG = (OSCsense/LENGTHsense+1)*100/2,
    OSCantisense_AVG = (OSCantisense/LENGTHantisense+1)*100/2,
    DM6sense_AVG = (DM6sense/LENGTHsense+1)*100/2,
    DM6antisense_AVG = (DM6antisense/LENGTHantisense+1)*100/2,
    sRNA_antisense_AVG = sRNA_antisense/LENGTHantisense,
    sRNA_sense_AVG = sRNA_sense/LENGTHsense
  )%>%
  #mutate away NaN by -1 in the data frame
  mutate_all(~replace_na(., 0))%>%
  select(-LENGTHsense,-LENGTHantisense)%>%
  arrange(OSCantisense)%>%
  select(CHR,contains("antisense"))


#plot heatmap visualizing the overlap between the different TE coverages

# Create the bar plot
ordered_chr <- Y %>%
  arrange(OSCantisense_AVG) %>%
  .$CHR %>%
  unique() %>%
  factor(levels = .)

bar_plot <- Y%>%
  select(-contains("AVG"),CHR)%>%
  mutate(CHR = factor(CHR, levels = levels(ordered_chr))) %>%
  select(CHR,contains("sRNA"))%>%
  pivot_longer(cols = -c(CHR), names_to = "STRAND", values_to = "VALUE")%>%
  ggplot( aes(x=CHR, y=abs(VALUE))) +
  geom_bar(stat="identity", width = 0.5, fill = "grey", color = "black", size = 0.1) +
  coord_flip() +
  theme_cowplot(14) +
  facet_wrap(~STRAND) +
  theme(plot.margin = margin(l = 0, r = 0, b = 0, t = 0))+
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    axis.title.y = element_blank()
  )

piRNA_heatmap = Y%>%
  select(-contains("AVG"),CHR)%>%
  mutate(CHR = factor(CHR, levels = levels(ordered_chr))) %>%
  select(CHR,contains("sRNA"))%>%
  pivot_longer(cols = -c(CHR), names_to = "STRAND", values_to = "VALUE")%>%
  ggplot(aes(x=STRAND, y=CHR, fill=log10(abs(VALUE))))+
  geom_tile(color="black")+
  scale_fill_acton( reverse = TRUE, range=c(0,1))+
  theme_bw()+
  theme(
        axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())

DGEfiltered = DGE%>%
  filter(TE %in% Y$CHR)%>%
  select(TE, log2FoldChange)%>%
  mutate(TE = factor(TE, levels = levels(ordered_chr)))

maxVAL = DGEfiltered %>%
  summarize(maxVAL = max(log2FoldChange))%>%
  pull(maxVAL)

TE_heatmap = DGEfiltered%>%
  filter(TE %in% Y$CHR)%>%
  select(TE, log2FoldChange)%>%
  mutate(TE = factor(TE, levels = levels(ordered_chr))) %>%
  ggplot(aes(x=1, y=TE, fill=log2FoldChange))+
  geom_tile(color="black")+
  scale_fill_acton( reverse = TRUE, limits=c(0, maxVAL), oob = scales::squish)+
  theme_bw()+
  theme(
        axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())


# Create the heatmap
heatmap <- Y %>%
  mutate(CHR = factor(CHR, levels = levels(ordered_chr))) %>%
  select(contains("AVG"),CHR)%>%
  select(  -contains("sRNA"))%>%
  pivot_longer(cols = -c(CHR), names_to = "TE", values_to = "overlap")%>%
  mutate(STRAND = case_when(
    str_detect(TE,"antisense") ~ "antisense",
    TRUE ~ "sense"
  ))%>%
  mutate(TE = str_remove(TE,"antisense"),
         TE = str_remove(TE,"sense"))%>%
  ggplot(aes(x=TE, y=CHR, fill=overlap))+
  geom_tile(color="black")+
  geom_text(aes(label = round(overlap, digits=2)), size = 2, color = "black") +
  scale_fill_lajolla( reverse = TRUE, range=c(0,1))+
  facet_wrap(~STRAND)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.text.y = element_text(size = 6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "left")


# Combine the plots
combined_plot <- heatmap + bar_plot + plot_layout(widths = c(1,1))

# Print the combined plot
print(combined_plot)
ggsave("TEcoverage_heatmap.fullTE.pdf", combined_plot, width = 8, height = 8)

#########################################################################################
#only analyze regions covered by sRNAs 

#determmine overlap between srna and TE coverage
X=TABLE %>%
  filter(CHR %in% max_coverage$CHR)%>%
  group_by(CHR)%>%
  mutate(
    OSCsense = case_when(
      (sRNA.sense.bin == 1 & OSC.sense.bin == 1) ~ 1,
      (sRNA.sense.bin == 1 & OSC.sense.bin == 0) ~ -1,
                         TRUE ~ 0),
    OSCantisense = case_when(
      (sRNA.antisense.bin == 1 & OSC.antisense.bin == 1) ~ 1,
      (sRNA.antisense.bin == 1 & OSC.antisense.bin == 0) ~ -1,
      TRUE ~ 0),
    DM6sense = case_when(
      (sRNA.sense.bin == 1 & dm6.sense.bin == 1) ~ 1,
      (sRNA.sense.bin == 1 & dm6.sense.bin == 0) ~ -1,
      TRUE ~ 0),
    DM6antisense = case_when(
      (sRNA.antisense.bin == 1 & dm6.antisense.bin == 1) ~ 1,
      (sRNA.antisense.bin == 1 & dm6.antisense.bin == 0) ~ -1,
      TRUE ~ 0)
  )

#deterrmine length of sRNA covered region
sRNAcovered = X %>%
  pivot_longer(cols = contains("TEcoverage.sRNA"), names_to = "TYPE", values_to = "coverage" )%>%
  filter(coverage >CUTOFF)%>%
  group_by(CHR,TYPE)%>%
  summarize(
    LENGTH = n()
  )

#summarize 
Y = X %>% 
  group_by(CHR)%>%
  summarize(
    OSCsense = sum(OSCsense),
    OSCantisense = sum(OSCantisense),
    DM6sense = sum(DM6sense),
    DM6antisense = sum(DM6antisense),
    LENGTHsense=sum(sRNA.sense.bin),
    LENGTHantisense=sum(sRNA.antisense.bin)
  )%>%
  mutate(
    OSCsense = OSCsense/LENGTHsense,
    OSCantisense = OSCantisense/LENGTHantisense,
    DM6sense = DM6sense/LENGTHsense,
    DM6antisense = DM6antisense/LENGTHantisense
  )%>%
  #mutate away NaN by -1 in the data frame
  mutate_all(~replace_na(., 0))%>%
  select(-LENGTHsense,-LENGTHantisense)



Y %>%
  # select(-contains("antisense"))%>%
  select(CHR,contains("antisense"))%>%
  arrange(OSCantisense) %>%
  mutate(CHR = factor(CHR, levels = unique(CHR)))%>%
  pivot_longer(cols = -c(CHR), names_to = "TE", values_to = "overlap")%>%
  filter(CHR==currTE)

#plot heatmap visualizing the overlap between the different TE coverages
p = Y %>%
  arrange(OSCantisense) %>%
  pivot_longer(cols = -c(CHR), names_to = "TE", values_to = "overlap")%>%
  mutate(STRAND = case_when(
    str_detect(TE,"antisense") ~ "antisense",
    TRUE ~ "sense"
  ))%>%
  #remove sense/antisense from the TE column
  mutate(TE = str_remove(TE,"antisense"),
         TE = str_remove(TE,"sense"))%>%
  # select(-contains("antisense"))%>%
  mutate(CHR = factor(CHR, levels = unique(CHR)))%>%
  ggplot(aes(x=TE, y=CHR, fill=overlap))+
    geom_tile()+
    scale_fill_viridis_c()+
    facet_wrap(~STRAND)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          axis.text.y = element_text(size = 6))

ggsave("TEcoverage_heatmap.sRNA_covered_regions.pdf", p, width = 5, height = 8)



#plot coverage of selected TE
for(currTE in med_coverage$CHR){
  print(currTE)
  plotTABLE = TABLE %>%
    filter(CHR==currTE)%>%
    select(CHR,POS,contains("sRNA"),contains(".bin"))

  minVAL= plotTABLE %>%
    summarize(minVAL = min(TEcoverage.sRNA.antisense.bedgraph))%>%
    pull(minVAL)
  maxVAL= plotTABLE %>%
    summarize(maxVAL = max(TEcoverage.sRNA.sense.bedgraph))%>%
    pull(maxVAL)

  plotTABLE = plotTABLE %>%
    mutate(POS = as.numeric(POS),
           OSC.antisense.bin = OSC.antisense.bin*minVAL*1.05,
           OSC.sense.bin = OSC.sense.bin*maxVAL*1.5,
           dm6.antisense.bin = dm6.antisense.bin*minVAL*1.05,
           dm6.sense.bin = dm6.sense.bin*maxVAL*1.5
    )%>%
    arrange(POS)

  collapse_blocks <- function(df) {
    df %>%
      arrange(POS) %>%
      group_by(CHR, TYPE, coverage) %>%
      mutate(
        block_id = cumsum(c(1, diff(POS) != 1))
      ) %>%
      group_by(CHR, TYPE, coverage, block_id) %>%
      summarise(
        POS_start = min(POS),
        POS_end = max(POS),
        .groups = "drop"
      ) %>%
      select(-block_id)
  }
  
  
  # Modified plotting code
  plotTABLE_pivoted <- plotTABLE %>%
    pivot_longer(cols = -c(CHR, POS), names_to = "TYPE", values_to = "coverage")
  
  # Separate and collapse area data (OSC)
  area_data <- plotTABLE_pivoted %>%
    filter(str_detect(TYPE, "OSC")) %>%
    collapse_blocks()
  
  # Keep line data (sRNA) as is
  line_data <- plotTABLE_pivoted %>%
    filter(str_detect(TYPE, "sRNA"))
  
  # Updated plot using geom_rect for blocks
  p <- ggplot() +
    geom_rect(data = area_data %>% filter(str_detect(TYPE, "OSC.antisense.bin")),
              aes(xmin = POS_start, xmax = POS_end, ymin = 0, ymax = coverage),
              fill = "red", alpha = 0.2) +
    geom_rect(data = area_data %>% filter(str_detect(TYPE, "OSC.sense.bin")),
              aes(xmin = POS_start, xmax = POS_end, ymin = 0, ymax = coverage),
              fill = "cyan", alpha = 0.2) +
    geom_line(data = line_data %>% filter(str_detect(TYPE, "sRNA.sense")),
              aes(x = POS, y = coverage), linewidth = 0.1) +
    geom_line(data = line_data %>% filter(str_detect(TYPE, "sRNA.antisense")),
              aes(x = POS, y = coverage), linewidth = 0.1) +
    theme_bw()
  
  ggsave(paste0("TEplots/",currTE,".OSC.pdf"), p,width = 8, height = 4)
  
  # Separate and collapse area data (OSC)
  area_data <- plotTABLE_pivoted %>%
    filter(str_detect(TYPE, "dm6")) %>%
    collapse_blocks()
  
  # Keep line data (sRNA) as is
  line_data <- plotTABLE_pivoted %>%
    filter(str_detect(TYPE, "sRNA"))
  
  p <- ggplot() +
    geom_rect(data = area_data %>% filter(str_detect(TYPE, "dm6.antisense.bin")),
              aes(xmin = POS_start, xmax = POS_end, ymin = 0, ymax = coverage),
              fill = "red", alpha = 0.2) +
    geom_rect(data = area_data %>% filter(str_detect(TYPE, "dm6.sense.bin")),
              aes(xmin = POS_start, xmax = POS_end, ymin = 0, ymax = coverage),
              fill = "cyan", alpha = 0.2) +
    geom_line(data = line_data %>% filter(str_detect(TYPE, "sRNA.sense")),
              aes(x = POS, y = coverage), linewidth = 0.1) +
    geom_line(data = line_data %>% filter(str_detect(TYPE, "sRNA.antisense")),
              aes(x = POS, y = coverage), linewidth = 0.1) +
    theme_bw()
  
  
  ggsave(paste0("TEplots/",currTE,".dm6.pdf"), p,width = 8, height = 4)
}



RAW = read_tsv("OSCs_siPIWI.TEhist") %>% 
  filter(TE=="springer")

RAW %>% 
  ggplot(aes(x=POS, y=VAL))+
  geom_line()
