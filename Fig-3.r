library(tidyverse)
library(cowplot)
library(ggalluvial)
library(dplyr)
library(ggforce)
library(plotly)
library(scales)
library(ggbeeswarm)
library(khroma)
library(ggh4x)
library(tidyverse)
library(jsonlite)
library(ggrastr)
theme_set(theme_cowplot())

setwd("~OSC-paper//Fig3/")


###############################################################################################
###############################################################################################
###############################################################################################
# ###########################   plot  SV summary   ##############################
###############################################################################################
###############################################################################################
###############################################################################################

RAW = read_tsv("TEsummary.OSC_r1.01.txt", col_names = TRUE)%>%
  mutate(
    TEclass = case_when(
      TE == 'noTE' ~ 'noTE',
      TRUE ~ 'TE'
    )
  )

###############################################################################################
#length profile
p = RAW %>%
  select(TEclass, ID, TYPE)%>%
  group_by(TEclass, TYPE) %>%
  summarise(COUNT = n()) %>%
  ggplot(aes(axis1=TYPE, axis2=TEclass, y = COUNT, fill = TYPE)) +
  geom_alluvium() +
  geom_stratum() +
  scale_y_continuous(breaks = seq(0, 5500, by = 1000)) +
  scale_x_discrete(expand = c(.1, .1))+
  geom_text(stat = 'stratum', aes(label = after_stat(stratum)), size = 3, na.rm = TRUE) +
  scale_fill_manual(
    values = c('insertion' = '#081626', 'deletion' = '#e69f00'),
    name = 'TE type'
  ) +
  labs(
    title = 'TE analysis of SVs in OSC_r1.01',
    x = 'TE class',
    y = 'Number of SVs'
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = 'bottom',
    plot.title = element_text(hjust = 0.5, size = 16, face = 'bold')
  ) 
p
ggsave(p, filename = 'TEsummary.OSC_r1.01.pdf', width = 4, height = 8, dpi = 300)


###############################################################################################
#length histogram including TE splitup
p = RAW %>%
  mutate(
    SVlength = abs(SVlength),
    SVlength = ifelse(SVlength > 12000, 12000, SVlength),
    ALPHA = case_when(
      TEclass == 'noTE' ~ 0.3,
      TRUE ~ 1
    )
  ) %>%
  ggplot(aes(x = SVlength, fill = TYPE, alpha=TEclass)) +
  # Insertions (positive y)
  geom_histogram(
    data = ~filter(., TYPE == 'insertion'),
    bins = 50, color = 'black', linewidth = 0.1, position = 'stack'
  ) +
  # Deletions (negative y)
  geom_histogram(
    data = ~filter(., TYPE == 'deletion'),
    aes(y = -after_stat(count)),
    bins = 50, color = 'black', linewidth = 0.1, position = 'stack'
  ) +
  scale_fill_manual(
    values = c('insertion' = '#081626', 'deletion' = '#e69f00'),
    name = 'TE type'
  ) +
  scale_alpha_manual(values=c('noTE'=0.3, 'TE'=1)) +
  scale_x_log10() +
  theme_cowplot(14) +
  theme(
    panel.border = element_rect(color = 'black', fill = NA, linewidth = 0.3),
  ) +
  labs(
    x = 'SV length (capped at 12 kb)',
    y = 'Count (insertions positive, deletions negative)'
  )
p
ggsave(p, filename = 'SV-histogram.OSC_r1.01.pdf', width = 16, height = 4, dpi = 300)



###############################################################################################
#plot changes per TE including deregulation potential 
for (i in c("full", "filtered")){

  plotTABLE = RAW %>% 
    filter(TE != 'noTE', SVlength > 1000) %>%
    group_by(TE,TYPE) %>%
    summarise(COUNT = n()) %>%
    mutate(
      COUNT = case_when(
        TYPE == 'deletion' ~ -COUNT,
        TRUE ~ COUNT
      )
    )%>%
    group_by(TE) %>%
    mutate(ORDERVAL = max(ifelse(TYPE == 'insertion', COUNT, 0))) %>%
    ungroup() %>%
    mutate(TE = fct_reorder(TE, ORDERVAL, .desc = FALSE))
  
  #filter TEs with <5 SVs
  if(i == "filtered"){
    plotTABLE = plotTABLE %>%
      group_by(TE) %>%
      filter(sum(abs(COUNT)) >= 10) %>%
      ungroup()
  }
  
  p = plotTABLE %>%
    ggplot(aes(x = TE, y = COUNT, fill=TYPE)) +
    geom_bar(stat = 'identity', width = 0.5, position = 'stack') +
    coord_flip() +
    scale_fill_manual(
      values = c('insertion' = '#081626', 'deletion' = '#e69f00'),
      name = 'TE type'
    ) +
    theme_cowplot(14) +
    theme(
      panel.border = element_rect(color = 'black', fill = NA, linewidth = 0.3),
      axis.title.y = element_blank(),
      legend.position.inside = c(0.8, 0.8)
    ) +
    labs(
      x = 'number of SVs containing indicated TE'
    )
  p
  
  DGE = read_tsv('RNAseq_DGE.txt')%>%
    filter(str_detect(gene_id,"TE:"), !str_detect(gene_id,"_AS")) %>%
    separate(gene_id, into = c("X","TE"), sep = ":")
  
  #add order val into DGE based on TE
  DGE = DGE %>%
    left_join(
      plotTABLE %>%
        select(TE, ORDERVAL) %>%
        distinct(),
      by = 'TE'
    ) %>%
    filter(! is.na(ORDERVAL))%>%
    mutate(TE = fct_reorder(TE, ORDERVAL, .desc = FALSE))%>%
    select(TE, log2FoldChange)
    
  
  TE_heatmap = DGE%>%
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
  TE_heatmap
  
  combined_plot <- p + TE_heatmap + plot_layout(widths = c(10,1))
  combined_plot
  
  HEIGHT=plotTABLE$TE%>%
    unique()%>%
    length()
  
  HEIGHT=HEIGHT*2/10
  print(HEIGHT)
   ggsave(combined_plot, filename = paste("TE-splitup.OSC_r1.01.",i,".pdf",sep=""), width = 10, height = HEIGHT, dpi = 300)
}


###############################################################################################
###############################################################################################
###############################################################################################
# ###########################   plot ChIP heatmap analysis  ##############################
###############################################################################################
###############################################################################################
###############################################################################################

RAW <- read_tsv('TEsummary.OSC_r1.01.txt', col_names = TRUE ) %>%
  mutate(
    TEclass = case_when(
      TE == 'noTE' ~ 'noTE',
      TRUE ~ 'TE'
    )
  )


#######
#read H3K9 matrix from plot-step including clusterd categories


# path to your plotHeatmap matrix
FILE = "TEs_H3K9me3_matrix.plotHeatmap_Sienski-2012.txt"

#Parse JSON header for metadata
FIRSTline = readr::read_lines(FILE, n_max = 1)
stopifnot(startsWith(FIRSTline, "@"))
META = jsonlite::fromJSON(sub("^@", "", FIRSTline))

SAMPLElabels = META$sample_labels               # e.g. c("siLuc_H3K9me3","siPiwi_H3K9me3")
SAMPLEbounds = META$sample_boundaries           # e.g. c(0, 1001, 2002)  (zero-based)
GROUPlabels = META$group_labels                # e.g. c("cluster_1","cluster_2")
GROUPbounds = META$group_boundaries            # zero-based row boundaries

BINsize    <- META[["bin size"]][1]            # 10
US    <- META$upstream[1]                 # 5000
DS  <- META$downstream[1]               # 5000
nBINS       <- US/BINsize + 1 + DS/BINsize  # 1001
BINlabels  <- - (US/BINsize) : (DS/BINsize)  # -500..500

#Read the numeric part (skip JSON line)
# Columns: chrom, start, end, name, score, strand, then bins...
RAW <- readr::read_tsv(
  FILE, skip = 1, col_names = FALSE, progress = FALSE,
  col_types = cols(
    X1 = col_character(),  # chrom
    X2 = col_double(),     # start
    X3 = col_double(),     # end
    X4 = col_character(),  # name
    X5 = col_double(),     # score
    X6 = col_character(),  # strand
    .default = col_double()
  )
)

#Locate the bin columns per sample (convert deepTools' 0-based boundaries to 1-based + offset for first 6 columns)
#For sample k: columns are (6 + samp_bounds[k] + 1) ... (6 + samp_bounds[k+1])
SAMPLEbins_cols <- map2(
  SAMPLEbounds[-length(SAMPLEbounds)],
  SAMPLEbounds[-1],
  ~ seq.int(6 + .x + 1, 6 + .y)
)
names(SAMPLEbins_cols) <- SAMPLElabels

#Build per-sample long tables with explicit bin labels
make_long_for_sample <- function(SAMPLEname) {
  cols <- SAMPLEbins_cols[[SAMPLEname]]
  stopifnot(length(cols) == nBINS)
  
  # set bin labels as column names for this sample’s bin matrix
  vals <- RAW[, cols, drop = FALSE]
  colnames(vals) <- as.character(BINlabels)
  
  tibble(
    CHR  = RAW$X1,
    START  = as.integer(RAW$X2),
    END    = as.integer(RAW$X3),
    NAME   = RAW$X4,
    STRAND = RAW$X6
  ) %>%
    bind_cols(vals) %>%
    pivot_longer(
      cols = all_of(as.character(BINlabels)),
      names_to = "BIN",
      values_to = SAMPLEname
    ) %>%
    mutate(BIN = as.integer(BIN))  # -500..500
}

LONGlist = lapply(SAMPLElabels, make_long_for_sample)

# Join the samples side-by-side by region keys + bin
TABLE = reduce(
  LONGlist,
  ~ full_join(.x, .y, by = c("CHR","START","END","NAME","STRAND","BIN"))
)

# --- 5) Attach k-means cluster per region using group_boundaries (deepTools uses 0-based row indexing)
row_zero_idx = seq_len(nrow(RAW)) - 1L
GROUPindex <- findInterval(row_zero_idx, GROUPbounds, rightmost.closed = TRUE)
REGIONcluster = tibble(
  CHR  = RAW$X1,
  START  = as.integer(RAW$X2),
  END    = as.integer(RAW$X3),
  NAME   = RAW$X4,
  STRAND = RAW$X6,
  CLUSTER = factor(GROUPlabels[pmax(1, pmin(GROUPindex, length(GROUPlabels)))],
                   levels = GROUPlabels)
)

# Add cluster to the long tibble
plotTABLE <- TABLE %>%
  filter(! str_detect(NAME,"random"), !str_detect(NAME, "stellate")) %>%
  left_join(REGIONcluster, by = c("CHR","START","END","NAME","STRAND")) %>%
  arrange(CLUSTER)%>%
  mutate(
    POS = BIN * 10
  )%>%
  # mutate(
  #   POS = case_when(
  #     STRAND == "-" ~ POS*-1,
  #     TRUE ~ POS
  #   )
  # )%>%
  select(-c(CHR,START,END,STRAND,BIN)) %>%
  pivot_longer(
    cols = -c(NAME,CLUSTER,POS),
    names_to = "SAMPLE",
    values_to = "COVERAGE"
  )%>%
  separate(NAME, into = c("TE","ID","STRANDoverlap"), sep = ":!:", convert = TRUE, remove  = FALSE)%>%
  group_by(NAME,SAMPLE) %>%                                 # per TE (row)
  mutate(row_mean = mean(COVERAGE, na.rm = TRUE)) %>% 
  ungroup() %>%
  arrange(CLUSTER, row_mean) %>%                        # sort by mean
  mutate(row_order = factor(NAME, levels = unique(NAME)))%>%
  {}


#plot as heatmap
p1 = plotTABLE %>%
  ggplot(aes(x = POS, y = row_order, fill = log10(COVERAGE+0.1))) +
  geom_tile_rast() +
  scale_fill_lajolla(limits = c(0.1, 2), oob = scales::squish) +
  facet_grid2(CLUSTER~SAMPLE, scales = "free_y",  space = "free_y" ) +
  geom_vline(xintercept = 0,  color = "black", size = 0.3)+
  theme_cowplot()+
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

#summary line-plot
p2 = plotTABLE %>%
  group_by(CLUSTER,SAMPLE,POS)%>%
  summarise(MEAN = mean(COVERAGE, na.rm = TRUE),
            SD = sd(COVERAGE, na.rm = TRUE),
            N = n(),
            SE = SD/sqrt(N),
            CI_lower = MEAN - qt(1 - (0.05 / 2), N - 1) * SE,
            CI_upper = MEAN + qt(1 - (0.05 / 2), N - 1) * SE
  ) %>%
  ggplot(aes(x = POS, y = MEAN, color = CLUSTER, fill = CLUSTER)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2, color = NA) +
  geom_vline(xintercept = 0,  color = "black", size = 0.3)+
  facet_row(~SAMPLE) +
  theme_cowplot()+
  theme(
    # axis.text.y = element_blank(),
    # axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.text.x = element_blank()
  )

#arrange plots into one
pp = plot_grid(p2,p1, ncol = 1, rel_heights = c(1,3), align = "v", axis = "lr")
pp
ggsave(pp, filename = "TE_H3K9.pdf", width = 8, height = 20)  



###############################################################################################
###############################################################################################
###############################################################################################
# ###########################   plot TE category analysis  ##############################
###############################################################################################
###############################################################################################
###############################################################################################
#create per TE plot

TEtable = plotTABLE %>%
  select(NAME,CLUSTER, TE) %>% 
  unique()%>%
  group_by(TE)%>%
  mutate(TEcount = n()) %>%
  filter(TEcount > 10) %>%
  group_by(TE,CLUSTER,TEcount)%>%
  summarise(CLcount = n()) %>%
  {}


regTEs1=c("412","mdg1","gypsy5","Tabor",'Stalker2',"297","gtwin","Quasimodo","gypsy","blood")
regTEs2=c("Stalker","Idefix","NOF","Repbase_Transib1","GATE","HMS-Beagle2","X-element","17.6","Fw2","Repbase_Chouto","rover")
regTEs3=c("Max-element","ZAM","opus","Bari1","Cr1a","TART-A","gypsy10_new","F-element")

regTEs1=c(
  "gypsy5",
  "412",
  "mdg1",
  "297",
  "gtwin",
  "Tabor",
  "gypsy",
  "Stalker",
  "Stalker2",
  "NOF",
  "blood",
  "Idefix",
  "Quasimodo",
  "Repbase_Transib1",
  "HMS-Beagle2",
  "17.6",
  "Repbase_FB4",
  "rover"
)
regTEs2=""

TEtable = TEtable %>%
  mutate(
    CLASS = case_when(
      TE %in% regTEs1 ~ "Cat1",
      TE %in% regTEs2 ~ "Cat1",
      TRUE ~ 'non-regulated'
    )
  )

p=TEtable %>% 
  ggplot(aes(x=reorder(TE,-TEcount), y=CLcount,fill=CLUSTER))+
  geom_bar(stat='identity', position = 'stack')+
  facet_grid2(~CLASS, scales = "free",  space = "free_x")+
  theme_cowplot(14)+
  theme(
    axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1),
    strip.background = element_blank(),
    legend.position = "none"
  )
p
ggsave(p, filename = "TE_H3K9_barplot.pdf", width = 5, height = 4)


TEtable = plotTABLE %>%
  filter(TE=="copia")%>%
  mutate(
    CLUSTER= case_when(
      CLUSTER=="cluster_2" ~ "cluster_1",
      TRUE ~ CLUSTER
    )
  )%>%
  select(NAME,CLUSTER, TE,STRANDoverlap) %>% 
  unique()%>%
  group_by(TE,STRANDoverlap)%>%
  mutate(TEcount = n()) %>%
  filter(TEcount > 10) %>%
  group_by(TE,CLUSTER,TEcount,STRANDoverlap)%>%
  summarise(CLcount = n()) %>%
  {}

p = TEtable %>%
  group_by(CLUSTER) %>%
  mutate(Percent = CLcount / sum(CLcount) * 100) %>%
  ggplot(aes(x = CLUSTER, y = Percent, fill = STRANDoverlap)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(y = "Percent")+
  theme_cowplot(14)
p

ggsave(p, filename = "TE_copia_strand_barplot.pdf", width = 4, height = 4)


###############################################################################################
###############################################################################################
###############################################################################################
# ###########################   plot TE small RNA coverage  ##############################
###############################################################################################
###############################################################################################
###############################################################################################
#TE histograms

RAWsense = read_tsv("TE_hist_sense.man.txt") %>% 
  filter(TE %in% c("copia","gypsy"))%>%
  select(-contains("siArmi"))%>%
  mutate(
    STRAND="sense"
  )

RAWantisense = read_tsv("TE_hist_antisense.man.txt") %>% 
  filter(TE %in% c("copia","gypsy"))%>%
  select(-contains("siArmi"))%>%
  mutate(
    STRAND="antisense" 
  )

RAW = rbind(RAWsense, RAWantisense) %>%
  pivot_longer(
    cols = -c(TE,POS,STRAND),
    names_to = "SAMPLE",
    values_to = "COVERAGE"
  )%>%
  group_by(TE,STRAND,POS)%>%
  summarise(MEAN = mean(COVERAGE, na.rm = TRUE),
            SD = sd(COVERAGE, na.rm = TRUE),
            N = n(),
            SE = SD/sqrt(N),
            CI_lower = MEAN - qt(1 - (0.05 / 2), N - 1) * SE,
            CI_upper = MEAN + qt(1 - (0.05 / 2), N - 1) * SE
  ) 


p = RAW %>%
  ggplot(aes(x = POS, y = MEAN, color = STRAND, fill = STRAND)) +
  geom_line(size = 0.2) +
  geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.2, color = NA) +
  geom_hline(yintercept = 0,  color = "black", size = 0.3)+
  scale_y_continuous(limits = c(-10000,15000))+
  scale_color_manual(values=c("sense"="black","antisense"="black"))+
  scale_fill_manual(values=c("sense"="black","antisense"="black"))+
  facet_col(~TE) +
  theme_cowplot()+
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.title.x = element_blank(),
    # axis.ticks.x = element_blank(),
    # axis.line.x = element_blank(),
    # axis.text.x = element_blank()
  )
p
ggsave(p, filename = "TE_hist.sRNA.pdf", width = 4.5, height = 6)

