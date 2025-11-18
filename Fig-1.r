library(tidyverse)
library(dplyr)
library(ggforce)
library(cowplot)
library(plotly)
library(scales)
library(ggbeeswarm)
library(ggrastr)
library(ggh4x)
setwd("~OSC-paper//Fig1/")


###############################################################################################
###############################################################################################
###############################################################################################
# ###########################   plot  SNPs along the chromosome  ##############################
###############################################################################################
###############################################################################################
###############################################################################################

RAWsnp = read_tsv("OSC_r1.01_SNPs_along_chromosome.txt") %>% 
  filter(CHR %in% c("2L_RagTag","2R_RagTag", "3L_RagTag", "3R_RagTag", "X_RagTag", "4_RagTag"), unique100nt == 2 ) %>%
  separate_rows(VAL, sep = ",", convert = TRUE) %>%
  group_by(CHR, POS, POSend,unique100nt) %>%
  summarise(VAL = max(VAL), .groups = "drop")
  {}

p = RAWsnp %>% 
  ggplot(aes(x=as.numeric(POS),y=as.numeric(VAL)))+
  geom_point_rast(size=0.02, alpha=0.1)+
  facet_grid2(~CHR, scales = "free",  space = "free_x")+
  theme_cowplot(14)+
  theme(
    panel.border = element_rect(color="black", fill=NA, linewidth=0.3)
  )
p

ggsave("SNPs_along_chromosome.pdf", p, width=20, height=4)

SNPdata = RAWsnp %>% 
  mutate(
    lengthCLASS = "SNP",
    ALPHA=0.01,
    SIZE=0.1,
    TYPE="SNP",
    TE="noTE",
    TEtype="noTE"
  )%>%
  select(CHR, POS, VAL, lengthCLASS,ALPHA, SIZE,TYPE,TE,TEtype)

###############################################################################################
###############################################################################################
###############################################################################################
# ###########################   plot  SVs along the chromosome  ##############################
###############################################################################################
###############################################################################################
###############################################################################################

RAW = read_tsv("OSC_r1.01_SVs_along_chromosome.txt") %>% 
  filter(CHR %in% c("2L_RagTag","2R_RagTag", "3L_RagTag", "3R_RagTag", "X_RagTag", "4_RagTag") ) %>%
  separate_rows(VAL, sep = ",", convert = TRUE) %>%
  mutate(
    lengthCLASS = case_when(
      abs(SVlength) < 5000 ~ "short",
      TRUE ~ "long"
    ),
    ALPHA=0.5,
    SIZE=0.2,
    TEtype = case_when(
      TE == "noTE" ~ "noTE",
      TE== "NA" ~ "NA",
      TRUE ~ "TE"
    ),
    
  )%>%
{}

#merge tables
MERGED = rbind(RAW %>% select(CHR, POS, VAL, lengthCLASS, SIZE,ALPHA,TYPE,TE,TEtype), SNPdata)

p = MERGED %>% 
  ggplot(aes(x=as.numeric(POS),y=as.numeric(VAL), alpha=ALPHA, size=SIZE))+
  geom_point_rast(aes(color=TEtype))+
  facet_grid2(lengthCLASS~CHR, scales = "free",  space = "free_x")+
  scale_alpha_identity()+
  scale_size_identity()+
  scale_color_manual(values=c("noTE"="black", "NA"="white", "TE"="#cc79a7"))+
  theme_cowplot(14)+
  theme(
    panel.border = element_rect(color="black", fill=NA, linewidth=0.3)
  )
p

ggsave("variations_along_chromosome.pdf", p, width=16, height=6)

###############################################################################################
###############################################################################################
###############################################################################################
# ###########################   plot  length-histogram  ##############################
###############################################################################################
###############################################################################################
###############################################################################################

#length histogram
p=RAW %>% 
  mutate(
    SVlength = abs(SVlength),
    SVlength = case_when(
      SVlength > 12000 ~ 12000,
      TRUE ~ SVlength
    )
  )%>%
  filter(TE!="NA")%>%
  ggplot(aes(x=SVlength, fill=TEtype))+
  geom_histogram(bins=50, color="black", linewidth=0.1)+
  scale_fill_manual(values=c("noTE"="black", "NA"="white", "TE"="red"))+
  theme_cowplot(6)+
  theme(
    panel.border = element_rect(color="black", fill=NA, linewidth=0.3),
  )
p
ggsave("SV_length_histogram.pdf", p, width=8, height=3)

ggplotly(p, tooltip = c("text"))

###############################################################################################
###############################################################################################
###############################################################################################
# ###########################   calculate read and SNP coverage  ##############################
###############################################################################################
###############################################################################################
###############################################################################################

#smoothed Illumina coverage
RAW = read_tsv("Illumina_coveraged.smoothed.bg",
                     col_names = c("CHR", "START", "END", "VAL"),
                     col_types = "ciid")%>%
  filter(CHR %in% c("2L_RagTag","2R_RagTag", "3L_RagTag", "3R_RagTag", "X_RagTag")) 
  

# compute a mid-point for plotting
RAW <- RAW %>%
  mutate(MID = (START + END) / 2)

# plot smoothed histogram/curve
ggplot(RAW, aes(x = MID, y = VAL)) +
  geom_line(color = "steelblue") +
  geom_smooth(se = FALSE, span = 0.00001, color = "red") +
  theme_minimal() +
  scale_y_continuous(limits=c(0,150)) +
  facet_wrap(~CHR, scales="free_y")+
  labs(x = "Genomic position", y = "Signal", 
       title = "Smoothed bedGraph track")


#####################################################################################
#SNPS PER KB

RAW = read_tsv("OSC_r1.01_SNPs_per_1kb.txt")

RAW %>%
  ggplot(aes(x=1, y=het_count))+
  geom_quasirandom(size=0.3)+
  scale_y_log10()


RAW = read_tsv("dm6_SNPs_per_1kb.txt")%>%
  mutate(
    totalSNP = het_count+hom_count
  )

