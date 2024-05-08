library(qqman)
library(ggplot2)
library(dplyr)
library(ggrepel)

setwd("/home/durwa004/stock526/Aim3/GWAS/")  # Set the working directory

# Read PLINK association results
gwa2 = read.table("SCD_CaseControl.pheno.QC.assoc.assoc", header = TRUE)

# Calculate cumulative positions for each chromosome
don <- gwa2 %>%
  group_by(CHR) %>%
  summarise(chr_len=max(BP)) %>%
  mutate(tot=cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
  select(-chr_len) %>%
  left_join(gwa, ., by=c("CHR"="CHR")) %>%
  arrange(CHR, BP) %>%
  mutate(BPcum=BP+tot)

# Calculate axis positions for labeling
axisdf = don %>%
  group_by(CHR) %>%
  summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

# Set significance thresholds
et_sig <- -log10(0.05/8682578)
gw_sig <- -log10(0.05/17182003)

# Create Manhattan plot for PLINK analysis
MP_plink <- ggplot(don, aes(x=BPcum, y=-log10(P))) +
  geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +  
  geom_hline(yintercept = et_sig, linetype = "dashed", color = "blue") + 
  geom_hline(yintercept = gw_sig, linetype = "dashed", color = "red") +
  scale_x_continuous(label = axisdf$CHR, breaks = axisdf$center) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  labs(title = "Plink",
       x = "Genomic Position",
       y = "-log10(p-value)")

# Save the plot as a PDF file
ggsave("Manhattan_plot_plink.pdf", MP_plink, height = 5, width = 10)








library(qqman)
library(ggplot2)
library(dplyr)
library(ggrepel)

setwd("/home/durwa004/stock526/Aim3/GWAS/")  # Set the working directory

# Read PLINK association results
gwa2 = read.table("SCD_CaseControl.pheno.QC.assoc.assoc", header = TRUE)

# Calculate cumulative positions for each chromosome
don <- gwa2 %>%
  group_by(CHR) %>%
  summarise(chr_len=max(BP)) %>%
  mutate(tot=cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
  select(-chr_len) %>%
  left_join(gwa, ., by=c("CHR"="CHR")) %>%
  arrange(CHR, BP) %>%
  mutate(BPcum=BP+tot)

# Calculate axis positions for labeling
axisdf = don %>%
  group_by(CHR) %>%
  summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

# Set significance thresholds
et_sig <- -log10(0.05/8682578)
gw_sig <- -log10(0.05/17182003)

# Create Manhattan plot for PLINK analysis
MP_plink <- ggplot(don, aes(x=BPcum, y=-log10(P))) +
  geom_point(aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +  
  geom_hline(yintercept = et_sig, linetype = "dashed", color = "blue") + 
  geom_hline(yintercept = gw_sig, linetype = "dashed", color = "red") +
  scale_x_continuous(label = axisdf$CHR, breaks = axisdf$center) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  labs(title = "Plink",
       x = "Genomic Position",
       y = "-log10(p-value)")

# Save the plot as a PDF file
ggsave("Manhattan_plot_plink.pdf", MP_plink, height = 5, width = 10)
