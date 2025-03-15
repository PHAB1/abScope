suppressWarnings(library(dplyr))
suppressWarnings(library(alakazam))

# args - Input and output paths
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("É necessário fornecer dois argumentos: [input_fasta] [output_fasta]")
}

df_inp_path <- args[1]
out_dir <- args[2]

# import AIRR
df <- airr::read_rearrangement(df_inp_path) %>%
  filter(!is.na(v_call), !is.na(j_call), !is.na(junction_length)) %>%
  mutate(v_call = sub(",.*", "", v_call),
         j_call = sub(",.*", "", j_call)) %>%
  mutate(v_call = sub("\\*.*", "", v_call),
         j_call = sub("\\*.*", "", j_call))

# Render Plotly

# Assign sorted levels and subset to all IGHV in subject_id/ sample_id
df <- countGenes(df, gene="v_call", mode="gene")

df$seq_freq <- df$seq_freq*100
vdjBarplot_g1 <- ggplot(df, aes_string(x = "gene", y = "seq_freq")) +
  theme_bw() +
  ggtitle("IGH-V-D-J Usage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ylab("Percent of repertoire") +
  xlab("") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  #scale_fill_brewer(palette = "Set1") +
  geom_bar(stat = "identity", position = "dodge", width = 0.7)
      
#plotly::ggplotly(vdjBarplot_g1)
# save pdf by gene file
ggsave(paste0(out_dir,"vdjBarplot_ByGene.pdf"), plot = vdjBarplot_g1, width = 8, height = 6)

# By family
df <- countGenes(df, gene="v_call", mode="family") 
df <- data_summary(df, varname="seq_freq")

df$seq_freq <- df$seq_freq*100
df$sd <- df$sd*100
df$se <- df$se*100

pdf(paste0(out_dir,"vdjBarplot_ByFamily.pdf"), width = 8, height = 6)
ggplot(df, aes(x=gene, y=seq_freq)) + 
  geom_bar(stat="identity", color="black",
    position=position_dodge()) +
  geom_errorbar(aes(ymin=seq_freq, ymax=seq_freq+se), width=.2,
    position=position_dodge(.9)) + 
  theme_bw() +
  ggtitle("IGH-V-D-J Usage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ylab("Percent of repertoire") +
  xlab("") +
  scale_y_continuous(labels = scales::percent_format(scale = 1))
dev.off()
