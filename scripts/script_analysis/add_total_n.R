library(dplyr)
library(readr)

args.df <- "pellet_VL_R3.tsv"
args.c1 <- "../igcicle/pellet_VL_R3/1_igCicle.tsv"
args.c2 <- "../igcicle/pellet_VL_R3/2_igCicle.tsv"

df <- read.csv(args.df,sep="\t")
c1 <- read.csv(args.c1,sep="\t")
c2 <- read.csv(args.c2,sep="\t")

add_total_n <- function(df, c1,c2) {
    all_ori <- filter(rbind(c1,c2), fwr1!="",cdr1_aa!="", fwr2!="",cdr2_aa!="",fwr3!="",cdr3_aa!="",fwr4!="")
    df$total_n <- nrow(all_ori)-2
    write.table(df, paste0(strsplit(args.df, split = ".", fixed = TRUE)[[1]][1], "_total_n.tsv"),sep="\t", quote = FALSE )
}

add_total_n(df,c1,c2)
