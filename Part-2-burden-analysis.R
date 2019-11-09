# first set current working as working directory
library("ggplot2")
library("readr")
library("ggrepel")
options(warn=-1)
rm(list = ls())

# +\- bin definition
bin <- 4 # Meaning 9aa window -4<--1-->+4
# Locating all files for the analysis 
files_clinvar <- list.files(path="db/",pattern="*clinvar-hgmd\\.binary", full.names=T, recursive=FALSE)
files_gnomad <- list.files(path="db/",pattern="*gnomad\\.binary", full.names=T, recursive=FALSE)

# Main loop
for(i in 1:length(files_clinvar)){
  # Loading single family file
  data_c <- read_delim(files_clinvar[i], "\t", escape_double = FALSE, trim_ws = TRUE)
  data_g <- read_delim(files_gnomad[i], "\t", escape_double = FALSE, trim_ws = TRUE)
  family.name <- gsub("gnomad\\.binary","",files_gnomad[i])
  family.name.final <- gsub("db/","",family.name)
  
  cont <- 0
  # Defining Columns
  members <- ((ncol(data_c)-4)/2)
  data_f <- data_c[,c(1:(members+2))]
  data_f$Adj_bin_count <- ""
  data_f$DM_adj_bin_count <- ""
  data_f$sum <- ""
  data_f$mpr <- ""
  data_f$`Gene:Disease` <- data_c$`Gene:Disease`
  data_f$fisher.p <- ""
  data_f$or <- ""
  data_f$ci1 <- ""
  data_f$ci2 <- ""
  data_f$adj.p <- ""
  data_f$log.adj.p <- ""
  data_f$fmpr <- "Neutral"
  data_f$proxy <- "NA"
  # Creating counters
  aa_hit_total_c <- sum(data_c[,(ncol(data_c)-1)])
  aa_hit_total_g <- sum(data_g[,(ncol(data_c)-1)])
  # Main bin count
  for (b in seq((1+bin), (nrow(data_c)-bin), by = (1+bin))) {
    cont = cont + 1
    aa_hit_bin_c <- sum(data_c[c((b-bin):(b+bin)),(ncol(data_c)-1)])
    aa_hit_bin_g <- sum(data_g[c((b-bin):(b+bin)),(ncol(data_c)-1)])
    aa_hit_outbin_c <- aa_hit_total_c - aa_hit_bin_c
    aa_hit_outbin_g <- aa_hit_total_g - aa_hit_bin_g
    # Contingency table and fisher test 
    cont_table_g = matrix(c(aa_hit_bin_c, aa_hit_bin_g, aa_hit_outbin_c, aa_hit_outbin_g), nrow = 2, dimnames = list(c("clinvar", "gnomad"), c("Hit-in", "Hit-out")))
    f_g <- fisher.test(cont_table_g, conf.level = 0.95, alternative = "greater")
    ci_g <- as.vector(f_g$conf.int)
    data_f$Adj_bin_count[b] <- aa_hit_bin_g/members
    data_f$DM_adj_bin_count[b] <- aa_hit_bin_c/members
    data_f$sum[b] <- ((aa_hit_bin_g/members) - (aa_hit_bin_c/members))
    data_f$mpr[b] <- ifelse((aa_hit_bin_g/members) < (aa_hit_bin_c/members), ((aa_hit_bin_g/members) - (aa_hit_bin_c/members)), "Neutral") 
    data_f$fisher.p[b] <- f_g$p.value
    data_f$or[b] <- f_g$estimate  
    data_f$ci1[b] <- ci_g[1]
    data_f$ci2[b] <- ci_g[2]
    data_f$proxy[b] <- "proxy"
  }
  # Multiple Testing adjustment
  pg <- as.vector(t(data_f$fisher.p))
  data_f$adj.p <- p.adjust(pg, method="bonferroni")
  data_f$log.adj.p <- log10(data_f$adj.p)
  # Extracting main stats
  bin_f <- 0
  bin_or <- 0
  cont2 <- 0
  for (v in seq((1+bin), (nrow(data_g)-bin), by = (1+bin))) {  
    if (data_f$adj.p[v]<0.05) {
      resta <- v - bin_f
      if (resta == (1+bin) & bin_or > as.numeric(data_f$or[v])) {
        start <- v
      } else {
        start <- v- bin
      }
      for (a in start:(v+bin)) {
        data_f$fmpr[a] <- "PER"
        data_f$fisher.p[a]     <- data_f$fisher.p[v]
        data_f$or[a]           <- data_f$or[v]
        data_f$ci1[a]          <- data_f$ci1[v]
        data_f$ci2[a]          <- data_f$ci2[v]
        data_f$adj.p[a]        <- data_f$adj.p[v]
        data_f$log.adj.p[a]        <- data_f$log.adj.p[v]
      }
      bin_f <- v
      bin_or <- as.numeric(data_f$or[v])
      data_f$proxy[v] <- "proxy-tag"
      cont2 = cont2 + 1
    } # loop the aa inside the bin
  } # loop the bin
  
  # Combining overlapping PER bins 
  data_f$mpr.tag <- "NA"
  data_f$ini <- 0
  data_f$fin <- 0
  data_f$size <- 0
  if (cont2 > 0) {
    w1<-subset(data_f, subset = proxy == "proxy-tag", select = c(Index, or, ci1, adj.p))
    w1$tags <- "1"
    w1$ini <- 0
    w1$fin <- 0
    tag <- 1
    if (nrow(w1)>1) {
      for (m in seq(2,(nrow(w1)))) {
        dif <- w1$Index[m] - w1$Index[m-1]
        if(dif > (1+bin)) {
          tag = tag + 1
          w1$tags[m] <- tag
        } else {
          w1$tags[m] <- tag     
        }
      }
      w2 <- w1[order(-as.numeric(w1$or)),]   
      w3 <- w2[!duplicated(w2[,c('tags')]),] 
      mpr <- w3[order(as.numeric(w3$Index)),]
      mpr$ini <- 0
      mpr$fin <- 0
      mpr$size <- 0
      mpr$name <- ""
      # Identifing Singular PERs
      for (t in seq(1,tag)) {
        p <- subset(w1, subset = tags == t, select = c(Index))
        mpr$ini[t] <- (min(p)-bin)
        mpr$fin[t] <- (max(p)+bin)
        mpr$size[t] <- (mpr$fin[t] - mpr$ini[t] + 1)
        mpr$name[t] <- paste0("PER",t,sep = "")
        pos <- mpr$Index[t]
        data_f$mpr.tag[pos] <- paste("PER",t,sep = "")
        data_f$ini[pos] <- mpr$ini[t]
        data_f$fin[pos] <- mpr$fin[t]
        data_f$size[pos] <- mpr$fin[t] - mpr$ini[t] + 1
      }
    } else if (nrow(w1) == 1) {
      pos = w1$Index[1]
      data_f$mpr.tag[pos] <- paste("PER1")
      data_f$ini[pos] <- w1$Index[1] - bin
      data_f$fin[pos] <- w1$Index[1] + bin
      data_f$size[pos] <- 2*bin + 1
    }
  }
  # Writing Results
  rm(p,w1,w2,w3,dif,m,t,tag,pg,aa_hit_bin_c,aa_hit_bin_g,aa_hit_outbin_c,aa_hit_outbin_g,bin_f,bin_or,cont2,cont)
  data_f <- data_f[, !(colnames(data_f) %in% c("sum","mpr","foo"))]
  write.table(x = data_f, paste0(getwd(),family.name.final,"bin9.stats"), sep="\t", quote=F, row.names=F, col.names=T)
} #loop de archivos

############# PLOT

# Load
files <- list.files(pattern="*\\.bin9\\.stats$", full.names=T, recursive=FALSE)

for(i in 1:length(files)){
  # Main loop
  #file <- gsub("./","", files[i])
  file <- files[i]
  familia <- gsub(".bin9.stats","", file)
  ## BURDEN
  dat<- read_delim(paste0(file), "\t", 
                   escape_double = FALSE, trim_ws = TRUE) #choose the Info
  {
    p1<-ggplot(dat, aes(x=Index)) +
      # Colors
      geom_line(aes(y=Adj_bin_count, col="gnomAD"),na.rm = T) +
      geom_line(aes(y=DM_adj_bin_count, col="Pathogenic"),na.rm = T) +
      geom_area(aes(y=Adj_bin_count), fill="turquoise3", alpha=0.5) +
      geom_area(aes(y=DM_adj_bin_count), fill="magenta3", alpha=0.35) +
      # Fixed values
      geom_text(aes( 0, 0, label = "Adjusted p-value:", vjust = 1, hjust = 0), size = 3, colour = "grey") +
      geom_text(aes( 0, -1.30103, label = "0.05", vjust = 1, hjust = 0), size = 3, colour = "grey") +
      geom_text(aes( 0, -2, label = "0.01", vjust = 1, hjust = 0), size = 3, colour = "grey") +
      geom_text(aes( 0, -3, label = "0.001", vjust = 1, hjust = 0), size = 3, colour = "grey") +
      geom_text(aes( 0, -4, label = "0.0001", vjust = 1, hjust = 0), size = 3, colour = "grey") +
      geom_text(aes( 0, -5, label = "0.00001", vjust = 1, hjust = 0), size = 3, colour = "grey") +
      geom_text(aes( 0, -6, label = "0.000001", vjust = 1, hjust = 0), size = 3, colour = "grey") +
      
      geom_hline(aes(yintercept=(-1), alpha=0.25), colour="darkgrey", linetype="dashed", show.legend = F)  +
      geom_hline(aes(yintercept=(-2), alpha=0.25), colour="darkgrey", linetype="dashed", show.legend = F)  +
      geom_hline(aes(yintercept=(-3), alpha=0.25), colour="darkgrey", linetype="dashed", show.legend = F)  +
      geom_hline(aes(yintercept=(-4), alpha=0.25), colour="darkgrey", linetype="dashed", show.legend = F)  +
      geom_hline(aes(yintercept=(-5), alpha=0.25), colour="darkgrey", linetype="dashed", show.legend = F)  +
      geom_hline(aes(yintercept=(-6), alpha=0.25), colour="darkgrey", linetype="dashed", show.legend = F)  +
      
      scale_y_continuous(breaks = c(0, seq(0, max(dat$Adj_bin_count, na.rm = T), by = 2))) +
      scale_color_manual(name=c("General Population (gnomAD)","Patient variants (ClinVar/HGMD)","Pathogenic enriched region (PER)"),
                         breaks=c("gnomAD","Pathogenic","PER"),
                         values=c("turquoise3","magenta3","indianred1")) +
      theme_bw() + theme(panel.border = element_rect(fill = NA, colour = "darkgrey", size = 0.1),
                         legend.box.background = element_rect(),
                         legend.title = element_blank(),
                         legend.position = c(.99, .99),
                         legend.justification = c("right", "top"),
                         legend.box.just = "right",
                         legend.margin = margin(6, 6, 6, 6),
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank()) +
      labs(title=paste0("Missense burden analysis: ",familia),
           x = "Aligned aminoacid sequence",
           y = "Missense burden") +
      coord_cartesian(ylim = c((min(dat$log.adj.p,na.rm = T)-0.25),
                               max(dat$Adj_bin_count,na.rm = T)+0.25), expand = FALSE) +
      
      # p-values
      geom_line(aes(y=(ifelse(fmpr=="PER", log.adj.p,0)), col="PER"), na.rm = T) +
      geom_bar(aes(y=(ifelse(fmpr=="PER", log.adj.p,0)),x=dat$Index), fill="indianred1", stat="identity") +
      
      # PERs labels
      geom_label_repel(aes(y=dat$log.adj.p, x = Index, label=mpr.tag),
                       nudge_y = -1, size=2, fontface=  "plain", segment.color = "grey",
                       arrow = arrow(length = unit(0.01, 'npc'), ends = "first")) +
      geom_hline(aes(yintercept=0), colour="darkgrey", linetype="solid")
  }
  ggsave(filename = paste0(getwd(),"/",file,".pdf"),plot = p1,device = "pdf")
}
print("Done!")




