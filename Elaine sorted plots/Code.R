# start with cyto_norm.cts.txt and control_norm.cts.txt

control_cts <- read.table("ctrl_norm.cts.txt", sep="\t", header = TRUE, as.is = TRUE)
cyto_cts <- read.table("cyto_norm.cts.txt", sep="\t", header = TRUE, as.is = TRUE)

#Index Celltype positions
Aria_index = c(3, 7, 11, 15)
Aria_S2_index = c(2, 6, 10, 14)
PA_index = c(4, 8, 12, 16)
TYTO_S2_index = c(5, 9, 13, 17)

#Initialise output DF
Control_Means = NULL

#Loop over rows (CONTROL)
for (i in 1:nrow(control_cts)) {
	row <- control_cts[i,]
  
	gene <- as.character(row$row.names) 
	Aria <- rowMeans(subset(row, select = c(Aria_index)))
	Aria_S2 <- rowMeans(subset(row, select = c(Aria_S2_index)))
	PA <- rowMeans(subset(row, select = c(PA_index)))
	TYTO_S2 <- rowMeans(subset(row, select = c(TYTO_S2_index)))
  
	Control_Means = rbind(Control_Means, data.frame(gene, Aria, Aria_S2, PA, TYTO_S2))
}

write.table(Control_Means, file = "Control_Means.txt", sep="\t", , quote=F, row.names = F)

#Initialise output DF
Cyto_Means = NULL

#Loop over rows (CYTO)
for (i in 1:nrow(cyto_cts)) {
	row <- cyto_cts[i,]
  
	gene <- as.character(row$row.names) 
	Aria <- rowMeans(subset(row, select = c(Aria_index)))
	Aria_S2 <- rowMeans(subset(row, select = c(Aria_S2_index)))
	PA <- rowMeans(subset(row, select = c(PA_index)))
	TYTO_S2 <- rowMeans(subset(row, select = c(TYTO_S2_index)))
  
	Cyto_Means = rbind(Cyto_Means, data.frame(gene, Aria, Aria_S2, PA, TYTO_S2))
}

write.table(Cyto_Means, file = "Cyto_Means.txt", sep="\t", quote=F, row.names = F)

# Mean values for celltypes in Control + Cyto completed. 

Control_Means <- read.table("Control_Means.txt", sep="\t", header = TRUE, as.is = TRUE)
Cyto_Means <- read.table("Cyto_Means.txt", sep="\t",  header = TRUE, as.is = TRUE)

vec = c()
for (i in 1:nrow(Control_Means)){
	row <- Control_Means[i,]
	if ((row$Aria > row$Aria_S2 && row$Aria > row$TYTO_S2) && (row$PA > row$Aria_S2 && row$PA > row$TYTO_S2)){
	gene <- as.character(row$gene)
	vec <- c(vec, gene)
	}
}
write.table(vec, file="Control_S2_negative_list.txt", quote=F, row.names = F)

#Return genes where S2 positive > S2 negative (CYTO)
vec = c()
for (i in 1:nrow(Cyto_Means)){
	row <- Cyto_Means[i,]
	if ((row$Aria > row$Aria_S2 && row$Aria > row$TYTO_S2) && (row$PA > row$Aria_S2 && row$PA > row$TYTO_S2)){
	gene <- as.character(row$gene)
	vec <- c(vec, gene)
	}
}
write.table(vec, file="Cyto_S2_negative_list.txt", quote=F, row.names = F)

#Return genes where S2 positive > S2 negative (Control)
vec = c()
for (i in 1:nrow(Control_Means)){
	row <- Control_Means[i,]
	if ((row$Aria_S2 > row$Aria & row$Aria_S2 > row$PA) && (row$TYTO_S2 > row$Aria & row$TYTO_S2 > row$PA)){
	gene <- as.character(row$gene)
	vec <- c(vec, gene)
	}
}
write.table(vec, file="Control_S2_positive_list.txt", quote=F, row.names = F)


#Return genes where S2 positive > S2 negative (Cyto)
vec = c()
for (i in 1:nrow(Cyto_Means)){
	row <- Cyto_Means[i,]
	if ((row$Aria_S2 > row$Aria & row$Aria_S2 > row$PA) && (row$TYTO_S2 > row$Aria & row$TYTO_S2 > row$PA)){
	gene <- as.character(row$gene)
	vec <- c(vec, gene)
	}
}
write.table(vec, file="Cyto_S2_positive_list.txt", quote=F, row.names = F)

# Now have DF containing gene names only of S2+ > S2- etc..

#subset original df's using 4 output files from previous step

Control_S2_p <- read.table("Control_S2_positive_list.txt", header = T, as.is = T)
Control_S2_n <- read.table("Control_S2_negative_list.txt", header = T, as.is = T)

Cyto_S2_p <- read.table("Cyto_S2_positive_list.txt", header = T, as.is = T)
Cyto_S2_n <- read.table("Cyto_S2_negative_list.txt", header = T, as.is = T)

subset_ctrl_S2_p <- control_cts[control_cts$row.names %in% Control_S2_p$x,]
subset_ctrl_S2_n <- control_cts[control_cts$row.names %in% Control_S2_n$x,]

subset_cyto_S2_p <- cyto_cts[cyto_cts$row.names %in% Cyto_S2_p$x,]
subset_cyto_S2_n <- cyto_cts[cyto_cts$row.names %in% Cyto_S2_n$x,]

# change 'setwd()' and DF in for loop accordingly:

setwd("/Users/barrydigby/Desktop/Galway_Genomics/Fun/Cyto_S2_negative")
gg_df = NULL

for (i in 1:nrow(subset_cyto_S2_n)){
  row <- subset_cyto_S2_n[i,]
  gene <- as.character(row$row.names)
  Aria <- subset(row, select = c(Aria_index))
  Aria_S2 <- subset(row, select = c(Aria_S2_index))
  PA <- subset(row, select = c(PA_index))
  Tyto_S2 <- subset(row, select = c(TYTO_S2_index))
  
  gg_df <- setNames(data.frame("Aria", t(Aria)), c("celltype", "expr"))
  x <- setNames(data.frame("Aria S2", t(Aria_S2)), c("celltype", "expr"))
  y <- setNames(data.frame("PA", t(PA)), c("celltype", "expr"))
  z <- setNames(data.frame("Tyto S2", t(Tyto_S2)), c("celltype", "expr"))
	
	
  #get summary statistics for ggplot:
  #load dplyr, unload plyr!
  df.summary <- gg_df %>%
  group_by(celltype) %>%
  summarise(
    sd = sd(expr),
    expr = mean(expr))
	
  #need to calc sd in gg_df to set ylim height (zoom in on high vals)
  max <- max((gg_df$expr) + 0.5*(gg_df$sd))  
  min <- min(gg_df$expr)
  
  pdf(paste0(gene,".pdf"), height = 3, width = 3)
  
  plot <- ggplot(gg_df, aes(celltype, expr)) +
          geom_bar(stat = "identity", data = df.summary,
                   fill = NA, color = "black", width = 0.5, size = 1.0) +
          geom_jitter(position = position_jitter(0.0),
                      color = "black") + 
          geom_errorbar(aes(ymin = expr, ymax = expr+sd),
                        data = df.summary, width = 0.25, size=0.2) +
          theme_pubr() +
          coord_cartesian(ylim = c(min, max)) +
          scale_y_continuous(expand = expand_scale(mult = c(0, .1))) +
          labs(title = paste(gene), y = "Normalized Expression", x = "cell type") +
          theme(
                plot.title = element_text(size=12, face="bold"),
                axis.title.x = element_text(size=12, face="bold"),
                axis.title.y = element_text(size=12, face="bold")
        )

  print(plot)
  dev.off()

# coord_cartesian required or else ylim = c(min,max) does not work. 
#scale_y_continuous to remove white space at bottom of bar & x-axis
# the ggbarplot worked really well, but was not able to dictate error bars like ggplot was:

gg_df <- rbind(gg_df, x, y, z)
	
  max <- max((gg_df$expr) + 0.5*(gg_df$sd)) 
  min <- min((gg_df$expr)- 0.5*(gg_df$sd))
  
  pdf(paste0(gene,".pdf"), height = 3, width = 3)
  plot <- ggbarplot(gg_df, x = "celltype", y = "expr",
            add = c("mean_sd", "jitter"),
            color = "black",
            ylim = c(min, max),
            position = position_dodge(0.2),
            ggtheme = theme_bw(),
            font.main = "bold",
            font.x = "bold",
            font.y = "bold",
            xlab = "cell type",
            ylab = "Normalized Expression",
            title = paste(gene))
  print(plot)
  dev.off()
}
