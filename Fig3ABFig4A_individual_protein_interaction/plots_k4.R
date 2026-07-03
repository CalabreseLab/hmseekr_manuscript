
##############
# get the top5RBP list and plot whether they are sig or not
# for p 0.05 0.01 and 0.001
# plot XNM all together
###############################


setwd("XNM")
library(dplyr)
library(tidyr)
library(circlize)
library(ComplexHeatmap)

filelist<-dir(path='./',pattern='XNM_hmseekrhits_paired_unspliced_pval_CKnormed_p.*.csv',full.names = F)

for (fname in filelist) {
  
  df<-read.csv(fname,header=T)
  
  df$group<-substr(df$query, 1, 1)
  df$group<-factor(df$group,levels=c('X','N','M'))
  
  gcat<-df$group
  
  dfrbp<-df[,c(1,3)]
  
  
  # Split the 'top5' column into multiple columns
  dfrbp <- dfrbp %>%
    mutate(top5_split = strsplit(as.character(top5RBP), ";")) %>%
    unnest_wider(top5_split, names_sep = "_") %>%
    replace(is.na(.), "")  # Replace NAs with empty strings
  
  # Rename columns to desired names
  colnames(dfrbp) <- c("query","top5RBP", "RBP1", "RBP2", "RBP3", "RBP4", "RBP5")
  
  dfrbp<-as.matrix(dfrbp[,c(3:7)])
  
  dfpval<-df[,c(1,4)]
  
  # Split the 'top5' column into multiple columns
  # use adjpval
  dfpval <- dfpval %>%
    mutate(pval_split = strsplit(as.character(top5RBP_adjpval), ";")) %>%
    unnest_wider(pval_split, names_sep = "_") %>%
    mutate(across(starts_with("pval_split"), ~ as.numeric(na_if(., "NA"))))  # Convert to numeric
  
  # Rename columns to desired names
  colnames(dfpval) <- c("query","top5RBP_pval", "pval1", "pval2", "pval3", "pval4", "pval5")
  
  dfpval<-as.data.frame(dfpval)
  
  dfpval$query<-gsub('IST|ISTint|EAT1chk|ALAT1chk','',dfpval$query)
  
  rownames(dfpval)<-dfpval$query
  
  dfpval<-as.matrix(dfpval[,c(3:7)])
  
  
  ###############
  # set min pval
  xsmall<-1e-310
  
  col_fun <- function(x) {
    
    ifelse(x > log10(0.05),
           "white",
           # If it is ≤ 0.05, we interpolate from 0 to 0.05
           circlize::colorRamp2(breaks = c(log10(xsmall), log10(0.05)),
                                colors = c("#ff0000", "#ffd8d8"))(x))
  }
  
  dfpval[which(dfpval<(xsmall))]<-xsmall
  dfpvallog10<-log10(dfpval)
  
  # Create the heatmap with p-values and label blocks with RBP names
  xht<-Heatmap(
    dfpvallog10,
    name = "adj_p\n",
    col = col_fun,  
    show_heatmap_legend = TRUE,  # Hide legend since it's binary
    cluster_rows = FALSE,         # Disable clustering if not needed
    cluster_columns = FALSE,
    show_column_names = FALSE,
    row_names_side = "left",
    column_names_gp = gpar(fontsize = 30),
    row_names_gp = gpar(fontsize = 30, fontface = 3),
    row_split = gcat,
    row_gap = unit(15, "mm"),
    row_title = NULL,
    heatmap_legend_param = list(
      at = c(log10(xsmall), log10(0.05)),              
      labels = c("0", "0.05"),
      border = "black",
      labels_gp = gpar(fontsize = 30),
      title_gp = gpar(fontsize = 30),
      legend_direction = "vertical",
      legend_height = unit(10, "cm"),
      grid_width = unit(1, "cm")
    ),
    cell_fun = function(j, i, x, y, width, height, fill) {
      value <- dfpvallog10[i, j]
      
      # Set fill color based on col_fun, cap values above 0.05 as white
      fill_color <- ifelse(is.na(value), "white", col_fun(value))
      
      # Draw the rectangle with the appropriate fill color
      grid.rect(x, y, width, height, gp = gpar(fill = fill_color, col = "#e0e0e0", lwd = 3))
      
      # Add RBP label in each cell, skip label for NA values
      if (!is.na(value)) {
        grid.text(dfrbp[i, j], x, y, gp = gpar(fontsize = 30))
      }
    }
  )
  
  fp<-gsub('XNM_hmseekrhits_paired_unspliced_pval_CKnormed_','',fname)
  fp<-gsub('.csv','',fp,fixed=T)
  
  pdf(paste0('XNM_top5RBPs_heatmap_unsplice_',fp,'.pdf'),width=15,height=28)
  draw(xht,heatmap_legend_side = "right",padding = unit(c(3, 3, 3, 3), "mm"))
  dev.off()
  
}



###################################
# plot total hit counts for XNM all together

setwd("XNM")
library(dplyr)
library(tidyr)
library(circlize)
library(ComplexHeatmap)

t05<-read.csv('XNM_hmseekrhits_paired_unspliced_pval_CKnormed_p05.csv',header=T)
t01<-read.csv('XNM_hmseekrhits_paired_unspliced_pval_CKnormed_p01.csv',header=T)
t001<-read.csv('XNM_hmseekrhits_paired_unspliced_pval_CKnormed_p001.csv',header=T)

sum(t05$query!=t01$query)
sum(t05$query!=t001$query)

comb<-data.frame(query=t05$query,
                 p0.05=t05$hit_totalcount,
                 p0.01=t01$hit_totalcount,
                 p0.001=t001$hit_totalcount)

comb$group<-substr(comb$query, 1, 1)

comb$query<-gsub('IST|ISTint|EAT1chk|ALAT1chk','',comb$query)

xheaders <- c('X1','XrA','XrF','X2','XrB1','X3','XrB2',
              'X4','X5','XrD','X6','X7','X8','X9',
              'XrE','X10','X11','X12','X13','X14') 

nheaders<-paste0('N',c(1:22))
mheaders<-paste0('M',c(1:9))

xnm<-c(xheaders,nheaders,mheaders)

comb_plot<- comb %>% 
  mutate(query=factor(query,levels=xnm)) %>%
  arrange(query)

comb_plot$group<-factor(comb_plot$group,levels=c('X','N','M'))
gcat<-comb_plot$group
comb_plot$group<-NULL

combdata<-as.matrix(comb_plot[,c(2:4)])

rownames(combdata)<-comb_plot$query

hbreaks<-c(0, 2000, 20000, 200000)
hlabels<-c('0', '2,000', '20,000', '200,000')

col_fun <- colorRamp2(
  breaks = hbreaks,
  colors = c("white","#77dd1c","#FCFC03","#dd1c77")
)


# Create the heatmap with p-values and label blocks with RBP names
xnm<-Heatmap(
  combdata,
  name = "totalhits\n",
  col = col_fun,  
  show_heatmap_legend = TRUE,  # Hide legend since it's binary
  cluster_rows = FALSE,         # Disable clustering if not needed
  cluster_columns = FALSE,
  show_column_names = FALSE,
  column_names_side = "top",
  row_names_side = "left",
  column_names_gp = gpar(fontsize = 30),
  column_names_rot = 0,
  column_names_centered = TRUE,
  row_names_gp = gpar(fontsize = 30,fontface = 3),
  row_split = gcat,
  row_gap = unit(15, "mm"),
  row_title = NULL,
  heatmap_legend_param = list(
    at = hbreaks,              
    labels = hlabels,
    border = "black",
    labels_gp = gpar(fontsize = 30),
    title_gp = gpar(fontsize = 30),
    legend_direction = "vertical",
    legend_height = unit(10, "cm"),
    grid_width = unit(1, "cm")
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    value <- combdata[i, j]
    
    # Set fill color based on col_fun, cap values above 0.05 as white
    fill_color <- col_fun(value)
    
    # Draw the rectangle with the appropriate fill color
    grid.rect(x, y, width, height, gp = gpar(fill = fill_color, col = "#e0e0e0", lwd = 3))
    
    # Add RBP label in each cell, skip label for NA values
    if (!is.na(value)) {
      label <- prettyNum(value, big.mark = ",", scientific = FALSE, preserve.width = "none")
      grid.text(label, x, y, gp = gpar(fontsize = 30))
    }
  }
)


pdf('XNM_totalcount_heatmap_unsplice_allps.pdf',width=9,height=28)
draw(xnm,heatmap_legend_side = "right",padding = unit(c(3, 3, 3, 3), "mm"))
dev.off()

###################################
# plot total transcript counts for XNM all together

setwd("XNM")
library(dplyr)
library(tidyr)
library(circlize)
library(ComplexHeatmap)

t05<-read.csv('XNM_hmseekrhits_paired_unspliced_pval_CKnormed_p05.csv',header=T)
t01<-read.csv('XNM_hmseekrhits_paired_unspliced_pval_CKnormed_p01.csv',header=T)
t001<-read.csv('XNM_hmseekrhits_paired_unspliced_pval_CKnormed_p001.csv',header=T)

sum(t05$query!=t01$query)
sum(t05$query!=t001$query)

comb<-data.frame(query=t05$query,
                 p0.05=t05$transcript_totalcount,
                 p0.01=t01$transcript_totalcount,
                 p0.001=t001$transcript_totalcount)

comb$group<-substr(comb$query, 1, 1)

comb$query<-gsub('IST|ISTint|EAT1chk|ALAT1chk','',comb$query)

xheaders <- c('X1','XrA','XrF','X2','XrB1','X3','XrB2',
              'X4','X5','XrD','X6','X7','X8','X9',
              'XrE','X10','X11','X12','X13','X14') 

nheaders<-paste0('N',c(1:22))
mheaders<-paste0('M',c(1:9))

xnm<-c(xheaders,nheaders,mheaders)

comb_plot<- comb %>% 
  mutate(query=factor(query,levels=xnm)) %>%
  arrange(query)

comb_plot$group<-factor(comb_plot$group,levels=c('X','N','M'))
gcat<-comb_plot$group
comb_plot$group<-NULL

combdata<-as.matrix(comb_plot[,c(2:4)])

rownames(combdata)<-comb_plot$query

hbreaks<-c(0, 1000, 3000, 6000)
hlabels<-c('0', '1,000', '3,000', '6,000')

col_fun <- colorRamp2(
  breaks = hbreaks,
  colors = c("white","#77dd1c","#FCFC03","#dd1c77")
)


# Create the heatmap with p-values and label blocks with RBP names
xnm<-Heatmap(
  combdata,
  name = "totalhits\n",
  col = col_fun,  
  show_heatmap_legend = TRUE,  # Hide legend since it's binary
  cluster_rows = FALSE,         # Disable clustering if not needed
  cluster_columns = FALSE,
  show_column_names = FALSE,
  column_names_side = "top",
  row_names_side = "left",
  column_names_gp = gpar(fontsize = 30),
  column_names_rot = 0,
  column_names_centered = TRUE,
  row_names_gp = gpar(fontsize = 30,fontface = 3),
  row_split = gcat,
  row_gap = unit(15, "mm"),
  row_title = NULL,
  heatmap_legend_param = list(
    at = hbreaks,              
    labels = hlabels,
    border = "black",
    labels_gp = gpar(fontsize = 30),
    title_gp = gpar(fontsize = 30),
    legend_direction = "vertical",
    legend_height = unit(10, "cm"),
    grid_width = unit(1, "cm")
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    value <- combdata[i, j]
    
    # Set fill color based on col_fun, cap values above 0.05 as white
    fill_color <- col_fun(value)
    
    # Draw the rectangle with the appropriate fill color
    grid.rect(x, y, width, height, gp = gpar(fill = fill_color, col = "#e0e0e0", lwd = 3))
    
    # Add RBP label in each cell, skip label for NA values
    if (!is.na(value)) {
      label <- prettyNum(value, big.mark = ",", scientific = FALSE, preserve.width = "none")
      grid.text(label, x, y, gp = gpar(fontsize = 30))
    }
  }
)


pdf('XNM_transcript_totalcount_heatmap_unsplice_allps.pdf',width=9,height=28)
draw(xnm,heatmap_legend_side = "right",padding = unit(c(3, 3, 3, 3), "mm"))
dev.off()

