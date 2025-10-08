library(limma)
library(scales)
#library(DataCombine)
library(colourpicker)
library(gridExtra)
library(gridGraphics)
library(ggpubr)

#load("U:/Marcus/New/all data seperate_joined/out.Rdata")
#setwd("U:/Marcus/New/all data seperate_joined")

col1 <- c("#228b22",#forestgreen
          "#186118", #dark green
          "#4f94cd", #steelblue3
          "#95bfe1", #light blue
          "#666666", #gray40
          "#dbdbdb", #gray86
          "#e12f0a", #red
          "#ed826c", #light red
          "#01e8fe", #aquamarine
          "#018b98", #dark aquamarine
          "#d601e1", #violet
          "#96019e", #dark violet
          "#fe8601", #orange
          "#be4e0f", #brown
          "#000000", #black
          "#ffffff"
) 



Microarrayplot <- function(Geo, Annotation, Gff, start, end, replicon, graph_size = 3,color=c(),
                           subgenes = FALSE, incom_genes ="_partial", min_micro = NA, 
                           max_micro=NA, filter_micro=TRUE, order_micro = NULL,
                           filter_reads =TRUE,
                           rename_graphs=NULL, line_visible=TRUE, arrow_body_height= 7, 
                           arrowhead_height=10, arrowhead_width=8, Gsize =1.2, 
                           msize= 4, Mwidth =7, ntlength =NA, lineSize = 1){
  if(is.null(graph_size)) graph_size <- 3
  if(is.null(subgenes)) subgenes <- TRUE
  if(is.null(line_visible)) line_visible <- TRUE
  if(is.null(incom_genes)) incom_genes <- "_partial"
  if(is.null(arrow_body_height)) arrow_body_height <- 7
  if(is.null(arrowhead_height)) arrowhead_height <- 10
  if(is.null(arrowhead_width)) arrowhead_width <- 8
  if(is.null(max_micro)) max_micro <- NA
  if(is.null(ntlength)) ntlength <- NA
  
  max1 <- max(Annotation$end[which(Annotation$Replicon==replicon)], Gff$end[which(Gff$Replicon==replicon)])
  end <- min(max1+500, end)
  
  
  N <- str_remove(colnames(Geo),  zero_or_more("[-_. ]")%R%one_or_more(DGT)%R%END)
  N1 <- unique(N[!is.na(N)])
  if(is.null(order_micro)) order_micro <- N1
  
  
  #filter
  q <- filter_micro
  if(is.null(q)){q<- N}
  Data <- Geo[,which(N %in% q)]
  N <- N1[which(N1 %in% q)]
  rename_graphs <- rename_graphs[which(N1 %in% q)]
  
  Data2 <- cbind(Data, Annotation)
  Data2$med <- (as.integer(Data2$start)+as.integer(Data2$end))/2
  
  #reduce data to interest
  D <- Data2[which(Data2$start > start & Data2$end < end & Data2$Replicon==replicon), ]
  
  if(nrow(D)==0){D[1,] <- rep(NA, ncol(D))}
  
  n <- 1:(which(colnames(D)=="Row")[1]-1)
  #U <- unique(str_sub(colnames(O1)[n],1,-3))
  Un <- ncol(D)
  for(i in N){
    n1 <- which(str_remove(colnames(D), zero_or_more("[-_. ]")%R%one_or_more(DGT)%R%END)==i)
    D[, ncol(D)+1] <- sapply(1:nrow(D), function(i2)
      log(mean(as.numeric(2^(D[i2,n1]))),base =2))
      #exp(mean(as.numeric(log(D[i2,n1])))))
  }
  colnames(D)[Un+1:length(N)] <- str_c(N,"m")
  
  
  g <- Un+1:length(N)
  U <- unique(D$Sequence)
  for(i in U){
    n1 <- which(D$Sequence ==i)
    if(length(n1)>1){
      
      for(i2 in g){
        D[n1[1], i2] <- exp(mean(as.numeric(log(D[n1, i2]))))
      }
      D <- D[-n1[2:length(n1)],]
    }
  }
  
  
  
  D2 <- D[,which(colnames(D)=="Row"):ncol(D)]
  Op <- D2 %>% pivot_longer(cols = which(str_remove(colnames(D2), "m"%R%END) %in% N), names_to = "run", 
                           values_to = "expression")
  Op$run <- str_remove(Op$run, "m"%R%END)
  
  si <- 13
  
  if(!is.logical(Gff$orientation)){
    test <- c("+", "fwd", "fw", 1, "forward", "Forward", "f","TRUE")
    Gff$orientation <- Gff$orientation %in% test
  }
  
  
  Op$group <- str_c(Op$run,"___", Op$Name)
  
  #Op$end <- droplevels(Op$end)
  Op$num <- 1:nrow(Op)
  #Op$end <- as.numeric(as.character(Op$end))
  Op$run <-  factor(Op$run, levels = order_micro)
  
  Op <- Op[order(match(Op$run, order_micro)), ]
  Op$group <-  factor(Op$group, levels = unique(Op$group))
  
  
  if(is.na(min_micro)){min_micro <- floor(min(Op$expression))-1}
  if(is.na(max_micro)){max_micro <- ceiling(max(Op$expression))+2}
  
  #col <- as.character(sapply(N, function(i) color[which(N1 %in% i)]))
  col <- color
  r <- which(Op$orientation=="+" )
  f <- which(Op$orientation=="-" )
  
  
  
  p1 <- ggplot(Op[f,], aes( color = run))+
    geom_line(aes(x=as.integer(med), group = group, y = expression),
              linewidth = lineSize)+ 
    geom_segment(aes(x = as.integer(start), y = expression, 
                     xend = as.integer(end), yend = expression, group =num),
                 linewidth = lineSize*1.5) +
    theme_bw() +                                                         #deletes background
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.line.x.bottom=element_blank(),
          axis.line.y.left=element_line(linewidth=1),
          axis.text=element_text(size=si, face="bold"),
          legend.box = "horizontal",
          legend.text = element_text(size = 14),
          axis.title = element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y.left = element_line(linewidth=1),
          legend.background = element_blank(),
          plot.background = element_rect(fill='transparent', color=NA))+ #cosmetic
    scale_x_continuous(position="bottom",limits = c(start,end), expand = c(0, 0))+ 
    scale_y_continuous(expand = c(0, 0), limits = c(min_micro, max_micro))+
    labs(color = "")+
    scale_color_manual(values=col)+
    guides(color = guide_legend(override.aes = list(linewidth = 10)))
    
  
  
  p2 <- ggplot(Op[r,], aes( color = run))+
    geom_line(aes(x=as.integer(med), group= group, y = -expression),
              linewidth = lineSize, show.legend = FALSE) + 
    geom_segment(aes(x = as.integer(start), y = -expression, 
                     xend = as.integer(end), yend = -expression, group =num),
                 linewidth = lineSize*1.5, show.legend = FALSE)+
    theme_bw() +                                                         #deletes background
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.line.x.bottom=element_line(linewidth=1),
          axis.line.y.left=element_line(linewidth=1),
          axis.text=element_text(size=si, face="bold"),
          axis.title = element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y.left = element_line(linewidth=1),
          plot.background = element_rect(fill='transparent', color=NA))+ #cosmetic
    scale_x_continuous(position="top",limits = c(start,end), expand = c(0, 0))+ 
    scale_y_continuous(expand = c(0, 0), limits = c(-max_micro, -min_micro)) +
    #labs(color = "")+                                                     #no Name for legend
    scale_color_manual(values=col)       
  
  
  
  
  
  
  
  g0 <- which(((Gff$start>=start & Gff$end<=end) |
                 (Gff$start<start & Gff$end>start) |
                 (Gff$start<end & Gff$end>end)) &
                Gff$Replicon==replicon)#get all relevant Gff-elements
  
  lb <- arrow_body_height #width of gene-arrow
  lh <-arrowhead_height
  lw <-arrowhead_width
  
  if(length(g0)>0){
    G <- Gff[g0,]
    
    s <- which((G$start<start & !G$orientation) |
                 G$end>end & G$orientation)         #which Gff-elements would have an arrow and over the edge
    s1 <- which(G$start<start | G$end>end)             #which Gff-elements are too long
    
    if(length(s1)>0){ G$gene[s1] <- str_c(G$gene[s1], incom_genes) }  #rename elements that are to long
    
    
    G$start <- sapply(G$start, function(i) max(i, start))   #Gff-start is mininmal start of area of Interest
    G$end <- sapply(G$end, function(i) min(i, end))         #Gff-end is maxinmal end of area of Interest
    
    G$start1<-NA                  #get second category for Gff-elemnts with different arrow 
    G$end1<-NA
    G$startall <- G$start         #get category for all Gff-elements for labeling
    G$endall <- G$end
    if(length(s)>0){
      G$start1[s] <- G$start[s]   
      G$end1[s] <- G$end[s]
      G$start[s] <- NA
      G$end[s] <- NA
    }
    
    
    gf <- which(G$orientation) #which Gff-element is forward
    gr <- which(!G$orientation)#which Gff-element is reverse
  }else{
    G <- data.frame(matrix(NA,ncol=ncol(Gff)))
    colnames(G) <- colnames(Gff)
    gf <-NULL
    gr <-NULL
  }
  
  if(length(unique(G$lane))>1){
    g <- unique(G$lane)
    l1 <- length(unique(G$lane))
    r1 <-1:(l1*2)+nrow(G)
    gf <- c(gf,1:(l1)+nrow(G))
    gr <- c(gr,1:l1+nrow(G)+l1)
    G[r1,] <-NA
    G$orientation[r1] <- c(rep(TRUE,l1),rep(FALSE,l1))
    G$lane[r1] <- rep(g,2)
    G$Replicon[r1] <- G$Replicon[1]
    
  }
  
  if(is.na(ntlength)){
    G$gene[(G$endall -G$startall)/(end-start)*Mwidth < strwidth(G$gene, units="inches")*1.5 ] <-""
  }else{
    G$gene[abs(G$endall -G$startall) < ntlength ] <-""
  }
  
  
  if(length(gf)==0){     #in case there is no forward Gff-element...
    
    p3 <- ggplot() +     #write empty plot
      geom_gene_arrow()
    
  }else{
    col1 <- unique(G$color[gf])
    names(col1)<- unique(G$color[gf])
    
    p3 <- ggplot(G[gf,], aes(y = as.character(lane))) +
      geom_gene_arrow(aes(forward = orientation, fill=color,             #write Gene-map
                          xmin = as.integer(start), xmax = as.integer(end)),
                      arrowhead_height = unit(lh, "mm"), 
                      arrowhead_width = unit(max(0,lw), "mm"), 
                      arrow_body_height = unit(lb, "mm"))+     #first Gff-elements inside totally in area with regular arrowheads
      geom_gene_arrow(aes(forward = orientation, xmin = as.integer(start1), 
                          xmax = as.integer(end1), fill=color),
                      arrowhead_height = unit(lb*0.6, "mm"), 
                      arrowhead_width = unit(lb*0.3, "mm"), 
                      arrow_body_height = unit(lb, "mm"))+    #secound Gff-elements with ending outside area with irregular arrowheads
      scale_fill_manual(values=col1) 
  }
  
  if(line_visible){
    p3 <- p3 + theme_genes()
  }else{
    p3 <- p3 + theme_void()
  }
  
  p3 <- p3 + scale_x_continuous(limits = c(start,end), expand = c(0, 0))+
    theme(axis.line.x.bottom=element_line(linewidth=1),
          axis.text=element_text(size=si, face="bold"), #x-axis is plotted
          axis.title = element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.position="none",
          plot.background = element_rect(fill='transparent', color=NA))
  
  if(length(gf)>0){ 
    if(any(!is.na(G$from[gf]))){ 
      if(subgenes){
        p3 <- p3 + geom_subgene_arrow(aes(xmin = as.integer(startall), xmax = as.integer(endall), 
                                          y = as.character(lane), xsubmin = as.integer(from),
                                          xsubmax = as.integer(to)), color="black",
                                      fill=G$subcolor[gf],
                                      arrowhead_height = unit(lb, "mm"), 
                                      arrowhead_width = unit(0, "mm"), 
                                      arrow_body_height = unit(lb, "mm"))
      }}
    p3 <- p3 +  geom_text(aes(x = (startall+endall)/2, label = gene, fontface = label_type),
                          size = msize)
  }
  
  
  
  if(length(gr)==0){                                  #in case there is no reverse Gff-element...
    p4 <- ggplot() +geom_gene_arrow()#write empty plot
    
  }else{
    col1 <- unique(G$color[gr])
    names(col1)<- unique(G$color[gr])
    p4 <- ggplot(G[gr,], aes(y = as.character(10-as.integer(lane)))) + #write Gene-map (reverse)
      geom_gene_arrow(aes(forward = orientation, xmin = as.integer(start),            
                          xmax = as.integer(end), fill = color),
                      arrowhead_height = unit(lh, "mm"), 
                      arrowhead_width = unit(lw, "mm"), 
                      arrow_body_height = unit(lb, "mm")) +    #first Gff-elements inside totally in area with regular arrowheads
      geom_gene_arrow(aes(forward = orientation, xmin = as.integer(start1), 
                          xmax = as.integer(end1), fill = color),
                      arrowhead_height = unit(lb*0.6, "mm"), 
                      arrowhead_width = unit(lb*0.3, "mm"), 
                      arrow_body_height = unit(lb, "mm"))+    #secound Gff-elements with ending outside area with irregular arrowheads
      scale_fill_manual(values=col1)     #colorization
  }
  
  if(line_visible){
    p4 <- p4+ theme_genes()
  }else{
    p4 <- p4+ theme_void()
  }
  
  p4 <- p4+ 
    scale_x_continuous(limits = c(start,end), expand = c(0, 0))+         #x-axis range is exactly as the lineplots
    theme(axis.text=element_text(size=12, face="bold"),    #text of axis
          axis.title = element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.position="none",                          #no legend
          axis.ticks.x = element_blank(),                  
          axis.text.x = element_blank(),
          plot.background = element_rect(fill='transparent', color=NA))+
    easy_remove_x_axis()    #no x-axis
  
  if(length(gr)>0){
    if(any(!is.na(G$from[gr]))){ 
      if(subgenes){
        p4 <- p4 + geom_subgene_arrow(aes(xmin = as.integer(startall), xmax = as.integer(endall), 
                                          y = as.character(10-as.integer(lane)), xsubmin = as.integer(from),
                                          xsubmax = as.integer(to)), color="black",
                                      fill=G$subcolor[gr],
                                      arrowhead_height = unit(lb, "mm"), 
                                      arrowhead_width = unit(0, "mm"),
                                      arrow_body_height=unit(lb, "mm"))
      }
    }
    p4 <- p4+  geom_text(aes(x = (startall+endall)/2, label = gene, fontface = label_type),
                         size = msize)
  }
  
  
  p_lab <- ggplot() +               #write empty plot, with text
    annotate(geom = "text", x = 0, y = 0, label = "                             number of reads", 
             fontface = "bold", angle = 90, size=6) +
    coord_cartesian(clip = "off")+
    theme_void()                    #no background
  
  p_xlab <- ggplot() +               #write empty plot, with text
    annotate(geom = "text", x = mean(c(start,end)), y = 0, label = "position", fontface = "bold", size=4) +
    coord_cartesian(clip = "off")+
    theme_void() 
  
  plegend <- as_ggplot(get_legend(p1))
  
  ap <- p2 %>% 
    insert_top(p4, height=1/graph_size) %>% 
    insert_top(p3, height = 1/graph_size) %>%
    insert_top(p1+ theme(legend.position = "none"), height = 1)   #combine plots with each other, graph_size free variable
  ap <- ap %>% #insert_left(p_lab, width = 0.1)%>%  #combined complete plots with y-axis title
    insert_bottom(p_xlab, height = 0.15)
  
  ap_grob <- grid.grabExpr(print(ap))
  
  ap <- grid.arrange(ap_grob, plegend, ncol = 2, widths=c(1, 0.2))
  
  
  return(ap) #return plot
  
}









load_Gff3 <- function(input){
  G <- ape::read.gff(input)
  G$ID <- str_extract(G$attributes,one_or_more(WRD)%R%DOT%R%DGT)
  G$Gen <- str_extract(G$attributes,"gene-"%R%one_or_more(WRD))
  G$Gen <- str_sub(G$Gen,6,-1)
  G$strand <- sapply(G$strand, function(i) c(1,-1)[which(c("+","-")==i)])
  b <- which(!duplicated(G[,c(4,5,7,11)]) & !is.na(G$Gen))
  G <- G[b,c(1,11,3:5,7)]
  colnames(G) <- c("Replicon", "gene", "type", "start", "end", "orientation")
  
  h <- which(G$start>G$end)#checks values low to "start"; high to "end"
  if(length(h)>0){
    g1 <- sapply(h, function(i) min(G$start[i],G$end[i])) 
    g2 <- sapply(h, function(i) max(G$start[i],G$end[i])) 
    G$start[h] <- g1               #gives the lower value to column "start"
    G$end[h] <- g2                 #gives the higher value to column "end"
  }
  
  G$color <- "#4682B4"
  G$color[which(G$type !="gene")] <- "#D6DA18"
  G$from <- NA
  G$to <- NA
  G$subcolor <- NA
  G$label_type <- "bold"
  G$label_type[which(G$type =="gene")] <- "bold.italic"
  G$lane <- 1
  return(G)
}


load_Gff <- function(input){
  G <- ape::read.gff(input)
  G$gene <- str_match(G$attributes,"ID=\\s*(.*?)\\s*;")[,2]
  colnames(G)[which(colnames(G)=="strand")] <- "orientation"
  colnames(G)[which(colnames(G)=="seqid")] <- "Replicon"
  
  h <- which(G$start>G$end)
  if(length(h)>0){
    g1 <- sapply(h, function(i) min(G$start[i],G$end[i]))
    g2 <- sapply(h, function(i) max(G$start[i],G$end[i]))
    G$start[h] <- g1
    G$end[h] <- g2
  }
  
  G$color <- "#4682B4"
  G$color[which(G$type !="gene")] <- "#D6DA18"
  G$from <- NA
  G$to <- NA
  G$subcolor <- NA
  G$label_type <- "bold"
  G$label_type[which(G$type =="gene")] <- "bold.italic"
  G$lane <- 1
  return(G)
}

load_csv <- function(input){
  G <- read.csv(input, sep=";")
  if(colnames(G)[1]=="X" & all(G[,1]==1:nrow(G))){
    G <- G[,2:ncol(G)]
  }
  o <- colnames(G)
  
  h <- which(G$start>G$end)
  if(length(h)>0){
    g1 <- sapply(h, function(i) min(G$start[i],G$end[i]))
    g2 <- sapply(h, function(i) max(G$start[i],G$end[i]))
    G$start[h] <- g1
    G$end[h] <- g2
  }
  
  if(!any(o %in% "type")){ G$type <- "unknown" }
  if(!any(o %in% "color")){ 
    G$color <- "#4682B4"
    G$color[which(G$type !="gene")] <- "#D6DA18"}
  G$color <- str_replace_all(G$color,SPC,"")
  G$color[is.na(G$color)] <- "#D6DA18"
  
  if(!any(o %in% "from")){ G$from <- NA }
  if(!any(o %in% "to")){ G$to <- NA }
  if(!any(o %in% "subcolor")){ G$subcolor <- NA }
  if(!any(o %in% "label_type")){ G$label_type <- "bold.italic" 
  }else{ G$label_type[is.na(G$label_type)] <- "bold.italic" }
  if(!any(o %in% "lane")){ G$lane <- 1 
  }else{ G$lane[is.na(G$lane)] <- 1 }
  
  return(G)
  
  
}



