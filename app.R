library(shiny)
library(shinydashboard)
library(keras)
library(DT)  # For displaying a data table
library(magrittr) 
library(Matrix) 
library(shinyhelper) 
library(data.table)
library(tensorflow)
library(reticulate)
library(hdf5r) 
library(ggdendro) 
library(gridExtra) 
library(shinythemes)
library(ggplot2)
library(ggrepel)
library(bslib)
library(shinyjs)
library(Seurat)
library(rmarkdown)
library(shinyalert)
library(stringr)
library(isoband)

genes_neural <- read.csv("genes_neural.csv")
genes_neural <- genes_neural$x

my_matrix_filled_merged <- read.csv("my_matrix_filled_merged.csv")
rownames(my_matrix_filled_merged) <- my_matrix_filled_merged$X
my_matrix_filled_merged <- my_matrix_filled_merged[,-1]

light <- bs_theme(bootswatch = "cerulean")
dark <- bs_theme(bootswatch = "cyborg")
    
sc1conf = readRDS("sc1conf.rds")
sc1def  = readRDS("sc1def.rds")
sc1gene = readRDS("sc1gene.rds")
sc1meta = readRDS("sc1meta.rds")

sc1conf$UI[sc1conf$UI == "orig.ident"] <- "Multiple Myeloma Patient"
sc1conf$ID[sc1conf$ID == "orig.ident"] <- "Multiple Myeloma Patient"
sc1conf$UI[sc1conf$UI == "sample_id"] <- "Treatment Response"
sc1conf$ID[sc1conf$ID == "sample_id"] <- "Treatment Response"
sc1conf$UI[sc1conf$UI == "seurat_clusters"] <- "Seurat Predicted Clusters"
sc1conf$ID[sc1conf$ID == "seurat_clusters"] <- "Seurat Predicted Clusters"

sc1def$meta1 <- "Multiple Myeloma Patient"
sc1def$meta2 <- "Seurat Predicted Clusters"
sc1def$grp1 <- "Multiple Myeloma Patient"
sc1def$grp2 <- "Seurat Predicted Clusters"

colnames(sc1meta)[2] <- "Multiple Myeloma Patient"
colnames(sc1meta)[6] <- "Seurat Predicted Clusters"
colnames(sc1meta)[7] <- "Treatment Response"

sc1def$genes <- as.character(sc1def$genes)
sc1def$genes <- genes_neural[c(1:5, 66:70)]


### Useful stuff 
# Colour palette, etc
#####
cList = list(c("grey85","#FFF7EC","#FEE8C8","#FDD49E","#FDBB84", 
               "#FC8D59","#EF6548","#D7301F","#B30000","#7F0000"), 
             c("#4575B4","#74ADD1","#ABD9E9","#E0F3F8","#FFFFBF", 
               "#FEE090","#FDAE61","#F46D43","#D73027")[c(1,1:9,9)], 
             c("#FDE725","#AADC32","#5DC863","#27AD81","#21908C", 
               "#2C728E","#3B528B","#472D7B","#440154")) 
names(cList) = c("White-Red", "Blue-Yellow-Red", "Yellow-Green-Purple") 

# Panel sizes 
pList = c("400px", "600px", "800px") 
names(pList) = c("Small", "Medium", "Large") 
pList2 = c("500px", "700px", "900px") 
names(pList2) = c("Small", "Medium", "Large") 
pList3 = c("600px", "800px", "1000px") 
names(pList3) = c("Small", "Medium", "Large") 
sList = c(18,24,30) 
names(sList) = c("Small", "Medium", "Large") 
lList = c(5,6,7) 
names(lList) = c("Small", "Medium", "Large") 

# Function to extract legend 
g_legend <- function(a.gplot){  
    tmp <- ggplot_gtable(ggplot_build(a.gplot))  
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")  
    legend <- tmp$grobs[[leg]]  
    legend 
}  

# Plot theme 
sctheme <- function(base_size = 24, XYval = TRUE, Xang = 0, XjusH = 0.5){ 
    oupTheme = theme( 
        text =             element_text(size = base_size, family = "Helvetica"), 
        panel.background = element_rect(fill = "white", colour = NA), 
        axis.line =   element_line(colour = "black"), 
        axis.ticks =  element_line(colour = "black", size = base_size / 20), 
        axis.title =  element_text(face = "bold"), 
        axis.text =   element_text(size = base_size), 
        axis.text.x = element_text(angle = Xang, hjust = XjusH), 
        legend.position = "bottom", 
        legend.key =      element_rect(colour = NA, fill = NA) 
    ) 
    if(!XYval){ 
        oupTheme = oupTheme + theme( 
            axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
            axis.text.y = element_blank(), axis.ticks.y = element_blank()) 
    } 
    return(oupTheme) 
} 

### Common plotting functions 
# Plot cell information on dimred 
scDRcell <- function(inpConf, inpMeta, inpdrX, inpdrY, inp1, inpsub1, inpsub2, 
                     inpsiz, inpcol, inpord, inpfsz, inpasp, inptxt, inplab){ 
    if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
    # Prepare ggData 
    ggData = inpMeta[, c(inpConf[UI == inpdrX]$ID, inpConf[UI == inpdrY]$ID, 
                         inpConf[UI == inp1]$ID, inpConf[UI == inpsub1]$ID),  
                     with = FALSE] 
    colnames(ggData) = c("X", "Y", "val", "sub") 
    rat = (max(ggData$X) - min(ggData$X)) / (max(ggData$Y) - min(ggData$Y)) 
    bgCells = FALSE 
    if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
        bgCells = TRUE 
        ggData2 = ggData[!sub %in% inpsub2] 
        ggData = ggData[sub %in% inpsub2] 
    } 
    if(inpord == "Max-1st"){ 
        ggData = ggData[order(val)] 
    } else if(inpord == "Min-1st"){ 
        ggData = ggData[order(-val)] 
    } else if(inpord == "Random"){ 
        ggData = ggData[sample(nrow(ggData))] 
    } 
    
    # Do factoring if required 
    if(!is.na(inpConf[UI == inp1]$fCL)){ 
        ggCol = strsplit(inpConf[UI == inp1]$fCL, "\\|")[[1]] 
        names(ggCol) = levels(ggData$val) 
        ggLvl = levels(ggData$val)[levels(ggData$val) %in% unique(ggData$val)] 
        ggData$val = factor(ggData$val, levels = ggLvl) 
        ggCol = ggCol[ggLvl] 
    } 
    
    # Actual ggplot 
    ggOut = ggplot(ggData, aes(X, Y, color = val)) 
    if(bgCells){ 
        ggOut = ggOut + 
            geom_point(data = ggData2, color = "snow2", size = inpsiz, shape = 16) 
    } 
    ggOut = ggOut + 
        geom_point(size = inpsiz, shape = 16) + xlab(inpdrX) + ylab(inpdrY) + 
        sctheme(base_size = sList[inpfsz], XYval = inptxt) 
    if(is.na(inpConf[UI == inp1]$fCL)){ 
        ggOut = ggOut + scale_color_gradientn("", colours = cList[[inpcol]]) + 
            guides(color = guide_colorbar(barwidth = 15)) 
    } else { 
        sListX = min(nchar(paste0(levels(ggData$val), collapse = "")), 200) 
        sListX = 0.75 * (sList - (1.5 * floor(sListX/50))) 
        ggOut = ggOut + scale_color_manual("", values = ggCol) + 
            guides(color = guide_legend(override.aes = list(size = 5),  
                                        nrow = inpConf[UI == inp1]$fRow)) + 
            theme(legend.text = element_text(size = sListX[inpfsz])) 
        if(inplab){ 
            ggData3 = ggData[, .(X = mean(X), Y = mean(Y)), by = "val"] 
            lListX = min(nchar(paste0(ggData3$val, collapse = "")), 200) 
            lListX = lList - (0.25 * floor(lListX/50)) 
            ggOut = ggOut + 
                geom_text_repel(data = ggData3, aes(X, Y, label = val), 
                                color = "grey10", bg.color = "grey95", bg.r = 0.15, 
                                size = lListX[inpfsz], seed = 42) 
        } 
    } 
    if(inpasp == "Square") { 
        ggOut = ggOut + coord_fixed(ratio = rat) 
    } else if(inpasp == "Fixed") { 
        ggOut = ggOut + coord_fixed() 
    } 
    return(ggOut) 
} 

scDRnum <- function(inpConf, inpMeta, inp1, inp2, inpsub1, inpsub2, 
                    inpH5, inpGene, inpsplt){ 
    if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
    # Prepare ggData 
    ggData = inpMeta[, c(inpConf[UI == inp1]$ID, inpConf[UI == inpsub1]$ID), 
                     with = FALSE] 
    colnames(ggData) = c("group", "sub") 
    h5file <- H5File$new(inpH5, mode = "r") 
    h5data <- h5file[["grp"]][["data"]] 
    ggData$val2 = h5data$read(args = list(inpGene[inp2], quote(expr=))) 
    ggData[val2 < 0]$val2 = 0 
    h5file$close_all() 
    if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
        ggData = ggData[sub %in% inpsub2] 
    } 
    
    # Split inp1 if necessary 
    if(is.na(inpConf[UI == inp1]$fCL)){ 
        if(inpsplt == "Quartile"){nBk = 4} 
        if(inpsplt == "Decile"){nBk = 10} 
        ggData$group = cut(ggData$group, breaks = nBk) 
    } 
    
    # Actual data.table 
    ggData$express = FALSE 
    ggData[val2 > 0]$express = TRUE 
    ggData1 = ggData[express == TRUE, .(nExpress = .N), by = "group"] 
    ggData = ggData[, .(nCells = .N), by = "group"] 
    ggData = ggData1[ggData, on = "group"] 
    ggData = ggData[, c("group", "nCells", "nExpress"), with = FALSE] 
    ggData[is.na(nExpress)]$nExpress = 0 
    ggData$pctExpress = 100 * ggData$nExpress / ggData$nCells 
    ggData = ggData[order(group)] 
    colnames(ggData)[3] = paste0(colnames(ggData)[3], "_", inp2) 
    return(ggData) 
} 
# Plot gene expression on dimred 
scDRgene <- function(inpConf, inpMeta, inpdrX, inpdrY, inp1, inpsub1, inpsub2, 
                     inpH5, inpGene, 
                     inpsiz, inpcol, inpord, inpfsz, inpasp, inptxt){ 
    if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
    # Prepare ggData 
    ggData = inpMeta[, c(inpConf[UI == inpdrX]$ID, inpConf[UI == inpdrY]$ID, 
                         inpConf[UI == inpsub1]$ID),  
                     with = FALSE] 
    colnames(ggData) = c("X", "Y", "sub") 
    rat = (max(ggData$X) - min(ggData$X)) / (max(ggData$Y) - min(ggData$Y)) 
    
    h5file <- H5File$new(inpH5, mode = "r") 
    h5data <- h5file[["grp"]][["data"]] 
    ggData$val = h5data$read(args = list(inpGene[inp1], quote(expr=))) 
    ggData[val < 0]$val = 0 
    h5file$close_all() 
    bgCells = FALSE 
    if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
        bgCells = TRUE 
        ggData2 = ggData[!sub %in% inpsub2] 
        ggData = ggData[sub %in% inpsub2] 
    } 
    if(inpord == "Max-1st"){ 
        ggData = ggData[order(val)] 
    } else if(inpord == "Min-1st"){ 
        ggData = ggData[order(-val)] 
    } else if(inpord == "Random"){ 
        ggData = ggData[sample(nrow(ggData))] 
    } 
    
    # Actual ggplot 
    ggOut = ggplot(ggData, aes(X, Y, color = val)) 
    if(bgCells){ 
        ggOut = ggOut + 
            geom_point(data = ggData2, color = "snow2", size = inpsiz, shape = 16) 
    } 
    ggOut = ggOut + 
        geom_point(size = inpsiz, shape = 16) + xlab(inpdrX) + ylab(inpdrY) + 
        sctheme(base_size = sList[inpfsz], XYval = inptxt) +  
        scale_color_gradientn(inp1, colours = cList[[inpcol]]) + 
        guides(color = guide_colorbar(barwidth = 15)) 
    if(inpasp == "Square") { 
        ggOut = ggOut + coord_fixed(ratio = rat) 
    } else if(inpasp == "Fixed") { 
        ggOut = ggOut + coord_fixed() 
    } 
    return(ggOut) 
} 

# Plot gene coexpression on dimred 
bilinear <- function(x,y,xy,Q11,Q21,Q12,Q22){ 
    oup = (xy-x)*(xy-y)*Q11 + x*(xy-y)*Q21 + (xy-x)*y*Q12 + x*y*Q22 
    oup = oup / (xy*xy) 
    return(oup) 
} 
scDRcoex <- function(inpConf, inpMeta, inpdrX, inpdrY, inp1, inp2, 
                     inpsub1, inpsub2, inpH5, inpGene, 
                     inpsiz, inpcol, inpord, inpfsz, inpasp, inptxt){ 
    if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
    # Prepare ggData 
    ggData = inpMeta[, c(inpConf[UI == inpdrX]$ID, inpConf[UI == inpdrY]$ID, 
                         inpConf[UI == inpsub1]$ID),  
                     with = FALSE] 
    colnames(ggData) = c("X", "Y", "sub") 
    rat = (max(ggData$X) - min(ggData$X)) / (max(ggData$Y) - min(ggData$Y)) 
    
    h5file <- H5File$new(inpH5, mode = "r") 
    h5data <- h5file[["grp"]][["data"]] 
    ggData$val1 = h5data$read(args = list(inpGene[inp1], quote(expr=))) 
    ggData[val1 < 0]$val1 = 0 
    ggData$val2 = h5data$read(args = list(inpGene[inp2], quote(expr=))) 
    ggData[val2 < 0]$val2 = 0 
    h5file$close_all() 
    bgCells = FALSE 
    if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
        bgCells = TRUE 
        ggData2 = ggData[!sub %in% inpsub2] 
        ggData = ggData[sub %in% inpsub2] 
    } 
    
    # Generate coex color palette 
    cInp = strsplit(inpcol, "; ")[[1]] 
    if(cInp[1] == "Red (Gene1)"){ 
        c10 = c(255,0,0) 
    } else if(cInp[1] == "Orange (Gene1)"){ 
        c10 = c(255,140,0) 
    } else { 
        c10 = c(0,255,0) 
    } 
    if(cInp[2] == "Green (Gene2)"){ 
        c01 = c(0,255,0) 
    } else { 
        c01 = c(0,0,255) 
    } 
    c00 = c(217,217,217) ; c11 = c10 + c01 
    nGrid = 16; nPad = 2; nTot = nGrid + nPad * 2 
    gg = data.table(v1 = rep(0:nTot,nTot+1), v2 = sort(rep(0:nTot,nTot+1))) 
    gg$vv1 = gg$v1 - nPad ; gg[vv1 < 0]$vv1 = 0; gg[vv1 > nGrid]$vv1 = nGrid 
    gg$vv2 = gg$v2 - nPad ; gg[vv2 < 0]$vv2 = 0; gg[vv2 > nGrid]$vv2 = nGrid 
    gg$cR = bilinear(gg$vv1, gg$vv2, nGrid, c00[1], c10[1], c01[1], c11[1]) 
    gg$cG = bilinear(gg$vv1, gg$vv2, nGrid, c00[2], c10[2], c01[2], c11[2]) 
    gg$cB = bilinear(gg$vv1, gg$vv2, nGrid, c00[3], c10[3], c01[3], c11[3]) 
    gg$cMix = rgb(gg$cR, gg$cG, gg$cB, maxColorValue = 255) 
    gg = gg[, c("v1", "v2", "cMix")] 
    
    # Map colours 
    ggData$v1 = round(nTot * ggData$val1 / max(ggData$val1)) 
    ggData$v2 = round(nTot * ggData$val2 / max(ggData$val2)) 
    ggData$v0 = ggData$v1 + ggData$v2 
    ggData = gg[ggData, on = c("v1", "v2")] 
    if(inpord == "Max-1st"){ 
        ggData = ggData[order(v0)] 
    } else if(inpord == "Min-1st"){ 
        ggData = ggData[order(-v0)] 
    } else if(inpord == "Random"){ 
        ggData = ggData[sample(nrow(ggData))] 
    } 
    
    # Actual ggplot 
    ggOut = ggplot(ggData, aes(X, Y)) 
    if(bgCells){ 
        ggOut = ggOut + 
            geom_point(data = ggData2, color = "snow2", size = inpsiz, shape = 16) 
    } 
    ggOut = ggOut + 
        geom_point(size = inpsiz, shape = 16, color = ggData$cMix) + 
        xlab(inpdrX) + ylab(inpdrY) + 
        sctheme(base_size = sList[inpfsz], XYval = inptxt) + 
        scale_color_gradientn(inp1, colours = cList[[1]]) + 
        guides(color = guide_colorbar(barwidth = 15)) 
    if(inpasp == "Square") { 
        ggOut = ggOut + coord_fixed(ratio = rat) 
    } else if(inpasp == "Fixed") { 
        ggOut = ggOut + coord_fixed() 
    } 
    return(ggOut) 
} 

scDRcoexLeg <- function(inp1, inp2, inpcol, inpfsz){ 
    # Generate coex color palette 
    cInp = strsplit(inpcol, "; ")[[1]] 
    if(cInp[1] == "Red (Gene1)"){ 
        c10 = c(255,0,0) 
    } else if(cInp[1] == "Orange (Gene1)"){ 
        c10 = c(255,140,0) 
    } else { 
        c10 = c(0,255,0) 
    } 
    if(cInp[2] == "Green (Gene2)"){ 
        c01 = c(0,255,0) 
    } else { 
        c01 = c(0,0,255) 
    } 
    c00 = c(217,217,217) ; c11 = c10 + c01 
    nGrid = 16; nPad = 2; nTot = nGrid + nPad * 2 
    gg = data.table(v1 = rep(0:nTot,nTot+1), v2 = sort(rep(0:nTot,nTot+1))) 
    gg$vv1 = gg$v1 - nPad ; gg[vv1 < 0]$vv1 = 0; gg[vv1 > nGrid]$vv1 = nGrid 
    gg$vv2 = gg$v2 - nPad ; gg[vv2 < 0]$vv2 = 0; gg[vv2 > nGrid]$vv2 = nGrid 
    gg$cR = bilinear(gg$vv1, gg$vv2, nGrid, c00[1], c10[1], c01[1], c11[1]) 
    gg$cG = bilinear(gg$vv1, gg$vv2, nGrid, c00[2], c10[2], c01[2], c11[2]) 
    gg$cB = bilinear(gg$vv1, gg$vv2, nGrid, c00[3], c10[3], c01[3], c11[3]) 
    gg$cMix = rgb(gg$cR, gg$cG, gg$cB, maxColorValue = 255) 
    gg = gg[, c("v1", "v2", "cMix")] 
    
    # Actual ggplot 
    ggOut = ggplot(gg, aes(v1, v2)) + 
        geom_tile(fill = gg$cMix) + 
        xlab(inp1) + ylab(inp2) + coord_fixed(ratio = 1) + 
        scale_x_continuous(breaks = c(0, nTot), label = c("low", "high")) + 
        scale_y_continuous(breaks = c(0, nTot), label = c("low", "high")) + 
        sctheme(base_size = sList[inpfsz], XYval = TRUE) 
    return(ggOut) 
} 

scDRcoexNum <- function(inpConf, inpMeta, inp1, inp2, 
                        inpsub1, inpsub2, inpH5, inpGene){ 
    if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
    # Prepare ggData 
    ggData = inpMeta[, c(inpConf[UI == inpsub1]$ID), with = FALSE] 
    colnames(ggData) = c("sub") 
    h5file <- H5File$new(inpH5, mode = "r") 
    h5data <- h5file[["grp"]][["data"]] 
    ggData$val1 = h5data$read(args = list(inpGene[inp1], quote(expr=))) 
    ggData[val1 < 0]$val1 = 0 
    ggData$val2 = h5data$read(args = list(inpGene[inp2], quote(expr=))) 
    ggData[val2 < 0]$val2 = 0 
    h5file$close_all() 
    if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
        ggData = ggData[sub %in% inpsub2] 
    } 
    
    # Actual data.table 
    ggData$express = "none" 
    ggData[val1 > 0]$express = inp1 
    ggData[val2 > 0]$express = inp2 
    ggData[val1 > 0 & val2 > 0]$express = "both" 
    ggData$express = factor(ggData$express, levels = unique(c("both", inp1, inp2, "none"))) 
    ggData = ggData[, .(nCells = .N), by = "express"] 
    ggData$percent = 100 * ggData$nCells / sum(ggData$nCells) 
    ggData = ggData[order(express)] 
    colnames(ggData)[1] = "expression > 0" 
    return(ggData) 
} 

# Plot violin / boxplot 
scVioBox <- function(inpConf, inpMeta, inp1, inp2, 
                     inpsub1, inpsub2, inpH5, inpGene, 
                     inptyp, inppts, inpsiz, inpfsz){ 
    if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
    # Prepare ggData 
    ggData = inpMeta[, c(inpConf[UI == inp1]$ID, inpConf[UI == inpsub1]$ID), 
                     with = FALSE] 
    colnames(ggData) = c("X", "sub") 
    
    # Load in either cell meta or gene expr
    if(inp2 %in% inpConf$UI){ 
        ggData$val = inpMeta[[inpConf[UI == inp2]$ID]] 
    } else { 
        h5file <- H5File$new(inpH5, mode = "r") 
        h5data <- h5file[["grp"]][["data"]] 
        ggData$val = h5data$read(args = list(inpGene[inp2], quote(expr=))) 
        ggData[val < 0]$val = 0 
        set.seed(42) 
        tmpNoise = rnorm(length(ggData$val)) * diff(range(ggData$val)) / 1000 
        ggData$val = ggData$val + tmpNoise 
        h5file$close_all() 
    } 
    if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
        ggData = ggData[sub %in% inpsub2] 
    } 
    
    # Do factoring 
    ggCol = strsplit(inpConf[UI == inp1]$fCL, "\\|")[[1]] 
    names(ggCol) = levels(ggData$X) 
    ggLvl = levels(ggData$X)[levels(ggData$X) %in% unique(ggData$X)] 
    ggData$X = factor(ggData$X, levels = ggLvl) 
    ggCol = ggCol[ggLvl] 
    
    # Actual ggplot 
    if(inptyp == "violin"){ 
        ggOut = ggplot(ggData, aes(X, val, fill = X)) + geom_violin(scale = "width") 
    } else { 
        ggOut = ggplot(ggData, aes(X, val, fill = X)) + geom_boxplot() 
    } 
    if(inppts){ 
        ggOut = ggOut + geom_jitter(size = inpsiz, shape = 16) 
    } 
    ggOut = ggOut + xlab(inp1) + ylab(inp2) + 
        sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +  
        scale_fill_manual("", values = ggCol) +
        theme(legend.position = "none")
    return(ggOut) 
} 

# Plot proportion plot 
scProp <- function(inpConf, inpMeta, inp1, inp2, inpsub1, inpsub2, 
                   inptyp, inpflp, inpfsz){ 
    if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
    # Prepare ggData 
    ggData = inpMeta[, c(inpConf[UI == inp1]$ID, inpConf[UI == inp2]$ID, 
                         inpConf[UI == inpsub1]$ID),  
                     with = FALSE] 
    colnames(ggData) = c("X", "grp", "sub") 
    if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
        ggData = ggData[sub %in% inpsub2] 
    } 
    ggData = ggData[, .(nCells = .N), by = c("X", "grp")] 
    ggData = ggData[, {tot = sum(nCells) 
    .SD[,.(pctCells = 100 * sum(nCells) / tot, 
           nCells = nCells), by = "grp"]}, by = "X"] 
    
    # Do factoring 
    ggCol = strsplit(inpConf[UI == inp2]$fCL, "\\|")[[1]] 
    names(ggCol) = levels(ggData$grp) 
    ggLvl = levels(ggData$grp)[levels(ggData$grp) %in% unique(ggData$grp)] 
    ggData$grp = factor(ggData$grp, levels = ggLvl) 
    ggCol = ggCol[ggLvl] 
    
    # Actual ggplot 
    if(inptyp == "Proportion"){ 
        ggOut = ggplot(ggData, aes(X, pctCells, fill = grp)) + 
            geom_col() + ylab("Cell Proportion (%)") 
    } else { 
        ggOut = ggplot(ggData, aes(X, nCells, fill = grp)) + 
            geom_col() + ylab("Number of Cells") 
    } 
    if(inpflp){ 
        ggOut = ggOut + coord_flip() 
    } 
    ggOut = ggOut + xlab(inp1) + 
        sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +  
        scale_fill_manual("", values = ggCol) + 
        theme(legend.position = "right") 
    return(ggOut) 
} 

# Get gene list 
scGeneList <- function(inp, inpGene){ 
    geneList = data.table(gene = unique(trimws(strsplit(inp, ",|;|
")[[1]])), 
                          present = TRUE) 
    geneList[!gene %in% names(inpGene)]$present = FALSE 
    return(geneList) 
} 

# Plot gene expression bubbleplot / heatmap 
scBubbHeat <- function(inpConf, inpMeta, inp, inpGrp, inpPlt, 
                       inpsub1, inpsub2, inpH5, inpGene, inpScl, inpRow, inpCol, 
                       inpcols, inpfsz, save = FALSE){ 
    if(is.null(inpsub1)){inpsub1 = inpConf$UI[1]} 
    # Identify genes that are in our dataset 
    geneList = scGeneList(inp, inpGene) 
    geneList = geneList[present == TRUE] 
    shiny::validate(need(nrow(geneList) <= 50, "More than 50 genes to plot! Please reduce the gene list!")) 
    shiny::validate(need(nrow(geneList) > 1, "Please input at least 2 genes to plot!")) 
    
    # Prepare ggData 
    h5file <- H5File$new(inpH5, mode = "r") 
    h5data <- h5file[["grp"]][["data"]] 
    ggData = data.table() 
    for(iGene in geneList$gene){ 
        tmp = inpMeta[, c("sampleID", inpConf[UI == inpsub1]$ID), with = FALSE] 
        colnames(tmp) = c("sampleID", "sub") 
        tmp$grpBy = inpMeta[[inpConf[UI == inpGrp]$ID]] 
        tmp$geneName = iGene 
        tmp$val = h5data$read(args = list(inpGene[iGene], quote(expr=))) 
        ggData = rbindlist(list(ggData, tmp)) 
    } 
    h5file$close_all() 
    if(length(inpsub2) != 0 & length(inpsub2) != nlevels(ggData$sub)){ 
        ggData = ggData[sub %in% inpsub2] 
    } 
    shiny::validate(need(uniqueN(ggData$grpBy) > 1, "Only 1 group present, unable to plot!")) 
    
    # Aggregate 
    ggData$val = expm1(ggData$val) 
    ggData = ggData[, .(val = mean(val), prop = sum(val>0) / length(sampleID)), 
                    by = c("geneName", "grpBy")] 
    ggData$val = log1p(ggData$val) 
    
    # Scale if required 
    colRange = range(ggData$val) 
    if(inpScl){ 
        ggData[, val:= scale(val), keyby = "geneName"] 
        colRange = c(-max(abs(range(ggData$val))), max(abs(range(ggData$val)))) 
    } 
    
    # hclust row/col if necessary 
    ggMat = dcast.data.table(ggData, geneName~grpBy, value.var = "val") 
    tmp = ggMat$geneName 
    ggMat = as.matrix(ggMat[, -1]) 
    rownames(ggMat) = tmp 
    if(inpRow){ 
        hcRow = dendro_data(as.dendrogram(hclust(dist(ggMat)))) 
        ggRow = ggplot() + coord_flip() + 
            geom_segment(data = hcRow$segments, aes(x=x,y=y,xend=xend,yend=yend)) + 
            scale_y_continuous(breaks = rep(0, uniqueN(ggData$grpBy)), 
                               labels = unique(ggData$grpBy), expand = c(0, 0)) + 
            scale_x_continuous(breaks = seq_along(hcRow$labels$label), 
                               labels = hcRow$labels$label, expand = c(0, 0.5)) + 
            sctheme(base_size = sList[inpfsz]) + 
            theme(axis.title = element_blank(), axis.line = element_blank(), 
                  axis.ticks = element_blank(), axis.text.y = element_blank(), 
                  axis.text.x = element_text(color="white", angle = 45, hjust = 1)) 
        ggData$geneName = factor(ggData$geneName, levels = hcRow$labels$label) 
    } else { 
        ggData$geneName = factor(ggData$geneName, levels = rev(geneList$gene)) 
    } 
    if(inpCol){ 
        hcCol = dendro_data(as.dendrogram(hclust(dist(t(ggMat))))) 
        ggCol = ggplot() + 
            geom_segment(data = hcCol$segments, aes(x=x,y=y,xend=xend,yend=yend)) + 
            scale_x_continuous(breaks = seq_along(hcCol$labels$label), 
                               labels = hcCol$labels$label, expand = c(0.05, 0)) + 
            scale_y_continuous(breaks = rep(0, uniqueN(ggData$geneName)), 
                               labels = unique(ggData$geneName), expand=c(0,0)) + 
            sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) + 
            theme(axis.title = element_blank(), axis.line = element_blank(), 
                  axis.ticks = element_blank(), axis.text.x = element_blank(), 
                  axis.text.y = element_text(color = "white")) 
        ggData$grpBy = factor(ggData$grpBy, levels = hcCol$labels$label) 
    } 
    
    # Actual plot according to plottype 
    if(inpPlt == "Bubbleplot"){ 
        # Bubbleplot 
        ggOut = ggplot(ggData, aes(grpBy, geneName, color = val, size = prop)) + 
            geom_point() +  
            sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) +  
            scale_x_discrete(expand = c(0.05, 0)) +  
            scale_y_discrete(expand = c(0, 0.5)) + 
            scale_size_continuous("proportion", range = c(0, 8), 
                                  limits = c(0, 1), breaks = c(0.00,0.25,0.50,0.75,1.00)) + 
            scale_color_gradientn("expression", limits = colRange, colours = cList[[inpcols]]) + 
            guides(color = guide_colorbar(barwidth = 15)) + 
            theme(axis.title = element_blank(), legend.box = "vertical") 
    } else { 
        # Heatmap 
        ggOut = ggplot(ggData, aes(grpBy, geneName, fill = val)) + 
            geom_tile() +  
            sctheme(base_size = sList[inpfsz], Xang = 45, XjusH = 1) + 
            scale_x_discrete(expand = c(0.05, 0)) +  
            scale_y_discrete(expand = c(0, 0.5)) + 
            scale_fill_gradientn("expression", limits = colRange, colours = cList[[inpcols]]) + 
            guides(fill = guide_colorbar(barwidth = 15)) + 
            theme(axis.title = element_blank()) 
    } 
    
    # Final tidy 
    ggLeg = g_legend(ggOut) 
    ggOut = ggOut + theme(legend.position = "none") 
    if(!save){ 
        if(inpRow & inpCol){ggOut =  
            grid.arrange(ggOut, ggLeg, ggCol, ggRow, widths = c(7,1), heights = c(1,7,2),  
                         layout_matrix = rbind(c(3,NA),c(1,4),c(2,NA)))  
        } else if(inpRow){ggOut =  
            grid.arrange(ggOut, ggLeg, ggRow, widths = c(7,1), heights = c(7,2),  
                         layout_matrix = rbind(c(1,3),c(2,NA)))  
        } else if(inpCol){ggOut =  
            grid.arrange(ggOut, ggLeg, ggCol, heights = c(1,7,2),  
                         layout_matrix = rbind(c(3),c(1),c(2)))  
        } else {ggOut =  
            grid.arrange(ggOut, ggLeg, heights = c(7,2),  
                         layout_matrix = rbind(c(1),c(2)))  
        }  
    } else { 
        if(inpRow & inpCol){ggOut =  
            arrangeGrob(ggOut, ggLeg, ggCol, ggRow, widths = c(7,1), heights = c(1,7,2),  
                        layout_matrix = rbind(c(3,NA),c(1,4),c(2,NA)))  
        } else if(inpRow){ggOut =  
            arrangeGrob(ggOut, ggLeg, ggRow, widths = c(7,1), heights = c(7,2),  
                        layout_matrix = rbind(c(1,3),c(2,NA)))  
        } else if(inpCol){ggOut =  
            arrangeGrob(ggOut, ggLeg, ggCol, heights = c(1,7,2),  
                        layout_matrix = rbind(c(3),c(1),c(2)))  
        } else {ggOut =  
            arrangeGrob(ggOut, ggLeg, heights = c(7,2),  
                        layout_matrix = rbind(c(1),c(2)))  
        }  
    } 
    return(ggOut) 
} 
#####


# Load or define your neural network model here
model_NN <- keras::load_model_hdf5("my_model_MM_prueba_2_all.h5")

# Define the UI
ui <- navbarPage(theme = light, checkboxInput("dark_mode", "Dark mode"), collapsable = FALSE, title = "Multiple Myeloma App",position = "static-top",
                tabsetPanel(
                            tabPanel(HTML("Single-cell RNA-seq Workflow"),
                                     mainPanel( uiOutput("pdfview"))
                            ),
                            tabPanel(HTML("Single-cell RNA-seq Analyses in Multiple Myeloma"),
                                     style = "margin-bottom: 4.5rem;",
                                        # Custom CSS to increase the sidebar height
                                        tags$style(HTML(".sidebar { height: 200vh; overflow-y: auto; }")),
                                        tabsetPanel(
                                                  # Tabpanel 1  ##### 
                                                    tabPanel(
                                                        HTML("CellInfo vs GeneExpr"),
                                                        h4("Cell information vs gene expression on reduced dimensions"),
                                                        "In this tab, users can visualise both cell information and gene ",
                                                        "expression side-by-side on low-dimensional representions.",
                                                        br(),br(),
                                                        fluidRow(
                                                            column(
                                                                3, h4("Dimension Reduction"),
                                                                fluidRow(
                                                                    column(
                                                                        12, selectInput("sc1a1drX", "X-axis:", choices = sc1conf[dimred == TRUE]$UI,
                                                                                        selected = sc1def$dimred[1]),
                                                                        selectInput("sc1a1drY", "Y-axis:", choices = sc1conf[dimred == TRUE]$UI,
                                                                                    selected = sc1def$dimred[2]))
                                                                )
                                                            ), # End of column (6 space)
                                                            column(
                                                                3, actionButton("sc1a1togL", "Toggle to subset cells"),
                                                                conditionalPanel(
                                                                    condition = "input.sc1a1togL % 2 == 1",
                                                                    selectInput("sc1a1sub1", "Cell information to subset:",
                                                                                choices = sc1conf[grp == TRUE]$UI[-2],
                                                                                selected = sc1def$grp1),
                                                                    uiOutput("sc1a1sub1.ui"),
                                                                    actionButton("sc1a1sub1all", "Select all groups", class = "btn btn-primary"),
                                                                    actionButton("sc1a1sub1non", "Deselect all groups", class = "btn btn-primary")
                                                                )
                                                            ), # End of column (6 space)
                                                            column(
                                                                6, actionButton("sc1a1tog0", "Toggle graphics controls"),
                                                                conditionalPanel(
                                                                    condition = "input.sc1a1tog0 % 2 == 1",
                                                                    fluidRow(
                                                                        column(
                                                                            6, sliderInput("sc1a1siz", "Point size:",
                                                                                           min = 0, max = 4, value = 1.25, step = 0.25),
                                                                            radioButtons("sc1a1psz", "Plot size:",
                                                                                         choices = c("Small", "Medium", "Large"),
                                                                                         selected = "Medium", inline = TRUE),
                                                                            radioButtons("sc1a1fsz", "Font size:",
                                                                                         choices = c("Small", "Medium", "Large"),
                                                                                         selected = "Medium", inline = TRUE)
                                                                        ),
                                                                        column(
                                                                            6, radioButtons("sc1a1asp", "Aspect ratio:",
                                                                                            choices = c("Square", "Fixed", "Free"),
                                                                                            selected = "Square", inline = TRUE),
                                                                            checkboxInput("sc1a1txt", "Show axis text", value = FALSE)
                                                                        )
                                                                    )
                                                                )
                                                            )  # End of column (6 space)
                                                        ),   # End of fluidRow (4 space)
                                                        fluidRow(
                                                            column(
                                                                6, style="border-right: 2px solid black", h4("Cell information"),
                                                                fluidRow(
                                                                    column(
                                                                        6, selectInput("sc1a1inp1", "Cell information:",
                                                                                       choices = sc1conf$UI[c(1,5,6)],
                                                                                       selected = sc1def$meta1) %>%
                                                                            helper(type = "inline", size = "m", fade = TRUE,
                                                                                   title = "Cell information to colour cells by",
                                                                                   content = c("Select cell information to colour cells",
                                                                                               "- Categorical covariates have a fixed colour palette",
                                                                                               paste0("- Continuous covariates are coloured in a ",
                                                                                                      "Blue-Yellow-Red colour scheme, which can be ",
                                                                                                      "changed in the plot controls")))
                                                                    ),
                                                                    column(
                                                                        6, actionButton("sc1a1tog1", "Toggle plot controls"),
                                                                        conditionalPanel(
                                                                            condition = "input.sc1a1tog1 % 2 == 1",
                                                                            radioButtons("sc1a1col1", "Colour (Continuous data):",
                                                                                         choices = c("White-Red","Blue-Yellow-Red","Yellow-Green-Purple"),
                                                                                         selected = "Blue-Yellow-Red"),
                                                                            radioButtons("sc1a1ord1", "Plot order:",
                                                                                         choices = c("Max-1st", "Min-1st", "Original", "Random"),
                                                                                         selected = "Original", inline = TRUE),
                                                                            checkboxInput("sc1a1lab1", "Show cell info labels", value = TRUE)
                                                                        )
                                                                    )
                                                                ),
                                                                fluidRow(column(12, uiOutput("sc1a1oup1.ui"))),
                                                                downloadButton("sc1a1oup1.pdf", "Download PDF"),
                                                                downloadButton("sc1a1oup1.png", "Download PNG"), br(),
                                                                div(style="display:inline-block",
                                                                    numericInput("sc1a1oup1.h", "PDF / PNG height:", width = "138px",
                                                                                 min = 4, max = 20, value = 6, step = 0.5)),
                                                                div(style="display:inline-block",
                                                                    numericInput("sc1a1oup1.w", "PDF / PNG width:", width = "138px",
                                                                                 min = 4, max = 20, value = 8, step = 0.5)), br(),
                                                                actionButton("sc1a1tog9", "Toggle to show cell numbers / statistics"),
                                                                conditionalPanel(
                                                                    condition = "input.sc1a1tog9 % 2 == 1",
                                                                    h4("Cell numbers / statistics"),
                                                                    radioButtons("sc1a1splt", "Split continuous cell info into:",
                                                                                 choices = c("Quartile", "Decile"),
                                                                                 selected = "Decile", inline = TRUE),
                                                                    dataTableOutput("sc1a1.dt")
                                                                )
                                                            ), # End of column (6 space)
                                                            column(6, 
                                                                   h4("Gene expression"),
                                                                    fluidRow(
                                                                        column(6, 
                                                                               selectInput("sc1a1inp2", "Gene name:", choices=NULL) %>%
                                                                                helper(type = "inline", size = "m", fade = TRUE,
                                                                                       title = "Gene expression to colour cells by",
                                                                                       content = c("Select gene to colour cells by gene expression",
                                                                                                   paste0("- Gene expression are coloured in a ",
                                                                                                          "White-Red colour scheme which can be ",
                                                                                                          "changed in the plot controls")))
                                                                        ),
                                                                        column(6, 
                                                                               actionButton("sc1a1tog2", "Toggle plot controls"),
                                                                                conditionalPanel(
                                                                                    condition = "input.sc1a1tog2 % 2 == 1",
                                                                                    radioButtons("sc1a1col2", "Colour:",
                                                                                                 choices = c("White-Red","Blue-Yellow-Red","Yellow-Green-Purple"),
                                                                                                 selected = "White-Red"),
                                                                                    radioButtons("sc1a1ord2", "Plot order:",
                                                                                                 choices = c("Max-1st", "Min-1st", "Original", "Random"),
                                                                                                 selected = "Max-1st", inline = TRUE)
                                                                                )
                                                                        )
                                                                ),
                                                                fluidRow(column(12, uiOutput("sc1a1oup2.ui"))),
                                                                downloadButton("sc1a1oup2.pdf", "Download PDF"),
                                                                downloadButton("sc1a1oup2.png", "Download PNG"), br(),
                                                                div(style="display:inline-block",
                                                                    numericInput("sc1a1oup2.h", "PDF / PNG height:", width = "138px",
                                                                                 min = 4, max = 20, value = 6, step = 0.5)),
                                                                div(style="display:inline-block",
                                                                    numericInput("sc1a1oup2.w", "PDF / PNG width:", width = "138px",
                                                                                 min = 4, max = 20, value = 8, step = 0.5))
                                                            )  # End of column (6 space)
                                                        )    # End of fluidRow (4 space)
                                                    ),     # End of tab (2 space) TABPANEL1
                                                  # Tabpanel 2  ##### 
                                                    ### Tab1.a2: cellInfo vs cellInfo on dimRed
                                                    tabPanel(
                                                        HTML("CellInfo vs CellInfo"),
                                                        h4("Cell information vs cell information on dimension reduction"),
                                                        "In this tab, users can visualise two cell informations side-by-side ",
                                                        "on low-dimensional representions.",
                                                        br(),
                                                        br(),
                                                        fluidRow(
                                                            column(3, 
                                                                   h4("Dimension Reduction"),
                                                                   selectInput("sc1a2drX", "X-axis:", choices = sc1conf[dimred == TRUE]$UI, selected = sc1def$dimred[1]),
                                                                   selectInput("sc1a2drY", "Y-axis:", choices = sc1conf[dimred == TRUE]$UI, selected = sc1def$dimred[2])
                                                            ), # End of column (6 space)
                                                            column(3, 
                                                                   actionButton("sc1a2togL", "Toggle to subset cells"),
                                                                    conditionalPanel(condition = "input.sc1a2togL % 2 == 1",
                                                                        selectInput("sc1a2sub1", "Cell information to subset:",
                                                                                choices = sc1conf[grp == TRUE]$UI[-2],
                                                                                selected = sc1def$grp1),
                                                                        uiOutput("sc1a2sub1.ui"),
                                                                        actionButton("sc1a2sub1all", "Select all groups", class = "btn btn-primary"),
                                                                        actionButton("sc1a2sub1non", "Deselect all groups", class = "btn btn-primary")
                                                                )
                                                            ), # End of column (6 space)
                                                            column(6, 
                                                                   actionButton("sc1a2tog0", "Toggle graphics controls"),
                                                                    conditionalPanel(
                                                                        condition = "input.sc1a2tog0 % 2 == 1",
                                                                        fluidRow(
                                                                            column(6, 
                                                                                    sliderInput("sc1a2siz", "Point size:",
                                                                                               min = 0, max = 4, value = 1.25, step = 0.25),
                                                                                    radioButtons("sc1a2psz", "Plot size:",
                                                                                             choices = c("Small", "Medium", "Large"),
                                                                                             selected = "Medium", inline = TRUE),
                                                                                    radioButtons("sc1a2fsz", "Font size:",
                                                                                             choices = c("Small", "Medium", "Large"),
                                                                                             selected = "Medium", inline = TRUE)
                                                                            ),
                                                                            column(6, 
                                                                                   radioButtons("sc1a2asp", "Aspect ratio:",
                                                                                                choices = c("Square", "Fixed", "Free"),
                                                                                                selected = "Square", inline = TRUE),
                                                                                    checkboxInput("sc1a2txt", "Show axis text", value = FALSE)
                                                                            )
                                                                        )
                                                                )
                                                            )  # End of column (6 space)
                                                        ),   # End of fluidRow (4 space)
                                                        fluidRow(
                                                            column(
                                                                6, style="border-right: 2px solid black", h4("Cell information 1"),
                                                                fluidRow(
                                                                    column(
                                                                        6, selectInput("sc1a2inp1", "Cell information:",
                                                                                       choices = sc1conf$UI[c(1,5,6)],
                                                                                       selected = sc1def$meta1) %>%
                                                                            helper(type = "inline", size = "m", fade = TRUE,
                                                                                   title = "Cell information to colour cells by",
                                                                                   content = c("Select cell information to colour cells",
                                                                                               "- Categorical covariates have a fixed colour palette",
                                                                                               paste0("- Continuous covariates are coloured in a ",
                                                                                                      "Blue-Yellow-Red colour scheme, which can be ",
                                                                                                      "changed in the plot controls")))
                                                                    ),
                                                                    column(
                                                                        6, actionButton("sc1a2tog1", "Toggle plot controls"),
                                                                        conditionalPanel(
                                                                            condition = "input.sc1a2tog1 % 2 == 1",
                                                                            radioButtons("sc1a2col1", "Colour (Continuous data):",
                                                                                         choices = c("White-Red", "Blue-Yellow-Red",
                                                                                                     "Yellow-Green-Purple"),
                                                                                         selected = "Blue-Yellow-Red"),
                                                                            radioButtons("sc1a2ord1", "Plot order:",
                                                                                         choices = c("Max-1st", "Min-1st", "Original", "Random"),
                                                                                         selected = "Original", inline = TRUE),
                                                                            checkboxInput("sc1a2lab1", "Show cell info labels", value = TRUE)
                                                                        )
                                                                    )
                                                                ),
                                                                fluidRow(column(12, uiOutput("sc1a2oup1.ui"))),
                                                                downloadButton("sc1a2oup1.pdf", "Download PDF"),
                                                                downloadButton("sc1a2oup1.png", "Download PNG"), br(),
                                                                div(style="display:inline-block",
                                                                    numericInput("sc1a2oup1.h", "PDF / PNG height:", width = "138px",
                                                                                 min = 4, max = 20, value = 6, step = 0.5)),
                                                                div(style="display:inline-block",
                                                                    numericInput("sc1a2oup1.w", "PDF / PNG width:", width = "138px",
                                                                                 min = 4, max = 20, value = 8, step = 0.5))
                                                            ), # End of column (6 space)
                                                            column(
                                                                6, h4("Cell information 2"),
                                                                fluidRow(
                                                                    column(
                                                                        6, selectInput("sc1a2inp2", "Cell information:",
                                                                                       choices = sc1conf$UI[c(1,5,6)],
                                                                                       selected = sc1def$meta2) %>%
                                                                            helper(type = "inline", size = "m", fade = TRUE,
                                                                                   title = "Cell information to colour cells by",
                                                                                   content = c("Select cell information to colour cells",
                                                                                               "- Categorical covariates have a fixed colour palette",
                                                                                               paste0("- Continuous covariates are coloured in a ",
                                                                                                      "Blue-Yellow-Red colour scheme, which can be ",
                                                                                                      "changed in the plot controls")))
                                                                    ),
                                                                    column(
                                                                        6, actionButton("sc1a2tog2", "Toggle plot controls"),
                                                                        conditionalPanel(
                                                                            condition = "input.sc1a2tog2 % 2 == 1",
                                                                            radioButtons("sc1a2col2", "Colour (Continuous data):",
                                                                                         choices = c("White-Red", "Blue-Yellow-Red",
                                                                                                     "Yellow-Green-Purple"),
                                                                                         selected = "Blue-Yellow-Red"),
                                                                            radioButtons("sc1a2ord2", "Plot order:",
                                                                                         choices = c("Max-1st", "Min-1st", "Original", "Random"),
                                                                                         selected = "Original", inline = TRUE),
                                                                            checkboxInput("sc1a2lab2", "Show cell info labels", value = TRUE)
                                                                        )
                                                                    )
                                                                ),
                                                                fluidRow(column(12, uiOutput("sc1a2oup2.ui"))),
                                                                downloadButton("sc1a2oup2.pdf", "Download PDF"),
                                                                downloadButton("sc1a2oup2.png", "Download PNG"), br(),
                                                                div(style="display:inline-block",
                                                                    numericInput("sc1a2oup2.h", "PDF / PNG height:", width = "138px",
                                                                                 min = 4, max = 20, value = 6, step = 0.5)),
                                                                div(style="display:inline-block",
                                                                    numericInput("sc1a2oup2.w", "PDF / PNG width:", width = "138px",
                                                                                 min = 4, max = 20, value = 8, step = 0.5))
                                                            )  # End of column (6 space)
                                                        )    # End of fluidRow (4 space)
                                                    ),     # End of tab (2 space) TABPANEL2
                                                  # Tabpanel 3 ##### 
                                                    ### Tab1.a3: geneExpr vs geneExpr on dimRed
                                                    tabPanel(
                                                        HTML("GeneExpr vs GeneExpr"),
                                                        h4("Gene expression vs gene expression on dimension reduction"),
                                                        "In this tab, users can visualise two gene expressions side-by-side ",
                                                        "on low-dimensional representions.",
                                                        br(),
                                                        br(),
                                                        fluidRow(
                                                            column(3, 
                                                                   h4("Dimension Reduction"),
                                                                   selectInput("sc1a3drX", "X-axis:", choices = sc1conf[dimred == TRUE]$UI,
                                                                                        selected = sc1def$dimred[1]),
                                                                    selectInput("sc1a3drY", "Y-axis:", choices = sc1conf[dimred == TRUE]$UI,
                                                                                    selected = sc1def$dimred[2])
                                                            ), # End of column (6 space)
                                                            column(3, 
                                                                   actionButton("sc1a3togL", "Toggle to subset cells"),
                                                                    conditionalPanel(
                                                                        condition = "input.sc1a3togL % 2 == 1",
                                                                        selectInput("sc1a3sub1", "Cell information to subset:",
                                                                                choices = sc1conf[grp == TRUE]$UI[-2],
                                                                                selected = sc1def$grp1),
                                                                        uiOutput("sc1a3sub1.ui"),
                                                                        actionButton("sc1a3sub1all", "Select all groups", class = "btn btn-primary"),
                                                                        actionButton("sc1a3sub1non", "Deselect all groups", class = "btn btn-primary")
                                                                    )
                                                            ), # End of column (6 space)
                                                            column(6, 
                                                                   actionButton("sc1a3tog0", "Toggle graphics controls"),
                                                                    conditionalPanel(
                                                                        condition = "input.sc1a3tog0 % 2 == 1",
                                                                    fluidRow(
                                                                        column(6, 
                                                                               sliderInput("sc1a3siz", "Point size:",
                                                                                           min = 0, max = 4, value = 1.25, step = 0.25),
                                                                                radioButtons("sc1a3psz", "Plot size:",
                                                                                         choices = c("Small", "Medium", "Large"),
                                                                                         selected = "Medium", inline = TRUE),
                                                                                radioButtons("sc1a3fsz", "Font size:",
                                                                                         choices = c("Small", "Medium", "Large"),
                                                                                         selected = "Medium", inline = TRUE)
                                                                        ),
                                                                        column(6, 
                                                                               radioButtons("sc1a3asp", "Aspect ratio:",
                                                                                            choices = c("Square", "Fixed", "Free"),
                                                                                            selected = "Square", inline = TRUE),
                                                                                checkboxInput("sc1a3txt", "Show axis text", value = FALSE)
                                                                        )
                                                                    )
                                                                )
                                                            )  # End of column (6 space)
                                                        ),   # End of fluidRow (4 space)
                                                        fluidRow(
                                                            column(6, 
                                                                   style="border-right: 2px solid black", h4("Gene expression 1"),
                                                                    fluidRow(
                                                                        column(6, 
                                                                               selectInput("sc1a3inp1", "Gene name:", choices=NULL) %>%
                                                                                helper(type = "inline", size = "m", fade = TRUE,
                                                                                       title = "Gene expression to colour cells by",
                                                                                       content = c("Select gene to colour cells by gene expression",
                                                                                                   paste0("- Gene expression are coloured in a ",
                                                                                                          "White-Red colour scheme which can be ",
                                                                                                          "changed in the plot controls")))
                                                                        ),
                                                                        column(6, 
                                                                               actionButton("sc1a3tog1", "Toggle plot controls"),
                                                                                conditionalPanel(
                                                                                    condition = "input.sc1a3tog1 % 2 == 1",
                                                                                    radioButtons("sc1a3col1", "Colour:",
                                                                                         choices = c("White-Red", "Blue-Yellow-Red",
                                                                                                     "Yellow-Green-Purple"),
                                                                                         selected = "White-Red"),
                                                                                    radioButtons("sc1a3ord1", "Plot order:",
                                                                                         choices = c("Max-1st", "Min-1st", "Original", "Random"),
                                                                                         selected = "Max-1st", inline = TRUE)
                                                                                )
                                                                        )
                                                                    ),
                                                                    fluidRow(column(12, uiOutput("sc1a3oup1.ui"))),
                                                                    downloadButton("sc1a3oup1.pdf", "Download PDF"),
                                                                    downloadButton("sc1a3oup1.png", "Download PNG"), br(),
                                                                    div(style="display:inline-block",
                                                                        numericInput("sc1a3oup1.h", "PDF / PNG height:", width = "138px",
                                                                                 min = 4, max = 20, value = 6, step = 0.5)),
                                                                    div(style="display:inline-block",
                                                                        numericInput("sc1a3oup1.w", "PDF / PNG width:", width = "138px",
                                                                                 min = 4, max = 20, value = 8, step = 0.5))
                                                            ), # End of column (6 space)
                                                            column(6, 
                                                                   h4("Gene expression 2"),
                                                                    fluidRow(
                                                                        column(6, 
                                                                               selectInput("sc1a3inp2", "Gene name:", choices=NULL) %>%
                                                                            helper(type = "inline", size = "m", fade = TRUE,
                                                                                   title = "Gene expression to colour cells by",
                                                                                   content = c("Select gene to colour cells by gene expression",
                                                                                               paste0("- Gene expression are coloured in a ",
                                                                                                      "White-Red colour scheme which can be ",
                                                                                                      "changed in the plot controls")))
                                                                        ),
                                                                        column(6, 
                                                                               actionButton("sc1a3tog2", "Toggle plot controls"),
                                                                                conditionalPanel(
                                                                                 condition = "input.sc1a3tog2 % 2 == 1",
                                                                                radioButtons("sc1a3col2", "Colour:",
                                                                                         choices = c("White-Red", "Blue-Yellow-Red",
                                                                                                     "Yellow-Green-Purple"),
                                                                                         selected = "White-Red"),
                                                                                radioButtons("sc1a3ord2", "Plot order:",
                                                                                             choices = c("Max-1st", "Min-1st", "Original", "Random"),
                                                                                         selected = "Max-1st", inline = TRUE)
                                                                                )
                                                                        )
                                                                    ),
                                                                    fluidRow(column(12, uiOutput("sc1a3oup2.ui"))),
                                                                    downloadButton("sc1a3oup2.pdf", "Download PDF"),
                                                                    downloadButton("sc1a3oup2.png", "Download PNG"), br(),
                                                                    div(style="display:inline-block",
                                                                        numericInput("sc1a3oup2.h", "PDF / PNG height:", width = "138px",
                                                                                 min = 4, max = 20, value = 6, step = 0.5)),
                                                                    div(style="display:inline-block",
                                                                        numericInput("sc1a3oup2.w", "PDF / PNG width:", width = "138px",
                                                                                 min = 4, max = 20, value = 8, step = 0.5))
                                                            )  # End of column (6 space)
                                                        )    # End of fluidRow (4 space)
                                                    ),     # End of tab (2 space) TABPANEL3
                                                  # Tabpanel 4 ##### 
                                                    ### Tab1.b2: Gene coexpression plot
                                                    tabPanel(
                                                        HTML("Gene coexpression"),
                                                        h4("Coexpression of two genes on reduced dimensions"),
                                                        "In this tab, users can visualise the coexpression of two genes ",
                                                        "on low-dimensional representions.",
                                                        br(),br(),
                                                        fluidRow(
                                                            column(
                                                                3, h4("Dimension Reduction"),
                                                                fluidRow(
                                                                    column(
                                                                        12, selectInput("sc1b2drX", "X-axis:", choices = sc1conf[dimred == TRUE]$UI,
                                                                                        selected = sc1def$dimred[1]),
                                                                        selectInput("sc1b2drY", "Y-axis:", choices = sc1conf[dimred == TRUE]$UI,
                                                                                    selected = sc1def$dimred[2]))
                                                                )
                                                            ), # End of column (6 space)
                                                            column(
                                                                3, actionButton("sc1b2togL", "Toggle to subset cells"),
                                                                conditionalPanel(
                                                                    condition = "input.sc1b2togL % 2 == 1",
                                                                    selectInput("sc1b2sub1", "Cell information to subset:",
                                                                                choices = sc1conf[grp == TRUE]$UI[-2],
                                                                                selected = sc1def$grp1),
                                                                    uiOutput("sc1b2sub1.ui"),
                                                                    actionButton("sc1b2sub1all", "Select all groups", class = "btn btn-primary"),
                                                                    actionButton("sc1b2sub1non", "Deselect all groups", class = "btn btn-primary")
                                                                )
                                                            ), # End of column (6 space)
                                                            column(
                                                                6, actionButton("sc1b2tog0", "Toggle graphics controls"),
                                                                conditionalPanel(
                                                                    condition = "input.sc1b2tog0 % 2 == 1",
                                                                    fluidRow(
                                                                        column(
                                                                            6, sliderInput("sc1b2siz", "Point size:",
                                                                                           min = 0, max = 4, value = 1.25, step = 0.25),
                                                                            radioButtons("sc1b2psz", "Plot size:",
                                                                                         choices = c("Small", "Medium", "Large"),
                                                                                         selected = "Medium", inline = TRUE),
                                                                            radioButtons("sc1b2fsz", "Font size:",
                                                                                         choices = c("Small", "Medium", "Large"),
                                                                                         selected = "Medium", inline = TRUE)
                                                                        ),
                                                                        column(
                                                                            6, radioButtons("sc1b2asp", "Aspect ratio:",
                                                                                            choices = c("Square", "Fixed", "Free"),
                                                                                            selected = "Square", inline = TRUE),
                                                                            checkboxInput("sc1b2txt", "Show axis text", value = FALSE)
                                                                        )
                                                                    )
                                                                )
                                                            )  # End of column (6 space)
                                                        ),   # End of fluidRow (4 space)
                                                        fluidRow(
                                                            column(
                                                                3, style="border-right: 2px solid black", h4("Gene Expression"),
                                                                selectInput("sc1b2inp1", "Gene 1:", choices=NULL) %>%
                                                                    helper(type = "inline", size = "m", fade = TRUE,
                                                                           title = "Gene expression to colour cells by",
                                                                           content = c("Select gene to colour cells by gene expression",
                                                                                       paste0("- Gene expression are coloured in a ",
                                                                                              "White-Red colour scheme which can be ",
                                                                                              "changed in the plot controls"))),
                                                                selectInput("sc1b2inp2", "Gene 2:", choices=NULL) %>%
                                                                    helper(type = "inline", size = "m", fade = TRUE,
                                                                           title = "Gene expression to colour cells by",
                                                                           content = c("Select gene to colour cells by gene expression",
                                                                                       paste0("- Gene expression are coloured in a ",
                                                                                              "White-Blue colour scheme which can be ",
                                                                                              "changed in the plot controls"))),
                                                                actionButton("sc1b2tog1", "Toggle plot controls"),
                                                                conditionalPanel(
                                                                    condition = "input.sc1b2tog1 % 2 == 1",
                                                                    radioButtons("sc1b2col1", "Colour:",
                                                                                 choices = c("Red (Gene1); Blue (Gene2)",
                                                                                             "Orange (Gene1); Blue (Gene2)",
                                                                                             "Red (Gene1); Green (Gene2)",
                                                                                             "Green (Gene1); Blue (Gene2)"),
                                                                                 selected = "Red (Gene1); Blue (Gene2)"),
                                                                    radioButtons("sc1b2ord1", "Plot order:",
                                                                                 choices = c("Max-1st", "Min-1st", "Original", "Random"),
                                                                                 selected = "Max-1st", inline = TRUE)
                                                                )
                                                            ), # End of column (6 space)
                                                            column(
                                                                6, style="border-right: 2px solid black",
                                                                uiOutput("sc1b2oup1.ui"),
                                                                downloadButton("sc1b2oup1.pdf", "Download PDF"),
                                                                downloadButton("sc1b2oup1.png", "Download PNG"), br(),
                                                                div(style="display:inline-block",
                                                                    numericInput("sc1b2oup1.h", "PDF / PNG height:", width = "138px",
                                                                                 min = 4, max = 20, value = 8, step = 0.5)),
                                                                div(style="display:inline-block",
                                                                    numericInput("sc1b2oup1.w", "PDF / PNG width:", width = "138px",
                                                                                 min = 4, max = 20, value = 10, step = 0.5))
                                                            ), # End of column (6 space)
                                                            column(
                                                                3, uiOutput("sc1b2oup2.ui"),
                                                                downloadButton("sc1b2oup2.pdf", "Download PDF"),
                                                                downloadButton("sc1b2oup2.png", "Download PNG"),
                                                                br(), h4("Cell numbers"),
                                                                dataTableOutput("sc1b2.dt")
                                                            )  # End of column (6 space)
                                                        )    # End of fluidRow (4 space)
                                                    ),     # End of tab (2 space) TABPANEL4
                                                  # Tabpanel 5 ##### 
                                                    ### Tab1.c1: violinplot / boxplot
                                                    tabPanel(
                                                        HTML("Violinplot / Boxplot"),
                                                        h4("Cell information / gene expression violin plot / box plot"),
                                                        "In this tab, users can visualise the gene expression or continuous cell information ",
                                                        "(e.g. Number of UMIs / module score) across groups of cells (e.g. libary / clusters).",
                                                        br(),br(),
                                                        fluidRow(
                                                            column(
                                                                3, style="border-right: 2px solid black",
                                                                selectInput("sc1c1inp1", "Cell information (X-axis):",
                                                                            choices = sc1conf[grp == TRUE]$UI[-2],
                                                                            selected = sc1def$grp1) %>%
                                                                    helper(type = "inline", size = "m", fade = TRUE,
                                                                           title = "Cell information to group cells by",
                                                                           content = c("Select categorical cell information to group cells by",
                                                                                       "- Single cells are grouped by this categorical covariate",
                                                                                       "- Plotted as the X-axis of the violin plot / box plot")),
                                                                selectInput("sc1c1inp2", "Cell Info / Gene name (Y-axis):", choices=NULL) %>%
                                                                    helper(type = "inline", size = "m", fade = TRUE,
                                                                           title = "Cell Info / Gene to plot",
                                                                           content = c("Select cell info / gene to plot on Y-axis",
                                                                                       "- Can be continuous cell information (e.g. nUMIs / scores)",
                                                                                       "- Can also be gene expression")),
                                                                radioButtons("sc1c1typ", "Plot type:",
                                                                             choices = c("violin", "boxplot"),
                                                                             selected = "violin", inline = TRUE),
                                                                checkboxInput("sc1c1pts", "Show data points", value = FALSE),
                                                                actionButton("sc1c1togL", "Toggle to subset cells"),
                                                                conditionalPanel(
                                                                    condition = "input.sc1c1togL % 2 == 1",
                                                                    selectInput("sc1c1sub1", "Cell information to subset:",
                                                                                choices = sc1conf[grp == TRUE]$UI[-2],
                                                                                selected = sc1def$grp1),
                                                                    uiOutput("sc1c1sub1.ui"),
                                                                    actionButton("sc1c1sub1all", "Select all groups", class = "btn btn-primary"),
                                                                    actionButton("sc1c1sub1non", "Deselect all groups", class = "btn btn-primary")
                                                                ), br(), br(),
                                                                actionButton("sc1c1tog", "Toggle graphics controls"),
                                                                conditionalPanel(
                                                                    condition = "input.sc1c1tog % 2 == 1",
                                                                    sliderInput("sc1c1siz", "Data point size:",
                                                                                min = 0, max = 4, value = 1.25, step = 0.25),
                                                                    radioButtons("sc1c1psz", "Plot size:",
                                                                                 choices = c("Small", "Medium", "Large"),
                                                                                 selected = "Medium", inline = TRUE),
                                                                    radioButtons("sc1c1fsz", "Font size:",
                                                                                 choices = c("Small", "Medium", "Large"),
                                                                                 selected = "Medium", inline = TRUE))
                                                            ), # End of column (6 space)
                                                            column(9, uiOutput("sc1c1oup.ui"),
                                                                   downloadButton("sc1c1oup.pdf", "Download PDF"),
                                                                   downloadButton("sc1c1oup.png", "Download PNG"), br(),
                                                                   div(style="display:inline-block",
                                                                       numericInput("sc1c1oup.h", "PDF / PNG height:", width = "138px",
                                                                                    min = 4, max = 20, value = 8, step = 0.5)),
                                                                   div(style="display:inline-block",
                                                                       numericInput("sc1c1oup.w", "PDF / PNG width:", width = "138px",
                                                                                    min = 4, max = 20, value = 10, step = 0.5))
                                                            )  # End of column (6 space)
                                                        )    # End of fluidRow (4 space)
                                                    ),     # End of tab (2 space) 
                                                  # Tabpanel 6 ##### 
                                                    ### Tab1.c2: Proportion plot
                                                    tabPanel(
                                                        HTML("Proportion plot"),
                                                        h4("Proportion / cell numbers across different cell information"),
                                                        "In this tab, users can visualise the composition of single cells based on one discrete ",
                                                        "cell information across another discrete cell information. ",
                                                        "Usage examples include the library or cellcycle composition across clusters.",
                                                        br(),br(),
                                                        fluidRow(
                                                            column(
                                                                3, style="border-right: 2px solid black",
                                                                selectInput("sc1c2inp1", "Cell information to plot (X-axis):",
                                                                            choices = sc1conf[grp == TRUE]$UI[-2],
                                                                            selected = sc1def$grp2) %>%
                                                                    helper(type = "inline", size = "m", fade = TRUE,
                                                                           title = "Cell information to plot cells by",
                                                                           content = c("Select categorical cell information to plot cells by",
                                                                                       "- Plotted as the X-axis of the proportion plot")),
                                                                selectInput("sc1c2inp2", "Cell information to group / colour by:",
                                                                            choices = sc1conf[grp == TRUE]$UI[-2],
                                                                            selected = sc1def$grp1) %>%
                                                                    helper(type = "inline", size = "m", fade = TRUE,
                                                                           title = "Cell information to group / colour cells by",
                                                                           content = c("Select categorical cell information to group / colour cells by",
                                                                                       "- Proportion / cell numbers are shown in different colours")),
                                                                radioButtons("sc1c2typ", "Plot value:",
                                                                             choices = c("Proportion", "CellNumbers"),
                                                                             selected = "Proportion", inline = TRUE),
                                                                checkboxInput("sc1c2flp", "Flip X/Y", value = FALSE),
                                                                actionButton("sc1c2togL", "Toggle to subset cells"),
                                                                conditionalPanel(
                                                                    condition = "input.sc1c2togL % 2 == 1",
                                                                    selectInput("sc1c2sub1", "Cell information to subset:",
                                                                                choices = sc1conf[grp == TRUE]$UI[-2],
                                                                                selected = sc1def$grp1),
                                                                    uiOutput("sc1c2sub1.ui"),
                                                                    actionButton("sc1c2sub1all", "Select all groups", class = "btn btn-primary"),
                                                                    actionButton("sc1c2sub1non", "Deselect all groups", class = "btn btn-primary")
                                                                ), br(), br(),
                                                                actionButton("sc1c2tog", "Toggle graphics controls"),
                                                                conditionalPanel(
                                                                    condition = "input.sc1c2tog % 2 == 1",
                                                                    radioButtons("sc1c2psz", "Plot size:",
                                                                                 choices = c("Small", "Medium", "Large"),
                                                                                 selected = "Medium", inline = TRUE),
                                                                    radioButtons("sc1c2fsz", "Font size:",
                                                                                 choices = c("Small", "Medium", "Large"),
                                                                                 selected = "Medium", inline = TRUE))
                                                            ), # End of column (6 space)
                                                            column(9, uiOutput("sc1c2oup.ui"),
                                                                   downloadButton("sc1c2oup.pdf", "Download PDF"),
                                                                   downloadButton("sc1c2oup.png", "Download PNG"), br(),
                                                                   div(style="display:inline-block",
                                                                       numericInput("sc1c2oup.h", "PDF / PNG height:", width = "138px",
                                                                                    min = 4, max = 20, value = 8, step = 0.5)),
                                                                   div(style="display:inline-block",
                                                                       numericInput("sc1c2oup.w", "PDF / PNG width:", width = "138px",
                                                                                    min = 4, max = 20, value = 10, step = 0.5))
                                                            )  # End of column (6 space)
                                                        )    # End of fluidRow (4 space)
                                                    ),     # End of tab (2 space) TABPANEL6 
                                                  # Tabpanel 7 ##### 
                                                    ### Tab1.d1: Multiple gene expr
                                                    tabPanel(
                                                        HTML("Bubbleplot / Heatmap"),
                                                        h4("Gene expression bubbleplot / heatmap"),
                                                        "In this tab, users can visualise the gene expression patterns of ",
                                                        "multiple genes grouped by categorical cell information (e.g. library / cluster).", 
                                                        br(),
                                                        "The normalised expression are averaged, log-transformed and then plotted.",
                                                        br(),
                                                        br(),
                                                        fluidRow(
                                                            column(3, 
                                                                   style="border-right: 2px solid black",
                                                                    textAreaInput("sc1d1inp", 
                                                                                 HTML("List of gene names <br /> (Max 50 genes, separated <br /> by , or ; or newline):"),
                                                                                height = "200px",
                                                                                value = paste0(sc1def$genes, collapse = ", ")) %>%
                                                                    helper(type = "inline", size = "m", fade = TRUE,
                                                                           title = "List of genes to plot on bubbleplot / heatmap",
                                                                           content = c("Input genes to plot",
                                                                                       "- Maximum 50 genes (due to ploting space limitations)",
                                                                                       "- Genes should be separated by comma, semicolon or newline")),
                                                                    selectInput("sc1d1grp", "Group by:",
                                                                                choices = sc1conf[grp == TRUE]$UI[-2],
                                                                                selected = sc1conf[grp == TRUE]$UI[1]) %>%
                                                                    helper(type = "inline", size = "m", fade = TRUE,
                                                                               title = "Cell information to group cells by",
                                                                               content = c("Select categorical cell information to group cells by",
                                                                                           "- Single cells are grouped by this categorical covariate",
                                                                                           "- Plotted as the X-axis of the bubbleplot / heatmap")),
                                                                    radioButtons("sc1d1plt", "Plot type:",
                                                                                 choices = c("Bubbleplot", "Heatmap"),
                                                                                 selected = "Bubbleplot", inline = TRUE),
                                                                    checkboxInput("sc1d1scl", "Scale gene expression", value = TRUE),
                                                                    checkboxInput("sc1d1row", "Cluster rows (genes)", value = TRUE),
                                                                    checkboxInput("sc1d1col", "Cluster columns (samples)", value = FALSE),
                                                                    br(),
                                                                    actionButton("sc1d1togL", "Toggle to subset cells"),
                                                                    conditionalPanel(
                                                                        condition = "input.sc1d1togL % 2 == 1",
                                                                        selectInput("sc1d1sub1", "Cell information to subset:",
                                                                                    choices = sc1conf[grp == TRUE]$UI[-2],
                                                                                    selected = sc1def$grp1),
                                                                        uiOutput("sc1d1sub1.ui"),
                                                                        actionButton("sc1d1sub1all", "Select all groups", class = "btn btn-primary"),
                                                                        actionButton("sc1d1sub1non", "Deselect all groups", class = "btn btn-primary")
                                                                    ), 
                                                                    br(), 
                                                                    br(),
                                                                    actionButton("sc1d1tog", "Toggle graphics controls"),
                                                                    conditionalPanel(
                                                                        condition = "input.sc1d1tog % 2 == 1",
                                                                        radioButtons("sc1d1cols", "Colour scheme:",
                                                                                     choices = c("White-Red", "Blue-Yellow-Red",
                                                                                                 "Yellow-Green-Purple"),
                                                                                     selected = "Blue-Yellow-Red"),
                                                                        radioButtons("sc1d1psz", "Plot size:",
                                                                                     choices = c("Small", "Medium", "Large"),
                                                                                     selected = "Medium", inline = TRUE),
                                                                        radioButtons("sc1d1fsz", "Font size:",
                                                                                     choices = c("Small", "Medium", "Large"),
                                                                                     selected = "Medium", inline = TRUE)
                                                                    )
                                                            ), # End of column (6 space)
                                                            column(9, 
                                                                   h4(htmlOutput("sc1d1oupTxt")),
                                                                   uiOutput("sc1d1oup.ui"),
                                                                   downloadButton("sc1d1oup.pdf", "Download PDF"),
                                                                   downloadButton("sc1d1oup.png", "Download PNG"), 
                                                                   br(),
                                                                   div(style="display:inline-block",
                                                                       numericInput("sc1d1oup.h", "PDF / PNG height:", width = "138px",
                                                                                    min = 4, max = 20, value = 10, step = 0.5)),
                                                                   div(style="display:inline-block",
                                                                       numericInput("sc1d1oup.w", "PDF / PNG width:", width = "138px",
                                                                                    min = 4, max = 20, value = 10, step = 0.5))
                                                            )  # End of column (6 space)
                                                        )    # End of fluidRow (4 space)
                                                    )      # End of tab (2 space) TABPANEL7 
                                                               #####
                                        ) # tabsetpanel 2 ends
                            ), # end tabpanel 2
                            tabPanel(HTML("Neural Network for Treatment Prediction in Multiple Myeloma"),
                                     style = "margin-bottom: 4.5rem; width: 100%;",
                                     # Main panel for displaying outputs ----
                                     sidebarLayout(
                                         sidebarPanel(
                                             useShinyjs(),
                                             style = "width: 125%;",
                                             tabsetPanel(
                                                 tabPanel("RNA-seq Data",
                                                          # Content for RNA-seq Data tab
                                                          tags$img(
                                                              src = "single-cell_RNA_seq_protocol.png",
                                                              height = "auto",
                                                              width = "100%",
                                                              alt = "Something went wrong",
                                                              style = "margin-top: 20px; margin-bottom: 20px;"
                                                          ),
                                                          tags$div(
                                                              style = "background-color:#f5f5f5; border-left: 3px solid #2a78d6; padding: 8px 12px; margin-bottom: 12px; font-size: 12px; color: #333;",
                                                              tags$strong("Data privacy: "),
                                                              "Uploaded files are processed in memory for this session only, to generate a prediction. ",
                                                              "They are not stored, logged, or shared with third parties. Please do not upload data containing direct patient identifiers."
                                                          ),
                                                          downloadButton("downloadExampleData_bulk", "Download example data") %>%
                                                              helper(
                                                                  type = "inline", size = "m", fade = TRUE,
                                                                  title = "Example data",
                                                                  content = "A raw bulk RNA-seq expression matrix for 6 real patients from an independent, publicly available external
                                                                             validation cohort (GEO accession GSE159426, PAD regimen), used in the associated study. Upload it here to
                                                                             try the full preprocessing workflow before using your own data."),
                                                          br(), br(),
                                                          fileInput("RNA_seq", label = "Input an Expression Raw Matrix from RNA-seq (csv)",
                                                                    accept = c('text/csv', 'text/comma-separated-values', 'text/plain', '.csv')
                                                                    ),
                                                          actionButton("showData_raw", "Show Head of processed RNA-seq Data Uploaded"),
                                                          actionButton("removeRawTable", "Remove RNA-seq Raw Data"),
                                                          downloadButton("downloadData", "Download Processed RNA-seq Data for the next step"),
                                                          br(), br(), br()
                                                 ),
                                                 tabPanel("Single-cell RNA-seq Data",
                                                          # Content for Single-cell RNA-seq Data tab
                                                          tags$img(
                                                              src = "RNA_seq_protocol.png",
                                                              height = "auto",
                                                              width = "100%",
                                                              alt = "Something went wrong",
                                                              style = "margin-top: 20px; margin-bottom: 20px;"
                                                          ),
                                                          tags$div(
                                                              style = "background-color:#f5f5f5; border-left: 3px solid #2a78d6; padding: 8px 12px; margin-bottom: 12px; font-size: 12px; color: #333;",
                                                              tags$strong("Data privacy: "),
                                                              "Uploaded files are processed in memory for this session only, to generate a prediction. ",
                                                              "They are not stored, logged, or shared with third parties. Please do not upload data containing direct patient identifiers."
                                                          ),
                                                          downloadButton("downloadExampleData_sc", "Download example data") %>%
                                                              helper(
                                                                  type = "inline", size = "m", fade = TRUE,
                                                                  title = "Example data",
                                                                  content = "The same 6-patient expression matrix as in the RNA-seq Data tab (GEO accession GSE159426), usable here to
                                                                             try the single-cell-style preprocessing workflow before using your own data."),
                                                          br(), br(),
                                                          fileInput("Single_cell", label = "Input an Expression Raw Matrix from Single-cell RNA-seq (csv)",
                                                                    accept = c('text/csv', 'text/comma-separated-values', 'text/plain', '.csv')),
                                                          actionButton("showData_raw_sc", "Show Head of processed scRNA-seq Data Uploaded") %>%
                                                              helper(
                                                                  type = "inline", size = "m", fade = TRUE,
                                                                  title = "Information about Data Visualization",
                                                                  content = "This process may take some time, please wait"),
                                                          actionButton("removeRawTable_sc", "Remove Raw Data"),
                                                          downloadButton("downloadData_sc", "Download Processed scRNA-seq Data for the next step"),
                                                          br(), br(), br()
                                                 ),
                                                 tabPanel("Neural Network Prediction",
                                                          # Content for Neural Network Prediction tab
                                                          tags$div(
                                                              style = "background-color:#f5f5f5; border-left: 3px solid #2a78d6; padding: 8px 12px; margin-bottom: 12px; font-size: 12px; color: #333;",
                                                              tags$strong("Data privacy: "),
                                                              "Uploaded files are processed in memory for this session only, to generate a prediction. ",
                                                              "They are not stored, logged, or shared with third parties. Please do not upload data containing direct patient identifiers."
                                                          ),
                                                          downloadButton("downloadExampleData_prediction", "Download example data (known-outcome patients)") %>%
                                                              helper(
                                                                  type = "inline", size = "m", fade = TRUE,
                                                                  title = "Example data",
                                                                  content = "Already-preprocessed data (skip the two tabs on the left) for 2 real patients from the external
                                                                             validation cohort GSE159426, with known clinical outcome: MM104 (optimal responder) and MM116
                                                                             (suboptimal responder). Upload this file and click 'Run Neural Network' to see whether the
                                                                             predictions match these known labels."),
                                                          br(), br(),
                                                          fileInput("patient_data_expression", label = "Input the previous Downloaded Data for Neural Network Prediction",
                                                                    accept = c('text/csv', 'text/comma-separated-values', 'text/plain', '.csv')),
                                                          actionButton("runNetwork", "Run Neural Network") %>%
                                                              helper(
                                                                  type = "inline", size = "m", fade = TRUE,
                                                                  title = "Information about Neural Network Fundamentals",
                                                                  content = "This Neural Network was designed to predict response to: 
                                                                  Bortezomib, Melphalan and Prednisolone (VMP); 
                                                                  Bortezomib, Thalidomide and Dexamethasone (VTD);
                                                                  Bortezomib, Adriamycin and Dexamethasone (PAD)"),
                                                          actionButton("showData", "Show Data Uploaded"),
                                                          actionButton("removeTable", "Remove Table"),
                                                          tags$img(
                                                              src = "Screenshot_neural.png",
                                                              height = "auto",
                                                              width = "100%",
                                                              alt = "Something went wrong",
                                                              style = "margin-top: 20px;"
                                                          )
                                                        )
                                             )
                                         ),
                                         mainPanel(
                                             div(
                                                 verbatimTextOutput(outputId = "prediction") %>%
                                                 helper(
                                                     type = "inline", size = "m", fade = TRUE,
                                                     title = "How to interpret your result",
                                                     content = "For each sample, the network outputs a probability of belonging to the OPTIMAL-response group and a
                                                                probability of belonging to the SUBOPTIMAL-response group (the two always sum to 100%).
                                                                The predicted label is whichever group has the higher probability, unless that probability is
                                                                below 70%, in which case the sample is reported as UNDETERMINED rather than forcing a low-confidence call.
                                                                OPTIMAL indicates the model expects a favourable response to first-line bortezomib-based therapy
                                                                (VMP / VTD / PAD); SUBOPTIMAL indicates an expected poor response. This is a research-use prediction
                                                                tool, not a diagnostic device, and should not be used alone to guide clinical decisions."),
                                                 style = "max-width: 100%;"
                                                 ),
                                             div(
                                                 id = "rawTableContainer",
                                                 DTOutput('table_raw', width = "100%"),
                                                 style = "margin-top: 40px;"
                                             ),
                                             div(
                                                 id = "rawTableContainer_sc",
                                                 DTOutput('table_raw_sc', width = "100%"),
                                                 style = "margin-top: 40px;"
                                             ),
                                             div(
                                                 id = "tableContainer",
                                                 DTOutput('table', width = "100%"),
                                                 style = "margin-top: 40px;"
                                             ),
                                             style = "margin-top: 20px; margin-left: 160px; max-width: 50%;"  # Adjust the max-width as needed
                                         )
                                     )
                            )
                ), # tabsetpanel end
                # Static footer
                tags$div(
                    style = "background-color: #f5f5f5; padding: 10px; text-align: center; position: fixed; bottom: 0; left: 0; right: 0; min-height: 2.5rem; width: 100%;",
                    HTML(paste("<strong>Enrique de la Rosa Morón, Bioinformatics & Functional Genomics</strong>", "<br><span style='font-size: 10px; padding-top: 5px;'>Ouyang et
                               al. ShinyCell: Simple and sharable visualisation of single-cell gene expression data. Bioinformatics, doi:10.1093/bioinformatics/btab209</span>",
                               "<br><span style='font-size: 10px; color:#666;'>This application uses only strictly necessary session cookies required for it to function
                               (no third-party tracking or analytics cookies are set).</span>"))
                )
)


#############################################################################################################
#############################################################################################################
#############################################################################################################
#############################################################################################################

# Define server  ----
#####
server <- function(input, output, session) {
#####    
    shinyjs::useShinyjs()
    
    options(shiny.maxRequestSize = 3000 * 1024^2)
    
    observe(session$setCurrentTheme(
        if (isTRUE(input$dark_mode)) dark else light
    ))
    
    output$pdfview <- renderUI({
        tags$iframe(style="height:1000px; width:150%", src="single_cell_script_tutorial.html")
    })
    
#####
    ######
    # Function to reset the file input buttons
    resetFileInputs <- function() {
        shinyjs::reset("RNA_seq")
        shinyjs::reset("Single_cell")
        shinyjs::reset("patient_data_expression")
    }
    
    # Initial visibility of tables
    showRawTable <- FALSE
    showRawTable_sc <- FALSE
    showTable <- FALSE
    
    # Logic to control datatable visibility
    observeEvent(input$removeRawTable, {
        showRawTable <<- FALSE
        shinyjs::hide("rawTableContainer")
        resetFileInputs()
      #  session$reload()  # Reload the session
    })
    
    observeEvent(input$removeRawTable_sc, {
        showRawTable_sc <<- FALSE
        shinyjs::hide("rawTableContainer_sc")
        resetFileInputs()
       # session$reload()  # Reload the session
    })
    
    observeEvent(input$removeTable, {
        showTable <<- FALSE
        shinyjs::hide("tableContainer")
        resetFileInputs()
        #session$reload()  # Reload the session
    })
    
    observeEvent(input$RNA_seq, {
        showRawTable <<- TRUE
        shinyjs::show("rawTableContainer")
    })
    
    observeEvent(input$Single_cell, {
        showRawTable_sc <<- TRUE
        shinyjs::show("rawTableContainer_sc")
    })
    
    observeEvent(input$patient_data_expression, {
        showTable <<- TRUE
        shinyjs::show("tableContainer")
    })
    
    observeEvent(input$showData_raw_sc, {
        # Show a modal when the button is pressed
        shinyalert("Wait please!", "Data Processing is in progress...")
    })
    
    observeEvent(input$showData_raw, {
        # Show a modal when the button is pressed
        shinyalert("Wait please!", "Data Processing is in progress...")
    })
    
    UploadandprocessData <- function(inputfile, genes_neural, my_matrix_filled_merged) {
        data_raw <- read.csv(inputfile$datapath)
        rownames(data_raw) <- data_raw[,1]
        data_raw <- data_raw[,-1]
        data_raw <- as.matrix(data_raw)
        
        missing_genes <- setdiff(genes_neural, rownames(data_raw))
        
        if (length(missing_genes) > 0) {
            to_add <- matrix(data = 0, nrow = length(missing_genes), ncol = ncol(data_raw))
            data_raw <- rbind(data_raw, to_add)
            rownames(data_raw)[(nrow(data_raw) - length(missing_genes) + 1):nrow(data_raw)] <- missing_genes
        } 
        
        
        data_raw <- data_raw[genes_neural, ]
        data_raw <- cbind(my_matrix_filled_merged, data_raw)
        data_raw <- CreateSeuratObject(data_raw)
        data_raw <- NormalizeData(data_raw)
        data_raw <- FindVariableFeatures(data_raw, selection.method = "vst", nfeatures = 38651)
        data_raw <- ScaleData(data_raw, features = rownames(data_raw), scale.max = 10)
        
        return(data_raw)
    }
    
    data_matrix_raw <- reactive({
        req(input$RNA_seq)
        data_raw <- UploadandprocessData(input$RNA_seq, genes_neural, my_matrix_filled_merged)
    })
    
    output$table_raw <- renderDT({
        req(input$RNA_seq, input$showData_raw > 0)
        data_raw <- data_matrix_raw()

        if(ncol(data_raw[,19:ncol(data_raw)]) < 4){
            head(as.matrix(data_raw@assays$RNA@scale.data[, 19:ncol(data_raw)]))
        } else {
            head(as.matrix(data_raw@assays$RNA@scale.data[, 19:21]))
        }
    })
    
    processData <- function(inputfile, genes_neural, my_matrix_filled_merged) {
        data_raw_sc <- read.csv(inputfile$datapath)
        rownames(data_raw_sc) <- data_raw_sc[,1]
        data_raw_sc <- data_raw_sc[,-1]
        data_raw_sc <- as.matrix(data_raw_sc)
        
        missing_genes <- setdiff(genes_neural, rownames(data_raw_sc))
        
        if (length(missing_genes) > 0) {
            to_add <- matrix(data = 0, nrow = length(missing_genes), ncol = ncol(data_raw_sc))
            data_raw_sc <- rbind(data_raw_sc, to_add)
            rownames(data_raw_sc)[(nrow(data_raw_sc) - length(missing_genes) + 1):nrow(data_raw_sc)] <- missing_genes
        } 
        
        
        data_raw_sc <- data_raw_sc[genes_neural, ]
        data_raw_sc <- rowMeans(data_raw_sc)
        data_raw_sc <- cbind(my_matrix_filled_merged, data_raw_sc)
        data_raw_sc <- CreateSeuratObject(data_raw_sc)
        data_raw_sc <- NormalizeData(data_raw_sc)
        data_raw_sc <- FindVariableFeatures(data_raw_sc, selection.method = "vst", nfeatures = 38651)
        data_raw_sc <- ScaleData(data_raw_sc, features = rownames(data_raw_sc), scale.max = 10)
        
        return(data_raw_sc)
    }
    
    data_matrix_raw_sc <- reactive({
        req(input$Single_cell)
        Sys.sleep(5)
        data_raw_sc <- processData(input$Single_cell, genes_neural, my_matrix_filled_merged)
    })
    
    output$table_raw_sc <- renderDT({
        req(input$Single_cell, input$showData_raw_sc > 0)
        data_raw_sc <- data_matrix_raw_sc()
        
        head(as.matrix(data_raw_sc@assays$RNA@scale.data[1:5, 19]))
    })
    
    data_matrix <- reactive({
        req(input$patient_data_expression)
        data <- read.csv(input$patient_data_expression$datapath)
        rownames(data) <- data[,1]
        data <- data[,-1]
    })
    
    output$table <- renderDT({
        req(input$patient_data_expression)
        data <- read.csv(input$patient_data_expression$datapath)
        rownames(data) <- data[,1]
        data <- data[,-1]
        data <- as.matrix(data)
        
        if (req(input$showData > 0)) {
            if (nrow(data) < 4) {
                t(data[nrow(data),1:4])
            }
            else{
                t(data[1:3,1:4])
            } 
        } else {
            NULL
        }
    })
    #####

    #####
    output$prediction <- renderText({
        if (is.null(data_matrix())) {
            return("Please upload a valid expression matrix.")
        }
        
        mtx <- as.matrix(data_matrix())
        # Check if the "Run Neural Network" button is clicked
        run_clicked <- input$runNetwork
        
        if (run_clicked == 0) {
            return("Click 'Run Neural Network' to make predictions.")
        } else {
            # Make predictions using the loaded model
            prediction <- predict(model_NN, mtx)
            prediction_result <- c()
            # browser()
            for (i in seq(1:nrow(prediction))) {
                if (max(prediction[i,]) == prediction[[i,1]]) {
                    if(max(prediction[i,]) < 0.7) {
                    prediction_text <- paste0("undetermined (", round(prediction[[i,1]] * 100, 2), "%)")    
                    } else {
                    prediction_text <- paste0("optimal (", round(prediction[[i,1]] * 100, 2), "%)")
                    }
                } else {
                    if(max(prediction[i,]) < 0.7){
                    prediction_text <- paste0("undetermined (", round(prediction[[i,2]] * 100, 2), "%)")
                    } else {
                    prediction_text <- paste0("suboptimal (", round(prediction[[i,2]] * 100, 2), "%)") 
                    }    
                } 
                prediction_result <- rbind(prediction_result, prediction_text)
            }

            # Return the prediction
            # Join the lines with line breaks after a specific character count (e.g., 50 characters)
            line_break_width <- 80
            response_strings <- paste0("The predicted treatment response of the given patient ", seq(1:nrow(prediction)), " is ", toupper(prediction_result))
            response_strings_wrapped <- str_wrap(response_strings, width = line_break_width)
            paste(response_strings_wrapped, collapse = "\n")
        }
    })

    # Downloadable csv of selected dataset ----
    output$downloadData <- downloadHandler(
        filename = function() {
            paste("preprocessed_dataset_for_neural_network.csv", sep = "")
        },
        content = function(file) {
            data_raw_seurat <- data_matrix_raw()
            write.csv(as.matrix(t(data_raw_seurat@assays$RNA@scale.data[, 19:ncol(data_raw_seurat)])), file, row.names = TRUE)
        }
    )
    
    output$downloadData_sc <- downloadHandler(
        filename = function() {
            paste("preprocessed_dataset_for_neural_network.csv", sep = "")
        },
        content = function(file) {
            data_raw_seurat <- data_matrix_raw_sc()
            write.csv(as.matrix(t(data_raw_seurat@assays$RNA@scale.data[, 19:ncol(data_raw_seurat)])), file, row.names = TRUE)
        }
    )

    # Example data, so first-time users can try the full workflow before uploading their own files.
    # Real, public, external-validation-cohort data (GSE159426, 6 patients, PAD regimen).
    output$downloadExampleData_bulk <- downloadHandler(
        filename = function() {
            "example_raw_expression_matrix_GSE159426.csv"
        },
        content = function(file) {
            file.copy("www/example_raw_expression_matrix.csv", file)
        }
    )

    output$downloadExampleData_sc <- downloadHandler(
        filename = function() {
            "example_raw_expression_matrix_GSE159426.csv"
        },
        content = function(file) {
            file.copy("www/example_raw_expression_matrix.csv", file)
        }
    )

    # Already-preprocessed example (2 patients, known clinical outcome) for the
    # Neural Network Prediction tab, so users can validate predictions directly.
    output$downloadExampleData_prediction <- downloadHandler(
        filename = function() {
            "example_prediction_ready_MM104optimal_MM116suboptimal.csv"
        },
        content = function(file) {
            file.copy("www/example_prediction_ready.csv", file)
        }
    )
    
    
    #####
    #####
    ### For all tags and Server-side selectize 
    observe_helpers() 
    optCrt="{ option_create: function(data,escape) {return('<div class=\"create\"><strong>' + '</strong></div>');} }" 
    updateSelectizeInput(session, "sc1a1inp2", choices = names(sc1gene)[which(names(sc1gene) %in% genes_neural)], server = TRUE, 
                         selected = names(sc1gene)[566], options = list( 
                             maxOptions = 77, create = TRUE, persist = TRUE, render = I(optCrt))) 
    updateSelectizeInput(session, "sc1a3inp1", choices = names(sc1gene)[which(names(sc1gene) %in% genes_neural)], server = TRUE, 
                         selected = names(sc1gene)[566], options = list( 
                             maxOptions = 77, create = TRUE, persist = TRUE, render = I(optCrt))) 
    updateSelectizeInput(session, "sc1a3inp2", choices = names(sc1gene)[which(names(sc1gene) %in% genes_neural)], server = TRUE, 
                         selected = names(sc1gene)[19310], options = list( 
                             maxOptions = 77, create = TRUE, persist = TRUE, render = I(optCrt))) 
    updateSelectizeInput(session, "sc1b2inp1", choices = names(sc1gene)[which(names(sc1gene) %in% genes_neural)], server = TRUE, 
                         selected = names(sc1gene)[566], options = list( 
                             maxOptions = 77, create = TRUE, persist = TRUE, render = I(optCrt))) 
    updateSelectizeInput(session, "sc1b2inp2", choices = names(sc1gene)[which(names(sc1gene) %in% genes_neural)], server = TRUE, 
                         selected = names(sc1gene)[19310], options = list( 
                             maxOptions = 77, create = TRUE, persist = TRUE, render = I(optCrt))) 
    updateSelectizeInput(session, "sc1c1inp2", server = TRUE, 
                         choices = c(sc1conf[is.na(fID)]$UI,names(sc1gene)[which(names(sc1gene) %in% genes_neural)]), 
                         selected = sc1conf[is.na(fID)]$UI[1], options = list( 
                             maxOptions = length(sc1conf[is.na(fID)]$UI) + 3, 
                             create = TRUE, persist = TRUE, render = I(optCrt))) 
    
    
    ### Plots for tab a1 
    output$sc1a1sub1.ui <- renderUI({ 
        sub = strsplit(sc1conf[UI == input$sc1a1sub1]$fID, "\\|")[[1]] 
        checkboxGroupInput("sc1a1sub2", "Select which cells to show", inline = TRUE, 
                           choices = sub, selected = sub) 
    }) 
    observeEvent(input$sc1a1sub1non, { 
        sub = strsplit(sc1conf[UI == input$sc1a1sub1]$fID, "\\|")[[1]] 
        updateCheckboxGroupInput(session, inputId = "sc1a1sub2", label = "Select which cells to show", 
                                 choices = sub, selected = NULL, inline = TRUE) 
    }) 
    observeEvent(input$sc1a1sub1all, { 
        sub = strsplit(sc1conf[UI == input$sc1a1sub1]$fID, "\\|")[[1]] 
        updateCheckboxGroupInput(session, inputId = "sc1a1sub2", label = "Select which cells to show", 
                                 choices = sub, selected = sub, inline = TRUE) 
    }) 
    output$sc1a1oup1 <- renderPlot({ 
        scDRcell(sc1conf, sc1meta, input$sc1a1drX, input$sc1a1drY, input$sc1a1inp1,  
                 input$sc1a1sub1, input$sc1a1sub2, 
                 input$sc1a1siz, input$sc1a1col1, input$sc1a1ord1, 
                 input$sc1a1fsz, input$sc1a1asp, input$sc1a1txt, input$sc1a1lab1) 
    }) 
    output$sc1a1oup1.ui <- renderUI({ 
        plotOutput("sc1a1oup1", height = pList[input$sc1a1psz]) 
    }) 
    output$sc1a1oup1.pdf <- downloadHandler( 
        filename = function() { paste0("sc1",input$sc1a1drX,"_",input$sc1a1drY,"_",  
                                       input$sc1a1inp1,".pdf") }, 
        content = function(file) { ggsave( 
            file, device = "pdf", height = input$sc1a1oup1.h, width = input$sc1a1oup1.w, useDingbats = FALSE, 
            plot = scDRcell(sc1conf, sc1meta, input$sc1a1drX, input$sc1a1drY, input$sc1a1inp1,   
                            input$sc1a1sub1, input$sc1a1sub2, 
                            input$sc1a1siz, input$sc1a1col1, input$sc1a1ord1,  
                            input$sc1a1fsz, input$sc1a1asp, input$sc1a1txt, input$sc1a1lab1) ) 
        }) 
    output$sc1a1oup1.png <- downloadHandler( 
        filename = function() { paste0("sc1",input$sc1a1drX,"_",input$sc1a1drY,"_",  
                                       input$sc1a1inp1,".png") }, 
        content = function(file) { ggsave( 
            file, device = "png", height = input$sc1a1oup1.h, width = input$sc1a1oup1.w, 
            plot = scDRcell(sc1conf, sc1meta, input$sc1a1drX, input$sc1a1drY, input$sc1a1inp1,   
                            input$sc1a1sub1, input$sc1a1sub2, 
                            input$sc1a1siz, input$sc1a1col1, input$sc1a1ord1,  
                            input$sc1a1fsz, input$sc1a1asp, input$sc1a1txt, input$sc1a1lab1) ) 
        }) 
    output$sc1a1.dt <- renderDataTable({ 
        ggData = scDRnum(sc1conf, sc1meta, input$sc1a1inp1, input$sc1a1inp2, 
                         input$sc1a1sub1, input$sc1a1sub2, 
                         "sc1gexpr.h5", sc1gene, input$sc1a1splt) 
        datatable(ggData, rownames = FALSE, extensions = "Buttons", 
                  options = list(pageLength = -1, dom = "tB", buttons = c("copy", "csv", "excel"))) %>% 
            formatRound(columns = c("pctExpress"), digits = 2) 
    }) 
    
    output$sc1a1oup2 <- renderPlot({ 
        scDRgene(sc1conf, sc1meta, input$sc1a1drX, input$sc1a1drY, input$sc1a1inp2,  
                 input$sc1a1sub1, input$sc1a1sub2, 
                 "sc1gexpr.h5", sc1gene, 
                 input$sc1a1siz, input$sc1a1col2, input$sc1a1ord2, 
                 input$sc1a1fsz, input$sc1a1asp, input$sc1a1txt) 
    }) 
    output$sc1a1oup2.ui <- renderUI({ 
        plotOutput("sc1a1oup2", height = pList[input$sc1a1psz]) 
    }) 
    output$sc1a1oup2.pdf <- downloadHandler( 
        filename = function() { paste0("sc1",input$sc1a1drX,"_",input$sc1a1drY,"_",  
                                       input$sc1a1inp2,".pdf") }, 
        content = function(file) { ggsave( 
            file, device = "pdf", height = input$sc1a1oup2.h, width = input$sc1a1oup2.w, useDingbats = FALSE, 
            plot = scDRgene(sc1conf, sc1meta, input$sc1a1drX, input$sc1a1drY, input$sc1a1inp2,  
                            input$sc1a1sub1, input$sc1a1sub2, 
                            "sc1gexpr.h5", sc1gene, 
                            input$sc1a1siz, input$sc1a1col2, input$sc1a1ord2, 
                            input$sc1a1fsz, input$sc1a1asp, input$sc1a1txt) ) 
        }) 
    output$sc1a1oup2.png <- downloadHandler( 
        filename = function() { paste0("sc1",input$sc1a1drX,"_",input$sc1a1drY,"_",  
                                       input$sc1a1inp2,".png") }, 
        content = function(file) { ggsave( 
            file, device = "png", height = input$sc1a1oup2.h, width = input$sc1a1oup2.w, 
            plot = scDRgene(sc1conf, sc1meta, input$sc1a1drX, input$sc1a1drY, input$sc1a1inp2,  
                            input$sc1a1sub1, input$sc1a1sub2, 
                            "sc1gexpr.h5", sc1gene, 
                            input$sc1a1siz, input$sc1a1col2, input$sc1a1ord2, 
                            input$sc1a1fsz, input$sc1a1asp, input$sc1a1txt) ) 
        }) 
    
    
    ### Plots for tab a2 
    output$sc1a2sub1.ui <- renderUI({ 
        sub = strsplit(sc1conf[UI == input$sc1a2sub1]$fID, "\\|")[[1]] 
        checkboxGroupInput("sc1a2sub2", "Select which cells to show", inline = TRUE, 
                           choices = sub, selected = sub) 
    }) 
    observeEvent(input$sc1a2sub1non, { 
        sub = strsplit(sc1conf[UI == input$sc1a2sub1]$fID, "\\|")[[1]] 
        updateCheckboxGroupInput(session, inputId = "sc1a2sub2", label = "Select which cells to show", 
                                 choices = sub, selected = NULL, inline = TRUE) 
    }) 
    observeEvent(input$sc1a2sub1all, { 
        sub = strsplit(sc1conf[UI == input$sc1a2sub1]$fID, "\\|")[[1]] 
        updateCheckboxGroupInput(session, inputId = "sc1a2sub2", label = "Select which cells to show", 
                                 choices = sub, selected = sub, inline = TRUE) 
    }) 
    output$sc1a2oup1 <- renderPlot({ 
        scDRcell(sc1conf, sc1meta, input$sc1a2drX, input$sc1a2drY, input$sc1a2inp1,  
                 input$sc1a2sub1, input$sc1a2sub2, 
                 input$sc1a2siz, input$sc1a2col1, input$sc1a2ord1, 
                 input$sc1a2fsz, input$sc1a2asp, input$sc1a2txt, input$sc1a2lab1) 
    }) 
    output$sc1a2oup1.ui <- renderUI({ 
        plotOutput("sc1a2oup1", height = pList[input$sc1a2psz]) 
    }) 
    output$sc1a2oup1.pdf <- downloadHandler( 
        filename = function() { paste0("sc1",input$sc1a2drX,"_",input$sc1a2drY,"_",  
                                       input$sc1a2inp1,".pdf") }, 
        content = function(file) { ggsave( 
            file, device = "pdf", height = input$sc1a2oup1.h, width = input$sc1a2oup1.w, useDingbats = FALSE, 
            plot = scDRcell(sc1conf, sc1meta, input$sc1a2drX, input$sc1a2drY, input$sc1a2inp1,   
                            input$sc1a2sub1, input$sc1a2sub2, 
                            input$sc1a2siz, input$sc1a2col1, input$sc1a2ord1,  
                            input$sc1a2fsz, input$sc1a2asp, input$sc1a2txt, input$sc1a2lab1) ) 
        }) 
    output$sc1a2oup1.png <- downloadHandler( 
        filename = function() { paste0("sc1",input$sc1a2drX,"_",input$sc1a2drY,"_",  
                                       input$sc1a2inp1,".png") }, 
        content = function(file) { ggsave( 
            file, device = "png", height = input$sc1a2oup1.h, width = input$sc1a2oup1.w, 
            plot = scDRcell(sc1conf, sc1meta, input$sc1a2drX, input$sc1a2drY, input$sc1a2inp1,   
                            input$sc1a2sub1, input$sc1a2sub2, 
                            input$sc1a2siz, input$sc1a2col1, input$sc1a2ord1,  
                            input$sc1a2fsz, input$sc1a2asp, input$sc1a2txt, input$sc1a2lab1) ) 
        }) 
    
    output$sc1a2oup2 <- renderPlot({ 
        scDRcell(sc1conf, sc1meta, input$sc1a2drX, input$sc1a2drY, input$sc1a2inp2,  
                 input$sc1a2sub1, input$sc1a2sub2, 
                 input$sc1a2siz, input$sc1a2col2, input$sc1a2ord2, 
                 input$sc1a2fsz, input$sc1a2asp, input$sc1a2txt, input$sc1a2lab2) 
    }) 
    output$sc1a2oup2.ui <- renderUI({ 
        plotOutput("sc1a2oup2", height = pList[input$sc1a2psz]) 
    }) 
    output$sc1a2oup2.pdf <- downloadHandler( 
        filename = function() { paste0("sc1",input$sc1a2drX,"_",input$sc1a2drY,"_",  
                                       input$sc1a2inp2,".pdf") }, 
        content = function(file) { ggsave( 
            file, device = "pdf", height = input$sc1a2oup2.h, width = input$sc1a2oup2.w, useDingbats = FALSE, 
            plot = scDRcell(sc1conf, sc1meta, input$sc1a2drX, input$sc1a2drY, input$sc1a2inp2,   
                            input$sc1a2sub1, input$sc1a2sub2, 
                            input$sc1a2siz, input$sc1a2col2, input$sc1a2ord2,  
                            input$sc1a2fsz, input$sc1a2asp, input$sc1a2txt, input$sc1a2lab2) ) 
        }) 
    output$sc1a2oup2.png <- downloadHandler( 
        filename = function() { paste0("sc1",input$sc1a2drX,"_",input$sc1a2drY,"_",  
                                       input$sc1a2inp2,".png") }, 
        content = function(file) { ggsave( 
            file, device = "png", height = input$sc1a2oup2.h, width = input$sc1a2oup2.w, 
            plot = scDRcell(sc1conf, sc1meta, input$sc1a2drX, input$sc1a2drY, input$sc1a2inp2,   
                            input$sc1a2sub1, input$sc1a2sub2, 
                            input$sc1a2siz, input$sc1a2col2, input$sc1a2ord2,  
                            input$sc1a2fsz, input$sc1a2asp, input$sc1a2txt, input$sc1a2lab2) ) 
        }) 
    
    
    ### Plots for tab a3 
    output$sc1a3sub1.ui <- renderUI({ 
        sub = strsplit(sc1conf[UI == input$sc1a3sub1]$fID, "\\|")[[1]] 
        checkboxGroupInput("sc1a3sub2", "Select which cells to show", inline = TRUE, 
                           choices = sub, selected = sub) 
    }) 
    observeEvent(input$sc1a3sub1non, { 
        sub = strsplit(sc1conf[UI == input$sc1a3sub1]$fID, "\\|")[[1]] 
        updateCheckboxGroupInput(session, inputId = "sc1a3sub2", label = "Select which cells to show", 
                                 choices = sub, selected = NULL, inline = TRUE) 
    }) 
    observeEvent(input$sc1a3sub1all, { 
        sub = strsplit(sc1conf[UI == input$sc1a3sub1]$fID, "\\|")[[1]] 
        updateCheckboxGroupInput(session, inputId = "sc1a3sub2", label = "Select which cells to show", 
                                 choices = sub, selected = sub, inline = TRUE) 
    }) 
    output$sc1a3oup1 <- renderPlot({ 
        scDRgene(sc1conf, sc1meta, input$sc1a3drX, input$sc1a3drY, input$sc1a3inp1,  
                 input$sc1a3sub1, input$sc1a3sub2, 
                 "sc1gexpr.h5", sc1gene, 
                 input$sc1a3siz, input$sc1a3col1, input$sc1a3ord1, 
                 input$sc1a3fsz, input$sc1a3asp, input$sc1a3txt) 
    }) 
    output$sc1a3oup1.ui <- renderUI({ 
        plotOutput("sc1a3oup1", height = pList[input$sc1a3psz]) 
    }) 
    output$sc1a3oup1.pdf <- downloadHandler( 
        filename = function() { paste0("sc1",input$sc1a3drX,"_",input$sc1a3drY,"_",  
                                       input$sc1a3inp1,".pdf") }, 
        content = function(file) { ggsave( 
            file, device = "pdf", height = input$sc1a3oup1.h, width = input$sc1a3oup1.w, useDingbats = FALSE, 
            plot = scDRgene(sc1conf, sc1meta, input$sc1a3drX, input$sc1a3drY, input$sc1a3inp1,  
                            input$sc1a3sub1, input$sc1a3sub2, 
                            "sc1gexpr.h5", sc1gene, 
                            input$sc1a3siz, input$sc1a3col1, input$sc1a3ord1, 
                            input$sc1a3fsz, input$sc1a3asp, input$sc1a3txt) ) 
        }) 
    output$sc1a3oup1.png <- downloadHandler( 
        filename = function() { paste0("sc1",input$sc1a3drX,"_",input$sc1a3drY,"_",  
                                       input$sc1a3inp1,".png") }, 
        content = function(file) { ggsave( 
            file, device = "png", height = input$sc1a3oup1.h, width = input$sc1a3oup1.w, 
            plot = scDRgene(sc1conf, sc1meta, input$sc1a3drX, input$sc1a3drY, input$sc1a3inp1,  
                            input$sc1a3sub1, input$sc1a3sub2, 
                            "sc1gexpr.h5", sc1gene, 
                            input$sc1a3siz, input$sc1a3col1, input$sc1a3ord1, 
                            input$sc1a3fsz, input$sc1a3asp, input$sc1a3txt) ) 
        }) 
    
    output$sc1a3oup2 <- renderPlot({ 
        scDRgene(sc1conf, sc1meta, input$sc1a3drX, input$sc1a3drY, input$sc1a3inp2,  
                 input$sc1a3sub1, input$sc1a3sub2, 
                 "sc1gexpr.h5", sc1gene, 
                 input$sc1a3siz, input$sc1a3col2, input$sc1a3ord2, 
                 input$sc1a3fsz, input$sc1a3asp, input$sc1a3txt) 
    }) 
    output$sc1a3oup2.ui <- renderUI({ 
        plotOutput("sc1a3oup2", height = pList[input$sc1a3psz]) 
    }) 
    output$sc1a3oup2.pdf <- downloadHandler( 
        filename = function() { paste0("sc1",input$sc1a3drX,"_",input$sc1a3drY,"_",  
                                       input$sc1a3inp2,".pdf") }, 
        content = function(file) { ggsave( 
            file, device = "pdf", height = input$sc1a3oup2.h, width = input$sc1a3oup2.w, useDingbats = FALSE, 
            plot = scDRgene(sc1conf, sc1meta, input$sc1a3drX, input$sc1a3drY, input$sc1a3inp2,  
                            input$sc1a3sub1, input$sc1a3sub2, 
                            "sc1gexpr.h5", sc1gene, 
                            input$sc1a3siz, input$sc1a3col2, input$sc1a3ord2, 
                            input$sc1a3fsz, input$sc1a3asp, input$sc1a3txt) ) 
        }) 
    output$sc1a3oup2.png <- downloadHandler( 
        filename = function() { paste0("sc1",input$sc1a3drX,"_",input$sc1a3drY,"_",  
                                       input$sc1a3inp2,".png") }, 
        content = function(file) { ggsave( 
            file, device = "png", height = input$sc1a3oup2.h, width = input$sc1a3oup2.w, 
            plot = scDRgene(sc1conf, sc1meta, input$sc1a3drX, input$sc1a3drY, input$sc1a3inp2,  
                            input$sc1a3sub1, input$sc1a3sub2, 
                            "sc1gexpr.h5", sc1gene, 
                            input$sc1a3siz, input$sc1a3col2, input$sc1a3ord2, 
                            input$sc1a3fsz, input$sc1a3asp, input$sc1a3txt) ) 
        }) 
    
    
    ### Plots for tab b2 
    output$sc1b2sub1.ui <- renderUI({ 
        sub = strsplit(sc1conf[UI == input$sc1b2sub1]$fID, "\\|")[[1]] 
        checkboxGroupInput("sc1b2sub2", "Select which cells to show", inline = TRUE, 
                           choices = sub, selected = sub) 
    }) 
    observeEvent(input$sc1b2sub1non, { 
        sub = strsplit(sc1conf[UI == input$sc1b2sub1]$fID, "\\|")[[1]] 
        updateCheckboxGroupInput(session, inputId = "sc1b2sub2", label = "Select which cells to show", 
                                 choices = sub, selected = NULL, inline = TRUE) 
    }) 
    observeEvent(input$sc1b2sub1all, { 
        sub = strsplit(sc1conf[UI == input$sc1b2sub1]$fID, "\\|")[[1]] 
        updateCheckboxGroupInput(session, inputId = "sc1b2sub2", label = "Select which cells to show", 
                                 choices = sub, selected = sub, inline = TRUE) 
    }) 
    output$sc1b2oup1 <- renderPlot({ 
        scDRcoex(sc1conf, sc1meta, input$sc1b2drX, input$sc1b2drY,   
                 input$sc1b2inp1, input$sc1b2inp2, input$sc1b2sub1, input$sc1b2sub2, 
                 "sc1gexpr.h5", sc1gene, 
                 input$sc1b2siz, input$sc1b2col1, input$sc1b2ord1, 
                 input$sc1b2fsz, input$sc1b2asp, input$sc1b2txt) 
    }) 
    output$sc1b2oup1.ui <- renderUI({ 
        plotOutput("sc1b2oup1", height = pList2[input$sc1b2psz]) 
    }) 
    output$sc1b2oup1.pdf <- downloadHandler( 
        filename = function() { paste0("sc1",input$sc1b2drX,"_",input$sc1b2drY,"_",  
                                       input$sc1b2inp1,"_",input$sc1b2inp2,".pdf") }, 
        content = function(file) { ggsave( 
            file, device = "pdf", height = input$sc1b2oup1.h, width = input$sc1b2oup1.w, useDingbats = FALSE, 
            plot = scDRcoex(sc1conf, sc1meta, input$sc1b2drX, input$sc1b2drY,  
                            input$sc1b2inp1, input$sc1b2inp2, input$sc1b2sub1, input$sc1b2sub2, 
                            "sc1gexpr.h5", sc1gene, 
                            input$sc1b2siz, input$sc1b2col1, input$sc1b2ord1, 
                            input$sc1b2fsz, input$sc1b2asp, input$sc1b2txt) ) 
        }) 
    output$sc1b2oup1.png <- downloadHandler( 
        filename = function() { paste0("sc1",input$sc1b2drX,"_",input$sc1b2drY,"_",  
                                       input$sc1b2inp1,"_",input$sc1b2inp2,".png") }, 
        content = function(file) { ggsave( 
            file, device = "png", height = input$sc1b2oup1.h, width = input$sc1b2oup1.w, 
            plot = scDRcoex(sc1conf, sc1meta, input$sc1b2drX, input$sc1b2drY,  
                            input$sc1b2inp1, input$sc1b2inp2, input$sc1b2sub1, input$sc1b2sub2, 
                            "sc1gexpr.h5", sc1gene, 
                            input$sc1b2siz, input$sc1b2col1, input$sc1b2ord1, 
                            input$sc1b2fsz, input$sc1b2asp, input$sc1b2txt) ) 
        }) 
    output$sc1b2oup2 <- renderPlot({ 
        scDRcoexLeg(input$sc1b2inp1, input$sc1b2inp2, input$sc1b2col1, input$sc1b2fsz) 
    }) 
    output$sc1b2oup2.ui <- renderUI({ 
        plotOutput("sc1b2oup2", height = "300px") 
    }) 
    output$sc1b2oup2.pdf <- downloadHandler( 
        filename = function() { paste0("sc1",input$sc1b2drX,"_",input$sc1b2drY,"_",  
                                       input$sc1b2inp1,"_",input$sc1b2inp2,"_leg.pdf") }, 
        content = function(file) { ggsave( 
            file, device = "pdf", height = 3, width = 4, useDingbats = FALSE, 
            plot = scDRcoexLeg(input$sc1b2inp1, input$sc1b2inp2, input$sc1b2col1, input$sc1b2fsz) ) 
        }) 
    output$sc1b2oup2.png <- downloadHandler( 
        filename = function() { paste0("sc1",input$sc1b2drX,"_",input$sc1b2drY,"_",  
                                       input$sc1b2inp1,"_",input$sc1b2inp2,"_leg.png") }, 
        content = function(file) { ggsave( 
            file, device = "png", height = 3, width = 4, 
            plot = scDRcoexLeg(input$sc1b2inp1, input$sc1b2inp2, input$sc1b2col1, input$sc1b2fsz) ) 
        }) 
    output$sc1b2.dt <- renderDataTable({ 
        ggData = scDRcoexNum(sc1conf, sc1meta, input$sc1b2inp1, input$sc1b2inp2, 
                             input$sc1b2sub1, input$sc1b2sub2, "sc1gexpr.h5", sc1gene) 
        datatable(ggData, rownames = FALSE, extensions = "Buttons", 
                  options = list(pageLength = -1, dom = "tB", buttons = c("copy", "csv", "excel"))) %>% 
            formatRound(columns = c("percent"), digits = 2) 
    }) 
    
    
    ### Plots for tab c1 
    output$sc1c1sub1.ui <- renderUI({ 
        sub = strsplit(sc1conf[UI == input$sc1c1sub1]$fID, "\\|")[[1]] 
        checkboxGroupInput("sc1c1sub2", "Select which cells to show", inline = TRUE, 
                           choices = sub, selected = sub) 
    }) 
    observeEvent(input$sc1c1sub1non, { 
        sub = strsplit(sc1conf[UI == input$sc1c1sub1]$fID, "\\|")[[1]] 
        updateCheckboxGroupInput(session, inputId = "sc1c1sub2", label = "Select which cells to show", 
                                 choices = sub, selected = NULL, inline = TRUE) 
    }) 
    observeEvent(input$sc1c1sub1all, { 
        sub = strsplit(sc1conf[UI == input$sc1c1sub1]$fID, "\\|")[[1]] 
        updateCheckboxGroupInput(session, inputId = "sc1c1sub2", label = "Select which cells to show", 
                                 choices = sub, selected = sub, inline = TRUE) 
    }) 
    output$sc1c1oup <- renderPlot({ 
        scVioBox(sc1conf, sc1meta, input$sc1c1inp1, input$sc1c1inp2, 
                 input$sc1c1sub1, input$sc1c1sub2, 
                 "sc1gexpr.h5", sc1gene, input$sc1c1typ, input$sc1c1pts, 
                 input$sc1c1siz, input$sc1c1fsz) 
    }) 
    output$sc1c1oup.ui <- renderUI({ 
        plotOutput("sc1c1oup", height = pList2[input$sc1c1psz]) 
    }) 
    output$sc1c1oup.pdf <- downloadHandler( 
        filename = function() { paste0("sc1",input$sc1c1typ,"_",input$sc1c1inp1,"_",  
                                       input$sc1c1inp2,".pdf") }, 
        content = function(file) { ggsave( 
            file, device = "pdf", height = input$sc1c1oup.h, width = input$sc1c1oup.w, useDingbats = FALSE, 
            plot = scVioBox(sc1conf, sc1meta, input$sc1c1inp1, input$sc1c1inp2, 
                            input$sc1c1sub1, input$sc1c1sub2, 
                            "sc1gexpr.h5", sc1gene, input$sc1c1typ, input$sc1c1pts, 
                            input$sc1c1siz, input$sc1c1fsz) ) 
        }) 
    output$sc1c1oup.png <- downloadHandler( 
        filename = function() { paste0("sc1",input$sc1c1typ,"_",input$sc1c1inp1,"_",  
                                       input$sc1c1inp2,".png") }, 
        content = function(file) { ggsave( 
            file, device = "png", height = input$sc1c1oup.h, width = input$sc1c1oup.w, 
            plot = scVioBox(sc1conf, sc1meta, input$sc1c1inp1, input$sc1c1inp2, 
                            input$sc1c1sub1, input$sc1c1sub2, 
                            "sc1gexpr.h5", sc1gene, input$sc1c1typ, input$sc1c1pts, 
                            input$sc1c1siz, input$sc1c1fsz) ) 
        }) 
    
    
    ### Plots for tab c2 
    output$sc1c2sub1.ui <- renderUI({ 
        sub = strsplit(sc1conf[UI == input$sc1c2sub1]$fID, "\\|")[[1]] 
        checkboxGroupInput("sc1c2sub2", "Select which cells to show", inline = TRUE, 
                           choices = sub, selected = sub) 
    }) 
    observeEvent(input$sc1c2sub1non, { 
        sub = strsplit(sc1conf[UI == input$sc1c2sub1]$fID, "\\|")[[1]] 
        updateCheckboxGroupInput(session, inputId = "sc1c2sub2", label = "Select which cells to show", 
                                 choices = sub, selected = NULL, inline = TRUE) 
    }) 
    observeEvent(input$sc1c2sub1all, { 
        sub = strsplit(sc1conf[UI == input$sc1c2sub1]$fID, "\\|")[[1]] 
        updateCheckboxGroupInput(session, inputId = "sc1c2sub2", label = "Select which cells to show", 
                                 choices = sub, selected = sub, inline = TRUE) 
    }) 
    output$sc1c2oup <- renderPlot({ 
        scProp(sc1conf, sc1meta, input$sc1c2inp1, input$sc1c2inp2,  
               input$sc1c2sub1, input$sc1c2sub2, 
               input$sc1c2typ, input$sc1c2flp, input$sc1c2fsz) 
    }) 
    output$sc1c2oup.ui <- renderUI({ 
        plotOutput("sc1c2oup", height = pList2[input$sc1c2psz]) 
    }) 
    output$sc1c2oup.pdf <- downloadHandler( 
        filename = function() { paste0("sc1",input$sc1c2typ,"_",input$sc1c2inp1,"_",  
                                       input$sc1c2inp2,".pdf") }, 
        content = function(file) { ggsave( 
            file, device = "pdf", height = input$sc1c2oup.h, width = input$sc1c2oup.w, useDingbats = FALSE, 
            plot = scProp(sc1conf, sc1meta, input$sc1c2inp1, input$sc1c2inp2,  
                          input$sc1c2sub1, input$sc1c2sub2, 
                          input$sc1c2typ, input$sc1c2flp, input$sc1c2fsz) ) 
        }) 
    output$sc1c2oup.png <- downloadHandler( 
        filename = function() { paste0("sc1",input$sc1c2typ,"_",input$sc1c2inp1,"_",  
                                       input$sc1c2inp2,".png") }, 
        content = function(file) { ggsave( 
            file, device = "png", height = input$sc1c2oup.h, width = input$sc1c2oup.w, 
            plot = scProp(sc1conf, sc1meta, input$sc1c2inp1, input$sc1c2inp2,  
                          input$sc1c2sub1, input$sc1c2sub2, 
                          input$sc1c2typ, input$sc1c2flp, input$sc1c2fsz) ) 
        }) 
    
    
    ### Plots for tab d1 
    output$sc1d1sub1.ui <- renderUI({ 
        sub = strsplit(sc1conf[UI == input$sc1d1sub1]$fID, "\\|")[[1]] 
        checkboxGroupInput("sc1d1sub2", "Select which cells to show", inline = TRUE, 
                           choices = sub, selected = sub) 
    }) 
    observeEvent(input$sc1d1sub1non, { 
        sub = strsplit(sc1conf[UI == input$sc1d1sub1]$fID, "\\|")[[1]] 
        updateCheckboxGroupInput(session, inputId = "sc1d1sub2", label = "Select which cells to show", 
                                 choices = sub, selected = NULL, inline = TRUE) 
    }) 
    observeEvent(input$sc1d1sub1all, { 
        sub = strsplit(sc1conf[UI == input$sc1d1sub1]$fID, "\\|")[[1]] 
        updateCheckboxGroupInput(session, inputId = "sc1d1sub2", label = "Select which cells to show", 
                                 choices = sub, selected = sub, inline = TRUE) 
    }) 
    output$sc1d1oupTxt <- renderUI({ 
        geneList = scGeneList(input$sc1d1inp, sc1gene) 
        if(nrow(geneList) > 50){ 
            HTML("More than 50 input genes! Please reduce the gene list!") 
        } else { 
            oup = paste0(nrow(geneList[present == TRUE]), " genes OK and will be plotted") 
            if(nrow(geneList[present == FALSE]) > 0){ 
                oup = paste0(oup, "<br/>", 
                             nrow(geneList[present == FALSE]), " genes not found (", 
                             paste0(geneList[present == FALSE]$gene, collapse = ", "), ")") 
            } 
            HTML(oup) 
        } 
    }) 
    output$sc1d1oup <- renderPlot({ 
        scBubbHeat(sc1conf, sc1meta, input$sc1d1inp, input$sc1d1grp, input$sc1d1plt, 
                   input$sc1d1sub1, input$sc1d1sub2, "sc1gexpr.h5", sc1gene, 
                   input$sc1d1scl, input$sc1d1row, input$sc1d1col, 
                   input$sc1d1cols, input$sc1d1fsz) 
    }) 
    output$sc1d1oup.ui <- renderUI({ 
        plotOutput("sc1d1oup", height = pList3[input$sc1d1psz]) 
    }) 
    output$sc1d1oup.pdf <- downloadHandler( 
        filename = function() { paste0("sc1",input$sc1d1plt,"_",input$sc1d1grp,".pdf") }, 
        content = function(file) { ggsave( 
            file, device = "pdf", height = input$sc1d1oup.h, width = input$sc1d1oup.w, 
            plot = scBubbHeat(sc1conf, sc1meta, input$sc1d1inp, input$sc1d1grp, input$sc1d1plt, 
                              input$sc1d1sub1, input$sc1d1sub2, "sc1gexpr.h5", sc1gene, 
                              input$sc1d1scl, input$sc1d1row, input$sc1d1col, 
                              input$sc1d1cols, input$sc1d1fsz, save = TRUE) ) 
        }) 
    output$sc1d1oup.png <- downloadHandler( 
        filename = function() { paste0("sc1",input$sc1d1plt,"_",input$sc1d1grp,".png") }, 
        content = function(file) { ggsave( 
            file, device = "png", height = input$sc1d1oup.h, width = input$sc1d1oup.w, 
            plot = scBubbHeat(sc1conf, sc1meta, input$sc1d1inp, input$sc1d1grp, input$sc1d1plt, 
                              input$sc1d1sub1, input$sc1d1sub2, "sc1gexpr.h5", sc1gene, 
                              input$sc1d1scl, input$sc1d1row, input$sc1d1col, 
                              input$sc1d1cols, input$sc1d1fsz, save = TRUE) ) 
        })    
    

}
#####

#####

shinyApp(ui, server)


