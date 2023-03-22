



yy<-function(){
  t1<-1
  
  t2<-2
}






plot.microNetProps <- function(x,
                               layout = "spring",
                               sameLayout = FALSE,
                               layoutGroup = NULL,
                               repulsion = 1,
                               groupNames = NULL,
                               groupsChanged = FALSE,
                               labels = NULL,
                               shortenLabels = "intelligent",
                               labelLength = 6L,
                               labelPattern = c(5,"'",3),
                               charToRm = NULL,
                               labelScale = TRUE,
                               labelFont = 1,
                               #nodes
                               nodeFilter = "none",
                               nodeFilterPar = NULL,
                               rmSingles = "none",
                               nodeSize = "fix",
                               nodeSizeSpread = 3,
                               nodeColor = "cluster",
                               colorVec = NULL,
                               featVecCol = NULL,
                               sameClustCol = TRUE,
                               sameColThresh = 2L,
                               nodeShape = NULL,
                               featVecShape = NULL,
                               nodeTransp = 60,
                               borderWidth = 1,
                               borderCol = "gray80",
                               #hubs
                               highlightHubs = TRUE,
                               hubTransp = NULL,
                               hubLabelFont = NULL,
                               hubBorderWidth = NULL,
                               hubBorderCol = "black",
                               #edges
                               edgeFilter = "none",
                               edgeFilterPar = NULL,
                               edgeInvisFilter = "none",
                               edgeInvisPar = NULL,
                               edgeWidth = 1,
                               colorNegAsso = TRUE,
                               posCol = NULL,
                               negCol = NULL,
                               cut = NULL,
                               edgeTranspLow = 0,
                               edgeTranspHigh = 0,
                               cexNodes = 1,
                               cexHubs = 1.2,
                               cexLabels = 1,
                               cexTitle = 1.2,
                               showTitle = NULL,
                               title1 = NULL,
                               title2 = NULL,
                               mar = c(1, 3, 3, 3),
                               ...){

  inputArgs <- c(as.list(environment()), list(...))

  outputArgs <- except_plot_networks(inputArgs)
  for(i in 1:length(outputArgs)){
    assign(names(outputArgs)[i], outputArgs[[i]])
  }

  xgroups <- x$input$groups
  isempty1 <- x$isempty$isempty1
  isempty2 <- x$isempty$isempty2

  # if twoNets==TRUE, two networks are plotted
  twoNets <- x$input$twoNets

  if(twoNets & !is.null(xgroups)){
    deletespace <- function(name){
      if(grepl(" ", strsplit(name, "")[[1]])[1]){
        gname <- strsplit(name, "")[[1]]
        name <- paste(gname[-1], collapse = "")
      }
      name
    }
    xgroups[1] <- deletespace(xgroups[1])
    xgroups[2] <- deletespace(xgroups[2])
  }

  #=============================================================================
  opar <- par()

  adja1_orig <- x$input$adjaMat1

  adja1 <- filter_edges(adja = adja1_orig, edgeFilter = edgeFilter,
                            edgeFilterPar = edgeFilterPar)

  adja1List <- filter_nodes(adja = adja1, nodeFilter = nodeFilter,
                           nodeFilterPar = nodeFilterPar, layout = layout,
                           degree = x$centralities$degree1,
                           between = x$centralities$between1,
                           close = x$centralities$close1,
                           eigen = x$centralities$eigenv1,
                           cluster = x$clustering$clust1)
  adja1 <- adja1List$adja
  keep1 <- adja1List$keep

  if(twoNets){
    adja2_orig <- x$input$adjaMat2

    adja2 <- filter_edges(adja = adja2_orig, edgeFilter = edgeFilter,
                         edgeFilterPar = edgeFilterPar)

    adja2List <- filter_nodes(adja = adja2_orig, nodeFilter = nodeFilter,
                             nodeFilterPar = nodeFilterPar, layout = layout,
                             degree = x$centralities$degree2,
                             between = x$centralities$between2,
                             close = x$centralities$close2,
                             eigen = x$centralities$eigenv2,
                             cluster = x$clustering$clust2)
    adja2 <- adja2List$adja
    keep2 <- adja2List$keep

    if(length(keep1) == 0 & length(keep2) == 0){
      stop("No nodes remaining in both networks after node filtering.")
    }

    if(sameLayout){
      keep.tmp <- union(keep1, keep2)
      keep <- colnames(adja1)[colnames(adja1) %in% keep.tmp]
      adja1 <- adja1[keep, keep]
      adja2 <- adja2[keep, keep]
    } else{
      adja1 <- adja1[keep1, keep1]
      adja2 <- adja2[keep2, keep2]
    }

  } else{
    if(length(keep1) == 0){
      stop("No nodes remaining after node filtering.")
    }
    adja1 <- adja1[keep1, keep1]
    adja2 <- NULL
    adja2_orig <- NULL
  }


  if(rmSingles != "none"){
    adja1.tmp <- adja1
    diag(adja1.tmp) <- 0
    zeros1 <- sapply(1:nrow(adja1.tmp), function(i){ all(adja1.tmp[i, ] == 0) })
    names(zeros1) <- colnames(adja1.tmp)

    if(any(zeros1 == TRUE)){
      torm1 <- which(zeros1 == TRUE)
    } else{
      torm1 <- NULL
    }

    if(twoNets){
      adja2.tmp <- adja2
      diag(adja2.tmp) <- 0
      zeros2 <- sapply(1:nrow(adja2.tmp), function(i){ all(adja2.tmp[i, ] == 0) })
      names(zeros2) <- colnames(adja2.tmp)

      if(any(zeros1 == TRUE)){
        torm2 <- which(zeros2 == TRUE)
      } else{
        torm2 <- NULL
      }

      torm_both <- intersect(torm1, torm2)

      if(rmSingles == "inboth" & length(torm_both) != 0){
        adja1 <- adja1.tmp[-torm_both, -torm_both]
        adja2 <- adja2.tmp[-torm_both, -torm_both]

      } else if(rmSingles == "all"){
        if(sameLayout){
          if(length(torm_both) != 0){
            warning('Argument "rmSingles" has been set to "inboth" due to same layout.')
            adja1 <- adja1.tmp[-torm_both, -torm_both]
            adja2 <- adja2.tmp[-torm_both, -torm_both]
          }
        } else{
          if(length(torm1) != 0) adja1 <- adja1.tmp[-torm1, -torm1]
          if(length(torm2) != 0) adja2 <- adja2.tmp[-torm2, -torm2]

        }
      }

      isempty1 <- all(adja1[upper.tri(adja1)] == 0)
      isempty2 <- all(adja2[upper.tri(adja2)] == 0)

    } else{
      if(length(torm1) != 0) adja1 <- adja1.tmp[-torm1, -torm1]
      isempty1 <- all(adja1[upper.tri(adja1)] == 0)
      isempty2 <- NULL
      if(isempty1) stop("Network is empty and cannot be plotted.")
    }
  }

  kept1 <- which(colnames(adja1_orig) %in% colnames(adja1))
  kept2 <- which(colnames(adja2_orig) %in% colnames(adja2))

  #==========================================================================
  # node colors

  if(nodeColor == "cluster"){
    clust1 <- x$clustering$clust1
    clust2 <- x$clustering$clust2

    if(!(is.null(clust1) & is.null(clust2))){
      clustcolors <- get_clust_cols(clust1 = clust1, clust2 = clust2,
                                    adja1 = adja1, adja2 = adja2,
                                    kept1 = kept1, kept2 = kept2,
                                    isempty1 = isempty1, isempty2 = isempty2,
                                    colorVec = colorVec,
                                    sameClustCol = sameClustCol,
                                    sameColThresh = sameColThresh,
                                    twoNets = twoNets)
      nodecol1 <- clustcolors$clustcol1
      nodecol2 <- clustcolors$clustcol2

    } else{
      warning('No clusterings returned from "netAnalyze".')
      nodecol1 <- rep("grey40", ncol(adja1))
      if(twoNets){
        nodecol2 <- rep("grey40", ncol(adja2))
      } else{
        nodecol2 <- NULL
      }
    }

  } else if(nodeColor == "feature"){
    featVecCol <- as.factor(featVecCol)
    stopifnot(all(colnames(adja1) %in% names(featVecCol)))

    if(is.null(colorVec)){
      cols <- rainbow(length(levels(featVecCol)))
    } else{
      if(length(colorVec) < length(levels(featVecCol))){
        stop(paste("Argument 'featVecCol' has", length(levels(featVecCol)),
                   "levels but 'colorVec' has only length ", length(colorVec)))
      }
      cols <- colorVec
    }

    feature1 <- featVecCol[colnames(adja1)]
    nodecol1 <- character(length(feature1))
    for(i in seq_along(levels(feature1))){
      nodecol1[feature1 == levels(feature1)[i]] <- cols[i]
    }
    nodecol2 <- NULL

    if(twoNets){
      stopifnot(all(colnames(adja2) %in% names(featVecCol)))

      feature2 <- featVecCol[colnames(adja2)]
      nodecol2 <- character(length(feature2))
      for(i in seq_along(levels(feature2))){
        nodecol2[feature2 == levels(feature2)[i]] <- cols[i]
      }
    }

  } else if(nodeColor == "colorVec"){
    stopifnot(all(colnames(adja1) %in% names(colorVec)))
    nodecol1 <- colorVec[colnames(adja1)]
    nodecol2 <- NULL

    if(twoNets){
      stopifnot(all(colnames(adja2) %in% names(colorVec)))
      nodecol2 <- colorVec[colnames(adja2)]
    }

  } else if(is.character(nodeColor)){
    nodecol1 <- rep(nodeColor, ncol(adja1))
    nodecol2 <- NULL
    if(twoNets){
      nodecol2 <- rep(nodeColor, ncol(adja2))
    }

  } else{
    nodecol1 <- rep("grey40", ncol(adja1))
    if(twoNets){
      nodecol2 <- rep("grey40", ncol(adja2))
    }

  }

  if(nodeTransp > 0){
    nodecol1 <- col_to_transp(nodecol1, nodeTransp)
    if(twoNets){
      nodecol2 <- col_to_transp(nodecol2, nodeTransp)
    }
  }

  #==========================================================================
  # define node shapes


  if(is.null(featVecShape)){
    if(is.null(nodeShape)){
      nodeShape1 <- nodeShape2 <- "circle"
    } else{
      nodeShape2 <- nodeShape2 <- nodeShape[1]
    }

  } else{
    stopifnot(all(colnames(adja1) %in% names(featVecShape)))

    shape_fact <- as.factor(featVecShape)
    nlevel <- length(levels(shape_fact))

    if(nlevel > 4){
      stop("Argument 'featVecShape' may at most have four factor levels.")
    }

    if(is.null(nodeShape)){
      nodeShape <- c("circle", "square", "triangle", "diamond")
    } else{
      if(length(nodeShape) < nlevel){
        stop(paste("Argument 'nodeShape' must have length", nlevel,
                   "because 'featVecShape' has", nlevel, "levels."))
      }
    }

    shapeVec <- as.character(featVecShape)
    names(shapeVec) <- names(featVecShape)

    for(i in 1:nlevel){
      shapeVec[featVecShape == levels(shape_fact)[i]] <- nodeShape[i]
    }

    nodeShape1 <- shapeVec[colnames(adja1)]

    if(twoNets){
      stopifnot(all(colnames(adja2) %in% names(featVecShape)))
      nodeShape2 <- shapeVec[colnames(adja2)]
    }
  }



  #==========================================================================
  # define size of vertices

  if(!is.numeric(nodeSize)){
    nodeSize1 <- get_node_size(nodeSize = nodeSize, nodeSizeSpread = nodeSizeSpread,
                          adja = adja1, countMat = x$input$countMat1,
                          normCounts = x$input$normCounts1, kept = kept1,
                          cexNodes = cexNodes, cexHubs = cexHubs,
                          hubs = x$hubs$hubs1, highlightHubs = highlightHubs,
                          degree = x$centralities$degree1,
                          between = x$centralities$between1,
                          close = x$centralities$close1,
                          eigen = x$centralities$eigenv1)
    if(twoNets){
      nodeSize2 <- get_node_size(nodeSize = nodeSize, nodeSizeSpread = nodeSizeSpread,
                            adja = adja2, countMat = x$input$countMat2,
                            normCounts = x$input$normCounts2, kept = kept2,
                            cexNodes = cexNodes, cexHubs = cexHubs,
                            hubs = x$hubs$hubs2, highlightHubs = highlightHubs,
                            degree = x$centralities$degree2,
                            between = x$centralities$between2,
                            close = x$centralities$close2,
                            eigen = x$centralities$eigenv2)
    }
  }

  #===============================================
  # define border colors and fonts

  border1 <- rep(borderCol, nrow(adja1))
  labelFont1 <- rep(labelFont, nrow(adja1))
  borderWidth1 <- rep(borderWidth, nrow(adja1))

  if(highlightHubs){
    if(is.null(hubLabelFont)){
      hubLabelFont <- labelFont * 2
    }
    if(is.null(hubBorderWidth)){
      hubBorderWidth <- borderWidth * 2
    }
    if(is.null(hubTransp)){
      hubTransp <- nodeTransp / 2
    }

    hubs1 <- x$hubs$hubs1
    hubs1 <- hubs1[hubs1 %in% colnames(adja1)]

    border1[match(hubs1, rownames(adja1))] <- hubBorderCol
    labelFont1[match(hubs1, rownames(adja1))] <- hubLabelFont
    borderWidth1[match(hubs1, rownames(adja1))] <- hubBorderWidth

    if(nodeTransp != hubTransp){
      hubcol1 <- col_to_transp(nodecol1, hubTransp)
      nodecol1[rownames(adja1) %in% hubs1] <- hubcol1[rownames(adja1) %in% hubs1]
    }

  }

  if(twoNets){
    border2 <- rep(borderCol, nrow(adja2))
    labelFont2 <- rep(labelFont, nrow(adja2))
    borderWidth2 <- rep(borderWidth, nrow(adja2))

    if(highlightHubs){
      hubs2 <- x$hubs$hubs2
      hubs2 <- hubs2[hubs2 %in% colnames(adja2)]

      border2[match(hubs2, rownames(adja2))] <- hubBorderCol
      labelFont2[match(hubs2, rownames(adja2))] <- hubLabelFont
      borderWidth2[match(hubs2, rownames(adja2))] <- hubBorderWidth

      if(nodeTransp != hubTransp){
        hubcol2 <- col_to_transp(nodecol2, hubTransp)
        nodecol2[rownames(adja2) %in% hubs2] <- hubcol2[rownames(adja2) %in% hubs2]
      }
    }
  }

  #==========================================================================

  if(showTitle){
    if(twoNets){
      if(is.null(title1) || is.null(title2)){
        if(is.null(groupNames)){
          main1 = paste0("group '" , xgroups[1], "'" )
          main2 = paste0("group '" , xgroups[2], "'" )
        } else{
          stopifnot(length(groupNames) == 2)

          main1 = as.character(groupNames[1])
          main2 = as.character(groupNames[2])
        }
      } else{
        main1 <- title1
        main2 <- title2
      }
    } else{
      if(is.null(title1)){
        main1 <- ""
        warning('Argument "showTitle" is TRUE, but no title was defined.')
      } else{
        main1 <- title1
      }
    }

  }


  #==========================================================================
  # define cut parameter for qgraph

  if(!is.null(cut) & (length(cut) == 1)){
    stopifnot(is.numeric(cut))
    cut1 <- cut2 <- cut
  } else if(!is.null(cut) & (length(cut) == 1)){
    stopifnot(is.numeric(cut))
    cut1 <- cut[1]
    cut2 <- cut[2]
  } else{
    weights1 <- adja1[upper.tri(adja1)]
    weights1 <- weights1[abs(weights1) > 0]
    nNodes1 <- ncol(adja1)

    if(twoNets){
      weights2 <- adja2[upper.tri(adja2)]
      weights2 <- weights2[abs(weights2) > 0]
      nNodes2 <- ncol(adja2)
    } else{
      weights2 <- 0
      nNodes2 <- 1
    }

    if (length(weights1) > (3 * nNodes1) || length(weights2) > (3 * nNodes2)) {
      cut1 <- max(sort(abs(weights1), decreasing = TRUE)[2 * nNodes1],
                  quantile(abs(weights1), 0.75))
      cut2 <- ifelse(twoNets, max(sort(abs(weights2), decreasing = TRUE)[2 * nNodes2],
                                  quantile(abs(weights2), 0.75)), 0)


    } else if (length(weights1) > 1 || length(weights2) > 1){
      cut1 <- quantile(abs(weights1), 0.75)
      cut2 <- ifelse(twoNets, quantile(abs(weights2), 0.75), 0)

    } else {
      cut1 <- 0
      cut2 <- 0
    }
  }

  #==========================================================================
  # define edge colors

  if(is.null(posCol)){
    poscol1 <- "#009900"
    poscol2 <- "darkgreen"
  } else{
    if(length(posCol) == 2){
      poscol1 <- posCol[1]
      poscol2 <- posCol[2]
    } else{
      poscol1 <- poscol2 <- posCol
    }

  }

  if(is.null(negCol)){
    negcol1 <- "red"
    negcol2 <- "#BF0000"
  } else if(negCol == FALSE){
    negcol1 <- poscol1
    negcol2 <- poscol2
  } else{
    if(length(negCol) == 2){
      negcol1 <- negCol[1]
      negcol2 <- negCol[2]
    } else{
      negcol1 <- negcol2 <- negCol
    }
  }

  if(edgeTranspLow > 0){  # transparency for values below cut
    poscol1 <- col_to_transp(poscol1, edgeTranspLow)
    negcol1 <- col_to_transp(negcol1, edgeTranspLow)
  }

  if(edgeTranspHigh > 0){  # transparency for values above cut
    poscol2 <- col_to_transp(poscol2, edgeTranspHigh)
    negcol2 <- col_to_transp(negcol2, edgeTranspHigh)
  }

  if(!is.null(x$input$assoEst1)){
    assoMat1 <- x$input$assoEst1
  } else{
    assoMat1 <- x$input$dissScale1
  }


  ###############################
  if(x$input$assoType == "dissimilarity"){
    assoMat1 <- x$input$dissEst1
  } else{
    assoMat1 <- x$input$assoEst1
  }

  ###############################
  colmat1 <- matrix(NA, ncol(assoMat1), ncol(assoMat1), dimnames = dimnames(assoMat1))


  if(colorNegAsso){
    colmat1[assoMat1 < 0] <- negcol1
    colmat1[assoMat1 >= 0] <- poscol1

    colmat1 <- colmat1[kept1, kept1]

    colmat1[abs(adja1) >= cut1 & colmat1 == poscol1] <- poscol2
    colmat1[abs(adja1) >= cut1 & colmat1 == negcol1] <- negcol2
  } else{
    colmat1 <- poscol1
    colmat1 <- colmat1[kept1, kept1]
    colmat1[abs(adja1) >= cut1 & colmat1 == poscol1] <- poscol2
  }

  #--------------------------------------------
  # define layout

  diag(adja1) <- 0

  if(is.matrix(layout)){
    lay1 <- layout

  } else if(is.function(layout)){
    gr <- graph_from_adjacency_matrix(adja1, mode = "undirected",
                                      weighted = TRUE, diag = FALSE)
    lay1 <- layout(gr)
    rownames(lay1) <- colnames(adja1)
  } else{
    lay1 <- qgraph(adja1, color = nodecol1, layout = layout, vsize = nodeSize1,
                   labels = colnames(adja1),label.scale = labelScale,
                   border.color = border1, border.width = borderWidth1,
                   label.font = labelFont1, label.cex = cexLabels,
                   edge.width = edgeWidth, repulsion = repulsion, cut = cut1,
                   DoNotPlot = TRUE, shape = nodeShape1, ...)$layout
    rownames(lay1) <- rownames(adja1)
  }
  lay2 <- NULL

  #--------------------------------------------
  # define labels


  if(is.null(labels)){

    adja.tmp <- rename_taxa(adja1, toRename = "both", shortenLabels = shortenLabels,
                        labelLength = labelLength, labelPattern = labelPattern,
                        charToRm = charToRm)
    labels1 <- rownames(adja.tmp)

  } else if(is.logical(labels)){
    labels1 <- labels
  } else if(is.list(labels)){
    labels1 <- labels[[1]]
  } else if(is.vector(labels)){
    labels1 <- labels
  } else{
    stop("Argument 'labels' must be either a vector with a label for each node,
         or a list of length 2 naming the nodes in each network.")
  }

  #--------------------------------------------
  # filter edges (without influencing the layout)


  if(edgeInvisFilter != "none"){
    adja1 <- filter_edges(adja1, edgeFilter = edgeInvisFilter,
                         edgeFilterPar = edgeInvisPar)
  }

  #==========================================================================
  if(twoNets){
    # define edge colors

    if(x$input$assoType == "dissimilarity"){
      assoMat2 <- x$input$dissEst2
    } else{
      assoMat2 <- x$input$assoEst2
    }

    colmat2 <- matrix(NA, ncol(assoMat2), ncol(assoMat2),
                      dimnames = dimnames(assoMat2))

    if(colorNegAsso){
      colmat2[assoMat2 < 0] <- negcol1
      colmat2[assoMat2 >= 0] <- poscol1

      colmat2 <- colmat2[kept2, kept2]

      colmat2[abs(adja2) >= cut2 & colmat2 == poscol1] <- poscol2
      colmat2[abs(adja2) >= cut2 & colmat2 == negcol1] <- negcol2
    } else{
      colmat2 <- poscol1
      colmat2 <- colmat2[kept2, kept2]
      colmat2[abs(adja2) >= cut2 & colmat2 == poscol1] <- poscol2
    }

    #--------------------------------------------
    # define layout

    diag(adja2) <- 0

    if(is.matrix(layout)){

      lay2 <- layout

    } else if(is.function(layout)){
      gr <- graph_from_adjacency_matrix(adja2, mode = "undirected",
                                        weighted = TRUE, diag = FALSE)
      lay2 <- layout(gr)
      rownames(lay2) <- colnames(adja2)

    } else{
      lay2 <- qgraph(adja2, color = nodecol2, layout = layout, vsize = nodeSize2,
                     labels = colnames(adja2), label.scale = labelScale,
                     border.color = border2, border.width = borderWidth2,
                     label.font = labelFont2, label.cex = cexLabels,
                     edge.width = edgeWidth,  repulsion = repulsion, cut = cut2,
                     DoNotPlot = TRUE, shape = nodeShape2, ...)$layout
      rownames(lay2) <- rownames(adja2)

    }
    print("shuchu\n")
    print(colnames(adja2))
    #--------------------------------------------
    # rename taxa
    if(is.null(labels)){
      adja.tmp <- rename_taxa(adja2, toRename = "both", shortenLabels = shortenLabels,
                              labelLength = labelLength, labelPattern = labelPattern,
                              charToRm = charToRm)
      labels2 <- rownames(adja.tmp)

    } else if(is.logical(labels)){
      labels2 <- labels
    } else if(is.list(labels)){
      labels2 <- labels[[2]]
    } else if(is.vector(labels)){
      labels2 <- labels
    }

    #--------------------------------------------
    # filter edges (without influencing the layout)

    if(edgeInvisFilter != "none"){
      adja2 <- filter_edges(adja2, edgeFilter = edgeInvisFilter,
                           edgeFilterPar = edgeInvisPar)
    }

  }


  #==========================================================================
  ### plot network(s)


  if(twoNets){
      par(mfrow = c(1,2))


    if(sameLayout){
      if(layoutGroup == 1){
        lay2 <- lay1
      } else{
        lay1 <- lay2
      }
    }

    if(groupsChanged){
      q2 <- qgraph(adja2, color = nodecol2, layout = lay2, vsize = nodeSize2,
                   labels = labels2, label.scale = labelScale,
                   border.color = border2, border.width = borderWidth2,
                   label.font = labelFont2, label.cex = cexLabels,
                   edge.width = edgeWidth, edge.color = colmat2,
                   repulsion = repulsion, cut = cut2,
                   mar = mar, shape = nodeShape2, ...)
      if(showTitle) title(main = list(main2, cex = cexTitle))
    }

    q1 <- qgraph(adja1, color = nodecol1, layout = lay1, vsize = nodeSize1,
                 labels = labels1,label.scale = labelScale,
                 border.color = border1, border.width = borderWidth1,
                 label.font = labelFont1, label.cex = cexLabels,
                 edge.width = edgeWidth, edge.color = colmat1,
                 repulsion = repulsion, cut = cut1,
                 mar = mar,  shape = nodeShape1, ...)
    if(showTitle) title(main = list(main1, cex = cexTitle))

    if(!groupsChanged){
      q2 <-  qgraph(adja2, color = nodecol2, layout = lay2, vsize = nodeSize2,
                    labels = labels2, label.scale = labelScale,
                    border.color = border2, border.width = borderWidth2,
                    label.font = labelFont2, label.cex = cexLabels,
                    edge.width = edgeWidth, edge.color = colmat2,
                    repulsion = repulsion, cut = cut2,
                    mar = mar,  shape = nodeShape2, ...)
      if(showTitle) title(main = list(main2, cex = cexTitle))
    }

    #---------------------------------------------
  } else{

    q1 <- qgraph(adja1, color = nodecol1, layout = lay1, vsize = nodeSize1,
                 labels = labels1, label.scale = labelScale,
                 border.color = border1, border.width = borderWidth1,
                 label.font = labelFont1, label.cex = cexLabels,
                 edge.width = edgeWidth, edge.color = colmat1,
                 repulsion = repulsion, cut = cut1,
                 mar = mar,  shape = nodeShape1, ...)
    if(showTitle) title(main = list(main1, cex = cexTitle))
    q2 <- NULL
    labels2 <- NULL
  }

  layout.out <- list(layout1 = lay1, layout2 = lay2)
  nodecol.out <- list(nodecol1 = nodecol1, nodecol2 = nodecol2)
  labels.out <- list(labels1 = labels1, labels2 = labels2)


  output <- list(q1 = q1, q2 = q2, layout = layout.out,
                 nodecolor = nodecol.out, labels = labels.out)

  #------------------------------------------------------------------
  #------------------------------------------------------------------

  par(mfrow = opar$mfrow, mar = opar$mar)

  invisible(output)


}
