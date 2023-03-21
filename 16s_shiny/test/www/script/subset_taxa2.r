subset_taxa2<-function (physeq, phyloseq2) # large physeq,small physeq
{
  if (is.null(tax_table(physeq))) {
    cat("Nothing subset. No taxonomyTable in physeq.\n")
    return(physeq)
  }
  else {
    if (taxa_are_rows(physeq)) {
      oldMA <- as(tax_table(physeq), "matrix")
      oldDF <- data.frame(oldMA)
      newDF<-oldDF[rownames(as.data.frame(otu_table(phyloseq2))),]
      
      newMA <- as(newDF, "matrix")
    }
    else {
      oldMA <- as(tax_table(physeq), "matrix")
      oldDF <- data.frame(oldMA)
      newDF<-oldDF[rownames(as.data.frame(t(otu_table(phyloseq2)))),]
      
      newMA <- as(newDF, "matrix")
    }

    if (inherits(physeq, "taxonomyTable")) {
      return(tax_table(newMA))
    }
    else {
      tax_table(physeq) <- tax_table(newMA)
      return(physeq)
      
    }
  }
  
}
