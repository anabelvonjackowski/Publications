#source: https://github.com/joey711/phyloseq/issues/654
## Plot composition at taxaRank2 level within taxa taxaSet1 at taxaRank1 level
## Restricts plot to numberOfTaxa taxa
plot_composition <- function(physeq, taxaRank1 = "Phylum", taxaSet1 = "Proteobacteria",
                             taxaRank2 = "Family", numberOfTaxa = 9, fill = NULL,
                             x = "Sample", y = "Abundance", facet_grid = NULL) {
  ## Args:
  ## - physeq: phyloseq class object
  ## - taxaRank1: taxonomic level in which to do the first subsetting
  ## - taxaSet1: subset of level taxaRank1 to use
  ## - taxaRank2: taxonomic level used to agglomerate
  ## - numberOfTaxa: number of (top) taxa to keep at level taxaRank2
  ##
  ## Returns:
  ## - ggplot2 graphics
  ggdata <- ggformat(physeq, taxaRank1, taxaSet1, taxaRank2, numberOfTaxa)
  p <- ggplot(ggdata, aes_string(x = x, y = y, fill = fill, color = fill, group = "Sample"))
  ## Manually change color scale to assign grey to "Unknown" (if any)
  if (!is.null(fill) && any(c("Unknown", "Other") %in% unique(ggdata[, fill]))) {
    ranks <- as.character(unique(ggdata[, fill]))
    ranks <- ranks[ ! ranks %in% c("Unknown", "Other")]
    colvals <- c(gg_color_hue(length(ranks)), "grey45", "black")
    names(colvals) <- c(ranks, "Unknown", "Other")
    ## Now add the manually re-scaled layer with Unassigned as grey
    p <- p + scale_fill_manual(values=colvals) + scale_color_manual(values = colvals)
    
  }
  p <- p + geom_bar(stat = "identity", position = "stack")
  if ( !is.null(facet_grid)) {
    p <- p + facet_grid(facets = facet_grid, scales = "free_x")
  }
  p <- p + theme(axis.text.x=element_text(angle=90), axis.title.x=element_blank())
  p <- p + ggtitle(paste("Composition within", taxaSet1, "(", numberOfTaxa, "top", taxaRank2, ")"))
  return(p)
}
