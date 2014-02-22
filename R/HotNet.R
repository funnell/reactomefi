#' NCI MAF to genescore
#'
#' Converts a NCI Mutation Annotation File to a two column data.frame where
#'  the first column contains gene names and the second column contains scores
#'  representing the proportion of samples in which a gene is mutated.
#'
#' @param maf NCI Mutation Annotation File.
#' @return data.frame Two column data.frame containing gene names and scores.
maf2genescores <- function(maf) {
    allowed.types <- c("frame_shift_del", "in_frame_del", "in_frame_ins",
                       "frame_shift_ins", "nonsense_mutation",
                       "nonstop_mutation", "missense_mutation",
                       "splice_site_ins", "splice_site_snp", "splice_site_snp",
                       "splice_site_del", "stopcodon_dnp", "splice_site",
                       "stop_codon_del", "init_met_del")
    maf <- subset(maf, tolower(Variant_Classification) %in% allowed.types)

    maf <- maf[c("Hugo_Symbol", "Tumor_Sample_Barcode")]
    maf[, "Tumor_Sample_Barcode"] <- substr(maf$Tumor_Sample_Barcode, 1, 12)
    maf <- maf[!duplicated(maf), ]

    num.samples <- length(unique(maf$Tumor_Sample_Barcode)) 
    genes <- unique(as.character(maf$Hugo_Symbol))
    genescores <- data.frame(gene=genes, score=rep(0, length(genes)))
    for (gene in genes) {
        gene.samples <- subset(maf, Hugo_Symbol == gene, "Tumor_Sample_Barcode")
        num.gene.samples <- nrow(gene.samples)
        score <- num.gene.samples / num.samples
        genescores[genescores["gene"] == gene, "score"] <- score
    }
    return(genescores)
}

#' @rdname service-methods
#' @aliases service,HotNet-method
setMethod("service", signature("HotNet"), function(object) {
    return(object@service)
})

#' @rdname modules-methods
#' @aliases modules,HotNet-method
setMethod("modules", signature("HotNet"), function(object) {
    return(object@modules)
})

#' @rdname subnet-methods
#' @aliases subnet,HotNet-method
setMethod("subnet", signature("HotNet"), function(object, fdr = 0.05) {
    mask <- sapply(modules(object), function(x) x$fdr <= fdr)
    sub.mods <- modules(object)[mask]
    if (length(sub.mods) == 0) {
        return(new("ReactomeFINetwork", service = service(object)))
    }

    mod.genes <- data.frame()
    for (i in seq(length(sub.mods))) {
        genes.df <- data.frame(gene = sub.mods[[i]]$genes, module = i,
                               stringsAsFactors=F)
        mod.genes <- rbind(mod.genes, genes.df)
    }
    network <- new("ReactomeFINetwork", service = service(object))
    network <- build(network, mod.genes$gene)
    modules(network) <- mod.genes
    return(network)
})
