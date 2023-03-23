# input, topic is the cistopicobject@binarized.cistopics, format should be "Topic1"
# return a data frame containing the enriched GO terms, their associated genes, and the fold enrichment and adjusted p-values.

# It takes a topic number as input and performs gene ontology (GO) enrichment analysis on the genes associated with that topic.
# It first extracts the binary vector of peak assignments for that topic from the cisTopicObject object,
# and uses it to retrieve the Ensembl gene IDs associated with those peaks.

# It then uses the enrichGO function from the clusterProfiler package to perform GO enrichment analysis on those genes
# using the "biological process" ontology.

# Finally, it calculates the fold enrichment for each significantly enriched GO term

# GO_topic_anno <- function(topic, cisTopicObject){
#   topic.list <- cisTopicObject@binarized.cisTopics
#   topic.1 <- topic.list[[topic]]
#   # get peaks and genes coord data frame
#   pg <- as.data.frame(cisTopicObject@region.data[["ENSEMBL"]], cisTopicObject@region.names)
#   colnames(pg) <- c("ENSEMBL")
#   geneID <- pg[rownames(topic.1), ]
#   anno.topic <- as.data.frame(geneID, rownames(topic.1))
#
#   # Run GO enrichment analysis
#   ego2 <- enrichGO(gene = anno.topic$geneID,
#                    OrgDb = org.Dr.eg.db,
#                    keyType = 'ENSEMBL',
#                    ont = "BP",
#                    pAdjustMethod = "BH",
#                    pvalueCutoff  = 0.01,
#                    qvalueCutoff  = 0.05)
#
#   godata <- as_tibble(ego2@result)
#   godata <- godata %>% filter(p.adjust < 0.05) %>% select(ID, Description, GeneRatio, BgRatio, p.adjust, geneID)
#
#
#   # Calculate enrichment score
#   results <- vector()
#   results.1 <- vector()
#
#   # Looping through each value in the "values" column of the data frame
#   for (i in 1:nrow(godata)) {
#     # Converting the value to a numerical value and performing division
#     result <- eval(parse(text = godata$GeneRatio[i]))/1
#     result.1 <- eval(parse(text = godata$BgRatio[i]))/1
#
#     # Converting the result to a double and appending it to the results vector
#     results <- c(results, as.double(result))
#     results.1 <- c(results.1, as.double(result.1))
#   }
#   godata$fold.enrichment <- results/results.1
#
#   godata$topic <- paste0(topic)
#
#   return(godata)
# }


# input a godata object (output of the GO_topic_anno function), a cisTopicObject, a topic number, and a Seurat object (seurat).
# returns the top motifs that are enriched in the given topic based on their p-value and fold enrichment.

# First extracts the Ensembl IDs of the genes that are enriched in the given topic using the godata object
# It then matches the Ensembl IDs to peak IDs in the cisTopicObject
# It uses these peak IDs to perform motif analysis on the Seurat object using the FindMotifs function.
# find_enriched_motifs <- function(godata, cisTopicObject, topic, seurat) {
#   gene <- godata$geneID
#   ensembl <- strsplit(gene, "/")
#   pg <- as.data.frame(cisTopicObject@region.data[["ENSEMBL"]], cisTopicObject@region.names)
#   colnames(pg) <- c("ENSEMBL")
#   topic.list <- cisTopicObject@binarized.cisTopics
#   topic.1 <- topic.list[[topic]]
#   geneID <- pg[rownames(topic.1), ]
#   anno.topic <- as.data.frame(geneID, rownames(topic.1))
#   anno.topic$peaks <- rownames(anno.topic)
#   peak.mf <- anno.topic[anno.topic$geneID %in% ensembl[[3]], "peaks" ]
#   peak.mf <- rownames(anno.topic)
#   topic.peak <- gsub(":", "-", peak.mf)
#   DefaultAssay(seurat) <- "peaks"
#   # need to setup background peaks
#   open.peaks <- AccessiblePeaks(seurat)
#   meta.feature <- GetAssayData(seurat, assay = "peaks", slot = "meta.features")
#   peaks.matched <- MatchRegionStats(
#     meta.feature = meta.feature[open.peaks, ],
#     query.feature = meta.feature[topic.peak, ],
#     n = 50000
#   )
#   enriched.motifs.topic <- FindMotifs(
#     object = seurat,
#     features = topic.peak,
#     background = peaks.matched
#
#   )
#   top.enriched.motifs.topic <- enriched.motifs.topic[enriched.motifs.topic$p.adjust < 0.05, ]
#   top.enriched.motifs.topic$log <- log10(top.enriched.motifs.topic$fold.enrichment)
#   return(top.enriched.motifs.topic)
# }
