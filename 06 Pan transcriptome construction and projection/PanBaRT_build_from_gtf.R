#!/usr/bin/env Rscript

# PanBaRT_build_from_gtf.R
# usage: Rscript PanBaRT_build_from_gtf.R <input_gtf> <prefix> [chimeric_tolerance=0.05]
# Note: input_gtf is the result of minimap2 --> genePred-->gtf, mapping the transcripts of GsRTD to the gtf of linear pan

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 2) stop("Usage: Rscript PanBaRT_build_from_gtf.R <input_gtf> <prefix> [chimeric_tolerance]")
gtf_file <- args[1]
prefix <- args[2]
chimeric_tolerance <- ifelse(length(args)>=3, as.numeric(args[3]), 0.05)

library(rtracklayer)
library(GenomicRanges)
library(IRanges)
library(stringr)
library(dplyr)

message("Importing GTF: ", gtf_file)
gr_all <- import(gtf_file) # GTF expected to contain exon features

# ---------------------------------------------------------------------
# Helper functions (implementations of authors' helpers)
# ---------------------------------------------------------------------
# getMultiExonTrans: return GRanges (exon entries) for transcripts with >1 exon
getMultiExonTrans <- function(gr){
  exons <- gr[gr$type %in% c("exon","CDS","utr","exon")] # keep exon-like entries
  if(length(exons)==0) return(GRanges())
  tx_to_exons <- split(exons, exons$transcript_id)
  multi <- tx_to_exons[sapply(tx_to_exons, length) > 1]
  if(length(multi)==0) return(GRanges())
  # reduce each transcript's exons to canonical representation
  reduced <- unlist(reduce(multi))
  # we need transcript_id metadata: attach transcript_id from list names
  # After reduce, widths may combine; set transcript_id according to grouping
  # create GRanges where names are transcript ids
  ids <- rep(names(multi), sapply(multi, function(x) 1))
  # But reduce above loses grouping; easier approach: for each transcript create reduced ranges and set metadata
  grlist <- GRangesList()
  for(tid in names(multi)){
    r <- reduce(multi[[tid]])
    r$transcript_id <- tid
    r$type <- "exon"
    r$gene_id <- NA
    grlist[[tid]] <- r
  }
  res <- unlist(grlist)
  return(res)
}

# getMonoExonTrans: return GRanges for transcripts with exactly 1 exon
getMonoExonTrans <- function(gr){
  exons <- gr[gr$type %in% c("exon","CDS","utr","exon")]
  tx_to_exons <- split(exons, exons$transcript_id)
  mono <- tx_to_exons[sapply(tx_to_exons, length) == 1]
  if(length(mono)==0) return(GRanges())
  grlist <- GRangesList()
  for(tid in names(mono)){
    r <- mono[[tid]]
    # keep metadata
    r$type <- "exon"
    r$transcript_id <- tid
    r$gene_id <- NA
    grlist[[tid]] <- r
  }
  res <- unlist(grlist)
  return(res)
}

# getIntronChain: create a chain string for each multi-exon transcript
getIntronChain <- function(multi_trans){
  # multi_trans is GRanges with $transcript_id and exon ranges (reduced)
  txlist <- split(multi_trans, multi_trans$transcript_id)
  chain <- sapply(names(txlist), function(tid){
    exs <- txlist[[tid]]
    exs <- exs[order(start(exs))]
    if(length(exs) <= 1) return(NA_character_)
    introns <- IRanges(start=end(exs)[-length(exs)] + 1, end=start(exs)[-1] - 1)
    paste0(as.character(seqnames(exs)[1]),":", paste0(start(introns),"-", end(introns), collapse=";"))
  }, USE.NAMES = TRUE)
  return(chain)
}

# getOverlapRange: for mono-exon transcripts, get merged loci (reduced)
getOverlapRange <- function(mono_trans){
  if(length(mono_trans)==0) return(GRanges())
  grlist <- split(mono_trans, mono_trans$group) # group must exist
  # But if no group preassigned, just reduce by seqnames
  reduced <- reduce(mono_trans)
  # create artificial gene_id for each reduced region
  reduced$gene_id <- paste0("mono_locus_", seq_along(reduced))
  return(reduced)
}

# getTransRange: return a GRanges of transcript ranges (start=min exon start, end=max exon end)
getTransRange <- function(gr_trans){
  # gr_trans is a GRanges of exons with $transcript_id
  txlist <- split(gr_trans, gr_trans$transcript_id)
  rngs <- lapply(names(txlist), function(tid){
    x <- txlist[[tid]]
    GRanges(seqnames = seqnames(x)[1],
            ranges = IRanges(start = min(start(x)), end = max(end(x))),
            strand = unique(as.character(strand(x))),
            transcript_id = tid)
  })
  res <- do.call(c, rngs)
  return(res)
}

# assignGeneID: simple assignment based on seqnames + start order
assignGeneID <- function(gr){
  # gr is a GRanges where each record corresponds to an exon of a "collapsed transcript"
  # We will define genes by reducing gr by overlap and then assign unique IDs per locus
  tranges <- getTransRange(gr)
  merged <- reduce(tranges)
  # For each transcript, find which merged locus it overlaps
  hits <- findOverlaps(tranges, merged)
  map <- data.frame(transcript_id = tranges$transcript_id[queryHits(hits)],
                    locus = subjectHits(hits), stringsAsFactors = FALSE)
  # create gene names
  genes <- paste0("gene_", seq_along(merged))
  mapping <- setNames(genes[map$locus], map$transcript_id)
  # assign gene_id back to gr (per exon)
  # Build a vector mapping per transcript id
  gr$gene_id <- mapping[gr$transcript_id]
  return(gr)
}

# getIntronicGene: transcripts entirely within intron region of another transcript
getIntronicGene <- function(query, subject){
  # query,subject: gr where metadata transcript_id exists
  # idea: build introns for multi-exon subject transcripts and test if query transcript ranges are within those introns entirely
  subj_multi <- getMultiExonTrans(subject)
  if(length(subj_multi)==0) return(GRanges())
  # build intron ranges per subject transcript
  txlist <- split(subj_multi, subj_multi$transcript_id)
  intron_gr <- GRangesList()
  for(tid in names(txlist)){
    exs <- txlist[[tid]][order(start(txlist[[tid]]))]
    if(length(exs) <= 1) next
    introns <- IRanges(start = end(exs)[-length(exs)] + 1, end = start(exs)[-1] - 1)
    if(length(introns)>0){
      g <- GRanges(seqnames = seqnames(exs)[1], ranges = introns, strand = strand(exs)[1])
      g$parent_tx <- tid
      intron_gr[[tid]] <- g
    }
  }
  intron_all <- unlist(intron_gr)
  if(length(intron_all)==0) return(GRanges())

  q_ranges <- getTransRange(query)
  # check if q_ranges are within any intron of intron_all
  hits <- findOverlaps(q_ranges, intron_all, type = "within")
  if(length(hits) == 0) return(GRanges())
  sel_q <- q_ranges[unique(queryHits(hits))]
  return(sel_q)
}

psetdiff <- function(gr, chimeric_list){
  # gr: GRanges of exons for many transcripts, with $transcript_id metadata
  # chimeric_list: a GRangesList (keyed by transcript_id) of intervals to remove for each transcript
  require(GenomicRanges)

  out_list <- GRangesList()
  txlist <- split(gr, gr$transcript_id)

  for(tid in names(txlist)){
    tx_gr <- txlist[[tid]]

    # get the to-remove intervals for this transcript
    # chimeric_list may be a GRangesList or ordinary list; handle both
    torem <- NULL
    if(!is.null(chimeric_list[[tid]])) {
      torem <- chimeric_list[[tid]]
    } else if(!is.null(chimeric_list[tid])) {
      torem <- chimeric_list[tid][[1]]
    }

    # If nothing to remove, keep original
    if(is.null(torem) || length(torem) == 0){
      out_list[[tid]] <- tx_gr
      next
    }

    # Ensure 'torem' is a GRanges (unlist if it's a GRangesList)
    if(inherits(torem, "GRangesList") || inherits(torem, "CompressedGRangesList")){
      torem <- unlist(torem, use.names = FALSE)
    }
    if(!inherits(torem, "GRanges")){
      # try to coerce, otherwise skip removal but warn
      warning(sprintf("psetdiff: torem for %s is not GRanges (class=%s). Skipping removal for this transcript.", tid, class(torem)[1]))
      out_list[[tid]] <- tx_gr
      next
    }

    # Use GenomicRanges::setdiff explicitly to avoid dispatch to base::setdiff
    new <- GenomicRanges::setdiff(tx_gr, torem)

    if(length(new) > 0){
      new$transcript_id <- tid
      out_list[[tid]] <- new
    } else {
      # if nothing remains after cutting, we skip (no exons left)
      # optionally create zero-length placeholder? we skip here.
    }
  }

  # unlist and return; preserve GRanges type
  if(length(out_list) == 0) return(GRanges())
  res <- unlist(out_list, use.names = FALSE)
  return(res)
}
# ---- end psetdiff ----

# ---------------------------------------------------------------------
# Step 1: Processing of multi-exon and single-exon transcripts
# ---------------------------------------------------------------------
message("Step1: split multi-exon and mono-exon transcripts and collapse identical intron-chains")

gr <- gr_all
# keep only exons
gr <- gr[gr$type == "exon"]

# Multi-exon transcripts
multi_trans <- getMultiExonTrans(gr)
message("Multi-exon transcripts (reduced): ", length(unique(multi_trans$transcript_id)))

intron_chain <- getIntronChain(multi_trans)
# attach group key
multi_trans$group <- intron_chain[multi_trans$transcript_id]
# If NA (some how), set group to transcript id
multi_trans$group[is.na(multi_trans$group)] <- multi_trans$transcript_id[is.na(multi_trans$group)]
multi_trans_group <- split(multi_trans$transcript_id, multi_trans$group)
multi_trans_group <- lapply(multi_trans_group, unique)

# collapse each group: reduce ranges per group => one "representative" transcript region
multi_split <- split(multi_trans, multi_trans$group)
multi_split <- endoapply(multi_split, function(x) reduce(x))
multi_trans_collapsed <- unlist(multi_split)
multi_trans_collapsed$type <- 'exon'
multi_trans_collapsed$gene_id <- names(multi_trans_collapsed)
multi_trans_collapsed$transcript_id <- names(multi_trans_collapsed)
names(multi_trans_collapsed) <- NULL

# Mono-exon transcripts
mono_trans <- getMonoExonTrans(gr)
message("Mono-exon transcripts: ", length(unique(mono_trans$transcript_id)))
# get loci by reducing mono exons
mono_loci <- reduce(mono_trans)
# build mapping: find which mono exon lies within which locus (within True)
hits <- findOverlaps(mono_trans, mono_loci, type = "within")
if(length(hits)>0){
  mono_trans <- mono_trans[queryHits(hits)]
  mono_trans$group <- paste0("mono_locus_", subjectHits(hits))
} else {
  # if none within (rare), group by each transcript
  mono_trans$group <- mono_trans$transcript_id
}
mono_trans_group <- split(mono_trans$transcript_id, mono_trans$group)
mono_trans_group <- lapply(mono_trans_group, unique)

mono_split <- split(mono_trans, mono_trans$group)
mono_split <- endoapply(mono_split, function(x) reduce(x))
mono_trans_collapsed <- unlist(mono_split)
mono_trans_collapsed$type <- 'exon'
mono_trans_collapsed$gene_id <- names(mono_trans_collapsed)
mono_trans_collapsed$transcript_id <- names(mono_trans_collapsed)
names(mono_trans_collapsed) <- NULL

# combine
trans_group <- c(mono_trans_group, multi_trans_group)
trans_group <- sapply(trans_group, function(x) paste0(unique(x), collapse = ';'))
names(trans_group) <- paste0('trans-', names(trans_group))

gr_collapsed <- c(multi_trans_collapsed, mono_trans_collapsed)
gr_collapsed$transcript_id <- paste0('trans-', gr_collapsed$transcript_id)
gr_collapsed$source_trans <- trans_group[gr_collapsed$transcript_id]
gr_collapsed <- sort(gr_collapsed, by = ~ seqnames + start + end)

# ---------------------------------------------------------------------
# Step 2: refine genes (chimeric & intronic detection) and assign new gene ids
# ---------------------------------------------------------------------
gr <- gr_collapsed
gr$gene_id0 <- gr$gene_id
gr$transcript_id0 <- gr$transcript_id
gr$observation <- NA

trans_range <- getTransRange(gr)
hits <- findOverlaps(query = trans_range, subject = trans_range)
idx <- which(queryHits(hits) != subjectHits(hits))
hits <- hits[idx,]

message('Refine chimeric genes')

gr1 <- trans_range[queryHits(hits)]
gr2 <- trans_range[subjectHits(hits)]
overlaps <- pintersect(gr1,gr2)
overlaps$pair1 <- gr1$transcript_id
overlaps$pair2 <- gr2$transcript_id

p1 <- width(overlaps)/width(gr1)
p2 <- width(overlaps)/width(gr2)
idx <- which(p1 < chimeric_tolerance & p2 < chimeric_tolerance)

if(length(idx) > 0){
  chemric <- overlaps[idx]

  # Group by pair1, reduce() each transcript separately and label the transcript_id, then merge
  pair1_split <- split(chemric, chemric$pair1)
  pair1_list <- GRangesList()
  for(tid in names(pair1_split)){
    # reduce this transcript's pieces
    r <- reduce(pair1_split[[tid]])
    if(length(r) > 0){
      r$transcript_id <- tid
      pair1_list[[tid]] <- r
    }
  }
  pair1_combined <- if(length(pair1_list) > 0) unlist(pair1_list, use.names = FALSE) else GRanges()

  # Group by pair2 and handle in a similar way
  pair2_split <- split(chemric, chemric$pair2)
  pair2_list <- GRangesList()
  for(tid in names(pair2_split)){
    r <- reduce(pair2_split[[tid]])
    if(length(r) > 0){
      r$transcript_id <- tid
      pair2_list[[tid]] <- r
    }
  }
  pair2_combined <- if(length(pair2_list) > 0) unlist(pair2_list, use.names = FALSE) else GRanges()

  chemric_combined <- c(pair1_combined, pair2_combined)
  if(length(chemric_combined) == 0){
    message("No chimeric intervals after reduction.")
  } else {
    chemric_combined <- sort(chemric_combined, by = ~ seqnames + start + end)

    gr_chemric <- gr[gr$transcript_id %in% chemric_combined$transcript_id]
    gr_no_chemric <- gr[!(gr$transcript_id %in% chemric_combined$transcript_id)]

    chemric_by_tx <- split(chemric_combined, chemric_combined$transcript_id)

    keep_ids <- unique(gr_chemric$transcript_id)
    chemric_by_tx <- chemric_by_tx[names(chemric_by_tx) %in% keep_ids]

    gr_chemric_cut <- psetdiff(gr_chemric, chemric_by_tx)
    gr_chemric_cut <- unlist(split(gr_chemric_cut, gr_chemric_cut$transcript_id), use.names = FALSE)
    if(length(gr_chemric_cut) > 0){
      gr_chemric_cut$transcript_id <- gr_chemric_cut$transcript_id
      gr_chemric_cut$gene_id <- gr_chemric_cut$transcript_id
      gr_chemric_cut$type <- 'exon'
    }

    gr_cut <- c(gr_chemric_cut, gr_no_chemric)
    gr_cut <- sort(gr_cut, by = ~ seqnames + start + end)
    gr_cut <- assignGeneID(gr = gr_cut)
    mapping <- data.frame(transcript_id=gr_cut$transcript_id,
                          gene_id=gr_cut$gene_id)
    mapping <- unique(mapping)
    rownames(mapping) <- mapping$transcript_id
    gr$gene_id <- mapping[gr$transcript_id,'gene_id']
  }
}
## ------------------ End chimeric fix ------------------

# Intronic genes
message('Refine intronic genes...')
gr_intronic <- getIntronicGene(query = gr, subject = gr)
if(length(gr_intronic) > 0){
  gr_intronic <- assignGeneID(gr_intronic)
  gr_intronic$gene_id <- paste0('intronic-', gr_intronic$gene_id)
  mapping_i <- data.frame(transcript_id = gr_intronic$transcript_id,
                          gene_id = gr_intronic$gene_id, stringsAsFactors = FALSE)
  mapping_i <- unique(mapping_i)
  rownames(mapping_i) <- mapping_i$transcript_id
  idx_i <- which(gr$transcript_id %in% gr_intronic$transcript_id)
  gr$observation[idx_i] <- 'intronic'
  gr$gene_id[idx_i] <- mapping_i[gr$transcript_id[idx_i],'gene_id']
}

# Generate new gene ids: prefix + seqname + zero-padded index
message('Generate new gene ids...')
gr <- sort(gr, by = ~ seqnames + start + end)
genes <- unique(gr$gene_id)
pad_n <- nchar(as.character(length(genes)))
genes_n <- stringr::str_pad(string = 1:length(genes), width = pad_n, pad = '0')
names(genes_n) <- genes
gene_id_new <- paste0(prefix, seqnames(gr), genes_n[gr$gene_id])
gr$gene_id <- gene_id_new

mapping <- data.frame(transcript_id = gr$transcript_id,
                      gene_id = gene_id_new, stringsAsFactors = FALSE)
mapping <- unique(mapping)
mapping <- mapping[order(mapping$gene_id),]
n <- sequence(rle(mapping$gene_id)$lengths)
transcript_id_new <- paste0(mapping$gene_id, '.', n)
names(transcript_id_new) <- mapping$transcript_id
gr$transcript_id <- transcript_id_new[gr$transcript_id]

# Annotate observations (multi/mono/exonic)
gr_multi <- getMultiExonTrans(gr)
idx_multi <- which(gr$transcript_id %in% gr_multi$transcript_id)
if(length(idx_multi) > 0) gr$observation[idx_multi] <- paste0(gr$observation[idx_multi], ';multi-exon')

gr_mono <- getMonoExonTrans(gr)
idx_mono <- which(gr$transcript_id %in% gr_mono$transcript_id)
if(length(idx_mono) > 0) gr$observation[idx_mono] <- paste0(gr$observation[idx_mono], ';mono-exon')

# exonic transcripts: mono exon transcripts that are within multi-exon transcripts
trans_range_mono <- getTransRange(gr_mono)
hits_within <- findOverlaps(trans_range_mono, gr, type = 'within')
if(length(hits_within)>0){
  mapping_exonic <- data.frame(from = trans_range_mono$transcript_id[queryHits(hits_within)],
                               to = gr$transcript_id[subjectHits(hits_within)],
                               stringsAsFactors = FALSE)
  mapping_exonic <- unique(mapping_exonic)
  trans_in_exonic <- unique(mapping_exonic$from[mapping_exonic$from != mapping_exonic$to])
  idx_exonic <- which(gr$transcript_id %in% trans_in_exonic)
  if(length(idx_exonic)>0) gr$observation[idx_exonic] <- paste0(gr$observation[idx_exonic], ';exonic-trans')
}

gr$observation <- gsub('NA;','', gr$observation)

# export final PanBaRT GTF and mapping
library(rtracklayer)
export(gr, paste0(prefix, "_PanBaRT.gtf"))

# mapping2: source transcripts -> PanBaRT transcripts mapping
# reconstruct mapping2 using gr$source_trans (which maps to original transcripts)
mapping2 <- data.frame(PanBaRT_transcript_id = gr$transcript_id,
                       PanBaRT_gene_id = gr$gene_id,
                       source_transcripts = gr$source_trans,
                       observation = gr$observation,
                       seqnames = as.character(seqnames(gr)),
                       start = start(gr),
                       end = end(gr),
                       strand = as.character(strand(gr)),
                       stringsAsFactors = FALSE)
write.table(mapping2, file = paste0(prefix,"_GsRTD_and_PanBaRT_match.tsv"),
            row.names = FALSE, sep = "\t", quote = FALSE)

message("PanBaRT build finished. Output files:")
message(paste0(prefix, "_PanBaRT.gtf"))
message(paste0(prefix, "_GsRTD_and_PanBaRT_match.tsv"))
