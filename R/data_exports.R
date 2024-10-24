#' Documentation for the data

#' Default databases by organism
#' 
#' @format Named character vector
#' @source GEO number
"defaultDbNames"

#' Annotations of Motifs to Transcription Factors (TFs)
#' 
#' @rdname motifAnnotations
#' @name motifAnnotations
#' @aliases motifAnnotations_hgnc motifAnnotations_mgi motifAnnotations_dmel
#' @title Annotations of Motifs to TFs
#' @description 
#' This dataset contains annotations linking motifs to transcription factors (TFs) 
#' for different organisms (Human, Mouse, Fly) based on motif collection version 9.
#' 
#' \bold{Version 9}: 
#' Annotations for motif collection version 9 ('mc9nr').
#' \itemize{
#'     \item{Human: motifAnnotations_hgnc_v9 (\emph{motifs-v9-nr.hgnc-m0.001-o0.0.tbl})}
#'     \item{Mouse: motifAnnotations_mgi_v9 (\emph{motifs-v9-nr.mgi-m0.001-o0.0.tbl})}
#'     \item{Fly: motifAnnotations_dmel_v9 (\emph{motifs-v9-nr.flybase-m0.001-o0.0.tbl})}
#' }
#' 
#' These objects are designed for use with RcisTarget without modification, 
#' but users can explore them to get more information about the motifs. 
#' The original files are available from \code{https://resources.aertslab.org/cistarget/motif2tf/}.
#' 
#' \bold{Columns:}
#' \itemize{
#'  \item{\bold{motif}: }{Motif ID.}
#'  \item{\bold{TF}: }{Transcription factor or inferred gene.}
#'  \item{\bold{directAnnotation}, \bold{inferred_Orthology}, \bold{inferred_MotifSimil}: }{Boolean values indicating whether the motif is annotated to the TF 
#'      directly ("directAnnotation") or inferred by orthology ("inferred_Orthology") or motif similarity ("inferred_MotifSimil").}
#'  \item{\bold{Description}: }{Description of the annotation source.}
#'  \item{\bold{annotationSource}: }{Source of the annotation as a factor. Levels: directAnnotation, 
#'      inferredBy_Orthology, inferredBy_MotifSimilarity, inferredBy_MotifSimilarity_n_Orthology.}
#' }
#' 
#' @docType data
#' @seealso importAnnotations, RcisTarget
#' @keywords datasets
NULL


