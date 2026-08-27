#' Run Cell-Cell Communication analysis
#' 
#' This function acts as a unified entry point for different cell-cell communication (CCI) 
#' methods (e.g., CellChat). It delegates to the specific backend method and stores a 
#' standardized `cci_interactions` data frame in the `SingleCellExperiment` object.
#'
#' @title RunCCI
#' @param object A SingleCellExperiment object.
#' @param method Character, the CCI method to use. Currently supports `"CellChat"`.
#' @param keep_object Logical, whether to store the original backend object (e.g., CellChat object) 
#' in the metadata. Defaults to `FALSE` to prevent the `SingleCellExperiment` object from becoming too large.
#' @param ... Arguments passed to the specific backend (e.g., `RunCellChat`).
#' @return A SingleCellExperiment object with CCI results stored in the metadata.
#' @export
RunCCI <- function(object, method = c("CellChat", "CellPhoneDB", "NicheNet"), keep_object = FALSE, ...) {
    method <- match.arg(method)
    if (method == "CellChat") {
        # Force return="sce" for RunCCI to maintain standard SCE workflow
        args <- list(...)
        if (is.null(args$return)) {
            args$return <- "sce"
        }
        return(do.call(RunCellChat, c(list(sce = object, keep_object = keep_object), args)))
    } else if (method == "CellPhoneDB") {
        return(RunCellPhoneDB(object, ...))
    } else if (method == "NicheNet") {
        return(RunNicheNet(object, ...))
    }
}

#' Run CellChat communication analysis
#' 
#' https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat-vignette.html#c-starting-from-a-singlecellexperiment-object
#' 
#' @title RunCellChat
#' @param sce A SingleCellExperiment object.
#' @param group Grouping variable used for communication analysis. If `NULL`,
#' sclet resolves it from `ActiveIdent(sce)`. When the active identity is
#' `colLabels`, sclet uses a temporary metadata column derived from
#' `SingleCellExperiment::colLabels(sce)`.
#' @param assay_name the name of assay. If NULL, sclet resolves it from the
#' selected `layer`.
#' @param layer layer used for communication analysis. If NULL, use
#' `DefaultLayer(sce)`.
#' @param species one of human, mouse
#' @param db_item one of  "all", "except Non-protein","Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact" and "Non-protein Signaling"
#' @param type one of  "triMean", "truncatedMean", "thresholdedMean", "median"
#' @param trim the fraction (0 to 0.25) of observations to be trimmed from each end of x before the mean is computed
#' @param min.cells min cells to filter
#' @param name communication record id. Defaults to `"cellchat"`.
#' @param return one of `"cellchat"`, `"sce"` or `"both"`
#' @param keep_object Logical, whether to store the complete CellChat object in the SCE metadata. Defaults to FALSE.
#' @return by default a CellChat object; if `return = "sce"`, returns the
#' updated SingleCellExperiment; if `return = "both"`, returns a list with both
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment assayNames
#' @export
RunCellChat <- function(sce, group = NULL,
                        assay_name = NULL,
                        layer = NULL,
                        species = "human",
                        db_item = c("Secreted Signaling"),
                        type = "triMean",
                        trim = 0.1,
                        min.cells = 10,
                        name = "cellchat",
                        return = c("cellchat", "sce", "both"),
                        keep_object = FALSE) {

    if (!requireNamespace("CellChat", quietly = TRUE)) {
        stop("Package 'CellChat' is needed for this function to work. Please install it.")
    }

    species <- match.arg(species, c("human", "mouse", "zebrafish"))
    db_item <- match.arg(db_item, c("all", "except Non-protein",
                                    "Secreted Signaling", "ECM-Receptor", 
                                    "Cell-Cell Contact" , "Non-protein Signaling"))
    return <- match.arg(return)
    type <- match.arg(type, c("triMean", "truncatedMean", "thresholdedMean", "median"))

    if(! "logcounts" %in% assayNames(sce)){
        sce <- NormalizeData(sce)
    }

    source <- sclet_resolve_expression_source(
        object = sce,
        layer = layer,
        assay = assay_name,
        prefer_nonscaled = TRUE,
        context = "CellChat analysis"
    )
    assay_name <- source$assay
    data <- assay(sce, assay_name)
    meta <- as.data.frame(colData(sce))
    resolved_group <- sclet_resolve_group_column(sce, meta = meta, group = group)
    meta <- resolved_group$meta
    group <- resolved_group$group

    cellchat_obj <- CellChat::createCellChat(object = data, 
                                             meta = meta, 
                                             group.by = group)


    db <- switch(species,
        human = CellChat::CellChatDB.human,
        mouse = CellChat::CellChatDB.mouse,
        zebrafish = CellChat::CellChatDB.zebrafish
    )

    if (db_item == "all"){
        db.use <- db
    } else if(db_item == "except Non-protein"){
        db.use <- CellChat::subsetDB(db)
    } else {
        db.use <- CellChat::subsetDB(db, search = db_item, key = "annotation")
    }


    cellchat_obj@DB <- db.use
    cellchat_obj <- CellChat::subsetData(cellchat_obj)

    cellchat_obj <- CellChat::identifyOverExpressedGenes(cellchat_obj)
    cellchat_obj <- CellChat::identifyOverExpressedInteractions(cellchat_obj)


    cellchat_obj <- CellChat::computeCommunProb(cellchat_obj, type = type)
    cellchat_obj <- CellChat::filterCommunication(cellchat_obj, min.cells = min.cells)

    cellchat_obj <- CellChat::computeCommunProbPathway(cellchat_obj)
    cellchat_obj <- CellChat::aggregateNet(cellchat_obj)

    # Extract standard CCI interactions dataframe
    df <- CellChat::subsetCommunication(cellchat_obj)
    cci_interactions <- data.frame(
        source_group = df$source,
        target_group = df$target,
        ligand = df$ligand,
        receptor = df$receptor,
        pathway = df$pathway_name,
        probability = df$prob,
        p_value = df$pval,
        method = "CellChat",
        stringsAsFactors = FALSE
    )

    sce <- sclet_set_analysis(
        sce,
        "cellchat",
        list(
            id = name,
            method = "CellChat",
            object = if (keep_object) cellchat_obj else NULL,
            interactions = cci_interactions,
            group = group,
            assay = assay_name,
            layer = source$layer,
            species = species,
            db_item = db_item,
            params = list(
                type = type,
                trim = trim,
                min.cells = min.cells
            )
        )
    )
    sce <- sclet_set_analysis_state(
        object = sce,
        type = "communication",
        id = name,
        method = "CellChat",
        inputs = list(
            assay = assay_name,
            layer = source$layer,
            group = group
        ),
        artifacts = list(
            analysis_key = "cellchat"
        ),
        params = list(
            species = species,
            db_item = db_item,
            type = type,
            trim = trim,
            min.cells = min.cells
        ),
        summary = list(
            n_groups = length(unique(meta[[group]]))
        ),
        active = TRUE
    )
    sce <- sclet_log_command(
        sce,
        "RunCellChat",
        params = list(
            group = group,
            assay_name = assay_name,
            layer = source$layer,
            species = species,
            db_item = db_item,
            type = type,
            trim = trim,
            min.cells = min.cells,
            name = name
        ),
        outputs = list(
            analysis = "communication",
            communication = name
        )
    )

    if (identical(return, "sce")) {
        return(sce)
    }
    if (identical(return, "both")) {
        return(list(
            sce = sce,
            cellchat = cellchat_obj
        ))
    }

    attr(cellchat_obj, "sce") <- sce
    return(cellchat_obj)
}

#' Run CellPhoneDB communication analysis
#' 
#' @title RunCellPhoneDB
#' @param sce A SingleCellExperiment object.
#' @param group Grouping variable used for communication analysis. If `NULL`,
#' sclet resolves it from `ActiveIdent(sce)`.
#' @param assay_name the name of assay. If NULL, sclet resolves it from the selected `layer`.
#' @param layer layer used for communication analysis. If NULL, use `DefaultLayer(sce)`.
#' @param species one of human, mouse. CellPhoneDB is primarily designed for human, but sclet can map mouse genes.
#' @param name communication record id. Defaults to `"cellphonedb"`.
#' @param threads Number of threads to use. Defaults to 4.
#' @param iterations Number of iterations for statistical analysis. Defaults to 1000.
#' @param pvalue P-value threshold. Defaults to 0.05.
#' @param ... Additional parameters.
#' @return A SingleCellExperiment object with updated CCI state.
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment assayNames
#' @importFrom basilisk basiliskStart basiliskRun basiliskStop
#' @export
RunCellPhoneDB <- function(sce, group = NULL,
                           assay_name = NULL,
                           layer = NULL,
                           species = "human",
                           name = "cellphonedb",
                           threads = 4,
                           iterations = 1000,
                           pvalue = 0.05,
                           ...) {

    species <- match.arg(species, c("human", "mouse"))
    
    if(! "logcounts" %in% assayNames(sce)){
        sce <- NormalizeData(sce)
    }

    source <- sclet_resolve_expression_source(
        object = sce,
        layer = layer,
        assay = assay_name,
        prefer_nonscaled = TRUE,
        context = "CellPhoneDB analysis"
    )
    assay_name <- source$assay
    
    meta <- as.data.frame(colData(sce))
    resolved_group <- sclet_resolve_group_column(sce, meta = meta, group = group)
    meta <- resolved_group$meta
    group <- resolved_group$group

    # Prepare inputs for Python
    # CellPhoneDB requires a meta.txt (cell, cell_type) and counts.txt (gene, cell)
    # We will pass these via basilisk and anndata
    
    mat <- assay(sce, assay_name)
    gene_names <- rownames(mat)
    
    # If mouse, we ideally need ortholog mapping, but for simplicity we will just uppercase for now
    # Note: A robust ortholog mapping should be implemented here for production
    if (species == "mouse") {
        gene_names <- toupper(gene_names)
    }
    
    cell_meta <- data.frame(
        cell = colnames(mat),
        cell_type = as.character(meta[[group]]),
        stringsAsFactors = FALSE
    )
    
    # Execute Python CellPhoneDB via basilisk
    message("Running CellPhoneDB via basilisk...")
    message("(First run will take several minutes to set up the Python environment)")
    res_list <- sclet_basilisk_run(
        sclet_cellphonedb_env,
        function(mat, gene_names, cell_meta, threads, iterations, pvalue) {
        if (!requireNamespace("reticulate", quietly = TRUE)) {
            stop("Package 'reticulate' is needed for this function to work. Please install it.")
        }
        
        # Using tempdir for CPDB outputs
        out_dir <- tempfile(pattern = "cpdb_out")
        dir.create(out_dir)
        out_dir <- normalizePath(out_dir, winslash = "/", mustWork = FALSE)
        
        write_counts_csv <- function(x, genes, file) {
            con <- file(file, open = "wt")
            on.exit(close(con), add = TRUE)

            writeLines(
                paste(c("", colnames(x)), collapse = ","),
                con = con
            )

            for (i in seq_len(nrow(x))) {
                values <- as.vector(x[i, , drop = TRUE])
                writeLines(
                    paste(c(genes[i], values), collapse = ","),
                    con = con
                )
            }
        }

        # Write inputs
        counts_file <- file.path(out_dir, "counts.csv")
        meta_file <- file.path(out_dir, "meta.csv")

        # CellPhoneDB consumes a dense CSV file, but we can stream rows to disk
        # instead of materializing the full expression matrix in memory first.
        write_counts_csv(mat, gene_names, counts_file)
        utils::write.csv(cell_meta, meta_file, row.names = FALSE, quote = FALSE)
        
        # Python execution
        reticulate::py_run_string(paste0("
import os
import pandas as pd
from cellphonedb.utils import db_utils
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

cpdb_out_dir = '", out_dir, "'
meta_file = os.path.join(cpdb_out_dir, 'meta.csv')
counts_file = os.path.join(cpdb_out_dir, 'counts.csv')

cpdb_dir = os.path.join(os.path.expanduser('~'), '.cpdb')
db_version = 'v5.0.0'
db_file_path = os.path.join(cpdb_dir, 'cellphonedb.zip')

if not os.path.exists(db_file_path):
    print('Downloading CellPhoneDB database...')
    db_utils.download_database(cpdb_dir, db_version)

cpdb_statistical_analysis_method.call(
    cpdb_file_path=db_file_path,
    meta_file_path=meta_file,
    counts_file_path=counts_file,
    counts_data='hgnc_symbol',
    output_path=cpdb_out_dir,
    threads=", threads, ",
    iterations=", iterations, ",
    pvalue=", pvalue, ",
    output_suffix='sclet'
)

# Read results
pvals_path = os.path.join(cpdb_out_dir, 'statistical_analysis_pvalues_sclet.txt')
means_path = os.path.join(cpdb_out_dir, 'statistical_analysis_means_sclet.txt')

if os.path.exists(pvals_path) and os.path.exists(means_path):
    pvals_df = pd.read_csv(pvals_path, sep='\t')
    means_df = pd.read_csv(means_path, sep='\t')
else:
    pvals_df = pd.DataFrame()
    means_df = pd.DataFrame()
        "))
        
        # Retrieve results
        pvals <- reticulate::py$pvals_df
        means <- reticulate::py$means_df
        
        list(pvalues = pvals, means = means)
    },
        mat = mat,
        gene_names = gene_names,
        cell_meta = cell_meta,
        threads = threads,
        iterations = iterations,
        pvalue = pvalue
    )
    
    # Parse results into standard format
    pvals <- res_list$pvalues
    means <- res_list$means
    
    if (is.null(pvals) || nrow(pvals) == 0) {
        stop("CellPhoneDB failed to produce valid output.")
    }
    
    # Extract cell pair columns
    pair_cols <- grep("\\|", colnames(pvals), value = TRUE)
    
    interactions_list <- list()
    for (pair in pair_cols) {
        cells <- strsplit(pair, "\\|")[[1]]
        source_cell <- cells[1]
        target_cell <- cells[2]
        
        sub_pvals <- pvals[[pair]]
        sub_means <- means[[pair]]
        
        valid_idx <- which(sub_pvals <= pvalue)
        if (length(valid_idx) > 0) {
            
            ligand_val <- if ("gene_a" %in% names(pvals)) {
                ifelse(is.na(pvals$gene_a[valid_idx]) | pvals$gene_a[valid_idx] == "", pvals$partner_a[valid_idx], pvals$gene_a[valid_idx])
            } else {
                pvals$interacting_pair[valid_idx]
            }
            
            receptor_val <- if ("gene_b" %in% names(pvals)) {
                ifelse(is.na(pvals$gene_b[valid_idx]) | pvals$gene_b[valid_idx] == "", pvals$partner_b[valid_idx], pvals$gene_b[valid_idx])
            } else {
                pvals$interacting_pair[valid_idx]
            }
            
            interactions_list[[pair]] <- data.frame(
                source_group = source_cell,
                target_group = target_cell,
                ligand = ligand_val,
                receptor = receptor_val,
                pathway = pvals$interacting_pair[valid_idx],
                probability = sub_means[valid_idx],
                p_value = sub_pvals[valid_idx],
                method = "CellPhoneDB",
                stringsAsFactors = FALSE
            )
        }
    }
    
    if (length(interactions_list) > 0) {
        cci_interactions <- do.call(rbind, interactions_list)
        rownames(cci_interactions) <- NULL
    } else {
        cci_interactions <- data.frame()
    }
    
    # Set state
    sce <- sclet_set_analysis(
        sce,
        "cellphonedb",
        list(
            id = name,
            method = "CellPhoneDB",
            object = NULL,
            interactions = cci_interactions,
            group = group,
            assay = assay_name,
            layer = source$layer,
            species = species,
            params = list(
                threads = threads,
                iterations = iterations,
                pvalue = pvalue
            )
        )
    )
    
    sce <- sclet_set_analysis_state(
        object = sce,
        type = "communication",
        id = name,
        method = "CellPhoneDB",
        inputs = list(
            assay = assay_name,
            layer = source$layer,
            group = group
        ),
        artifacts = list(
            analysis_key = "cellphonedb"
        ),
        params = list(
            species = species,
            threads = threads,
            iterations = iterations,
            pvalue = pvalue
        ),
        summary = list(
            n_groups = length(unique(meta[[group]]))
        ),
        active = TRUE
    )
    
    return(sce)
}

#' Run NicheNet communication analysis
#' 
#' @title RunNicheNet
#' @param sce A SingleCellExperiment object.
#' @param group Grouping variable used for communication analysis. If `NULL`, sclet resolves it from `ActiveIdent(sce)`.
#' @param receiver Character, name of the receiver cell group.
#' @param sender Character or character vector, name(s) of sender cell groups. Defaults to `"all"`.
#' @param genes_oi Character vector of genes of interest in the receiver. If `NULL`, the function will use `scran::findMarkers` to find DEGs for the receiver cluster vs others.
#' @param ligand_target_matrix NicheNet ligand-target matrix. Required.
#' @param lr_network NicheNet ligand-receptor network. Required.
#' @param weighted_networks NicheNet weighted networks. Required.
#' @param assay_name the name of assay. If NULL, sclet resolves it from the selected `layer`.
#' @param layer layer used for communication analysis. If NULL, use `DefaultLayer(sce)`.
#' @param name communication record id. Defaults to `"nichenet"`.
#' @param top_n_ligands Number of top ligands to extract. Defaults to 20.
#' @param exprs_pct Numeric, threshold for percentage of expressing cells to consider a gene expressed. Defaults to 0.10.
#' @return A SingleCellExperiment object with updated CCI state.
#' @importFrom SummarizedExperiment assay
#' @importFrom SummarizedExperiment assayNames
#' @importFrom utils head
#' @export
RunNicheNet <- function(sce, group = NULL,
                        receiver,
                        sender = "all",
                        genes_oi = NULL,
                        ligand_target_matrix,
                        lr_network,
                        weighted_networks,
                        assay_name = NULL,
                        layer = NULL,
                        name = "nichenet",
                        top_n_ligands = 20,
                        exprs_pct = 0.10) {

    if (!requireNamespace("nichenetr", quietly = TRUE)) {
        stop("Package 'nichenetr' is needed for this function. Please install it.")
    }

    if (missing(ligand_target_matrix) || missing(lr_network) || missing(weighted_networks)) {
        stop("NicheNet requires prior models: ligand_target_matrix, lr_network, and weighted_networks. Please provide them.")
    }

    if(! "logcounts" %in% assayNames(sce)){
        sce <- NormalizeData(sce)
    }

    source <- sclet_resolve_expression_source(
        object = sce,
        layer = layer,
        assay = assay_name,
        prefer_nonscaled = TRUE,
        context = "NicheNet analysis"
    )
    assay_name <- source$assay

    meta <- as.data.frame(colData(sce))
    resolved_group <- sclet_resolve_group_column(sce, meta = meta, group = group)
    meta <- resolved_group$meta
    group <- resolved_group$group

    exprs_mat <- assay(sce, assay_name)
    
    # Identify sender and receiver cells
    if (!receiver %in% meta[[group]]) {
        stop("Receiver cluster not found in group column.")
    }
    receiver_cells <- colnames(sce)[meta[[group]] == receiver]
    
    if (length(sender) == 1 && sender == "all") {
        sender_cells <- colnames(sce)
    } else {
        missing_senders <- setdiff(sender, meta[[group]])
        if (length(missing_senders) > 0) {
            stop("Sender clusters not found: ", paste(missing_senders, collapse = ", "))
        }
        sender_cells <- colnames(sce)[meta[[group]] %in% sender]
    }

    # Define expressed genes
    expressed_genes_receiver <- rownames(exprs_mat)[sclet_matrix_rowMeans(exprs_mat[, receiver_cells, drop=FALSE] > 0) >= as.numeric(exprs_pct)]
    expressed_genes_sender <- rownames(exprs_mat)[sclet_matrix_rowMeans(exprs_mat[, sender_cells, drop=FALSE] > 0) >= as.numeric(exprs_pct)]
    background_expressed_genes <- expressed_genes_receiver

    # Define genes of interest
    if (is.null(genes_oi)) {
        if (!requireNamespace("scran", quietly = TRUE)) {
            stop("Package 'scran' is needed to automatically find marker genes. Please install it or provide 'genes_oi'.")
        }
        message("Automatically finding marker genes for receiver: ", receiver)
        markers <- scran::findMarkers(exprs_mat, meta[[group]], pval.type = "any")
        receiver_markers <- markers[[receiver]]
        genes_oi <- rownames(receiver_markers)[receiver_markers$FDR < 0.05 & receiver_markers$summary.logFC > 0.25]
    }
    
    if (length(genes_oi) == 0) {
        stop("No genes of interest found. NicheNet requires a valid gene set of interest.")
    }

    # Extract ligands and receptors
    ligands <- lr_network$from
    receptors <- lr_network$to
    
    expressed_ligands <- intersect(ligands, expressed_genes_sender)
    expressed_receptors <- intersect(receptors, expressed_genes_receiver)
    
    lr_network_expressed <- lr_network[lr_network$from %in% expressed_ligands & lr_network$to %in% expressed_receptors, ]
    potential_ligands <- unique(lr_network_expressed$from)
    
    if (length(potential_ligands) == 0) {
        stop("No potential ligands expressed in sender cells matching receptors in receiver cells.")
    }

    # Predict ligand activities
    message("Predicting ligand activities...")
    ligand_activities <- nichenetr::predict_ligand_activities(
        geneset = genes_oi,
        background_expressed_genes = background_expressed_genes,
        ligand_target_matrix = ligand_target_matrix,
        potential_ligands = potential_ligands
    )
    
    # Sort and get top ligands
    ligand_activities <- ligand_activities[order(ligand_activities$pearson, decreasing = TRUE), ]
    top_ligands <- head(ligand_activities$test_ligand, top_n_ligands)

    # Get ligand-receptor links
    message("Getting weighted ligand-receptor links...")
    lr_links <- nichenetr::get_weighted_ligand_receptor_links(
        best_upstream_ligands = top_ligands,
        expressed_receptors = expressed_receptors,
        lr_network = lr_network,
        weighted_networks_lr_sig = weighted_networks$lr_sig
    )
    
    # Format into standard interactions dataframe
    if (is.null(lr_links) || nrow(lr_links) == 0) {
        warning("NicheNet returned empty ligand-receptor matrix.")
        cci_interactions <- data.frame()
    } else {
        # lr_links usually has format: from, to, weight
        df <- as.data.frame(lr_links)
        
        # Merge with pearson score from activities
        df <- merge(df, ligand_activities[, c("test_ligand", "pearson")], by.x = "from", by.y = "test_ligand", all.x = TRUE)
        
        sender_str <- paste(sender, collapse = ",")
        
        cci_interactions <- data.frame(
            source_group = sender_str,
            target_group = receiver,
            ligand = as.character(df$from),
            receptor = as.character(df$to),
            pathway = NA_character_,
            probability = df$pearson, # Use Pearson correlation as probability/score
            p_value = NA_real_,
            method = "NicheNet",
            stringsAsFactors = FALSE
        )
        
        # Sort by pearson score
        cci_interactions <- cci_interactions[order(cci_interactions$probability, decreasing = TRUE), ]
        rownames(cci_interactions) <- NULL
    }

    # Set state
    sce <- sclet_set_analysis(
        sce,
        "nichenet",
        list(
            id = name,
            method = "NicheNet",
            object = NULL,
            interactions = cci_interactions,
            group = group,
            assay = assay_name,
            layer = source$layer,
            params = list(
                receiver = receiver,
                sender = sender,
                top_n_ligands = top_n_ligands,
                exprs_pct = exprs_pct
            )
        )
    )
    
    sce <- sclet_set_analysis_state(
        object = sce,
        type = "communication",
        id = name,
        method = "NicheNet",
        inputs = list(
            assay = assay_name,
            layer = source$layer,
            group = group
        ),
        artifacts = list(
            analysis_key = "nichenet"
        ),
        params = list(
            receiver = receiver,
            sender = sender,
            top_n_ligands = top_n_ligands
        ),
        summary = list(
            n_interactions = nrow(cci_interactions)
        ),
        active = TRUE
    )
    
    return(sce)
}
