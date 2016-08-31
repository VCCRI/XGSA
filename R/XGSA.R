
#' @title XGSA: Cross Species Gene Set Analysis
#'
#' @description XGSA is a statistical method for cross-species gene set analysis. This package contains the neccessary R code to run XGSA as described in "XGSA: a statistical method for cross-species gene set analysis". XGSA is implemented within a helpful framework to speed up cross-species analyses.
#'
#'
#' @author Djordje Djordjevic <D.djordjevic@victorchang.edu.au>
#'
#' @docType package
#' @name xgsa-package
#' @aliases xgsa XGSA
#'
#' @examples
#' ## Do a sample analyses using cardiac data.
#'
#'
NULL


#
# Copyright Victor Chang Cardiac Research Institute 2016
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
##########################
#
#
# Users should also check the names of functions loaded by this script
#
# Sometimes you might get temporarily disconnected from the BioMart web server with an error message like:
#
# "Error in value[[3L]](cond) :
#  Request to BioMart web service failed. Verify if you are still connected to the internet.  Alternatively the BioMart web service is temporarily down."
#
# If this happens just try to re run the command that failed.
#
# XGSA depends on two libraries, slam and biomaRt
# To install these execute the following commands:
#
# source("http://bioconductor.org/biocLite.R")
# biocLite("biomaRt")
# install.packages("slam")
#
# Processing the Gene Ontology also depends on the packages AnnotationDBI and igraph
# install.packages("AnnotationDBI")
# install.packages("igraph")
#
##########################

require(slam)
require(biomaRt)
require(AnnotationDbi)
require(igraph)
require(GO.db)
require(graph)

##########################
# When Ensembl move servers (rarely) or it is down (less rare) these variable may need to be changed, for example to an archived version

#' @title xgsa_params
#' @rdname xgsa_params
#' @name xgsa_params
#' @export biomartID
#' @export hostID
#' @export supportedSpecies
#'
biomartID <- "ENSEMBL_MART_ENSEMBL"

#hostID <-  "www.ensembl.org"
hostID <-  "asia.ensembl.org"
#hostID <-  "dec2015.archive.ensembl.org"

#' @title find_supported_datasets
#'
#' @description
#' This function determines which species are currently supported by Ensembl
#'
#' @rdname find_supported_datasets
#' @name find_supported_datasets
#'
#' @details
#' This function determines which species are currently supported by Ensembl
#'
#' @return This function returns a vector or organism names in the form 'hsapiens'
#'
#'
#' @importFrom biomaRt useMart listDatasets
#'
#' @export
#'
#' @examples
#' supported_species <- find_supported_datasets()
#' head(supported_species)


find_supported_datasets <- function(default=TRUE){
  ensembl <- useMart(biomartID, dataset='hsapiens_gene_ensembl', host=hostID)
  sets <- listDatasets(ensembl)
  organisms <- gsub("([a-z]*)_.*","\\1",sets[,1],"_")
  return(organisms)
}

##########################
# This is stored globally at the start of the session to avoid having to be recalled
supportedSpecies <- find_supported_datasets()



##########################
# This function executes the paired Fisher's exact tests which form the basis of XGSA
# It returns a single numeric p-value
# Input parameters are the Ensembl gene IDs in the two species, corresponding to the rows and columns of the homology matrix, the homology matrix itself, maximum and minimum gene set sizes to be considered, and the two gene universes to consider.


#' @title paired_fishers_exact_tests
#'
#' @description
#' This function executes the paired Fisher's exact tests which form the basis of XGSA.
#'
#' @rdname paired_fishers_exact_tests
#' @name paired_fishers_exact_tests
#'
#' @details
#' This function executes the paired Fisher's exact tests which form the basis of XGSA.
#' @return This function returns a single numeric p-value
#'
#' @param row_genes Ensembl gene IDs from species A
#' @param col_genes Ensembl gene IDs from species B
#' @param homology_matrix The homology matrix for species A and B
#' @param min Minimum gene set sizes to be considered, default is 5
#' @param max Maximum gene set sizes to be considered, default is 500
#' @param rowuniverse Gene universe in species A
#' @param coluniverse Gene universe in species B
#'
#' @importFrom slam row_sums col_sums
#'
#' @export
#'
#' @examples
#' Used within XGSA framework

paired_fishers_exact_tests<-function(row_genes, col_genes, homology_matrix, min=5, max=500, rowuniverse, coluniverse){

  #only test those genes that actually exist in the matrix
  row_genes <- unique(as.character(row_genes[row_genes%in%rownames(homology_matrix)]))
  col_genes <- unique(as.character(col_genes[col_genes%in%colnames(homology_matrix)]))

  # figure out the correct universe size
  rowuniverse <- unique(as.character(rowuniverse[rowuniverse%in%rownames(homology_matrix)]))
  coluniverse <- unique(as.character(coluniverse[coluniverse%in%colnames(homology_matrix)]))

  rowuniverse <- intersect(rowuniverse, rownames(homology_matrix)[which(row_sums(homology_matrix[,coluniverse, drop=FALSE])>0)])
  coluniverse <- intersect(coluniverse, colnames(homology_matrix)[which(col_sums(homology_matrix[rowuniverse,, drop=FALSE])>0)])


  #check size of gene sets
  if((length(row_genes) < min) | (length(col_genes) < min) | (length(row_genes) > max) | (length(col_genes) > max)) { return(NULL) }

  #share orthology information by predicting gene set genes
  predicted_col_genes <- colnames(homology_matrix)[which(col_sums(homology_matrix[row_genes,, drop=FALSE])>0)]
  predicted_row_genes <- rownames(homology_matrix)[which(row_sums(homology_matrix[,col_genes, drop=FALSE])>0)]

  rowhits <- intersect(row_genes, predicted_row_genes)
  colhits <- intersect(col_genes, predicted_col_genes)

  row.p <- fisher.test(matrix(c(length(rowhits),length(predicted_row_genes)-length(rowhits),length(row_genes)-length(rowhits),length(rowuniverse)-length(unique(union(row_genes, predicted_row_genes)))),2,2), alternative='greater')$p.value
  col.p <- fisher.test(matrix(c(length(colhits),length(predicted_col_genes)-length(colhits),length(col_genes)-length(colhits),length(coluniverse)-length(unique(union(col_genes, predicted_col_genes)))),2,2), alternative='greater')$p.value


  return(max(row.p, col.p))
}



##########################
# This function returns the genes that overlap between the two sets
# It returns a list with two vectors of gene IDs, matching thw row and column species of the homology table

#' @title get_overlap_genes
#'
#' @description
#' This function returns the genes that overlap between the two sets.
#'
#' @rdname get_overlap_genes
#' @name get_overlap_genes
#'
#' @details
#' This function returns the genes that overlap between the two sets.
#' @return This function returns a list with two vectors of gene IDs, matching the row and column species of the homology table.
#'
#' @param row_genes Ensembl gene IDs from species A
#' @param col_genes Ensembl gene IDs from species B
#' @param homology_matrix The homology matrix for species A and B
#' @param min Minimum gene set sizes to be considered, default is 5
#' @param max Maximum gene set sizes to be considered, default is 500
#'
#' @importFrom slam row_sums col_sums
#'
#' @export
#'
#' @examples
#' Used within XGSA framework
get_overlap_genes<-function(row_genes, col_genes, homology_matrix, min=5, max=500){

  #only test those genes that actually exist in the matrix
  row_genes <- unique(as.character(row_genes[row_genes%in%rownames(homology_matrix)]))
  col_genes <- unique(as.character(col_genes[col_genes%in%colnames(homology_matrix)]))

  #check size of gene sets
  #if((length(row_genes) < min) | (length(col_genes) < min) | (length(row_genes) > max) | (length(col_genes) > max)) { return(NULL) }

  #share orthology information by predicting gene set genes
  predicted_col_genes <- colnames(homology_matrix)[which(col_sums(homology_matrix[row_genes,, drop=FALSE])>0)]
  predicted_row_genes <- rownames(homology_matrix)[which(row_sums(homology_matrix[,col_genes, drop=FALSE])>0)]

  rowhits <- intersect(row_genes, predicted_row_genes)
  colhits <- intersect(col_genes, predicted_col_genes)

  overlap <- list(rowhits, colhits)
  names(overlap) <- c("rowhits","colhits")
  return(overlap)
}



##########################
# This function generates a sparse matrix from a table with two or three columns: ID ID Value(optional)
# It uses the simple_triplet_matrix structure from the package 'slam'


#' @title generate_sparse_matrix
#'
#' @description
#' This function generates a sparse matrix from a table with two or three columns: ID ID Value(optional).
#' It uses the simple_triplet_matrix structure from the package 'slam'.
#'
#' @rdname generate_sparse_matrix
#' @name generate_sparse_matrix
#'
#' @details
#' This function generates a sparse matrix from a table with two or three columns: ID ID Value(optional).
#' It uses the simple_triplet_matrix structure from the package 'slam'.
#' @return This function return a sparse matrix
#'
#' @param homology_table Table of homology mapping between two species, as retrieved by the function get_homology_matrix()
#'
#' @importFrom slam simple_triplet_matrix
#'
#' @export
#'
#' @examples
#' Not used trivially
#'

generate_sparse_matrix <- function(homology_table){
  #require(slam)
  # Compare only genes with homology mapping
  homology_table <- unique(homology_table)[!unique(homology_table)[,2]=="",]

  if(ncol(homology_table)==3){
    #hms <- sparseMatrix(i = match(homology_table[,1], unique(homology_table[,1]) ), j = match(homology_table[,2], unique(homology_table[,2])), x = homology_table[,3])
    hms <- simple_triplet_matrix(i = match(homology_table[,1], unique(homology_table[,1]) ), j = match(homology_table[,2], unique(homology_table[,2])), v = homology_table[,3])
  } else if(ncol(homology_table)==2){
    hms <- simple_triplet_matrix(i = match(homology_table[,1], unique(homology_table[,1]) ), j = match(homology_table[,2], unique(homology_table[,2])), v = rep(1,nrow(homology_table)))
    #hms <- sparseMatrix(i = match(homology_table[,1], unique(homology_table[,1]) ), j = match(homology_table[,2], unique(homology_table[,2])))
  } else {
    print("ERROR: incorrect number of columns in homology table ")
    return(NULL)
  }
  rownames(hms) = unique(homology_table[,1])
  colnames(hms) = unique(homology_table[,2])
  return(hms)
}
##########################


##########################
# generates the best reciprocal hits matrix
# note this can still give not perfect 1-1 mapping if homologs have identical sequence conservation in each direction


#' @title generate_best_reciprocal_hits_matrix
#'
#' @description
#' This function generates a homology matrix based on the best reciprocal hits method.
#' Note: this can still give not perfect 1-1 mapping if homologs have identical sequence conservation in each direction
#'
#' @rdname generate_best_reciprocal_hits_matrix
#' @name generate_best_reciprocal_hits_matrix
#'
#' @details
#' This function generates a homology matrix based on the best reciprocal hits method.
#' Note: this can still give not perfect 1-1 mapping if homologs have identical sequence conservation in each direction
#' @return This function returns a sparse matrix
#'
#' @param hms1 Homology matrix with homology scores with species 1 as rows
#' @param hms2 Homology matrix with homology scores with species 2 as rows
#' @param homology_table Table of homology mapping between two species (species 1 in column 1 and species 2 in column 2)
#'
#' @importFrom slam simple_triplet_matrix
#'
#' @export
#'
#' @examples
#' Not used trivially
#'
generate_best_reciprocal_hits_matrix <- function(hms1, hms2, homology_table){

  # Need To Add: if there are no values in the matrices return failure
  # Compare only genes with homology mapping
  homology_table <- unique(homology_table)[!unique(homology_table)[,2]=="",]

  reciprocal.include <- vector()
  for(i in 1:nrow(homology_table)){
    h <- homology_table[i,1]
    z <- homology_table[i,2]
    if ((z %in% names(which(hms1[h,]==max(hms1[h,]))))&(h %in% names(which(hms2[z,]==max(hms2[z,]))))){
      reciprocal.include <- append(reciprocal.include, 1)
    } else {
      reciprocal.include <- append(reciprocal.include, 0)
    }
  }

  sub_homology_table <- homology_table[reciprocal.include==1,]

  hms_BRH <- generate_sparse_matrix(sub_homology_table)
  return(hms_BRH)
}
###########################

##########################
# use ENSEMBL to generate homology matrix
# input parameteres are strings like 'hsapiens' and 'drerio' and a number indicating which species to get sequence identity wrt.

#' @title get_homology_table
#'
#' @description
#' This function uses information from Ensembl to generate a homology table between two species.
#'
#' @rdname get_homology_table
#' @name get_homology_table
#'
#' @details
#' This function uses information from Ensembl to generate a homology table between two species.
#' @return This function return a data frame with 2 or 3 columns, representing the Ensembl gene IDs in species 1 and 2, and the sequence identity score if requested.
#'
#' @param species1 Species 1 name, in the form 'hsapiens'
#' @param species2 Species 2 name, in the form 'hsapiens'
#' @param sequence_identity_reference Flag of which sequence identity value to retrieve. 0 (default) returns no value); 1 returns the sequence of identity of genes in species 2 compared to species 1; 2 returns the sequence of identity of genes in species 1 compared to species 2.
#'
#' @importFrom biomaRt useMart getBM
#'
#' @export
#'
#' @examples
#' human_zebrafish_homology_table <- get_homology_table('hsapiens', 'drerio')
#' head(human_zebrafish_homology_table)
#'

get_homology_table <- function(species1, species2, sequence_identity_reference = 0){
  #require(biomaRt)
  dataset_name <- paste(species1, "_gene_ensembl", sep="")
  homolog_attribute <- paste(species2, "_homolog_ensembl_gene", sep="")
  ensembl <- useMart(biomartID, dataset=dataset_name, host=hostID)

  if(sequence_identity_reference == 1){
    homology_table <- getBM(attributes=c("ensembl_gene_id", homolog_attribute, paste(species2, "_homolog_perc_id", sep="")), mart = ensembl)
  } else if (sequence_identity_reference == 2){
    homology_table <- getBM(attributes=c("ensembl_gene_id", homolog_attribute, paste(species2, "_homolog_perc_id_r1", sep="")), mart = ensembl)
  } else {
    homology_table <- getBM(attributes=c("ensembl_gene_id", homolog_attribute), mart = ensembl)
  }
  return(homology_table)
}

##########################
# Returns a data frame with one column containing all Ensembl gene IDs for the given species
# Input parameter is the species string (eg. 'hsapiens')

#' @title get_ENSEMBL_gene_list
#'
#' @description
#' This helper function returns a data frame with one column containing all Ensembl gene IDs for the given species
#'
#' @rdname get_ENSEMBL_gene_list
#' @name get_ENSEMBL_gene_list
#'
#' @details This helper function returns a data frame with one column containing all Ensembl gene IDs for the given species
#' @return This helper function returns a data frame with one column containing all Ensembl gene IDs for the given species
#'
#' @param species Species name in the form 'hsapiens'
#'
#' @importFrom biomaRt useMart getBM
#'
#' @export
#'
#' @examples
#' human_ensembl_IDs <- get_ENSEMBL_gene_list('hsapiens')
#' head(human_ensembl_IDs)
#'

get_ENSEMBL_gene_list <- function(species){
  require(biomaRt)
  dataset_name <- paste(species, "_gene_ensembl", sep="")
  ensembl <- useMart(biomartID, dataset=dataset_name, host=hostID)
  gene.list <- getBM(attributes=c("ensembl_gene_id"), mart = ensembl)
  return(gene.list)
}

##########################
# Useful to convert from transcript IDs to gene IDs
# Returns a data frame with 2 columns
# Input parameter is the species string (eg. 'hsapiens')

#' @title get_ENSEMBL_trasncript_table
#'
#' @description
#' This helper function returns a data frame with Ensembl gene IDs in the first column and Ensembl transcript IDs in the second column. This can be useful when gene expression data has been mapped to transcripts for example.
#'
#' @rdname get_ENSEMBL_trasncript_table
#' @name get_ENSEMBL_trasncript_table
#'
#' @details This helper function returns a data frame with Ensembl gene IDs in the first column and Ensembl transcript IDs in the second column. This can be useful when gene expression data has been mapped to transcripts for example.
#' @return This helper function returns a data frame with Ensembl gene IDs in the first column and Ensembl transcript IDs in the second column. This can be useful when gene expression data has been mapped to transcripts for example.
#'
#' @param species Species name in the form 'hsapiens'
#'
#' @importFrom biomaRt useMart getBM
#'
#' @export
#'
#' @examples
#' transcript_table <- get_ENSEMBL_trasncript_table('hsapiens')
#' head(transcript_table)
#'

get_ENSEMBL_trasncript_table <- function(species){
  require(biomaRt)
  dataset_name <- paste(species, "_gene_ensembl", sep="")
  ensembl <- useMart(biomartID, dataset=dataset_name, host=hostID)
  transcript.table <- getBM(attributes=c("ensembl_gene_id","ensembl_transcript_id"), mart = ensembl)
  return(transcript.table)
}

##########################
# Retrieves the mapping  between ENSEMBL ids and gene symbols
# Returns a data frame with 2 columns
# Input parameter is the species string (eg. 'hsapiens')

#' @title get_ENSEMBL_symbol_map
#'
#' @description
#' This helper function retrieves the mapping  between ENSEMBL gene IDs and common gene symbols
#'
#' @rdname get_ENSEMBL_symbol_map
#' @name get_ENSEMBL_symbol_map
#'
#' @details
#' This helper function retrieves the mapping  between ENSEMBL gene IDs and common gene symbols
#' @return This function returns a data frame with Ensembl gene IDs in the first column and common gene symbols in the second column.
#'
#' @param species Species name in the form 'hsapiens'
#'
#' @importFrom biomaRt useMart getBM
#'
#' @export
#'
#' @examples
#' human_symbol_map <- get_ENSEMBL_symbol_map('hsapiens')
#' head(human_symbol_map)

get_ENSEMBL_symbol_map <- function(species){
  require(biomaRt)
  dataset_name <- paste(species, "_gene_ensembl", sep="")
  ensembl <- useMart(biomartID, dataset=dataset_name, host=hostID)
  ENSEMBL_symbol_map <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"), mart = ensembl)
  return(ENSEMBL_symbol_map)
}

##########################
# Retrieves the GO term names and IDs of the latest Gene Ontology annotations from ENSEMBL
# Returns a data frame containing GO terms and term names
# Input parameter is the species string (eg. 'hsapiens')

#' @title get_GO_names
#'
#' @description
#' This function retrieves the GO term names and IDs of the latest Gene Ontology annotations from Ensembl
#'
#' @rdname get_GO_names
#' @name get_GO_names
#'
#' @details
#' This function retrieves the GO term names and IDs of the latest Gene Ontology annotations from Ensembl
#' @return This function returns a data frame containing GO terms in column 1 and term names in column 2
#'
#' @param species Species name in the form 'hsapiens'
#'
#' @importFrom biomaRt useMart getBM
#'
#' @export
#'
#' @examples
#' go_names <- get_GO_names('hsapiens')
#' head(go_names)

get_GO_names <- function(species){
  require(biomaRt)
  dataset_name <- paste(species, "_gene_ensembl", sep="")
  ensembl <- useMart(biomartID, dataset=dataset_name, host=hostID)
  GOTerms <- getBM(mart=ensembl, attributes=c('go_id','name_1006'))
  return(GOTerms)
}

##########################
# Retrieves the latest Gene Ontology annotations from ENSEMBL
# Returns a data frame containing gene ids, associated go terms, GO evidence codes, and GO ontology type (BP, MF, CC)
# Input parameter is the species string (eg. 'hsapiens')

#' @title get_GO_mappings
#'
#' @description
#' This function retrieves the latest Gene Ontology annotations from Ensembl
#'
#' @rdname get_GO_mappings
#' @name get_GO_mappings
#'
#' @details
#' This function retrieves the latest Gene Ontology annotations from Ensembl
#' @return This function returns a data frame containing Ensembl gene IDs, annotated go terms, GO evidence codes and GO ontology type (BP, MF, CC), in columns 1:4 respectively.
#'
#' @param species Species name in the form 'hsapiens'
#'
#' @importFrom biomaRt useMart getBM
#'
#' @export
#'
#' @examples
#' human_GO_mappings <- get_GO_mappings('hsapiens')
#' head(human_GO_mappings)

get_GO_mappings <- function(species){
  require(biomaRt)
  dataset_name <- paste(species, "_gene_ensembl", sep="")
  ensembl <- useMart(biomartID, dataset=dataset_name, host=hostID)
  EG2GO <- getBM(mart=ensembl, attributes=c('ensembl_gene_id','go_id','go_linkage_type','namespace_1003'))
  return(EG2GO)
}

##########################
# Retrieves GO mappings with only the supplied evidence codes and from only the supplied ontologies and converts the result into a list
# Returns a list of GO terms with constituent ENESMBL gene IDS
# Input parameters are the species string (eg. 'hsapiens'), a character vector of evidence codes and a character vector of ontology codes (eg. 'BP')
# Defaults to non-computational evidence

#' @title get_GO_list_from_ontologies_with_evidence_codes
#'
#' @description
#' This function retrieves GO mappings with only the supplied evidence codes and from only the supplied ontologies and converts the result into a list.
#'
#' @rdname get_GO_list_from_ontologies_with_evidence_codes
#' @name get_GO_list_from_ontologies_with_evidence_codes
#'
#' @details
#' This function retrieves GO mappings with only the supplied evidence codes and from only the supplied ontologies and converts the result into a list.
#' @return This function returns a named list of GO terms with annotated Ensembl gene IDs within each element.
#'
#' @param species Species name in the form 'hsapiens'
#' @param evidence.codes A character vector of evidence codes, defaults to 'direct' evidence = c("EXP","IDA","IPI","IMP","IGI","IEP", "TAS", "IC"). A description of GO evidence codes can be found here: http://geneontology.org/page/guide-go-evidence-codes
#' @param ontologies A character vector of the desired ontologies to query, defaults to all three of them - i.e. c("biological_process", "molecular_function", "cellular_component")
#'
#' @export
#'
#' @examples
#' Human_trim_GO <- get_GO_list_from_ontologies_with_evidence_codes('hsapiens')
#' summary(Human_trim_GO[1:5])

get_GO_list_from_ontologies_with_evidence_codes <- function(species, evidence.codes=  c("EXP","IDA","IPI","IMP","IGI","IEP", "TAS", "IC"), ontologies = c("biological_process", "molecular_function", "cellular_component")){
  GOMap <- get_GO_mappings(species)
  subGOMap <- GOMap[(GOMap$go_linkage_type %in% evidence.codes)&(GOMap$namespace_1003 %in% ontologies),]
  GOTerms <- split(subGOMap$ensembl_gene_id, subGOMap$go_id)
  return(GOTerms)
}


##########################
# Collapses all child GO term annotations to their parent terms, resulting in the full redundant GO
# Returns a list of GO terms with constituent gene IDS
# Input parameter is the list of GO term annotations produced by "get_GO_list_from_ontologies_with_evidence_codes" for example

#' @title collapseGO
#'
#' @description
#' This function collapses all child GO term annotations to their parent terms, resulting in the full redundant GO lists.
#'
#' @rdname collapseGO
#' @name collapseGO
#'
#' @details
#' This function collapses all child GO term annotations to their parent terms, resulting in the full redundant GO lists.
#' @return This function returns a named list of GO terms with annotated Ensembl gene IDs within each element.
#'
#' @param GO_Terms A non-redundant list of GO term annotations produced by "get_GO_list_from_ontologies_with_evidence_codes" for example.
#'
#' @import igraph
#' @import graph
#' @import GO.db
#' @import AnnotationDbi
#'
#' @export
#'
#' @examples
#' Human_trim_GO <- get_GO_list_from_ontologies_with_evidence_codes('hsapiens')
#' summary(Human_trim_GO[1:5])
#'
#' Human_full_GO <- collapseGO(Human_trim_GO)
#' summary(Human_full_GO[1:5])

collapseGO <- function(GO_Terms){

  suppressMessages(require(igraph))
  require(AnnotationDbi)

  bp <- AnnotationDbi::makeGOGraph(ont = "bp")
  mf <- AnnotationDbi::makeGOGraph(ont = "mf")
  cc <- AnnotationDbi::makeGOGraph(ont = "cc")

  go.bp <- igraph.from.graphNEL(bp)
  go.mf <- igraph.from.graphNEL(mf)
  go.cc <- igraph.from.graphNEL(cc)

  go.g <- igraph::union(go.bp, go.mf, go.cc)

  new_GO_Terms <- lapply(names(GO_Terms), function(G){
    if(!G%in%names(V(go.g))){
      return(unique(na.omit(unlist(GO_Terms[c(G)]))))
    }
    child_terms <- names(subcomponent(go.g, G, mode = "in"))
    return(unique(na.omit(unlist(GO_Terms[c(child_terms,G)]))))
  })

  names(new_GO_Terms) <- names(GO_Terms)
  return(new_GO_Terms)
}


##########################
# Retrieves GO mappings with only the supplied evidence codes and from only the supplied ontologies and converts the result into a list
# Collapses all child term annotations to their parent terms, resulting in the full redundant GO
# Returns a list of GO terms with constituent ENESMBL gene IDS
# Input parameters are the species string (eg. 'hsapiens'), a character vector of evidence codes and a character vector of ontology codes (eg. 'BP')
# Defaults to non-computational evidence

#' @title get_GO
#'
#' @description
#' This function retrieves GO mappings with only the supplied evidence codes and from only the supplied ontologies and converts the result into a list. It retrieves the full and redundant GO mappings.
#'
#' @rdname get_GO
#' @name get_GO
#'
#' @details
#' This function retrieves GO mappings with only the supplied evidence codes and from only the supplied ontologies and converts the result into a list. It retrieves the full and redundant GO mappings.
#' @return This function returns a named list of GO terms with annotated Ensembl gene IDs within each element.
#'
#' @param species Species name in the form 'hsapiens'
#' @param evidence.codes A character vector of evidence codes, defaults to 'direct' evidence = c("EXP","IDA","IPI","IMP","IGI","IEP", "TAS", "IC"). A description of GO evidence codes can be found here: http://geneontology.org/page/guide-go-evidence-codes
#' @param ontologies A character vector of the desired ontologies to query, defaults to all three of them - i.e. c("biological_process", "molecular_function", "cellular_component")
#'
#' @export
#'
#' @examples
#' Human_full_GO <- get_GO('hsapiens')

get_GO <- function(species, evidence.codes=  c("EXP","IDA","IPI","IMP","IGI","IEP", "TAS", "IC"), ontologies = c("biological_process", "molecular_function", "cellular_component")){
  go <- get_GO_list_from_ontologies_with_evidence_codes(species, evidence.codes, ontologies)
  print("retrieved GO")
  c.go <- collapseGO(go)
  return(c.go)
}


##########################
# Returns a list with two elements (one for each species).
# Each element contains a vector of gene IDs in the species that have more than 1 homolog in the other species
# Input parameters are the two species strings (eg. 'hsapiens', 'drerio')

#' @title get_complex_genes
#'
#' @description
#' This helper function retrieves complex genes between two species.
#'
#' @rdname get_complex_genes
#' @name get_complex_genes
#'
#' @details
#' This helper function retrieves complex genes between two species.
#' @return This function returns a list with two elements (one for each species), each containing a vector of gene IDs in that species that have more than 1 homolog in the other species
#'
#' @param species1 Species 1 name in the form 'hsapiens'
#' @param species2 Species 2 name in the form 'hsapiens'
#'
#' @export
#'
#' @examples
#' human_zebrafish_complex_genes <- get_complex_genes('hsapiens','drerio')

get_complex_genes <- function(species1, species2){
  hm <- get_homology_matrix(species1, species2)
  hm.boolean <- hm > 0
  complex.cols <- which(col_sums(hm.boolean) > 1)
  complex.rows <- which(row_sums(hm.boolean) > 1)
  species1.complex.genes <- rownames(hm.boolean)[complex.rows]
  species2.complex.genes <- colnames(hm.boolean)[complex.cols]
  complex.genes <- list(species1.complex.genes, species2.complex.genes)
  names(complex.genes) <- c(species1,species2)
  return(complex.genes)
}

##########################
# Calculates the complexity score (Chi) of the gene sets in an XGSA data set with respect to a second species
# Returns a list of complexity scores, one score for each gene set
# Input parameters are the XGSA data set object and the other species string (eg. 'drerio')

#' @title calculate_genesetlist_complexity
#'
#' @description
#' This helper function calculates the complexity score of the gene sets in an XGSA data set with respect to a second species. The complexity score is defined in the publication, as the fraction of gene set that has complex homology to another species.
#'
#' @rdname calculate_genesetlist_complexity
#' @name calculate_genesetlist_complexity
#'
#' @details
#' This function calculates the complexity score of the gene sets in an XGSA data set with respect to a second species. The complexity score is defined in the publication, as the fraction of gene set that has complex homology to another species.
#' @return This function returns a list of complexity scores, one score for each gene set
#'
#' @param dataset An XGSA data set created by 'new_XGSA_dataset'
#' @param species2 Species name in the form 'hsapiens', naming the species to compare the dataset against for gene set complexity.
#'
#' @export
#'
#' @examples
#' test_mouse_genes <- c("ENSMUSG00000033837", "ENSMUSG00000031965", "ENSMUSG00000053110", "ENSMUSG00000030557", "ENSMUSG00000051159", "ENSMUSG00000040289", "ENSMUSG00000035458", "ENSMUSG00000028780", "ENSMUSG00000062327", "ENSMUSG00000037868")
#' test_mouse_dataset <- new_XGSA_dataset(species = 'mmusculus', data = list(test_genes = test_mouse_genes), name = 'Test Mouse Genes')
#' mouse_zebrafish_complexity <- calculate_genesetlist_complexity(test_mouse_dataset,'drerio')
#' print(mouse_zebrafish_complexity)

calculate_genesetlist_complexity <- function(dataset, species2){
  #check species aren't the same
  species1 <- dataset$species
  if( species1 == species2 ){
    print("ERROR: no complex homology between the same species")
    return(NULL)
  }
  complex.genes <- get_complex_genes(species1, species2)
  data1 <- dataset$data
  data1.complexity <- lapply(data1, function(X){
    return(length(intersect(X, complex.genes[[species1]])) / length(X))
  })
  return(data1.complexity)
}

##########################
# Calculates an alternative complexity score (not described in the paper) of the gene sets in an XGSA data set with respect to a second species
# Returns a list of complexity scores, one score for each gene set
# Input parameters are the XGSA data set object and the other species string (eg. 'drerio')

#' @title calculate_alternate_genesetlist_complexity
#'
#' @description
#' This function calculates an alternative complexity score (not described in the paper) of the gene sets in an XGSA data set with respect to a second species. The alternate complexity score is the total number of homology edges between a gene set and another species divides by the total number of genes in a gene set. As such it ranges from 1 to infinity.
#'
#' @rdname calculate_alternate_genesetlist_complexity
#' @name calculate_alternate_genesetlist_complexity
#'
#' @details
#' This function calculates an alternative complexity score (not described in the paper) of the gene sets in an XGSA data set with respect to a second species. The alternate complexity score is the total number of homology edges between a gene set and another species divides by the total number of genes in a gene set. As such it ranges from 1 to infinity.
#' @return This function returns a list of complexity scores, one score for each gene set
#'
#' @param dataset An XGSA data set created by 'new_XGSA_dataset'
#' @param species2 Species name in the form 'hsapiens', naming the species to compare the dataset against for gene set complexity.
#'
#' @export
#'
#' @examples
#' test_mouse_genes <- c("ENSMUSG00000033837", "ENSMUSG00000031965", "ENSMUSG00000053110", "ENSMUSG00000030557", "ENSMUSG00000051159", "ENSMUSG00000040289", "ENSMUSG00000035458", "ENSMUSG00000028780", "ENSMUSG00000062327", "ENSMUSG00000037868")
#' test_mouse_dataset <- new_XGSA_dataset(species = 'mmusculus', data = list(test_genes = test_mouse_genes), name = 'Test Mouse Genes')
#' mouse_zebrafish_complexity <- calculate_alternate_genesetlist_complexity(test_mouse_dataset,'drerio')
#' print(mouse_zebrafish_complexity)

calculate_alternate_genesetlist_complexity <- function(dataset, species2){
  #check species aren't the same
  species1 <- dataset$species
  if( species1 == species2 ){
    print("ERROR: no complex homology between the same species")
    return(NULL)
  }
  hm <- get_homology_matrix(species1, species2)
  data1 <- dataset$data
  data1.complexity <- lapply(data1, function(X){
    return(sum(row_sums(hm[X[X%in%rownames(hm)],]))/sum(X%in%rownames(hm)))
  })
  return(data1.complexity)
}

##########################
# Creates a homology matrix between two species by retrieving the homology table and then converting it into a sparse matrix
# Returns a sparse matrix
# Input parameters are the two species strings (eg. 'hsapiens', 'drerio') and a boolean of whether to return the sequence identity values or not

#' @title get_homology_matrix
#'
#' @description
#' This function creates a homology matrix between two species by retrieving the homology table and then converting it into a sparse matrix.
#'
#' @rdname get_homology_matrix
#' @name get_homology_matrix
#'
#' @details
#' This function creates a homology matrix between two species by retrieving the homology table and then converting it into a sparse matrix.
#' @return This function returns a sparse matrix where rows represent the genes in species 1 and columns the genes in species 2.
#'
#' @param species1 Species 1 name, in the form 'hsapiens'
#' @param species2 Species 2 name, in the form 'hsapiens'
#' @param seq.identity Boolean of whether or not  to retrieve sequence identity information, defaults to FALSE
#'
#' @export
#'
#' @examples
#' human_zebrafish_matix <- get_homology_matrix('hsapiens','drerio')
#' print(data.matrix(human_zebrafish_matix[1:5,1:5]))

get_homology_matrix <- function(species1, species2, seq.identity = FALSE){
  if(seq.identity == FALSE){
    if(!exists("homology.matrix.list", envir = .GlobalEnv)){
      assign("homology.matrix.list", list(), envir = .GlobalEnv)
    }
    direct <- !is.null(homology.matrix.list[[species1]][[species2]])
    indirect <- !is.null(homology.matrix.list[[species2]][[species1]])
    if(direct){
      # homology matrix already genearted and stored
      return(homology.matrix.list[[species1]][[species2]])
    } else if(indirect){
      # inverse homology matrix already generated and stored
      return(t(homology.matrix.list[[species2]][[species1]]))
    } else {
      if(species1 == species2){

        gene.list <- get_ENSEMBL_gene_list(species1)
        hm <- generate_sparse_matrix(cbind(gene.list, gene.list))

        #DO a same species analysis

      } else {

        # create a new homology matrix and return it
        ht <- get_homology_table(species1, species2)
        hm <- generate_sparse_matrix(ht)
      }
      if(species1 %in% names(homology.matrix.list)){
        homology.matrix.list[[species1]][[species2]] <<- hm
      } else {
        homology.matrix.list[[species1]] <<- list()
        homology.matrix.list[[species1]][[species2]] <<- hm
      }
      return(homology.matrix.list[[species1]][[species2]])
    }
  } else {
    if(!exists("homology.seqid.matrix.list", envir = .GlobalEnv)){
      assign("homology.seqid.matrix.list", list(), envir = .GlobalEnv)
    }
    direct <- !is.null(homology.seqid.matrix.list[[species1]][[species2]])
    if(direct){
      # homology matrix already genearted and stored
      return(homology.seqid.matrix.list[[species1]][[species2]])
    } else {
      if(species1 == species2){
        gene.list <- get_ENSEMBL_gene_list(species1)
        hm <- generate_sparse_matrix(cbind(gene.list, gene.list))
      } else {
        # create a new homology matrix and return it
        ht <- get_homology_table(species1, species2, 1)
        hm <- generate_sparse_matrix(ht)
      }
      if(species1 %in% names(homology.seqid.matrix.list)){
        homology.seqid.matrix.list[[species1]][[species2]] <<- hm
      } else {
        homology.seqid.matrix.list[[species1]] <<- list()
        homology.seqid.matrix.list[[species1]][[species2]] <<- hm
      }
      return(homology.seqid.matrix.list[[species1]][[species2]])
    }
  }
}



##########################
# Creates an XGSA data set object
# Input parameters are the species string (eg. 'hsapiens'), a list of named gene sets containing ENSEMBL IDs, the data type (currently only 'genesetlist' is supported), the name of the dataset (eg. 'HumanGeneOntology') and a character vector containing the ENSEMBL gene ID universe for this data set.
# Returns a list with 5 elements describing the data set

#' @title new_XGSA_dataset
#'
#' @description
#' This function creates an XGSA data set object that is required for performing tests.
#'
#'
#' @rdname new_XGSA_dataset
#' @name new_XGSA_dataset
#'
#' @details
#' This function creates an XGSA data set object that is required for performing tests.
#'
#' @return This function returns a list with 5 elements describing the data set, the species, gene sets, data type, name of data set and the applicable gene universe.
#'
#' @param species Species name, in the form 'hsapiens'
#' @param data A named list of gene sets where each element contains a character vector of Ensembl gene IDs. This can be a named list of only a single element containing gene IDs.
#' @param type The type of data provided, currently only 'genesetlist' is supported.
#' @param name A name for the data set / collection of gene sets
#' @param universe The gene universe for the gene sets in data. They must all have the same gene universe. If no gene universe is provided (default) all Ensembl gene IDs for the indicated species will be used as the gene universe.
#'
#' @export
#'
#' @examples
#' test_mouse_genes <- c("ENSMUSG00000033837", "ENSMUSG00000031965", "ENSMUSG00000053110", "ENSMUSG00000030557", "ENSMUSG00000051159", "ENSMUSG00000040289", "ENSMUSG00000035458", "ENSMUSG00000028780", "ENSMUSG00000062327", "ENSMUSG00000037868")
#' test_mouse_dataset <- new_XGSA_dataset(species = 'mmusculus', data = list(test_genes = test_mouse_genes), name = 'Test Mouse Genes')
#' print(summary(test_mouse_dataset))

new_XGSA_dataset <- function(species, data, type = "genesetlist", name, universe="empty"){
  #check species
  if(!species %in% supportedSpecies){
    print("ERROR: Species not supported")
    return(NULL)
  }

  #check type
  if(!type%in%c("geneset","genesetlist")){
    print("ERROR: Type must be geneset or genesetlist")
    return(NULL)
  }

  #check data
  if(type=="geneset"){
    if(!sum(is(data)[1:2] == c("character","vector") ) == 2){
      print("ERROR: Data is in wrong format")
      return(NULL)
    }
  }  else if(type=="genesetlist"){
    if(!(sum(summary(data)[,"Mode"]=="character") == length(data))){
      print("ERROR: Data is in wrong format")
      return(NULL)
    }
  }
  if(universe[1] == "empty"){
    print("Retrieving gene universe...")
    universe <- unlist(get_ENSEMBL_gene_list(species))
  }

  dataset <- list(species, data, type, name, universe)
  names(dataset) <- c("species", "data", "type", "name", "universe")
  return(dataset)
}

##########################
# Checks that a data set is valid

#' @title check_XGSA_dataset
#'
#' @description
#' This function checks that an XGSA data set is valid
#'
#' @rdname check_XGSA_dataset
#' @name check_XGSA_dataset
#'
#' @details
#' This function checks that an XGSA data set is valid
#'
#' @return This function returns a boolean.
#'
#' @param dataset An XGSA data set
#'
#' @export
#'
#' @examples
#' test_mouse_genes <- c("ENSMUSG00000033837", "ENSMUSG00000031965", "ENSMUSG00000053110", "ENSMUSG00000030557", "ENSMUSG00000051159", "ENSMUSG00000040289", "ENSMUSG00000035458", "ENSMUSG00000028780", "ENSMUSG00000062327", "ENSMUSG00000037868")
#' test_mouse_dataset <- new_XGSA_dataset(species = 'mmusculus', data = list(test_genes = test_mouse_genes), name = 'Test Mouse Genes')
#' check_XGSA_dataset(test_mouse_dataset)


check_XGSA_dataset <- function(dataset){
  spec.check <- is.character(dataset$species)
  data.check <- is.character(dataset$data) | is.list(dataset$data)
  type.check <- dataset$type %in% c("geneset","genesetlist")
  if(spec.check & data.check & type.check){
    return(TRUE)
  } else {
    print("Data is wrong")
    return(FALSE)
  }
}

##########################
# Runs an XGSA test on two data sets
# Input parameters are the two data sets, the test type (currently only "Fisher" is supported) and the minimum and maximum gene set sizes to consider.
# Returns a complex structure with P values and gene overlaps.

#' @title run_XGSA_test
#'
#' @description
#' This function runs an XGSA test on two XGSA data sets.
#'
#' @rdname run_XGSA_test
#' @name run_XGSA_test
#'
#' @details
#' This function runs an XGSA test on two XGSA data sets.
#'
#' @return
#' This function returns a named list with one elements for each comparison. Each list element contains another list with two named elements; 'pvals' - the p-value; 'genes' - the overlapping genes from each species.
#'
#' @param dataset1 A XGSA data set
#' @param dataset2 A XGSA data set
#' @param test The type of statistical test to perform, currently only 'Fisher' (default) is supported
#' @param min Minimum gene set sizes to be considered, default is 5
#' @param max Maximum gene set sizes to be considered, default is 500
#'
#' @export
#'
#' @examples
#' ##########################
#' # Step 1) Create sample mouse gene set.
#' #####
#' # In this example we will extract a list of mouse cardiac development genes from a previous published study.
#' cardiac.perturbation.data <- read.table("http://cardiaccode.victorchang.edu.au/data/Cardiaccode_in_vivo_evidence_2014_06_10.txt", sep=" ", header=TRUE, quote=" ")
#' mouse.cardiac.genes <- unique(unlist(cardiac.perturbation.data[cardiac.perturbation.data$Species == "Mus musculus", c("Regulator","Target")]))
#' # XGSA harnesses the Ensembl homology pipeline, and so we need to convert all of our data sets to Ensembl gene IDs
#' # We will convert the mouse gene symbols to Ensembl IDs using the XGSA helper function "get_ENSEMBL_symbol_map".
#' mouse.ensembl.symbol.map <- get_ENSEMBL_symbol_map(species = 'mmusculus')
#' mouse.cardiac.ensembl.symbols <- mouse.ensembl.symbol.map$ensembl_gene_id[mouse.ensembl.symbol.map$external_gene_name %in% mouse.cardiac.genes]
#'
#' # Now we have a list of Ensembl IDs for mouse cardiac genes, we will turn it into an XGSA data set.
#' # Note that the input data MUST be a named list, which allows for multiple gene sets in the same data set.
#' # As we don't have a defined gene universe for this data set we will use all mouse Ensembl IDs that have an external gene symbol.
#' mouse.data <- new_XGSA_dataset(species = 'mmusculus', data = list(mouseCardiacGenes = mouse.cardiac.ensembl.symbols), type = 'genesetlist', name = 'MouseCardiacGenes', universe = unique(mouse.ensembl.symbol.map$ensembl_gene_id))
#'
#' ##########################
#' # Step 2) Create reference zebrafish GO set.
#' #####
#' # In this example we will compare to the zebrafish Gene Ontology using "direct" evidence only - this means the annotations are NOT transferred between species.
#' # We will use another XGSA helper function to retrieve the latest Gene Ontology information from Ensembl "get_GO_list_from_ontologies_with_evidence_codes".
#' # The gene universe we will use is all ofthe zebrafish biological process genes that we are testing.
#' zebrafish.GO <- get_GO('drerio', ontologies = "biological_process")
#' zebrafish.GO <- zebrafish.GO[lapply(zebrafish.GO, length) > 10 & lapply(zebrafish.GO, length) < 500]
#' zebrafish.GO.data <- new_XGSA_dataset(species = "drerio", data = zebrafish.GO, type = 'genesetlist', name = "ZebrafishGO", universe = unique(unlist(zebrafish.GO)))
#'
#' ##########################
#' # Step 3) The test!
#' #####
#' # Now we can compare the mouse cardiac genes to the zebrafish gene ontology.
#' mouse.cardiac.vs.zebrafish.GO.results <- run_XGSA_test(mouse.data, zebrafish.GO.data)
#'
#' ##########################
#' # Step 4) Examining the results.
#' #####
#' # We need to separate the pvalues and the overlapping gene IDs, because XGSA returns both.
#' # The p.values are stored in the first element of each result, and the overlapping genes are stored in the second element.
#' resulting.pvals <- lapply(mouse.cardiac.vs.zebrafish.GO.results, function(X){ X[[1]] })
#' resulting.overlap.genes <- lapply(mouse.cardiac.vs.zebrafish.GO.results, function(X){ X[[2]] })
#'
#' # Now we perform Benjamini Hochberg multiple hypothesis testing correction to the pvalues.
#' adjusted.pvals <- p.adjust(unlist(resulting.pvals), method = "BH")
#'
#' # We need to make the names of our results interpretable for humans, so we extract the GO Term IDs
#' names(adjusted.pvals) <- unlist(lapply(strsplit(names(adjusted.pvals) ," "), function(X){return(X[[2]])}))
#' # We can use another XGSA helper function to find out the GO term names.
#' zebrafish.GO.names <- get_GO_names('drerio')
#' # And finally we get interpretable names
#' names(adjusted.pvals) <- zebrafish.GO.names[match( names(adjusted.pvals), zebrafish.GO.names$go_id),"name_1006"]
#'
#' significant.GO.Terms <- adjusted.pvals[which(adjusted.pvals < 0.05)]
#'
#' # Now let's look at the 10 most significant GO term results
#' head(sort(significant.GO.Terms),10)
#'
#' # Zebraish cardiac development terms are significantly enriched in mouse cardiac development genes, and vice-versa.

run_XGSA_test <- function(dataset1, dataset2, test="fisher", min=5, max=500){
  # check data sets
  if(!(check_XGSA_dataset(dataset1) & check_XGSA_dataset(dataset2))){
    print("ERROR: Data is wrong")
    return(NA)
  }

  name1 <- dataset1$name
  name2 <- dataset2$name

  species1 <- dataset1$species
  species2 <- dataset2$species

  data1 <- dataset1$data
  data2 <- dataset2$data

  type1 <- dataset1$type
  type2 <- dataset2$type

  universe1 <- dataset1$universe
  universe2 <- dataset2$universe

  hm <- get_homology_matrix(species1, species2)

  if(type1 == "genesetlist" & type2 == "genesetlist"){
    if(test == "fisher"){
      results <- lapply(data1, function(X){
        results2 <- lapply(data2, function(Y){
          res <- paired_fishers_exact_tests(X, Y, hm, min=min, max=max, rowuniverse=universe1, coluniverse=universe2)
          if(!is.null(res)){
            return(res)
          }else{
            return(1)
          }})
        names(results2) <- names(data2)
        overlap2 <- lapply(data2, function(Y){
          return(get_overlap_genes(X, Y, hm, min=min, max=max))
        })
        names(overlap2) <- names(data2)
        return(list(pvals = results2, genes = overlap2))
      })
      names(results) <- names(data1)
      return(results)
    } else {
      print("No other tests implemented yet")
      return(NA)
    }
  } else {
    print("Other comparisons not implemented yet")
    return(NA)
  }
}


