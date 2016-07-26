# XGSA.R
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
# This file contains the neccessary R code to run XGSA as described in "XGSA: an unbiased statistical method for cross-species gene set analysis"
# XGSA is implemented in a framework that will create several global variables and objects to improve performance.
# It is important that users do not have their own objects named the following:
#
# host.ID
# biomart.ID
# homology.matrix.list
# supported.species
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
# XGSA depends on two libraries, Matrix and biomaRt
# To install these execute the following commands:
#
# source("http://bioconductor.org/biocLite.R")
# biocLite("biomaRt")
# install.packages("Matrix")
# 
##########################

require(Matrix)
require(biomaRt)

##########################
# When Ensembl move servers (rarely) or it is down (less rare) these variable may need to be changed, for example to an archived version
biomart.ID <- "ENSEMBL_MART_ENSEMBL"
host.ID <-  "www.ensembl.org"
#host.ID <-  "dec2015.archive.ensembl.org"


find_supported_datasets <- function(default=TRUE){
  ensembl <- useMart(biomart.ID, dataset='hsapiens_gene_ensembl', host=host.ID)
  sets <- listDatasets(ensembl)
  organisms <- gsub("([a-z]*)_.*","\\1",sets[,1],"_")
  return(organisms)
}

##########################
# This is stored globally at the start of the session to avoid having to be recalled
supported.species <- find_supported_datasets()

##########################
# This function executes the paired Fisher's exact tests which form the basis of XGSA
# It returns a single numeric p-value
# Input parameters are the Ensembl gene IDs in the two species, corresponding to the rows and columns of the homology matrix, the homology matrix itself, maximum and minimum gene set sizes to be considered, and the two gene universes to consider.
paired_fishers_exact_tests<-function(row_genes, col_genes, homology_matrix, min=5, max=500, rowuniverse, coluniverse){
    
    #only test those genes that actually exist in the matrix
    row_genes <- unique(as.character(row_genes[row_genes%in%rownames(homology_matrix)]))
    col_genes <- unique(as.character(col_genes[col_genes%in%colnames(homology_matrix)]))

    # figure out the correct universe size
    rowuniverse <- as.character(rowuniverse[rowuniverse%in%rownames(homology_matrix)])
    coluniverse <- as.character(coluniverse[coluniverse%in%colnames(homology_matrix)])
    
    rowuniverse <- intersect(rowuniverse, rownames(homology_matrix)[which(rowSums(homology_matrix[,coluniverse, drop=FALSE])>0)])
    coluniverse <- intersect(coluniverse, colnames(homology_matrix)[which(colSums(homology_matrix[rowuniverse,, drop=FALSE])>0)])
    
    
    #check size of gene sets
    if((length(row_genes) < min) | (length(col_genes) < min) | (length(row_genes) > max) | (length(col_genes) > max)) { return(NULL) }
    
    #share orthology information by predicting gene set genes
    predicted_col_genes <- colnames(homology_matrix)[which(colSums(homology_matrix[row_genes,, drop=FALSE])>0)]
    predicted_row_genes <- rownames(homology_matrix)[which(rowSums(homology_matrix[,col_genes, drop=FALSE])>0)]
    
    rowhits <- intersect(row_genes, predicted_row_genes)
    colhits <- intersect(col_genes, predicted_col_genes)
    
    row.p <- fisher.test(matrix(c(length(rowhits),length(predicted_row_genes)-length(rowhits),length(row_genes)-length(rowhits),length(rowuniverse)-length(unique(union(row_genes, predicted_row_genes)))),2,2), alternative='greater')$p.value 
    col.p <- fisher.test(matrix(c(length(colhits),length(predicted_col_genes)-length(colhits),length(col_genes)-length(colhits),length(coluniverse)-length(unique(union(col_genes, predicted_col_genes)))),2,2), alternative='greater')$p.value 
    
    
    return(max(row.p, col.p))
}


##########################
# This function returns the genes that overlap between the two sets
# It returns a list with two vectors of gene IDs, matching thw row and column species of the homology table
# Input parameters are the Ensembl gene IDs in the two species corresponding to the rows and columns of the homology matrix, the homology matrix itself, maximum and minimum gene set sizes to be considered.
get_overlap_genes<-function(row_genes, col_genes, homology_matrix, min=5, max=500){
  
  #only test those genes that actually exist in the matrix
  row_genes <- as.character(row_genes[row_genes%in%rownames(homology_matrix)])
  col_genes <- as.character(col_genes[col_genes%in%colnames(homology_matrix)])
  
  #check size of gene sets
  #if((length(row_genes) < min) | (length(col_genes) < min) | (length(row_genes) > max) | (length(col_genes) > max)) { return(NULL) }
  
  #share orthology information by predicting gene set genes
  predicted_col_genes <- colnames(homology_matrix)[which(colSums(homology_matrix[row_genes,, drop=FALSE])>0)]
  predicted_row_genes <- rownames(homology_matrix)[which(rowSums(homology_matrix[,col_genes, drop=FALSE])>0)]
  
  rowhits <- intersect(row_genes, predicted_row_genes)
  colhits <- intersect(col_genes, predicted_col_genes)
  
  overlap <- list(rowhits, colhits)
  names(overlap) <- c("rowhits","colhits")
  return(overlap)
}

#loading functions

##########################
#this function generates a sparse matrix from a table with two or three columns: ID ID Value(optional)
generate_sparse_matrix <- function(homology_table){
  # Compare only genes with homology mapping
  homology_table <- unique(homology_table)[!unique(homology_table)[,2]=="",]
  
  if(ncol(homology_table)==3){
    hms <- sparseMatrix(i = match(homology_table[,1], unique(homology_table[,1]) ), j = match(homology_table[,2], unique(homology_table[,2])), x = homology_table[,3])
  } else if(ncol(homology_table)==2){
    hms <- sparseMatrix(i = match(homology_table[,1], unique(homology_table[,1]) ), j = match(homology_table[,2], unique(homology_table[,2])))
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
get_homology_table <- function(species1, species2, sequence_identity_reference = 0){
  require(biomaRt)
  dataset_name <- paste(species1, "_gene_ensembl", sep="")
  homolog_attribute <- paste(species2, "_homolog_ensembl_gene", sep="") 
  ensembl <- useMart(biomart.ID, dataset=dataset_name, host=host.ID)
  
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
get_ENSEMBL_gene_list <- function(species1){
  require(biomaRt)
  dataset_name <- paste(species1, "_gene_ensembl", sep="")
  ensembl <- useMart(biomart.ID, dataset=dataset_name, host=host.ID)
  gene.list <- getBM(attributes=c("ensembl_gene_id"), mart = ensembl)
  return(gene.list)  
}

##########################
# Useful to convert from transcript IDs to gene IDs
# Returns a data frame with 2 columns
# Input parameter is the species string (eg. 'hsapiens')
get_ENSEMBL_trasncript_table <- function(species){
  require(biomaRt)
  dataset_name <- paste(species, "_gene_ensembl", sep="")
  ensembl <- useMart(biomart.ID, dataset=dataset_name, host=host.ID)
  transcript.table <- getBM(attributes=c("ensembl_gene_id","ensembl_transcript_id"), mart = ensembl)
  return(transcript.table)  
}

##########################
# Retrieves the mapping  between ENSEMBL ids and gene symbols
# Returns a data frame with 2 columns
# Input parameter is the species string (eg. 'hsapiens')
get_ENSEMBL_symbol_map <- function(species){
  require(biomaRt)
  dataset_name <- paste(species, "_gene_ensembl", sep="")
  ensembl <- useMart(biomart.ID, dataset=dataset_name, host=host.ID)
  ENSEMBL_symbol_map <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"), mart = ensembl)
  return(ENSEMBL_symbol_map)
}

##########################
# Retrieves the GO term names and IDs of the latest Gene Ontology annotations from ENSEMBL
# Returns a data frame containing GO terms and term names
# Input parameter is the species string (eg. 'hsapiens')
get_GO_names <- function(species){
  require(biomaRt)
  dataset_name <- paste(species, "_gene_ensembl", sep="")
  ensembl <- useMart(biomart.ID, dataset=dataset_name, host=host.ID)
  GOTerms <- getBM(mart=ensembl, attributes=c('go_id','name_1006'))
  return(GOTerms)
}

##########################
# Retrieves the latest Gene Ontology annotations from ENSEMBL
# Returns a data frame containing gene ids, associated go terms, GO evidence codes, and GO ontology type (BP, MF, CC)
# Input parameter is the species string (eg. 'hsapiens')
get_GO_mappings <- function(species){
  require(biomaRt)
  dataset_name <- paste(species, "_gene_ensembl", sep="")
  ensembl <- useMart(biomart.ID, dataset=dataset_name, host=host.ID)
  EG2GO <- getBM(mart=ensembl, attributes=c('ensembl_gene_id','go_id','go_linkage_type','namespace_1003'))
  return(EG2GO)
} 

##########################
# Retrieves GO mappings and converts the result into a list
# Returns a list of GO terms with constituent ENESMBL gene IDS
# Input parameter is the species string (eg. 'hsapiens')
get_GO_list <- function(species){
  GOMap <- get_GO_mappings(species)
  GOTerms <- split(GOMap$ensembl_gene_id, GOMap$go_id)
  return(GOTerms)
}

##########################
# Retrieves GO mappings with only the supplied evidence codes and converts the result into a list
# Returns a list of GO terms with constituent ENESMBL gene IDS
# Input parameters are the species string (eg. 'hsapiens') and a character vector of evidence codes
# Defaults to non-computational evidence 
get_GO_list_with_evidence_codes <- function(species, evidence.codes = c("EXP","IDA","IPI","IMP","IGI","IEP", "TAS", "IC")){
  GOMap <- get_GO_mappings(species)
  subGOMap <- GOMap[GOMap$go_linkage_type %in% evidence.codes,]
  GOTerms <- split(subGOMap$ensembl_gene_id, subGOMap$go_id)
  return(GOTerms)
}

##########################
# Retrieves GO mappings with only the supplied evidence codes and from only the supplied ontologies and converts the result into a list
# Returns a list of GO terms with constituent ENESMBL gene IDS
# Input parameters are the species string (eg. 'hsapiens'), a character vector of evidence codes and a character vector of ontology codes (eg. 'BP')
# Defaults to non-computational evidence 
get_GO_list_from_ontologies_with_evidence_codes <- function(species, evidence.codes=  c("EXP","IDA","IPI","IMP","IGI","IEP", "TAS", "IC"), ontologies){
  GOMap <- get_GO_mappings(species)
  subGOMap <- GOMap[(GOMap$go_linkage_type %in% evidence.codes)&(GOMap$namespace_1003 %in% ontologies),]
  GOTerms <- split(subGOMap$ensembl_gene_id, subGOMap$go_id)
  return(GOTerms)
}

##########################
# Returns a list with two elements (one for each species). 
# Each element contains a vector of gene IDs in the species that have more than 1 homolog in the other species
# Input parameters are the two species strings (eg. 'hsapiens', 'drerio') 
get_complex_genes <- function(species1, species2){
  hm <- get_homology_matrix(species1, species2)
  hm.boolean <- hm > 0
  complex.cols <- which(colSums(hm.boolean) > 1)
  complex.rows <- which(rowSums(hm.boolean) > 1)
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
calculate_genesetlist_complexity <- function(dataset1, species2){
  #check species aren't the same
  species1 <- dataset1$species
  if( species1 == species2 ){
    print("ERROR: no complex homology between the same species")
    return(NULL)
  }
  complex.genes <- get_complex_genes(species1, species2)
  data1 <- dataset1$data
  data1.complexity <- lapply(data1, function(X){
    return(length(intersect(X, complex.genes[[species1]])) / length(X))
  })
  return(data1.complexity)
}

##########################
# Calculates an alternative complexity score (not described in the paper) of the gene sets in an XGSA data set with respect to a second species
# Returns a list of complexity scores, one score for each gene set
# Input parameters are the XGSA data set object and the other species string (eg. 'drerio') 
calculate_alternate_genesetlist_complexity <- function(dataset1, species2){
  #check species aren't the same
  species1 <- dataset1$species
  if( species1 == species2 ){
    print("ERROR: no complex homology between the same species")
    return(NULL)
  }
  hm <- get_homology_matrix(species1, species2)
  data1 <- dataset1$data
  data1.complexity <- lapply(data1, function(X){
    return(sum(rowSums(hm[X[X%in%rownames(hm)],]))/sum(X%in%rownames(hm)))
  })
  return(data1.complexity)
}

##########################
# Creates a homology matrix between two species by retrieving the homology table and then converting it into a sparse matrix
# Returns a sparse matrix
# Input parameters are the two species strings (eg. 'hsapiens', 'drerio') and a boolean of whether to return the sequence identity values or not 
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
      # create a new homology matrix and return it
      ht <- get_homology_table(species1, species2)
      hm <- generate_sparse_matrix(ht)
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
      # create a new homology matrix and return it
      ht <- get_homology_table(species1, species2, 1)
      hm <- generate_sparse_matrix(ht)
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

new_XGSA_dataset <- function(species, data, type, name, universe="empty"){
  #check species
  if(!species %in% supported.species){
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
  
  
  # check if matrix exists
  
  if(species1 == species2){
    
    gene.list <- get_ENSEMBL_gene_list(species1)
    hm <- generate_sparse_matrix(cbind(gene.list, gene.list))
    
    # DO a same species analysis
    
  } else {
    
    hm <- get_homology_matrix(species1, species2)
    
  }
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
        return(list(results2,overlap2))
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


