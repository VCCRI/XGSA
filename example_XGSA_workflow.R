# example_XGSA_workflow.R
##########################
#
# This file contains R code to run an example XGSA analysis. We will compare mouse cardiac genes with the Zebrafish Gene Ontology. 
# Make sure you have installed the packages "biomaRt" and "Matrix", the commands to do so are below.
#
# solving issues with libcurl and XML dependencies:
# sudo yum install curl-develsudo yum install curl-devel
# sudo apt-get install libcurl4-openssl-dev libxml2-dev
#
# solving issues with https from unix "In file(filename, "r", encoding = encoding) : unsupported URL scheme"
# library(RCurl)
#    eval( expr = 
#                 parse( text = getURL("https://raw.githubusercontent.com/VCCRI/XGSA/master/XGSA.R",
#                                                             ssl.verifypeer=FALSE) ))
#
#
#
# source("http://bioconductor.org/biocLite.R")
# biocLite("biomaRt")
# install.packages("Matrix")
#
##########################
# Step 1) load XGSA functions
#####
source("https://raw.githubusercontent.com/VCCRI/XGSA/master/XGSA.R")

##########################
# Step 2) Create sample mouse gene set.
#####
# In this example we will extract a list of mouse cardiac development genes from a previous published study.
cardiac.perturbation.data <- read.table("http://cardiaccode.victorchang.edu.au/data/Cardiaccode_in_vivo_evidence_2014_06_10.txt", sep="\t", header=TRUE, quote="\"")
mouse.cardiac.genes <- unique(unlist(cardiac.perturbation.data[cardiac.perturbation.data$Species == "Mus musculus", c("Regulator","Target")]))

# XGSA harnesses the Ensembl homology pipeline, and so we need to convert all of our data sets to Ensembl gene IDs
# We will convert the mouse gene symbols to Ensembl IDs using the XGSA helper function "get_ENSEMBL_symbol_map".
mouse.ensembl.symbol.map <- get_ENSEMBL_symbol_map(species = 'mmusculus')
mouse.cardiac.ensembl.symbols <- mouse.ensembl.symbol.map$ensembl_gene_id[mouse.ensembl.symbol.map$external_gene_name %in% mouse.cardiac.genes]

# Now we have a list of Ensembl IDs for mouse cardiac genes, we will turn it into an XGSA data set.
# Note that the input data MUST be a named list, which allows for multiple gene sets in the same data set.
# As we don't have a defined gene universe for this data set we will use all mouse Ensembl IDs that have an external gene symbol.
mouse.data <- new_XGSA_dataset(species = 'mmusculus', data = list(mouseCardiacGenes = mouse.cardiac.ensembl.symbols), type = 'genesetlist', name = 'MouseCardiacGenes', universe = unique(mouse.ensembl.symbol.map$ensembl_gene_id))

##########################
# Step 3) Create reference zebrafish GO set.
#####
# In this example we will compare to the zebrafish Gene Ontology using "direct" evidence only - this means the annotations are NOT transferred between species.
# We will use another XGSA helper function to retrieve the latest Gene Ontology information from Ensembl "get_GO_list_from_ontologies_with_evidence_codes".
# The gene universe we will use is all ofthe zebrafish biological process genes that we are testing. 
zebrafish.GO <- get_GO('drerio', ontologies = "biological_process")
zebrafish.GO <- zebrafish.GO[lapply(zebrafish.GO, length) > 10 & lapply(zebrafish.GO, length) < 500]
zebrafish.GO.data <- new_XGSA_dataset(species = "drerio", data = zebrafish.GO, type = 'genesetlist', name = "ZebrafishGO", universe = unique(unlist(zebrafish.GO)))

##########################
# Step 4) The test!
#####
# Now we can compare the mouse cardiac genes to the zebrafish gene ontology.
mouse.cardiac.vs.zebrafish.GO.results <- run_XGSA_test(mouse.data, zebrafish.GO.data)

##########################
# Step 5) Examining the results.
#####
# We need to separate the pvalues and the overlapping gene IDs, because XGSA returns both.
# The p.values are stored in the first element of each result, and the overlapping genes are stored in the second element.
resulting.pvals <- lapply(mouse.cardiac.vs.zebrafish.GO.results, function(X){ X[[1]] })
resulting.overlap.genes <- lapply(mouse.cardiac.vs.zebrafish.GO.results, function(X){ X[[2]] })

# Now we perform Benjamini Hochberg multiple hypothesis testing correction to the pvalues.
adjusted.pvals <- p.adjust(unlist(resulting.pvals), method = "BH")

# We need to make the names of our results interpretable for humans, so we extract the GO Term IDs
names(adjusted.pvals) <- unlist(lapply(strsplit(names(adjusted.pvals) ,"\\."), function(X){return(X[[2]])}))
# We can use another XGSA helper function to find out the GO term names.
zebrafish.GO.names <- get_GO_names('drerio')
# And finally we get interpretable names
names(adjusted.pvals) <- zebrafish.GO.names[match( names(adjusted.pvals), zebrafish.GO.names$go_id),"name_1006"]

significant.GO.Terms <- adjusted.pvals[which(adjusted.pvals < 0.05)]

# Now let's look at the 10 most significant GO term results
head(sort(significant.GO.Terms),10)
#cell fate commitment                             mesoderm development 
#8.363504e-09                                     2.305281e-07 
#cell fate specification                 embryonic heart tube development 
#2.734734e-07                                     9.498289e-07 
#cardiocyte differentiation                         diencephalon development 
#4.236265e-06                                     4.506317e-06 
#heart morphogenesis                             endoderm development 
#9.151835e-06                                     3.238483e-05 
#heart development regulation of endodermal cell fate specification 
#3.238483e-05                                     3.238483e-05 
#
# Zebraish cardiac development terms are significantly enriched in mouse cardiac development genes, and vice-versa.
#
# Please now try your own analysis! Find some gene sets (make them a decent size, 10 - 500 is nice) from any of the supported species and compare them to each other or a reference set like the gene ontology.
# To find supported species check:
supported.species
#[1] "oanatinus"         "cporcellus"        "gaculeatus"        "itridecemlineatus" "lafricana"         "choffmanni"       
#[7] "csavignyi"         "fcatus"            "rnorvegicus"       "psinensis"         "cjacchus"          "ttruncatus"       
#[13] "scerevisiae"       "celegans"          "csabaeus"          "oniloticus"        "amexicanus"        "trubripes"        
#[19] "pmarinus"          "eeuropaeus"        "falbicollis"       "etelfairi"         "cintestinalis"     "ptroglodytes"     
#[25] "nleucogenys"       "sscrofa"           "ocuniculus"        "dnovemcinctus"     "pcapensis"         "tguttata"         
#[31] "mlucifugus"        "hsapiens"          "pformosa"          "tbelangeri"        "mfuro"             "ggallus"          
#[37] "xtropicalis"       "ecaballus"         "pabelii"           "drerio"            "xmaculatus"        "tnigroviridis"    
#[43] "lchalumnae"        "amelanoleuca"      "mmulatta"          "pvampyrus"         "panubis"           "mdomestica"       
#[49] "acarolinensis"     "vpacos"            "tsyrichta"         "ogarnettii"        "dmelanogaster"     "mmurinus"         
#[55] "loculatus"         "olatipes"          "oprinceps"         "ggorilla"          "dordii"            "oaries"           
#[61] "mmusculus"         "mgallopavo"        "gmorhua"           "saraneus"          "aplatyrhynchos"    "sharrisii"        
#[67] "meugenii"          "btaurus"           "cfamiliaris" 
