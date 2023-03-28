# Author: Ben Ringham
# Licence: Unlicenced

options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings. This causes an error for me.
# suppressMessages(loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8"))

options(stringAsfactors = FALSE, useFancyQuotes = FALSE)

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("-i", "--id_file"), type="character", help="Path to file with IDs to convert"),
  make_option(c("-c","--col_num"), type="integer", help="Column identifier that contains IDs."),
  make_option(c("-o","--organism"), type="character",help="Which organism the IDs belong to"),
  make_option(c("-t","--id_type"),type="character",help="Type of the incoming IDs"),
  make_option(c("-f","--out_file"),type="character",help="Path to output file (will be created)."),
  make_option(c("-z","--header"),type="character",help="If the input file has a header"),
  make_option(c("-a","--acknowledge_limit"),type="character",help="Continue with more than 500 genes? True or false."),
  make_option(c("-m","--merge"),type="character",help="Merge output with input file? True or false.")
)

parser <- OptionParser(usage = "%prog [options]", option_list=option_list)

args = parse_args(parser)

# Vars:
id_file <- args$id_file
id_type <- args$id_type
organism <- args$organism
out_file <- args$out_file
col_num <- as.numeric(args$col_num)
file_has_header <- as.logical(args$header)
acknowledge_limit <- as.logical(args$acknowledge_limit)
LIMIT = 500
MERGE <- as.logical(args$merge)


# Adapted from annotateMyIds
if(organism == "Hs"){
    suppressPackageStartupMessages(library(org.Hs.eg.db))
    db <- org.Hs.eg.db
} else if (organism == "Mm"){
    suppressPackageStartupMessages(library(org.Mm.eg.db))
    db <- org.Mm.eg.db
} else if (organism == "Dm"){
    suppressPackageStartupMessages(library(org.Dm.eg.db))
    db <- org.Dm.eg.db
} else if (organism == "Dr"){
    suppressPackageStartupMessages(library(org.Dr.eg.db))
    db <- org.Dr.eg.db
} else if (organism == "Rn"){
    suppressPackageStartupMessages(library(org.Rn.eg.db))
    db <- org.Rn.eg.db
} else if (organism == "At"){
    suppressPackageStartupMessages(library(org.At.tair.db))
    db <- org.At.tair.db
} else if (organism == "Gg"){
    suppressPackageStartupMessages(library(org.Gg.eg.db))
    db <- org.Gg.eg.db
} else {
    cat(paste("Organism type not supported", organism))
}
  suppressPackageStartupMessages(library(jsonlite))

gene_symbols <- as.character(read.table(id_file, sep="\t", header=file_has_header)[, col_num], quote="")

suppressMessages(entrez <- select(db, keys=gene_symbols, keytype=id_type, columns=c("ENTREZID")))

gene_ids <- entrez[complete.cases(entrez),]$ENTREZID
chunk_size <- 100
cn <- ceiling(length(gene_ids) / chunk_size)

if(!acknowledge_limit & length(gene_ids) > LIMIT) {
  cat("You have more than", LIMIT, "gene ids\n")
  cat("This will result in ", cn*chunk_size, " gene descriptions being fetched.\n")
  cat("If you want to fetch data for all ", length(gene_ids), " genes, please set acknowledge_limit to FALSE.\n")
  cat("Exiting.\n")
  TOOMANYGENES()
}

result_df <- data.frame(ENTREZ = character(), SYMBOL = character(), NAME = character(), SUMMARY = character(), CHR = character(), stringsAsFactors = FALSE)

for (i in seq_len(cn)) {

  chunk_genes <- gene_ids[(chunk_size*(i-1)+1):min(chunk_size*i,length(gene_ids))]

  gids <- paste(chunk_genes, collapse = ",")

  url <- paste0('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=', gids, '&retmode=json')

  tryCatch({
    cat("Fetching data for chunk ", i, " of ", cn, "\n")
    retdata <- jsonlite::fromJSON(url)

    result <- lapply(chunk_genes, function(g) {
  if(as.character(g) %in% names(retdata$result)) {
    nomenclaturesymbol <- retdata$result[[as.character(g)]]$nomenclaturesymbol
    summary <- retdata$result[[as.character(g)]]$summary
    nomenclaturename <- retdata$result[[as.character(g)]]$nomenclaturename
    maploc <- retdata$result[[as.character(g)]]$maplocation
  } else {
    nomenclaturesymbol <- ''
    summary <- ''
    nomenclaturename <- ''
    maploc <- ''
  }
  c(g, nomenclaturesymbol, nomenclaturename, summary, maploc)
})
    
    # Create a data frame from the results
      result_df_chunk <- data.frame(do.call(rbind, result), stringsAsFactors = FALSE)

    colnames(result_df_chunk) <- c("ENTREZ", "SYMBOL","NAME", "SUMMARY", "CHR")
  }, error = function(e) {
    cat("Error fetching data for chunk ", i, ": ", e$message, "\n")
    # JSONERROR()
  })

    # Append the new results to the existing data frame
  result_df <- rbind(result_df, result_df_chunk)
}

if(MERGE) {
  original <- read.table(id_file, sep="\t", header=TRUE)
  colnames(original)[col_num] <- "SYMBOL"
  merged <- merge(original, result_df, by = "SYMBOL", all.x = TRUE) # merge the two data frames by gene_id
  write.table(merged, out_file, row.names = FALSE, col.names = TRUE, sep="\t")
} else {
    write.table(result_df, out_file, row.names = FALSE, col.names = TRUE, sep="\t")
}
