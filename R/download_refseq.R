{
  taxonomy_link <- "https://raw.githubusercontent.com/Xudong-Han/MetaScope/master/data/taxonomy.txt"
  taxonomy.table <- read.table(taxonomy_link, header = T, sep = "\t")
  taxon_list <- c(taxonomy.table$superkingdom, taxonomy.table$kingdom,
                  taxonomy.table$phylum, taxonomy.table$class, taxonomy.table$order,
                  taxonomy.table$family, taxonomy.table$genus, taxonomy.table$species,
                  taxonomy.table$strain)
  taxon_list <- taxon_list[!is.na(taxon_list)]
  ## check if user provided a vaild taxon
  if (!taxon %in% taxon_list){
    stop("Your input is not a valid taxon")
  }
  ## get the rank of the input taxon
  classification.table <- classification(get_uid(taxon,messages=F)[[1]],
                                         db = 'ncbi')[[1]]
  rank_input <- classification.table$rank[nrow(classification.table)]
  message(taxon," is a ",rank_input)
  
  ## get the ncbi scinames of children species or strains
  children_list <- get_children(taxon,rank_input,data = taxonomy.table)
  
  ## get the parent taxon in superkingdom rank

  ## download the assembly summary refseq table from ncbi which includes genome download link
  refseq_link <- "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary.txt"
  refseq_table <- read.table(refseq_link, header = T, sep = "\t",
                             comment.char = "",
                             quote = "", skip = 1)
  ## filter the table, keep the lines with species or strains of input
  species_table <- refseq_table[which(refseq_table$organism_name %in% children_list),]
  
  ## Reduce the table size based on reference or represenative
  if (representative) {
    reference <- TRUE
  }
  if (representative & reference) {
    species_table <- species_table[species_table$refseq_category %in%
                                     c("reference genome",
                                       "representative genome"), ]
  } else {
    if (!representative & reference) {
      species_table <- species_table[species_table$refseq_category ==
                                       "reference genome", ]
    } 
  }
  total_genomes <- nrow(species_table)
  if (total_genomes == 0){
    message("No availabe genome for ",taxon)
  } else{
    message(paste("Downloading", total_genomes, taxon, "genomes from RefSeq"))
    
    ## delete existing genome files and combined fasta--make these user-defined
    
    download_dir <- paste(taxon, "refseq_download", sep = "_")
    # remove any existing files/directories
    unlink(download_dir, recursive = TRUE, force = TRUE)
    if (compress) {
      combined_fasta <- paste(taxon, "fasta.gz", sep = ".")
      combined_fasta_patho <- paste(taxon, "pathoscope.fasta.gz", sep = ".")
    } else {
      combined_fasta <- paste(taxon, "fasta", sep = ".")
      combined_fasta_patho <- paste(taxon, "pathoscope.fasta", sep = ".")
    }
    
    tryCatch({
      suppressWarnings(file.remove(combined_fasta))
      suppressWarnings(file.remove(combined_fasta_patho))
    }) # start with a new combined file
    
    ## Download the genome
    for (i in 1:nrow(species_table)) {
      tryCatch({
        if (i%%10 == 0) {
          message(paste("Number of Genomes Downloaded: ", i, "/",
                        total_genomes, " (", round(100 * i/total_genomes, 2),
                        "%)", sep = ""))
        }
        
        ## Download the genome
        genome_file <- paste(basename(as.character(species_table[i, ]$ftp_path)),
                             "genomic.fna.gz", sep = "_")
        location <- paste(species_table[i, ]$ftp_path, genome_file, sep = "/")
        destination <- paste(download_dir, genome_file, sep = "/")
        if (!dir.exists(download_dir)) {
          dir.create(download_dir)
        }
        download.file(location, destination)
        
        ## read in the genome
        ref <- Biostrings::readDNAStringSet(destination)
        
        ## write to file
        Biostrings::writeXStringSet(ref, combined_fasta, append = T,
                                    compress = compress)
        
        ## format for pathoscope and write to file
        accession <- NULL
        for (j in strsplit(names(ref), " ")) {
          accession <- c(accession, j[1])
        }
        names(ref) <- paste("ti|", species_table[i, ]$taxid, "|org|",
                            gsub(" ", "_", species_table[i, ]$organism_name),
                            "|accession|", accession, sep = "")
        Biostrings::writeXStringSet(ref, combined_fasta_patho, append = T,
                                    compress = compress)
        
        ## delete intermediate download files
        unlink(download_dir, recursive = TRUE, force = TRUE)
      }, error = function(e) {
        cat("ERROR :", conditionMessage(e), "\n")
      })
    }
    # Remove pathoscope file if unwanted
    if (!patho_out) file.remove(combined_fasta_patho)
    
    # Ensure removal of intermediate folder of files
    unlink("viral_refseq_download", recursive = TRUE, force = TRUE)
    
    message("DONE! Downloaded ", i, " genomes to ", combined_fasta)
  }
}
