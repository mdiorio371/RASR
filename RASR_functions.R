
list.of.packages <-
  c(
    "tidyverse", "data.table", "tidyr", 
    "lubridate", "ggplot2",
    "curl", "RColorBrewer", "ape",
    "pbmcapply", "igraph", "viridis",
    "gridExtra", "gridExtra", "scales", 
    "knitr", "RCurl", "kableExtra", "boot", 
    "car", "modelr", "broom", "ggpubr", "cowplot",
    "ggridges", "R.utils", "MASS", ""
  )
new.packages <- 
  list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0) install.packages(new.packages)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}


library(tidyverse)
library(data.table)
library(pbmcapply)
library(tidyr)
library(lubridate) 
library(Biostrings)
library(ggplot2)
library(curl)
library(RColorBrewer)
library(ape)
library(igraph)
library(viridis)  
library(gridExtra)
library(scales)
library(knitr)
library(RCurl)
library(kableExtra)
library(boot)
library(car)
library(modelr)
library(broom)
#install.packages("gggenomes")
library(gggenomes)
library(ggpubr)
library(cowplot)
#install.packages("ggplotify")
library(ggplotify)
library(ggtree)
library(ggridges)
library(R.utils)
library(MASS)
library(TOSTER)


get_ncbi_table <- 
  function(db_dir = "db"){
    if (!dir.exists(db_dir)){dir.create(db_dir)}
    ncbi_dir <- 
      "ftp://ftp.ncbi.nlm.nih.gov/genomes"
    # refseq
    refseq_file <- 
      sprintf(
        "%s/refseq/bacteria/assembly_summary.txt",
        ncbi_dir
      )
    
    # genbank genomes 
    genbank_file <- 
      sprintf(
        "%s/genbank/bacteria/assembly_summary.txt",
        ncbi_dir
      )
    refseq <- 
      fread(
        refseq_file, 
        stringsAsFactors = F, 
        quote = "", 
        select = c(1,5,6,7,8,11,12,15,16,17,18,19,20)
      ) %>% 
      filter(
        assembly_level =="Complete Genome" &
          version_status =="latest"
      ) %>%
      `colnames<-`(c("assembly", colnames(.)[-1])) %>%
      as_tibble %>%
      mutate(db = "refseq")
    ## Access the Genbank genomes that are not in refseq
    genbank <- 
      fread(
        genbank_file, 
        stringsAsFactors = F, 
        quote = "", 
        select = c(1,5,6,7,8,11,12,15,16,17,18,19,20)
      ) %>% 
      filter(
        assembly_level =="Complete Genome" &
          version_status =="latest"
      ) %>%
      `colnames<-`(c("assembly", colnames(.)[-1])) %>%
      as_tibble %>%
      filter(paired_asm_comp != "identical") %>%
      mutate(db = "genbank")
    
    # refseq & genbank combined  
    bac_summary <- 
      bind_rows(refseq, genbank) %>%
      mutate(
        fna_path = 
          sprintf(
            "%s/%s_genomic.fna.gz",
            ftp_path, gsub(".*/", "", ftp_path)
          )
      )
    
    # bac_summary %>%
    #   mutate(
    #     seq_rel_date = as.Date(seq_rel_date)      
    #   ) %>%
    #   filter(seq_rel_date <= "2021-01-26")
    
    # get taxonomy information
    taxa_file <- 
      "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz"
    
    taxa_con <- 
      curl_download(
        taxa_file,
        destfile = "db/new_taxdump.tar.gz"
      )
    ### extract the ranked lineage table
    untar(taxa_con, exdir = db_dir, files = "rankedlineage.dmp")
    ## read it in
    taxa_table <- 
      read_tsv(
        "db/rankedlineage.dmp", 
        col_names = 
          c(
            "taxid", "name", "s", "g", "f", "o", "c","p", "k", "d"
          ),
        col_types=("i-c-c-c-c-c-c-c-c-c-")
      ) %>% 
      filter(d=="Bacteria")
    ## Combinte the tables
    bac_taxa <- 
      left_join(bac_summary, taxa_table, by = "taxid") %>% 
      mutate(
        p2 = 
          ifelse(
            p=="Proteobacteria", 
            sub("proteobacteria", "",c), 
            p
          )
      ) %>%
      dplyr::select(-c(version_status, assembly_level, k, d)) %>%
      mutate(
        species = 
          unlist(
            pbmclapply(
              organism_name,
              function(x){
                spec_out <-
                  strsplit(x, " ")[[1]][1:2] %>% paste(collapse = " ")
                spec_out_2 <-
                  strsplit(spec_out, " ")[[1]][2]
                if (spec_out_2=="sp."){
                  spec_out <-
                    x
                  #strsplit(x, " ")[[1]][1:3] %>% paste(collapse = " ")
                }
                return(spec_out)
              }
            )
          )
      )
    ## note sys time 
    return(bac_taxa)
  }

get_report_table <-
  function(org_table, core_number=7){
    
    report_details <-
      pbmclapply(
        (org_table %>% pull(ftp_path)),
        function(ftp_path){
          #tester <- list()
          # for (i in (1:length((org_table %>% pull(ftp_path))))){
          # ftp_path <-
          #    (org_table %>% pull(ftp_path))[1]
          if (ftp_path=="na" ){
            report_details <- 
              tibble(
                asm_name = NA,
                org_name = NA,
                submitter = NA,
                accession = NA, 
                len = NA, 
                strain = NA,
                date_sub = NA,
                technology = NA
              )
          } else {
            rep_file <- 
              sprintf(
                "%s/%s_assembly_report.txt",
                ftp_path, gsub(".*/", "", ftp_path)
              ) 
            rep_read <-
              tryCatch(readLines(rep_file))
            #tester[[i]] <- rep_read
            
            #fread(., fill=TRUE, stringsAsFactors = F)
            strain <- 
              rep_read[grepl("# Infraspecific name:", rep_read)] %>% 
              sub(".*strain=", "", .)
            if (length(strain)==0){
              strain <- NA
            }
            technology <- 
              rep_read[grepl("# Sequencing", rep_read)] %>% 
              sub(".*technology:", "", .) %>% 
              trimws
            if (length(technology)==0){
              technology <- NA
            }
            out_tib <- 
              tibble(
                asm_name = 
                  sub(".*: ", "", rep_read[1]) %>%
                  trimws,
                org_name = 
                  sub(".*: ", "", rep_read[2]) %>%
                  trimws,
                submitter = 
                  rep_read[grepl("Submitter", rep_read)] %>% 
                  sub(".*Submitter:", "", .) %>% 
                  trimws,
                accession = 
                  rep_read[grepl("Chromosome", rep_read)] %>% 
                  sub(".*=\\t", "", .) %>% 
                  sub("\t.*", "", .),
                len = 
                  rep_read[grepl("Chromosome", rep_read)] %>% 
                  sub(".*Assembly\\t", "", .) %>%
                  sub("\t.*", "", .),
                strain = 
                  strain,
                date_sub = 
                  rep_read[grepl("Date", rep_read)] %>% 
                  sub(".*Date:", "", .) %>% 
                  trimws,
                technology = 
                  technology
                #ftp_path =
                # ftp_path
              )
          }
          
          return(out_tib)
        }, mc.cores = core_number
      ) %>%
      bind_rows()
    # report_details[
    #     which((lapply(report_details, nrow) %>% 
    #                unlist)>1)
    # ] 
    return(report_details)
  }

generate_dirs <- 
  function(phylum_i, genus_i, org_i){
    sync_dir <- 
      sprintf(
        "db/synced/%s/%s_%s",
        phylum_i, genus_i, org_i
      ) 
    gff_dir <- 
      sprintf(
        "db/gff/%s/%s_%s",
        phylum_i, genus_i, org_i
      ) 
    
    identity_dir <- 
      sprintf(
        "identity/%s/%s_%s",
        phylum_i, genus_i, org_i
      ) 
    
    nucmer_dir <- 
      sprintf(
        "nucmer/%s/%s_%s",
        phylum_i, genus_i, org_i
      ) 
    r_dir <- 
      sprintf(
        "r_vars/%s/%s_%s",
        phylum_i, genus_i, org_i
      ) 
    
    results_dir <- 
      sprintf("results/%s/%s_%s",
              phylum_i, genus_i, org_i
      ) 
    r_dir <- 
      sprintf("r_vars/%s/%s_%s",
              phylum_i, genus_i, org_i
      ) 
    
    if (!dir.exists(results_dir)){dir.create(results_dir, recursive = T)}
    if (!dir.exists(r_dir)){dir.create(r_dir, recursive = T)}
    if (!dir.exists(identity_dir)){dir.create(identity_dir, recursive = T)}
    if (!dir.exists(sync_dir)){dir.create(sync_dir, recursive = T)}
    if (!dir.exists(gff_dir)){dir.create(gff_dir, recursive = T)}
    if (!dir.exists(nucmer_dir)){dir.create(nucmer_dir, recursive = T)}
  }

download_gff <- 
  function(report_details, gff_dir, core_number=7){
    gff_table <- 
      pbmclapply(
        1:nrow(report_details),
        function(i){
          ## get the gff path
          asm_name <- 
            pull(report_details, asm_name)[i]
          gff_path <- 
            (report_details %>% 
               pull(ftp_path))[i] %>%
            sprintf(
              "%s/%s_genomic.gff.gz",
              ., gsub(".*/", "", .)
            )
          ## Set the download path
          gff_outfile <- 
            sprintf(
              "%s/%s",
              gff_dir,
              basename(gff_path)
            )
          # if it's there, don't downlaod it
          if (file.exists(gff_outfile)){
            gff_tib <- 
              tibble(
                gff = basename(gff_outfile),
                gff_exists = T,
                asm_name = asm_name
              )
          } else{ # if not, try downloading it
            test_downlaod <- 
              try(
                curl_download(
                  url = 
                    gff_path,
                  destfile = 
                    gff_outfile
                ), 
                silent = T
              )
            gff_tib <- 
              tibble(
                gff = basename(gff_outfile),
                gff_exists = !grepl("Error", test_downlaod),
                asm_name = asm_name
              )
          }
          return(gff_tib)
        }, mc.cores = core_number
      ) %>% 
      bind_rows()
    return(gff_table)
  }


sync_with_gff <- 
  function(
    full_species_table, sync_dir, 
    gff_dir, synced_gff_dir, core_number=7
    ){
    ## check for dnaAs
    #####
    sync_report <- 
      pbmclapply(
        1:nrow(full_species_table),
        function(i){
          #for (i in 1:nrow(report_details)){
          ## get the gff path
          gff_path <- 
            sprintf(
              "%s/%s_genomic.gff.gz",
              gff_dir,
              basename((full_species_table %>% pull(ftp_path))[i])
            )
          
          ## Set the download path
          gff_outfile <- 
            sprintf(
              "%s/%s",
              gff_dir,
              basename(gff_path)
            )
          genome_acc <- 
            (full_species_table %>% pull(accession))[i]
          out_file <- 
            sprintf(
              "%s/%s.txt",
              sync_dir, genome_acc
            )
          gff <- 
            read.gff(gff_outfile)
          
          ## Check if it's been synced
          if (file.exists(out_file)) {
            out_string <- 
              sprintf(
                "%s already synced",
                genome_acc
              )
          } else if (nrow(gff)==0){
            out_string <- 
              sprintf(
                "%s has no GFF file",
                genome_acc
              )
          } else{
            #for (i in 1:length(fna_paths)){
            asm_path <- 
              (full_species_table %>% pull(ftp_path))[i]
            genome_file <- 
              sprintf(
                "%s/%s_genomic.fna.gz",
                asm_path, 
                gsub(".*/", "", asm_path)
              )
            genome_seq <- 
              readDNAStringSet(genome_file)
            genome_seq <- 
              genome_seq[
                !(grepl("plasmid", names(genome_seq),  ignore.case = T))
                ]
            genome_seq <- 
              genome_seq[
                !(grepl("phage", names(genome_seq),  ignore.case = T))
                ]
            # take the bigger sequence ie the main one
            genome_seq <- 
              genome_seq[
                width(genome_seq)==max(width(genome_seq))
                ]
            genome_length <- 
              width(genome_seq)

            dnaA_vect <- 
              gff %>% 
              as_tibble() %>% 
              filter(
                (grepl(
                  "gene=dnaA", attributes, ignore.case = T
                ) |
                  grepl(
                    "product=chromosome replication initiator DnaA", 
                    attributes, ignore.case = T
                  ) |
                  grepl(
                    "product=chromosomal replication initiator protein DnaA", 
                    attributes, ignore.case = T
                  )) &
                  type == "CDS"
              ) %>% as.data.frame()
            
            if (nrow(dnaA_vect)>1){
              if (all(pull(dnaA_vect, strand)=="+")){
                dnaA_vect <- 
                  dnaA_vect %>%
                  filter(start ==min(start))
              } else if (all(pull(dnaA_vect, strand)=="-")) {
                dnaA_vect <- 
                  dnaA_vect %>%
                  filter(start ==max(start))
              } else {
                dnaA_vect <- 
                  dnaA_vect %>% dplyr::slice(1)
              }
            }
            
            if (nrow(dnaA_vect)==0) {
              out_string <- 
                sprintf(
                  "%s has no dnaA annotation",
                  genome_acc
                )
            } else {
              gff_ori_strand <- 
                pull(dnaA_vect, strand)
              ## syncronize to the dnaA position, reverse comp if necessary
              if (gff_ori_strand=="+"){
                ori_start <- 
                  pull(dnaA_vect, start)
                out_seq <- 
                  xscat(
                    subseq(
                      genome_seq, 
                      start = ori_start, end = genome_length),
                    subseq(
                      genome_seq, 
                      start = 1, end = (ori_start-1))
                  ) %>%
                  `names<-`(genome_acc)
              } else {
                ori_start <- 
                  genome_length - pull(dnaA_vect, end) + 1
                out_seq <- 
                  xscat(
                    subseq(
                      reverseComplement(genome_seq), 
                      start = ori_start, end = genome_length
                    ),
                    subseq(
                      reverseComplement(genome_seq), 
                      start = 1, end = (ori_start-1)
                    )
                  ) %>%
                  `names<-`(genome_acc)
                
              }
              writeXStringSet(
                out_seq, out_file
              )
              out_string <- 
                sprintf(
                  "%s synced",
                  genome_acc
                )
            }
          }
          ## the dnaA locus
          #print(i)}
          return(out_string)
        }, mc.cores = core_number
      ) %>% unlist
    return(sync_report)
  }



#get_genomes(report_details, genome_dir = "unsynced", gff_dir)

kmer_counts <- 
  function(g_name, repeats_dir, sync_dir){
    #g_name = genome_vector[1]
    k_count_com <- 
      sprintf(
        "%s -in=%s/%s.txt -out=%s/%s_31mer.txt",
        "kmercountexact.sh fastadump=f rcomp=t mincount=2",
        sync_dir, g_name, repeats_dir, g_name
      )
    k_out <- 
      sprintf("%s/%s_31mer.txt", repeats_dir, g_name)
    if (!file.exists(k_out)){
      system(k_count_com, ignore.stderr = T)
    }
    if (file.size(k_out)==0){
      kmer_count = 0
    } else {
      kmer_count <- 
        fread(k_out) %>% #filter(V2>1) %>%
        pull(V2) %>% sum
    }
    
    return(kmer_count)
  }



mash_commands <- 
  function(genus, org, identity_dir, sync_dir){
    sketchfile <- 
      sprintf(
        "mash sketch -p 8 -o %s/mash.msh -s 10000 '%s'/*.txt ",
        identity_dir, sync_dir
      )
    mash_lines <- 
      sprintf(
        "mash dist -p 8 %s/mash.msh %s/mash.msh > %s/%s_%s_mash.txt",
        identity_dir, identity_dir, 
        identity_dir, genus, org
      )
    return(capture.output(cat(sketchfile, mash_lines, sep = "\n")))
  }

identity_sld <- 
  function(
    identity_dir, genus_i, org_i, report_details, filetype = "mash",
    clust_method = "complete", sync_dir
  ){
    
    id_table <-
      fread(
        sprintf(
          "%s/%s_%s_%s.txt", 
          identity_dir, genus_i, org_i, filetype
        ), stringsAsFactors = F, select = c(1:3)
      ) %>%
      mutate(V1 = sub(".txt", "", basename(V1))) %>%
      mutate(V2 = sub(".txt", "", basename(V2))) %>%
      `colnames<-`(c("ref", "qry", "identity")) %>%
      mutate(identity = as.numeric(identity)) %>%
      as_tibble
    
    if (filetype=="mash"){
      id_table <- 
        id_table %>%
        mutate(
          identity = (1-identity)*100
        )
    }
    
    id_table_wide <- 
      id_table %>%
      dplyr::select(c(ref, qry, identity)) %>%
      spread(., qry, identity) %>%
      as.data.frame() %>%
      column_to_rownames("ref")
    
    ## for clustering purposes
    id_dist <- 
      100-as.matrix(id_table_wide) 
    id_dist[upper.tri(id_dist, diag = T)]<-NA
    
    id_dist <- 
      as.dist(id_dist)
    id_clust <- 
      hclust(id_dist, method = clust_method)
    clust_order <-
      id_clust$labels[id_clust$order]
    
    id_table_long <- 
      id_table_wide %>%
      rownames_to_column("ref") %>%
      gather(key = qry,value = "identity", na.rm = FALSE, -c(ref)) %>%
      group_by(ref) %>% 
      mutate(identity_mean = mean(identity)) %>%
      arrange(desc(identity_mean, ref, qry, identity)) %>%
      ungroup
    
    seq_lens <- 
      report_details %>% 
      dplyr::select(c(accession, len))
    
    if (any(is.na(pull(seq_lens, len)))){
      
    }
    if (
      any(
        !(
          (pull(id_table_long, ref) %>% unique) %in% 
          pull(seq_lens, accession)
        )
      )
    ){
      accs <- 
        (pull(id_table_long, ref) %>% unique)[
          !(
            (pull(id_table_long, ref) %>% unique) %in% 
              pull(seq_lens, accession)
          )
          ]
      accs_tib <- 
        lapply(accs, function(x){
          tibble(
            accession = x,
            len = width(readDNAStringSet(sprintf("%s/%s.txt", sync_dir, x)))
          )
        }) %>% bind_rows
      seq_lens <- 
        seq_lens %>%
        bind_rows(accs_tib)
    }
    
    
    
    id_sld <- 
      id_table_long[,c(1:3)] %>%
      left_join(
        ., rename(seq_lens, accession = "ref"), 
        by = "ref") %>%
      rename(len = "reflen") %>%
      left_join(
        ., rename(seq_lens, accession = "qry"), 
        by = "qry") %>%
      rename(len = "qrylen") %>%
      mutate(
        ld = abs(reflen-qrylen),
        sld = 100*(abs(reflen-qrylen)/mean(c(reflen, qrylen))),
        identity_sld = identity-sld
      ) %>%
      mutate(
        ref = factor(ref, levels = clust_order),
        qry = factor(qry, levels = clust_order)
      )
    return(id_sld)
  }

nucmer_maxmatch_commmand <- 
  function(ref, qry, outfile){
    sprintf(
      "nucmer --maxmatch --prefix=%s %s %s",
      outfile, ref, qry
      )
  }

generate_nucmer_commands <- 
  function(genome_matrix, nucmer_dir, sync_dir){
    nuc_commands <- 
      apply(
        genome_matrix,
        1,
        function(x){
          #ref <- grep(x[1], synced_files, value = T)
          #qry <- grep(x[2], synced_files, x[2], value = T)
          
          nuc_params <-
            "--mum --maxgap=500 --mincluster=100"
          nucmer <-
            sprintf(
              "nucmer %s --prefix=%s/%s_v_%s",
              nuc_params,
              nucmer_dir,
              x[1], x[2]
            )
          nucmer_command <- 
            sprintf(
              "then %s %s/%s.txt %s/%s.txt",
              nucmer,
              sync_dir, x[1], 
              sync_dir, x[2]
            )
          filter_out <- 
            sprintf(
              "%s/%s_v_%s_filtered.delta",
              nucmer_dir,
              x[1], 
              x[2]
            )
          
          nuc_check <- 
            sprintf(
              "if [ ! -e %s ]",
              filter_out
            )
          
          nucmer_filter <-
            sprintf(
              "delta-filter -r -q %s/%s_v_%s.delta > %s",
              nucmer_dir,
              x[1], x[2],
              filter_out
            )
          ## to save space, remove the original delta file
          nuc_remove <- 
            sprintf(
              "rm %s/%s_v_%s.delta",
              nucmer_dir,
              x[1], x[2]
            )
          
          return(
            capture.output(
              cat(
                nuc_check,
                nucmer_command, 
                nucmer_filter, 
                nuc_remove,
                "fi",
                sep = "; "
              )
            )
          )
          
        }
      )
    return(nuc_commands)
    
  }

one_vs_all <- 
  function(genome_vector, ref = "random", sync_dir, nucmer_dir){
    if(ref=="random") {
      set.seed(1234)
      ref <- 
        genome_vector[sample(length(genome_vector), 1)]
    } 
    pmat <- 
      tibble(
        ref = ref,
        qry = genome_vector[genome_vector!=ref]
      ) %>% 
      as.matrix()
    
    message(sprintf("Aligning 1 vs %i sequences", nrow(pmat)))
    ##
    nuc_params <-
      "--mum --maxgap=500 --mincluster=100"
    
    if (!dir.exists(nucmer_dir)){dir.create(nucmer_dir, recursive = T)}
    
    nuc_commands <- 
      generate_nucmer_commands(pmat, nucmer_dir, sync_dir)
    
    ## aligning all sequences
    pbmclapply(nuc_commands, system, mc.cores = 7)
  }

all_vs_all <- 
  function(genome_tibble, sync_dir, nucmer_dir){
    # ani2$ref
    # ani2$qry
    pmat <- 
      tibble(
        ref = genome_tibble$ref,
        qry = genome_tibble$qry
      ) %>% 
      as.matrix()
    
    ##
    nuc_params <-
      "--mum --maxgap=500 --mincluster=100"
    
    if (!dir.exists(nucmer_dir)){dir.create(nucmer_dir, recursive = T)}
    
    nuc_commands <- 
      generate_nucmer_commands(pmat, nucmer_dir, sync_dir)
    
    ## aligning all sequences
    pbmclapply(nuc_commands, system, mc.cores = 7)
  }

read_delta <- 
  function(delta_path){
    delta_lines <- 
      (readLines(delta_path) %>% strsplit(., " "))[-c(1,2)]
    id_lines <- 
      delta_lines[lengths(delta_lines)==4][[1]]
    delta_alignments <- 
      delta_lines[lengths(delta_lines)==7] %>%
      unlist %>%
      matrix(., ncol = 7, byrow = T) 
    
    delta_alignments <- 
      delta_alignments %>%
      apply(., 2, as.numeric) %>%
      matrix(., ncol = 7) %>%
      `colnames<-`(c("rs", "re", "qs", "qe", "error", "e2", "zero")) %>% 
      as_tibble %>%
      dplyr::select(c(1:5)) %>%
      mutate(
        strand = ifelse(qe-qs > 0, '+', '-'),
        rid = strsplit(id_lines, " ")[[1]] %>% sub(">", "", .),
        qid = strsplit(id_lines, " ")[[2]],
        rlen = strsplit(id_lines, " ")[[3]] %>% as.numeric,
        qlen = strsplit(id_lines, " ")[[4]] %>% as.numeric,
        rcov = abs(re-rs+1),
        qcov = abs(qe-qs+1),
        perc_error = (round(100*error/ pmax(rcov, qcov), 2))
      ) %>% 
      mutate(
        meanlen = ceiling(rowMeans(. [, c("rcov", "qcov")])),
        xmid = rlen/2,
        X_dist = 
          ifelse(
            strand=="+",
            abs((rs-qs)/sqrt(2)),
            abs(qs + rs - qlen)/sqrt(2)
          ),
        X_dist_weight = 
          X_dist * (meanlen/mean(c(rlen, qlen)))
      )
    
    return(delta_alignments)
  }

colinear_delta <- 
  function(delta_table){
    #   delta_tibbles[[1]] # %>% plot_delta
    if (is.character(delta_table)){
      delta_table <- read_delta(delta_table)
    }
    
    delta_table <- 
      delta_table %>%
      mutate(
        qs2 = ifelse(strand=="+", qs, qe),
        qe2 = ifelse(strand=="+", qe, qs)
      )
    
    ref <- 
      delta_table %>%
      pull(rid) %>% 
      unique
    qry <- 
      delta_table %>%
      pull(qid) %>%
      unique
    ## seqlens
    rlen <- 
      delta_table %>%
      pull(rlen) %>%
      unique
    qlen <- 
      delta_table %>%
      pull(qlen) %>%
      unique
    
    ref_cov <- 
      IRanges(delta_table$rs, delta_table$re) %>% 
      reduce() %>%
      width %>% sum 
    qry_cov <- 
      IRanges(delta_table$qs2, delta_table$qe2) %>% 
      reduce() %>%
      width %>% sum 
    cov_prop <- 
      round(
        mean(c(ref_cov, qry_cov)) / 
          mean(c(rlen, qlen)), 
        3
      )
    if (cov_prop<0.8){
      out_logical = FALSE
    } else {
      fd <- 
        delta_table %>% 
        filter_delta(., minlen = 5e4)
      out_logical <- 
        ifelse(
          nrow(filter(fd, strand=="-"))==0,
          TRUE, FALSE
        )
    }
    return(out_logical)
    
  }

generate_nucmer_commands <- 
  function(genome_matrix, nucmer_dir, sync_dir){
    nuc_commands <- 
      apply(
        genome_matrix,
        1,
        function(x){
          #ref <- grep(x[1], synced_files, value = T)
          #qry <- grep(x[2], synced_files, x[2], value = T)
          
          nuc_params <-
            "--mum --maxgap=500 --mincluster=100"
          nucmer <-
            sprintf(
              "nucmer %s --prefix=%s/%s_v_%s",
              nuc_params,
              nucmer_dir,
              x[1], x[2]
            )
          nucmer_command <- 
            sprintf(
              "then %s %s/%s.txt %s/%s.txt",
              nucmer,
              sync_dir, x[1], 
              sync_dir, x[2]
            )
          filter_out <- 
            sprintf(
              "%s/%s_v_%s_filtered.delta",
              nucmer_dir,
              x[1], 
              x[2]
            )
          
          nuc_check <- 
            sprintf(
              "if [ ! -e %s ]",
              filter_out
            )
          
          nucmer_filter <-
            sprintf(
              "delta-filter -r -q %s/%s_v_%s.delta > %s",
              nucmer_dir,
              x[1], x[2],
              filter_out
            )
          ## to save space, remove the original delta file
          nuc_remove <- 
            sprintf(
              "rm %s/%s_v_%s.delta",
              nucmer_dir,
              x[1], x[2]
            )
          
          return(
            capture.output(
              cat(
                nuc_check,
                nucmer_command, 
                nucmer_filter, 
                nuc_remove,
                "fi",
                sep = "; "
              )
            )
          )
          
        }
      )
    return(nuc_commands)
    
  }

filter_delta <- 
  function(
    delta_table, contig_summary = FALSE, 
    maxgap = 1e4, minlen = 1e4, X_dist_diff = 5e4
  ){
    if (is.character(delta_table)){
      delta_table <- read_delta(delta_table)
    }
    
    out_tibble <- 
      delta_table %>%
      #dtf %>%
      group_by(strand) %>% 
      mutate(
        qry_gapsize = 
          qs - lag(qe, default = qs[1]),
        ref_gapsize = 
          rs - lag(re, default = rs[1]),
        X_diff = 
          X_dist-lag(X_dist, default = 0)
      ) %>%
      ungroup() %>%
      mutate(
        qry_gapsize = 
          ifelse(
            strand=="+",
            qry_gapsize,
            qry_gapsize*-1
          ),
        qry_gaps_up = 
          ifelse(
            (qry_gapsize) < maxgap,
            0,1
          ),
        qry_gaps_up = 
          cumsum(qry_gaps_up),
        qry_gaps_down = 
          ifelse(
            (qry_gapsize) > -maxgap,
            0,1
          ),
        qry_gaps_down = 
          cumsum(qry_gaps_down),
        ref_gaps = 
          ifelse(
            abs(ref_gapsize) < maxgap,
            0,1
          ),
        ref_gaps = 
          cumsum(ref_gaps),
        XDD = ## new
          X_dist-lag(X_dist, default = 0),
        XDF = 
          ifelse(
            abs(XDD) < X_dist_diff,
            0,1
          ),
        XDF_diff = 
          cumsum(XDF),
        new_contigs = ## new
          rleid(strand, qry_gaps_up, qry_gaps_down, ref_gaps, XDF_diff)
      ) %>%
      group_by(new_contigs, strand) %>% 
      dplyr::summarise(
        "X_dist" = mean(X_dist),
        #weighted.mean(x = X_dist, w =  meanlen/max(re)),
        "meanlen" = sum(meanlen),
        "rs" = min(rs),
        "re" = max(re),
        "qs" = 
          unique(ifelse(
            strand=="+",
            min(qs),
            max(qs)
          )),
        "qe" = 
          unique(ifelse(
            strand=="+",
            max(qe),
            min(qe)
          )),
        #"xmid" = unique(xmid),
        "rid" = unique(rid),
        "qid" = unique(qid),
        "slope" = (qe-qs)/(re-rs),
        "rlen" = unique(rlen),
        "qlen" = unique(qlen), 
        .groups = "keep"
      ) %>%
      ungroup %>%
      filter(
        meanlen > minlen
      ) %>%
      #group_by(strand) %>% 
      mutate(
        qry_gapsize = 
          qs - lag(qe, default = qs[1]),
        ref_gapsize = 
          rs - lag(re, default = rs[1]),
        XDD = 
          X_dist-lag(X_dist, default = 0),
        XDF = 
          ifelse(
            abs(XDD) < X_dist_diff,
            0,1
          ),
        XDF_diff = 
          cumsum(XDF),
        new_contigs = 
          rleid(strand, XDF_diff)
      ) %>%
      ungroup() %>%
      group_by(new_contigs, strand) %>% 
      dplyr::summarise(
        "X_dist" =  mean(X_dist),
        #weighted.mean(X_dist, meanlen),
        # "meanlen" = sum(meanlen),
        "rs" = min(rs),
        "re" = max(re),
        "qs" = 
          unique(ifelse(
            strand=="+",
            min(qs),
            max(qs)
          )),
        "qe" = 
          unique(ifelse(
            strand=="+",
            max(qe),
            min(qe)
          )),
        #"xmid" = unique(xmid),
        "rid" = unique(rid),
        "qid" = unique(qid),
        "slope" = (qe-qs)/(re-rs),
        "rlen" = unique(rlen),
        "qlen" = unique(qlen), 
        "meanlen" = mean(c(abs(qe-qs), abs(re-rs))),
        .groups = "keep"
      ) %>% 
      ungroup
    
    if (contig_summary){
      out_tibble <- 
        out_tibble %>%
        mutate(
          c_contig = rleid(strand)
        ) %>%
        #mutate(refmax = max(re)) %>%
        group_by(c_contig, strand) %>%
        dplyr::summarise(
          "rs" = min(rs),
          "re" = max(re),
          "qs" = 
            unique(ifelse(
              strand=="+",
              min(qs),
              max(qs)
            )),
          "qe" = 
            unique(ifelse(
              strand=="+",
              max(qe),
              min(qe)
            )
            ),
          "rlen" = unique(rlen),
          "qlen" = unique(qlen),
          "slope" = mean(slope),
          "X_dist" = mean(X_dist),
          "rid" = unique(rid),
          "qid" = unique(qid), 
          .groups = "keep"
          #"xmid" = unique(xmid),
          #"c_length" = sum(total_length),
          #"X_distance" = mean(X_distance),
          #refmax = unique(refmax)
        ) %>%
        ungroup()
    }
    
    return(out_tibble)
  }

Cluster_Genomes <- 
  function(ref, genome_vector, report_details, nucmer_dir, r_backup){
    clustering_results <- 
      list()
    j=1
    remaining_genomes <- genome_vector
    while (length(remaining_genomes) >0 ){
      iteration = 
        j
      if (length(remaining_genomes)==1){
        clustering_results[[j]] <- remaining_genomes
        break
      } else if ( length(remaining_genomes)==0) { break }
      
      message(sprintf("Iteration %i", iteration))
      
      set.seed(1234)
      sampled_seq <- 
        sample(length(remaining_genomes), 1)
      if (j==1){
        sampled_seq <- 
          which(remaining_genomes==ref)
      }
      pmat <- 
        tibble(
          ref = remaining_genomes[sampled_seq],
          qry = remaining_genomes[-sampled_seq]
        ) %>% 
        as.matrix()
      
      message(sprintf("Aligning 1 vs %i sequences", nrow(pmat)))
      
      ##
      nuc_params <-
        "--mum --maxgap=500 --mincluster=100"
      iter_nuc_dir <- #"nucmer/s_enterica/test"
        sprintf(
          "%s/i%i_1",
          nucmer_dir, 
          iteration
        )
      if (!dir.exists(iter_nuc_dir)){dir.create(iter_nuc_dir, recursive = T)}
      
      nuc_commands <- 
        generate_nucmer_commands(pmat, iter_nuc_dir, sync_dir)
      
      ## aligning all sequences
      # for (nc in nuc_commands){
      #   system(nc)
      # }
      
      pbmclapply(nuc_commands, system, mc.cores = 7)
      
      # Clustering
      message(
        sprintf(
          "Clustering..."
        )
      )
      col_tib <- 
        pbmclapply(
          1:nrow(pmat), 
          function(i) {
            #for(i in 1:nrow(pmat)){
            if(
              (sprintf(
                "%s/%s_v_%s_filtered.delta",
                iter_nuc_dir, pmat[i,1], pmat[i,2]
              ) %>% read_delta() %>% nrow)==0
            ){
              out_tib <- 
                tibble(
                  ref = pmat[i,1],
                  qry = pmat[i,2],
                  colinear = FALSE
                )
            } else{
              out_tib <-
                tibble(
                  ref = pmat[i,1],
                  qry = pmat[i,2],
                  colinear = 
                    colinear_delta(
                      read_delta(
                        sprintf(
                          "%s/%s_v_%s_filtered.delta",
                          iter_nuc_dir, pmat[i,1], pmat[i,2]
                        )
                      ) 
                    )
                )
            }
            
            return(out_tib)
          }, mc.cores = 7
        ) %>% 
        bind_rows()
      
      colinear_genomes <- 
        col_tib %>%
        filter(colinear) %>%
        pull(qry)
      #plot_delta(paste0(iter_nuc_dir, "/NZ_CP019411.1_NZ_CP006693.1_filtered.delta"))
      clustering_results[[j]] <- 
        c(remaining_genomes[sampled_seq], colinear_genomes)
      
      remaining_genomes <- 
        remaining_genomes[!(remaining_genomes %in% clustering_results[[j]])]
      
      message(
        sprintf(
          "%i sequences clustered",
          length(clustering_results[[j]])
        )
      )
      # backup_dir <- 
      #   sprintf("%s/clus", r_backup, org)
      # if(!dir.exists(backup_dir)){dir.create(backup_dir, recursive = T)}
      saveRDS(
        clustering_results, 
        sprintf(sprintf("%s/clustering_results_%i.rds", r_backup, j))
      )
      j=j+1
      
    }
    ## order high to low
    clustering_results_final <- 
      clustering_results[order(sapply(clustering_results,length),decreasing=T)]
    saveRDS(
      clustering_results_final, 
      sprintf("%s/clustering_results.rds", r_backup)
    )
    return(clustering_results_final)
  }

low_cov_or_inverted <- 
  function(delta_table){
    # delta_table <-
    #   delta_tibbles[[1]] # %>% plot_delta
    if (is.character(delta_table)){
      delta_table <- read_delta(delta_table)
    }
    
    ref <- 
      delta_table %>%
      pull(rid) %>% 
      unique
    qry <- 
      delta_table %>%
      pull(qid) %>%
      unique
    ## seqlens
    rlen <- 
      delta_table %>%
      pull(rlen) %>%
      unique
    qlen <- 
      delta_table %>%
      pull(qlen) %>%
      unique
    
    ref_cov <- 
      IRanges(delta_table$rs, delta_table$re) %>% 
      reduce() %>%
      width %>% sum 
    ref_cov_prop <- 
      round(ref_cov/rlen, 3)
    fd <- 
      delta_table %>% 
      filter_delta()#., minlen = 5e4)
    
    lci <- 
      tibble(
        low_coverage = (ref_cov_prop<0.8),
        inversions = nrow(filter(fd, strand=="-"))!=0
      )
    return(lci)
  }
LETTERS702 <- c(LETTERS, sapply(LETTERS, function(x) paste0(x, LETTERS)))

get_cluster_table <- 
  function(clustered, mash_tib, sync_dir, nucmer_dir){
    cluster_tib <- 
      lapply(
        1:length(clustered),
        function(i){
          tibble(
            accession = clustered[[i]],
            clust = i
          )
        }
      ) %>%
      bind_rows() %>%
      add_count(clust, name = "clust_count") %>%
      mutate(
        clust = LETTERS702[clust]
      )
    
    clust_identity <- 
      inner_join(
        mash_tib,
        (cluster_tib[,c(1:2)] %>% `colnames<-`(c("ref", "ref_clust"))),
        by = "ref"
      ) %>%
      inner_join(
        ., 
        (cluster_tib[,c(1:2)] %>% `colnames<-`(c("qry", "qry_clust"))), 
        by = "qry"
      ) 
    clust_identity_diff <- 
      clust_identity %>%
      filter(
        ref_clust != qry_clust, 
        ref_clust !="LOC", 
        qry_clust !="LOC"
      ) %>%
      group_by(ref_clust, qry_clust) %>%
      filter(identity_sld==max(identity_sld)) %>% 
      arrange(
        match(ref_clust, c(LETTERS702, "LOC")), 
        match(qry_clust,  c(LETTERS702, "LOC"))
      ) %>%
      mutate(
        clusts =
          paste(
            sort(
              c(as.character(ref_clust)[1], as.character(qry_clust)[1])
            ),
            collapse = ""
          )
      ) %>%
      ungroup %>%
      group_by(clusts) %>%
      filter(identity_sld==max(identity_sld)) %>%
      #arrange(clusts) %>%
      filter(!duplicated(clusts)) %>%
      arrange(
        match(ref_clust, c(LETTERS702, "LOC")), 
        match(qry_clust,  c(LETTERS702, "LOC"))
      )
    # arrange(ref_clust, qry_clust)
    
    ava <- 
      all_vs_all(
        clust_identity_diff[,1:2], 
        sync_dir, 
        nucmer_dir = sprintf("%s/ava", nucmer_dir)
      )
    
    zero_size <- 
      which(
        lapply(
          list.files(sprintf("%s/ava", nucmer_dir), full.names = T),
          function(fs){
            file.size(fs)==0
          }
        ) %>% unlist
      )
    
    
    if (length(zero_size)>0){
      zfiles <- 
        list.files(sprintf("%s/ava", nucmer_dir), full.names = T)[zero_size]
      file.remove(zfiles)
      redone_command <- 
        generate_nucmer_commands(
          clust_identity_diff[zero_size,1:2], nucmer_dir, sync_dir
        ) 
      lapply(redone_command, system)
      
      ava <- 
        all_vs_all(
          clust_identity_diff[,1:2], 
          sync_dir, 
          nucmer_dir = sprintf("%s/ava", nucmer_dir)
        )
    }
    
    
    delta_tibbles <- 
      pbmclapply(
        1:nrow(clust_identity_diff),
        function(ij){
          #for (i in 1:nrow(clust_identity_diff)){
          delta_file <- 
            (list.files(
              sprintf("%s/ava", nucmer_dir), 
              full.names = T
            ) %>% 
              grep(
                pattern = 
                  paste0(pull(clust_identity_diff, ref)[ij], "_"), 
                value = T
              ) %>% 
              grep(
                pattern = 
                  paste0(pull(clust_identity_diff, qry)[ij], "_"), 
                value = T
              ))[1] 
          if (file.size(delta_file)==0){
            out_tib <- NULL
          } else {
            out_tib <- 
              read_delta(delta_file) #%>% low_cov_or_inverted()
          }
          return(out_tib)
          
        }
      )
    
    # delta_tibbles = delta_tibbles[-which(sapply(delta_tibbles, is.null))]
    # which(lapply(delta_tibbles, is.null) %>% unlist)
    #for (ij in 1:length(delta_tibbles)){low_cov_or_inverted(delta_tibbles[[ij]])}
    #delta_tibbles
    
    clust_identity_diff_2 <- 
      clust_identity_diff %>%
      bind_cols(
        pbmclapply(delta_tibbles, low_cov_or_inverted, mc.cores = 1) %>% bind_rows()
      ) %>%
      arrange(match(ref_clust, LETTERS702), match(qry_clust, LETTERS702)) %>%
      mutate(
        old_clust =
          qry_clust,
        #substr(clusts, start = 2, stop = 2),
        LOC = 
          (!inversions )#& low_coverage),
        
      )
    #clust_identity_diff_2 %>% filter(LOC)
    
    #for (dt in delta_tibbles){low_cov_or_inverted(dt)}
    #clust_identity_diff_2 %>% filter(!low_coverage & !inversions)
    if (any(clust_identity_diff_2$LOC)){
      cluster_tib <- 
        cluster_tib %>%
        mutate(
          clust2 = 
            ifelse(
              (clust %in% filter(clust_identity_diff_2, LOC)$old_clust),
              "LOC",
              clust
            )
        ) %>% 
        arrange(match(clust2,  c(LETTERS702, "LOC"))) %>%
        mutate(
          clustrle = 
            rleid(clust2), 
          clust = 
            ifelse(
              clust2=="LOC",
              "LOC",
              LETTERS702[clustrle]
            )
        ) %>% 
        dplyr::select(c(1:2)) %>%
        add_count(clust, name = "clust_count")
    } 
    clust_identity <- 
      suppressWarnings(
        inner_join(
          mash_tib,
          (cluster_tib[,c(1:2)] %>% `colnames<-`(c("ref", "ref_clust"))),
          by = "ref"
        ) %>%
          inner_join(
            ., 
            (cluster_tib[,c(1:2)] %>% `colnames<-`(c("qry", "qry_clust"))), 
            by = "qry"
          ) %>%
          mutate(
            ref_clust = 
              factor(
                ref_clust, 
                levels = 
                  c(LETTERS702, "LOC")
              ),
            qry_clust = 
              factor(
                qry_clust, 
                levels = levels(ref_clust)
              )
          ) %>%
          arrange(
            match(ref_clust, c(LETTERS702, "LOC")), 
            match(qry_clust,  c(LETTERS702, "LOC"))
          )
      )
    
    clust_identity_diff <- 
      clust_identity %>%
      filter(
        ref_clust != qry_clust, 
        ref_clust !="LOC", 
        qry_clust !="LOC"
      ) %>%
      group_by(ref_clust, qry_clust) %>%
      filter(identity_sld==max(identity_sld)) %>% 
      mutate(
        clusts = 
          paste(
            sort(
              c(as.character(ref_clust)[1], as.character(qry_clust)[1])
            ), 
            collapse = ""
          )
      ) %>%
      ungroup %>%
      group_by(clusts) %>%
      filter(identity_sld==max(identity_sld)) %>% 
      arrange(clusts) %>% 
      filter(!duplicated(clusts)) %>%
      arrange(
        match(ref_clust, c(LETTERS702, "LOC")), 
        match(qry_clust,  c(LETTERS702, "LOC"))
      )
    
    ava <- 
      all_vs_all(
        clust_identity_diff[,1:2], 
        sync_dir, 
        nucmer_dir = sprintf("%s/ava", nucmer_dir)
      )
    
    #clust_identity_diff_2
    return(cluster_tib)
    
  }


delta_midpoints <- 
  function(delta_table){
    if (is.character(delta_table)){
      delta_table <- read_delta(delta_table)
      #delta_table <- delta_tibbles[[1]]
    }
    ## report midpoint if the inversion is at least 1 % of the length
    filtered_table <- 
      delta_table %>%
      filter_delta() %>% 
      filter(
        meanlen/mean(c(delta_table$rlen[1], delta_table$qlen[1])) > 0.01
      ) %>%
      filter(strand =="-") %>% 
      mutate(
        ref_mp = ceiling((rs+re)/2),
        qry_mp = ceiling((qs+qe)/2)
      ) 
    out_tib <-
      dplyr::select(filtered_table, ref_mp, qry_mp, new_contigs) %>%
      group_by(new_contigs) %>%
      mutate(
        mp_mean = mean(c(ref_mp, qry_mp))
      ) %>% 
      ungroup %>%
      mutate(
        mean_len =
          mean(c(pull(delta_table, rlen)[1], pull(delta_table, qlen)[1])),
        ref_prop = ref_mp/mean_len,
        qry_prop = qry_mp/mean_len,
        mean_prop = mp_mean/mean_len
      )
    return(out_tib)
  }


plot_delta <- 
  function(
    delta_table, gtitle = "", 
    xlb = NULL, ylb = NULL,
    add_lines = TRUE
  ){
    
    if (is.character(delta_table)){
      delta_table <- read_delta(delta_table)
    }
    #delta_table <- delta_tibbles[[1]]
    myColors <- brewer.pal(5,"Set1")[c(1,2)]
    names(myColors) <- c("-", "+")
    colScale <- scale_colour_manual(name = "grp",values = myColors)
    
    p <- 
      ggplot(
        delta_table, 
        aes(
          x=rs, xend=re, y=qs, yend=qe, 
          colour= factor(strand, levels = c("-", "+"))
        ) 
      )+
      geom_segment(alpha=1, size = 1) + 
      #geom_point(alpha=.5, size = 1) + 
      theme_classic() + 
      #facet_grid(qid~., scales='free', space='free', switch='y') +
      theme(
        legend.position= "none",#c(.99,.4), 
        legend.justification=c(1,0),
        #axis.text.y=element_blank(), 
        #axis.ticks.y=element_blank(),
        plot.title = element_text(hjust = 0.5)
      ) +
      ggtitle(gtitle) +
      # geom_vline(
      #     xintercept = unique(pull(delta_table, xmid)),
      #     color="darkseagreen2",
      #     linetype="dashed", size=.25) +
      # geom_abline(
      #     intercept = 
      #         2*unique(pull(delta_table, xmid)), 
      #     slope = -1, 
      #     color="darkseagreen2",
      #     linetype="dashed", size=.25) +
      # geom_abline(
    #   intercept = 0, slope = 1, color="darkseagreen2",
    #   linetype="dashed", size=.5) +
    colScale +
      scale_x_continuous(limits=c(1, delta_table$rlen[1])) +
      scale_y_continuous(limits=c(1, delta_table$qlen[1]))
    
    #scale_colour_brewer(brewer.pal(n = 8, name = "Set1")[c(1,2)])
    #scale_colour_brewer(palette='Set1')
    
    if (is.null(xlb) & is.null(ylb)){
      p <- 
        p +             
        xlab(unique(pull(delta_table, rid))) + 
        ylab(unique(pull(delta_table, qid))) 
    } else{
      p <-
        p +
        xlab(xlb) + 
        ylab(ylb) 
    }
    if (add_lines){
      p <- 
        p +
        geom_abline(
          intercept = pull(delta_tibbles[[1]], qlen)[1],
          slope = -1, color="grey25",
          linetype="dashed", size=.1) +
        geom_abline(
          intercept = 0,
          slope = 1, color="grey25",
          linetype="dashed", size=.1)
    }
    return(p)
  }


plot_delta5 <- 
  function(
    delta_table, gtitle = "", 
    xlb = NULL, ylb = NULL,
    add_lines = TRUE
  ){
    
    if (is.character(delta_table)){
      delta_table <- read_delta(delta_table)
      
    }
    myColors <- brewer.pal(5,"Set1")[c(1,2)]
    names(myColors) <- c("-", "+")
    colScale <- scale_colour_manual(name = "grp",values = myColors)
    fdt <-
      filter_delta(delta_table)
    
    p <- 
      ggplot(
        delta_table, 
        aes(
          x=rs, xend=re, y=qs, yend=qe, 
          colour= factor(strand, levels = c("-", "+"))
        ) 
      )+
      geom_segment(alpha=1, size = 2) + 
      theme_classic() + 
      theme(
        legend.position= "none",#c(.99,.4), 
        legend.justification=c(1,0),
        axis.text=element_blank(), 
        axis.ticks=element_blank(),
        plot.title = element_text(hjust = 0.35),
        axis.title = element_text(size=20)
      ) +

    geom_hline(
      yintercept = fdt$qe[1], 
      color=(brewer.pal(6, "Set2"))[2],
      linetype="dashed", size=.35
    ) +
      geom_hline(
        yintercept = fdt$qs[2], 
        color=(brewer.pal(6, "Set2"))[2],
        linetype="dashed", size=.35
      ) +
      geom_vline(
        xintercept = fdt$re[1], 
        color=(brewer.pal(6, "Set2"))[2],
        linetype="dashed", size=.35
      ) +
      geom_vline(
        xintercept = fdt$re[2], 
        color=(brewer.pal(6, "Set2"))[2],
        linetype="dashed", size=.35
      )+
      geom_abline(
        intercept =
          unique(pull(delta_table, rlen)),
        slope = 
          -unique(pull(delta_table, qlen))/(unique(pull(delta_table, rlen)))
        ,
        color="darkgrey",
        linetype="dashed", size=.35) +
      geom_abline(
        intercept = 0, slope = 1, 
        color="darkgrey",
        linetype="dashed", size=.35) +
      colScale +
      scale_x_continuous(limits=c(1, delta_table$rlen[1]), expand = c(0,0)) +
      scale_y_continuous(limits=c(1, delta_table$qlen[1]), expand = c(0,0))
    p
    #scale_colour_brewer(brewer.pal(n = 8, name = "Set1")[c(1,2)])
    #scale_colour_brewer(palette='Set1')
    
    if (is.null(xlb) & is.null(ylb)){
      p <- 
        p +             
        xlab(unique(pull(delta_table, rid))) + 
        ylab(unique(pull(delta_table, qid))) 
    } else{
      p <-
        p +
        xlab(xlb) + 
        ylab(ylb) 
    }
    return(p)
  }
