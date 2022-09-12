#Full clustering pipeline


species_list <- 
  ncbi_table_clade_count %>% 
  ungroup %>%
  arrange(desc(spec_count))


# species_name <-
#   "Yersinia_pestis"

species_name <-
  species_list$species2[11]
finished_species <- 
  list.dirs("results")

which(species_list$species2==species_name)
#67
#80
#45
spec2_list
#for (species_name in rev(species_list$species2)){
for (species_name in spec2_list[75]){
  #tryCatch({
    ####### directory stuff
    genus_name <- 
      strsplit(species_name, "_")[[1]][1]
    
    clade_name <- 
      filter(ncbi_table_clade_count, species2==species_name) %>%
      pull(p2) %>% gsub(" ", "_", .)
    
    
    species_variable_names <-
      c(
        "sync_dir", "gff_dir", "identity_dir", 
        "nucmer_dir", "r_vars_dir", "results_dir", 
        "repeats_dir", "invs_dir"
      )
    
    species_dirs <- 
      gsub("_dir", "", species_variable_names)
    
    for (var_name in species_variable_names){
      dir_name <- 
        gsub("_dir", "", var_name)
      full_dir_name <- 
        sprintf("%s/%s/%s/%s", dir_name, clade_name, genus_name, species_name)
      
      assign(
        var_name, 
        full_dir_name
      )
      if (!dir.exists(full_dir_name)){dir.create(full_dir_name, recursive = T)}
    }
    
    ############# details from the original table
    
    species_table <- 
      ncbi_table %>%
      filter(species== sub("_", " ", species_name))
      #filter(name== sub("_", " ", species_name))
    
    ## Get additional information from the assembly report
    report_details <-
      get_report_table(species_table) %>%
      mutate(
        len = as.numeric(len),
        accession = make.unique(accession)
      ) %>% 
      distinct %>%
      dplyr::select(-c(submitter))

    ## Filter for contigs 
    if (nrow(report_details)>nrow(species_table)){
      report_details <- 
        report_details %>%
        group_by(asm_name) %>%
        filter(!is.na(len)) %>%
        filter(len == max(len, na.rm = T)) %>% ungroup()
    }
    
    report_details <- 
      report_details %>% #pull(len) %>% density %>% plot 
      filter(len < mean(len, na.rm = T)+1e6 & len > mean(len, na.rm = T)-1e6) 
    
    ## join the report table with the filtered table
    full_species_table <-   
      left_join(report_details, species_table, by = "asm_name") %>%
      mutate(len = as.numeric(len)) %>%
      filter(!duplicated(asm_name)) %>%
      dplyr::select(
        -c(
          submitter, refseq_category, s, 
          gbrs_paired_asm, paired_asm_comp, species_taxid
        )
      )
    ###### Synchronization
    gff_table <- 
      download_gff(full_species_table, gff_dir)
    
    gffs_available <- 
      pull(gff_table, gff_exists) %>% sum
    
    sprintf(
      "%i of %i GFF records downloaded",
      gffs_available, nrow(gff_table)
    )
    
    ## Synchronize genomes using their dnaA annotations from the GFF file
    sync_report <-
      sync_with_gff(
        full_species_table,
        sync_dir, 
        gff_dir
      )
    synced_table <- 
      tibble(
        Identified_genomes = nrow(full_species_table),
        Synced_genomes = length(list.files(sync_dir))
      )
    synced_table
    #### Repetitive quantification
    accession_list <- 
      list.files(sync_dir) %>%
      sub(".txt", "", .)
    
    if (synced_table$Synced_genomes<10){
      write_csv(report_details, sprintf("%s/report_summary.csv", results_dir))
      next
    }
    
    ## Repetitive sequnece information
    possible_kmers <-
      pbmclapply(
        accession_list,
        function(x){
          #for (x in accession_list){
          filter(full_species_table, accession == x) %>%
            dplyr::select(accession, len) %>%
            mutate( 
              kmers_31 =
                kmer_counts(g_name = x, repeats_dir, sync_dir),
              possible_31_mers =
                len-31+1,
              kmer_prop =
                kmers_31/possible_31_mers
            )
        }, mc.cores = 6
      ) %>%
      bind_rows %>%
      arrange(kmer_prop) %>%
      mutate(kmer_perc = 100*kmer_prop)
    
    report_summary <- 
      left_join(report_details, possible_kmers[,c(1,6)], by = "accession") 
    
    
    ##### Pairwise identity comparisons
    if (!
        file.exists(
          sprintf(
            "%s/%s_%s_mash.txt", identity_dir, genus_name, species_name
          )
        )
    ){
      mash_bash <- 
        mash_commands(genus_name, species_name, identity_dir, sync_dir)
      system(mash_bash[1], ignore.stderr = T)
      system(mash_bash[2])
    }
    
    
    
    ## get the genome with the most closest relatives
    mash_tib <-
      identity_sld(
        identity_dir, genus_name, species_name, report_details, filetype = "mash", 
        sync_dir = sync_dir
      )
    most_related <- 
      (mash_tib %>% 
         mutate(ref = as.character(ref)) %>% 
         group_by(ref) %>% 
         summarise(mean_id = mean(identity_sld), .groups = "keep") %>%
         arrange(desc(mean_id)) %>%
         pull(ref))[1]
  


  cluster_file <- 
    sprintf("%s/clustering_results.rds", r_vars_dir)
  if (file.exists(cluster_file)){
    clustered <- 
      readRDS(cluster_file)
    if (length(unlist(clustered))!=length(accession_list)){
      clustered <- 
        Cluster_Genomes(
          ref = most_related, 
          genome_vector = accession_list, 
          report_details = report_details,
          nucmer_dir = nucmer_dir, 
          r_backup = r_vars_dir
        )
    } else {next}
    
  } else {
    clustered <- 
      Cluster_Genomes(
        ref = most_related, 
        genome_vector = accession_list, 
        report_details = report_details,
        nucmer_dir = nucmer_dir, 
        r_backup = r_vars_dir
      )
  }
  saveRDS(clustered, cluster_file)
  
  cluster_tib <-
    suppressWarnings(
      get_cluster_table(clustered, mash_tib, sync_dir, nucmer_dir) 
    )
  
  cluster_count <- 
    cluster_tib %>% 
    filter(clust !="LOC") %>% pull(clust) %>% unique %>% length
  
  if (cluster_count<2){
    write_csv(
      left_join(
        report_details, cluster_tib, by = "accession"
                ), 
      sprintf("%s/report_summary.csv", results_dir))
    next
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
                LETTERS702[1:(cluster_count)]
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
  
  ## within vs outside clust identity
  # clust_identity_within <- 
  #   clust_identity %>%
  #   filter(
  #     ref_clust !="LOC", 
  #     qry_clust !="LOC",
  #     ref != qry
  #     ) %>%
  #   group_by(ref_clust, qry_clust) %>%
  #   summarise(mean_sld = mean(identity_sld), comparisons = n())
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

  
  ### comparison tables between high identity representatives of clusts
  delta_tibbles <-
    pbmclapply(
      1:nrow(clust_identity_diff),
      function(i){
        #for (i in 1:nrow(clust_identity_diff)){
        (list.files(
          sprintf("%s/ava", nucmer_dir),
          full.names = T
        ) %>%
          grep(
            pattern =
              (pull(clust_identity_diff, ref)[i]), 
            value = T
          ) %>%
          grep(
            pattern =
              (pull(clust_identity_diff, qry)[i]), 
            value = T
          ))[1] %>%
          read_delta 
      }
    )
  all_midpoints <-
    pbmclapply(
      delta_tibbles,
      function(x){
        tibble(
          ref = pull(x, rid)[1],
          qry = pull(x, qid)[1],
          midpoint =
            delta_midpoints(x) %>% 
            pull(mean_prop)# %>%
          #mean
        )
      }, mc.cores = 6
    ) %>%
    bind_rows() %>%
    add_count(midpoint, name = "mpc") %>%
    mutate(
      mpc_prop = mpc/sum(mpc)
    ) %>%
    mutate(mid_dist = midpoint-0.5) 
  
  inversions_per_comparison <-
    pbmclapply(
      delta_tibbles,
      function(x){
        filter_delta(x) %>%
          group_by(rid, qid) %>%
          summarise(
            inversions = sum(strand=="-")
          ) %>%
          `colnames<-`(c("ref", "qry", "inversions"))
      }, mc.cores = 6
    ) %>%
    bind_rows()

  clust_identity_diff_2 <- 
    left_join(
      ungroup(clust_identity_diff),
      inversions_per_comparison,
      by = c("ref", "qry")
    )
  write_csv(
    clust_identity_diff_2, 
    sprintf("%s/clust_identity_diff.csv", results_dir)
  )  
  write_csv(
    inversions_per_comparison, 
    sprintf("%s/inversions_per_comparison.csv", results_dir)
  )
  
  
  inv_tib_long <- 
    lapply(
      1:nrow(clust_identity_diff),
      function(i){
        delta_invs <- 
          delta_tibbles[[i]] %>%
          filter_delta() %>%
          filter(meanlen >5e4) %>%
          mutate(inv_cs = rleid(strand)) %>%
          filter(
            strand =="-",
            !duplicated(inv_cs)
          ) %>% 
          nrow
        out_tib <- 
          clust_identity_diff[i,] %>% 
          ungroup %>% 
          dplyr::select(ref_clust, qry_clust) %>%
          mutate(inversions = delta_invs)
      }
    ) %>%
    bind_rows %>% 
    bind_rows(
      ., 
      tibble(
        ref_clust =  .$qry_clust,
        qry_clust = .$ref_clust,
        inversions = .$inversions
      )
    ) %>%
    mutate(
      inversions = 
        ifelse(
          ref_clust==qry_clust,
          0,
          ifelse(
            inversions==0,
            1,
            inversions
          )
        )
    )
  inv_tib_wide <- 
    inv_tib_long %>%
    spread(key = qry_clust, value = inversions) %>%
    replace(is.na(.), 0)
  
  write_csv(inv_tib_wide, sprintf("%s/inversion_table.csv", results_dir))
  report_summary <- 
    left_join(cluster_tib, report_summary, by = "accession") 
  
  write_csv(report_summary, sprintf("%s/report_summary.csv", results_dir))
  write_csv(inv_tib_wide,  sprintf("%s/inv_tib_wide.csv", results_dir))
  
  midpoint_summary <- 
    all_midpoints
  
  inv_adjacency <- 
    inv_tib_long %>%
    mutate(inversions = ifelse(inversions==1, 1, 0)) %>%
    group_by(ref_clust) %>%
    spread(key = qry_clust, value = inversions) %>%
    column_to_rownames(var = "ref_clust") %>% 
    as.matrix() %>%
    replace(is.na(.), 0)
  
  cluster_tib_count <- 
    cluster_tib %>%
    count(clust) %>%
    arrange(desc(n)) %>%
    mutate(size_v = rescale(n, to=c(15, 40))) %>%
    filter(clust != "LOC")
  
  inv_adjacency <- 
    graph_from_adjacency_matrix(inv_adjacency, mode = "undirected")
  V(inv_adjacency)$size = 
    cluster_tib_count$size_v
  png(
    filename = sprintf("%s/pathway.png", results_dir), 
    width = 300, height = 225, units='mm', res = 300
  )
  plot(inv_adjacency, vertex.color="lightblue", edge.width= 3)
  dev.off()
  
  saveRDS(all_midpoints, sprintf("%s/all_midpoints.rds", r_vars_dir))
  saveRDS(clust_identity_diff, sprintf("%s/clust_identity_diff.rds", r_vars_dir))
  
### saving images
  if (length(delta_tibbles)<200){
    delta_plots <- 
      pbmclapply(
        #which(grepl("A", pull(clust_ani_diff, clusts))),
        1:nrow(clust_identity_diff),
        function(i){
          #for (i in 1:nrow(clust_ani_diff)){
          #for (i in which(grepl("A", pull(clust_ani_diff, clusts)))){
          plot_delta(
            delta_tibbles[[i]]
          ) +
            xlab(pull(clust_identity_diff, ref_clust)[i]) +
            ylab(pull(clust_identity_diff, qry_clust)[i]) +
            theme(
              axis.text = element_blank(), 
              axis.ticks = element_blank()
            )
        }
      )
  }
  if (length(delta_tibbles)>50 & length(delta_tibbles)<200){
    plot_chunks <- 
      split(
        1:length(delta_plots), ceiling(seq_along(1:length(delta_plots))/25)
      )
    #pb = txtProgressBar(min = 1, max =length(plot_chunks), initial = 1)
    for (i in 1:length(plot_chunks)){
      ava_plots <- 
        do.call(
          "grid.arrange", 
          c(
            delta_plots[plot_chunks[[i]]], 
            ncol=5
          )
        )
      ggsave(sprintf("%s/ava_plots_%i.png", results_dir, i), ava_plots)
      #setTxtProgressBar(pb,i)
    }
    
  } else if (length(delta_tibbles)>=200) {
    delta_plots <- 
      pbmclapply(
        #which(grepl("A", pull(clust_ani_diff, clusts))),
        1:(cluster_count),
        function(i){
          plot_delta(
            delta_tibbles[[i]]
          ) +
            xlab(pull(clust_identity_diff, ref_clust)[i]) +
            ylab(pull(clust_identity_diff, qry_clust)[i]) +
            theme(
              axis.text = element_blank(), 
              axis.ticks = element_blank()
            )
        }
      )
    
    if (cluster_count >25){
      plot_chunks <- 
        split(
          1:cluster_count, ceiling(seq_along(1:cluster_count)/25)
        )
      for (i in 1:length(plot_chunks)){
        ava_plots <-
          do.call(
            "grid.arrange",
            c(
              delta_plots[plot_chunks[[i]]],
              ncol=5
            )
          )
        ggsave(sprintf("%s/ava_plots_%i.png", results_dir, i), ava_plots)
        #setTxtProgressBar(pb,i)
      }
      #pb = txtProgressBar(min = 1, max =length(plot_chunks), initial = 1)
    } else{
      ava_plots <- 
        do.call(
          "grid.arrange", 
          c(
            delta_plots,#[1:cluster_count], 
            ncol=5
          )
        )
      ggsave(sprintf("%s/ava_plots_A.png", results_dir, i), ava_plots)
    }
    
    
    
  } else {
    ava_plots <- 
      do.call(
        "grid.arrange", 
        c(delta_plots, ncol=floor(sqrt(length(delta_plots))))
      )
    ggsave(sprintf("%s/ava_plots.png", results_dir), ava_plots)
  }
  rm(ava_plots)

  mean_mps <- 
    ggplot(
      all_midpoints
    )+
    stat_density(
      aes(x=mid_dist, y=..scaled..), 
      geom = "line", position = "identity", trim = F, adjust=0.1
    ) +
    annotate(geom="text", x=0, y=0, label="|",
             color="black") +
    annotate(geom="text", x=0, y=-0.09, label="dnaA",
             color="black", size = 5) +
    #geom_bar(aes(x=midpoint, y=..scaled..), stat="identity", width = 0.001) +
    # geom_text(aes(x=0, y=-0, label = "|")) +
    # geom_text(aes(x=0, y=-.09, label = "dnaA", size = 18)) +
    xlim(-0.5,0.5) +
    theme(
      legend.position = "none" ,
      plot.title = element_text(hjust = 0.5, size = 22),
      axis.text = element_blank() ,
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      plot.margin=unit(c(0,0,0,0),"mm")
    ) +
    ggtitle(sprintf("")) +
    scale_y_continuous(limits = c(-1, 1)) +
    geom_hline(yintercept = 0) +
    coord_polar(start = pi) 
  ggsave(sprintf("%s/mean_mps.png", results_dir), mean_mps)
    
  
  # ggplot(
  #     all_midpoints
  #   ) +
  #   stat_density(
  #     aes(x=mid_dist, y=..scaled..), 
  #     geom = "line", position = "identity", trim = F, adjust=0.1
  #   ) +
  #   geom_bar(aes(x=midpoint, y=mpc_prop), stat="identity", width = 0.001) +
  #   geom_text(aes(x=0, y=-0, label = "|")) +
  #   geom_text(aes(x=0, y=-.09, label = "OriC")) +
  #   xlim(-.5,.5) +
  #   theme(
  #     legend.position = "none" , 
  #     plot.title = element_text(hjust = 0.5, size = 14),
  #     axis.text = element_blank() ,
  #     axis.ticks = element_blank(),
  #     axis.title = element_blank()
  #   ) +
  #   ggtitle(
  #     sprintf(
  #       "%s \n Inversion Midpoint Proportion", 
  #       species_name
  #       )
  #   ) +
  #   scale_y_continuous(limits = c(-1, 1)) +
  #   geom_hline(yintercept = 0) +
  #   coord_polar() 
  
  identity_sld_hm <- 
    ggplot(
      data = clust_identity, 
      aes(
        x= factor(ref, unique(ref)),
        y= factor(qry, rev(unique(ref))),
        fill = identity_sld
      )
    ) + 
    geom_tile() + 
    theme_classic() +
    theme(
      axis.text.x = element_text(angle=45, hjust = 1),
      panel.background=element_rect(fill="white"),
      plot.title = element_text(size=20,face="bold", hjust = 0.5),
      axis.title = element_blank()
    ) + 
    ggtitle(sprintf("%s length adjusted identity", species_name)) + 
    theme(legend.title=element_text(face="bold", size=2)) + 
    scale_fill_viridis(option = "B") 
  
  if (length(accession_list)>50){
    identity_sld_hm <- 
      identity_sld_hm +
      theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()
      )
  }
  if (cluster_count<=20){
    identity_sld_hm <- 
      identity_sld_hm +
      facet_grid(
        rows = vars(qry_clust), cols = vars(ref_clust), 
        scales = "free", space = "free", switch = "y")
    
  }
  
  ggsave(sprintf("%s/identity_sld_hm.png", results_dir), identity_sld_hm)
  rm(identity_sld_hm)
  
  summary_tib <- 
    tibble(
      genoforms = 
        length(unique(pull(filter(cluster_tib, clust !="LOC"), clust))),
      clusters = 
        length(unique(pull(cluster_tib, clust))),
      mean_midpoint = 
        all_midpoints %>% pull(midpoint) %>% mean(.) %>% round(., 2),
      mean_kmer_prop = mean(pull(possible_kmers, kmer_prop)),
      kmer_prop_SD = sd(pull(possible_kmers, kmer_prop))
    )
  write_csv(summary_tib, sprintf("%s/out_summary.csv", results_dir))
  # print(
  #   sprintf(
  #     "%s done, genome %i", 
  #     species_name, 
  #     which(species_list$species2==species_name)
  #     )
  #   )
  print(
    sprintf(
      "%s done, genome %i of %i", 
      species_name, 
      which(spec2_list==species_name),
      length(spec2_list)
    )
  )
  spec2_list
  print("10 sec nap X)")
  Sys.sleep(10)
  print("Ok :)")
  #}, error=function(e){} )
}
  ## finish
  