## for 2 genes
expr_plot_2G <- function(t1, t2, filename, which_cell, ylimit1, ylimit2)
{
  local_cell <- subset(local, subset = cell_type == which_cell)
  gene_var <- c(grep(t1, rownames(local_cell), value = TRUE),
                grep(t2, rownames(local_cell), value = TRUE))
  
  p <- list()
  for (i in gene_var) 
  {
    gene1<- FetchData(local_cell, vars = i)
    colnames(gene1) <- "gene"
    df <- data.frame(local_cell$disease_id, local_cell$disease_state, gene1$gene)
    colnames(df) <- c("Samples", "Type", "Gene")
    
    df<-setDT(df)[ , .(Gene = mean(Gene)), by = Samples]
    
    df <- add_column(df, Type = stri_replace_all_regex(str = df$Samples,
                                                       pattern = c('^C.*','^A.*', '^T.*'), 
                                                       replacement = disease_order,
                                                       vectorize_all = FALSE), .after = "Samples") 
    
    p[[i]] <- ggviolin(df, x = "Type", y = "Gene", fill = "Type",
                       palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                       add=c("boxplot", "jitter"),add.params = list(fill="white"),
                       order = disease_order,
                       shape = "Type", #size = 0.1,
                       ylab = i, xlab = FALSE) #, xlab = "Groups")
  }
  p1 <- ggpar(p[[1]], ylim = ylimit1)
  p2 <- ggpar(p[[2]], ylim = ylimit2)
  
  # arrange the three plots in a single row
  prow <- plot_grid( p1 + theme(legend.position="none"),
                     p2 + theme(legend.position="none"),
                     align = 'vh',
                     hjust = -1,
                     nrow = 1
  )
  legend_b <- get_legend(p1 + theme(legend.position="top"))
  
  # add the legend underneath the row we made earlier. Give it 10% of the height
  # of one plot (via rel_heights).
  p <- plot_grid( prow, legend_b, ncol = 1, rel_heights = c(1, .2))
  return(p)
}

## for 3 genes
expr_plot_3G <- function(t1, t2, t3, filename, which_cell, ylimit1, ylimit2, ylimit3)
{
  local_cell <- subset(local, subset = cell_type == which_cell)
  gene_var <- c(grep(t1, rownames(local_cell), value = TRUE),
                grep(t2, rownames(local_cell), value = TRUE),
                grep(t3, rownames(local_cell), value = TRUE))
  
  p <- list()
  for (i in gene_var) 
  {
    gene1<- FetchData(local_cell, vars = i)
    colnames(gene1) <- "gene"
    df <- data.frame(local_cell$disease_id, local_cell$disease_state, gene1$gene)
    colnames(df) <- c("Samples", "Type", "Gene")
    # table(df$local_beta.sample_id)
    
    df<-setDT(df)[ , .(Gene = mean(Gene)), by = Samples]
    
    df <- add_column(df, Type = stri_replace_all_regex(str = df$Samples,
                                                       pattern = c('^C.*','^A.*', '^T.*'), 
                                                       replacement = disease_order,
                                                       vectorize_all = FALSE), .after = "Samples") 
    
    p[[i]] <- ggviolin(df, x = "Type", y = "Gene", fill = "Type",
                       palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                       add=c("boxplot", "jitter"),add.params = list(fill="white"),
                       order = disease_order,
                       shape = "Type", #size = 0.1,
                       ylab = i, xlab = FALSE) #, xlab = "Groups")
  }
  p1 <- ggpar(p[[1]], ylim = ylimit1)
  p2 <- ggpar(p[[2]], ylim = ylimit2)
  p3 <- ggpar(p[[3]], ylim = ylimit3)
  
  # arrange the three plots in a single row
  prow <- plot_grid( p1 + theme(legend.position="none"),
                     p2 + theme(legend.position="none"),
                     p3 + theme(legend.position="none"),
                     align = 'vh',
                     # labels = c("A", "B", "C"),
                     hjust = -1,
                     nrow = 1
  )
  legend_b <- get_legend(p1 + theme(legend.position="top"))
  
  # add the legend underneath the row we made earlier. Give it 10% of the height
  # of one plot (via rel_heights).
  p <- plot_grid( prow, legend_b, ncol = 1, rel_heights = c(1, .2))
  return(p)
}

## for 4 genes
expr_plot_4G <- function(t1, t2, t3, t4, filename, which_cell, ylimit1, ylimit2, ylimit3, ylimit4)
{
  local_cell <- subset(local, subset = cell_type == which_cell)
  gene_var <- c(grep(t1, rownames(local_cell), value = TRUE),
                grep(t2, rownames(local_cell), value = TRUE),
                grep(t3, rownames(local_cell), value = TRUE),
                grep(t4, rownames(local_cell), value = TRUE))
  
  p <- list()
  for (i in gene_var) 
  {
    gene1<- FetchData(local_cell, vars = i)
    colnames(gene1) <- "gene"
    df <- data.frame(local_cell$disease_id, local_cell$disease_state, gene1$gene)
    colnames(df) <- c("Samples", "Type", "Gene")
    # table(df$local_beta.sample_id)
    
    df<-setDT(df)[ , .(Gene = mean(Gene)), by = Samples]
    
    df <- add_column(df, Type = stri_replace_all_regex(str = df$Samples,
                                                       pattern = c('^C.*','^A.*', '^T.*'), 
                                                       replacement = disease_order,
                                                       vectorize_all = FALSE), .after = "Samples") 
    
    p[[i]] <- ggviolin(df, x = "Type", y = "Gene", fill = "Type",
                       palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                       add=c("boxplot", "jitter"),add.params = list(fill="white"),
                       order = disease_order,
                       shape = "Type", #size = 0.1,
                       ylab = i, xlab = FALSE) #, xlab = "Groups")
  }
  p1 <- ggpar(p[[1]], ylim = ylimit1)
  p2 <- ggpar(p[[2]], ylim = ylimit2)
  p3 <- ggpar(p[[3]], ylim = ylimit3)
  p4 <- ggpar(p[[4]], ylim = ylimit4)
  # arrange the three plots in a single row
  prow <- plot_grid( p1 + theme(legend.position="none"),
                     p2 + theme(legend.position="none"),
                     p3 + theme(legend.position="none"),
                     p4 + theme(legend.position="none"),
                     align = 'vh',
                     # labels = c("A", "B", "C"),
                     hjust = -1,
                     nrow = 1
  )
  legend_b <- get_legend(p1 + theme(legend.position="top"))
  
  # add the legend underneath the row we made earlier. Give it 10% of the height
  # of one plot (via rel_heights).
  p <- plot_grid( prow, legend_b, ncol = 1, rel_heights = c(1, .2))
  return(p)
}

## for 5 genes
expr_plot_5G <- function(t1, t2, t3, t4, t5, filename, which_cell, ylimit1, ylimit2, ylimit3, ylimit4, ylimit5)
{
  local_cell <- subset(local, subset = cell_type == which_cell)
  gene_var <- c(grep(t1, rownames(local_cell), value = TRUE),
                grep(t2, rownames(local_cell), value = TRUE),
                grep(t3, rownames(local_cell), value = TRUE),
                grep(t4, rownames(local_cell), value = TRUE),
                grep(t5, rownames(local_cell), value = TRUE))
  
  p <- list()
  for (i in gene_var)
  {
    gene1<- FetchData(local_cell, vars = i)
    colnames(gene1) <- "gene"
    df <- data.frame(local_cell$disease_id, local_cell$disease_state, gene1$gene)
    colnames(df) <- c("Samples", "Type", "Gene")
    # table(df$local_beta.sample_id)
    
    df<-setDT(df)[ , .(Gene = mean(Gene)), by = Samples]
    
    df <- add_column(df, Type = stri_replace_all_regex(str = df$Samples,
                                                       pattern = c('^C.*','^A.*', '^T.*'),
                                                       replacement = disease_order,
                                                       vectorize_all = FALSE), .after = "Samples")
    
    p[[i]] <- ggviolin(df, x = "Type", y = "Gene", fill = "Type",
                       palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                       add=c("boxplot", "jitter"),add.params = list(fill="white"),
                       order = disease_order,
                       shape = "Type", #size = 0.1,
                       ylab = i, xlab = FALSE) #, xlab = "Groups")
  }
  p1 <- ggpar(p[[1]], ylim = ylimit1)
  p2 <- ggpar(p[[2]], ylim = ylimit2)
  p3 <- ggpar(p[[3]], ylim = ylimit3)
  p4 <- ggpar(p[[4]], ylim = ylimit4)
  p5 <- ggpar(p[[5]], ylim = ylimit5)
  
  # arrange the three plots in a single row
  prow <- plot_grid( p1 + theme(legend.position="none"),
                     p2 + theme(legend.position="none"),
                     p3 + theme(legend.position="none"),
                     p4 + theme(legend.position="none"),
                     p5 + theme(legend.position="none"),
                     align = 'vh',
                     # labels = c("A", "B", "C"),
                     hjust = -1,
                     nrow = 1
  )
  legend_b <- get_legend(p1 + theme(legend.position="top"))
  
  # add the legend underneath the row we made earlier. Give it 10% of the height
  # of one plot (via rel_heights).
  p <- plot_grid( prow, legend_b, ncol = 1, rel_heights = c(1, .2))
  return(p)
}

#################################
## for 5 genes MC
expr_plot_5genes_MC <- function(t1, t2, t3, t4, t5, filename, local_cond, local_cond_filtered, disease_order_MC,
                                which_cell, ylimit1, ylimit2, ylimit3, ylimit4, ylimit5)
{
  local_cell <- subset(local_cond, subset = cell_type == which_cell)
  gene_var <- c(grep(t1, rownames(local_cell), value = TRUE),
                grep(t2, rownames(local_cell), value = TRUE),
                grep(t3, rownames(local_cell), value = TRUE),
                grep(t4, rownames(local_cell), value = TRUE),
                grep(t5, rownames(local_cell), value = TRUE))
  local_cell_filtered <- subset(local_cond_filtered, subset = cell_type == which_cell)
  
  p <- list()
  for (i in gene_var)
  {
    gene1<- FetchData(local_cell, vars = i)
    colnames(gene1) <- "gene"
    df1 <- data.frame(local_cell$disease_id, local_cell$disease_state, gene1$gene)
    colnames(df1) <- c("Samples", "Type", "Gene")
    
    gene2<- FetchData(local_cell_filtered, vars = i)
    colnames(gene2) <- "gene"
    df2 <- data.frame(local_cell_filtered$disease_id, local_cell_filtered$disease_state, gene2$gene)
    colnames(df2) <- c("Samples", "Type", "Gene")
    df2$Samples <- str_replace(df2$Samples, "AAB_", "MC-AAB_")
    df2$Type <- str_replace(df2$Type, "AAB", "MC-AAB")
    df <- rbind(df1,df2)
    df<-setDT(df)[ , .(Gene = mean(Gene)), by = Samples]
    
    df <- add_column(df, Type = stri_replace_all_regex(str = df$Samples,
                                                       pattern = c('^C.*','^M.*', '^T.*'), 
                                                       replacement = disease_order_MC,
                                                       vectorize_all = FALSE), .after = "Samples") 
    
    p[[i]] <- ggviolin(df, x = "Type", y = "Gene", fill = "Type",
                       palette = c("#00AFBB", "#E7B800", "#FC4E07"), # "#FC4E07",
                       add=c("boxplot", "jitter"),add.params = list(fill="white"),
                       order = disease_order_MC,
                       shape = "Type", #size = 0.1,
                       ylab = i, xlab = FALSE) #+ scale_x_discrete(labels=c("Control" = "Control", "AAB" = "MC-AAB",
    #"T1D" = "T1D")) + theme(legend.position = "none")
  }
  p1 <- ggpar(p[[1]], ylim = ylimit1)
  p2 <- ggpar(p[[2]], ylim = ylimit2)
  p3 <- ggpar(p[[3]], ylim = ylimit3)
  p4 <- ggpar(p[[4]], ylim = ylimit4)
  p5 <- ggpar(p[[5]], ylim = ylimit5)
  
  # arrange the three plots in a single row
  prow <- plot_grid( p1 + theme(legend.position="none"),
                     p2 + theme(legend.position="none"),
                     p3 + theme(legend.position="none"),
                     p4 + theme(legend.position="none"),
                     p5 + theme(legend.position="none"),
                     align = 'vh',
                     hjust = -1,
                     nrow = 1
  )
  legend_b <- get_legend(p1 + theme(legend.position="top"))
  
  # add the legend underneath the row we made earlier. Give it 10% of the height
  # of one plot (via rel_heights).
  p <- plot_grid(prow, ncol = 1, rel_heights = c(1, .2), legend_b) 
  return(p)
}

## for 4 genes MC
expr_plot_4genes_MC <- function(t1, t2, t3, t4, filename, local_cond, local_cond_filtered, disease_order_MC,
                                which_cell, ylimit1, ylimit2, ylimit3, ylimit4)
{
  local_cell <- subset(local_cond, subset = cell_type == which_cell)
  gene_var <- c(grep(t1, rownames(local_cell), value = TRUE),
                grep(t2, rownames(local_cell), value = TRUE),
                grep(t3, rownames(local_cell), value = TRUE),
                grep(t4, rownames(local_cell), value = TRUE))
  local_cell_filtered <- subset(local_cond_filtered, subset = cell_type == which_cell)
  
  p <- list()
  for (i in gene_var)
  {
    gene1<- FetchData(local_cell, vars = i)
    colnames(gene1) <- "gene"
    df1 <- data.frame(local_cell$disease_id, local_cell$disease_state, gene1$gene)
    colnames(df1) <- c("Samples", "Type", "Gene")
    
    gene2<- FetchData(local_cell_filtered, vars = i)
    colnames(gene2) <- "gene"
    df2 <- data.frame(local_cell_filtered$disease_id, local_cell_filtered$disease_state, gene2$gene)
    colnames(df2) <- c("Samples", "Type", "Gene")
    df2$Samples <- str_replace(df2$Samples, "AAB_", "MC-AAB_")
    df2$Type <- str_replace(df2$Type, "AAB", "MC-AAB")
    df <- rbind(df1,df2)
    df<-setDT(df)[ , .(Gene = mean(Gene)), by = Samples]
    
    df <- add_column(df, Type = stri_replace_all_regex(str = df$Samples,
                                                       pattern = c('^C.*','^M.*', '^T.*'), 
                                                       replacement = disease_order_MC,
                                                       vectorize_all = FALSE), .after = "Samples") 
    
    p[[i]] <- ggviolin(df, x = "Type", y = "Gene", fill = "Type",
                       palette = c("#00AFBB", "#E7B800", "#FC4E07"), # "#FC4E07",
                       add=c("boxplot", "jitter"),add.params = list(fill="white"),
                       order = disease_order_MC,
                       shape = "Type", #size = 0.1,
                       ylab = i, xlab = FALSE) #+ scale_x_discrete(labels=c("Control" = "Control", "AAB" = "MC-AAB",
    #"T1D" = "T1D")) + theme(legend.position = "none")
  }
  p1 <- ggpar(p[[1]], ylim = ylimit1)
  p2 <- ggpar(p[[2]], ylim = ylimit2)
  p3 <- ggpar(p[[3]], ylim = ylimit3)
  p4 <- ggpar(p[[4]], ylim = ylimit4)
  
  # arrange the three plots in a single row
  prow <- plot_grid( p1 + theme(legend.position="none"),
                     p2 + theme(legend.position="none"),
                     p3 + theme(legend.position="none"),
                     p4 + theme(legend.position="none"),
                     align = 'vh',
                     hjust = -1,
                     nrow = 1
  )
  legend_b <- get_legend(p1 + theme(legend.position="top"))
  
  # add the legend underneath the row we made earlier. Give it 10% of the height
  # of one plot (via rel_heights).
  p <- plot_grid(prow, ncol = 1, rel_heights = c(1, .2), legend_b) 
  return(p)
}


## for 3 genes MC
expr_plot_3genes_MC <- function(t1, t2, t3, filename, local_cond, local_cond_filtered, disease_order_MC,
                                which_cell, ylimit1, ylimit2, ylimit3)
{
  local_cell <- subset(local_cond, subset = cell_type == which_cell)
  gene_var <- c(grep(t1, rownames(local_cell), value = TRUE),
                grep(t2, rownames(local_cell), value = TRUE),
                grep(t3, rownames(local_cell), value = TRUE))
  local_cell_filtered <- subset(local_cond_filtered, subset = cell_type == which_cell)
  
  p <- list()
  for (i in gene_var)
  {
    gene1<- FetchData(local_cell, vars = i)
    colnames(gene1) <- "gene"
    df1 <- data.frame(local_cell$disease_id, local_cell$disease_state, gene1$gene)
    colnames(df1) <- c("Samples", "Type", "Gene")
    
    gene2<- FetchData(local_cell_filtered, vars = i)
    colnames(gene2) <- "gene"
    df2 <- data.frame(local_cell_filtered$disease_id, local_cell_filtered$disease_state, gene2$gene)
    colnames(df2) <- c("Samples", "Type", "Gene")
    df2$Samples <- str_replace(df2$Samples, "AAB_", "MC-AAB_")
    df2$Type <- str_replace(df2$Type, "AAB", "MC-AAB")
    df <- rbind(df1,df2)
    df<-setDT(df)[ , .(Gene = mean(Gene)), by = Samples]
    
    df <- add_column(df, Type = stri_replace_all_regex(str = df$Samples,
                                                       pattern = c('^C.*','^M.*', '^T.*'), 
                                                       replacement = disease_order_MC,
                                                       vectorize_all = FALSE), .after = "Samples") 
    
    p[[i]] <- ggviolin(df, x = "Type", y = "Gene", fill = "Type",
                       palette = c("#00AFBB", "#E7B800", "#FC4E07"), # "#FC4E07",
                       add=c("boxplot", "jitter"),add.params = list(fill="white"),
                       order = disease_order_MC,
                       shape = "Type", #size = 0.1,
                       ylab = i, xlab = FALSE) #+ scale_x_discrete(labels=c("Control" = "Control", "AAB" = "MC-AAB",
    #"T1D" = "T1D")) + theme(legend.position = "none")
  }
  p1 <- ggpar(p[[1]], ylim = ylimit1)
  p2 <- ggpar(p[[2]], ylim = ylimit2)
  p3 <- ggpar(p[[3]], ylim = ylimit3)
  
  # arrange the three plots in a single row
  prow <- plot_grid( p1 + theme(legend.position="none"),
                     p2 + theme(legend.position="none"),
                     p3 + theme(legend.position="none"),
                     align = 'vh',
                     hjust = -1,
                     nrow = 1
  )
  legend_b <- get_legend(p1 + theme(legend.position="top"))
  
  # add the legend underneath the row we made earlier. Give it 10% of the height
  # of one plot (via rel_heights).
  p <- plot_grid(prow, ncol = 1, rel_heights = c(1, .2), legend_b) 
  return(p)
}

## for 2 genes MC
expr_plot_2genes_MC <- function(t1, t2, filename, local_cond, local_cond_filtered, disease_order_MC,
                                which_cell, ylimit1, ylimit2)
{
  local_cell <- subset(local_cond, subset = cell_type == which_cell)
  gene_var <- c(grep(t1, rownames(local_cell), value = TRUE),
                grep(t2, rownames(local_cell), value = TRUE))
  local_cell_filtered <- subset(local_cond_filtered, subset = cell_type == which_cell)
  
  p <- list()
  for (i in gene_var)
  {
    gene1<- FetchData(local_cell, vars = i)
    colnames(gene1) <- "gene"
    df1 <- data.frame(local_cell$disease_id, local_cell$disease_state, gene1$gene)
    colnames(df1) <- c("Samples", "Type", "Gene")
    
    gene2<- FetchData(local_cell_filtered, vars = i)
    colnames(gene2) <- "gene"
    df2 <- data.frame(local_cell_filtered$disease_id, local_cell_filtered$disease_state, gene2$gene)
    colnames(df2) <- c("Samples", "Type", "Gene")
    df2$Samples <- str_replace(df2$Samples, "AAB_", "MC-AAB_")
    df2$Type <- str_replace(df2$Type, "AAB", "MC-AAB")
    df <- rbind(df1,df2)
    df<-setDT(df)[ , .(Gene = mean(Gene)), by = Samples]
    
    df <- add_column(df, Type = stri_replace_all_regex(str = df$Samples,
                                                       pattern = c('^C.*','^M.*', '^T.*'), 
                                                       replacement = disease_order_MC,
                                                       vectorize_all = FALSE), .after = "Samples") 
    
    p[[i]] <- ggviolin(df, x = "Type", y = "Gene", fill = "Type",
                       palette = c("#00AFBB", "#E7B800", "#FC4E07"), # "#FC4E07",
                       add=c("boxplot", "jitter"),add.params = list(fill="white"),
                       order = disease_order_MC,
                       shape = "Type", #size = 0.1,
                       ylab = i, xlab = FALSE) #+ scale_x_discrete(labels=c("Control" = "Control", "AAB" = "MC-AAB",
    #"T1D" = "T1D")) + theme(legend.position = "none")
  }
  p1 <- ggpar(p[[1]], ylim = ylimit1)
  p2 <- ggpar(p[[2]], ylim = ylimit2)
  
  # arrange the three plots in a single row
  prow <- plot_grid( p1 + theme(legend.position="none"),
                     p2 + theme(legend.position="none"),
                     align = 'vh',
                     hjust = -1,
                     nrow = 1
  )
  legend_b <- get_legend(p1 + theme(legend.position="top"))
  
  # add the legend underneath the row we made earlier. Give it 10% of the height
  # of one plot (via rel_heights).
  p <- plot_grid(prow, ncol = 1, rel_heights = c(1, .2), legend_b) 
  return(p)
}

###########################################################################################################
# ## for 5 genes
expr_plot_5genes_all <- function(t1, t2, t3, t4, t5, filename, local, disease_order,
                                 which_cell, ylimit1, ylimit2, ylimit3, ylimit4, ylimit5)
{
  local_cell <- subset(local, subset = cell_type == which_cell)
  gene_var <- c(grep(t1, rownames(local_cell), value = TRUE),
                grep(t2, rownames(local_cell), value = TRUE),
                grep(t3, rownames(local_cell), value = TRUE),
                grep(t4, rownames(local_cell), value = TRUE),
                grep(t5, rownames(local_cell), value = TRUE))

  p <- list()
  for (i in gene_var)
  {
    gene1<- FetchData(local_cell, vars = i)
    colnames(gene1) <- "gene"
    df <- data.frame(local_cell$disease_id, local_cell$disease_state, gene1$gene)
    colnames(df) <- c("Samples", "Type", "Gene")

    df<-setDT(df)[ , .(Gene = mean(Gene)), by = Samples]

    df <- add_column(df, Type = stri_replace_all_regex(str = df$Samples,
                                                       pattern = c('^C.*','^A.*', '^T.*'),
                                                       replacement = disease_order,
                                                       vectorize_all = FALSE), .after = "Samples")

    p[[i]] <- ggviolin(df, x = "Type", y = "Gene", fill = "Type",
                       palette = c("#00AFBB", "#E7B800", "#FC4E07"), #
                       add=c("boxplot", "jitter"),add.params = list(fill="white"),
                       order = disease_order,
                       shape = "Type", #size = 0.1,
                       ylab = i, xlab = FALSE) #, xlab = "Groups")
  }
  p1 <- ggpar(p[[1]], ylim = ylimit1)
  p2 <- ggpar(p[[2]], ylim = ylimit2)
  p3 <- ggpar(p[[3]], ylim = ylimit3)
  p4 <- ggpar(p[[4]], ylim = ylimit4)
  p5 <- ggpar(p[[5]], ylim = ylimit5)

  # arrange the three plots in a single row
  prow <- plot_grid( p1 + theme(legend.position="none"),
                     p2 + theme(legend.position="none"),
                     p3 + theme(legend.position="none"),
                     p4 + theme(legend.position="none"),
                     p5 + theme(legend.position="none"),
                     align = 'vh',
                     hjust = -1,
                     nrow = 1
  )
  legend_b <- get_legend(p1 + theme(legend.position="top"))

  # add the legend underneath the row we made earlier. Give it 10% of the height
  # of one plot (via rel_heights).
  p <- plot_grid( prow, legend_b, ncol = 1, rel_heights = c(1, .2))
  return(p)
}

## for 4 genes
expr_plot_4genes_all <- function(t1, t2, t3, t4, filename, local, disease_order,
                                 which_cell, ylimit1, ylimit2, ylimit3, ylimit4)
{
  local_cell <- subset(local, subset = cell_type == which_cell)
  gene_var <- c(grep(t1, rownames(local_cell), value = TRUE),
                grep(t2, rownames(local_cell), value = TRUE),
                grep(t3, rownames(local_cell), value = TRUE),
                grep(t4, rownames(local_cell), value = TRUE))

  p <- list()
  for (i in gene_var)
  {
    gene1<- FetchData(local_cell, vars = i)
    colnames(gene1) <- "gene"
    df <- data.frame(local_cell$disease_id, local_cell$disease_state, gene1$gene)
    colnames(df) <- c("Samples", "Type", "Gene")

    df<-setDT(df)[ , .(Gene = mean(Gene)), by = Samples]

    df <- add_column(df, Type = stri_replace_all_regex(str = df$Samples,
                                                       pattern = c('^C.*','^A.*', '^T.*'),
                                                       replacement = disease_order,
                                                       vectorize_all = FALSE), .after = "Samples")

    p[[i]] <- ggviolin(df, x = "Type", y = "Gene", fill = "Type",
                       palette = c("#00AFBB", "#E7B800", "#FC4E07"), #
                       add=c("boxplot", "jitter"),add.params = list(fill="white"),
                       order = disease_order,
                       shape = "Type", #size = 0.1,
                       ylab = i, xlab = FALSE) #, xlab = "Groups")
  }
  p1 <- ggpar(p[[1]], ylim = ylimit1)
  p2 <- ggpar(p[[2]], ylim = ylimit2)
  p3 <- ggpar(p[[3]], ylim = ylimit3)
  p4 <- ggpar(p[[4]], ylim = ylimit4)

  # arrange the three plots in a single row
  prow <- plot_grid( p1 + theme(legend.position="none"),
                     p2 + theme(legend.position="none"),
                     p3 + theme(legend.position="none"),
                     p4 + theme(legend.position="none"),
                     align = 'vh',
                     hjust = -1,
                     nrow = 1
  )
  legend_b <- get_legend(p1 + theme(legend.position="top"))

  # add the legend underneath the row we made earlier. Give it 10% of the height
  # of one plot (via rel_heights).
  p <- plot_grid( prow, legend_b, ncol = 1, rel_heights = c(1, .2))
  return(p)
}

## for top 3 genes
expr_plot_3genes_all <- function(t1, t2, t3, filename, local, disease_order,
                                 which_cell, ylimit1, ylimit2, ylimit3)
{
  local_cell <- subset(local, subset = cell_type == which_cell)
  gene_var <- c(grep(t1, rownames(local_cell), value = TRUE),
                grep(t2, rownames(local_cell), value = TRUE),
                grep(t3, rownames(local_cell), value = TRUE))

  p <- list()
  for (i in gene_var)
  {
    gene1<- FetchData(local_cell, vars = i)
    colnames(gene1) <- "gene"
    df <- data.frame(local_cell$disease_id, local_cell$disease_state, gene1$gene)
    colnames(df) <- c("Samples", "Type", "Gene")

    df<-setDT(df)[ , .(Gene = mean(Gene)), by = Samples]

    df <- add_column(df, Type = stri_replace_all_regex(str = df$Samples,
                                                       pattern = c('^C.*','^A.*', '^T.*'),
                                                       replacement = disease_order,
                                                       vectorize_all = FALSE), .after = "Samples")

    p[[i]] <- ggviolin(df, x = "Type", y = "Gene", fill = "Type",
                       palette = c("#00AFBB", "#E7B800", "#FC4E07"), #
                       add=c("boxplot", "jitter"),add.params = list(fill="white"),
                       order = disease_order,
                       shape = "Type", #size = 0.1,
                       ylab = i, xlab = FALSE) #, xlab = "Groups")
  }
  p1 <- ggpar(p[[1]], ylim = ylimit1)
  p2 <- ggpar(p[[2]], ylim = ylimit2)
  p3 <- ggpar(p[[3]], ylim = ylimit3)

  # arrange the three plots in a single row
  prow <- plot_grid( p1 + theme(legend.position="none"),
                     p2 + theme(legend.position="none"),
                     p3 + theme(legend.position="none"),
                     align = 'vh',
                     hjust = -1,
                     nrow = 1
  )
  legend_b <- get_legend(p1 + theme(legend.position="top"))

  # add the legend underneath the row we made earlier. Give it 10% of the height
  # of one plot (via rel_heights).
  p <- plot_grid( prow, legend_b, ncol = 1, rel_heights = c(1, .2))
  return(p)
}


## for top 2 genes
expr_plot_2genes_all <- function(t1, t2, filename, local, disease_order,
                                 which_cell, ylimit1, ylimit2)
{
  local_cell <- subset(local, subset = cell_type == which_cell)
  gene_var <- c(grep(t1, rownames(local_cell), value = TRUE),
                grep(t2, rownames(local_cell), value = TRUE))
  
  p <- list()
  for (i in gene_var)
  {
    gene1<- FetchData(local_cell, vars = i)
    colnames(gene1) <- "gene"
    df <- data.frame(local_cell$disease_id, local_cell$disease_state, gene1$gene)
    colnames(df) <- c("Samples", "Type", "Gene")
    
    df<-setDT(df)[ , .(Gene = mean(Gene)), by = Samples]
    
    df <- add_column(df, Type = stri_replace_all_regex(str = df$Samples,
                                                       pattern = c('^C.*','^A.*', '^T.*'),
                                                       replacement = disease_order,
                                                       vectorize_all = FALSE), .after = "Samples")
    
    p[[i]] <- ggviolin(df, x = "Type", y = "Gene", fill = "Type",
                       palette = c("#00AFBB", "#E7B800", "#FC4E07"), #
                       add=c("boxplot", "jitter"),add.params = list(fill="white"),
                       order = disease_order,
                       shape = "Type", #size = 0.1,
                       ylab = i, xlab = FALSE) #, xlab = "Groups")
  }
  p1 <- ggpar(p[[1]], ylim = ylimit1)
  p2 <- ggpar(p[[2]], ylim = ylimit2)
  
  # arrange the three plots in a single row
  prow <- plot_grid( p1 + theme(legend.position="none"),
                     p2 + theme(legend.position="none"),
                     align = 'vh',
                     hjust = -1,
                     nrow = 1
  )
  legend_b <- get_legend(p1 + theme(legend.position="top"))
  
  # add the legend underneath the row we made earlier. Give it 10% of the height
  # of one plot (via rel_heights).
  p <- plot_grid( prow, legend_b, ncol = 1, rel_heights = c(1, .2))
  return(p)
}
## for top 1 genes
expr_plot_1gene_all <- function(t1, filename, local, disease_order,
                                which_cell, ylimit1)
{
  local_cell <- subset(local, subset = cell_type == which_cell)
  gene_var <- c(grep(t1, rownames(local_cell), value = TRUE))

  p <- list()
  for (i in gene_var)
  {
    gene1<- FetchData(local_cell, vars = i)
    colnames(gene1) <- "gene"
    df <- data.frame(local_cell$disease_id, local_cell$disease_state, gene1$gene)
    colnames(df) <- c("Samples", "Type", "Gene")

    df<-setDT(df)[ , .(Gene = mean(Gene)), by = Samples]

    df <- add_column(df, Type = stri_replace_all_regex(str = df$Samples,
                                                       pattern = c('^C.*','^A.*', '^T.*'),
                                                       replacement = disease_order,
                                                       vectorize_all = FALSE), .after = "Samples")

    p[[i]] <- ggviolin(df, x = "Type", y = "Gene", fill = "Type",
                       palette = c("#00AFBB", "#E7B800", "#FC4E07"), #
                       add=c("boxplot", "jitter"),add.params = list(fill="white"),
                       order = disease_order,
                       shape = "Type", #size = 0.1,
                       ylab = i, xlab = FALSE) #, xlab = "Groups")
  }
  p1 <- ggpar(p[[1]], ylim = ylimit1)

  # arrange the three plots in a single row
  prow <- plot_grid( p1 + theme(legend.position="none"),
                     align = 'vh',
                     hjust = -1,
                     nrow = 1
  )
  legend_b <- get_legend(p1 + theme(legend.position="top"))

  # add the legend underneath the row we made earlier. Give it 10% of the height
  # of one plot (via rel_heights).
  p <- plot_grid( prow, legend_b, ncol = 1, rel_heights = c(1, .2))
  return(p)
}
################################################################################
## for 6 genes combined
expr_plot_6genes_combined <- function(t1, t2, t3, t4, t5, t6, filename, local, local_cond_filtered, disease_order_combined,
                                      which_cell, ylimit1, ylimit2, ylimit3, ylimit4, ylimit5, ylimit6)
{
  local_cell <- subset(local, subset = cell_type == which_cell)
  gene_var <- c(grep(t1, rownames(local_cell), value = TRUE),
                grep(t2, rownames(local_cell), value = TRUE),
                grep(t3, rownames(local_cell), value = TRUE),
                grep(t4, rownames(local_cell), value = TRUE),
                grep(t5, rownames(local_cell), value = TRUE),
                grep(t6, rownames(local_cell), value = TRUE))
  local_cell_filtered <- subset(local_cond_filtered, subset = cell_type == which_cell)
  
  p <- list()
  for (i in gene_var)
  {
    gene1<- FetchData(local_cell, vars = i)
    colnames(gene1) <- "gene"
    df1 <- data.frame(local_cell$disease_id, local_cell$disease_state, gene1$gene)
    colnames(df1) <- c("Samples", "Type", "Gene")
    
    gene2<- FetchData(local_cell_filtered, vars = i)
    colnames(gene2) <- "gene"
    df2 <- data.frame(local_cell_filtered$disease_id, local_cell_filtered$disease_state, gene2$gene)
    colnames(df2) <- c("Samples", "Type", "Gene")
    df2$Samples <- str_replace(df2$Samples, "AAB_", "MC-AAB_")
    df2$Type <- str_replace(df2$Type, "AAB", "MC-AAB")
    df <- rbind(df1,df2)
    df<-setDT(df)[ , .(Gene = mean(Gene)), by = Samples]
    
    df <- add_column(df, Type = stri_replace_all_regex(str = df$Samples,
                                                       pattern = c('^C.*', '^A.*', '^M.*', '^T.*'), 
                                                       replacement = disease_order_combined,
                                                       vectorize_all = FALSE), .after = "Samples") 
    
    p[[i]] <- ggviolin(df, x = "Type", y = "Gene", fill = "Type",
                       palette = c('#21918c', '#440154', '#FC4E07', '#fde725'), # "#FC4E07",
                       #palette = c("#00AFBB", "#E7B800", "#3CB371", "#FC4E07"), # "#FC4E07",
                       add=c("boxplot", "jitter"),add.params = list(fill="white"),
                       order = disease_order_combined,
                       shape = "Type", #size = 0.1,
                       ylab = i, xlab = FALSE) #+ scale_x_discrete(labels=c("Control" = "Control", "AAB" = "MC-AAB",
    #"T1D" = "T1D")) + theme(legend.position = "none")
  }
  p1 <- ggpar(p[[1]], ylim = ylimit1)
  p2 <- ggpar(p[[2]], ylim = ylimit2)
  p3 <- ggpar(p[[3]], ylim = ylimit3)
  p4 <- ggpar(p[[4]], ylim = ylimit4)
  p5 <- ggpar(p[[5]], ylim = ylimit5)
  p6 <- ggpar(p[[6]], ylim = ylimit6)
  
  # arrange the three plots in a single row
  prow <- plot_grid( p1 + theme(legend.position="none"),
                     p2 + theme(legend.position="none"),
                     p3 + theme(legend.position="none"),
                     p4 + theme(legend.position="none"),
                     p5 + theme(legend.position="none"),
                     p6 + theme(legend.position="none"),
                     align = 'vh',
                     hjust = -1,
                     nrow = 1
  )
  legend_b <- get_legend(p1 + theme(legend.position="top"))
  
  # add the legend underneath the row we made earlier. Give it 10% of the height
  # of one plot (via rel_heights).
  p <- plot_grid(prow, ncol = 1, rel_heights = c(1, .2), legend_b) 
  return(p)
}

## for 5 genes combined
expr_plot_5genes_combined <- function(t1, t2, t3, t4, t5, filename, local, local_cond_filtered, disease_order_combined,
                                      which_cell, ylimit1, ylimit2, ylimit3, ylimit4, ylimit5)
{
  local_cell <- subset(local, subset = cell_type == which_cell)
  gene_var <- c(grep(t1, rownames(local_cell), value = TRUE),
                grep(t2, rownames(local_cell), value = TRUE),
                grep(t3, rownames(local_cell), value = TRUE),
                grep(t4, rownames(local_cell), value = TRUE),
                grep(t5, rownames(local_cell), value = TRUE))
  local_cell_filtered <- subset(local_cond_filtered, subset = cell_type == which_cell)
  
  p <- list()
  for (i in gene_var)
  {
    gene1<- FetchData(local_cell, vars = i)
    colnames(gene1) <- "gene"
    df1 <- data.frame(local_cell$disease_id, local_cell$disease_state, gene1$gene)
    colnames(df1) <- c("Samples", "Type", "Gene")
    
    gene2<- FetchData(local_cell_filtered, vars = i)
    colnames(gene2) <- "gene"
    df2 <- data.frame(local_cell_filtered$disease_id, local_cell_filtered$disease_state, gene2$gene)
    colnames(df2) <- c("Samples", "Type", "Gene")
    df2$Samples <- str_replace(df2$Samples, "AAB_", "MC-AAB_")
    df2$Type <- str_replace(df2$Type, "AAB", "MC-AAB")
    df <- rbind(df1,df2)
    df<-setDT(df)[ , .(Gene = mean(Gene)), by = Samples]
    
    df <- add_column(df, Type = stri_replace_all_regex(str = df$Samples,
                                                       pattern = c('^C.*', '^A.*', '^M.*', '^T.*'), 
                                                       replacement = disease_order_combined,
                                                       vectorize_all = FALSE), .after = "Samples") 
    
    p[[i]] <- ggviolin(df, x = "Type", y = "Gene", fill = "Type",
                       palette = c('#21918c', '#440154', '#FC4E07', '#fde725'), # "#FC4E07",
		       #palette = c("#00AFBB", "#E7B800", "#3CB371", "#FC4E07"), # "#FC4E07",
                       add=c("boxplot", "jitter"),add.params = list(fill="white"),
                       order = disease_order_combined,
                       shape = "Type", #size = 0.1,
                       ylab = i, xlab = FALSE) #+ scale_x_discrete(labels=c("Control" = "Control", "AAB" = "MC-AAB",
    #"T1D" = "T1D")) + theme(legend.position = "none")
  }
  p1 <- ggpar(p[[1]], ylim = ylimit1)
  p2 <- ggpar(p[[2]], ylim = ylimit2)
  p3 <- ggpar(p[[3]], ylim = ylimit3)
  p4 <- ggpar(p[[4]], ylim = ylimit4)
  p5 <- ggpar(p[[5]], ylim = ylimit5)
  
  # arrange the three plots in a single row
  prow <- plot_grid( p1 + theme(legend.position="none"),
                     p2 + theme(legend.position="none"),
                     p3 + theme(legend.position="none"),
                     p4 + theme(legend.position="none"),
                     p5 + theme(legend.position="none"),
                     align = 'vh',
                     hjust = -1,
                     nrow = 1
  )
  legend_b <- get_legend(p1 + theme(legend.position="top"))
  
  # add the legend underneath the row we made earlier. Give it 10% of the height
  # of one plot (via rel_heights).
  p <- plot_grid(prow, ncol = 1, rel_heights = c(1, .2), legend_b) 
  return(p)
}


## for 4 genes combined
expr_plot_4genes_combined <- function(t1, t2, t3, t4, filename, local, local_cond_filtered, disease_order_combined,
                                which_cell, ylimit1, ylimit2, ylimit3, ylimit4)
{
  local_cell <- subset(local, subset = cell_type == which_cell)
  gene_var <- c(grep(t1, rownames(local_cell), value = TRUE),
                grep(t2, rownames(local_cell), value = TRUE),
                grep(t3, rownames(local_cell), value = TRUE),
                grep(t4, rownames(local_cell), value = TRUE))
  local_cell_filtered <- subset(local_cond_filtered, subset = cell_type == which_cell)
  
  p <- list()
  for (i in gene_var)
  {
    gene1<- FetchData(local_cell, vars = i)
    colnames(gene1) <- "gene"
    df1 <- data.frame(local_cell$disease_id, local_cell$disease_state, gene1$gene)
    colnames(df1) <- c("Samples", "Type", "Gene")
    
    gene2<- FetchData(local_cell_filtered, vars = i)
    colnames(gene2) <- "gene"
    df2 <- data.frame(local_cell_filtered$disease_id, local_cell_filtered$disease_state, gene2$gene)
    colnames(df2) <- c("Samples", "Type", "Gene")
    df2$Samples <- str_replace(df2$Samples, "AAB_", "MC-AAB_")
    df2$Type <- str_replace(df2$Type, "AAB", "MC-AAB")
    df <- rbind(df1,df2)
    df<-setDT(df)[ , .(Gene = mean(Gene)), by = Samples]
    
    df <- add_column(df, Type = stri_replace_all_regex(str = df$Samples,
                                                       pattern = c('^C.*', '^A.*', '^M.*', '^T.*'), 
                                                       replacement = disease_order_combined,
                                                       vectorize_all = FALSE), .after = "Samples") 
    
    p[[i]] <- ggviolin(df, x = "Type", y = "Gene", fill = "Type",
                      palette = c('#21918c', '#440154', '#FC4E07', '#fde725'), # "#FC4E07",                       
		      #palette = c("#00AFBB", "#E7B800", "#3CB371", "#FC4E07"), # "#FC4E07",
                       add=c("boxplot", "jitter"),add.params = list(fill="white"),
                       order = disease_order_combined,
                       shape = "Type", #size = 0.1,
                       ylab = i, xlab = FALSE) #+ scale_x_discrete(labels=c("Control" = "Control", "AAB" = "MC-AAB",
    #"T1D" = "T1D")) + theme(legend.position = "none")
  }
  p1 <- ggpar(p[[1]], ylim = ylimit1)
  p2 <- ggpar(p[[2]], ylim = ylimit2)
  p3 <- ggpar(p[[3]], ylim = ylimit3)
  p4 <- ggpar(p[[4]], ylim = ylimit4)
  
  # arrange the three plots in a single row
  prow <- plot_grid( p1 + theme(legend.position="none"),
                     p2 + theme(legend.position="none"),
                     p3 + theme(legend.position="none"),
                     p4 + theme(legend.position="none"),
                     align = 'vh',
                     hjust = -1,
                     nrow = 1
  )
  legend_b <- get_legend(p1 + theme(legend.position="top"))
  
  # add the legend underneath the row we made earlier. Give it 10% of the height
  # of one plot (via rel_heights).
  p <- plot_grid(prow, ncol = 1, rel_heights = c(1, .2), legend_b) 
  return(p)
}

## for 4 genes combined
expr_plot_4genes_combined_wolimits <- function(t1, t2, t3, t4, filename, local, local_cond_filtered, disease_order_combined,
                                      which_cell)
{
  local_cell <- subset(local, subset = cell_type == which_cell)
  gene_var <- c(grep(t1, rownames(local_cell), value = TRUE),
                grep(t2, rownames(local_cell), value = TRUE),
                grep(t3, rownames(local_cell), value = TRUE),
                grep(t4, rownames(local_cell), value = TRUE))
  local_cell_filtered <- subset(local_cond_filtered, subset = cell_type == which_cell)
  
  p <- list()
  for (i in gene_var)
  {
    gene1<- FetchData(local_cell, vars = i)
    colnames(gene1) <- "gene"
    df1 <- data.frame(local_cell$disease_id, local_cell$disease_state, gene1$gene)
    colnames(df1) <- c("Samples", "Type", "Gene")
    
    gene2<- FetchData(local_cell_filtered, vars = i)
    colnames(gene2) <- "gene"
    df2 <- data.frame(local_cell_filtered$disease_id, local_cell_filtered$disease_state, gene2$gene)
    colnames(df2) <- c("Samples", "Type", "Gene")
    df2$Samples <- str_replace(df2$Samples, "AAB_", "MC-AAB_")
    df2$Type <- str_replace(df2$Type, "AAB", "MC-AAB")
    df <- rbind(df1,df2)
    df<-setDT(df)[ , .(Gene = mean(Gene)), by = Samples]
    
    df <- add_column(df, Type = stri_replace_all_regex(str = df$Samples,
                                                       pattern = c('^C.*', '^A.*', '^M.*', '^T.*'), 
                                                       replacement = disease_order_combined,
                                                       vectorize_all = FALSE), .after = "Samples") 
    
    p[[i]] <- ggviolin(df, x = "Type", y = "Gene", fill = "Type",
                      palette = c('#21918c', '#440154', '#FC4E07', '#fde725'), # "#FC4E07",                       
#palette = c("#00AFBB", "#E7B800", "#3CB371", "#FC4E07"), # "#FC4E07",
                       add=c("boxplot", "jitter"),add.params = list(fill="white"),
                       order = disease_order_combined,
                       shape = "Type", #size = 0.1,
                       ylab = i, xlab = FALSE) #+ scale_x_discrete(labels=c("Control" = "Control", "AAB" = "MC-AAB",
    #"T1D" = "T1D")) + theme(legend.position = "none")
  }
  # p1 <- ggpar(p[[1]], ylim = ylimit1)
  # p2 <- ggpar(p[[2]], ylim = ylimit2)
  # p3 <- ggpar(p[[3]], ylim = ylimit3)
  # p4 <- ggpar(p[[4]], ylim = ylimit4)
  
  p1 <- ggpar(p[[1]])
  p2 <- ggpar(p[[2]])
  p3 <- ggpar(p[[3]])
  p4 <- ggpar(p[[4]])
  
  # arrange the three plots in a single row
  prow <- plot_grid( p1 + theme(legend.position="none"),
                     p2 + theme(legend.position="none"),
                     p3 + theme(legend.position="none"),
                     p4 + theme(legend.position="none"),
                     align = 'vh',
                     hjust = -1,
                     nrow = 1
  )
  legend_b <- get_legend(p1 + theme(legend.position="top"))
  
  # add the legend underneath the row we made earlier. Give it 10% of the height
  # of one plot (via rel_heights).
  p <- plot_grid(prow, ncol = 1, rel_heights = c(1, .2), legend_b) 
  return(p)
}

## for 3 genes combined
expr_plot_3genes_combined <- function(t1, t2, t3, filename, local, local_cond_filtered, disease_order_combined,
                                      which_cell, ylimit1, ylimit2, ylimit3)
{
  local_cell <- subset(local, subset = cell_type == which_cell)
  gene_var <- c(grep(t1, rownames(local_cell), value = TRUE),
                grep(t2, rownames(local_cell), value = TRUE),
                grep(t3, rownames(local_cell), value = TRUE))
  local_cell_filtered <- subset(local_cond_filtered, subset = cell_type == which_cell)
  
  p <- list()
  for (i in gene_var)
  {
    gene1<- FetchData(local_cell, vars = i)
    colnames(gene1) <- "gene"
    df1 <- data.frame(local_cell$disease_id, local_cell$disease_state, gene1$gene)
    colnames(df1) <- c("Samples", "Type", "Gene")
    
    gene2<- FetchData(local_cell_filtered, vars = i)
    colnames(gene2) <- "gene"
    df2 <- data.frame(local_cell_filtered$disease_id, local_cell_filtered$disease_state, gene2$gene)
    colnames(df2) <- c("Samples", "Type", "Gene")
    df2$Samples <- str_replace(df2$Samples, "AAB_", "MC-AAB_")
    df2$Type <- str_replace(df2$Type, "AAB", "MC-AAB")
    df <- rbind(df1,df2)
    df<-setDT(df)[ , .(Gene = mean(Gene)), by = Samples]
    
    df <- add_column(df, Type = stri_replace_all_regex(str = df$Samples,
                                                       pattern = c('^C.*', '^A.*', '^M.*', '^T.*'), 
                                                       replacement = disease_order_combined,
                                                       vectorize_all = FALSE), .after = "Samples") 
    
    p[[i]] <- ggviolin(df, x = "Type", y = "Gene", fill = "Type",
                      palette = c('#21918c', '#440154', '#FC4E07', '#fde725'), # "#FC4E07", 
#                      palette = c("#00AFBB", "#E7B800", "#3CB371", "#FC4E07"), # "#FC4E07",
                       # palette = c("#00AFBB", "#E7B800", "#FC4E07", "#3CB371"), # "#FC4E07",
                       
                       add=c("boxplot", "jitter"),add.params = list(fill="white"),
                       order = disease_order_combined,
                       shape = "Type", #size = 0.1,
                       ylab = i, xlab = FALSE) #+ scale_x_discrete(labels=c("Control" = "Control", "AAB" = "MC-AAB",
    #"T1D" = "T1D")) + theme(legend.position = "none")
  }
  p1 <- ggpar(p[[1]], ylim = ylimit1)
  p2 <- ggpar(p[[2]], ylim = ylimit2)
  p3 <- ggpar(p[[3]], ylim = ylimit3)
  
  # arrange the three plots in a single row
  prow <- plot_grid( p1 + theme(legend.position="none"),
                     p2 + theme(legend.position="none"),
                     p3 + theme(legend.position="none"),
                     align = 'vh',
                     hjust = -1,
                     nrow = 1
  )
  legend_b <- get_legend(p1 + theme(legend.position="top"))
  
  # add the legend underneath the row we made earlier. Give it 10% of the height
  # of one plot (via rel_heights).
  p <- plot_grid(prow, ncol = 1, rel_heights = c(1, .2), legend_b) 
  return(p)
}


## for 2 genes combined
expr_plot_2genes_combined <- function(t1, t2, filename, local, local_cond_filtered, disease_order_combined,
                                      which_cell, ylimit1, ylimit2)
{
  local_cell <- subset(local, subset = cell_type == which_cell)
  gene_var <- c(grep(t1, rownames(local_cell), value = TRUE),
                grep(t2, rownames(local_cell), value = TRUE))
  local_cell_filtered <- subset(local_cond_filtered, subset = cell_type == which_cell)
  
  p <- list()
  for (i in gene_var)
  {
    gene1<- FetchData(local_cell, vars = i)
    colnames(gene1) <- "gene"
    df1 <- data.frame(local_cell$disease_id, local_cell$disease_state, gene1$gene)
    colnames(df1) <- c("Samples", "Type", "Gene")
    
    gene2<- FetchData(local_cell_filtered, vars = i)
    colnames(gene2) <- "gene"
    df2 <- data.frame(local_cell_filtered$disease_id, local_cell_filtered$disease_state, gene2$gene)
    colnames(df2) <- c("Samples", "Type", "Gene")
    df2$Samples <- str_replace(df2$Samples, "AAB_", "MC-AAB_")
    df2$Type <- str_replace(df2$Type, "AAB", "MC-AAB")
    df <- rbind(df1,df2)
    df<-setDT(df)[ , .(Gene = mean(Gene)), by = Samples]
    
    df <- add_column(df, Type = stri_replace_all_regex(str = df$Samples,
                                                       pattern = c('^C.*', '^A.*', '^M.*', '^T.*'), 
                                                       replacement = disease_order_combined,
                                                       vectorize_all = FALSE), .after = "Samples") 
    
    p[[i]] <- ggviolin(df, x = "Type", y = "Gene", fill = "Type",
                      palette = c('#21918c', '#440154', '#FC4E07', '#fde725'), # "#FC4E07", 
#                      palette = c("#00AFBB", "#E7B800", "#3CB371", "#FC4E07"), # "#FC4E07",
                       add=c("boxplot", "jitter"),add.params = list(fill="white"),
                       order = disease_order_combined,
                       shape = "Type", #size = 0.1,
                       ylab = i, xlab = FALSE) #+ scale_x_discrete(labels=c("Control" = "Control", "AAB" = "MC-AAB",
    #"T1D" = "T1D")) + theme(legend.position = "none")
  }
  p1 <- ggpar(p[[1]], ylim = ylimit1)
  p2 <- ggpar(p[[2]], ylim = ylimit2)
  
  # arrange the three plots in a single row
  prow <- plot_grid( p1 + theme(legend.position="none"),
                     p2 + theme(legend.position="none"),
                     align = 'vh',
                     hjust = -1,
                     nrow = 1
  )
  legend_b <- get_legend(p1 + theme(legend.position="top"))
  
  # add the legend underneath the row we made earlier. Give it 10% of the height
  # of one plot (via rel_heights).
  p <- plot_grid(prow, ncol = 1, rel_heights = c(1, .2), legend_b) 
  return(p)
}

