require("tidyverse")
library("openxlsx")
require("xtable")
require("RColorBrewer")
require("grid")
require("gtable")
require("viridis")
require("ggpubr")
require("paletteer")


#' f
#' A function mapping a combination of maintained genotypes to a color
#' @param x The genotypes maintained
#' @param contab A translation table mapping a combination of maintained genotypes to a specific color
#'
#' @return A vector of colors. white is NA value
#' @export
#'
#' @examples
f <- function(x,contab){
  
  xtab <- tibble(value=x)
  
  outtab <- as.data.frame(left_join(xtab,contab,by="value"))
  outtab <- outtab %>% 
    mutate(labels=replace_na(labels,"white"))
  
  return(outtab[,"labels"])
}

#' f_poly
#'
#'Functions mapping a polymorphism pattern to a color
#'
#' @param x The maintained genotypes/input to the function is a vector of length nb values of cP(1)
#' @param conpoly A translation table mapping polymorphism patterns to a color
#' @param colored Should the color really be outputted. If FALSE just the patterns are generated.
#'
#' @return
#' @export ret A cellcolor string for a latex table
#'
#' @examples
f_poly <- function(x,conpoly,colored=TRUE){
  # Define the possible patterns
  pp <- c("P","M")
  # Split the string to get single haplotypes
  comb_haplo <- strsplit(x,", ")
  
  # For each cH1/cP1 combination 
  out <- lapply(comb_haplo,function(comb_haplotypes){
    # Get each single locus for all maintained haplotypes
    # result will be a data.frame with m-haplotype columns, n loci rows 
    dat <- data.frame(sapply(comb_haplotypes, strsplit,""))
    # Convert to numeric
    dat <- apply(dat,2,as.numeric)
    # Calculate the rowSums
    pat <- rowSums(dat)
    
    # Check if the row sum is 0 (all haplotypes have 0 allele) or ncol(dat) (all haplotypes have effective allele)
    patcheck <- as.numeric(pat==ncol(dat) | pat==0)
    # Construct the polymorphism pattern
    patret <- paste(paste0(pp[patcheck+1],1:length(patcheck)),collapse=", ")
    # Calculate the color for the polymorphism pattern
    # Bit logic
    patval <- sum(patcheck * 2^(0:(nrow(dat)-1)) )

    if(colored){
      ret <- paste0('\\cellcolor{', conpoly[patval+1], "}", patret)
    }else{
      ret <- patret
    }
    return(ret)
    
  })
  
return(unlist(out))
 
}

#' rename_columns
#'
#' Function to rename column such that they have the correct subscripts and superscripts
#' The function currently is very specific and converts CH1 to c_H^{(1)} and CH2 to c_H^{(2)}
#'
#' @param cnames Columnames
#'
#' @return
#' @export
#'
#' @examples
rename_columns <- function(cnames){
  # Prepare latex output
  newcnames <- paste0("$",cnames,"$")
  newcnames <- gsub("CH1","c_H^{(1)}",newcnames)
  newcnames <- gsub("CP1","c_P^{(1)}",newcnames)
  return(newcnames)
}

#' generate_combres
#' 
#' Generates a Pivot-table from the summary table
#'
#' @param sum_dat Summary table having summaries for various combination of parameter values
#' @param species The species for which to create the summary table. Can be H1, H2 or P
#' @param nfrom The column from which to construct the names
#' @param vfrom The column from which to collect the values
#' @param ffrom The column which should be kept fixed
#' @param npref Prefix to add to the names columns
#'
#' @return sum_over A tibble in wide format
#' @export
#'
#' @examples
generate_combres <- function(sum_dat,species,nfrom,vfrom,ffrom,npref){
  sum_over <- sum_dat %>%  
    filter(Species==species) %>% 
    select(!!nfrom,!!ffrom,!!vfrom) %>% 
    pivot_wider(names_from=!!nfrom,values_from=!!vfrom,names_prefix = npref) %>% 
    arrange(!!ffrom)
  return(sum_over)
}


#' generate_nbgenos_tab
#'
#' @param sum_table  A summary table for different combinations of cH_1 and cH_2.
#' @param latex Prepare output for latex table?
#'
#' @return A colored summary latex table. Cells are colored according to the number of genotypes maintained
#' @export
#'
#' @examples
generate_nbgenos_tab <- function(sum_table,latex=TRUE){
  sum_table_dat <- as.data.frame(sum_table)
  sum_table_dat[names(sum_table_dat)[-1]] <- lapply(sum_table_dat[names(sum_table_dat)[-1]], function(x){
    sapply(strsplit(x,","),length)})
  if(latex){colnames(sum_table_dat) <- rename_columns(names(sum_table))}
  return(sum_table_dat)
}

#' generate_poly_tab
#'
#' @param sum_table  A summary table for different combinations of cH_1 and cH_2 in wide format
#' @param conpoly The translation table
#' @colored TRUE return table with latex colors?
#' @return
#' @export
#'
#' @examples
generate_poly_tab <- function(sum_table,conpoly,colored=TRUE){
  sum_table_dat <- as.data.frame(sum_table)
  sum_table_dat[names(sum_table_dat)[-1]] <- lapply(sum_table_dat[names(sum_table_dat)[-1]], function(x){
    f_poly(x,conpoly,colored)})
  if(colored){
    colnames(sum_table_dat) <- rename_columns(names(sum_table))
  }
  return(sum_table_dat)
}


combmaintaineddat <- function(maintained_unparsed){
  maintained <- strsplit(maintained_unparsed, ", ")
  maintained_tabs <- sapply(maintained,function(x){
    # Split string into single loci
    split_haplos <- strsplit(x,"")
  })
  # Combine the list of maintained haplotypes for different cost combinations into a single data frame
  # Each column gives the allelic state of a single maintained haplotype
  # Each row is a single locus. Hence the number of of rows should correspond to the number of loci in the simulation.
  maintained_dat <- do.call(cbind.data.frame,maintained_tabs)
  # Convert all columns from string to numeric
  maintained_dat <- data.frame(apply(maintained_dat,2,as.numeric,drop=F))
  return(maintained_dat)
}

rangemaintained <- function(sum_dat,haplocol){
  haplos <- sum_dat %>% select(!!haplocol) %>% pull()
  a <- lapply(haplos, function(x){
    maintained_dat <- combmaintaineddat(x)
    found_ranges <- colSums(maintained_dat)
    
    possible_ranges <- 0:3
    
    range_boolean <- possible_ranges%in%found_ranges
    range_boolean
    
    range_weight <- 2^(possible_ranges)
    
    range_summary <- sum(range_boolean*range_weight)
    
    return(range_summary)
  })
  return(unlist(a))
}

rgenes_maintained <- function(sum_dat_host, haplocol){
  haplos <- sum_dat_host |> 
    select(!!haplocol) |> pull()
  a <- lapply(haplos, function(x){
    maintained_dat <- combmaintaineddat(x)
    found_ranges <- colSums(maintained_dat)
    
    possible_ranges <- 0:3
    
    range_boolean <- possible_ranges%in%found_ranges
    range_boolean
    
    range_weight <- 2^(possible_ranges)
    
    range_summary <- sum(range_boolean*range_weight)
    
    R_mat <- as.matrix(t(combmaintaineddat(x)))
    summary_Ralleles_present <- (colSums(R_mat) > 0)
    Poly_vec <- apply(R_mat,2, function(x){length(unique(x))-1})
    
    ret_dat <- tibble(range_summary = range_summary, R_1 = summary_Ralleles_present[1], 
                      R_2 = summary_Ralleles_present[2],
                      R_3 = summary_Ralleles_present[3],
                      poly_1 = Poly_vec[1],
                      poly_2 = Poly_vec[2],
                      poly_3 = Poly_vec[3])
    return(ret_dat)
  })
  daten <- bind_rows(a)
  combined <- bind_cols(sum_dat_host, daten) |>
    mutate(R_alleles = case_when(
      Species == "H1" ~ R_2 + 2*R_3,
      Species == "H2" ~ R_1 + 2*R_3,
      TRUE ~ -999
      )
    ) |>
    mutate(poly = case_when(
      Species == "H1" ~ poly_2 + 2*poly_3,
      Species == "H2" ~ poly_1 + 2*poly_3,
      TRUE ~ -999
    ))
}


plot_maintained <- function(sum_dat, haplocol, seed_select)
{
  sum_dat_host <- sum_dat |> 
    filter(Species%in%c("H1", "H2"))
  combined <- rgenes_maintained(sum_dat_host, haplocol)
  combined <- combined |>
    mutate(poly = factor(as.character(poly), levels = c("0","1","2", "3"))) |>
    mutate(R_alleles = factor(as.character(R_alleles), levels = c("0","1","2", "3"))) |>
    mutate(range_summary = factor(as.character(range_summary), levels = c("1","2", "3","4", "5", "6", "7")))
  
  pmain_plt <- ggplot(combined |> filter(seed == seed_select), aes(CP1, CH1)) +
    geom_tile(aes(fill=factor(poly)), color = "gray80", show.legend = TRUE) + 
    facet_grid(cols = vars(phi), rows = vars(Species), labeller=label_bquote(cols = phi == .(phi))) + 
    scale_fill_manual(values = c("0" = "#B2DFDBFF", 
                                 "1" = "#00A087FF",
                                 "2" = "#64B5F6",
                                 "3" = "#3C5488FF"), 
                      labels=c("0"="none", "1"= "private only","2"="ancestral only","3"="both"), 
                      name = "", 
                      drop=FALSE) +
    ggtitle("Polymorphism maintained") +
    theme_bw()
  
  rgenemain_plt <- ggplot(combined, aes(CP1, CH1)) +
    geom_tile(aes(fill=factor(R_alleles)), color = "gray80", show.legend = TRUE) + 
    facet_grid(cols = vars(phi), rows = vars(Species)) + 
    facet_grid(cols = vars(phi), rows = vars(Species), labeller=label_bquote(cols = phi == .(phi))) + 
    scale_fill_manual(values = c("0" = "#B2DFDBFF", 
                                 "1" = "#00A087FF",
                                 "2" = "#64B5F6",
                                 "3" = "#3C5488FF"), labels=c("0"="none", "1"= "private only","2"="ancestral only","3"="both"), 
                      name = "", 
                      drop = FALSE) +
    ggtitle("R alleles maintained") + 
    theme_bw()
  
  range_main_plt <- ggplot(combined, aes(CP1, CH1)) +
    geom_tile(aes(fill=factor(range_summary)), color = "gray80", show.legend = TRUE) + 
    facet_grid(cols = vars(phi), rows = vars(Species)) + 
    facet_grid(cols = vars(phi), rows = vars(Species), labeller=label_bquote(cols = phi == .(phi))) + 
    scale_fill_manual(values = c("1" = "#FFE0B2FF",
                                 "2" = "#FFA726FF",
                                 "3" = "#F57C00FF",
                                 "4" = "#FF95A8FF",
                                 "5" = "#EC407AFF",
                                 "6" = "#C2185BFF",
                                 "7" = "#8A4198FF"),        
                      labels=c("1" = "0",
                               "2" = "1",
                               "3" = "0 + 1",
                               "4" = "2",
                               "5" = "0 + 2",
                               "6" = "1 + 2",
                               "7" = "0 + 1 + 2"), 
                      breaks=as.character(1:7),  name = "",
                      drop = FALSE) +
    ggtitle("Resistance allele ranges maintained") + 
    theme_bw()
  
  
  arplt <- ggarrange(pmain_plt, rgenemain_plt, range_main_plt, nrow =3 )
  return(arplt)
}


#' f_poly_long
#' 
#' 
#' This function is used to calculate the polymorphism patterns for a given species for each locus based on the maintained haplotypes
#' The observed polymorphism pattern is a binary vector of length n-loci
#' 0 indicates that the given locus is polymorphic, 1 indicated the given locus is monomorphic
#' Once the polymorphism pattern is generated, it is summarized by converting the binary number into a decimal number
#' i.e. if l=3 and all three loci are monomorphic the pattern will be 7, if l=3 and 101 the pattern will be 5.
#' 
#' @param sum_dat A summary data file in long format. Where each row is a combination of $c_H^{(1)}$,  $c_H^{(2)}$, $c_P^{(1)}$, $c_P^{(2)}$, $\phi$ and Species.
#' @param haplocol name of the column which contains the haplotypes maintained in the given species at the end of the simulation run for the given combination of parameters.
#'
#' @return
#' @export
#'
#' @examples
f_poly_long <- function(sum_dat,haplocol){
  
  # Define the possible patterns
  pp <- c("P","M")
  
  haplos <- sum_dat %>% select(!!haplocol) %>% pull()
  a <- lapply(haplos, function(x){
    
    maintained_dat <- combmaintaineddat(x)
    #maintained <- strsplit(x, ", ")
    #maintained_tabs <- sapply(maintained,function(x){
    #  # Split string into single loci
    #  split_haplos <- strsplit(x,"")
    #})
    # Combine the list of maintained haplotypes for different cost combinations into a single data frame
    # Each column gives the allelic state of a single maintained haplotype
    # Each row is a single locus. Hence the number of of rows should correspond to the number of loci in the simulation.
    #maintained_dat <- do.call(cbind.data.frame,maintained_tabs)
    # Convert all columns from string to numeric
    #maintained_dat <- apply(maintained_dat,2,as.numeric)
    # Calculate the rowSums
    pat <- rowSums(maintained_dat)
    
    # Check if the row sum is 0 (all haplotypes have 0 allele) or ncol(dat) (all haplotypes have effective allele)
    patcheck <- as.numeric(pat==ncol(maintained_dat) | pat==0)
    # Construct the polymorphism pattern
    patret <- paste(paste0(pp[patcheck+1],1:length(patcheck)),collapse=", ")
    # Calculate the color for the polymorphism pattern
    # Bit logic
    patval <- sum(patcheck * 2^(0:(nrow(maintained_dat)-1)) )
    return(paste0("pattern",patval))
  })
  return(unlist(a))
}

#' Title
#'
#' Generates a plot for different combinations of $c_H^{(1)}$,  $c_H^{(2)}$, $c_P^{(1)}$, $c_P^{(2)}$, $\phi$ for a given polymorphism pattern combination
#' @param sum_dat A summary data file in long format. Where each row is a combination of $c_H^{(1)}$,  $c_H^{(2)}$, $c_P^{(1)}$, $c_P^{(2)}$, $\phi$ and Species.
#' 
#' sum_dat needs to have the following columns: Species, CH1, CH2, CP1, CP2, shortsum (the haplotypes maintained at the end of the simulations separated by ", ") and phi
#'
#'Depends on function: f_poly_long
#'
#' @return
#' @export
#'
#' @examples
plotpolymorphism <- function(sum_dat,haplocol){
  
  sum_dat$pattern <- f_poly_long(sum_dat,haplocol)
  sum_dat$translated_pattern <- translate_pattern(sum_dat %>% pull(pattern))
  
  species_names <- c(
    'H1'="Host 1",
    'H2'="Host 2",
    'P'="Pathogen"
  )
  
  scale_fill_viridis_nodrop <- scale_fill_manual(
    values = viridis(8), breaks=c("M,M,M","P,M,M","M,P,M","P,P,M","M,M,P","P,M,P","M,P,P","P,P,P"))
  
  
  cH_2 <- sum_dat %>% pull(CH2) %>% unique()
  cP_2 <- sum_dat %>% pull(CP2) %>% unique()
  
  p <- ggplot(sum_dat, aes(CH1,CP1)) + 
    geom_tile(aes(fill=factor(translated_pattern))) + 
    guides(fill=guide_legend(title="Polymorphism pattern")) +
    #facet_grid(Species ~ phi,labeller=label_bquote(cols = phi == .(phi),rows=.(species_names))) +
    facet_grid(Species ~ phi,labeller=label_bquote(cols = phi == .(phi))) +
    #facet_grid(Species ~ phi,labeller=labeller(
    #  Species = species_names
    #)) + 
    xlab(expression({c[H]}^{(1)})) + 
    ylab(expression({c[P]}^{(1)})) + 
    #scale_fill_viridis_d() +
    scale_fill_viridis_nodrop + 
    ggtitle(bquote(list({c[H]^{(2)}}==.(cH_2),{c[P]^{(2)}}==.(cP_2))))
  
  
  # Labels 
  labelR = "Observed patterns for species:"
  labelT = "Proportion of host species 1"
  
  # Get the ggplot grob
  z <- ggplotGrob(p)
  
  # Get the positions of the strips in the gtable: t = top, l = left, ...
  posR <- subset(z$layout, grepl("strip-r", name), select = t:r)
  posT <- subset(z$layout, grepl("strip-t", name), select = t:r)
  
  # Add a new column to the right of current right strips, 
  # and a new row on top of current top strips
  width <- z$widths[max(posR$r)]    # width of current right strips
  height <- z$heights[min(posT$t)]  # height of current top strips
  
  z <- gtable_add_cols(z, width, max(posR$r))  
  z <- gtable_add_rows(z, height, min(posT$t)-1)
  
  # Construct the new strip grobs
  stripR <- gTree(name = "Strip_right", children = gList(
    rectGrob(gp = gpar(col = NA, fill = "grey85")),
    textGrob(labelR, rot = -90, gp = gpar(fontsize = 8.8, col = "grey10"))))
  
  stripT <- gTree(name = "Strip_top", children = gList(
    rectGrob(gp = gpar(col = NA, fill = "grey85")),
    textGrob(labelT, gp = gpar(fontsize = 8.8, col = "grey10"))))
  
  # Position the grobs in the gtable
  z <- gtable_add_grob(z, stripR, t = min(posR$t)+1, l = max(posR$r) + 1, b = max(posR$b)+1, name = "strip-right")
  z <- gtable_add_grob(z, stripT, t = min(posT$t), l = min(posT$l), r = max(posT$r), name = "strip-top")
  
  # Add small gaps between strips
  z <- gtable_add_cols(z, unit(1/5, "line"), max(posR$r))
  z <- gtable_add_rows(z, unit(1/5, "line"), min(posT$t))
  
  
  return(z)
  ## Draw it
  #grid.newpage()
  #grid.draw(z)
  
}


plotpolymorphism_seed <- function(sum_dat,haplocol){
  
  sum_dat$pattern <- f_poly_long(sum_dat,haplocol)
  sum_dat$translated_pattern <- translate_pattern(sum_dat %>% pull(pattern))
  
  species_names <- c(
    'H1'="Host 1",
    'H2'="Host 2",
    'P'="Pathogen"
  )
  
  scale_fill_viridis_nodrop <- scale_fill_manual(
    values = viridis(8), breaks=c("M,M,M","P,M,M","M,P,M","P,P,M","M,M,P","P,M,P","M,P,P","P,P,P"))
  
  
  cH_2 <- sum_dat %>% pull(CH2) %>% unique()
  cP_2 <- sum_dat %>% pull(CP2) %>% unique()
  
  sum_dat_across <- sum_dat |> 
    group_by(phi, Species, CH1, CP1) |>
    summarize(nb_patterns = length(unique(translated_pattern)),
              sum_pat = paste0(table(translated_pattern), collapse = "_"),
              names_freq = paste0(names(table(translated_pattern)), collapse = "_"),
              most_freq_name = names(table(translated_pattern))[which(table(translated_pattern)==max(table(translated_pattern)))[1]],
              most_freq = as.numeric(table(translated_pattern)[which(table(translated_pattern)==max(table(translated_pattern)))[1]]),
              cat_initial = case_when(
                most_freq > 44 ~ "",
                TRUE ~ "*"
              )
              )
    
  
    p <- ggplot(sum_dat_across |> filter(!(Species == "H2" & phi == 1)), aes(CH1,CP1)) + 
    geom_tile(aes(fill=factor(most_freq_name)), color = "black") + 
    guides(fill=guide_legend(title="Polymorphism pattern")) +
    #facet_grid(Species ~ phi,labeller=label_bquote(cols = phi == .(phi),rows=.(species_names))) +
    facet_grid(Species ~ phi,labeller=label_bquote(cols = phi == .(phi))) +
    #facet_grid(Species ~ phi,labeller=labeller(
    #  Species = species_names
    #)) + 
    xlab(expression({c[H]}^{(1)})) + 
    ylab(expression({c[P]}^{(1)})) + 
    geom_text(aes(label = cat_initial), color = "black", size = 5) + 
    #scale_fill_viridis_d() +
    scale_fill_viridis_nodrop + 
    theme_bw() + 
    ggtitle(bquote(list({c[H]^{(2)}}==.(cH_2),{c[P]^{(2)}}==.(cP_2))))
  
  
  # Labels 
  labelR = "Observed patterns for species:"
  labelT = "Proportion of host species 1"
  
  # Get the ggplot grob
  z <- ggplotGrob(p)
  
  # Get the positions of the strips in the gtable: t = top, l = left, ...
  posR <- subset(z$layout, grepl("strip-r", name), select = t:r)
  posT <- subset(z$layout, grepl("strip-t", name), select = t:r)
  
  # Add a new column to the right of current right strips, 
  # and a new row on top of current top strips
  width <- z$widths[max(posR$r)]    # width of current right strips
  height <- z$heights[min(posT$t)]  # height of current top strips
  
  z <- gtable_add_cols(z, width, max(posR$r))  
  z <- gtable_add_rows(z, height, min(posT$t)-1)
  
  # Construct the new strip grobs
  stripR <- gTree(name = "Strip_right", children = gList(
    rectGrob(gp = gpar(col = NA, fill = "grey85")),
    textGrob(labelR, rot = -90, gp = gpar(fontsize = 8.8, col = "grey10"))))
  
  stripT <- gTree(name = "Strip_top", children = gList(
    rectGrob(gp = gpar(col = NA, fill = "grey85")),
    textGrob(labelT, gp = gpar(fontsize = 8.8, col = "grey10"))))
  
  # Position the grobs in the gtable
  z <- gtable_add_grob(z, stripR, t = min(posR$t)+1, l = max(posR$r) + 1, b = max(posR$b)+1, name = "strip-right")
  z <- gtable_add_grob(z, stripT, t = min(posT$t), l = min(posT$l), r = max(posT$r), name = "strip-top")
  
  # Add small gaps between strips
  z <- gtable_add_cols(z, unit(1/5, "line"), max(posR$r))
  z <- gtable_add_rows(z, unit(1/5, "line"), min(posT$t))
  
  
  return(z)
  ## Draw it
  #grid.newpage()
  #grid.draw(z)
  
}


#' Title
#' Function to plot an overview of the ranges maintained in the species at the final timepoint. Requires sum_dat the summary file for various combinations of CH1, CH2, phi and species.
#' @param sum_dat 
#' @param haplocol Name of the column which has the maintained haplotypes in 100, 000, 000 format
#'
#' @return
#' @export
#'
#' @examples
plotranges <- function(sum_dat,haplocol){
  
  sum_dat$ranges <- rangemaintained(sum_dat,haplocol)
  sum_dat$translated_ranges <- translate_rangesums(sum_dat %>% pull(ranges))
  
  species_names <- c(
    'H1'="Host 1",
    'H2'="Host 2",
    'P'="Pathogen"
  )
  
  
  # Set a manual scale
  scale_fill_viridis_nodrop <- scale_fill_manual(
      values = plasma(16), breaks=sort(c("0","1","0,1","2","0,2","1,2","0,1,2","3","0,3","1,3","0,1,3","2,3","0,2,3","1,2,3","0,1,2,3")))

  
  cH_2 <- sum_dat %>% pull(CH2) %>% unique()
  cP_2 <- sum_dat %>% pull(CP2) %>% unique()
  
  p <- ggplot(sum_dat, aes(CH1,CP1)) + 
    geom_tile(aes(fill=factor(translated_ranges))) + 
    guides(fill=guide_legend(title="Ranges maintained")) +
    #facet_grid(Species ~ phi,labeller=label_bquote(cols = phi == .(phi),rows=.(species_names)))  + 
    facet_grid(Species ~ phi,labeller=label_bquote(cols = phi == .(phi)))  + 
    theme_bw() + 
    xlab(expression({c[H]}^{(1)})) + 
    ylab(expression({c[P]}^{(1)})) + 
    #scale_fill_viridis_d(option="F")
    scale_fill_viridis_nodrop +
    ggtitle(bquote(list({c[H]^{(2)}}==.(cH_2),{c[P]^{(2)}}==.(cP_2))))
  

  
  
  # Labels 
  labelR = "Observed patterns for species:"
  labelT = "Proportion of host species 1"
  
  # Get the ggplot grob
  z <- ggplotGrob(p)
  
  # Get the positions of the strips in the gtable: t = top, l = left, ...
  posR <- subset(z$layout, grepl("strip-r", name), select = t:r)
  posT <- subset(z$layout, grepl("strip-t", name), select = t:r)
  
  # Add a new column to the right of current right strips, 
  # and a new row on top of current top strips
  width <- z$widths[max(posR$r)]    # width of current right strips
  height <- z$heights[min(posT$t)]  # height of current top strips
  
  z <- gtable_add_cols(z, width, max(posR$r))  
  z <- gtable_add_rows(z, height, min(posT$t)-1)
  
  # Construct the new strip grobs
  stripR <- gTree(name = "Strip_right", children = gList(
    rectGrob(gp = gpar(col = NA, fill = "grey85")),
    textGrob(labelR, rot = -90, gp = gpar(fontsize = 8.8, col = "grey10"))))
  
  stripT <- gTree(name = "Strip_top", children = gList(
    rectGrob(gp = gpar(col = NA, fill = "grey85")),
    textGrob(labelT, gp = gpar(fontsize = 8.8, col = "grey10"))))
  
  # Position the grobs in the gtable
  z <- gtable_add_grob(z, stripR, t = min(posR$t)+1, l = max(posR$r) + 1, b = max(posR$b)+1, name = "strip-right")
  z <- gtable_add_grob(z, stripT, t = min(posT$t), l = min(posT$l), r = max(posT$r), name = "strip-top")
  
  # Add small gaps between strips
  z <- gtable_add_cols(z, unit(1/5, "line"), max(posR$r))
  z <- gtable_add_rows(z, unit(1/5, "line"), min(posT$t))
  
  
  return(z)
  ## Draw it
  #grid.newpage()
  #grid.draw(z)
  
}


plotranges_seed <- function(sum_dat,haplocol){
  
  sum_dat$ranges <- rangemaintained(sum_dat,haplocol)
  sum_dat$translated_ranges <- translate_rangesums(sum_dat %>% pull(ranges))
  
  species_names <- c(
    'H1'="Host 1",
    'H2'="Host 2",
    'P'="Pathogen"
  )
  
  
  # Set a manual scale
  scale_fill_viridis_nodrop <- scale_fill_manual(
    values = plasma(16), breaks=sort(c("0","1","0,1","2","0,2","1,2","0,1,2","3","0,3","1,3","0,1,3","2,3","0,2,3","1,2,3","0,1,2,3")))
  
  
  cH_2 <- sum_dat %>% pull(CH2) %>% unique()
  cP_2 <- sum_dat %>% pull(CP2) %>% unique()
  

  sum_dat_across <- sum_dat |> 
    group_by(phi, Species, CH1, CP1) |>
    summarize(nb_patterns = length(unique(translated_ranges)),
              sum_pat = paste0(table(translated_ranges), collapse = "_"),
              names_freq = paste0(names(table(translated_ranges)), collapse = "_"),
              most_freq_name = names(table(translated_ranges))[which(table(translated_ranges)==max(table(translated_ranges)))[1]],
              most_freq = as.numeric(table(translated_ranges)[which(table(translated_ranges)==max(table(translated_ranges)))[1]]),
              cat_initial = case_when(
                most_freq > 44 ~ "",
                TRUE ~ "*"
              )
    )
  
  p <- ggplot(sum_dat_across |> filter(!(Species == "H2" & phi == 1)), aes(CH1,CP1)) + 
    geom_tile(aes(fill=factor(most_freq_name)), color = "grey80") + 
    guides(fill=guide_legend(title="Ranges maintained")) +
    #facet_grid(Species ~ phi,labeller=label_bquote(cols = phi == .(phi),rows=.(species_names)))  + 
    facet_grid(Species ~ phi,labeller=label_bquote(cols = phi == .(phi)))  + 
    geom_text(aes(label = cat_initial), color = "white") + 
    xlab(expression({c[H]}^{(1)})) + 
    ylab(expression({c[P]}^{(1)})) + 
    #geom_text(aes(label = most_freq), color = "white") + 
    #scale_fill_viridis_d(option="F")
    scale_fill_viridis_nodrop +
    theme_bw() + 
    ggtitle(bquote(list({c[H]^{(2)}}==.(cH_2),{c[P]^{(2)}}==.(cP_2))))
  
  
  
  
  # Labels 
  labelR = "Observed patterns for species:"
  labelT = "Proportion of host species 1"
  
  # Get the ggplot grob
  z <- ggplotGrob(p)
  
  # Get the positions of the strips in the gtable: t = top, l = left, ...
  posR <- subset(z$layout, grepl("strip-r", name), select = t:r)
  posT <- subset(z$layout, grepl("strip-t", name), select = t:r)
  
  # Add a new column to the right of current right strips, 
  # and a new row on top of current top strips
  width <- z$widths[max(posR$r)]    # width of current right strips
  height <- z$heights[min(posT$t)]  # height of current top strips
  
  z <- gtable_add_cols(z, width, max(posR$r))  
  z <- gtable_add_rows(z, height, min(posT$t)-1)
  
  # Construct the new strip grobs
  stripR <- gTree(name = "Strip_right", children = gList(
    rectGrob(gp = gpar(col = NA, fill = "grey85")),
    textGrob(labelR, rot = -90, gp = gpar(fontsize = 8.8, col = "grey10"))))
  
  stripT <- gTree(name = "Strip_top", children = gList(
    rectGrob(gp = gpar(col = NA, fill = "grey85")),
    textGrob(labelT, gp = gpar(fontsize = 8.8, col = "grey10"))))
  
  # Position the grobs in the gtable
  z <- gtable_add_grob(z, stripR, t = min(posR$t)+1, l = max(posR$r) + 1, b = max(posR$b)+1, name = "strip-right")
  z <- gtable_add_grob(z, stripT, t = min(posT$t), l = min(posT$l), r = max(posT$r), name = "strip-top")
  
  # Add small gaps between strips
  z <- gtable_add_cols(z, unit(1/5, "line"), max(posR$r))
  z <- gtable_add_rows(z, unit(1/5, "line"), min(posT$t))
  
  
  return(z)
  ## Draw it
  #grid.newpage()
  #grid.draw(z)

}


#' Title
#' Function to calculate the number of maintained haplotypes
#' Requires the maintained haplotypes for each combination of CH1, CH2, phi, species.
#' @param haplotypes 
#'
#' @return
#' @export
#'
#' @examples
getnbmaintained <- function(haplotypes){
  out_list <- lapply(haplotypes,function(h){
    length(strsplit(h,",")[[1]])
  })
  return(unlist(out_list))
}


#' Title
#' Function to plot the number of maintained haplotypes at the end of the simulation for CH1, CH2, phi, species combination.
#' @param sum_dat 
#' @param haplocol 
#'
#' @return
#' @export
#'
#' @examples
plotnbmaintained <- function(sum_dat,haplocol){
  sum_dat$nbmaintained <- getnbmaintained(sum_dat %>% pull(!!haplocol))

  species_names <- c(
    'H1'="Host 1",
    'H2'="Host 2",
    'P'="Pathogen"
  )
  
  cH_2 <- sum_dat %>% pull(CH2) %>% unique()
  cP_2 <- sum_dat %>% pull(CP2) %>% unique()
  
  scale_fill_viridis_nb <- scale_fill_manual(
    values = viridis(8), breaks=1:8)
  
  p <- ggplot(sum_dat, aes(CH1,CP1)) + 
    geom_tile(aes(fill=factor(nbmaintained))) + 
    guides(fill=guide_legend(title="Number of haplo-\ntypes maintained")) +
    #facet_grid(Species ~ phi,labeller=label_bquote(cols = phi == .(phi),rows=.(species_names))) +
    facet_grid(Species ~ phi,labeller=label_bquote(cols = phi == .(phi))) +
    #facet_grid(Species ~ phi,labeller=labeller(
    #  Species = species_names
    #)) + 
    xlab(expression({c[H]}^{(1)})) + 
    ylab(expression({c[P]}^{(1)})) + 
    scale_fill_viridis_nb +
    ggtitle(bquote(list({c[H]^{(2)}}==.(cH_2),{c[P]^{(2)}}==.(cP_2))))
  
  
  # Labels 
  labelR = "Observed patterns for species:"
  labelT = "Proportion of host species 1"
  
  # Get the ggplot grob
  z <- ggplotGrob(p)
  
  # Get the positions of the strips in the gtable: t = top, l = left, ...
  posR <- subset(z$layout, grepl("strip-r", name), select = t:r)
  posT <- subset(z$layout, grepl("strip-t", name), select = t:r)
  
  # Add a new column to the right of current right strips, 
  # and a new row on top of current top strips
  width <- z$widths[max(posR$r)]    # width of current right strips
  height <- z$heights[min(posT$t)]  # height of current top strips
  
  z <- gtable_add_cols(z, width, max(posR$r))  
  z <- gtable_add_rows(z, height, min(posT$t)-1)
  
  # Construct the new strip grobs
  stripR <- gTree(name = "Strip_right", children = gList(
    rectGrob(gp = gpar(col = NA, fill = "grey85")),
    textGrob(labelR, rot = -90, gp = gpar(fontsize = 8.8, col = "grey10"))))
  
  stripT <- gTree(name = "Strip_top", children = gList(
    rectGrob(gp = gpar(col = NA, fill = "grey85")),
    textGrob(labelT, gp = gpar(fontsize = 8.8, col = "grey10"))))
  
  # Position the grobs in the gtable
  z <- gtable_add_grob(z, stripR, t = min(posR$t)+1, l = max(posR$r) + 1, b = max(posR$b)+1, name = "strip-right")
  z <- gtable_add_grob(z, stripT, t = min(posT$t), l = min(posT$l), r = max(posT$r), name = "strip-top")
  
  # Add small gaps between strips
  z <- gtable_add_cols(z, unit(1/5, "line"), max(posR$r))
  z <- gtable_add_rows(z, unit(1/5, "line"), min(posT$t))
  
  
  return(z)
  ## Draw it
  #grid.newpage()
  #grid.draw(z)
}

plotnbmaintained_seed <- function(sum_dat,haplocol){
  
  sum_dat$nbmaintained <- getnbmaintained(sum_dat %>% pull(!!haplocol))

  species_names <- c(
    'H1'="Host 1",
    'H2'="Host 2",
    'P'="Pathogen"
  )
  
  cH_2 <- sum_dat %>% pull(CH2) %>% unique()
  cP_2 <- sum_dat %>% pull(CP2) %>% unique()
  
  scale_fill_viridis_nb <- scale_fill_manual(
    values = viridis(8), breaks=1:8)
  
  
  sum_dat_across <- sum_dat |> 
    group_by(phi, Species, CH1, CP1) |>
    summarize(nb_patterns = length(unique(nbmaintained)),
              sum_pat = paste0(table(nbmaintained), collapse = "_"),
              names_freq = paste0(names(table(nbmaintained)), collapse = "_"),
              most_freq_name = names(table(nbmaintained))[which(table(nbmaintained)==max(table(nbmaintained)))[1]],
              most_freq = as.numeric(table(nbmaintained)[which(table(nbmaintained)==max(table(nbmaintained)))[1]]),
              cat_initial = case_when(
                most_freq > 44 ~ "",
                TRUE ~ "*"
              )
    )
  
  p <- ggplot(sum_dat_across |> filter(!(Species == "H2" & phi == 1)), aes(CH1,CP1)) + 
    geom_tile(aes(fill=factor(most_freq_name)), color = "gray80") + 
    guides(fill=guide_legend(title="Number of haplo-\ntypes maintained")) +
    #facet_grid(Species ~ phi,labeller=label_bquote(cols = phi == .(phi),rows=.(species_names))) +
    facet_grid(Species ~ phi,labeller=label_bquote(cols = phi == .(phi))) +
    geom_text(aes(label = cat_initial), color = "white") + 
    #facet_grid(Species ~ phi,labeller=labeller(
    #  Species = species_names
    #)) + 
    xlab(expression({c[H]}^{(1)})) + 
    ylab(expression({c[P]}^{(1)})) + 
    scale_fill_viridis_nb +
    ggtitle(bquote(list({c[H]^{(2)}}==.(cH_2),{c[P]^{(2)}}==.(cP_2)))) + 
    theme_bw()
  
  
  # Labels 
  labelR = "Observed patterns for species:"
  labelT = "Proportion of host species 1"
  
  # Get the ggplot grob
  z <- ggplotGrob(p)
  
  # Get the positions of the strips in the gtable: t = top, l = left, ...
  posR <- subset(z$layout, grepl("strip-r", name), select = t:r)
  posT <- subset(z$layout, grepl("strip-t", name), select = t:r)
  
  # Add a new column to the right of current right strips, 
  # and a new row on top of current top strips
  width <- z$widths[max(posR$r)]    # width of current right strips
  height <- z$heights[min(posT$t)]  # height of current top strips
  
  z <- gtable_add_cols(z, width, max(posR$r))  
  z <- gtable_add_rows(z, height, min(posT$t)-1)
  
  # Construct the new strip grobs
  stripR <- gTree(name = "Strip_right", children = gList(
    rectGrob(gp = gpar(col = NA, fill = "grey85")),
    textGrob(labelR, rot = -90, gp = gpar(fontsize = 8.8, col = "grey10"))))
  
  stripT <- gTree(name = "Strip_top", children = gList(
    rectGrob(gp = gpar(col = NA, fill = "grey85")),
    textGrob(labelT, gp = gpar(fontsize = 8.8, col = "grey10"))))
  
  # Position the grobs in the gtable
  z <- gtable_add_grob(z, stripR, t = min(posR$t)+1, l = max(posR$r) + 1, b = max(posR$b)+1, name = "strip-right")
  z <- gtable_add_grob(z, stripT, t = min(posT$t), l = min(posT$l), r = max(posT$r), name = "strip-top")
  
  # Add small gaps between strips
  z <- gtable_add_cols(z, unit(1/5, "line"), max(posR$r))
  z <- gtable_add_rows(z, unit(1/5, "line"), min(posT$t))
  
  
  return(z)
  ## Draw it
  #grid.newpage()
  #grid.draw(z)
}


#' Title
#' Function to extract the final frequencies of the maintained haplotypes at the end of the simulation from the summary column of the summary file and to calculate the Simpson index. Note here it is calculated as 1-D, as D is the probability that two randomly picked entities are the same.
#' @param summarycol 
#'
#' @return
#' @export
#'
#' @examples
getsimpson <- function(summarycol){
  a <- lapply(summarycol,function(test){
    #https://stackoverflow.com/questions/24173194/remove-parentheses-and-text-within-from-strings-in-r
    freqs <-sapply(strsplit(strsplit(gsub("\\s*\\([^\\)]+\\)","",test),";")[[1]],": "),function(x){
      return(as.numeric(x[2]))
    }
    )
    return(1-sum(freqs^2))
  })
 return(unlist(a))
}


#' plotsimpson
#' Function to plot the final frequencies of the maintained haplotypes at the end of the simulation from the summary column of the summary file.
#' @param sum_dat 
#' @param summarycol 
#'
#' @return
#' @export
#'
#' @examples
plotsimpson <- function(sum_dat,summarycol){
  sum_dat$simpson <- getsimpson(sum_dat %>% pull(!!summarycol))
  
  species_names <- c(
    'H1'="Host 1",
    'H2'="Host 2",
    'P'="Pathogen"
  )
  
  
  cH_2 <- sum_dat %>% pull(CH2) %>% unique()
  cP_2 <- sum_dat %>% pull(CP2) %>% unique()
  
  # Set a manual scale
  scale_fill_viridis_nodrop <- scale_fill_manual(
    values = plasma(16), breaks=sort(c("0","1","0,1","2","0,2","1,2","0,1,2","3","0,3","1,3","0,1,3","2,3","0,2,3","1,2,3","0,1,2,3")))
  
  
  p <- ggplot(sum_dat, aes(CH1,CP1)) + 
    geom_tile(aes(fill=simpson)) + 
    guides(fill=guide_legend(title="Simpson: 1-D allele frequencies")) +
    #facet_grid(Species ~ phi,labeller=label_bquote(cols = phi == .(phi),rows=.(species_names)))  + 
    facet_grid(Species ~ phi,labeller=label_bquote(cols = phi == .(phi)))  + 
    xlab(expression({c[H]}^{(1)})) + 
    ylab(expression({c[P]}^{(1)})) +
    scale_fill_continuous(limits = c(0, 1)) +
    ggtitle(bquote(list({c[H]^{(2)}}==.(cH_2),{c[P]^{(2)}}==.(cP_2))))
  
  
  
  
  # Labels 
  labelR = "Observed patterns for species:"
  labelT = "Proportion of host species 1"
  
  # Get the ggplot grob
  z <- ggplotGrob(p)
  
  # Get the positions of the strips in the gtable: t = top, l = left, ...
  posR <- subset(z$layout, grepl("strip-r", name), select = t:r)
  posT <- subset(z$layout, grepl("strip-t", name), select = t:r)
  
  # Add a new column to the right of current right strips, 
  # and a new row on top of current top strips
  width <- z$widths[max(posR$r)]    # width of current right strips
  height <- z$heights[min(posT$t)]  # height of current top strips
  
  z <- gtable_add_cols(z, width, max(posR$r))  
  z <- gtable_add_rows(z, height, min(posT$t)-1)
  
  # Construct the new strip grobs
  stripR <- gTree(name = "Strip_right", children = gList(
    rectGrob(gp = gpar(col = NA, fill = "grey85")),
    textGrob(labelR, rot = -90, gp = gpar(fontsize = 8.8, col = "grey10"))))
  
  stripT <- gTree(name = "Strip_top", children = gList(
    rectGrob(gp = gpar(col = NA, fill = "grey85")),
    textGrob(labelT, gp = gpar(fontsize = 8.8, col = "grey10"))))
  
  # Position the grobs in the gtable
  z <- gtable_add_grob(z, stripR, t = min(posR$t)+1, l = max(posR$r) + 1, b = max(posR$b)+1, name = "strip-right")
  z <- gtable_add_grob(z, stripT, t = min(posT$t), l = min(posT$l), r = max(posT$r), name = "strip-top")
  
  # Add small gaps between strips
  z <- gtable_add_cols(z, unit(1/5, "line"), max(posR$r))
  z <- gtable_add_rows(z, unit(1/5, "line"), min(posT$t))
  
  
  return(z)
  ## Draw it
  #grid.newpage()
  #grid.draw(z)
  
}

#' translate_pattern 
#'
#' @param pattern The calculated polymorphism pattern of the form pattern0, pattern1, ... as a vector
#'
#' @return Return the pattern in X,X,X notation
#' @export
#'
#' @examples
translate_pattern <- function(pattern){
  dat <- tibble(pattern=pattern)
  dat_out <- dat %>% 
    mutate(translated_pattern=case_when(
      pattern == "pattern0" ~ "P,P,P",
      pattern == "pattern1" ~ "M,P,P",
      pattern == "pattern2" ~ "P,M,P",
      pattern == "pattern3" ~ "M,M,P",
      pattern == "pattern4" ~ "P,P,M",
      pattern == "pattern5" ~ "M,P,M",
      pattern == "pattern6" ~ "P,M,M",
      pattern == "pattern7" ~ "M,M,M"
    ))
  return(dat_out %>% pull(translated_pattern))
}

#' Title
#'
#' @param pat The ranges observed add the end of the simulation (Range= Number of 1-alleles) encoded into a binary retranslated to decimal
#'
#' @return 
#' @export
#'
#' @examples
translate_rangesums <- function(pat){
  dat <- tibble(pattern=pat)
  pattern_levels <- c("0","1","0,1","2","0,2","1,2","0,1,2","3","0,3","1,3","0,1,3","2,3","0,2,3","1,2,3","0,1,2,3")
  
  out <- pattern_levels[pat]
  out <- factor(out,levels=pattern_levels)
  return(out)
}

