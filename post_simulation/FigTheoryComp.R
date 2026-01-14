library(tidyverse)
library(deSolve)
library(reshape2)
library(nleqslv)
library(gridExtra)
require(ggimage)

GetQxy <- function(host_geno, path_geno, sigma) {
  d_xy <- host_geno * (1 - path_geno)
  d_xy <- sum(d_xy)
  Q_xy <- sigma^d_xy
  return(Q_xy)
}

GetB <- function(geno, omega, xi, host_log) {
  if(xi == 0) {
    B <- 1 - omega * sum(geno) / ifelse(host_log, 2, 3)
  } else {
    B <- 1 - omega * ((1 - exp(xi * sum(geno) / ifelse(host_log, 2, 3))) / (1 - exp(xi)))
  }
  return(B)
}

GetHostFitness <- function(B, Q, beta, path_state) {
  ret <- Q %*% path_state
  ret <- - beta * ret
  ret <- B * exp(ret)
  return(ret)
}

GetPathFitness <- function(B, Q, beta, phi, host_state) {
  phi_adj_state <- c(rep(phi, times = 4), rep((1-phi), times = 4)) * host_state
  ret <- t(Q) %*% phi_adj_state
  ret <- B * exp(beta * ret)
  return(ret)
}

GetHostFreqs <- function(host_range) {
  out_freqs <- rep(0, times = 4)
  if(host_range == 0) {
    out_freqs[1] <- 1
  } else if(host_range == 1) {
    out_freqs[2] <- 1/3
    out_freqs[3] <- 2/3
  } else if(host_range == 2) {
    out_freqs[4] <- 1
  }
  return(out_freqs)
}


GetPathFreqs <- function(path_range) {
  out_freqs <- rep(0, times = 8)
  if(path_range == 0) {
    out_freqs[1] <- 1
  } else if(path_range == 1) {
    out_freqs[2] <- 1/3
    out_freqs[3] <- 1/3
    out_freqs[4] <- 1/3
  } else if(path_range == 2) {
    out_freqs[5] <- 1/3
    out_freqs[6] <- 1/3
    out_freqs[7] <- 1/3
  } else if(path_range == 3) {
    out_freqs[8] <- 1
  }
  return(out_freqs)
}


GetNonInvRanges <- function(pars) {
  
  range_fitness <- data.frame()
  
  host_ranges <- c(0, 1, 2)
  path_ranges <- c(0, 1, 2, 3)
  
  for(h_range in host_ranges) {
    for(m_range in host_ranges) {
      for(p_range in path_ranges) {
        
        h_freq <- GetHostFreqs(h_range)
        m_freq <- GetHostFreqs(m_range)
        
        p_freq <- GetPathFreqs(p_range)
        
        host_fitness <- GetHostFitness(B = pars$B[1:8],
                                       Q = pars$Q,
                                       beta = pars$beta_H,
                                       path_state = p_freq)
        
        path_fitness <- GetPathFitness(B = pars$B[9:16],
                                       Q = pars$Q,
                                       phi = pars$phi,
                                       beta = pars$beta_P,
                                       host_state = c(h_freq, m_freq))
        
        h_fit <- host_fitness[1:4]
        m_fit <- host_fitness[5:8]
        
        h_max <- which(h_fit  >= 0.999 * max(h_fit))
        h_inv <- !prod(h_max == which(h_freq > 0))
        
        m_max <- which(m_fit  >= 0.999 * max(m_fit))
        m_inv <- !prod(m_max == which(m_freq > 0))      
        
        p_max <- which(path_fitness >= 0.999 * max(path_fitness))
        p_inv <- !prod(p_max == which(p_freq > 0))
        
        if(length(which(p_freq > 0)) %% length(p_max) != 0) {
          print("GASP")
          print(paste0("Host H Range: ", h_range,
                       " Host M Range: ", m_range,
                       " Path Range: ", p_range))
          print(length(which(p_freq > 0)))
          
          print(length(p_max))
          print(path_fitness)
          print(paste0("H Inv: ", h_inv,
                       " M Inv: ", m_inv,
                       " P Inv: ", p_inv))
        } else
          
          range_fitness <- rbind(range_fitness,
                                 data.frame(Host_H_Range = h_range,
                                            Host_M_Range = m_range,
                                            Path_Range = p_range,
                                            Host_H_Invasible = h_inv,
                                            Host_M_Invasible = m_inv,
                                            Path_Invasible = p_inv))
      }
    }
  }
  
  non_inv_ranges <- range_fitness %>%
    mutate(NonInvasible = (!Host_H_Invasible) * (!Host_M_Invasible) * (!Path_Invasible)) %>%
    filter(NonInvasible == 1)
  
  if(nrow(non_inv_ranges) == 0) {
    non_inv_ranges <- data.frame(Host_H_Range = NA,
                                 Host_M_Range = NA,
                                 Path_Range = NA,
                                 Host_H_Invasible = TRUE,
                                 Host_M_Invasible = TRUE,
                                 Path_Invasible = TRUE,
                                 NonInvasible = 0)
  }
  
  return(non_inv_ranges)
}

# integrates the dynamics using deSolve
IntegrateDynamics <- function(inistate, pars, endtime, timestep, fn){
  times <- seq(0, endtime, by = timestep)
  timeseries <- as.data.frame(ode(inistate, times, fn, pars))  
  return(timeseries)
}

BuildPars <- function(input_params) {
  pars <- with(input_params, {
    Q <- matrix(0, nrow = 8, ncol = 8)
    for(i in 1:nrow(Q)) {
      for(j in 1:ncol(Q)) {
        Q[i, j] <- GetQxy(host_geno = host_geno[i,], path_geno[j,], sigma)
      }
    }
    
    B <- rep(0, times = 16)
    for(i in 1:8) {
      B[i] <- GetB(geno = host_geno[i,], omega = omega_H, xi = xi_H, host_log = T)
    }
    
    for(i in 1:8) {
      B[i + 8] <- GetB(geno = path_geno[i,], omega = omega_P, xi = xi_P, host_log = F)
    }
    pars <- list(B = B, Q = Q, beta_H = beta_H, beta_P = beta_P, phi = phi)
    return(pars)
  })
  return(pars)
}


beta_H <- beta_M <- 1
beta_P <- 1

xi_H <- xi_M <- 3
xi_P <- 3

omega_H <- 0.1
omega_M <- 0.1
omega_P <- 0.1

phi <- seq(0.5, 0.8, length.out = 19)
sigma <- 0.85

in_pars <- crossing(beta_H = beta_H,
                    beta_P = beta_P,
                    xi_H = xi_H,
                    xi_P = xi_P,
                    omega_H = omega_H,
                    omega_P = omega_P,
                    phi = phi,
                    sigma = sigma)


ini_host <- c(rep(1 / 4, 4), rep(1/4, 4))
#ini_host <- runif(8)
#ini_host[1:4] <- ini_host[1:4] / sum(ini_host[1:4])
#ini_host[5:8] <- ini_host[5:8] / sum(ini_host[5:8])
names(ini_host) <- c("000", "001", "010", "011", "000", "001", "100", "101")
host_geno <- data.frame(Loc1 = as.numeric(substr(names(ini_host), 1, 1)),
                        Loc2 = as.numeric(substr(names(ini_host), 2, 2)),
                        Loc3 = as.numeric(substr(names(ini_host), 3, 3)))

ini_path <- rep(1 / 8, 8)
names(ini_path) <- c("000", "001", "010", "100", "011", "101", "110", "111")
path_geno <- data.frame(Loc1 = as.numeric(substr(names(ini_path), 1, 1)),
                        Loc2 = as.numeric(substr(names(ini_path), 2, 2)),
                        Loc3 = as.numeric(substr(names(ini_path), 3, 3)))
ini_state <- c(ini_host, ini_path)

end_time <- 1e5
time_step <- 1

time_cutoff <- end_time * 0.9

start_time <- Sys.time()

range_results <- data.frame()
comp_results <- data.frame()
for(i in 1:nrow(in_pars)) {
  cur_pars <- in_pars[i,]
  pars <- BuildPars(cur_pars)
  
  out <- IntegrateDynamics(ini_state, pars, end_time, time_step, fn = ReplDyn)
  colnames(out) <- paste0(c("",
                            rep("H", times = 4),
                            rep("M", times = 4),
                            rep("P", times = 8)),
                          colnames(out))
  
  final_state <- out %>%
    filter(time >= time_cutoff)
  final_state <- colMeans(final_state[,-1])
  final_state <- cbind(cur_pars, as.data.frame(t(final_state)))
  comp_results <- rbind(comp_results, final_state)
  
  #non_inv_ranges <- GetNonInvRanges(pars)
  #cur_ranges <- cbind(cur_pars, non_inv_ranges)
  #range_results <- rbind(range_results, cur_ranges)
}

print(Sys.time() - start_time)

plot_results <- comp_results %>%
  select(phi, starts_with("H"), starts_with("M"), starts_with("P")) %>%
  melt(id.vars = c("phi")) %>%
  mutate(value = ifelse(value < 1e-6, 0, value)) %>%
  mutate(value = ifelse(value > 0, value, NA)) %>%
  mutate(Genotype = substr(variable, 2, 4),
         Type = substr(variable, 1, 1)) %>%
  mutate(Type = case_when(Type == "H" ~ "Host H",
                          Type == "M" ~ "Host M",
                          Type == "P" ~ "Pathogen"))

hostH_results <- plot_results %>%
  filter(Type == "Host H" & (phi < 2/3 | Genotype == "000"))

nonunique_results <- plot_results %>%
  filter(Type == "Host H" & phi >= 2/3 & Genotype != "000")
nonunique_results$image <- "h1.png"

hostM_results <- plot_results %>%
  filter(Type == "Host M")

path_results <- plot_results %>%
  filter(phi >= 2/3 & variable == "P011")

path1_nonuniq <- plot_results %>%
  filter(Type == "Pathogen" & phi < 2/3)
path1_nonuniq$image <- "p1a.png"

path2_nonuniq <- plot_results %>%
  filter(Type == "Pathogen" & phi >= 2/3 & Genotype != "011")
path2_nonuniq$image <- "p2.png"

unique_results <- rbind(hostH_results, hostM_results, path_results)
nonunique_results <- rbind(nonunique_results, path1_nonuniq, path2_nonuniq)

plot_preds <- data.frame(Phi = seq(min(plot_results$phi), 
                                   max(plot_results$phi), length.out = 100)) %>%
  mutate(H001 = ifelse(Phi <= 2/3, (3 * Phi - 1) / (3 * Phi), NA),
         H010 = ifelse(Phi <= 2/3, 1 - H001, NA),
         M100 = ifelse(Phi <= 2/3, 1 / (3 - 3 * Phi), NA),
         M001 = ifelse(Phi <= 2/3, 1 - M100, NA),
         P011 = ifelse(Phi <= 2/3, 1/3, NA),
         P101 = ifelse(Phi <= 2/3, 1/3, NA),
         P110 = ifelse(Phi <= 2/3, 1/3, NA),
  ) %>%
  melt(id.vars = "Phi",
       variable.name = "Genotype",
       value.name = "value") %>%
  mutate(Type = substr(Genotype, 1, 1)) %>%
  mutate(Type = case_when(Type == "H" ~ "Host H",
                          Type == "M" ~ "Host M",
                          Type == "P" ~ "Pathogen")) %>%
  mutate(Genotype = substr(Genotype, 2, 4))


baseline_genos <- unique_results[1:8,]
baseline_genos$phi <- NA
baseline_genos$Genotype <- names(ini_path)

unique_results <- rbind(unique_results, baseline_genos)

plTheoryComp <- ggplot() +
  annotate("rect", xmin = 2/3, xmax = max(plot_results$phi), ymin = -Inf, ymax = Inf, 
           alpha = 0.3, fill = "lightblue") +
  geom_line(data = plot_preds, aes(x = Phi, y = value, color = Genotype), alpha = 0.5, size = 1) +
  geom_point(data = unique_results, aes(x = phi, y = value, color = Genotype, shape = Type),
             size = 4.5, alpha = 1) +
  geom_image(data = nonunique_results, aes(x = phi, y = value, image = image), size = 0.035) + # Adjust size as needed
  scale_color_manual(values =  c("000"="gray80",
                                 "100"="#D9565CFF",
                                 "010"="#F28A8AFF",
                                 "001"="brown4",
                                 "110"="#1BB6AFFF",
                                 "101"="#088BBEFF",
                                 "011"="#172869FF"
                                 #"111"="gray20"
                                )) +
  facet_grid( ~ Type) +
  theme_classic() +
  labs(x = expression("Relative proportion of host species H" ~ (phi)),
       y = "Frequency") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text = element_text(size=30),
        legend.title = element_text(size = 20),
        strip.text.x = element_text(size = 25),
        strip.text.y = element_text(size = 25),
        legend.text=element_text(size = 20),
        strip.background = element_blank(),
        panel.spacing = unit(2, "lines")) +
  guides(shape = "none")
plTheoryComp

jpeg("FigTheoryComp.jpeg", width = 16, height = 6, units = "in", res = 300)
plTheoryComp
dev.off()
