library(tidyverse)
library(deSolve)
library(reshape2)
library(nleqslv)
library(gridExtra)

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

omega_H <- seq(0.025, 0.3, length.out = 20)
omega_M <- 0
omega_P <- seq(0.025, 0.3, length.out = 20)

phi <- c(0.5, 0.6, 0.65, 0.7, 0.8, 0.9)
sigma <- 0.85

in_pars <- crossing(beta_H = beta_H,
                      beta_P = beta_P,
                      xi_H = xi_H,
                      xi_P = xi_P,
                      omega_H = omega_H,
                      omega_P = omega_P,
                      phi = phi,
                      sigma = sigma)

in_pars$omega_M <- in_pars$omega_H



ini_host <- c(rep(1 / 4, 4), rep(1/4, 4))
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
coex_results <- data.frame()
for(i in 1:nrow(in_pars)) {
  
  if(i %% 100 == 0) print(paste("Run:", i, "out of", nrow(in_pars)))
  
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
  coex_results <- rbind(coex_results, final_state)
  
  #non_inv_ranges <- GetNonInvRanges(pars)
  #cur_ranges <- cbind(cur_pars, non_inv_ranges)
  #range_results <- rbind(range_results, cur_ranges)
}

print(Sys.time() - start_time)

extinct_threshold <- 1e-6

plot_results <- coex_results %>%
  dplyr::select(omega_H, omega_P, starts_with("H"), starts_with("M"), starts_with("P")) %>%
  melt(id.vars = c("omega_H", "omega_P", "phi")) %>%
  mutate(value = ifelse(value < extinct_threshold, 0, value)) %>%
  mutate(value = ifelse(value <= 0, NA, value))

plStrainFreqs <- ggplot(plot_results, aes(x = omega_H, y = omega_P, fill = value)) +
  geom_tile(color = "white", lwd = 0.5) + # Add white borders
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  facet_grid(phi ~ variable) +
  theme_classic() +
  labs(x = "Host Cost (omega_H)",
       y = "Pathogen Cost (omega_P)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(text = element_text(size=30),
        legend.title = element_blank(),
        strip.text.x = element_text(size = 25),
        strip.text.y = element_text(size = 25),
        legend.text=element_text(size = 25))

plStrainFreqs

jpeg("FigCoexistingStrains.jpeg", width = 48, height = 24, units = "in", res = 300)
plStrainFreqs
dev.off()

presence_results <- cbind(coex_results[,1:9],
                          coex_results[,10:ncol(coex_results)] > 1e-6)

presence_results$State <- NA
for(i in 1:nrow(presence_results)) {
  h_pres <- presence_results[i, 10:13]
  m_pres <- presence_results[i, 14:17]
  p_pres <- presence_results[i, 18:25]
  
  h_state <- colnames(out)[which(h_pres == TRUE) + 1]
  h_state <- substr(h_state, 2, nchar(h_state))
  if(length(h_state) == 0) {
    h_state <- "Extinct"
  } else {
    h_state <- paste(h_state, collapse = "-")
  }
  
  m_state <- colnames(out)[which(m_pres == TRUE) + 5]
  m_state <- substr(m_state, 2, nchar(m_state))
  
  if(length(m_state) == 0) {
    m_state <- "Extinct"
  } else {
    m_state <- paste(m_state, collapse = "-")
  }
  
  p_state <- colnames(out)[which(p_pres == TRUE) + 9]
  p_state <- substr(p_state, 2, nchar(p_state))
  
  if(length(p_state) == 0) {
    p_state <- "Extinct"
  } else {
    p_state <- paste(p_state, collapse = "-")
  }
  
  presence_results$HState[i] <- h_state
  presence_results$MState[i] <- m_state
  presence_results$PState[i] <- p_state
}

plot_pres <- presence_results %>%
  dplyr::select(omega_H, omega_P, phi, HState, MState, PState) %>%
  melt(id.vars = c("omega_H", "omega_P", "phi")) %>%
  mutate(variable = case_when(variable == "HState" ~ "Host H",
                              variable == "MState" ~ "Host M",
                              variable == "PState" ~ "Pathogen"))

plStates <- ggplot(plot_pres, aes(x = omega_H, y = omega_P, fill = value)) +
  geom_tile() +
  scale_fill_viridis_d(option = "magma") +
  facet_grid(variable ~ phi, labeller = label_bquote(cols = phi == .(phi))) +
  theme_classic() +
  labs(x = expression("Maximum cost of resistance" ~ (Omega[H/M])),
       y = expression("Maximum cost of virulence" ~ (Omega[P])),
       fill = "Coexisting Genotypes") +
  theme(text = element_text(size=30),
        title = element_text(size = 30, face = "bold"),
        legend.title = element_text(size = 20),
        strip.text.x = element_text(size = 25),
        strip.text.y = element_text(size = 25),
        legend.text=element_text(size = 20),
        strip.background = element_blank(),
        legend.position = "right",
        panel.spacing = unit(0, "lines")) +
  ggtitle("(a)")
plStates

range_data <- plot_pres %>%
  separate_longer_delim(col = value, delim = "-") %>%
  mutate(Range = str_count(value, "1")) %>%
  group_by(omega_H, omega_P, phi, variable) %>%
  count(Range) %>%
  mutate(Range = paste0("Range", Range)) %>%
  pivot_wider(names_from = Range, values_from = n)

range_data$Range0 <- ifelse(is.na(range_data$Range0), NA, 0)
range_data$Range1 <- ifelse(is.na(range_data$Range1), NA, 1)
range_data$Range2 <- ifelse(is.na(range_data$Range2), NA, 2)
range_data$Range3 <- ifelse(is.na(range_data$Range3), NA, 3)

#range_data[is.na(range_data)] <- ""
range_data <- range_data %>%
  unite(AllRange,
        Range0, Range1, Range2, Range3,
        sep = " + ", na.rm = TRUE)


plRanges <- ggplot(range_data, aes(x = omega_H, y = omega_P, fill = AllRange)) +
  geom_tile() +
  scale_fill_manual(values = c("0" = "#AFB42B",
                               "1" = "#FFA726FF",
                               "0 + 1" = "#F57C00FF",
                               "2" = "#FF95A8FF",
                               "0 + 2" = "#EC407AFF",
                               "1 + 2" = "#C2185BFF",
                               "0 + 1 + 2" = "#8A4198FF",
                               "3" = "mediumorchid3",
                               "0 + 3" = "lightblue",
                               "1 + 3" = "darkseagreen2",
                               "0 + 1 + 3" = "limegreen",
                               "2 + 3" = "#008EA0FF",
                               "0 + 2 + 3" = "darkturquoise",
                               "1 + 2 + 3" = "#1A5355FF",
                               "0 + 1 + 2 + 3" = "black")) +
  facet_grid(variable ~ phi, labeller = label_bquote(cols = phi == .(phi))) +
  theme_classic() +
  labs(x = expression("Maximum cost of resistance" ~ (Omega[H/M])),
       y = expression("Maximum cost of virulence" ~ (Omega[P])),
       fill = "# of resistance (R) alleles\n or virulence (V) alleles/\n in maintained haplotypes") +
  theme(text = element_text(size=30),
        title = element_text(size = 30, face = "bold"),
        legend.title = element_text(size = 20),
        strip.text.x = element_text(size = 25),
        strip.text.y = element_text(size = 25),
        legend.text=element_text(size = 20),
        strip.background = element_blank(),
        legend.position = "right",
        panel.spacing = unit(0, "lines")) +
  ggtitle("(b)")
plRanges


jpeg("FigCoexistingStates.jpeg", width = 22, height = 24, units = "in", res = 300)
grid.arrange(plStates, plRanges, nrow = 2)
dev.off()
