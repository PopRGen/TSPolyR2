# Hanna Maerkle
# Plot to visualise genotype cost functions from Ashby and Boots 2017

require("tidyverse")

# Note function will only work if the vectors of all input variables have equal length

genotype_cost_shape <- function(cx1,cx2,grange){
  if(any(grange > 1)|any(grange < 0)) {stop("Range can only take values between 0 and 1")}
  cost <- 1-cx1*(1-exp(cx2*grange))/(1-exp(cx2))
  return(cost)
}

# Function for linear costs
genotype_cost_lin <- function(cx1,grange){
  cost <- 1 - cx1*grange
  return(cost)
}

# General function
genotype_cost <- function(cx1,cx2,grange){
  out <- numeric()
  noz <- which(cx2!=0)
  
  out[noz] <- genotype_cost_shape(cx1[noz],cx2[noz],grange[noz])
  out[-noz] <- genotype_cost_lin(cx1[-noz],grange[-noz])
  return(out)
}

# Defining cs1 and cs2
cs1 <- seq(0.1,0.3,by=0.1)
cs2 <- seq(-3,3,by=0.1)
grange <- seq(0,1,by=0.01)

dat <- expand.grid(cs1=cs1,cs2=cs2,grange=grange)

dat <- dat %>% 
  mutate(cost=genotype_cost(cs1,cs2,grange))

cx_names <- paste0("Omega==",cs1)
cx_names <- setNames(cx_names,as.character(cs1))


tplt <-ggplot(dat,aes(grange,cost,group=cs2)) + 
  facet_wrap(vars(cs1), labeller = labeller(cs1  = as_labeller(cx_names,  label_parsed))) 

tplt <- tplt + geom_line(aes(color=cs2)) + 
  theme_bw() + 
  theme(text = element_text(size=16), panel.grid.minor = element_blank()) +
  scale_color_gradient2(low = "deeppink4", mid = "white",
                        high = "deepskyblue3", name=expression(paste("Shape\nparameter ", Xi))) + 
  ylab("Relative genotype fitness") + 
  xlab("Proportion of loci which have resistance (R) or virulence (V) allele") + 
  theme(panel.spacing.x = unit(0.5, "cm")) 

tplt

pdf("Suppl1_fitness_functions.pdf",width=10, height = 3)
  print(tplt)
dev.off()


# The simple version

# Defining cs1 and cs2
cs1_simple <- 0.3 #seq(0.2,0.3,by=0.1)
cs2_simple <- seq(-3,3,by=3)
grange <- seq(0,1,by=0.01)

dat_simple <- expand.grid(cs1_simple=cs1_simple,cs2_simple=cs2_simple,grange=grange)

dat_simple <- dat_simple %>% 
  mutate(cost=genotype_cost(cs1_simple,cs2_simple,grange))

cx_names <- paste0("c[H]^(1)==",cs1_simple)
cx_names <- setNames(cx_names,as.character(cs1_simple))



intercept_dat <- data.frame(cs1_simple=unique(dat_simple[["cs1_simple"]])) %>%
  mutate(value=1-cs1_simple)

illranges <- seq(0,1,by=1/3)

illudat <- expand.grid(cs1_simple=cs1_simple,cs2_simple=cs2_simple,illranges=illranges) %>%
  mutate(cost=genotype_cost(cs1_simple,cs2_simple,illranges))

  
  
tplt_simple <-ggplot(dat_simple,aes(grange,cost,group=cs2_simple)) + 
  facet_wrap(vars(cs1_simple), labeller = labeller(cs1_simple  = as_labeller(cx_names,  label_parsed))) 

tplt_simple <- tplt_simple + 
  geom_line(aes(color=factor(cs2_simple)),linewidth=1.5) + 
  theme_bw() + 
  theme(panel.grid=element_blank(),text = element_text(size=16))

tplt_simple <- tplt_simple + 
  scale_color_manual(values=c("darkred","grey","deepskyblue3"), name=parse(text="c[H]^(2)")) + 
  ylab("Relative fitness") + 
  #ylab(expression(paste("Relative fitness (1-"~c[S],")"))) + 
  #xlab("Resistance range\n(# resistance alleles/# resistance loci)") 
  #xlab(expression(paste("Resistance range: ", frac(foo, bar))))
  xlab("Proportion of loci with R-alleles") + 
  geom_hline(data=intercept_dat,aes(yintercept = value),linetype=3) +
  geom_hline(yintercept = 1,linetype=3) +
  ylim(0,1) +
  geom_point(data=illudat,aes(illranges,cost,group=cs2_simple,color=factor(cs2_simple)),size=2)

  
tplt_simple

pdf("Genotype_fitness_cost_AshbyBoots2017_simple.pdf",height=4,width=6)
print(tplt_simple)
dev.off()


