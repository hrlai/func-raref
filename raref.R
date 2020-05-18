require(mvabund)
require(FD)
require(reshape)
require(dplyr)
require(foreach)
require(abind)

source("code/Chao_2014_FD_Hill_noBeta.R") # FD Hill number code from Chiu and Chao (2014)

data(antTraits)

ssm <- as.matrix(antTraits$abund)
rownames(ssm) <- rownames(antTraits$abund)

stm <- as.matrix(antTraits$traits[,c(1,5)])

# randomly assign origin to ant species
set.seed(101)
origin <- 
  data.frame(Species = colnames(ssm),
             Exotic = sample(c(0, 1), ncol(ssm), replace = TRUE))


# Species-species matrix filled with trait distances
# Using Gower's distance, but see Chiu & Chao (2014) PLoS ONE for possible downside (non-ultrametric)
get.D <- function(traits) {
  D <- gowdis(traits)
  if (is.euclid(D) == FALSE) D <- cailliez(D)
  D <- as.matrix(D)
  return(D)
}

D <- get.D(stm)


## First, decide on the extrapolation size
## look at the distribution of total abundance across sites
# hist(rowSums(ssm.all), breaks = 30)
# quantile(rowSums(ssm.all), seq(0, 1, 0.1))
## a few sites have <5 total abundance and these are pulling our lower limit down
## previously I didn't remove them, resulting in very low extrapolation number
## let's use the 10th percentile and remove extremely sparse sites
## so that we can extrapolate more and do a more "reliable" rarefaction--extrapolation
total.abun.lower.lim <- floor(as.numeric(quantile(rowSums(ssm), 0.1)))
raref.size <- 3 * total.abun.lower.lim
# hist(rowSums(ssm.all)); median(rowSums(ssm.all))
sparse.plots <- names(which(rowSums(ssm) < total.abun.lower.lim))
# remove sparse plots
ssm <- ssm[setdiff(rownames(ssm), sparse.plots), ]

n.raref <- 100

# species rarefaction/extrapolation
# Sample n = raref.size individuals from observed assemblages
doMC::registerDoMC(cores = 2)
acomb <- function(...) abind(..., along = 4)

SDFD.raref.raw <- 
  foreach (n = 1:n.raref, .combine = "acomb") %dopar% {
    tmp <- 
      apply(ssm, 1, function(x) {
        # sample individuals from each plot with replacement
        sample(rep(names(x), x), raref.size, replace = TRUE)
      })
    tmp <- 
      reshape::melt(tmp) %>% 
      rename(Plot = X2,
             Species = value) %>% 
      left_join(origin, by = "Species")
    
    # rarefied communities
    ssm.tmp.T <- with(tmp, tapply(Species, list(Species, Plot), length, default = 0))
    ssm.tmp.N <- with(tmp %>% filter(Exotic == 0), 
                      tapply(Species, list(Species, Plot), length, default = 0))
    ssm.tmp.E <- with(tmp %>% filter(Exotic == 1), 
                      tapply(Species, list(Species, Plot), length, default = 0))    
    ssm.tmp <- list(T = ssm.tmp.T, N = ssm.tmp.N, E = ssm.tmp.E)
    
    # calculate rarefied FD
    out <- 
      sapply(ssm.tmp, function(x) {
        # SD
        SD.raref.raw <- 
          cbind(SD0 = specnumber(t(x)),
                SD1 = exp(diversity(t(x), index = "shannon")),
                SD2 = diversity(t(x), index = "invsimpson"))
        # FD
        D.tmp <- D[rownames(x), rownames(x)]
        FD.raref.raw <- 
          do.call(cbind, 
                  lapply(c(0, 1, 2), function(q) {Func2014(D.tmp, x, q = q)$FuncD}))
        colnames(FD.raref.raw) <- paste0("FD", c(0, 1, 2))
        # combine SD and FD
        return(cbind(SD.raref.raw, FD.raref.raw))
      },
      simplify = "array")
    
    return(out)
  }

SDFD.raref <- 
  list(
    Mean = apply(SDFD.raref.raw, 1:3, function(x) mean(x[is.finite(x)], na.rm = TRUE)),
    LCI  = apply(SDFD.raref.raw, 1:3, function(x) quantile(x[is.finite(x)], probs = 0.025, na.rm = TRUE)),
    UCI  = apply(SDFD.raref.raw, 1:3, function(x) quantile(x[is.finite(x)], probs = 0.975, na.rm = TRUE)),
    SD   = apply(SDFD.raref.raw, 1:3, function(x) sd(x[is.finite(x)], na.rm = TRUE))
  )

