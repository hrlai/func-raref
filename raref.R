require(mvabund)
data(antTraits)

# randomly assign origin to ant species



# Species-species matrix filled with trait distances
# Using Gower's distance, but see Chiu & Chao (2014) PLoS ONE for possible downside (non-ultrametric)
D.all <- gowdis(stm.all)
if (is.euclid(D.all) == FALSE) D.all <- cailliez(D.all)
D.all <- as.matrix(D.all)
D <- lapply(ssm.list, function(mat) {
  temp <- D.all[colnames(mat), colnames(mat)]
  return(temp)
})

## First, decide on the extrapolation size
## look at the distribution of total abundance across sites
# hist(rowSums(ssm.all), breaks = 30)
# quantile(rowSums(ssm.all), seq(0, 1, 0.1))
## a few sites have <5 total abundance and these are pulling our lower limit down
## previously I didn't remove them, resulting in very low extrapolation number
## let's use the 10th percentile and remove extremely sparse sites
## so that we can extrapolate more and do a more "reliable" rarefaction--extrapolation
total.abun.lower.lim <- floor(as.numeric(quantile(rowSums(ssm.all), 0.1)))
raref.size <- 3 * total.abun.lower.lim
# hist(rowSums(ssm.all)); median(rowSums(ssm.all))
sparse.plots <- names(which(rowSums(ssm.all) < total.abun.lower.lim))
# remove sparse plots
ssm.list <- lapply(ssm.list, function(x) x[setdiff(rownames(x), sparse.plots), ])

# species rarefaction/extrapolation
# Sample n = raref.size individuals from observed assemblages
doMC::registerDoMC(cores = 2)
acomb <- function(...) abind(..., along = 4)

SDFD.raref.raw <- 
  foreach (n = 1:n.raref, .combine = "acomb") %dopar% {
    tmp <- 
      apply(ssm.list$T, 1, function(x) {
        # sample individuals from each plot with replacement
        sample(rep(names(x), x), raref.size, replace = TRUE)
      })
    tmp <- 
      reshape::melt(tmp) %>% 
      rename(Plot = X2,
             Species = value) %>% 
      left_join(select(use_spp, Species:Exotic))
    
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
        D.tmp <- D$T[rownames(x), rownames(x)]
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

