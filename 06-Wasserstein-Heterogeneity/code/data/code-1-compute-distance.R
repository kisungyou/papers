# COMPUTE PAIRWISE W2 DISTANCE AMONG ALL IMAGES PER DIGIT

# macro-level setting -----------------------------------------------------
# switch local vs. remote
if (interactive()){                       # placeholder
  myid = 0
} else {
  args = commandArgs(trailingOnly = TRUE) # get numeric scenario number
  myid = as.numeric(args[1])              # and make it as number
  myid = as.integer(myid)                 # total of 10 settings; 1-10
}

# pilot run with subset
use_subset = FALSE

# id conversion to the digit label by subtracting 1
current_label = as.integer(myid-1)

# setup -------------------------------------------------------------------
# load the library and dataset
library(T4transport)
data("digits")

# filename to save
save_name = file.path(getwd(),"distances",paste0("digit_",current_label,".RData"))
if (file.exists(save_name)) {
  message("Output already exists: ", save_name)
  quit(save = "no", status = 0)
}

# subset selection for pilot run
if (use_subset){
  # tempoarary small one
  old_digit = digits
  set.seed(10)
  idsub = sample(1:5000, 1000)
  digits = list()
  digits$image = old_digit$image[idsub]
  digits$label = old_digit$label[idsub]
}

# helper: compute one digit's pairwise distance matrix in W2 --------------
compute_digit_pdmat <- function(digit_label) {
  target_images <- digits$image[digits$label == digit_label]
  n_target <- length(target_images)
  
  # precompute measures ONCE (this is huge)
  measures <- vector("list", length=n_target)
  for (it in 1:n_target){
    measures[[it]] = T4transport::img2measure(target_images[[it]])
  }
  
  pdmat <- matrix(0, n_target, n_target)
  for (i in 1:(n_target - 1)) {
    mi <- measures[[i]]
    for (j in (i + 1):n_target) {
      mj <- measures[[j]]
      d <- T4transport::wasserstein(mi$support, mj$support,
                                    wx = mi$weight, wy = mj$weight)$distance
      pdmat[i, j] <- d
      pdmat[j, i] <- d
    }
  }
  return(pdmat)
}


# execution ---------------------------------------------------------------
pdmat = compute_digit_pdmat(current_label) # compute
save(pdmat, file=save_name) # save
