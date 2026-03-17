# GATHER ALL COMPUTED OBJECTS INTO A SINGLE DATAFILE


# setup -------------------------------------------------------------------
rm(list=ls())
library(rstudioapi)
library(T4transport)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load the original data --------------------------------------------------
# call the data
data(digits)

# empty list for saving the images and objects
list_images = list()
list_pdmats = list()


# iterate -----------------------------------------------------------------
for (i in 1:10){
  now_digit  = as.integer(i-1) # current digit
  load(file.path(getwd(),"distances",paste0("digit_",now_digit,".RData"))) # load
  
  list_images[[i]] = digits$image[digits$label == now_digit]
  list_pdmats[[i]] = pdmat
}

# save --------------------------------------------------------------------
save_name = paste0("per-digit-",nrow(pdmat),".RData")
save(list_images, list_pdmats, file=save_name)
