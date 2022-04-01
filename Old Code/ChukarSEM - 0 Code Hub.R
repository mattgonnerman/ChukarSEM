#Set Seed for Model Calibration
set.seed(2)

# Perform initial data cleanup and organization
source("./ChukarSEM - 1 Data Prep.R")

### Run initial model for Splines
# Model Code
source("./ChukarSEM - 2 NIMBLE Model.R")
# Prepare Data for Nimble
source("./ChukarSEM - 2a NIMBLE Prep.R")
# Run the Model
source("./ChukarSEM - 2b Run Model.R")
# Assess Model outputs
source("./ChukarSEM - 2c Assess Outputs.R")


### Run Final Model
# Model Code
source("./ChukarSEM - 3 NIMBLE Model.R")
# Prepare Data for Nimble
source("./ChukarSEM - 3a NIMBLE Prep.R")
# Run the Model
source("./ChukarSEM - 3b Run Model.R")
# Assess Model outputs
source("./ChukarSEM - 3c Assess Outputs.R")
