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