## README LPT model output analysis folder
Date last edited: 10/30/2020

### This folder contains R script files to analyze the model output from the LPT model 

NOTE: The main goal is to create one main script with each analysis as a modularized function
---

### Script names

- Post_Process_LPT_Model_Output.R
- mainPP_analysis.R
- Dispersal_and_Retention_Analyses.R
- SR_test.R
- ConnectivityMatrices.R

--- 

### Descriptions 

#### mainPP_Analysis.R

This R script will eventually become the main script that calls up the other analyses as modularized functions. 

#### Post_Process_LPT_Model_Output.R

This R script is used to make a first pass at the data to clean it. This script cleans the data for the following:
- unlists the data and converts it into a data.frame object
- removes any portion of a particles trajectory that leaves the spatial domain of the ROMS model. 


#### Dispersal_and_Retention_Analyses.R

This R script is used to create maps of the below metrics:
- Frequency of Occurrence (FO)
- Particle Import (IMP)
- Retention
- Average Distance Traveled in the Longitudinal and Latitudinal Direction (Dlat, Dlon)
- Depth vs. Time
- Latutide vs. Time
- Time spent in the model domain for the entire simulation


#### ConnectivityMartices.R

This R script is used to create connectivity matrices that are able to show if particles are being sourced and or sunk from the same or different region. 


#### SR_test.R

This R script is used to create a barplot showing the number of particles that return to the same region they are sourced from, but for only particles that settle in locations with maximum depts of 100 m. 

