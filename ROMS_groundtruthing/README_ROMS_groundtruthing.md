## README ROMS grountruthing folder
Date last edited: 10/30/2020

### This folder contains R script files used to compare and contrast the ROMS data to in situ oceanographic data such as water temperature, current and wind velocity, as well as comparing the 2 km and 250 m ROMS. 

NOTE: The main goal is to create one main script with each analysis as a modularized function. 
--- 

### Script names:
- Data_comparison_insitu_ROMS.R
- Particle_traj_comparison.R
- Plotting_Groundtruthing_Data.R

--- 

### Descriptions 

#### Data_comparison_insitu_ROMS.R 

This R script downloads and processes in situ oceanographic data such as water temperature, salinity, and current velocity to the same ROMS variables.


#### Particle_traj_comparison.R

This script extracts particles vertical velocity over the entire simualtion to compare if particles are experiencing different average/sum vertical velocities. 


#### Plotting_Groundtruthing_Data.R

This R script creates the figures to compare and contrast the ROMS and in situ data. 

