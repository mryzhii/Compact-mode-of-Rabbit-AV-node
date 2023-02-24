# Compact-mode-of-Rabbit-AV-node
This folder contains MATLAB source code for the Rabbit AV node model
from the paper:
"A compact multi-functional model of the rabbit atrioventricular 
node with dual pathways", Frontiers in Physiology, 14 (2023). 
DOI: 10.3389/fphys.2023.1126648

by M.Ryzhii (University of Aizu) and 
E.Ryzhii (Fukushima Medical University).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tested with MATLAB R2022b

rabbitAVN.m - main file
              Set 'Mode' to the corresponding value (see inside)

avn_data1.m - rabbit AVN model structure and parameters
avn_fun.m   - function for MATLAB ODE solver
avn_plot.m  - analyses ladder diagrams and creates plots
