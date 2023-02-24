# Compact-mode-of-Rabbit-AV-node
This folder contains MATLAB source code for the Rabbit AV node model
from the paper:<br>
"A compact multi-functional model of the rabbit atrioventricular 
node with dual pathways",<br>
Frontiers in Physiology, 14 (2023). <br>
DOI: 10.3389/fphys.2023.1126648

by M. Ryzhii (University of Aizu) and 
E. Ryzhii (Fukushima Medical University).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%<br>
Tested with MATLAB R2022b<br>

rabbitAVN.m - main file<br>
              Set 'Mode' to the corresponding value (see inside)<br>
avn_data1.m - rabbit AVN model structure and parameters<br>
avn_fun.m   - function for MATLAB ODE solver<br>
avn_plot.m  - analyses ladder diagrams and creates plots<br>
