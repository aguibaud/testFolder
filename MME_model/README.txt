A City-Scale Model of Fire Development in the Mediterranean and Middle Eastern Cities 
======================================================

Developed by:
Yonatan Shaham (yjonas83@gmail.com)
Itzhak Benenson (bennya@post.tau.ac.il)

License: GPLv3.0
see: http://www.gnu.org/licenses/gpl-3.0.txt


Installation
-------------

1. Download and install GAMA-platform: https://github.com/gama-platform/gama/wiki
2. Under the "Gama Project" window, right-click on "user models"->"import"->"existing project into workspace"
3. Select the MME model folder
4. Check "copy project into workspace"

*Files compatibility was tested for gama version 1.7RC on Windows 7
 
Running the Model
----------------------
Use the "MME model" file.
There are 2 experiments:
1. Regular display: 1 ignition in a random apartment, full GUI.
2. Serial_apts: going through all ground floor apartments, 1 ignition in each.

Parameters can be set for each experiment, see:
https://github.com/gama-platform/gama/wiki/DefiningParameters
https://github.com/gama-platform/gama/wiki/BatchExperiments

Running the data generation algorithm
-------------------------------------------
Use the "MME data generation algorithm" file.
Run the experiment named "running_the_algorithm".
One building will be processed in each simulation step.
Warnings are likely to appear and do not indicate a problem. They are created since the algorithm uses surplus polygons.
Output files prefix can be set using parameters window.

Input file should include the following fields:
1. BLDG_HT: building height in meters.
2. L: fuel load (wood kg/m2)
3. roof: 0 - roof in non-flammable, 1 - room is flammable

DO NOT USE NON-PROJECTED FILES (such as WGS:84)







