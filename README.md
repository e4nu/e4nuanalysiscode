# e4nu Analysis code
This repository contains the e4nu analysis code and utils to plot the relevant information
The following are brief descriptions of the code structure, various classes and configuration examples. It also provides with instructions for new e4nu analysers.
---------------
## Compile and run e4nu analyses
To configure the software:
1. Setup the enviromental variables in "e4nu_analysis_env.sh": `source e4nu_analysis_env.sh;`
2. Compile the code "make" from the E4NUANALYSIS directory: `cd $E4NUANALYSIS; make;`
4. If a previous installation exists, `cd $E4NUANALYSIS; make clean ; make;`
3. Run the main application, `./e4nuanalysis`

After each run, you will get a new root file containing information from the valid analised events. The information stored in this file, as well as it's name, can be configured from a configuration file (see Configuration section).
**NOTICE**: as of now, the code only works at the gpvms. 
---------------
## Directory Structure
The code can be found in different subdirectories, with the following structure: 
- **src** : source directory. Further described below.
- **data** : contains data files used to run analyses. There's two sub-directories, `AcceptanceMaps` and `FiducialsCorrections`, containing acceptance maps and fiducial corrections. As of now, only CLAS6 information is available. 
- **ConfigFiles** : it contains examples of configuration files. For instance, see `example_configuration.txt` 
- **PlottingScripts** : it contains examples of plotting scripts that can be used to plot the output of e4nu analysis. 

The structure of the src directory is as follows:
- **conf** : it contains most of the constants and configurables. Most of the constants depend on the configuration. It provides with the corresponding getter functions to access those. 
- **utils** : it provides with utils which can be used for any analysis. Some examples are `Fiducial.h`, `KinematicUtils.h` (with event kinematics definitions), and so on.  
- **physics** : it contains the implementation of `EventI`, and the corresponding specialitzations for data and MC, `CLAS6Event` and `MCEvent`. It also contains the `EventHolder` class and it's specialitzations. These are used to load the data from the input root files. 
- **analysis** : this directory contains the main classes used for analysis. The structure it's explained in the next section. 
- **apps** : contains the main executable
-  make : it contains some make files needed for compilation 
---------------
## E4nu Analysis Code
The e4nu analysis code is highly factorized into different classes, which inherit from each other. The main structure is depicted in the diagram below:
![e4nu diagram](https://github.com/e4nu/e4nuanalysiscode/blob/origin/Develop/RefactorizedCode/PlottingScripts/e4nuanalysis_diagram.png)

---------------
## User Guide
