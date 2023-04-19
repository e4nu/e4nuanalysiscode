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

NOTICE: as of now, the code only works at the gpvms. 

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

Each class has a specific role within the e4nu analysis:
- **ConfigueI**: it is responsible for the analysis configuration. The configuration is setup with an input txt file (see ConfigFiles/example_configuration.txt). All aspects of the analysis are configurable. It can be used to change the signal Topology definition, or turn off/on analysis cuts/requirements. 
- **BackgroundI**: it deals with all background substraction methods. 
- **AnalysisI**: it deals with analysis features common between data and MC. In particular:
  1. Requires valid event weights
  2. Electron kinematic cuts
  3. Cuts on Q2, W 
  4. Cooks the event - it removes all particles not specified in topology definition. For instance, in the case of a 1p0pi analysis (1p0pip0pim0pi00photons), it would remove neutrons, kaons, or other particles from the event. This simplifies the loops later in the analysis. 
- **MCCLAS6AnalysisI**: it deals with analysis features specific to MC data for CLAS6 analysis:
  1. Applies a minimum momentum cut on hadrons
  2. Smears hadrons kinematics
  3. Fiducials are taken care for
  4. Accepance weights are included 
- **CLAS6AnalysisI**: it deals with analysis features specific to CLAS6 data.
- **MCCLAS6StandardAnalysis** and **CLAS6StandardAnalysis**: they inherit from the interfaces. The standard classes simply call the `MCCLAS6AnalysisI::GetValidEvent(id)` functions. For analysis that differ from the standard one, a new class should be added with a new implementation of `GetValidEvent(id)`. The analysis is configured with the `AnalysisID` keyword. For now, only the standard analysis is available (analysis id of 0). 
- **E4NuAnalysis**: it is responsible to call either the MC or data objects according to the Configuration. It also deals with the **signal/bkg** selection. 

---------------
## User Guide
