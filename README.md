# e4nu Analysis code

This repository contains the e4nu analysis code and utils to plot the relevant information.
The following are brief descriptions of the code structure and various classes. Additional information on the analysis can be found on the e4nu Wiki.

---------------

## Compile and run e4nu analyses

To configure the software:
1. Setup the environmental variables according to the farm you're in. For the FNAL gpvm, `source e4nu_gpvm_env.sh;`. For the ifarm, `source e4nu_ifarm_env.sh;`.
2. Compile the code "make" from the E4NUANALYSIS directory: `cd $E4NUANALYSIS; make;`
4. If a previous installation exists, `cd $E4NUANALYSIS; make clean ; make;`
3. Run the main application, `./e4nuanalysis`. You can pass a different configuration file by passing the name to `./e4nuanalysis my_config.txt`. The configuration file must live in ConfigFiles.

After each run, you will get a new root file containing information from the valid analysed events. The information stored in this file, as well as it's name, can be configured from a configuration file (see Configuration guide Wiki).

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
- **physics** : it contains the implementation of `EventI`, and the corresponding specializations for data and MC, `CLAS6Event` and `MCEvent`. It also contains the `EventHolder` class and it's specializations. These are used to load the data from the input root files. 
- **analysis** : this directory contains the main classes used for analysis. The structure it's explained in the next section. 
- **apps** : contains the main executable, `e4nuanalysis.cxx`
-  make : it contains some make files needed for compilation

---------------

## E4nu Analysis Code

The analysis app [e4nuanalysis.cxx](https://github.com/e4nu/e4nuanalysiscode/blob/master/src/apps/e4nuanalysis.cxx) is the main executable. It's content is simple: it instiantates an E4NuAnalysis object (see below details on inheritance chain), configured with `ConfigFiles/example_configuration`. This object is then used to load the data from an input root file ([analysis -> LoadData()](https://github.com/e4nu/e4nuanalysiscode/blob/master/src/apps/e4nuanalysis.cxx#L28)), run the analysis and background subtraction methods, and finally, store the analized information in a TTree as well as Histograms.

The main analysis flow is as follows:
![e4nu flow](https://github.com/e4nu/e4nuanalysiscode/blob/master/PlottingScripts/e4nu_analysis_flow.png)

The e4nu analysis code is highly factorized into different classes, which inherit from each other. ConfigureI is the base class, E4NuAnalyisis is the top derived class. The main structure is depicted in the diagram below:
![e4nu diagram](https://github.com/e4nu/e4nuanalysiscode/blob/master/PlottingScripts/e4nuanalysis_diagram.png)
The (...) boxes indicate that new classes might be added to accomodate new analysis cuts, specific to a new analysis (see user gide section below).

Each class has a specific role within the e4nu analysis:
- **ConfigueI**: it is responsible for the analysis configuration. The configuration is setup with an input txt file (see ConfigFiles/example_configuration.txt). All aspects of the analysis are configurable. It can be used to change the signal Topology definition, or turn off/on analysis cuts/requirements. 
- **BackgroundI**: it deals with all background subtraction methods. More information on the available configuration setups can be found in Configuration Guide.
- **AnalysisI**: it deals with analysis features common between data and MC. In particular:
  1. Requires valid event weights
  2. Electron kinematic cuts
  3. Cuts on Q2, W
  4. Applies a minimum momentum cut on hadrons 
  5. [Cooks the event](https://github.com/e4nu/e4nuanalysiscode/blob/master/src/analysis/AnalysisI.cxx#L106) - it removes all particles not specified in topology definition. For instance, in the case of a 1p0pi analysis, it would remove neutrons, kaons, or other particles from the event. This simplifies the loops later in the analysis. 
- **MCCLAS6AnalysisI**: it deals with analysis features specific to MC data for CLAS6 analysis:
  1. Smears hadrons kinematics
  2. Fiducials are taken care for
  3. Acceptance weights are computed 
- **CLAS6AnalysisI**: it deals with analysis features specific to CLAS6 data.
- **MCCLAS6StandardAnalysis** and **CLAS6StandardAnalysis**: they inherit from the MCCLAS6AnalysisI and CLAS6AnalysisI interfaces. The standard classes are templates to facilitate the integration of new analysis by new users. For this reason, the standard classes simply call the `MCCLAS6AnalysisI::GetValidEvent(id)` or `CLAS6AnalysisI::GetValidEvent(id)` functions. For analysis that differ from the standard one, a new class should be added with a new implementation of `GetValidEvent(id)`, using these classes as templates. The analysis is configured with the `AnalysisID` keyword. For now, only the standard analysis is available (analysis id of 0). New analysis would require a new analysis ID.
- **E4NuAnalysis**: this is the main class used for analysis. It is responsible to call either the MC or data objects according to the Configuration. It also deals with the **signal/bkg** selection. An E4NuAnalysis object is defined in `e4nuanalysis.cxx` using a configuration file. 

