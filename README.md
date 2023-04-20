# e4nu Analysis code

This repository contains the e4nu analysis code and utils to plot the relevant information
The following are brief descriptions of the code structure, various classes and configuration examples. It also provides with instructions for new e4nu analysers.

---------------

## Compile and run e4nu analyses

To configure the software:
1. Setup the environmental variables in "e4nu_analysis_env.sh": `source e4nu_analysis_env.sh;`
2. Compile the code "make" from the E4NUANALYSIS directory: `cd $E4NUANALYSIS; make;`
4. If a previous installation exists, `cd $E4NUANALYSIS; make clean ; make;`
3. Run the main application, `./e4nuanalysis`

After each run, you will get a new root file containing information from the valid analysed events. The information stored in this file, as well as it's name, can be configured from a configuration file (see Configuration section).

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
- **apps** : contains the main executable
-  make : it contains some make files needed for compilation 

---------------

## E4nu Analysis Code

The main analysis flow is as follows:
![e4nu flow](https://github.com/e4nu/e4nuanalysiscode/blob/origin/Develop/RefactorizedCode/PlottingScripts/e4nu_analysis_flow.png)

The e4nu analysis code is highly factorized into different classes, which inherit from each other. The main structure is depicted in the diagram below:
![e4nu diagram](https://github.com/e4nu/e4nuanalysiscode/blob/origin/Develop/RefactorizedCode/PlottingScripts/e4nuanalysis_diagram.png)

Each class has a specific role within the e4nu analysis:
- **ConfigueI**: it is responsible for the analysis configuration. The configuration is setup with an input txt file (see ConfigFiles/example_configuration.txt). All aspects of the analysis are configurable. It can be used to change the signal Topology definition, or turn off/on analysis cuts/requirements. 
- **BackgroundI**: it deals with all background subtraction methods. 
- **AnalysisI**: it deals with analysis features common between data and MC. In particular:
  1. Requires valid event weights
  2. Electron kinematic cuts
  3. Cuts on Q2, W 
  4. Cooks the event - it removes all particles not specified in topology definition. For instance, in the case of a 1p0pi analysis, it would remove neutrons, kaons, or other particles from the event. This simplifies the loops later in the analysis. 
- **MCCLAS6AnalysisI**: it deals with analysis features specific to MC data for CLAS6 analysis:
  1. Applies a minimum momentum cut on hadrons
  2. Smears hadrons kinematics
  3. Fiducials are taken care for
  4. Acceptance weights are included 
- **CLAS6AnalysisI**: it deals with analysis features specific to CLAS6 data.
- **MCCLAS6StandardAnalysis** and **CLAS6StandardAnalysis**: they inherit from the interfaces. The standard classes simply call the `MCCLAS6AnalysisI::GetValidEvent(id)` functions. For analysis that differ from the standard one, a new class should be added with a new implementation of `GetValidEvent(id)`. The analysis is configured with the `AnalysisID` keyword. For now, only the standard analysis is available (analysis id of 0). 
- **E4NuAnalysis**: it is responsible to call either the MC or data objects according to the Configuration. It also deals with the **signal/bkg** selection. 

---------------

## User Guide

New users might want to include analysis features specific to their analysis. This section depicts the best coding practices for this goal in the e4nuanalysis software.

- The following classes ***musn't be modified***: BackgroundI, AnalysisI, CLAS6AnalysisI and MCCLAS6AnalysisI. These classes **are completely generic**. If the user doesn't want to use a specific cut or effect, it should simply be turned off using the configuration file.
- If additional features have to be included, the user should use either the CLAS6StandardAnalysis and CLAS6StandardAnalysis classes as templates for their new class (with a speci
fic name for their analysis). In addition, these classes must be accordingly included in `E4NuAnalysis`. A new analysis ID must be asigned to these classes.
- ConfigureI should not be modified, unless the user wants to add new configurables.

---------------

## Configuration Guide 

Run configurables:
- **EBeam**: beam energy. It should match the root file content
- **TargetPdg**: pdg of the target used in the run
- **NEvents**: number of events to run in your analysis
- **FirstEvent**: first event to start runing from. It can be used for parallelization

Analysis topology definition:
- **IsData**: used to inform the software whether the data is experimental (true) or not (false)
- **IsCLAS6Analysis**: bool set to true or false. For now it can only be true as CLAS12 is not available yet
- **Toplogy**: used to define the topology of the analysis. The format must be pdg1:multiplicity1,pdg2:multiplicity2 and so on

Background subtraction method configurables:
- **MaxBackgroundMultiplicity**: maximum background multiplicity to consider in your background substraction method
- **NRotations**: number of rotations used in the background substraction method

Analysis cuts: set to true or false to turn on or off
- **ApplyPhiOpeningAngle**
- **ApplyMomCut**
- **UsePhiThetaBand**
- **ApplyThetaSlice**
- **UseAllSectors** 
- **ApplyFiducial**
- **ApplyAccWeights**
- **ApplyReso**
- **ApplyMottWeight**
- **ApplyGoodSectorPhiSlice**
- **ApplyOutEMomCut**
- **ApplyQ2Cut**
- **ApplyWCut**
- **SubtractBkg**

Histogram configurables:
- **RangeList**: min1:max1,min2:max2,..,minN:maxN
- **ObservableList**: obs1,obs2,...,obsN
- **NBinsList**: nb1,nb2,...,nbN
- **NormalizeHists**: set to true to normalize from event distribuiton to cross section
- **DebugBkg**: add background plots for debugging

Input and output files configurables:
- **InputFile**: path to input root files with events to analize
- **OutputFile**: output root files with analised events and histograms
- **XSecFile**: path to xml file for MC normalization



