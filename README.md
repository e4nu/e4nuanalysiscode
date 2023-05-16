# e4nu Analysis code

This repository contains the e4nu analysis code and utils to plot the relevant information.
The following are brief descriptions of the code structure, various classes and configuration examples. It also provides with instructions for new e4nu analysers.

---------------

## Compile and run e4nu analyses

To configure the software:
1. Setup the environmental variables according to the farm you're in. For the FNAL gpvm, `source e4nu_gpvm_env.sh;`. For the ifarm, `source e4nu_ifarm_env.sh;`.
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

---------------

## Accessing Event Information

Data is loaded from root files to Event Holders. [EventHolderI](https://github.com/e4nu/e4nuanalysiscode/blob/master/src/physics/EventHolderI.h) is the base class. There is a derived class for CLAS6 data and one for MC data, [CLAS6EventHolder.h](https://github.com/e4nu/e4nuanalysiscode/blob/master/src/physics/CLAS6EventHolder.h) and [MCEventHolder.h](https://github.com/e4nu/e4nuanalysiscode/blob/master/src/physics/MCEventHolder.h) respectively. From the holders, you can get an event provided the event id. The corresponding CLAS6 and MC data holders are instiantated in `CLAS6AnalysisI.h` and `MCCLAS6AnalysisI.h`.  `E4NuAnalysis` deals with the loop over all events, see [Analyse(void)](https://github.com/e4nu/e4nuanalysiscode/blob/master/src/analysis/E4NuAnalysis.cxx#L76). This is the main function called in `e4nuanalyser.cxx`.

An event is defined with the [EventI class](https://github.com/e4nu/e4nuanalysiscode/blob/maser/src/physics/EventI.h). There's two specialitzations of this base class, one for MC events, [MCEvent](https://github.com/e4nu/e4nuanalysiscode/blob/master/src/physics/MCEvent.h), and one for CLAS6 data events, [CLAS6Event](https://github.com/e4nu/e4nuanalysiscode/blob/master/src/physics/CLAS6Event.h).
The main difference between `CLAS6Event` and `MCEvent` members is that `MCEvent` has true level GENIE information, such as the type of interaction, or before FSI kinematics. These are not accessible for `CLAS6Event`s, hence the different specialitzation. 

The base class, `EventI`, contains all the relevant information for analysis: 
- Event ID
- Target pdg code
- Incoming and Outcoming lepton pdg codes
- Initial lepton kinematics (TLorentzVector)
- Final lepton kinematics (TLorentzVector)
- Final hadrons kinematics: map<int,vector< TLorentzVector > >. The key is the particle pdg code. For each pdg, there's a vector containing all Lorentz Vector. For instance, if there's two pions, there will be two entries with pdg code 211. Pdg codes can easily be accessed using the [ParticleI members](https://github.com/e4nu/e4nuanalysiscode/blob/master/src/conf/ParticleI.h). For instance, to get the number of pions in our event, we can simply do `(event->GetFinalParticles4Mom())[conf::kPdgPiP].size()`. 
- The same variables as above are stored for the uncorrected kinematics (before we take into account smearing effects). These are only relevant for MC data. 
- Number of reconstructed hadrons (for protons, pions, and so on)
- Mott weight
- Event weight 
- Acceptance weight (only relevant for MC)
- Total Event weight : total weight = event weight * acceptance weight * mott weight 
- Event multiplicity : number of topology particles in event. For instance, for a 1p0pi analysis, a 1p0pi event has multiplicity 1. A 2p0pi event has multiplicity 2. It should only be used after the event has been cooked (this is handled in AnalysisI). 

The corresponding getter and setter functions are available to access and modify an event during the analysis. 

---------------

## User Guide
New analysis can be performed by ***simply modifying the configuration file***, where the topology definition can be adjusted to the new analysis. 

E4nu users might want to ***include analysis features specific to their analysis***. This section depicts the best coding practices for this goal in the e4nuanalysis software.

- **Musn't modify**: The following classes ***musn't be modified***: BackgroundI, AnalysisI, CLAS6StandardAnalysis, MCCLAS5StandardAnalysis, CLAS6AnalysisI and MCCLAS6AnalysisI. These classes **are completely generic**. If the user doesn't want to use a specific cut or effect, this should simply be turned off using the configuration file. 
- **To modify**: If additional features have to be included, the user should use either the MCCLAS6StandardAnalysis and CLAS6StandardAnalysis classes as templates for their new class (with a specific name for their analysis). For instance, if the user wants to write new features specific to a Transparency measurement, for instance, these should be specified in two new classes, such as MCCLAS6TransparencyAnalysis and CLAS6TransparencyAnalysis. In addition, the new classes must be included in `E4NuAnalysis`, as `E4NuAnalysis` will inherit from those. A new analysis ID must be asigned to these classes to be able to configure it with the input configuration file.
- **Can modify**: ConfigureI should not be modified, unless the user wants to add new configurables. The main executable structure, `e4nuanalysis.cxx`, musnt't be modified, as it calls the analysis functions in a specific order. However, the user might want to use a different configuration file. It is possible to change the configuration file name from [there](https://github.com/e4nu/e4nuanalysiscode/blob/721dd5f41d51f4827630165cb8f86bac7c127865/src/apps/e4nuanalysis.cxx#L25). 

---------------

## Configuration Guide 

***Run configurables***:
- **EBeam**: beam energy. It should match the root file content
- **TargetPdg**: pdg of the target used in the run
- **NEvents**: number of events to run in your analysis
- **FirstEvent**: first event to start runing from. It can be used for parallelization

***Analysis topology definition***:
- **IsData**: used to inform the software whether the data is experimental (true) or not (false)
- **IsCLAS6Analysis**: bool set to true or false. For now it can only be true as CLAS12 is not available yet
- **Toplogy**: used to define the topology of the analysis. The format must be pdg1:multiplicity1,pdg2:multiplicity2 and so on

***Background subtraction method configurables***:
- **MaxBackgroundMultiplicity**: maximum background multiplicity to consider in your background substraction method
- **NRotations**: number of rotations used in the background substraction method
- **SubtractBkg**: bool. If true, the background substraction method is used. 

***AnalysisI cuts***: set to true or false to turn on or off
- **ApplyPhiOpeningAngle**: see [line](https://github.com/e4nu/e4nuanalysiscode/blob/e1669032a67c265d7725fc78678ec6515b966580/src/analysis/AnalysisI.cxx#L68).
- **ApplyThetaSlice**: the limits for electron angle are defined [here](https://github.com/e4nu/e4nuanalysiscode/blob/e1669032a67c265d7725fc78678ec6515b966580/src/conf/AnalysisConstantsI.h#L24).
- **UseAllSectors**: if false, only some sectors are used, see [DetectorUtils](https://github.com/e4nu/e4nuanalysiscode/blob/e1669032a67c265d7725fc78678ec6515b966580/src/utils/DetectorUtils.cxx#L48) for more details.
- **ApplyFiducial**: used to turn off or on fiducials. Fiducials are used in the [MC analysis](https://github.com/e4nu/e4nuanalysiscode/blob/e1669032a67c265d7725fc78678ec6515b966580/src/analysis/MCCLAS6AnalysisI.cxx#L141) as well as the background subtraction method.
- **ApplyGoodSectorPhiSlice**: see [conf::GoodSectorPhiSlice(double phi)](https://github.com/e4nu/e4nuanalysiscode/blob/e1669032a67c265d7725fc78678ec6515b966580/src/conf/AnalysisCutsI.cxx#L42)
- **ApplyOutEMomCut**: the limits are defined [here](https://github.com/e4nu/e4nuanalysiscode/blob/e1669032a67c265d7725fc78678ec6515b966580/src/conf/AnalysisConstantsI.h#L15).
- **ApplyQ2Cut**: [cut on Q2](https://github.com/e4nu/e4nuanalysiscode/blob/e1669032a67c265d7725fc78678ec6515b966580/src/conf/AnalysisCutsI.cxx#L52) which depends on the beam energy. 
- **ApplyWCut**: [cut on W](https://github.com/e4nu/e4nuanalysiscode/blob/e1669032a67c265d7725fc78678ec6515b966580/src/conf/AnalysisCutsI.cxx#L67).
- **ApplyMomCut**: it considers a minimum momentum cut for hadrons. The cut depends on the pdg of the hadron and the beam energy. You can find the exact values [here](https://github.com/e4nu/e4nuanalysiscode/blob/e1669032a67c265d7725fc78678ec6515b966580/src/conf/AnalysisCutsI.cxx#L13), and it's [implementation](https://github.com/e4nu/e4nuanalysiscode/blob/882df6efff649773bdd17b25726fa962058b0141/src/analysis/AnalysisI.cxx#L121).

***MCCLAS6AnalysisI Cuts***: set to either true or false to turn on or off
- **ApplyAccWeights**: used to consider efficiency maps. The maps location is `$E4NUANALYSIS/data/AcceptanceMaps/CLAS6/`. They depend on the beam energy, target and hadron pdg (see [AccpetanceMapsI](https://github.com/e4nu/e4nuanalysiscode/blob/e1669032a67c265d7725fc78678ec6515b966580/src/conf/AccpetanceMapsI.cxx#L14)). It is implemented in [MCCLAS6AnalysisI](https://github.com/e4nu/e4nuanalysiscode/blob/e1669032a67c265d7725fc78678ec6515b966580/src/analysis/MCCLAS6AnalysisI.cxx#L175).
- **ApplyReso**: used to smear the particles momentum. It only depends on the hadron pdg. You can find the values used [here](https://github.com/e4nu/e4nuanalysiscode/blob/e1669032a67c265d7725fc78678ec6515b966580/src/conf/ParticleI.h#L31). See the implementation [here](https://github.com/e4nu/e4nuanalysiscode/blob/e1669032a67c265d7725fc78678ec6515b966580/src/analysis/MCCLAS6AnalysisI.cxx#L198).
- **ApplyMottWeight**: it scales the events by the [Mott Scaling](https://github.com/e4nu/e4nuanalysiscode/blob/e1669032a67c265d7725fc78678ec6515b966580/src/physics/MCEvent.cxx#L38) to correct for the different coupling

It is also possible to change the configuration to use GENIE information before FSI effects. To do so, simply do:
- **No FSI** true
To use only true signal events, 
- **TrueSignal** true 

***Histogram configurables***:
- **RangeList**: min1:max1,min2:max2,..,minN:maxN
- **ObservableList**: obs1,obs2,...,obsN
- **NBinsList**: nb1,nb2,...,nbN
- **NormalizeHists**: set to true to normalize from event distribution to cross section
- **DebugBkg**: add background plots for debugging

You can find the available observables [here](https://github.com/e4nu/e4nuanalysiscode/blob/e029793c6e445fe2179e42a30e3c55eeaf1af980/src/physics/EventI.cxx#L149)

***Input and output files configurables***:
- **InputFile**: path to input root files with events to analize
- **OutputFile**: output root files with analised events and histograms
- **XSecFile**: path to xml file for MC normalization
- **ComputeAccCorr**: if true, it also computes the true and reconstructed signal spectra. It is used for the acceptance correction calculation.


