# Main scripts
Several scripts are available. Each of them have different purposes. 
- `plot_e4nuanalysis`: it is used to plot the output of the `e4nuanalysis` output ROOT file in a common format. 
- `plot_radiative`: (e,e') events will radiate photons due to bremstrahlung and other QED effects. This effectively modifies the born cross section. Effectively, this can be taken into account with effective corrections that depend on the event kinematics and topology. This script plots the corrections as a function of different kinematic variables. It is used to study the relative importance and effect of radiative corrections. 
- `plot_bkg`: it stores the background plots in a nice format
- `plot_comparison`: compares several histograms at the event rate level