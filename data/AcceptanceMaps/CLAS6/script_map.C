{
gROOT->ProcessLine(".L ntuple.C+");
gROOT->ProcessLine("ntuple().Loop()");
gROOT->ProcessLine(".q");
};
