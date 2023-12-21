{
  gROOT->ProcessLine(".L treeProducer_simulation.cpp++");
  gROOT->ProcessLine("treeProducer_simulation().Loop()");
  gROOT->ProcessLine(".q");
};
