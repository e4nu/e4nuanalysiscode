import os, optparse
import ROOT

op = optparse.OptionParser(usage=__doc__)
op.add_option("--file", dest="file", help="ROOT files from which to read the Tree")
op.add_option("--output-name", dest="output_name", default="e4nu_analysis_comparison",help="Output name for the pdf file. Default %default")
op.add_option("--x-label", dest="xlabel", default="", help="X Axis label")
op.add_option("--y-label", dest="ylabel", default="", help="Y Axis label")
op.add_option("--plot-as-graph", dest="graph", default=False, action="store_true", help="Store histograms as graphs")
op.add_option("--is-MC", dest="MC",default=False, action="store_true", help="Is MC data")
opts, args = op.parse_args()

# Open the root file
file = ROOT.TFile.Open(opts.file)

canvas = ROOT.TCanvas("canvas", "Histogram Canvas", 800, 600)

tree_name = "CLAS6Tree"
if opts.MC : 
    tree_name = "MCTree"

tree = file.Get(tree_name)
hist_W = ROOT.TH1D("W", "W", 80,0.6,1.5)
hist_aT = ROOT.TH1D("AlphaT", "AlphaT", 80,0.,180)
hist_pimom = ROOT.TH1D("Pion_momenta", "Pion Momenta", 40,0.1,0.7)
hist_pitheta = ROOT.TH1D("Pion_angle", "Pion Theta", 40, 0, 150 ) 

# Loop over the entries in the tree
for entry in tree:
    TotWeight = entry.TotWeight 
    IsBkg = entry.IsBkg 
    EM = entry.EM
    QEL = entry.QEL
    RES = entry.RES
    MEC = entry.MEC
    DIS = entry.DIS
    TrueW = entry.TrueW
    pfl = entry.pfl
    pfl_phi = entry.pfl_phi
    pfl_theta = entry.pfl_theta
    pip_mom = entry.pip_mom
    pip_theta = entry.pip_theta
    alphaT = entry.AlphaT
    #hist_W.Fill(TrueW)
    #if TotWeight == 0 : continue
    hist_aT.Fill(alphaT)#,TotWeight)

#hist_W.Draw()

hist_aT.Draw("hist")
canvas.Print(opts.output_name+".pdf")
file.Close()
    
