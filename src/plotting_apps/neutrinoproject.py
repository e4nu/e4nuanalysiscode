from mpl_toolkits import mplot3d
import ROOT
from ROOT import TMath
from ROOT import gPad
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfFile, PdfPages
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import sys
import time
import os
from numpy import double
from fpdf import FPDF
from pypdf import PdfMerger
from ctypes import *


# Defining observable dictionary outside
observable_dict = {}
observable_dict["TotWeight"] = c_double(0)
observable_dict["RecoW"] = c_double(0)
observable_dict["ECal"] = c_double(0)
observable_dict["Recoq3"]= c_double(0)
observable_dict["pfl"] = c_double(0)
observable_dict["pfl_theta"] = c_double(0)
observable_dict["pfl_phi"] = c_double(0)
observable_dict["proton_mom"] = c_double(0)
observable_dict["proton_phi"] = c_double(0)
observable_dict["proton_theta"] = c_double(0)
observable_dict["HadAlphaT"] = c_double(0)
observable_dict["HadDeltaPT"] = c_double(0)
observable_dict["HadDeltaPhiT"] = c_double(0)
observable_dict["ElectronSector"] = int(-99)

def compute_event_rate( file_name, hist_name, tree_name, observable_x, observable_y, nbins_x, min_val_x, max_val_x, nbins_y, min_val_y, max_val_y ):
    # We want this function to be able to open any root file and loop over the events
    rootfile = ROOT.TFile.Open(file_name, "READ")
    if not rootfile:
        print("ERROR: the", file_name," does not exist.")
        return
    tree = rootfile.Get(tree_name)  # This will be different for data
    hist = ROOT.TH2D( hist_name, "", nbins_x, min_val_x, max_val_x, nbins_y, min_val_y, max_val_y)  # Need to be able to set title

    tree.SetBranchAddress( "TotWeight", observable_dict["TotWeight"] )
    tree.SetBranchAddress( "RecoW", observable_dict["RecoW"] )
    tree.SetBranchAddress( "ECal", observable_dict["ECal"] )
    tree.SetBranchAddress( "Recoq3", observable_dict["Recoq3"] )
    tree.SetBranchAddress( "pfl", observable_dict["pfl"] )
    tree.SetBranchAddress( "pfl_theta", observable_dict["pfl_theta"] )
    tree.SetBranchAddress( "pfl_phi", observable_dict["pfl_phi"] )
    tree.SetBranchAddress( "proton_mom", observable_dict["proton_mom"] )
    tree.SetBranchAddress( "proton_phi", observable_dict["proton_phi"] )
    tree.SetBranchAddress( "proton_theta", observable_dict["proton_theta"] )
    tree.SetBranchAddress( "HadAlphaT", observable_dict["HadAlphaT"] )
    tree.SetBranchAddress( "HadDeltaPT", observable_dict["HadDeltaPT"] )
    tree.SetBranchAddress( "HadDeltaPhiT", observable_dict["HadDeltaPhiT"] )

    i = 0
    while tree.GetEntry(i):
        content_x = 0
        content_y = 0
        w = observable_dict["TotWeight"]
        if observable_x[-4::1] == "_phi" or observable_x[-6::1] == "_theta":
            content_x = observable_dict[observable_x] * 180 / ROOT.TMath.Pi()
        else:
            content_x = observable_dict[observable_x]
        if observable_y[-4::1] == "_phi" or observable_y[-6::1] == "_theta":
            content_y = observable_dict[observable_y] * 180 / ROOT.TMath.Pi()
        else:
            content_y = observable_dict[observable_y]
        hist.Fill(observable_dict[observable_x], observable_dict[observable_y], w)
        i += 1
    return hist

def compute_acceptance_2d(file_names, output_name, observable_x, observable_y, nbins_x, min_val_x, max_val_x, nbins_y, min_val_y, max_val_y):

    dict_of_hists = {}
    hist_recoacc
    hist_trueacc

    for file_name in file_names:
        canvas = ROOT.TCanvas("canvas","canvas", 800, 600)
        canvas.cd()

        hist_recoacc = compute_event_rate( file_name+"_truereco.root", "Reco MC ACC", "MCCLAS6Tree", observable_x, observable_y, nbins_x, min_val_x, max_val_x, nbins_y, min_val_y, max_val_y )
        hist_trueacc = compute_event_rate( file_name+"_true.root", "True MC ACC", "MCCLAS6Tree", observable_x, observable_y, nbins_x, min_val_x, max_val_x, nbins_y, min_val_y, max_val_y )
        hists = [hist_recoacc, hist_trueacc]
        dict_of_hists[file_name] = hists

    output_name = output_name + "_acceptance_correction_" + observable_x + "_vs_" + observable_y
    output_file = ROOT.TFile((output_name + ".root"), "RECREATE")

    ratios_list = []

    for i in range(len(file_names)):
        ratio = dict_of_hists[file_names[i]][1].Clone()
        ratio.Divide(dict_of_hists[file_names[i]][0])
        ratio.SetName(f"Acceptance_Correction_Model_{i}")
        ratio.GetXaxis().SetTitle(observable_x)
        ratio.GetYaxis().SetTitle(observable_y)
        ratio.Draw("colz")
        canvas.SetName(f"Canvas_Acceptance_Correction_Model_{i}")
        canvas.Write()
        ratio.Write()
        ratios_list.append(ratio)

    n = len(ratios_list)
    average = ratios_list[0].Clone()
    for i in range(1,n):
        average.Add(ratios_list[0])
    average.Scale(1/float(n))
    ratio.SetName("Total_Acceptance_Correction")
    ratio.GetXaxis().SetTitle(observable_x)
    ratio.GetYaxis().SetTitle(observable_y)
    average.Draw("colz")
    canvas.Write()
    average.Write()

    return n

if __name__ == "__main__":
    file_name = ["/Users/juliatenavidal/Desktop/Postdoc/e4nu/AnalisedFiles/2024Generation/e4nuanalysis_1p0pi_G18_10a_Dipole_CFG_Q2_01_1GeV_eCarbon_NoRad"]
    output_name = "/Users/juliatenavidal/Desktop/Postdoc/e4nu/AnalisedFiles/clas6analysis_1p0pi_1GeV_"
    observable_x = "proton_mom"
    nbins_x = 30
    min_val_x = 0.6
    max_val_x = 2.3
    observable_y = "ECal"
    nbins_y = 30
    min_val_y = 0.6
    max_val_y = 2.3

    n = compute_acceptance_2d(file_name, output_name, observable_x, observable_y, nbins_x, min_val_x, max_val_x, nbins_y, min_val_y, max_val_y)
    print( n )
