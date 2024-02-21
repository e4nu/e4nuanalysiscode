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

def compute_acceptance_2d(file_names, output_name, observable_x, observable_y, nbins_x, min_val_x, max_val_x, nbins_y, min_val_y, max_val_y):

    dict_of_hists = {}

    for file_name in file_names:
        print (file_name)
        rootfile_mcrecoacc = ROOT.TFile.Open(file_name+"_truereco.root", "READ")
        rootfile_mctrueacc = ROOT.TFile.Open(file_name+"_true.root", "READ")
        canvas = ROOT.TCanvas("canvas","canvas", 800, 600)
        canvas.cd()

        if not rootfile_mcrecoacc:
            print("ERROR: the", file_name, "_truereco.root does not exist.")
            return
        if not rootfile_mctrueacc:
            print("ERROR: the", file_name, "_true.root does not exist.")
            return

        tree_mcrecoacc = rootfile_mcrecoacc.Get("MCCLAS6Tree")
        tree_mctrueacc = rootfile_mctrueacc.Get("MCCLAS6Tree")

        hist_recoacc = ROOT.TH2D("Reco MC ACC", "", nbins_x, min_val_x, max_val_x, nbins_y, min_val_y, max_val_y)
        hist_trueacc = ROOT.TH2D("True MC ACC", "", nbins_x, min_val_x, max_val_x, nbins_y, min_val_y, max_val_y)

        trees = [tree_mcrecoacc, tree_mctrueacc]
        hists = [hist_recoacc, hist_trueacc]

        for j in tqdm(range(len(trees))):
            content_x = 0
            content_y = 0

            trees[j].SetBranchAddress( "TotWeight", observable_dict["TotWeight"] )
            #observable_dict["IsBkg"] = trees[j].SetBranchAddress( "IsBkg", IsBkg )
            trees[j].SetBranchAddress( "RecoW", observable_dict["RecoW"] )
            trees[j].SetBranchAddress( "ECal", observable_dict["ECal"] )
            trees[j].SetBranchAddress( "Recoq3", observable_dict["Recoq3"] )
            trees[j].SetBranchAddress( "pfl", observable_dict["pfl"] )
            trees[j].SetBranchAddress( "pfl_theta", observable_dict["pfl_theta"] )
            trees[j].SetBranchAddress( "pfl_phi", observable_dict["pfl_phi"] )
            trees[j].SetBranchAddress( "proton_mom", observable_dict["proton_mom"] )
            trees[j].SetBranchAddress( "proton_phi", observable_dict["proton_phi"] )
            trees[j].SetBranchAddress( "proton_theta", observable_dict["proton_theta"] )
            trees[j].SetBranchAddress( "HadAlphaT", observable_dict["HadAlphaT"] )
            trees[j].SetBranchAddress( "HadDeltaPT", observable_dict["HadDeltaPT"] )
            trees[j].SetBranchAddress( "HadDeltaPhiT", observable_dict["HadDeltaPhiT"] )
            #trees[j].SetBranchAddress( "ElectronSector", observable_dict["ElectronSector"] )

            i = 0
            while trees[j].GetEntry(i):
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

                hists[j].Fill(observable_dict[observable_x], observable_dict[observable_y], w)
                print( observable_dict[observable_x], observable_dict[observable_y], w)
                i += 1
                # print(w, content_x, content_y)
                # hists[i].SetLineWidth(3)

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
    #average.Scale(1/float(n))
    ratio.SetName("Total_Acceptance_Correction")
    ratio.GetXaxis().SetTitle(observable_x)
    ratio.GetYaxis().SetTitle(observable_y)
    average.Draw("colz")
    canvas.Write()
    average.Write()

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

    compute_acceptance_2d(file_name, output_name, observable_x, observable_y, nbins_x, min_val_x, max_val_x, nbins_y, min_val_y, max_val_y)
