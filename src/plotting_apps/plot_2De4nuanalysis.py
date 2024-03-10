import ROOT
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
import sys
import time
import argparse
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
observable_dict["DataNormalization"] = c_double(0)
observable_dict["MCNormalization"] = c_double(0)


def compute_event_rate( file_name, hist_name, tree_name, observable_x, observable_y, nbins_x, min_val_x, max_val_x, nbins_y, min_val_y, max_val_y ):
    # We want this function to be able to open any root file and loop over the events
    hist = ROOT.TH2D( hist_name, "", nbins_x, min_val_x, max_val_x, nbins_y, min_val_y, max_val_y)  # Need to be able to set title
    rootfile = ROOT.TFile.Open(file_name, "READ")
    if not rootfile:
        print("ERROR: the", file_name," does not exist.")
        return
    tree = rootfile.Get(tree_name)  # This will be different for data

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

    Normalization_Factor = 0
    if tree_name == "CLAS6Tree":
        tree.SetBranchAddress( "DataNormalization", observable_dict["DataNormalization"] )
        Normalization_Factor = observable_dict["DataNormalization"]
    elif tree_name == "MCCLAS6Tree":
        tree.SetBranchAddress( "MCNormalization", observable_dict["MCNormalization"] )
        Normalization_Factor = observable_dict["MCNormalization"]

    i = 0
    while tree.GetEntry(i):
        content_x = 0
        content_y = 0
        w = observable_dict["TotWeight"]
        if observable_x[-4::1] == "_phi" or observable_x[-6::1] == "_theta":
            content_x = observable_dict[observable_x] #* 180 / ROOT.TMath.Pi()
        else:
            content_x = observable_dict[observable_x]
        if observable_y[-4::1] == "_phi" or observable_y[-6::1] == "_theta":
            content_y = observable_dict[observable_y] #* 180 / ROOT.TMath.Pi()
        else:
            content_y = observable_dict[observable_y]
        hist.Fill(content_x, content_y, w)  
        # hist.Fill(observable_dict[observable_x], observable_dict[observable_y], w)
        i += 1
    
    rootfile.Close()
    return hist, Normalization_Factor

def update_dict_of_hists(dict, file_name, hist_recoacc, hist_trueacc):

    hists = [hist_recoacc, hist_trueacc]
    dict[file_name] = hists

def format_hist(hist, title, observable_x, observable_y):
    hist.SetName(title)
    hist.SetTitle(title)
    hist.GetXaxis().SetTitle(observable_x)
    hist.GetYaxis().SetTitle(observable_y)
    

def compute_acceptance_2d(dictionary, output_name, observable_x, observable_y):

    average = list(dictionary.values())[0][0].Clone()
    average.Reset()
    
    output_name = output_name + "_acceptance_correction_" + observable_x + "_vs_" + observable_y
    output_file = ROOT.TFile((output_name + ".root"), "RECREATE")

    ratios_list = []

    i = 0
    for value in list(dictionary.values()):
        canvas = ROOT.TCanvas("canvas","canvas", 800, 600)
        canvas.cd()
        format_hist(value[0],f"True_Reco_Histogram_{i}", observable_x, observable_y)
        value[0].Draw("colz")
        canvas.SetName(f"Canvas_True_Reco_Histogram_{i}")
        canvas.Write()
        value[0].Write()

        format_hist(value[1],f"True_Histogram_{i}", observable_x, observable_y)
        value[1].Draw("colz")
        canvas.SetName(f"Canvas_True_Histogram_{i}")
        canvas.Write()
        value[1].Write()

        ratio = value[1].Clone()
        ratio.Divide(value[0])
        format_hist(ratio,f"Acceptance_Correction_Model_{i}", observable_x, observable_y)
        ratio.Draw("colz")
        canvas.SetName(f"Canvas_Acceptance_Correction_Model_{i}")
        canvas.Write()
        ratio.Write()
        ratios_list.append(ratio)
        i = i + 1

    n = len(ratios_list)
    for i in range(n):
        average.Add(ratios_list[0])
    average.Scale(1/float(n))
    ratio.SetName("Total_Acceptance_Correction")
    ratio.GetXaxis().SetTitle(observable_x)
    ratio.GetYaxis().SetTitle(observable_y)
    average.Draw("colz")
    canvas.Write()
    average.Write()

    output_file.Close()

    return average

def multiply(hist1, hist2):
    hist = hist1.Clone()
    hist.Multiply(hist2)
    return hist

def projection_hist(hist, projection_name, observable_x, start_y_bin = 1, end_y_bin = -1):
    hist_projection = hist.ProjectionX(projection_name, start_y_bin, end_y_bin)
    format_hist(hist_projection, projection_name, observable_x, "Event Rate")
    hist_projection.GetYaxis().SetTitle("Cross Section")
    hist_projection.Draw("colz")
    canvas.SetName("Canvas_" + projection_name)
    canvas.Write()
    hist_projection.Write()
    return hist_projection

def draw_print_two_histos(hist1, hist2, canvas_name, file):
    hist1.Draw("hist colz")
    hist2.Draw("SAME")
    canvas.SetName(canvas_name)
    canvas.Write()
    canvas.Print(file)


if __name__ == "__main__":
    # file_data = "C:/Users/ks202/Desktop/University/2024A/Physics_Laboratory_C/e4nuanalysis_1p0pianalysis_e_on_1000060120_2261MeV_3MaxBkgMult_clas6data.root"
    # file_names = ["C:/Users/ks202/Desktop/University/2024A/Physics_Laboratory_C/e4nuanalysis_1p0pi_GEM21_11a_Dipole_CFG_Q2_04_2GeV_eCarbon_NoRad"]
    # output_name = "C:/Users/ks202/Desktop/University/2024A/Physics_Laboratory_C/clas6analysis_1p0pi_1GeV_"
    # observable_x = "proton_mom"
    # nbins_x = 30
    # min_val_x = 0.25
    # max_val_x = 2
    # observable_y = "proton_theta"
    # nbins_y = 30
    # min_val_y = 0
    # max_val_y = 180

    # --data_file "C:/Users/ks202/Desktop/University/2024A/Physics_Laboratory_C/e4nuanalysis_1p0pianalysis_e_on_1000060120_2261MeV_3MaxBkgMult_clas6data.root"
    # --files_names "C:/Users/ks202/Desktop/University/2024A/Physics_Laboratory_C/e4nuanalysis_1p0pi_GEM21_11a_Dipole_CFG_Q2_04_2GeV_eCarbon_NoRad" 
    # --output_name "C:/Users/ks202/Desktop/University/2024A/Physics_Laboratory_C/clas6analysis_1p0pi_1GeV_"
    # --print_file "C:/Users/ks202/Desktop/University/2024A/Physics_Laboratory_C/Canvas_Cross_Section_Comparison.pdf"
    # --observable_x "proton_mom"
    # --nbins_x 30
    # --min_val_x 0.6
    # --max_val_x 2.3
    # --observable_y "ECal"
    # --nbins_y 30
    # --min_val_y 0.6
    # --max_val_y 2.3

    parser = argparse.ArgumentParser(usage=__doc__)
    #files
    parser.add_argument('--data_file', metavar = 'data_file', type = str, help = 'Please enter the file with the experimental data')
    parser.add_argument('--files_names', metavar = 'file_names', type = list, help = 'Please enter the list of Monte-Carlo model files')
    parser.add_argument('--output_name', metavar = 'output_name', type = str, help = 'Please enter the file directory into which you want to write the histograms')
    parser.add_argument('--print_file', metavar = 'print_file', type = str, help = 'Please enter the file into which you want to print the final histograms')
    #observable_x
    parser.add_argument('--observable_x', metavar = 'observable_x', type = str, help = 'Please enter the first observable that you want to use')
    parser.add_argument('--nbins_x', metavar = 'nbins_x', type = int, help = 'Please enter the number of bins for the first observable')
    parser.add_argument('--min_val_x', metavar = 'min_val_x', type = float, help = 'Please enter the maximum value of the first observable')
    parser.add_argument('--max_val_x', metavar = 'max_val_x', type = float, help = 'Please enter the minimum value of the first observable')
    #observable_y
    parser.add_argument('--observable_y', metavar = 'observable_y', type = str, help = 'Please enter the second observable that you want to use')
    parser.add_argument('--nbins_y', metavar = 'nbins_y', type = int, help = 'Please enter the number of bins for the second observable')
    parser.add_argument('--min_val_y', metavar = 'min_val_y', type = float, help = 'Please enter the maximum value of the second observable')
    parser.add_argument('--max_val_y', metavar = 'max_val_y', type = float, help = 'Please enter the minimum value of the second observable')

    args = parser.parse_args()

    file_data = args.data_file
    file_names = args.files_names
    output_name = args.output_name
    observable_x = args.observable_x
    nbins_x = args.nbins_x
    min_val_x = args.min_val_x
    max_val_x = args.max_val_x
    observable_y = args.observable_y
    nbins_y = args.nbins_y
    min_val_y = args.min_val_y
    max_val_y = args.max_val_y

    dictionary = {}

    for file_name in file_names:
        hist_recoacc, Norm_Reco = compute_event_rate( file_name+"_truereco.root", "Reco MC ACC", "MCCLAS6Tree", observable_x, observable_y, nbins_x, min_val_x, max_val_x, nbins_y, min_val_y, max_val_y )
        hist_recoacc.Scale(Norm_Reco, option = "width")
        hist_trueacc, Norm_True = compute_event_rate( file_name+"_true.root", "True MC ACC", "MCCLAS6Tree", observable_x, observable_y, nbins_x, min_val_x, max_val_x, nbins_y, min_val_y, max_val_y )
        hist_trueacc.Scale(Norm_True, option = "width")
        
        update_dict_of_hists(dictionary, file_name, hist_recoacc, hist_trueacc)

    acceptance = compute_acceptance_2d(dictionary, output_name, observable_x, observable_y)

    canvas = ROOT.TCanvas("canvas","canvas", 800, 600)
    canvas.cd()

    hist_data, Norm_Data = compute_event_rate( file_data, "CLAS6Data", "CLAS6Tree", observable_x, observable_y, nbins_x, min_val_x, max_val_x, nbins_y, min_val_y, max_val_y )
    hist_data.Scale(Norm_Data, option = "width")
    final_output_name = output_name + "_Acceptance_Corrected_" + observable_x + "_vs_" + observable_y
    final_output_file = ROOT.TFile((final_output_name + ".root"), "RECREATE")
    final_output_file.cd()

    format_hist(hist_data, "Data_Histogram", observable_x, observable_y)
    hist_data.Draw("colz")
    canvas.SetName(f"Data_Histogram")
    canvas.Write()
    hist_data.Write()

    corrected_hist = multiply(acceptance, hist_data)
    format_hist(corrected_hist, "Acceptance_Corrected_Histogram", observable_x, observable_y)
    corrected_hist.Draw("colz")
    canvas.SetName(f"Canvas_Acceptance_Corrected_Histogram")
    canvas.Write()
    corrected_hist.Write()

    values = list(dictionary.values())
    keys = list(dictionary.keys())

    start_y_bin = 10
    end_y_bin = 20

    hist_true_projections = {}
    for key in keys:
        i = keys.index(key)
        hist_true_mc_projection = projection_hist(dictionary[key][1], f"True_MC_Histogram_Projection_Model_{i}", observable_x, start_y_bin, end_y_bin)
        total_hist_true_mc_projection = projection_hist(dictionary[key][1], f"True_MC_Histogram_Single_Differential_Model_{i}", observable_x)
        hist_true_projections[key] = [hist_true_mc_projection, total_hist_true_mc_projection]

    # hist_data_projection = projection_hist(hist_data, "Data_Histogram_Projection", observable_x, start_y_bin, end_y_bin)

    corrected_hist_projection = projection_hist(corrected_hist, "Acceptance_Corrected_Histogram_Projection_Double_Differential", observable_x, start_y_bin, end_y_bin)

    total_corrected_hist_projection = projection_hist(corrected_hist, "Acceptance_Corrected_Histogram_Single_Differential", observable_x)

    # file = f"C:/Users/ks202/Desktop/University/2024A/Physics_Laboratory_C/Canvas_Cross_Section_Comparison_{observable_x}_vs_{observable_y}.pdf"
    file = args.print_file
    canvas.Print(file + "[")
    for key in keys:
        i = keys.index(key)
        draw_print_two_histos(hist_true_projections[key][0], corrected_hist_projection, f"Canvas_Second_Differential_Cross_Section_Comparison_Model_{i}", file)
        draw_print_two_histos(hist_true_projections[key][1], total_corrected_hist_projection, f"Canvas_Differential_Cross_Section_Comparison_Model_{i}", file)
    canvas.Print(file + "]")

    