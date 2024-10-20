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
observable_dict["E_{Cal}"] = c_double(0)
observable_dict["Efl"] = c_double(0)
observable_dict["Recoq3"]= c_double(0)
observable_dict["pfl"] = c_double(0)
observable_dict["pfl_{#theta}"] = c_double(0)
observable_dict["pfl_{#phi}"] = c_double(0)
observable_dict["p_{mom}"] = c_double(0)
observable_dict["p_{#phi}"] = c_double(0)
observable_dict["p_{#theta}"] = c_double(0)
observable_dict["Had#alpha_{T}"] = c_double(0)
observable_dict["P_{T}"] = c_double(0)
observable_dict["Had#delta#phi_{T}"] = c_double(0)
observable_dict["ElectronSector"] = int(-99)
observable_dict["DataNormalization"] = c_double(0)
observable_dict["MCNormalization"] = c_double(0)
observable_dict["QEL"] = c_bool(False)
observable_dict["RES"] = c_bool(False)
observable_dict["DIS"] = c_bool(False)
observable_dict["MEC"] = c_bool(False)
observable_dict["all"] = c_bool(True)

unit_dict = {}

unit_dict["E_{Cal}"] = "GeV"
unit_dict["Efl"] = "GeV"
unit_dict["pfl"] = "GeV/c"
unit_dict["pfl_{#theta}"] = "deg"
unit_dict["pfl_{#phi}"] = "deg"
unit_dict["p_{mom}"] = "GeV/c"
unit_dict["p_{#phi}"] = "deg"
unit_dict["p_{#theta}"] = "deg"
unit_dict["Had#alpha_{T}"] = "deg"
unit_dict["P_{T}"] = "GeV/c"
unit_dict["Had#delta#phi_{T}"] = "deg"


def compute_event_rate( file_name, hist_name, tree_name, observable_x, observable_y, nbins_x, min_val_x, max_val_x, nbins_y, min_val_y, max_val_y, hist_type = "all"):
    # We want this function to be able to open any root file and loop over the events
    hist = ROOT.TH2D( hist_name, "", nbins_x, min_val_x, max_val_x, nbins_y, min_val_y, max_val_y)  # Need to be able to set title
    rootfile = ROOT.TFile.Open(file_name, "READ")
    if not rootfile:
        print("ERROR: the", file_name," does not exist.")
        return
    tree = rootfile.Get(tree_name)  # This will be different for data

    tree.SetBranchAddress( "TotWeight", observable_dict["TotWeight"] )
    tree.SetBranchAddress( "RecoW", observable_dict["RecoW"] )
    tree.SetBranchAddress( "ECal", observable_dict["E_{Cal}"] )
    tree.SetBranchAddress( "Efl", observable_dict["Efl"] )
    tree.SetBranchAddress( "Recoq3", observable_dict["Recoq3"] )
    tree.SetBranchAddress( "pfl", observable_dict["pfl"] )
    tree.SetBranchAddress( "pfl_theta", observable_dict["pfl_{#theta}"] )
    tree.SetBranchAddress( "pfl_phi", observable_dict["pfl_{#phi}"] )
    tree.SetBranchAddress( "proton_mom", observable_dict["p_{mom}"] )
    tree.SetBranchAddress( "proton_phi", observable_dict["p_{#phi}"] )
    tree.SetBranchAddress( "proton_theta", observable_dict["p_{#theta}"] )
    tree.SetBranchAddress( "HadAlphaT", observable_dict["Had#alpha_{T}"] )
    tree.SetBranchAddress( "HadDeltaPT", observable_dict["P_{T}"] )
    tree.SetBranchAddress( "HadDeltaPhiT", observable_dict["Had#delta#phi_{T}"] )

    Normalization_Factor = 0
    if tree_name == "CLAS6Tree":
        tree.SetBranchAddress( "DataNormalization", observable_dict["DataNormalization"] )
        Normalization_Factor = observable_dict["DataNormalization"]
    elif tree_name == "MCCLAS6Tree":
        tree.SetBranchAddress( "MCNormalization", observable_dict["MCNormalization"] )
        Normalization_Factor = observable_dict["MCNormalization"]
        tree.SetBranchAddress( "QEL", observable_dict["QEL"])
        tree.SetBranchAddress( "RES", observable_dict["RES"])
        tree.SetBranchAddress( "DIS", observable_dict["DIS"])
        tree.SetBranchAddress( "MEC", observable_dict["MEC"])

    i = 0
    while tree.GetEntry(i):
        content_x = 0
        content_y = 0
        if observable_dict[hist_type]:
            w = observable_dict["TotWeight"]
            if observable_x[-7::1] == "_{#phi}" or observable_x[-9::1] == "_{#theta}":
                content_x = observable_dict[observable_x] #* 180 / ROOT.TMath.Pi()
            else:
                content_x = observable_dict[observable_x]
            if observable_y[-7::1] == "_{#phi}" or observable_y[-9::1] == "_{#theta}":
                content_y = observable_dict[observable_y] #* 180 / ROOT.TMath.Pi()
            else:
                content_y = observable_dict[observable_y]
            hist.Fill(content_x, content_y, w)
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

def cross_section(observable):
    cs = ""
    if observable[:6:1] == "#delta":
        cs = "#frac{d#sigma}{d%s}[#frac{#mub}{%s}]" %(observable[6::1], unit_dict[observable])
    else:
        cs = "#frac{d#sigma}{d%s}[#frac{#mub}{%s}]" %(observable, unit_dict[observable])
    return cs

def double_cross_section(observable_x, observable_y):
    cs = ""
    observable_new_x = observable_x
    observable_new_y = observable_y
    if observable_x[:6:1] == "#delta":
        observable_new_x = observable_x[:6:1]
    if observable_y[:6:1] == "#delta":
        observable_new_y = observable_y[:6:1]
    if unit_dict[observable_x] == unit_dict[observable_y]:
        cs = "#frac{d^{2}#sigma}{d%sd%s}[#frac{#mub}{%s^{2}}]" %(observable_new_x, observable_new_y, unit_dict[observable_x])
    else:
        cs = "#frac{d^{2}#sigma}{d%sd%s}[#frac{#mub}{%s*%s}]" %(observable_new_x, observable_new_y, unit_dict[observable_x], unit_dict[observable_y])
    return cs

def projection_hist(hist, projection_name, observable_x, observable_y, start_y_bin = 1, end_y_bin = -1):
    hist_projection = hist.ProjectionX(projection_name, start_y_bin, end_y_bin)
    if start_y_bin - end_y_bin == -1:
        format_hist(hist_projection, projection_name, f"{observable_x}[{unit_dict[observable_x]}]", double_cross_section(observable_x, observable_y))
    else:
        format_hist(hist_projection, projection_name, f"{observable_x}[{unit_dict[observable_x]}]", cross_section(observable_x))
    hist_projection.Draw("colz")
    canvas.SetName("Canvas_" + projection_name)
    canvas.Write()
    hist_projection.Write()
    return hist_projection

def draw_print_histos(hists, canvas_name, file):
    list_of_colors = [ROOT.kBlack, ROOT.kBlue, ROOT.kGreen, ROOT.kOrange, ROOT.kRed, ROOT.kBlack]
    legend_list = ["Total", "QEL", "RES", "DIS", "MEC", "Data"]
    hists[0].SetLineColor(list_of_colors[0])
    hists[0].SetMarkerStyle(8)
    hists[0].SetStats(0)
    hists[0].SetTitle(canvas_name)
    hists[0].Draw("hist colz C")
    for i in range(1,len(hists)-1):
        hists[i].SetLineColor(list_of_colors[i])
        hists[i].SetMarkerStyle(8)
        hists[i].SetStats(0)
        hists[i].Draw("hist SAME C")
    hists[-1].SetLineColor(list_of_colors[-1])
    hists[-1].Draw("SAME")
    # Top left corner:
    # leg = ROOT.TLegend(0.1,0.7,0.48,0.9)
    # Top right corner:
    leg = ROOT.TLegend(0.9,0.7,0.48,0.9)
    leg.SetBorderSize(0)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.03)
    for i in range(len(hists)-1):
        leg.AddEntry(hists[i], legend_list[i],"l")
    leg.AddEntry(hists[-1], legend_list[-1],"lep")
    leg.Draw()
    canvas.SetName(canvas_name)
    canvas.Write()
    canvas.Print(file)


if __name__ == "__main__":

    start_time = time.time()

    file_data = "C:/Users/ks202/Desktop/University/2024A/Physics_Laboratory_C/e4nuanalysis_1p0pianalysis_e_on_1000060120_2261MeV_3MaxBkgMult_clas6data.root"
    file_names = ["C:/Users/ks202/Desktop/University/2024A/Physics_Laboratory_C/e4nuanalysis_1p0pi_GEM21_11a_Dipole_CFG_Q2_04_2GeV_eCarbon_NoRad"]
    output_name = "C:/Users/ks202/Desktop/University/2024A/Physics_Laboratory_C/clas6analysis_1p0pi_1GeV_"
    observable_x = "P_{T}"
    nbins_x = 100
    min_val_x = 0.0
    max_val_x = 1.0
    observable_y = "E_{Cal}"
    nbins_y = 100
    min_val_y = 0.0
    max_val_y = 2.5
    start_y_bin = 90
    end_y_bin = 91
    start_projected_value_y = round(min_val_y + start_y_bin * (max_val_y - min_val_y)/nbins_y, 2)
    end_projected_value_y = round(min_val_y + end_y_bin * (max_val_y - min_val_y)/nbins_y, 2)

    # --data_file "C:/Users/ks202/Desktop/University/2024A/Physics_Laboratory_C/e4nuanalysis_1p0pianalysis_e_on_1000060120_2261MeV_3MaxBkgMult_clas6data.root"
    # --files_names "C:/Users/ks202/Desktop/University/2024A/Physics_Laboratory_C/e4nuanalysis_1p0pi_GEM21_11a_Dipole_CFG_Q2_04_2GeV_eCarbon_NoRad" 
    # --output_name "C:/Users/ks202/Desktop/University/2024A/Physics_Laboratory_C/clas6analysis_1p0pi_1GeV_"
    # --print_file "C:/Users/ks202/Desktop/University/2024A/Physics_Laboratory_C/Canvas_Cross_Section_Comparison.pdf"
    # --observable_x "p_{mom}"
    # --nbins_x 30
    # --min_val_x 0.6
    # --max_val_x 2.3
    # --observable_y "E_{Cal}"
    # --nbins_y 30
    # --min_val_y 0.6
    # --max_val_y 2.3

    # parser = argparse.ArgumentParser(usage=__doc__)
    # #files
    # parser.add_argument('--data_file', metavar = 'data_file', type = str, help = 'Please enter the file with the experimental data')
    # parser.add_argument('--files_names', metavar = 'file_names', type = list, help = 'Please enter the list of Monte-Carlo model files')
    # parser.add_argument('--output_name', metavar = 'output_name', type = str, help = 'Please enter the file directory into which you want to write the histograms')
    # parser.add_argument('--print_file', metavar = 'print_file', type = str, help = 'Please enter the file into which you want to print the final histograms')
    # #observable_x
    # parser.add_argument('--observable_x', metavar = 'observable_x', type = str, help = 'Please enter the first observable that you want to use')
    # parser.add_argument('--nbins_x', metavar = 'nbins_x', type = int, help = 'Please enter the number of bins for the first observable')
    # parser.add_argument('--min_val_x', metavar = 'min_val_x', type = float, help = 'Please enter the maximum value of the first observable')
    # parser.add_argument('--max_val_x', metavar = 'max_val_x', type = float, help = 'Please enter the minimum value of the first observable')
    # #observable_y
    # parser.add_argument('--observable_y', metavar = 'observable_y', type = str, help = 'Please enter the second observable that you want to use')
    # parser.add_argument('--nbins_y', metavar = 'nbins_y', type = int, help = 'Please enter the number of bins for the second observable')
    # parser.add_argument('--min_val_y', metavar = 'min_val_y', type = float, help = 'Please enter the maximum value of the second observable')
    # parser.add_argument('--max_val_y', metavar = 'max_val_y', type = float, help = 'Please enter the minimum value of the second observable')

    # args = parser.parse_args()

    # file_data = args.data_file
    # file_names = args.files_names
    # output_name = args.output_name
    # observable_x = args.observable_x
    # nbins_x = args.nbins_x
    # min_val_x = args.min_val_x
    # max_val_x = args.max_val_x
    # observable_y = args.observable_y
    # nbins_y = args.nbins_y
    # min_val_y = args.min_val_y
    # max_val_y = args.max_val_y

    dictionary = {}
    dict_dictionary = {}

    for file_name in file_names:
        hist_recoacc, Norm_Reco = compute_event_rate( file_name+"_truereco.root", "Reco MC ACC", "MCCLAS6Tree", observable_x, observable_y, nbins_x, min_val_x, max_val_x, nbins_y, min_val_y, max_val_y )
        hist_recoacc.Scale(Norm_Reco)
        dict_hist_trueacc = {}
        hist_trueacc, Norm_True = compute_event_rate( file_name+"_true.root", "True MC ACC", "MCCLAS6Tree", observable_x, observable_y, nbins_x, min_val_x, max_val_x, nbins_y, min_val_y, max_val_y )
        hist_trueacc_qel = compute_event_rate( file_name+"_true.root", "True MC ACC", "MCCLAS6Tree", observable_x, observable_y, nbins_x, min_val_x, max_val_x, nbins_y, min_val_y, max_val_y, "QEL" )[0]
        hist_trueacc_res = compute_event_rate( file_name+"_true.root", "True MC ACC", "MCCLAS6Tree", observable_x, observable_y, nbins_x, min_val_x, max_val_x, nbins_y, min_val_y, max_val_y, "RES" )[0]
        hist_trueacc_dis = compute_event_rate( file_name+"_true.root", "True MC ACC", "MCCLAS6Tree", observable_x, observable_y, nbins_x, min_val_x, max_val_x, nbins_y, min_val_y, max_val_y, "DIS" )[0]
        hist_trueacc_mec = compute_event_rate( file_name+"_true.root", "True MC ACC", "MCCLAS6Tree", observable_x, observable_y, nbins_x, min_val_x, max_val_x, nbins_y, min_val_y, max_val_y, "MEC" )[0]
        dict_hist_trueacc['hist'] = hist_trueacc
        dict_hist_trueacc['QEL'] = hist_trueacc_qel
        dict_hist_trueacc['RES'] = hist_trueacc_res
        dict_hist_trueacc['DIS'] = hist_trueacc_dis
        dict_hist_trueacc['MEC'] = hist_trueacc_mec
        for key in list(dict_hist_trueacc.keys()):
            dict_hist_trueacc[key].Scale(Norm_True)
        
        update_dict_of_hists(dictionary, file_name, hist_recoacc, hist_trueacc)
        dict_dictionary[file_name] = dict_hist_trueacc

    acceptance = compute_acceptance_2d(dictionary, output_name, observable_x, observable_y)

    canvas = ROOT.TCanvas("canvas","canvas", 800, 600)
    canvas.cd()

    hist_data, Norm_Data = compute_event_rate( file_data, "CLAS6Data", "CLAS6Tree", observable_x, observable_y, nbins_x, min_val_x, max_val_x, nbins_y, min_val_y, max_val_y )
    hist_data.Scale(Norm_Data)
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

    # values = list(dict_dictionary.values())
    keys = list(dict_dictionary.keys())

    hist_true_projections = {}
    for key in keys:
        i = keys.index(key)
        hist_true_mc_projection = projection_hist(dict_dictionary[key]['hist'], f"True_MC_Histogram_Projection_Model_{i}({min_val_x} < {observable_x} < {max_val_x},  {start_projected_value_y} < {observable_y} < {end_projected_value_y})", observable_x, observable_y, start_y_bin, end_y_bin)
        total_hist_true_mc_projection = projection_hist(dict_dictionary[key]['hist'], f"True_MC_Histogram_Single_Differential_Model_{i}({min_val_x} < {observable_x} < {max_val_x})", observable_x, observable_y)
        hist_true_mc_projection_qel = projection_hist(dict_dictionary[key]['QEL'], f"True_MC_Histogram_Projection_Model_{i}_QEL({min_val_x} < {observable_x} < {max_val_x},  {start_projected_value_y} < {observable_y} < {end_projected_value_y})", observable_x, observable_y, start_y_bin, end_y_bin)
        total_hist_true_mc_projection_qel = projection_hist(dict_dictionary[key]['QEL'], f"True_MC_Histogram_Single_Differential_Model_{i}_QEL({min_val_x} < {observable_x} < {max_val_x})", observable_x, observable_y)
        hist_true_mc_projection_res = projection_hist(dict_dictionary[key]['RES'], f"True_MC_Histogram_Projection_Model_{i}_RES({min_val_x} < {observable_x} < {max_val_x},  {start_projected_value_y} < {observable_y} < {end_projected_value_y})", observable_x, observable_y, start_y_bin, end_y_bin)
        total_hist_true_mc_projection_res = projection_hist(dict_dictionary[key]['RES'], f"True_MC_Histogram_Single_Differential_Model_{i}_RES({min_val_x} < {observable_x} < {max_val_x})", observable_x, observable_y)
        hist_true_mc_projection_dis = projection_hist(dict_dictionary[key]['DIS'], f"True_MC_Histogram_Projection_Model_{i}_DIS({min_val_x} < {observable_x} < {max_val_x},  {start_projected_value_y} < {observable_y} < {end_projected_value_y})", observable_x, observable_y, start_y_bin, end_y_bin)
        total_hist_true_mc_projection_dis = projection_hist(dict_dictionary[key]['DIS'], f"True_MC_Histogram_Single_Differential_Model_{i}_DIS({min_val_x} < {observable_x} < {max_val_x})", observable_x, observable_y)
        hist_true_mc_projection_mec = projection_hist(dict_dictionary[key]['MEC'], f"True_MC_Histogram_Projection_Model_{i}_MEC({min_val_x} < {observable_x} < {max_val_x},  {start_projected_value_y} < {observable_y} < {end_projected_value_y})", observable_x, observable_y, start_y_bin, end_y_bin)
        total_hist_true_mc_projection_mec = projection_hist(dict_dictionary[key]['MEC'], f"True_MC_Histogram_Single_Differential_Model_{i}_MEC({min_val_x} < {observable_x} < {max_val_x})", observable_x, observable_y)
        hist_true_projections[key] = [[hist_true_mc_projection, hist_true_mc_projection_qel, hist_true_mc_projection_res, hist_true_mc_projection_dis, hist_true_mc_projection_mec], [total_hist_true_mc_projection, total_hist_true_mc_projection_qel, total_hist_true_mc_projection_res, total_hist_true_mc_projection_dis, total_hist_true_mc_projection_mec]]

    # hist_data_projection = projection_hist(hist_data, "Data_Histogram_Projection", observable_x, start_y_bin, end_y_bin)

    corrected_hist_projection = projection_hist(corrected_hist, "Acceptance_Corrected_Histogram_Projection_Double_Differential", observable_x, observable_y, start_y_bin, end_y_bin)

    total_corrected_hist_projection = projection_hist(corrected_hist, "Acceptance_Corrected_Histogram_Single_Differential", observable_x, observable_y)

    file = f"C:/Users/ks202/Desktop/University/2024A/Physics_Laboratory_C/Cross_Section_Comparison_{observable_x}_vs_{observable_y}_{start_projected_value_y}_{end_projected_value_y}.pdf"
    # file = args.print_file
    canvas.Print(file + "[")
    for key in keys:
        i = keys.index(key)
        hist_true_projection = hist_true_projections[key][0]
        hist_true_projection.append(corrected_hist_projection)
        total_hist_true_projection = hist_true_projections[key][1]
        total_hist_true_projection.append(total_corrected_hist_projection)
        for i in range(2):
            for hist in hist_true_projections[key][i]:
                hist.Scale(1, option = "width")
        draw_print_histos(hist_true_projection, f"Canvas_Second_Differential_Cross_Section_Comparison_Model_{i}({min_val_x} < {observable_x} < {max_val_x}, {start_projected_value_y} < {observable_y} < {end_projected_value_y})", file)
        draw_print_histos(total_hist_true_projection, f"Canvas_Differential_Cross_Section_Comparison_Model_{i}({min_val_x} < {observable_x} < {max_val_x})", file)
    canvas.Print(file + "]")

    time_0 = time.time() - start_time
    print("Program runtime: %s minutes" % (round(time_0/60,2)))
    
