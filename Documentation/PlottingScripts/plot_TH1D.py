import os, optparse
import ROOT

def plot_TH1D_as_TGraph(root_file_names, output_file_name, histogram_names, legend_names, x_axis_labels, y_axis_labels, convert_graph):
    # Create a canvas
    canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)

    # Add a legend
    legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)

    files = []
    histograms = []
    graphs = []

    if len(root_file_names) == len(histogram_names) :     
        for i in range(len(root_file_names)) :

            x_axis_label = ''
            if len(x_axis_labels) == len(root_file_names) :
                x_axis_label = x_axis_labels[i]
            y_axis_label = ''
            if len(y_axis_labels) == len(root_file_names) :
                y_axis_label = y_axis_labels[i]
    
            # Open the root file
            #print(root_file_names[i])
            files.append(ROOT.TFile(root_file_names[i]))
    
            # Get the histogram
            histograms.append(files[i].Get(histogram_names[i]))
            histograms[i].SetLineColor(i+1)
            histograms[i].SetLineStyle(i+1)            
            histograms[i].GetXaxis().SetTitle(x_axis_label)
            histograms[i].GetXaxis().SetTitle(y_axis_label)
            # Set the histogram color and style
            # Convert the histogram to a TGraph                                                                                      
            if convert_graph: 
                graphs.append(ROOT.TGraph(histograms[i]))
                graphs[i].SetLineColor(i+1)
                graphs[i].SetLineStyle(i+1)
                graphs[i].GetXaxis().SetTitle(x_axis_label)
                graphs[i].GetXaxis().SetTitle(y_axis_label)
            
            if i == 0:
                if convert_graph : 
                    graphs[i].Draw("AC")
                else : 
                    histograms[i].Draw("hist")
            else:
                if convert_graph :
                    graphs[i].Draw("C")
                else :
                    histograms[i].Draw("hist same")

    for i, legend_name in enumerate(legend_names):
        legend.AddEntry(histograms[i], legend_name, "l")

    if len(legend_names) != 0 :
        legend.Draw()

    # Save the plot as a pdf file
    canvas.SaveAs(output_file_name+".pdf")
        

op = optparse.OptionParser(usage=__doc__)
op.add_option("--files", dest="files", help="ROOT files from which to read the histogram")
op.add_option("--output-name", dest="output_name", default="e4nu_analysis_comparison",help="Output name for the pdf file. Default %default")
op.add_option("--hist_names", dest="hists", help="Name of TH1Ds to plot")
op.add_option("--legend-list", dest="legends", default="", help="Name for the legend. Default none.")
op.add_option("--x-labels", dest="xlabels", default="", help="X Axis labels")
op.add_option("--y-labels", dest="ylabels", default="", help="Y Axis labels")
op.add_option("--plot-as-graph", dest="graph", default=False, action="store_true", help="Store histograms as graphs")
opts, args = op.parse_args()

req_file_list = (opts.files).split(',')
req_hist_list = (opts.hists).split(',')
req_legends = []
if opts.legends : 
    req_legends = (opts.legends).split(',')
req_xlabels = []
if opts.xlabels :
    req_xlabels = (opts.xlabels).split(',')
req_ylabels = []
if opts.ylabels :
    req_ylabels = (opts.ylabels).split(',')

plot_TH1D_as_TGraph(req_file_list,opts.output_name,req_hist_list,req_legends,req_xlabels,req_ylabels,opts.graph)
