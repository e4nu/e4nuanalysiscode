void MakeTheSmallPadPretty(TH1D* histo){

	histo->SetTitle("");
	histo->GetXaxis()->SetTitle("");
	histo->GetXaxis()->SetLabelSize(0.00);

	histo->GetYaxis()->SetTitle("(G-D) / D");
	histo->GetYaxis()->SetNdivisions(6);
	histo->GetYaxis()->SetLabelSize(0.15);
	histo->GetYaxis()->SetTitleSize(0.18);
	histo->GetYaxis()->SetTitleOffset(0.15);
};
