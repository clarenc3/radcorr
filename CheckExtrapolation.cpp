// Simple ROOT script for checking interpolation
// Run by root -l 'CheckExtrapolation.cpp'
void CheckExtrapolation() {
  TFile *f = new TFile("Elspectrum_muon_neutrino_merge_extend2.root");
  TGraph *g = (TGraph*)f->Get("Elspectrum_muon_neutrino_neutron_1GeV");

  // Now make a new Graph where we extrapolate using G
  TH1D *h = new TH1D("h", "h", 100, 0, 10);
  TH1D *hs = new TH1D("hs", "hs", 100, 0, 10);
  for (int i = 0; i < h->GetXaxis()->GetNbins(); ++i) {
    h->SetBinContent(i+1, g->Eval(h->GetXaxis()->GetBinCenter(i+1), 0, ""));
    hs->SetBinContent(i+1, g->Eval(h->GetXaxis()->GetBinCenter(i+1), 0, "S"));
  }

  TCanvas *canv = new TCanvas("canv", "canv", 1024, 1024);
  h->Draw();
  h->GetXaxis()->SetTitle("Q^{2} (GeV^{2}/c^{4})");
  h->GetYaxis()->SetTitle("d#sigma_{Radcorr}/d#sigma_{Tree}");

  hs->SetLineColor(kRed);
  hs->Draw("same");
  g->SetLineColor(kBlue);
  g->Draw("same");

  h->SetTitle("Linear interp.");
  hs->SetTitle("Spline interp.");
  g->SetTitle("Graph");
  g->SetLineWidth(2);
  g->SetLineColor(kGreen);

  g->SetMarkerSize(0);
  h->SetMarkerSize(0);
  hs->SetMarkerSize(0);

  g->SetFillStyle(0);
  h->SetFillStyle(0);
  hs->SetFillStyle(0);

  TLegend *leg = canv->BuildLegend(0.5, 0.5, 0.9, 0.9);
  leg->SetHeader("1 GeV, #nu_{#mu}");
}
