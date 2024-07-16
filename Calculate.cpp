#include <iostream>

#include "RadCorrCalculator.h"

// ROOT libraries for printing
#include "TH2D.h"
#include "TCanvas.h"
#include "TString.h"
#include "TStyle.h"

// Example of how the radiative corrections calculator can be used
// Here filling a histogram in four-momentum transfer (Q^2) and neutrino energy (E_nu) with the weight
// Using ROOT, printing to a multi-page PDF for muon neutrinos and muon anti-neutrinos

int main() {

  gStyle->SetOptStat(0);
  gStyle->SetNumberContours(255);

  bool linear = true; // Linear or TSpline3 interpolation. Use linear

  RadCorrCalc calc; // The calculator

  int nbinsenu = 300;
  double minenu = 0;
  double maxenu = 3;
  int nbinsq2 = 300;
  double minq2 = 0;
  double maxq2 = 3;
  double dEnu = (maxenu-minenu)/nbinsenu;
  double dQ2 = (maxq2-minq2)/nbinsq2;

  // Make the canvas
  TCanvas *canv = new TCanvas("canv", "canv", 1024, 1024);
  TString canvname = "calc_radcorr_new1pi_andextrap";
  if (linear) {
    calc.SetLinearInterp();
    canvname+="_linear";
  } else {
    calc.SetSplineInterp();
    canvname+="_spline";
  }
  //canvname+=Form("_test_final.pdf");
  canvname+=Form("_nueinc_newstyle.pdf");
  canv->Print(canvname+"[");
  // Make the TH2D
  TH2D *plot[4];

  // Loop over numu and numubar
  for (int type = 0; type < kNuEBar+1; ++type) {
    calc.SetNuType(NuType(type)); // Set the neutrino type, 0 = numu, 1 = numubar

    TString typestring;
    if      (type == 0) typestring = "#nu_{#mu}";
    else if (type == 1) typestring = "#bar{#nu}_{#mu}";
    else if (type == 2) typestring = "#nu_{e}";
    else if (type == 3) typestring = "#bar{#nu}_{e}";

    TString namestring = typestring;
    namestring.ReplaceAll("#","");
    namestring.ReplaceAll("{","");
    namestring.ReplaceAll("}","");

    plot[type] = new TH2D(Form("radcorr_weights_%s", namestring.Data()), Form("radcorr_weights_%s;E_{#nu} [GeV];Q^{2} [GeV^{2}]", typestring.Data()), nbinsenu, minenu, maxenu, nbinsq2, minq2, maxq2);

    for (int i = 0; i < nbinsenu; ++i) {
      double Enu = dEnu*(i+1);
      for (int j = 0; j < nbinsq2; ++j) {
        double Q2 = dQ2*(j+1);
        double weight = calc.CalcWeight(Enu, Q2);
        plot[type]->SetBinContent(i+1, j+1, weight);
      }
    }

    canv->cd();
    plot[type]->SetMaximum(1.5);
    plot[type]->SetMinimum(0.5);
    plot[type]->Draw("colz");
    canv->Print(canvname);
  }

  TString rootname = canvname;
  rootname.ReplaceAll(".pdf", ".root");
  TFile *output = new TFile(rootname, "recreate");
  plot[0]->Write("radcorr_weights_numu");
  plot[1]->Write("radcorr_weights_numubar");
  plot[2]->Write("radcorr_weights_nue");
  plot[3]->Write("radcorr_weights_nuebar");
  output->Close();

  canv->Print(canvname+"]");
}
