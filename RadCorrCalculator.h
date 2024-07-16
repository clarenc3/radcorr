#ifndef _RADCORRCALC_H_SEEN_
#define _RADCORRCALC_H_SEEN_

#include <iostream>
#include "TGraph.h"
#include "TFile.h"

//#define DEBUG

// Just a handy enum
enum NuType {kNuMu = 0, kNuMuBar = 1, kNuE = 2, kNuEBar = 3};

class RadCorrCalc {
  public:
    // Constructor and destructor
    RadCorrCalc();

    void SetupNuMu();
    void SetupNuE();
    void SetupGraphs();

    // Just delete the graphs
    ~RadCorrCalc() {
      // Delete the numu(bar) graphs
      for (int i = 0; i < kNuMuBar+1; ++i) {
        for (int j = 0; j < nEnuNuMu; ++j) {
          if (GraphsNuMu[i][j] != NULL) delete GraphsNuMu[i][j];
        }
        if (GraphsNuMu[i] != NULL) delete[] GraphsNuMu[i];
      }

      // Delete the nue(bar) graphs
      for (int i = 0; i < kNuMuBar+1; ++i) {
        for (int j = 0; j < nEnuNuE; ++j) {
          if (GraphsNuE[i][j] != NULL) delete GraphsNuE[i][j];
        }
        if (GraphsNuE[i] != NULL) delete[] GraphsNuE[i];
      }
    }

    // Calculate a weight using a fixed Enu, Q2 and neutrino type
    double CalcWeight(double Enu, double Q2);

    // Old method, mostly for documentation. Produces some discontinuities in the interpolation
    double CalcWeightOld(double Enu, double Q2);

    // Setters
    void SetLeptonMass(double input) { leptonmass = input; }; // In GeV

    // Set the neutrino type, and the lepton mass
    void SetNuType(NuType nu) { 
      nutype = nu; 
      if (nu == kNuMu || nu == kNuMuBar) {
        leptonmass = 0.10566;
        nEnu = nEnuNuMu;
        EnuRange = EnuRangeNuMu;
      } else if (nu == kNuE || nu == kNuEBar) {
        leptonmass = 0.511e-3;
        nEnu = nEnuNuE;
        EnuRange = EnuRangeNuE;
      } else {
        std::cerr << "Calculator does not support neutrino type " << nu << std::endl;
        std::cerr << "0 = numu, 1 = numubar, 2 = nue, 3 = nuebar" << std::endl;
        throw;
      }
    };

    void SetLinearInterp() { drawcmd = ""; };
    void SetSplineInterp() { drawcmd = "S"; };

    // Getters
    std::string GetSplineInterp() { return drawcmd; };
    double* GetEnuRange() { return EnuRange; };
    double* GetEnuRangeNuMu() { return EnuRangeNuMu; };
    double* GetEnuRangeNuE() { return EnuRangeNuE; };
    int GetNEnu() { return nEnu; };
    int GetNEnuNuMu() { return nEnuNuMu; };
    int GetNEnuNuE() { return nEnuNuE; };

    // Get the maximum Q2 for a given Enu and muon mass
    // Get the max Q2 for a given Enu to check interpolaton
    // CCQE only
    double GetQ2max(double Enu) {
      // Nucleon mass
      const double Mn = 0.93956542052;
      const double Mp = 0.93827208816;
      const double M = (Mn+Mp)/2.;
      const double M2 = M*M;
      const double leptonmass2 = leptonmass*leptonmass;
      const double Enu2 = Enu*Enu;

      double val = -(M+Enu)*leptonmass2+2*M*Enu2+Enu*sqrt(pow(2*M*Enu-leptonmass2, 2)-4*leptonmass2*M2);
      val /= (M+2*Enu);
      return val;
    }

  private:
    // Input TGraphs, first index is numu or numubar, second index is the fixed enu
    TGraph ***GraphsNuMu;
    TGraph ***GraphsNuE;
    // Have different number of pre-calculated points for numu and nue, so keep this pointing to the relevant one
    TGraph ***Graphs;

    // Enu Range that the inputs come in
    double EnuRangeNuMu[20];
    double EnuRangeNuE[9];
    double *EnuRange;

    // Just a hard-coded check to see if all the graphs are there
    int nEnuNuMu;
    int nEnuNuE;
    int nEnu;

    // Lepton mass in GeV
    double leptonmass;
    NuType nutype;

    // Interpolation method, spline or linear
    std::string drawcmd;

    bool CheckSetup(double &Enu, double &Q2);

    bool Checked;

    // Q2 tolerance, used in old calculation for checking Q2 range
    double Q2tol;
};

#endif
