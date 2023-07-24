#ifndef _RADCORRCALC_H_SEEN_
#define _RADCORRCALC_H_SEEN_

#include <iostream>
#include "TGraph.h"
#include "TFile.h"

// Just a handy enum
enum NuType {kNumu = 0, kNumuBar = 1};

class RadCorrCalc {
  public:
    // Constructor and destructor
    RadCorrCalc();
    ~RadCorrCalc();

    // Calculate a weight using a fixed Enu, Q2 and neutrino type
    double CalcWeight(double Enu, double Q2);

    // Old method, mostly for documentation. Produces some discontinuities in the interpolation
    double CalcWeightOld(double Enu, double Q2);

    // Setters
    void SetLeptonMass(double input) { leptonmass = input; }; // In GeV
    void SetNuType(NuType nu) { nutype = nu; };
    void SetLinearInterp() { drawcmd = ""; };
    void SetSplineInterp() { drawcmd = "S"; };

    // Getters
    std::string GetSplineInterp() { return drawcmd; };
    double* GetEnuRange() { return EnuRange; };
    int GetNEnu() { return nEnu; };
    // Get the maximum Q2 for a given Enu and muon mass
    double GetQ2max(double Enu);

  private:
    // Input TGraphs, first index is numu or numubar, second index is the fixed enu
    TGraph ***Graphs;

    // Enu Range that the inputs come in
    double EnuRange[20];
    // Just a hard-coded check to see if all the graphs are there
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
