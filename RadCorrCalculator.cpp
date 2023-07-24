#include "RadCorrCalculator.h"
//#define DEBUG
//

// Radcorr weight calculator
RadCorrCalc::RadCorrCalc() {
  TFile *fInputs[2];

  std::string fileloc = "inputs";
  // Two files for numu and numubar
  // 7 March 2023 use extended tables
  fInputs[kNumu] = new TFile((fileloc+"/Elspectrum_muon_neutrino_merge_extend2.root").c_str());
  fInputs[kNumuBar] = new TFile((fileloc+"/Elspectrum_muon_antineutrino_merge_extend2.root").c_str());

  if (!fInputs[kNumu]->IsOpen() || fInputs[kNumu]->IsZombie()) {
    std::cerr << "Input file for muon neutrinos does not exist" << std::endl;
    std::cerr << "Provided location: " << fInputs[kNumu]->GetName() << std::endl;
    throw;
  }

  if (!fInputs[kNumuBar]->IsOpen() || fInputs[kNumuBar]->IsZombie()) {
    std::cerr << "Input file for muon anti-neutrinos does not exist" << std::endl;
    std::cerr << "Provided location: " << fInputs[kNumuBar]->GetName() << std::endl;
    throw;
  }

  // The neutrino energy ranges
  EnuRange[0] = 0.2;
  EnuRange[1] = 0.3;
  EnuRange[2] = 0.5;
  EnuRange[3] = 0.6;
  EnuRange[4] = 0.75;
  EnuRange[5] = 0.9;
  EnuRange[6] = 1.0;
  EnuRange[7] = 2.0;
  EnuRange[8] = 3.0;
  EnuRange[9] = 5.0;
  EnuRange[10] = 10.;
  EnuRange[11] = 15.;
  EnuRange[12] = 20.;
  EnuRange[13] = 25.;
  EnuRange[14] = 30.;
  EnuRange[15] = 35.;
  EnuRange[16] = 40.;
  EnuRange[17] = 45.;
  EnuRange[18] = 50.;
  EnuRange[19] = 70.;

  // Hard-code the expected number of fixed Enu calculations to check against the size
  nEnu = 20;
  if (nEnu != sizeof(EnuRange)/sizeof(double)) {
    std::cerr << "The number of neutrino energies expected (" << nEnu << ") does not match the hard-coded range!" << std::endl;
    std::cerr << "Hard-coded range: " << sizeof(EnuRange)/sizeof(double) << std::endl;
    for (unsigned int i = 0; i < sizeof(EnuRange)/sizeof(double); ++i) {
      std::cerr << "EnuRange[" <<  i << "]=" << EnuRange[i] << std::endl;
    }
    throw;
  }

  // Check the number of keys in the file to make sure it matches the number of neutrino energies
  // Might as well be paranoid
  if (nEnu != fInputs[kNumu]->GetListOfKeys()->GetSize()) {
    std::cerr << "The number of neutrino energies expected (" << nEnu << ") does not match the entries in the TFile for muon neutrinos (" << fInputs[kNumu]->GetName() << ")" << std::endl;
    std::cerr << "Number of entries in file: " << fInputs[kNumu]->GetListOfKeys()->GetSize() << std::endl;
    throw;
  }
  
  if (nEnu != fInputs[kNumuBar]->GetListOfKeys()->GetSize()) {
    std::cerr << "The number of neutrino energies expected (" << nEnu << ") does not match the entries in the TFile for muon anti-neutrinos (" << fInputs[kNumuBar]->GetName() << ")" << std::endl;
    std::cerr << "Number of entries in file: " << fInputs[kNumuBar]->GetListOfKeys()->GetSize() << std::endl;
    throw;
  }

  // The spline base names
  std::string basename_numu = "Elspectrum_muon_neutrino_neutron";
  std::string basename_numub = "Elspectrum_muon_antineutrino_proton";

  Graphs = new TGraph**[kNumuBar+1];
  // Loop over neutrino and anti-neutrino
  for (int i = 0; i < kNumuBar+1; ++i) {
    Graphs[i] = new TGraph*[nEnu];
    std::string basename;
    if (i == kNumu) basename = basename_numu;
    else basename = basename_numub;
    Graphs[i][0] = (TGraph*)fInputs[i]->Get((basename+"_02GeV").c_str())->Clone();
    Graphs[i][1] = (TGraph*)fInputs[i]->Get((basename+"_03GeV").c_str())->Clone();
    Graphs[i][2] = (TGraph*)fInputs[i]->Get((basename+"_05GeV").c_str())->Clone();
    Graphs[i][3] = (TGraph*)fInputs[i]->Get((basename+"_06GeV").c_str())->Clone();
    Graphs[i][4] = (TGraph*)fInputs[i]->Get((basename+"_075GeV").c_str())->Clone();
    Graphs[i][5] = (TGraph*)fInputs[i]->Get((basename+"_09GeV").c_str())->Clone();
    Graphs[i][6] = (TGraph*)fInputs[i]->Get((basename+"_1GeV").c_str())->Clone();
    Graphs[i][7] = (TGraph*)fInputs[i]->Get((basename+"_2GeV").c_str())->Clone();
    Graphs[i][8] = (TGraph*)fInputs[i]->Get((basename+"_3GeV").c_str())->Clone();
    Graphs[i][9] = (TGraph*)fInputs[i]->Get((basename+"_5GeV").c_str())->Clone();
    Graphs[i][10] = (TGraph*)fInputs[i]->Get((basename+"_10GeV").c_str())->Clone();
    Graphs[i][11] = (TGraph*)fInputs[i]->Get((basename+"_15GeV").c_str())->Clone();
    Graphs[i][12] = (TGraph*)fInputs[i]->Get((basename+"_20GeV").c_str())->Clone();
    Graphs[i][13] = (TGraph*)fInputs[i]->Get((basename+"_25GeV").c_str())->Clone();
    Graphs[i][14] = (TGraph*)fInputs[i]->Get((basename+"_30GeV").c_str())->Clone();
    Graphs[i][15] = (TGraph*)fInputs[i]->Get((basename+"_35GeV").c_str())->Clone();
    Graphs[i][16] = (TGraph*)fInputs[i]->Get((basename+"_40GeV").c_str())->Clone();
    Graphs[i][17] = (TGraph*)fInputs[i]->Get((basename+"_45GeV").c_str())->Clone();
    Graphs[i][18] = (TGraph*)fInputs[i]->Get((basename+"_50GeV").c_str())->Clone();
    Graphs[i][19] = (TGraph*)fInputs[i]->Get((basename+"_70GeV").c_str())->Clone();
  }

  fInputs[0]->Close();
  fInputs[1]->Close();

  // Default the lepton mass to be something silly
  leptonmass = -999;
  nutype = NuType::kNumu;
  drawcmd = ""; // Use linear interpolation

  // Tolerance of Q2 in GeV2
  Q2tol = 1E-2;

  Checked = false;
}

bool RadCorrCalc::CheckSetup(double &Enu, double &Q2) {

  if (!Checked) {
    if (leptonmass < 0) {
      std::cerr << "Nonsensical leptonmass: " << leptonmass << std::endl;
      std::cerr << "Did you not call SetLeptonMass(double) before calling CalcWeight?" << std::endl;
      throw;
    }

    if (nutype < 0) {
      std::cerr << "Nonsensical nutype: " << nutype << std::endl;
      std::cerr << "Did you not call SetNuType(int) before calling CalcWeight?" << std::endl;
      throw;
    }
  }
  Checked = true;

  // Also check that this Q2 is allowed for this Enu
  double maxq2 = GetQ2max(Enu);
  if (Q2 > maxq2) {
#ifdef DEBUG
    std::cout << "Q2 " << Q2 << " is above allowed Q2 " << maxq2 << " for Enu " << Enu << std::endl;
    std::cout << "returning weight 1" << std::endl;
#endif
    return false;
  }

  // Check Enu is above minimum Enu
  // Could probably use the calculation for 0.2 here?
  if (Enu < EnuRange[0]) {
#ifdef DEBUG
    std::cout << "Enu less than minimum pre-calculated, return 1" << std::endl;
#endif
    return false;
  }

  // check non-zero Q2
  // Could also check against the minimum Q2?
  if (Q2 < 1E-5) {
#ifdef DEBUG
    std::cout << "Too low Q2, returning 1" << std::endl;
#endif
    return false;
  }

  // Check if neutrino energy is above maximum available
  // If so, proceed but with neutrino energy set to maximum
  if (Enu > EnuRange[nEnu-1]) {
    Enu = EnuRange[nEnu-1];
  } else if (Enu < EnuRange[0]) {
    Enu = EnuRange[0];
  }

  return true;
}

// Calculate the weight from the radiative correction in Q2
// Function of Enu and Q2, using linear interpolation
double RadCorrCalc::CalcWeight(double Enu, double Q2) {

  // Check if allowed Enu/Q2 combination, and class been setup
  if (!CheckSetup(Enu, Q2)) return 1.0;

  // Find the nearest graph point in Enu
  // Get the true neutrino energy and interpolate between nearest points
  // This is always the point below our provided Enu
  int nearest = 0;
  for (int j = 0; j < nEnu; ++j) {
    if (Enu > EnuRange[j]) nearest = j;
  }

  // Then get the maximum Q2 for each of the Enu points 
  // to make sure we're not extrapolating unphysically
  double Q2max_near = GetQ2max(EnuRange[nearest]);

  // Find the Q2 point before the boundary, see if it's close
  TGraph *g = Graphs[nutype][nearest]; // Get the TGraph for the lower point in Enu
  const double *x = g->GetX();
  const double *y = g->GetY();
  const int npoints = g->GetN();
  int jump_point = -1;
  double old_max = y[0];
  for (int j = 0; j < npoints; ++j) {
    if (x[j] > Q2max_near || fabs(y[j] - old_max) > 0.9) {
      jump_point = j;
      break;
    }
    old_max = y[j];
  }

  // The allowed Q2 is always going to be lower at the lower energy bin
  // If this is the case, evaluate the spline just before the drop off
  double low = 0;
  // Sometimes the Q2 will be right on the edge of allowed Q2 phase space
  if (jump_point > 0 && Q2 > x[jump_point-1]) {
    low = y[jump_point-1];
  } else {
    low = g->Eval(Q2, 0, drawcmd.c_str());
  }

  // This is the point above
  int nextbin = nearest+1;

  // The Enu might be beyond our last range
  if (Enu > EnuRange[nEnu-1]) nextbin = nearest;

  // Now also need to check the high Q2
  double Q2max_far = GetQ2max(EnuRange[nextbin]);
  // Find the point 
  TGraph *g_far = Graphs[nutype][nextbin];
  const int npoints_far = g->GetN();
  const double *x_far = g->GetX();
  const double *y_far = g->GetY();
  int jump_point_far = -1;
  for (int j = 0; j < npoints_far; ++j) {
    if (x_far[j] > Q2max_far) {
      jump_point_far = j;
      break;
    }
  }

  double high = 0;
  // Sometimes the Q2 will be right on the edge of allowed
  if (jump_point_far > 0 && Q2 > x_far[jump_point_far-1]) {
    high = y_far[jump_point_far-1];
  } else {
    high = g_far->Eval(Q2, 0, drawcmd.c_str());
  }

  // linear intepolation
  double weight = (high-low)*(Enu-EnuRange[nearest])/(EnuRange[nextbin]-EnuRange[nearest])+low;

  // Put in a weight cap
  if (weight > 10) weight = 10;
  if (weight < 0) weight = 0;

  return weight;
}

// Get the max Q2 for a given Enu to check interpolaton
// CCQE only
double RadCorrCalc::GetQ2max(double Enu) {
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

// Just delete the graphs
RadCorrCalc::~RadCorrCalc() {
  for (int i = 0; i < kNumuBar+1; ++i) {
    for (int j = 0; j < nEnu; ++j) {
      delete Graphs[i][j];
    }
    delete[] Graphs[i];
  }
}

// Calculate the weight from the radiative correction in Q2
// Function of Enu and Q2, using linear interpolation
double RadCorrCalc::CalcWeightOld(double Enu, double Q2) {

  if (leptonmass < 0) {
    std::cerr << "Nonsensical leptonmass: " << leptonmass << std::endl;
    std::cerr << "Did you not call SetLeptonMass(double) before calling CalcWeight?" << std::endl;
    throw;
  }

  if (nutype < 0) {
    std::cerr << "Nonsensical nutype: " << nutype << std::endl;
    std::cerr << "Did you not call SetNuType(int) before calling CalcWeight?" << std::endl;
    throw;
  }

  // Also check that this Q2 is allowed for this Enu
  double maxq2 = GetQ2max(Enu);
  if (Q2 > maxq2) {
#ifdef DEBUG
    std::cout << "Q2 " << Q2 << " is above allowed Q2 " << maxq2 << " for Enu " << Enu << std::endl;
    std::cout << "returning weight 1" << std::endl;
#endif
    return 1;
  }

  // Check Enu is above minimum Enu
  // Could probably use the calculation for 0.2 here?
  if (Enu < EnuRange[0]) {
#ifdef DEBUG
    std::cout << "Enu less than minimum pre-calculated, return 1" << std::endl;
#endif
    return 1;
  }

  // check non-zero Q2
  // Could also check against the minimum Q2?
  if (Q2 < 1E-5) {
#ifdef DEBUG
    std::cout << "Too low Q2, returning 1" << std::endl;
#endif
    return 1;
  }

  // Find the nearest graph point in Enu
  // Get the true neutrino energy and interpolate between nearest points
  // This is always the point below our provided Enu
  int nearest = 0;
  for (int j = 0; j < nEnu; ++j) {
    if (Enu > EnuRange[j]) nearest = j;
  }
  std::cout << "Enu: " << Enu << " " << " nearest: " << EnuRange[nearest] << " which has index " << nearest << " next: " << EnuRange[nearest+1] << std::endl;

  // Then get the maximum Q2 for each of the Enu points 
  // to make sure we're not extrapolating unphysically
  double Q2max_near = GetQ2max(EnuRange[nearest]);

  // The allowed Q2 is always going to be lower at the lower energy bin
  // If this is the case, evaluate the spline just before the drop off
  double low = 0;
  // Sometimes the Q2 will be right on the edge of allowed Q2 phase space
  // CWRET: maybe this is the problem?
  if (Q2 > Q2max_near-Q2tol) {
    TGraph *g = Graphs[nutype][nearest]; // Get the TGraph for the lower point in Enu
    const double *y = g->GetY();
    int npoints = g->GetN();
    // Try to find the spot where we're no longer dropping abruptly
    double prevy = y[0]; // Start at the first y value
                         // Which point does the jump occur at
    int jump = 0;
    double maximum_y = 0;
    // Find at what point the jump happens
    for (int j = 0; j < npoints; ++j) {
      // Find the maximum y[j]
      if (y[j] > maximum_y) {
        maximum_y = y[j];
      }
      // Look for jumps in the points
      if (prevy-y[j] > 0.9) {
        jump = j;
        break;
      }
      prevy = y[j];
    }
    low = y[jump-1];
    //low = maximum_y;
  } else {
    low = Graphs[nutype][nearest]->Eval(Q2, 0, drawcmd.c_str());
  }

  // This is the point above
  int nextbin = nearest+1;

  // The Enu might be beyond our last range
  if (Enu > EnuRange[nEnu-1]) nextbin = nearest;

  // Now also need to check the high Q2
  double Q2max_far = GetQ2max(EnuRange[nextbin]);
  double high = 0;
  // Sometimes the Q2 will be right on the edge of allowed
  if (Q2 > Q2max_far-Q2tol) {
    // Try to find the spot where we're no longer dropping abruptly
    TGraph *g = Graphs[nutype][nextbin];
    double *y = g->GetY();
    int npoints = g->GetN();
    double prevy = y[0];
    // Which point does the job occur at
    int jump = 0;
    double maximum_y = 0;
    // Find at what point the jump happens
    for (int j = 0; j < npoints; ++j) {
      // Find the maximum y[j]
      if (y[j] > maximum_y) {
        maximum_y = y[j];
      }
      // Look for jumps in the points
      if (prevy-y[j] > 0.9) {
        jump = j;
        break;
      }
      prevy = y[j];
    }
    high = y[jump-1];
    //high = maximum_y;
  } else {
    high = Graphs[nutype][nextbin]->Eval(Q2, 0, drawcmd.c_str());
  }

  // linear intepolation
  double weight = (high-low)*(Enu-EnuRange[nearest])/(EnuRange[nearest+1]-EnuRange[nearest])+low;

  // Put in a weight cap
  if (weight > 10) weight = 10;
  if (weight < 0) weight = 0;

  return weight;
}
