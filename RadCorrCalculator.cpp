#include "RadCorrCalculator.h"

//#define DEBUG
//#define DEBUG2

// Radcorr weight calculator
RadCorrCalc::RadCorrCalc() {

  // Setup numu objects
  SetupNuMu();
  // Setup nue objects
  SetupNuE();
  // Set up the final graphs
  SetupGraphs();

  // Default the lepton mass to be something silly
  leptonmass = -999;
  nutype = NuType::kNuMu;

  // Setup assuming numu
  SetNuType(nutype);

  drawcmd = ""; // Use linear interpolation

  // Tolerance of Q2 in GeV2
  Q2tol = 1E-2;

  // Has the setup been sanity checked?
  Checked = false;
}

void RadCorrCalc::SetupNuMu() {

  TFile *fInputs[2] = {NULL};

  std::string fileloc = "inputs";
  // Two files for numu and numubar
  // 7 March 2023 use extended tables
  fInputs[kNuMu] = new TFile((fileloc+"/Elspectrum_muon_neutrino_merge_extend2.root").c_str());
  fInputs[kNuMuBar] = new TFile((fileloc+"/Elspectrum_muon_antineutrino_merge_extend2.root").c_str());

  if (!fInputs[kNuMu]->IsOpen() || fInputs[kNuMu]->IsZombie()) {
    std::cerr << "Input file for muon neutrinos does not exist" << std::endl;
    std::cerr << "Provided location: " << fInputs[kNuMu]->GetName() << std::endl;
    throw;
  }

  if (!fInputs[kNuMuBar]->IsOpen() || fInputs[kNuMuBar]->IsZombie()) {
    std::cerr << "Input file for muon anti-neutrinos does not exist" << std::endl;
    std::cerr << "Provided location: " << fInputs[kNuMuBar]->GetName() << std::endl;
    throw;
  }
  // The neutrino energy ranges
  EnuRangeNuMu[0] = 0.2;
  EnuRangeNuMu[1] = 0.3;
  EnuRangeNuMu[2] = 0.5;
  EnuRangeNuMu[3] = 0.6;
  EnuRangeNuMu[4] = 0.75;
  EnuRangeNuMu[5] = 0.9;
  EnuRangeNuMu[6] = 1.0;
  EnuRangeNuMu[7] = 2.0;
  EnuRangeNuMu[8] = 3.0;
  EnuRangeNuMu[9] = 5.0;
  EnuRangeNuMu[10] = 10.;
  EnuRangeNuMu[11] = 15.;
  EnuRangeNuMu[12] = 20.;
  EnuRangeNuMu[13] = 25.;
  EnuRangeNuMu[14] = 30.;
  EnuRangeNuMu[15] = 35.;
  EnuRangeNuMu[16] = 40.;
  EnuRangeNuMu[17] = 45.;
  EnuRangeNuMu[18] = 50.;
  EnuRangeNuMu[19] = 70.;

  // Hard-code the expected number of fixed Enu calculations to check against the size
  nEnuNuMu = 20;
  if (nEnuNuMu != sizeof(EnuRangeNuMu)/sizeof(double)) {
    std::cerr << "The number of neutrino energies expected (" << nEnuNuMu << ") does not match the hard-coded range!" << std::endl;
    std::cerr << "Hard-coded range: " << sizeof(EnuRangeNuMu)/sizeof(double) << std::endl;
    for (unsigned int i = 0; i < sizeof(EnuRangeNuMu)/sizeof(double); ++i) {
      std::cerr << "EnuRangeNuMu[" <<  i << "]=" << EnuRangeNuMu[i] << std::endl;
    }
    throw;
  }

  // Check the number of keys in the file to make sure it matches the number of neutrino energies
  // Might as well be paranoid
  if (nEnuNuMu != fInputs[kNuMu]->GetListOfKeys()->GetSize()) {
    std::cerr << "The number of neutrino energies expected (" << nEnuNuMu << ") does not match the entries in the TFile for muon neutrinos (" << fInputs[kNuMu]->GetName() << ")" << std::endl;
    std::cerr << "Number of entries in file: " << fInputs[kNuMu]->GetListOfKeys()->GetSize() << std::endl;
    throw;
  }
  
  if (nEnuNuMu != fInputs[kNuMuBar]->GetListOfKeys()->GetSize()) {
    std::cerr << "The number of neutrino energies expected (" << nEnuNuMu << ") does not match the entries in the TFile for muon anti-neutrinos (" << fInputs[kNuMuBar]->GetName() << ")" << std::endl;
    std::cerr << "Number of entries in file: " << fInputs[kNuMuBar]->GetListOfKeys()->GetSize() << std::endl;
    throw;
  }

  // The spline base names
  std::string basename_numu = "Elspectrum_muon_neutrino_neutron";
  std::string basename_numub = "Elspectrum_muon_antineutrino_proton";

  GraphsNuMu = new TGraph**[kNuMuBar+1];
  // Loop over neutrino and anti-neutrino
  for (int i = 0; i < kNuMuBar+1; ++i) {
    GraphsNuMu[i] = new TGraph*[nEnuNuMu];
    std::string basename;
    if (i == kNuMu) basename = basename_numu;
    else basename = basename_numub;
    GraphsNuMu[i][0] = (TGraph*)fInputs[i]->Get((basename+"_02GeV").c_str())->Clone();
    GraphsNuMu[i][1] = (TGraph*)fInputs[i]->Get((basename+"_03GeV").c_str())->Clone();
    GraphsNuMu[i][2] = (TGraph*)fInputs[i]->Get((basename+"_05GeV").c_str())->Clone();
    GraphsNuMu[i][3] = (TGraph*)fInputs[i]->Get((basename+"_06GeV").c_str())->Clone();
    GraphsNuMu[i][4] = (TGraph*)fInputs[i]->Get((basename+"_075GeV").c_str())->Clone();
    GraphsNuMu[i][5] = (TGraph*)fInputs[i]->Get((basename+"_09GeV").c_str())->Clone();
    GraphsNuMu[i][6] = (TGraph*)fInputs[i]->Get((basename+"_1GeV").c_str())->Clone();
    GraphsNuMu[i][7] = (TGraph*)fInputs[i]->Get((basename+"_2GeV").c_str())->Clone();
    GraphsNuMu[i][8] = (TGraph*)fInputs[i]->Get((basename+"_3GeV").c_str())->Clone();
    GraphsNuMu[i][9] = (TGraph*)fInputs[i]->Get((basename+"_5GeV").c_str())->Clone();
    GraphsNuMu[i][10] = (TGraph*)fInputs[i]->Get((basename+"_10GeV").c_str())->Clone();
    GraphsNuMu[i][11] = (TGraph*)fInputs[i]->Get((basename+"_15GeV").c_str())->Clone();
    GraphsNuMu[i][12] = (TGraph*)fInputs[i]->Get((basename+"_20GeV").c_str())->Clone();
    GraphsNuMu[i][13] = (TGraph*)fInputs[i]->Get((basename+"_25GeV").c_str())->Clone();
    GraphsNuMu[i][14] = (TGraph*)fInputs[i]->Get((basename+"_30GeV").c_str())->Clone();
    GraphsNuMu[i][15] = (TGraph*)fInputs[i]->Get((basename+"_35GeV").c_str())->Clone();
    GraphsNuMu[i][16] = (TGraph*)fInputs[i]->Get((basename+"_40GeV").c_str())->Clone();
    GraphsNuMu[i][17] = (TGraph*)fInputs[i]->Get((basename+"_45GeV").c_str())->Clone();
    GraphsNuMu[i][18] = (TGraph*)fInputs[i]->Get((basename+"_50GeV").c_str())->Clone();
    GraphsNuMu[i][19] = (TGraph*)fInputs[i]->Get((basename+"_70GeV").c_str())->Clone();
  }

  fInputs[0]->Close();
  fInputs[1]->Close();
}

void RadCorrCalc::SetupNuE() {

  TFile *fInputs[2] = {NULL};

  std::string fileloc = "inputs";
  // Two files for numu and numubar
  // 7 March 2023 use extended tables
  fInputs[kNuMu] = new TFile((fileloc+"/EMspectrum_electron_neutrino_neutron_merge.root").c_str());
  fInputs[kNuMuBar] = new TFile((fileloc+"/EMspectrum_electron_antineutrino_proton_merge.root").c_str());

  if (!fInputs[kNuMu]->IsOpen() || fInputs[kNuMu]->IsZombie()) {
    std::cerr << "Input file for muon neutrinos does not exist" << std::endl;
    std::cerr << "Provided location: " << fInputs[kNuMu]->GetName() << std::endl;
    throw;
  }

  if (!fInputs[kNuMuBar]->IsOpen() || fInputs[kNuMuBar]->IsZombie()) {
    std::cerr << "Input file for muon anti-neutrinos does not exist" << std::endl;
    std::cerr << "Provided location: " << fInputs[kNuMuBar]->GetName() << std::endl;
    throw;
  }

  // The neutrino energy ranges
  EnuRangeNuE[0] = 0.2;
  EnuRangeNuE[1] = 0.3;
  EnuRangeNuE[2] = 0.5;
  EnuRangeNuE[3] = 1.0;
  EnuRangeNuE[4] = 2.0;
  EnuRangeNuE[5] = 3.0;
  EnuRangeNuE[6] = 5.0;
  EnuRangeNuE[7] = 10.0;
  EnuRangeNuE[8] = 20.0;

  // Hard-code the expected number of fixed Enu calculations to check against the size
  nEnuNuE = 9;
  if (nEnuNuE != sizeof(EnuRangeNuE)/sizeof(double)) {
    std::cerr << "The number of neutrino energies expected (" << nEnuNuE << ") does not match the hard-coded range!" << std::endl;
    std::cerr << "Hard-coded range: " << sizeof(EnuRangeNuE)/sizeof(double) << std::endl;
    for (unsigned int i = 0; i < sizeof(EnuRangeNuE)/sizeof(double); ++i) {
      std::cerr << "EnuRangeNuE[" <<  i << "]=" << EnuRangeNuE[i] << std::endl;
    }
    throw;
  }

  // Check the number of keys in the file to make sure it matches the number of neutrino energies
  // Might as well be paranoid
  if (nEnuNuE != fInputs[kNuMu]->GetListOfKeys()->GetSize()) {
    std::cerr << "The number of neutrino energies expected (" << nEnuNuE << ") does not match the entries in the TFile for electron neutrinos (" << fInputs[kNuMu]->GetName() << ")" << std::endl;
    std::cerr << "Number of entries in file: " << fInputs[kNuMu]->GetListOfKeys()->GetSize() << std::endl;
    throw;
  }
  
  if (nEnuNuE != fInputs[kNuMuBar]->GetListOfKeys()->GetSize()) {
    std::cerr << "The number of neutrino energies expected (" << nEnuNuE << ") does not match the entries in the TFile for electron anti-neutrinos (" << fInputs[kNuMuBar]->GetName() << ")" << std::endl;
    std::cerr << "Number of entries in file: " << fInputs[kNuMuBar]->GetListOfKeys()->GetSize() << std::endl;
    throw;
  }

  // The spline base names
  std::string basename_nue = "EMspectrum_electron_neutrino_neutron";
  std::string basename_nueb = "EMspectrum_electron_antineutrino_proton";

  GraphsNuE = new TGraph**[kNuMuBar+1];
  // Loop over neutrino and anti-neutrino
  for (int i = 0; i < kNuMuBar+1; ++i) {
    GraphsNuE[i] = new TGraph*[nEnuNuE];
    std::string basename;
    if (i == kNuMu) basename = basename_nue;
    else basename = basename_nueb;
    GraphsNuE[i][0] = (TGraph*)fInputs[i]->Get((basename+"_02GeV").c_str())->Clone();
    GraphsNuE[i][1] = (TGraph*)fInputs[i]->Get((basename+"_03GeV").c_str())->Clone();
    GraphsNuE[i][2] = (TGraph*)fInputs[i]->Get((basename+"_05GeV").c_str())->Clone();
    GraphsNuE[i][3] = (TGraph*)fInputs[i]->Get((basename+"_1GeV").c_str())->Clone();
    GraphsNuE[i][4] = (TGraph*)fInputs[i]->Get((basename+"_2GeV").c_str())->Clone();
    GraphsNuE[i][5] = (TGraph*)fInputs[i]->Get((basename+"_3GeV").c_str())->Clone();
    GraphsNuE[i][6] = (TGraph*)fInputs[i]->Get((basename+"_5GeV").c_str())->Clone();
    GraphsNuE[i][7] = (TGraph*)fInputs[i]->Get((basename+"_10GeV").c_str())->Clone();
    GraphsNuE[i][8] = (TGraph*)fInputs[i]->Get((basename+"_20GeV").c_str())->Clone();
  }
  fInputs[0]->Close();
  fInputs[1]->Close();
}

void RadCorrCalc::SetupGraphs() {
  Graphs = new TGraph**[kNuEBar+1];
  Graphs[0] = new TGraph*[nEnuNuMu];
  Graphs[1] = new TGraph*[nEnuNuMu];
  Graphs[2] = new TGraph*[nEnuNuE];
  Graphs[3] = new TGraph*[nEnuNuE];

  for (int i = 0; i < 2; ++i) {
    for (int j = 0; j < nEnuNuMu; ++j) {
      Graphs[i][j] = GraphsNuMu[i][j];
    }
  }
  for (int i = 2; i < 4; ++i) {
    for (int j = 0; j < nEnuNuE; ++j) {
      Graphs[i][j] = GraphsNuE[i-2][j];
    }
  }
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

  // Check if this Enu/Q2 combination is allowed, and that everything has been setup
  if (!CheckSetup(Enu, Q2)) return 1.0;

  // Find the nearest graph point in Enu
  // Get the true neutrino energy and interpolate between nearest points
  // This is always the point *below* provided Enu
  int nearest = 0;
  for (int j = 0; j < nEnu; ++j) {
    if (Enu > EnuRange[j]) nearest = j;
  }

  // Then get the maximum Q2 for each of the Enu points 
  // to make sure we're not extrapolating unphysically
  double Q2max_near = GetQ2max(EnuRange[nearest]);

  // Find the Q2 point before the boundary, see if it's close
  // Get the TGraph for the lower point in Enu
  TGraph *g = Graphs[nutype][nearest];
  const double *x = g->GetX();
  const double *y = g->GetY();
  const int npoints = g->GetN();

  // Find where the discontinuity happens in the input
  // This should be around the maximum Q2 for the fixed neutrino energy
  int q2_jump_point = -1;
  double old_max = y[0];
  for (int j = 0; j < npoints; ++j) {
    // If we've passed the maximum Q2 point
    if (x[j] > Q2max_near) {
      q2_jump_point = j-1;
      break;
    }
  }
  old_max = y[q2_jump_point];

#ifdef DEBUG2
  std::cout << "Enu: " << Enu << " Q2: " << Q2 << std::endl;
#endif
  double min_deltadelta = 100;
  double max_deltadelta = -999;
  double min_m1 = 100; // Closest to zero
  double min_p1 = 100; // Average delta minimum (most negative)
  double max_weight = 0;

  // Point of minimum first derivative
  int point_min_m1 = 0;
  // Point of minimum first derivative on N+1 point
  int point_min_p1 = 0;
  // Point of minimum second derivative
  int point_min_2 = 0;
  // Point of maximum second derivative
  int point_max_2 = 0;
  // Point of identified jump
  int point_jump = 0;
  // Point of maximum weight
  int point_max_weight = 0;

  // Loop over points, starting at second point going to second to last point
  // Needed to second derivative using N-1 and N+1th points
  for (int j = 1; j < npoints-1; ++j) {

    // Update maximum weight
    if (y[j] > max_weight) {
      max_weight = y[j];
      point_max_weight = j;
    }
    // First differentials (positive and negative direction)
    double deltap1 = (y[j+1]-y[j])/(x[j+1]-x[j]);
    double deltam1 = (y[j]-y[j-1])/(x[j]-x[j-1]);

    // Second differential
    double deltadelta = (deltap1-deltam1)/(x[j+1]-x[j-1]);

    // Find the smallest first derivative; look for the derivative relative the previous point, not the future point (comes next)
    if (fabs(deltam1) < fabs(min_m1) && y[j] > 0 && deltam1 > 0) {
      min_m1 = fabs(deltam1);
      point_min_m1 = j;
    }

    // Look at derivative in positive direction
    if (fabs(deltap1) < fabs(min_p1) && y[j] > 0 && deltap1 < 0) {
      min_p1 = fabs(deltap1);
      point_min_p1 = j;
    }

    // Find the smallest second derivative
    if (fabs(deltadelta) < fabs(min_deltadelta) && y[j] > 0) {
      min_deltadelta = fabs(deltadelta);
      point_min_2 = j;
    }

    // Find the largest second derivative
    if (fabs(deltadelta) > fabs(max_deltadelta) && y[j] > 0) {
      max_deltadelta = fabs(deltadelta);
      point_max_2 = j;
    }

    // If deltap1 and delta m1 have opposite sign, likely at a maximum (inflection)
    // This should also agree with the maximum second differential

#ifdef DEBUG2
    std::cout << "x[" << std::setw(2) << j << "]=" << std::setw(10) << x[j] << " y[" << std::setw(2) << j << "]=" << std::setw(10) << y[j] << ", deltap1=" << std::setw(15) << deltap1 << ", deltam1=" << std::setw(15) << deltam1 << ", deltadelta: " << std::setw(15) << deltadelta << std::endl;
#endif


    //if (fabs(deltadelta) > 0.1 && point_jump == 0 && j > point_max && fabs(deltap1) > 0.01 && fabs(deltam1) > 0.01) {
    //if (fabs(deltadelta) > 0.01 && j > point_max && deltap1 < -0.01 && fabs(deltam1) > 0.01) {

    // Find when the second order differential becomes larger than 0.01
    // and we've got negative slope at N+1 direction, and we had a noticeable slope in N-1 direction
    // Remove point_max check in case jump happens at maximum point
    if (fabs(deltadelta) > 0.01 && deltap1 < -0.01 && fabs(deltam1) > 0.01) {
      point_jump = j;
    }
  }

  // Use the identified point of the jump if identified
  if (point_jump > 0) {
    point_jump = point_min_p1;
  // If a jump point was not identified by the differentials, use the maximum weight
  } else { 
    point_jump = point_max_weight;
  }

  // Check if the identified jump point is above allowed Q2
  if (point_jump > q2_jump_point && q2_jump_point > 0) {
    point_jump = q2_jump_point;
  }
#ifdef DEBUG2
  std::cout << "*************************" << std::endl;
  std::cout << "Max Q2: " << Q2max_near << " at point " << q2_jump_point << std::endl;
  std::cout << std::setw(25) << "Smallest 2nd diff: " <<  std::setw(10) << min_deltadelta << ", at point " << point_min_2 << ", x=" << x[point_min_2] << std::endl;
  std::cout << std::setw(25) << "Largest 2nd diff: " <<  std::setw(10) << max_deltadelta << ", at point " << point_max_2 << ", x=" << x[point_max_2] << std::endl;
  std::cout << std::setw(25) << "Smallest 1st diff m1: " << std::setw(10) << min_m1 << ", at point " << point_min_m1 << ", x=" << x[point_min_m1] << std::endl;
  std::cout << std::setw(25) << "Smallest 1st diff p1: " << std::setw(10) << min_p1 << ", at point " << point_min_p1 << ", x=" << x[point_min_p1] << std::endl;
  std::cout << std::setw(25) << "Maximum: " << std::setw(10) << max_weight << ", at point " << point_max_weight << ", x=" << x[point_max_weight] << std::endl;
  std::cout << std::setw(25) << "Identified jump point: " << std::setw(10) << " " << ", at point " << point_jump << ", x=" << x[point_jump] << std::endl;
  std::cout << "*************************" << std::endl;
#endif


  // If Q2 is greater than where the jump happens, trace back
  //if (Q2 > x point_jump > 0) {
    //high = y[point_jump];

  // The allowed Q2 is always going to be lower at the lower energy bin
  // If this is the case, evaluate the spline just before the drop off
  double low = 0;
  // Sometimes the Q2 will be right on the edge of allowed Q2 phase space
  //if (jump_point > 0 && Q2 > x[jump_point-1]) {
    //low = y[jump_point-1];
  //} else if (jump_point == -1 && Q2 > x[jump_point-1]) {
    //low = y[jump_point-1];
  if (Q2 > x[point_jump] && point_jump > 0) {
    low = y[point_jump];
  // Jump point not successfully identified but Q2 is larger than allowed -> use last point
  //} else if (Q2 > x[npoints] && point_jump == 0) {
    //low = y[npoints-1];
    //std::cout << low << std::endl;
  } else {
    low = g->Eval(Q2, 0, drawcmd.c_str());
  }
#ifdef DEBUG2
  std::cout << "Enu near: " << EnuRange[nearest] << " Q2max near: " << Q2max_near << " weight: " << low << std::endl;
#endif

  // Now move on to the next point (the upper edge of Enu)

  // This is the point above
  int nextbin = nearest+1;

  // The Enu might be beyond our last range
  if (Enu > EnuRange[nEnu-1]) nextbin = nearest;

  // Now also need to check the high Q2
  double Q2max_far = GetQ2max(EnuRange[nextbin]);
  // Find the point 
  TGraph *g_far = Graphs[nutype][nextbin];
  const int npoints_far = g_far->GetN();
  const double *x_far = g_far->GetX();
  const double *y_far = g_far->GetY();
  int jump_point_far = -1;
  old_max = y_far[0];
  for (int j = 0; j < npoints_far; ++j) {
    if (x_far[j] >= Q2max_far) {
      jump_point_far = j;
      break;
    }
    old_max = y_far[j];
  }

  double high = 0;

  // Sometimes the Q2 will be right on the edge of allowed
  if (jump_point_far > 0 && Q2 > x_far[jump_point_far-1]) {
    high = y_far[jump_point_far-1];
    // Sometimes an event sits just on the Q2 boundary and the precalculated point is right there
    // Looks like the input is sometimes miscalculated here
    if (high == 0) high = y_far[jump_point_far-2];
    // If Q2 is larger than our maximum calculation and we haven't seen a jump yet
  } else if (jump_point_far == -1 && Q2 > x_far[npoints_far-1]) {
    high = y_far[npoints_far-1];
  } else {
    high = g_far->Eval(Q2, 0, drawcmd.c_str());
  }


  // linear intepolation
  double weight = (high-low)*(Enu-EnuRange[nearest])/(EnuRange[nextbin]-EnuRange[nearest])+low;

  // Put in a weight cap
  if (weight < 0) weight = 0;
#ifdef DEBUG2
  std::cout << "Enu far: " << EnuRange[nextbin] << " Q2max far: " << Q2max_far << " weight: " << high << std::endl;
  std::cout << "Overall weight: " << weight << std::endl;
#endif

  return weight;
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
  if (Enu < EnuRangeNuMu[0]) {
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
  for (int j = 0; j < nEnuNuMu; ++j) {
    if (Enu > EnuRangeNuMu[j]) nearest = j;
  }
  std::cout << "Enu: " << Enu << " " << " nearest: " << EnuRangeNuMu[nearest] << " which has index " << nearest << " next: " << EnuRangeNuMu[nearest+1] << std::endl;

  // Then get the maximum Q2 for each of the Enu points 
  // to make sure we're not extrapolating unphysically
  double Q2max_near = GetQ2max(EnuRangeNuMu[nearest]);

  // The allowed Q2 is always going to be lower at the lower energy bin
  // If this is the case, evaluate the spline just before the drop off
  double low = 0;
  // Sometimes the Q2 will be right on the edge of allowed Q2 phase space
  // CWRET: maybe this is the problem?
  if (Q2 > Q2max_near-Q2tol) {
    TGraph *g = GraphsNuMu[nutype][nearest]; // Get the TGraph for the lower point in Enu
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
    low = GraphsNuMu[nutype][nearest]->Eval(Q2, 0, drawcmd.c_str());
  }

  // This is the point above
  int nextbin = nearest+1;

  // The Enu might be beyond our last range
  if (Enu > EnuRangeNuMu[nEnuNuMu-1]) nextbin = nearest;

  // Now also need to check the high Q2
  double Q2max_far = GetQ2max(EnuRangeNuMu[nextbin]);
  double high = 0;
  // Sometimes the Q2 will be right on the edge of allowed
  if (Q2 > Q2max_far-Q2tol) {
    // Try to find the spot where we're no longer dropping abruptly
    TGraph *g = GraphsNuMu[nutype][nextbin];
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
    high = GraphsNuMu[nutype][nextbin]->Eval(Q2, 0, drawcmd.c_str());
  }

  // linear intepolation
  double weight = (high-low)*(Enu-EnuRangeNuMu[nearest])/(EnuRangeNuMu[nearest+1]-EnuRangeNuMu[nearest])+low;

  // Put in a weight cap
  if (weight > 10) weight = 10;
  if (weight < 0) weight = 0;

  return weight;
}
