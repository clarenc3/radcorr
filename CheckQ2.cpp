#include <iostream>
#include <iomanip>
#include "RadCorrCalculator.h"

// Simple application to check what the weight is if Q2 is right on the Q2_{max} point
int main() {
  RadCorrCalc calc;

  calc.SetNuType(NuType::kNuMuBar); // Set the neutrino type, 0 = numu, 1 = numubar

  for (int type = 0; type < kNuEBar+1; ++type) {
    calc.SetNuType(NuType(type)); // Set the neutrino type, 0 = numu, 1 = numubar
    std::cout << "-----------------------------" << std::endl;
    std::cout << "NuType: " << type << std::endl;
    int nEnu = calc.GetNEnu();
    double *EnuRange = calc.GetEnuRange();
    for (int i = 0; i < nEnu; ++i) {
      std::cout << "-----------------------------" << std::endl;
      double enu = EnuRange[i];
      double q2max = calc.GetQ2max(EnuRange[i]);
      double weight = calc.CalcWeight(enu, q2max);
      std::cout << "Enu: " << std::setw(5) << enu << " | Q2max: " << std::setw(8) << std::setprecision(5) << q2max << " | Weight: " << std::setw(8) << std::setprecision(5) << weight << std::endl;
    }
  }
}
