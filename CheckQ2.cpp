#include <iostream>
#include <iomanip>
#include "RadCorrCalculator.h"

// Simple application to check the weight given right on the Q2 maximum border
int main() {
  RadCorrCalc calc;
  calc.SetLeptonMass(0.10566);
  double *EnuRange = calc.GetEnuRange();
  int nEnu = calc.GetNEnu();
  for (int i = 0; i < nEnu; ++i) {
    std::cout << "***" << std::endl;
    double enu = EnuRange[i];
    double q2max = calc.GetQ2max(EnuRange[i]);
    double weight = calc.CalcWeight(enu, q2max);
    std::cout << "Enu: " << std::setw(5) << enu << " | Q2max: " << std::setw(8) << std::setprecision(5) << q2max << " | Weight: " << std::setw(8) << std::setprecision(5) << weight << std::endl;
  }
}
