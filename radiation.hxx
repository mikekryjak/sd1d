

#ifndef __RADIATION_H__
#define __RADIATION_H__

#include <field3d.hxx>
#include <bout_types.hxx>

#include <vector>
#include <cmath>

class RadiatedPower {
public:
  const Field3D power(const Field3D &Te, const Field3D &Ne, const Field3D &Ni);
  
  virtual BoutReal power(BoutReal Te, BoutReal ne, BoutReal ni) = 0;
  
private:
};

class InterpRadiatedPower : public RadiatedPower {
public:
  InterpRadiatedPower(const std::string &file);
  
  BoutReal power(BoutReal Te, BoutReal ne, BoutReal ni);
  
private:
  std::vector<BoutReal> te_array;  // Te in eV
  std::vector<BoutReal> p_array;   // Radiative loss rate in Watts m^3
};


/// Rates supplied by Eva Havlicova
class HydrogenRadiatedPower : public RadiatedPower {
public:
  BoutReal power(BoutReal Te, BoutReal ne, BoutReal ni);
  
  // Collision rate coefficient <sigma*v> [m3/s]
  BoutReal ionisation(BoutReal Te);
  
  //<sigma*v> [m3/s]
  BoutReal recombination(BoutReal n, BoutReal Te);
  
  // <sigma*v> [m3/s]
  BoutReal chargeExchange(BoutReal Te);
  
  // <sigma*v> [m3/s]
  BoutReal excitation(BoutReal Te);
  
private:
  
};

/*!
 * Hydrogen rates, fitted by Hannah Willett May 2015
 * University of York
 */
class UpdatedRadiatedPower : public RadiatedPower {
public:
  BoutReal power(BoutReal Te, BoutReal ne, BoutReal ni);  

  // Ionisation rate coefficient <sigma*v> [m3/s]
  BoutReal ionisation(BoutReal Ne, BoutReal T); 
  BoutReal ionisation_old(BoutReal T);
  
  // Recombination rate coefficient <sigma*v> [m3/s]
  BoutReal recombination(BoutReal n, BoutReal T);
  
  // Charge exchange rate coefficient <sigma*v> [m3/s]
  BoutReal chargeExchange(BoutReal Te);
  
  BoutReal excitation(BoutReal Ne, BoutReal Te);
  BoutReal excitation_old(BoutReal Te);
  
  // Yulin's neutral excited state population coefficients [Nn(H(n=x)) / Nn(H)]
  // Ratios of excited state to ground state populations.
  // NOTE THAT TEMPERATURE IS FIRST AND DENSITY SECOND OUTPUT
  BoutReal Channel_H_2_amjuel(BoutReal T,BoutReal Ne);
  BoutReal Channel_H_3_amjuel(BoutReal T,BoutReal Ne);
  BoutReal Channel_H_4_amjuel(BoutReal T,BoutReal Ne);
  BoutReal Channel_H_5_amjuel(BoutReal T,BoutReal Ne);
  BoutReal Channel_H_6_amjuel(BoutReal T,BoutReal Ne);
  
private:
  
};


/// Carbon in coronal equilibrium 
/// From I.H.Hutchinson Nucl. Fusion 34 (10) 1337 - 1348 (1994)
class HutchinsonCarbonRadiation : public RadiatedPower {
  BoutReal power(BoutReal Te, BoutReal ne, BoutReal ni) {
    return ne * ni * 2e-31*pow(Te/10., 3) / (1. + pow(Te/10., 4.5));
  }
};

/*
class PostJensen : public RadiatedPower {
public:
  BoutReal power(BoutReal Te, BoutReal ne, BoutReal ni) {
    if( (Te < Tmin) || (Te > Tmax) )
      return 0.0;
    
    return 0.0;
  }
protected:
  BoutReal Tmin, Tmax;
  BoutReal A[6];
  
  struct PJ_Data {
    const char* label; // Short name
    const char* name;  // Long name
    BoutReal Tmin, Tmax;
    BoutReal data[6];
  };
  
  static PJ_Data power_data[] = {
    {"C", "Carbon"},
    {0}
  };
};
*/

#endif // __RADIATION_H__
