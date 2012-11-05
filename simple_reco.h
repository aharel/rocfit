// Routines for Simple mtt Reconstruction

#ifndef simple_reco_h
#define simple_reco_h

#include "TLorentzVector.h"

class simple_reco {

 public :

  int _dbg;

  simple_reco( int debug_level = 0 );
  
  double solve_W( const TLorentzVector& lvl, const TLorentzVector& lvnu, double& qz0 ); // returns determinant. Only x and y components of lvnu are used.

  bool reco_nu( double mex, double mey, const TLorentzVector& lvl ); // returns whether reco succeeded

  bool reco_mtt( double mex, double mey, const TLorentzVector& lvl,  
		 const TLorentzVector& lvj1, const TLorentzVector& lvj2, const TLorentzVector& lvj3, const TLorentzVector& lvj4 ); // returns whether reco succeeded

  const TLorentzVector& nu() const; // resulting neutrino 4-vec
  const TLorentzVector& alt_nu() const; // the alternative, less likely, reco. result

  bool two_solutions() const; // whether we have two qz solution (in which case, alt_nu and alt_mtt are useful)
  double det() const; // determinant from solving the pZ equation (before any rescaling)
  double scale() const; // the resulting MET scaling
  double mtt() const; // the resulting ttbar mass
  double alt_mtt() const; // the alternative, less likeliy, reco. ttbar mass

  mutable double _El_sq, _met_sq, _R, _D, _a, _other_qz; // temp cache. private except for debug purposes.

 private:

  TLorentzVector _nu, _alt_nu;
  double _det, _scale;
  double _mtt, _alt_mtt;
  
};

#endif



