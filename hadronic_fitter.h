/*

C++ reimplementation of rocfit_kin.py's hadronic fitter using TMinuit as the underlying fitter.

Lots more documentation is there.

 */

#ifndef rocfit_hadronic_fitter_h
#define rocfit_hadronic_fitter_h

#include "TLorentzVector.h"
#include "TF1.h"
#include "TH1D.h"
#include "TMinuit.h"


class hadronic_fitter {

 public:
  
  int _dbg, _max_calls, _minuit_tolerance, _minuit_print_level;
  double _effs[3], _probs[3], _cums[3], _res[2]; // debug tracking
  bool _use_improve, _pre_simplex;

  bool _converged;

  hadronic_fitter(); // only the default constructor can be called without invoking ACLiC outside python

  void setup( int prob_track_level = 0, double top_mass = 172.0, double top_width = 0.65, double W_mass = 80.4, double W_width = 1.0425 ); 
  bool fit( const double* P, const double* Q, const double* H, const TH1& PITF, const TH1& QITF, const TH1& HITF,
	    const TF1& Peff, const TF1& Qeff, const TF1& Heff );
  bool fit( const TLorentzVector& P, const TLorentzVector& Q, const TLorentzVector& H, const TH1& PITF, const TH1& QITF, const TH1& HITF,
	    const TF1& Peff, const TF1& Qeff, const TF1& Heff );

  // Get (pre calculated) outputs
  double ratio( int index ) const { return ( index >= 0 && index < 3 ) ? _params[ index ] : 0.0 ; }
  int n_calls() const { return _minuit.fNfcn; }
  double MLL() const { return _MLL; }

  TLorentzVector fitT() const { return _T; }
  TLorentzVector fitW() const { return _W; }

  // This is the main workhorse function - so it does not check for errors and happily crashes to save time!
  double calc_MLL( double* par, bool track_prob = false ); // public since it is called by TMinuit

  void set_print_level( int lvl ) { _minuit.SetPrintLevel( lvl );}

  void use_migrad() { _minimization_command = "MIGRAD"; _pre_simplex = false; }
  void use_minimize() { _minimization_command = "MINI"; _pre_simplex = false; }
  void use_simplex() { _minimization_command = "SIMPLEX"; _pre_simplex = false; }
  void use_simplex_first() { _minimization_command = "MIGRAD"; _pre_simplex = true; }

  // ========================== public only for debug purposes ============================
  const TMinuit& minuit () const { return _minuit; }

  void fitted_val_err( int ip, double& val, double& err ) const;
  double fitted_val( int ip ) const;
  double fitted_err( int ip ) const;

 private:

  // setup
  int _prob_track_level;
  double _top_mass, _top_iwidth, _W_mass, _W_iwidth;
  TString _minimization_command;

  // pointers to inputs
  const TF1 *_Peff, *_Qeff, *_Heff;
  const TH1 *_PITF, *_QITF, *_HITF;

  TMinuit _minuit;

  TLorentzVector _genP, _genQ, _genH, _T, _W;
  TLorentzVector _obsP, _obsQ, _obsH;

  double _params[ 3 ];

  double _MLL;

  void update_gen( double current_ratios[] );

};

#endif
