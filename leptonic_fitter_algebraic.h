/*

C++ reimplementation of algebraic_leptonic_fit.py's leptonic fitter.

Lots more documentation is there and in rocfit_kin.

 */

#ifndef rocfit_leptonic_fitter_algebraic_h
#define rocfit_leptonic_fitter_algebraic_h

#include "TLorentzVector.h"
#include "TF1.h"
#include "TH1D.h"
#include "Math/Minimizer.h"
#include "Math/Functor.h"

#include <TMatrixD.h>

class leptonic_fitter_algebraic {

 public:
  
  int _dbg, _max_calls, _minimizer_print_level;
  double _eff, _probs[3], _cums[3], _tolerance; // probability tracking (for studying detailed performance)
  bool _converged, _swapped;

  leptonic_fitter_algebraic(); // only the default constructor can be called without invoking ACLiC outside python

  void setup( int prob_track_level = 0, double top_mass = 172.0, double top_width = 0.65, double W_mass = 80.4, double W_width = 1.0425 ); 

  bool fit( const double* B, const TH1& BITF, const TF1& Beff, 
	    const double* lep, 
	    double MEX, double MEY, const TF1& dnuPDF );

  bool fit( const TLorentzVector& B, const TH1& BITF, const TF1& Beff, 
	    const TLorentzVector& lep, 
	    double MEX, double MEY, const TF1& dnuPDF );

  double likeliest_scale( const TH1& ITF );

  // Get (pre calculated) outputs
  double param_value( int index = 0 ) const { return ( index >= 0 && index < 1 ) ? _params[ index ] : 0.0 ; }
  int n_calls() const { return _mini->NCalls(); }
  double MLL() const { return _MLL; }
  double penalty() const { return _penalty; }

  TLorentzVector fitNu() const { return _genN; }
  TLorentzVector fitT() const { return _T; }
  TLorentzVector fitW() const { return _W; }

  // This are the main workhorse function - so they do not check for errors and happily crashes to save time!
  double calc_MLL( const double* par, bool track_prob = false );

  void set_print_level( int lvl ) { _mini->SetPrintLevel( lvl );}

  // static methods that are used internally and may also be of some use externally
  static char* lv2str( const TLorentzVector& lv );
  static char* m2str( const TMatrixD& mat );
  static double cofactor( const TMatrixD& D, int i, int j ); // Assumes D is 3 by 3, as needed within the class.
  static TLorentzVector lv_with_mass( const TMatrixD& mat, double mass = 0. ); // Assumes mat is at least of size 3
  inline static TMatrixD rotation( int axis, double angle ) { return rotation( axis, TMath::Cos( angle ), TMath::Sin( angle ) ); }
  static TMatrixD rotation( int axis, double cos_angle, double sin_angle ); // not sure Lazy matrices will really gain here.

  // ========================== public only for debug purposes ============================
  const ROOT::Math::Minimizer& mini() const { return *_mini; }
  //const ROOT::Math::Functor& func() const ( return _functor; }

  TMatrixD _Ebl;
  TMatrixD _Emat;
  TMatrixD _dnu_mat; // V_0 - E(lab), i.e., the matrix that gives delta nu for each t-vec solution as \vec{dnu} = M_dnu \vec{t}
  TMatrixD _M_form; 
  TMatrixD _N, _NT;
  TMatrixD _P;
  TMatrixD _D;
  std::vector< double > _real_eigs;
  TMatrixD _tvec; // the (non-geometric) vector that generates the _nu_vec. Of the form (cos t, sin t, 1)
  TMatrixD _pen_nu_vec;
  TMatrixD _nu_vec;
  TMatrixD _dnu_vec;

  // Warning: out will *not* own its content array! It uses a static array --> the next call to swap_x_y will mess it up! (Yey ROOT!)
  void swap_x_y( const TMatrixD& in, TMatrixD& out ); // only for 3*3 matrices.

 private:

  //constants
  static double _Qnums[9], _unit_circle_nums[9];
  TMatrixD _Q, _unit_circle; // CINT doesn't like static TMatrixDs?

  // setup
  int _prob_track_level;
  double _top_mass, _W_mass, _TmW_m2, _W_m2 ;
  double _t_c0, _t_c1, _w_c0, _w_c1; // free (c0) and linear (c1) coefficiencts for parabola from log(norm(mass_term)). 

  // pointers to inputs
  const TF1 *_Beff;
  const TH1 *_BITF;
  // copies of inputs (since TF1::Integral isn't const-correct)
  TF1 _dnuPDF;

  ROOT::Math::Minimizer* _mini;
  ROOT::Math::Functor _functor;

  // event observables (per event setup)
  TLorentzVector _obsB, _obsL;
  double _obsX, _obsY;

  double _cos_lep_b, _sin_lep_b;
  double _lep_p;
  double _b_m2;
  double _x0, _y1;
  double _x1_0; // value of x1 for deltaB==0
  double _Z2_0; //   "   "  Z2  "      "
  double _denom;
  TMatrixD _R_T;
  TMatrixD _Nu;
  TMatrixD _invFlatVar;


  // ????? double _lepT2, _lepT, _lepPhi; // pre-calculated kinematics of the observed lepton

  // information on the current evaluation point
  TLorentzVector _genB, _genN, _T, _W; // genL == obsL, since tests without this constraint showed little correlation in sL

  double _params[ 1 ];
  double _MLL;

  double _deltaB;
  
  double _penalty;

  double solutions_MLL ( double sB, double dnux, double dnuy, bool track_prob = false ); // An underlying workhorse function, uses _genB as input
  void update_non_nu_gen( double sB ); // only for the final state objects that require it, i.e., only for the b-jet
  void update_nu_and_decay_chain( const TMatrixD& nu_vec );
  void update_nu_and_decay_chain( const TLorentzVector& nu );
  void update_nu_and_decay_chain( double nu_px, double nu_py, double nu_pz );
  void update_penalty( double Z2 );
  
  void real_eigenvalues( const TMatrixD& P ); // updates _real_eigs with the real eigenvalues of P.
  double met_distance( const TMatrixD& tvec ); // Returns the distance^2 (i.e., the chi^2) between the nu(tvec) and the measured met

  void add_intersections_with_circle( std::vector< TMatrixD > &vecs, // add them to this container
				      double S, double y1, // line is y=y1+Sx (appears in equations as y-y0=S(x-x0) and y1=y0-Sx0)
				      bool swap = false );  // x and y are swapped, unswap in solutions before booking them into vecs
};

#endif
