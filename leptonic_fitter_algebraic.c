/* 
   See header file for documentation
 */

#include "leptonic_fitter_algebraic.h"

#include "TLorentzVector.h"
#include "TH2D.h"
#include "TMath.h"
#include "TMatrixTBase.h"
#include "TMatrixDEigen.h"
#include "TVectorD.h"
#include "Math/Factory.h"

#include <iostream>
#include <stdexcept>
#include <vector>
#include <cassert>

using std::cout;
using std::cerr;
using std::endl;

static leptonic_fitter_algebraic* leptonic_fitter_algebraic_object = 0; // file scope global

double leptonic_fitter_algebraic::_Qnums[9] = { 0, -1, 0,     1, 0, 0,    0, 0, 0};
double leptonic_fitter_algebraic::_unit_circle_nums[9] = { 1, 0, 0,     0, 1, 0,    0, 0, -1};

// the function interface needed for minuit and used for the more modern interface as well
//____________________________________________________________________________________________________________________________________________

double leptonic_fitter_algebraic_function( const double *par ) // file scope global
{
  return leptonic_fitter_algebraic_object->calc_MLL( par );
}
//____________________________________________________________________________________________________________________________________________

double leptonic_fitter_algebraic::solutions_MLL( double sB, double dnux, double dnuy, bool track_prob )
{
  double Beff = _Beff->Eval( TMath::Log( _genB.E() ) );

  Int_t Bbin = _BITF->GetXaxis()->FindFixBin( sB );
  double BITF = _BITF->GetBinContent( Bbin );
  double Xprob = _dnuPDF.Eval( TMath::Abs( dnux ) );
  double Yprob = _dnuPDF.Eval( TMath::Abs( dnuy ) );

  double prob = Beff * BITF * Xprob * Yprob;

  double obsMLL = ( prob <= 0 ) ? 666.0 : - TMath::Log( prob );

  // Penalties from the invariant masses are calculated, when needed, elsewhere

  if( _dbg > 99 ) {
    int threshold = track_prob ? 150 : 300;
    if( _dbg >= threshold ) cout<<"DBG"<<threshold<<" lfa::MLL of "<<sB<<" "<<dnux<<" "<<dnuy<<" -> obsMLL: "<<obsMLL<<endl;
  }

  if( track_prob ) {
    _eff = Beff;
    _probs[ 0 ] = BITF;
    _probs[ 1 ] = Xprob;
    _probs[ 2 ] = Yprob;
    _cums[0] = _BITF->Integral( 0, Bbin );
    _cums[1] = _dnuPDF.Integral( 0, TMath::Abs( dnux ) );
    _cums[2] = _dnuPDF.Integral( 0, TMath::Abs( dnuy ) );
  }
  return obsMLL;
}
//____________________________________________________________________________________________________________________________________________

void leptonic_fitter_algebraic::update_non_nu_gen( double sB )
{
  _deltaB = sB - 1;
  _genB = TMath::Max( 1E-2, sB ) * _obsB;
}
//____________________________________________________________________________________________________________________________________________

void leptonic_fitter_algebraic::update_nu_and_decay_chain( double nu_px, double nu_py, double nu_pz )
{
  TLorentzVector nu;
  nu.SetXYZM( nu_px, nu_py, nu_pz, 0. );
  update_nu_and_decay_chain( nu );
}
// _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _

void leptonic_fitter_algebraic::update_nu_and_decay_chain( const TMatrixD& nu_vec )
{
  TLorentzVector nu = lv_with_mass( nu_vec );
  update_nu_and_decay_chain( nu );
}
// _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _

void leptonic_fitter_algebraic::update_nu_and_decay_chain( const TLorentzVector& nu )
{  
  _genN = nu;
  _W = nu + _obsL;
  _T = _W + _genB;
}
//____________________________________________________________________________________________________________________________________________

void leptonic_fitter_algebraic::update_penalty( double Z2 )
{
  if( Z2 < 0 ) { // no neutrino solutions - b and lepton are outside physically allowed range
    _pen_nu_vec = TMatrixD( _Emat, TMatrixD::kMult, _tvec );
    update_nu_and_decay_chain( _pen_nu_vec );
    double diff_W = _W.M() - _W_mass;
    double diff_t = _T.M() - _top_mass;
    _penalty = _t_c0 + _t_c1 * diff_t * diff_t   +   _w_c0 + _w_c1 * diff_W * diff_W;

    if( _dbg > 109 ) cout<<"DBG110 LFA penalty nu P: "<<lv2str(_genN)<<" -> W: "<<lv2str(_W)<<" -> t: "<<lv2str(_T)
			 <<" -> pens: "<<_w_c0 + _w_c1 * diff_W * diff_W<<", "<<_t_c0 + _t_c1 * diff_t * diff_t<<" -> "<<_penalty<<endl;
  } else {
    _penalty = 0.;
  }
}
//____________________________________________________________________________________________________________________________________________

double leptonic_fitter_algebraic::calc_MLL( const double* par, bool track_prob )
{
  double sB = par[0];
  if( _dbg > 219 ) cout<<"DBG220 lfa::calc for: "<<sB<<endl;
  update_non_nu_gen( sB );
  if( _dbg > 221 ) cout<<"DBG222 lfa::calc _genB: "<<lv2str(_genB)<<endl;

  double x1 = _x1_0 + 0.5*(_TmW_m2 - _b_m2*sB*sB) / (_denom*sB);
  double Z2 = _Z2_0 - 2*_x0*x1;
  double Z = (Z2 >= 0) ? TMath::Sqrt( Z2 ) : 0;
  if( _dbg > 224 ) cout<<"DBG225 lfa::calc Z2: "<<Z2<<endl;

  double entries[9] = { -Z*_y1/_x0,  0,  x1-_lep_p,     Z, 0, _y1,       0, Z, 0 };
  _Ebl.Use( 3, 3, entries );
  _Emat = _R_T;
  _Emat *= _Ebl;
    
  if( _dbg > 220 && Z2 >= 0 && _lep_p < 2 * _top_mass ) { // very energetic leptons can cause inaccuracies
    for( double iphi = 0; iphi < 4; ++iphi ) {
      double test_phi = iphi*TMath::Pi()/2;
      double tnums[3] = { TMath::Cos( test_phi ), TMath::Sin( test_phi ), 1. };
      TMatrixD test_tvec( 3, 1, tnums );
      TMatrixD test_nu( _Emat, TMatrixD::kMult, test_tvec );
      update_nu_and_decay_chain( test_nu );
      cout<<"DBG221 lfa::calc test_phi: "<<test_phi<<", test_nu: "<<lv2str(_genN)<<", -> W: "<<lv2str(_W)<<", t: "<<lv2str(_T)<<endl;
      assert( TMath::Abs( _W.M() - _W_mass ) < 1 && TMath::Abs( _T.M() - _top_mass ) < 1 );
    }
  }
  _dnu_mat = _Nu;
  _dnu_mat -= _Emat;

  bool tvec_set = false;
  if( Z > 0.002 ) { // ellipse is large enough to find the best point on it.
    _M_form.Transpose( _dnu_mat );
    _M_form *= _invFlatVar;
    _M_form *= _dnu_mat;
    if( _dbg > 221 ) cout<<"DBG222 lfa::calc _M_form: "<<m2str( _M_form )<<endl;

    // find the closest tvec
    std::vector< TMatrixD > extrimal_tvecs;
    _N = _M_form;
    _N *= _Q;
    _NT.Transpose( _N );
    _P = _N;
    _P += _NT;

    TMatrixD UP( _unit_circle, TMatrixD::kMult, _P );
    real_eigenvalues( UP );
    if( _dbg > 221 ) {
      cout<<"DBG222 lfa::calc _P: "<<m2str( _P )<<" -> the following "<<_real_eigs.size()<<" evals: ";
      for( unsigned int ire=0; ire < _real_eigs.size(); ++ire ) cout<<_real_eigs[ire]<<", ";
      cout<<endl;
    }
    if( _real_eigs.size() == 0 ) {
      cerr<<"ERROR: found no real eigenvectors of UP, which makes no sense."<<endl
	  <<"       Z: "<<Z<<" UP: "<<m2str( UP )<<endl;
      throw std::runtime_error( "leptonic_fitter_algebraic has numerical problems - found no real eigenvectors of UP.");
    }

    // find the degenrate matrix, if possible, choose one with intersecting lines
    double c22 = 1;
    for( unsigned int ire=0; ire < _real_eigs.size(); ++ire ) {
      _D = _unit_circle;
      _D *= - _real_eigs[ ire ];
      _D += _P;
    
      c22 = cofactor( _D, 2, 2 );
      if( _dbg > 221 ) cout<<"DBG222 lfa::calc e#: "<<ire<<" -> _D: "<<m2str( _D )<<" -> c22: "<<c22<<endl;
      if( c22 < 0 ) break;
    }

    // D is degenrate --> represents an equation which factorizes into two line equations.
    // Several cases which determine how we intersect these lines with the unit circle

    bool as_parallel_lines = c22 >= 0; // which means it's really 0, and positive value is due to numerical inaccuracies
                                       // so this is the first logical decision in the code. The actual first if is technical.
    if( _swapped ) { // will avoid D11=~0 in all but ~10-10 of the events
      static TMatrixD tmp;
      tmp.ResizeTo( _D );
      tmp =  _D;
      swap_x_y( tmp, _D );
    }
    double D11 = _D[1][1];
    if( D11 == 0 ) { // vertical lines, where y does not appear in the equations of the two lines
      cerr<<"ERROR: leptonic_fitter_algebraic does not handle the super-rare case of D00=D11==0..\n"
	  <<"ERROR: Estimated frequency is only 10^-10....\n"
	  <<"       Please use this as the first test case..."<<endl;

    } else { // normal case: D11 isn't 0, which means that y appears in the line equations

      if( as_parallel_lines ) {
	cout<<"Information: leptonic_fitter_algebraic encountered the \"0 probability\" case of parallel lines"<<endl;
	double c00 = cofactor( _D, 0, 0 );
	if( c00 > 0 ) {
	  cerr<<"ERROR: algorithm assumes for the \"0 probability\" case of positive c22, at least c00 is non-positive. D11!=0, Z: "
	      <<Z<<" c22: "<<c22<<" c00: "<<c00<<endl;
	  throw std::runtime_error( "leptonic_fitter_algebraic has numerical problems with c22 >= 0 and c00 > 0 (D11!=0)"); 
	}
	double sqrtNc00 = TMath::Sqrt( -c00 );
	double A = D11;
	double B = _D[0][1];
	double S = -B/A;
	// line's x0 is 0
	for( int isl = 0; isl < 2; ++isl ) {
	  int sign_of_line = isl ? 1 : -1;
	  double free_coeff = _D[1][2] + sign_of_line * sqrtNc00;
	  double y0 = -free_coeff / A;
	  add_intersections_with_circle( extrimal_tvecs, S, y0, _swapped );
	}
      } else { // normal case - intersecting lines

	double sqrtNc22 = TMath::Sqrt( -c22 );
	double line_y0 = cofactor( _D, 1, 2 ) / c22;
	double A = D11;
	double line_x0 = cofactor( _D, 0, 2 ) / c22;
	double B0 = _D[0][1];
	if( _dbg > 249 ) cout<<"DBG250 lfa::_calc has c22: "<<c22<<", line x0: "<<line_x0<<", y0: "<<line_y0<<", A: "<<A<<", B = "<<B0<<" +/- "<<sqrtNc22<<endl;
	for( int isl = 0; isl < 2; ++isl ) {
	  int sign_of_line = isl ? 1 : -1;
	  double B = B0 + sign_of_line * sqrtNc22;
	  double S = -B/A;
	  double y1 = line_y0 - S * line_x0;
	  if( _dbg > 299 ) cout<<"DBG300 lfa::_calc adding intersections for S: "<<S<<", y1: "<<y1<<endl;
	  add_intersections_with_circle( extrimal_tvecs, S, y1, _swapped );
	} // loop on lines
      } // if lines intersect or parallel
    } // if y appears in equations

    if( extrimal_tvecs.size() ) {
      unsigned int i_closest = 0;
      double min_distance = met_distance( extrimal_tvecs[ 0 ] );
      for( unsigned int it = 1; it < extrimal_tvecs.size(); ++it ) {
	double distance = met_distance( extrimal_tvecs[ it ] );
	if( distance < min_distance ) {
	  i_closest = it;
	  min_distance = distance;
	}
      }
      _tvec = extrimal_tvecs[ i_closest ];
      tvec_set = true;
    } else {
      cout<<"Information: leptonic_fitter_algebraic found no extrimal solutions for p_nu ellipse of size: "<<Z<<"GeV\n"
	  <<"             which is therefore neglected"<<endl;
    }
  } // if ellipse big enough to minimize on
  if( ! tvec_set ) {
    static double special_tnums[3] = { 0., 0., 1. };
    _tvec.SetMatrixArray( special_tnums );
  }
  if( _dbg > 219 ) cout<<"DBG220 chose _tvec: "<<m2str(_tvec)<<endl;

  update_penalty( Z2 );

  _dnu_vec = TMatrixD( _dnu_mat, TMatrixD::kMult, _tvec );
  const double *dnu_array = _dnu_vec.GetMatrixArray();
  _MLL = solutions_MLL( sB, dnu_array[0], dnu_array[1], track_prob ) + _penalty;

  if( _dbg > 219 ) {
    TMatrixD test_vec( _Emat, TMatrixD::kMult, _tvec );
    update_nu_and_decay_chain( test_vec );
    cout<<"DBG220 lfa::calc returning _MLL: "<<_MLL<<" <- Final_nu: "<<lv2str(_genN)<<", -> W: "<<lv2str(_W)<<", t: "<<lv2str(_T)<<" -> penalty: "<<_penalty<<endl;
    if( _penalty == 0 && _lep_p < 2 * _top_mass ) { // very energetic leptons can cause inaccuracies) 
      assert( TMath::Abs( _W.M() - _W_mass ) < 1 && TMath::Abs( _T.M() - _top_mass ) < 1 );
    }
  }
  return _MLL;
}
//____________________________________________________________________________________________________________________________________________

double leptonic_fitter_algebraic::met_distance( const TMatrixD& tvec )
{
  TMatrixD tmp( tvec, TMatrixD::kTransposeMult, _M_form );
  TMatrixD out( tmp, TMatrixD::kMult, tvec );
  return out[0][0];
}
//____________________________________________________________________________________________________________________________________________

leptonic_fitter_algebraic::leptonic_fitter_algebraic()
: _Ebl( 3, 3 ), _Emat( 3, 3 ), _dnu_mat( 3, 3 ), _M_form( 3, 3 ), _N( 3, 3 ), _NT( 3, 3 ), _P( 3, 3 ), _D( 3, 3 ), 
  _tvec( 3, 1 ), _pen_nu_vec( 3, 1 ), _nu_vec( 3, 1 ), _dnu_vec( 3, 1 ), 
  _Q( 3, 3, _Qnums ), _unit_circle( 3, 3, _unit_circle_nums ),
  _functor( &leptonic_fitter_algebraic_function, 1 ),
  _R_T( 3, 3 ), _Nu( 3, 3 ), _invFlatVar( 3, 3 )

{
  _dbg = _prob_track_level = 0;
  _minimizer_print_level = -999;
  _max_calls = 1000;
  _tolerance = 1.E-6;
  _MLL = 0.0;

  // set up the minimizer, a-la the ROOT example NumericalMinimization.C
  _mini = ROOT::Math::Factory::CreateMinimizer( "Minuit", "Simplex" );
  assert( 0 != _mini );
  _mini->SetMaxFunctionCalls( _max_calls ); // for Minuit/Minuit2 
  _mini->SetMaxIterations( _max_calls );    // for GSL 

  setup();
}
//____________________________________________________________________________________________________________________________________________

 void leptonic_fitter_algebraic::setup( int prob_track_level, double top_mass, double top_width, double W_mass, double W_width )
{
  _prob_track_level = prob_track_level;
  _top_mass = top_mass;
  _W_mass = W_mass;
  _W_m2 = W_mass * W_mass;
  _TmW_m2 = top_mass * top_mass - _W_m2;
  double tw = ( top_width > 0 ) ? top_width : 0.65 ;
  double ww = ( W_width > 0 ) ? W_width : 1.0425 ;
  _t_c0 = TMath::Log( tw * TMath::Sqrt( TMath::TwoPi() ) );
  _t_c1 = 0.5 / ( tw * tw );
  _w_c0 = TMath::Log( ww * TMath::Sqrt( TMath::TwoPi() ) );
  _w_c1 = 0.5 / ( ww * ww );
  _converged = _swapped = false;
}
//____________________________________________________________________________________________________________________________________________

bool leptonic_fitter_algebraic::fit( const double* B, const TH1& BITF, const TF1& Beff, 
				     const double* lep,
				     double MEX, double MEY, const TF1& dnuPDF )
{
  return fit( TLorentzVector( B ), BITF, Beff, TLorentzVector( lep ), MEX, MEY, dnuPDF );
}
//____________________________________________________________________________________________________________________________________________

bool leptonic_fitter_algebraic::fit( const TLorentzVector& B, const TH1& BITF, const TF1& Beff, 
				     const TLorentzVector& lep, 
				     double MEX, double MEY, const TF1& dnuPDF )
{
  if( _dbg > 19 ) cout<<"DBG20 Entered leptonic_fitter_algebraic::fit with B mass: "<<B.M()<<", l_m:"<<lep.M()<<", MET: "<<MEX<<" "<<MEY<<endl;
  if( B.M() <= 0 ) throw std::runtime_error( "leptonic_fitter_algebraic was given a b-jet with an illegal (non-positive) mass!"); 
  if( lep.M() < 0 ) throw std::runtime_error( "leptonic_fitter_algebraic was given a lepton with an illegal (negative) mass!"); 
  _converged = _swapped = false;
  _obsB = B;
  _obsL = lep;

  _BITF = &BITF;
  _Beff = &Beff;
  _dnuPDF = dnuPDF;

  _b_m2 = B.M2();

  double lep_b_angle = lep.Angle( B.Vect() );
  double cos_lep_b = TMath::Cos( lep_b_angle );
  double sin_lep_b = TMath::Sin( lep_b_angle );
  double b_p = B.P();
  double b_e = B.E();
  _denom = b_e - cos_lep_b * b_p;
  
  _lep_p = lep.P();
  _x0 = - _W_m2 / ( 2 * _lep_p );
  _y1 = - sin_lep_b * _x0 * b_p / _denom;
  _x1_0 = _x0 * b_e / _denom  -  _y1*_y1 / _x0;
  _Z2_0 = _x0*_x0 - _W_m2 - _y1*_y1;
  if( _dbg > 219 ) cout<<"DBG220 lfa updated lepton with: "<<lv2str( lep )<<" -> x0:"<<_x0<<", y1: "<<_y1<<", x1_0: "<<_x1_0<<", Z2_0: "<<_Z2_0<<endl;

  static double bnums[3];
  bnums[0] = B.X();
  bnums[1] = B.Y();
  bnums[2] = B.Z();
  TMatrixD bXYZ( 3, 1, bnums );
  _R_T = rotation( 2, lep.Phi() ); // R_z^T
  _R_T *= rotation( 1, lep.Theta() - 0.5*TMath::Pi() ); // R_z^T R_y^T
  TMatrixD rotation_vect( _R_T, TMatrixD::kTransposeMult, bXYZ ); // R_y R_z
  double* rotation_array = rotation_vect.GetMatrixArray();
  double phi_x = - TMath::ATan2( rotation_array[2], rotation_array[1] );
  if( _dbg > 99 ) cout<<"DBG100 lfa x rotation vector is:"<<rotation_array[0]<<" "<<rotation_array[1]<<" "<<rotation_array[2]<<" -> phi_x:"<<phi_x<<endl;
  _R_T *= rotation( 0, - phi_x ); // R_z^T R_y^T R_x^T

  // set up _Nu's non-zero elements so that \vec{nu} = Nu \vec{t} for any \vec{t} (since only t's 3nd component is used, and its always 1).
  _Nu[0][2] = MEX;
  _Nu[1][2] = MEY;

  double iVarMET = TMath::Power( TMath::Max( 1., dnuPDF.GetHistogram()->GetRMS() ), -2 );
  _invFlatVar[0][0] = _invFlatVar[1][1] = iVarMET; // set up the chi^2 distance with the right order of magnitude (generalizes to rotated covariance matrix)
  if( _dbg > 209 ) cout<<"DBG210 lfa "<<dnuPDF.GetName()<<" --> iVarMET:"<<iVarMET<<endl;

  // (re)define fit parameter, so all fits start off on an equal footing
  _mini->SetPrintLevel( _minimizer_print_level );
  _mini->Clear();
  _mini->SetFunction( _functor );
  leptonic_fitter_algebraic_object = this; // set the function in the functor pointing back to this object. Doubtfull that all this redirection is needed...
  _mini->SetTolerance( _tolerance );
  bool OK = _mini->SetLimitedVariable( 0, "sB", 1.0, 0.4, 0.1, 6.0 );
  //bool OK = _mini->SetVariable( 0, "sB", 1.0, 0.4 );
  if( ! OK ) {cerr<<"minimizer (@lfa) failed to SetVariable."<<endl; return false;}

  // define 1 sigma in terms of the function
  _mini->SetErrorDef( 0.5 ); // since this is a likelihood fit

  // do the minimization
  OK = _mini->Minimize(); 
  if( _dbg > 19 && ( ! OK || _dbg > 59 ) ) cout<<"DBG INFO: initial fit @lfa returned OK: "<<OK<<", has status: "<<_mini->Status()<<endl;

  _converged = OK; // use status somehow? depends on fitter?

  // read parameters
  const double *xs = _mini->X();
  for( int ip = 0; ip < 1; ++ip ) _params[ ip ] = xs[ ip ];

  // return all intermediate results to the minimum, in particular, the discriminant
  calc_MLL( _params, true );
  TMatrixD nu_vec( _Emat, TMatrixD::kMult, _tvec );
  update_nu_and_decay_chain( nu_vec );
  if( _dbg > 203 ) cout<<"DBG204 lfa finalized _genN: "<<lv2str(_genN)<<", _W: "<<lv2str(_W)<<", & _t: "<<lv2str(_T)<<endl;

  _MLL = _mini->MinValue();
  return true;
} 
//___________________________________________________________________________________________________________________________________________

double leptonic_fitter_algebraic::likeliest_scale( const TH1& ITF )
{
  Int_t ibin = ITF.GetMaximumBin();
  return ITF.GetBinCenter( ibin );
}
//___________________________________________________________________________________________________________________________________________
 
char* leptonic_fitter_algebraic::lv2str( const TLorentzVector& lv )
{
  return Form( "%f,%f,%f,%f", lv.X(), lv.Y(), lv.Z(), lv.M() );
}
//___________________________________________________________________________________________________________________________________________

char* leptonic_fitter_algebraic::m2str( const TMatrixD& mat )
{
  Int_t ncol = mat.GetNcols();
  Int_t nrow = mat.GetNrows();
  assert( ncol < 5 && nrow < 5 );
  TString out;
  for( int ic=0; ic < ncol; ++ic ) {
    for( int ir=0; ir < nrow; ++ir ) {
      out += Form( "%7g", mat[ir][ic] );
      if( ir < nrow - 1 ) out += ", ";
    }
    if( ic < ncol - 1 ) out += "  ;  ";
  }
  return Form( "%s", out.Data() );
}
//___________________________________________________________________________________________________________________________________________

void leptonic_fitter_algebraic::real_eigenvalues( const TMatrixD& mat )
{
  Int_t ncol = mat.GetNcols();
  assert( ncol < 10 && mat.GetNrows() == ncol );

  _real_eigs.clear();

  TMatrixDEigen eigen( mat );
  const TVectorD& e_res = eigen.GetEigenValuesRe();
  const TVectorD& e_ims = eigen.GetEigenValuesIm();
  for( int ir=0; ir < ncol; ++ir ) {
    if( e_ims[ ir ] == 0 ) _real_eigs.push_back( e_res[ ir ] );
  }
}
//___________________________________________________________________________________________________________________________________________

void leptonic_fitter_algebraic::add_intersections_with_circle( std::vector< TMatrixD > &vecs, double S, double y1, bool swap )
{
  static double line[3] = { 0, 0, 1 };
  double disc = 1 + S*S - y1*y1;
  if( disc < 0 ) return;
  double sqrtDisc = TMath::Sqrt( disc );
  double denom = 1+S*S;
  for( int ix = 0; ix < 2; ++ix ) {
    int disc_sign = ix ? 1 : -1;
    double sol_x = ( -S*y1 + disc_sign*sqrtDisc ) / denom;
    line[  swap ] = sol_x;
    line[ !swap ] = y1 + S*sol_x;
    vecs.push_back( TMatrixD( 3, 1, line ) );
  }
}
//___________________________________________________________________________________________________________________________________________

double leptonic_fitter_algebraic::cofactor( const TMatrixD& mat, int i, int j )
{
  assert( mat.GetNcols() == 3 && mat.GetNrows() == 3 );
  assert( i >= 0 && i < 3 && j >= 0 && j < 3 );

  int i0 = ( i      ) ? 0 : 1;
  int i1 = ( i == 2 ) ? 1 : 2;
  int j0 = ( j      ) ? 0 : 1;
  int j1 = ( j == 2 ) ? 1 : 2;

  int sign = ( ( i+j ) % 2 ) ? -1 : 1;

  return sign * ( mat[i0][j0] * mat[i1][j1] - mat[i0][j1] * mat[i1][j0] );
}
//___________________________________________________________________________________________________________________________________________

TLorentzVector leptonic_fitter_algebraic::lv_with_mass( const TMatrixD& mat, double mass )
{
  assert( mat.GetNcols() == 1 || mat.GetNrows() == 1 );
  TLorentzVector lv;
  const double *array = mat.GetMatrixArray();
  lv.SetXYZM( array[0], array[1], array[2], mass );
  return lv;
}
//___________________________________________________________________________________________________________________________________________

TMatrixD leptonic_fitter_algebraic::rotation( int axis, double cos, double sin )
{
  assert( axis >= 0 );
  static double elem[9];
  for( int ix=0; ix < 3; ++ix ) {
    for( int iy=0; iy < 3; ++iy ) {
      elem[ 3*ix + iy ] = ( ix == iy ) ? cos : 0;
    }
  }
  elem[ 3*((axis+1)%3)+(axis+2)%3 ] = -sin;
  elem[ 3*((axis+2)%3)+(axis+1)%3 ] = sin;
  elem[ 4*(axis%3) ] = 1.;
  TMatrixD out;
  out.Use( 3, 3, elem );
  return out;
}

//___________________________________________________________________________________________________________________________________________

void leptonic_fitter_algebraic::swap_x_y( const TMatrixD& in, TMatrixD& out )
{
  if( in.GetNoElements() != 9 ) 
    throw std::runtime_error( "ERROR! leptonic_fitter_algebraic::swap_x_y handles only 3 by 3 matrices");

  if( &in == &out )
    throw std::runtime_error( "ERROR! leptonic_fitter_algebraic::swap_x_y needs different in and out objects");

  static int remap[ 9 ]={ 4, 3, 5,   1, 0, 2,  7, 6, 8 };
  static double vec_out[9];

  const double *vec_in = in.GetMatrixArray();
  for( int ic=0;ic<9;++ic ) vec_out[ ic ] = vec_in[ remap[ ic ] ];
  out.Use( 3, 3, vec_out );
}
