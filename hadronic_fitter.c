/* 
   See header file for documentation
 */

#include "hadronic_fitter.h"

#include <TLorentzVector.h>
#include <TH2D.h>
#include <TMath.h>

#include <iostream>
#include <stdexcept>

using std::cout;
using std::cerr;
using std::endl;

static hadronic_fitter* hadronic_fitter_object = 0; // file scope global

// the function interface needed for minuit
//____________________________________________________________________________________________________________________________________________

void hadronic_fitter_function( int &, double *, double &f, double *par, int ) // file scope global
{
  f = hadronic_fitter_object->calc_MLL( par );
}

//____________________________________________________________________________________________________________________________________________

double hadronic_fitter::calc_MLL( double* par, bool track_prob )
{
  update_gen( par );

  double Peff = _Peff->Eval( TMath::Log( _genP.E() ) );
  double Qeff = _Qeff->Eval( TMath::Log( _genQ.E() ) );
  double Heff = _Heff->Eval( TMath::Log( _genH.E() ) );
  double eff = Peff * Qeff * Heff;

  Int_t Pbin = _PITF->GetXaxis()->FindFixBin( par[0] );
  Int_t Qbin = _QITF->GetXaxis()->FindFixBin( par[1] );
  Int_t Hbin = _HITF->GetXaxis()->FindFixBin( par[2] );
  
  double PITF = _PITF->GetBinContent( Pbin );
  double QITF = _QITF->GetBinContent( Qbin );
  double HITF = _HITF->GetBinContent( Hbin );
  double prob = PITF * QITF * HITF * eff;

  double obsMLL = ( prob <= 0 ) ? 666.0 : - TMath::Log( prob );

  double W_residual = ( _W.M() - _W_mass ) * _W_iwidth;
  double top_residual = ( _T.M() - _top_mass ) * _top_iwidth;
  double consMLL = 0.5 * ( W_residual * W_residual + top_residual * top_residual );
  if( _dbg > 199 ) cout<<"DBG200 hf::MLL of"<<par[0]<<" "<<par[1]<<" "<<par[2]<<" -> obsMLL: "<<obsMLL<<", consMLL:"<<consMLL<<" -> "
		       <<obsMLL+consMLL<<endl;

  if( track_prob ) {
    _effs[ 0 ] = Peff;
    _effs[ 1 ] = Qeff;
    _effs[ 2 ] = Heff;
    _probs[ 0 ] = PITF;
    _probs[ 1 ] = QITF;
    _probs[ 2 ] = HITF;
    _cums[0] = _PITF->Integral( 0, Pbin );
    _cums[1] = _QITF->Integral( 0, Qbin );
    _cums[2] = _HITF->Integral( 0, Hbin );
    _res[0] = W_residual;
    _res[1] = top_residual;
  }
  return obsMLL + consMLL;
}
//____________________________________________________________________________________________________________________________________________

void hadronic_fitter::update_gen( double current_ratios[] )
{
  double brats[3];
  for( int ip = 0; ip < 3; ++ip ) brats[ ip ] = TMath::Max( 1E-2, current_ratios[ ip ] );
  _genP = brats[0] * _obsP;
  _genQ = brats[1] * _obsQ;
  _genH = brats[2] * _obsH;

  _W = _genP + _genQ;
  _T = _W + _genH;
}
//____________________________________________________________________________________________________________________________________________

hadronic_fitter::hadronic_fitter()
: _minuit( 3 )
{
  _minuit_print_level = -999;
  _dbg = 0;
  _max_calls = 1000; // actual max calls three times larger due to MIGRAD,HESSE,MIGRAD pattern
  _minuit_tolerance = 1.E-6;
  _MLL = 0.0;
  _use_improve = true;
  _pre_simplex = false;
  use_minimize();
  setup();
}
//____________________________________________________________________________________________________________________________________________

void hadronic_fitter::setup( int prob_track_level, double top_mass, double top_width, double W_mass, double W_width )
{
  _prob_track_level = prob_track_level;
  _top_mass = top_mass;
  _top_iwidth = ( top_width > 0 ) ? ( 1 / top_width ) : 1.0 ;
  _W_mass = W_mass;
  _W_iwidth = ( W_width > 0 ) ? ( 1 / W_width ) : 1.0 ;
}
//____________________________________________________________________________________________________________________________________________

bool hadronic_fitter::fit( const double* P, const double* Q, const double* H, const TH1& PITF, const TH1& QITF, const TH1& HITF,
			   const TF1& Peff, const TF1& Qeff, const TF1& Heff )
{
  return fit( TLorentzVector( P ), TLorentzVector( Q ), TLorentzVector( H ), PITF, QITF, HITF, Peff, Qeff, Heff );
}
//____________________________________________________________________________________________________________________________________________

bool hadronic_fitter::fit( const TLorentzVector& P, const TLorentzVector& Q, const TLorentzVector& H, 
			   const TH1& PITF, const TH1& QITF, const TH1& HITF,
			   const TF1& Peff, const TF1& Qeff, const TF1& Heff )
{
  if( _dbg > 19 ) cout<<"DBG20 Entered hadronic_fitter::fit with P[0]: "<<std::flush<<P[0]<<endl;
  _converged = false;
  _obsP = P;
  _obsQ = Q;
  _obsH = H;

  _PITF = &PITF;
  _QITF = &QITF;
  _HITF = &HITF;
  _Peff = &Peff;
  _Qeff = &Qeff;
  _Heff = &Heff;

  // (re)define fit parameter, so all fits start off on an equal footing
  set_print_level( _minuit_print_level );
  _minuit.mncler();
  int error_code = 0;
  char letters[ 3 ] = { 'P', 'Q', 'H' };
  for( int ip = 0; ip < 3; ++ip ) {
    if( _dbg > 89 && _minuit_print_level < 0 ) cout<<"DBG90 Defining MINUIT param #"<<ip<<" which is "<<Form( "s%c", letters[ ip ] )<<endl;
    _minuit.mnparm( ip, Form( "s%c", letters[ ip ] ), 1.0, 0.1, 0.01, 100, error_code );
    if( error_code ) {cerr<<"minuit failed to define parameter #"<<ip<<", error_code: "<<error_code
			  <<endl; throw std::runtime_error( "Failed to define Minuit parameter" ); }
  }

  // set fitting function, making sure it points back to this object
  hadronic_fitter_object = this;
  _minuit.SetFCN( hadronic_fitter_function );
  if( error_code ) {cerr<<"minuit failed to SetFCN, error_code: "<<error_code<<endl; return false;}

  // define 1 sigma in terms of the function
  double arguments[ 10 ];
  arguments[ 0 ] = 0.5; // since this is a likelihood fit
  _minuit.mnexcm( "SET ERR", arguments, 1, error_code );
  if( error_code ) {cerr<<"minuit failed to define 1 sigma, error_code: "<<error_code<<endl; return false;}

  // minimization
  if( _pre_simplex ) {
    arguments[0] = _max_calls;
    arguments[1] = _minuit_tolerance;
    _minuit.mnexcm( "SIMPLEX", arguments, 2, error_code );
    if( _dbg > 29 && ( error_code || _dbg > 59 ) ) cout<<"DBG INFO: (hf, pre) simplex returned error_code: "<<error_code<<endl;
  }

  arguments[0] = _max_calls;
  arguments[1] = _minuit_tolerance;
  _minuit.mnexcm( _minimization_command, arguments, 2, error_code );
  if( _dbg > 19 && ( error_code || _dbg > 59 ) ) cout<<"DBG INFO: initial migrad returned error_code: "<<error_code<<endl;
  // the first minimization often reports an error due to difficulties determining the error matrix
  // this can be diagnosed when n_calls() < _max_calls, but otherwise it's not known whether the fitting was good enough.

  if( _use_improve ) {
    // improving the minimum
    // _minuit.mnexcm( "IMPROVE", arguments, 1, error_code ); this also calls mncuve, generating stupid error messages that can not be turned off
    _minuit.mnimpr();
    if( _minuit.fLnewmn ) { // "unexpectedly" found a new minimum
      if( _dbg > 19 ) cout<<"DBG INFO: improve found a new minimum"<<endl;
      arguments[0] = _max_calls;
      arguments[1] = _minuit_tolerance;
      _minuit.mnexcm( "MIGRAD", arguments, 2, error_code );
      if( _dbg > 19 && ( error_code || _dbg > 59 ) ) cout<<"DBG INFO: hf's post-improve migrad returned error_code: "<<error_code<<endl;
    }
  }

  // linear error estimation, just because (unfortunately) it sometimes improves the fit
  arguments[0] = _max_calls;
  _minuit.mnexcm( "HESSE" , arguments, 1, error_code);
  if( _dbg > 9 && ( error_code || _dbg > 59 ) ) cout<<"DBG INFO: hesse returned error_code: "<<error_code<<endl;

  // 2nd minization to check convergence
  _minuit.mnexcm( _minimization_command, arguments, 2, error_code );
  if( _dbg > 9 && ( ((error_code & 0x4) > 0 && _dbg > 39 ) || _dbg > 59 || (error_code - (error_code&0x4))>0 ) ) {
    cout<<"DBG INFO: final minimization returned error_code: "<<error_code<<" after "<<n_calls()<<" calls."<<endl;
  }
  _converged = ( error_code & 0x04 ) == 0;

  // read parameters
  for( int ip = 0; ip < 3; ++ip ) _params[ ip ] = fitted_val( ip );

  // return all intermediate results to the minimum (for debugging and performance studies) & ensure up to date W & top
  if( _prob_track_level > 0 ) 
    calc_MLL( _params, true );
  else 
    update_gen( _params );

  // some fit statistics output
  double fmin,fedm,errdef;
  int    npari, nparx, istat;
  _minuit.mnstat( fmin, fedm, errdef, npari, nparx, istat );
 
  _MLL = fmin;

  return true;
}
//____________________________________________________________________________________________________________________________________________

void hadronic_fitter::fitted_val_err( int ip, double& val, double& err ) const
{
  _minuit.GetParameter( ip, val, err );
}
// _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _

double hadronic_fitter::fitted_val( int ip ) const
{
  double val, err;
  fitted_val_err( ip, val, err );
  return val;
}
// _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _

double hadronic_fitter::fitted_err( int ip ) const
{
  double val, err;
  fitted_val_err( ip, val, err );
  return err;
}
