#include "simple_reco.h"
#include <TMath.h>

#include <iostream>

using std::cout;
using std::cerr;
using std::endl;

simple_reco::simple_reco( int debug_level )
{
  _dbg = debug_level;
  _det = _scale = _mtt = _alt_mtt = -6.66;
}


bool simple_reco::reco_mtt( double mex, double mey, const TLorentzVector& lvl, 
			    const TLorentzVector& lvj1, const TLorentzVector& lvj2, const TLorentzVector& lvj3, const TLorentzVector& lvj4 )
{
  _mtt = _alt_mtt = -6.66;

  bool nu_OK = reco_nu( mex, mey, lvl );
  if( ! nu_OK ) return false;

  TLorentzVector lvW;
  lvW = _nu + lvl;
  if( _dbg > 29 ) cout<<"DBG30 simple_reco::reco_mtt  lvW.M(): "<<lvW.M()<<endl;
  TLorentzVector lvTt( lvW );
  lvTt += lvj1 + lvj2 + lvj3 + lvj4;
  _mtt = lvTt.M();

  TLorentzVector lvoW;
  lvoW = _alt_nu + lvl;
  if( _dbg > 29 ) cout<<"DBG30 simple_reco::reco_mtt  lvoW.M(): "<<lvoW.M()<<endl;
  TLorentzVector lvoTt( lvoW );
  lvoTt += lvj1 + lvj2 + lvj3 + lvj4;
  _alt_mtt = lvoTt.M();
  return true;
}

bool simple_reco::reco_nu( double mex, double mey, const TLorentzVector& lvl )
{
  _det = _scale = -6.66;
  // Gave up on checking input lepton masses, since their accuracy in TopTrees is pathetic

  _nu.SetXYZM( mex, mey, 0., 0. );

  double qz = 0;
  _det = solve_W( lvl, _nu, qz );

  // Find neutrino 4-vector
  if( _det >= 0 ) {

    _nu.SetXYZM( mex, mey, qz, 0. );
    
  } else { // negative determinant - scaling the MET to satisfy the mass constraints (uses the mutables set by solve_W beforehand)

    double met_a = 4 * ( _R * _R  + _a * _met_sq );
    double met_b = 4 * _R * _D;
    double met_c = _D * _D;

    double met_det = met_b * met_b - 4 * met_a * met_c;
    if( _dbg > 29 ) cout<<"DBG30 simple_reco::reco_nu  met_a: "<<met_a<<", b: "<<met_b<<", c: "<<met_c<<", det: "<<met_det<<endl;

    if( met_det < 0 ) { 
      if( met_det < 1E-10 ) {cerr<<"ERROR! met_det<0 - can't solve for scaled MET."<<endl; return false; }
      met_det = 0;
    }
        
    if( met_det == 0 ) {
      _scale = - met_b / ( 2 * met_a );
    } else {
        
      double met1 = ( - met_b + TMath::Sqrt( met_det ) ) / ( 2 * met_a );
      double met2 = ( - met_b - TMath::Sqrt( met_det ) ) / ( 2 * met_a );

      if( met1 <= 0 && met2 <= 0 ) { cerr<<"ERROR! can't scale MET - both met1 & 2 are negative"<<endl; return false; }
      bool use1 = ( met2 <= 0 ) || ( met1 > 0 && ( TMath::Abs( TMath::Log( met2 ) ) > TMath::Abs( TMath::Log( met1 ) ) ) );
      if( _dbg > 29 ) cout<<"DBG30 simple_reco::reco_nu met1: "<<met1<<", 2: "<<met2<<" -> use1? "<<use1<<endl;
      
      _scale = use1 ? met1 : met2;
      if( _scale <= 0 ) { cerr<<"BUG! simple_reco chose a negative MET scale met1: "<<met1<<" met2: "<<met2<<endl; return false; }
    }
    _nu.SetXYZM( mex * _scale, mey * _scale, qz, 0. );
        
  } // if MET scaling needed

  // Prepare alternative reconstruction
  _alt_nu = _nu;
  if( _other_qz != qz ) _alt_nu.SetXYZM( mex, mey, _other_qz, 0. );
  if( _other_qz != qz && _scale >= 0 ) { cerr<<"BUG! two qz solution and MET rescaled?!"<<endl; return false; }

  if( _scale < 0 ) _scale = 1.0;

  return true;
}

double simple_reco::solve_W( const TLorentzVector& lvl, const TLorentzVector& lvnu, double& qz )
{
  if( _dbg > 105 ) { cout<<"DBG106 simple_reco::solve_W  lep M: "<<lvl.M()<<" 4vec: "; lvl.Print(); cout<<"nu M: "<<lvnu.M()<<" 4vec: "; lvnu.Print(); }
  // solve for leptonic W
  const double m_W = 80.399; // GeV

  double mx = lvnu.Px();
  double my = lvnu.Py();
  double lx = lvl.Px();
  double ly = lvl.Py();
  double lz = lvl.Pz();
  double lE = lvl.E();
  double lm2 = lvl.M2();

  // - quadratic equation for q_Z (neutrino momentum component)
  _El_sq = lE * lE;
  _met_sq = mx * mx + my * my;
  _R = lx * mx + ly * my;
  _D = m_W * m_W - lm2;
  _a = lz * lz - _El_sq;
  double qz_b = lz * ( _D + 2 * _R );
  double qz_c = _D * _D / 4.0 + _R * _R + _D * _R - _El_sq * _met_sq;

  double qz_det = qz_b * qz_b - 4 * _a * qz_c;
  if( _dbg > 105 ) { cout<<"DBG106 simple_reco::solve_W  _El_sq: "<<_El_sq<<", _met_sq: "<<_met_sq<<", _R: "<<_R<<", _D: "<<_D
                         <<"; _a: "<<_a<<", b: "<<qz_b<<" c: "<<qz_c<<" --> det: "<<qz_det<<endl; }

  if( qz_det > 0 ) {
    
    double qz1 = ( - qz_b + TMath::Sqrt( qz_det ) ) / ( 2 * _a );
    double qz2 = ( - qz_b - TMath::Sqrt( qz_det ) ) / ( 2 * _a );
    
    if( TMath::Abs(qz1) > TMath::Abs(qz2 ) ) qz = qz2;
    else                                     qz = qz1;

    _other_qz = qz1 + qz2 - qz;
    
  } else {
    
    qz = - qz_b / ( 2 * _a );
    _other_qz = qz;

  }

  return qz_det;
}

const TLorentzVector& simple_reco::nu() const
{
  return _nu;
}

const TLorentzVector& simple_reco::alt_nu() const
{
  return _alt_nu;
}

bool simple_reco::two_solutions() const
{
  return _det > 0;
}

double simple_reco::det() const
{
  return _det;
}

double simple_reco::scale() const
{
  return _scale;
}

double simple_reco::mtt() const
{
  return _mtt;
}

double simple_reco::alt_mtt() const
{
  return _alt_mtt;
}
