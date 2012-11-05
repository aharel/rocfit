#!/usr/bin/python2.7
#
# A kinematic fitter for top pair l+jet events using Type III inverse transfer functions
#
# This implementation uses C++ low-level fitters, which use TMinuit, including
# and an algebraic fitter based on Burt's new-fangled fitter for the leptonic top.
#
import math
import ConfigParser
import os
import itertools
import operator

import ROOT as R

from rocfit_util import *
import lorentzv as LV

R.gROOT.LoadMacro( 'hadronic_fitter.c+' )
R.gROOT.LoadMacro( 'leptonic_fitter_algebraic.c+' )
from ROOT import hadronic_fitter
from ROOT import leptonic_fitter_algebraic

#####################################
def all_matched( parton_matches ):
    '''Assumes caller defined parton_matches so that duplicates are impossible'''
    return not -1 in parton_matches
        
#####################################
def parton_region( eta ):
    return min( 4, int( abs( eta ) * 2 ) )

#####################################
def lne_bin( lne ):
    return max( 0, min( 15, int( ( lne - 2.5 ) * 4 ) ) )

#####################################
def met_bin( met ):
    return max( 0, min( 5, int( ( met - 20 ) /30. ) ) )

#####################################
class rocfit( object ):
    '''A kinematic fitter for top pair events in the lepton+jets channel, using low-level C++ code
    '''
    def __init__( self, config_path = 'rocfit.conf', do_hist = True, dbg = 1, hard_btags = False,
                  prob_track_level = 0 ):

        self.dbg = dbg
        self.prob_track_level = prob_track_level
        if dbg > 9:
            print 'Initializing a rocfit_algebraic fitter. do_hist:',do_hist,'hard_btags:',hard_btags
        if do_hist: # done first, before ROOT starts to change file contexts
            Nbx = Nby = 100
            lowX = 60
            highX = 360
            lowY = 0
            highY = 200
            self.h_premass = R.TH2D( 'hpremass', 'Pre-fit masses', Nbx, lowX, highX, Nby, lowY, highY )
            self.add_hist_vars( 'h_premass' )
            self.h_funcalls = R.TH1D( 'hcalls', '# of function calls', 1000, 0, 2000 )
            self.add_hist_vars( 'h_funcalls' )
            self.h_tagprob = R.TH1D( 'htagprob', 'log(tag) probability', 700, -14, 0 )
            self.add_hist_vars( 'h_tagprob' )
            self.h_lep_funcalls = R.TH1D( 'hlepcalls', '# of function calls in leptonic fitter', 1000, 0, 2000 )
            self.add_hist_vars( 'h_lep_funcalls' )
            if prob_track_level > 0 :
                for obj in ( 'P', 'Q', 'H', 'B' ) :
                    key = 'h_'+obj+'eff'
                    setattr( self, key, R.TH1D( 'h'+obj+'eff', obj+' eff', 200, 0, 1 ) )
                    self.add_hist_vars( key )
                for obj in ( 'P', 'Q', 'H', 'B', 'L', 'X', 'Y'  ) :
                    key = 'h_'+obj+'prob'
                    setattr( self, key, R.TH1D( 'h'+obj+'prob', obj+' prob', 200, 0, 1 ) )
                    self.add_hist_vars( key )
                    key = 'h_'+obj+'cum'
                    setattr( self, key, R.TH1D( 'h'+obj+'cum', obj+' cum. prob', 200, 0, 1 ) )
                    self.add_hist_vars( key )
                self.h_WH = R.TH1D( 'hWH', 'Hadronic W residual', 200, -5, 5 )
                self.add_hist_vars( 'h_WH' )
                self.h_TH = R.TH1D( 'hTH', 'Hadronic top residual', 200, -5, 5 )
                self.add_hist_vars( 'h_TH' )

        self.hard_btags = hard_btags

        config = ConfigParser.RawConfigParser()
        config.read
        config.read( config_path )
        print 'Opened config file \"'+config_path+'\" with sections:',config.sections()

        self.fcorr = config.getfloat( 'global', 'f_corr_met' )
        self.top_mass = config.getfloat( 'constants', 'top_mass' )
        self.W_mass = config.getfloat( 'constants', 'W_mass' )
        self.b_mass = config.getfloat( 'constants', 'b_mass' )
        self.widthT = config.getfloat( 'resolutions', 'eff_T_width' )
        self.widthW = config.getfloat( 'resolutions', 'eff_W_width' )
        
        
        try:
            btag_file_name = config.get( 'global', 'btag_file' )
            print 'Reading in btag distributions from file:',btag_file_name
            self.btag_file = R.TFile( btag_file_name ) # can not use temporary variable, since the TFile manages these histograms
            assert self.btag_file.IsOpen()

            self.b_tag_hist = self.btag_file.Get( 'hbtags' )
            self.c_tag_hist = self.btag_file.Get( 'hctags' )
            self.q_tag_hist = self.btag_file.Get( 'hqtags' )
            assert (not self.b_tag_hist is None) and (not self.c_tag_hist is None) and (not self.q_tag_hist is None)

        except ConfigParser.NoOptionError:
            if not self.hard_btags:
                print 'Config file doesn\'t contain tag rates, which are required since rocfit does not use hard b-tags'
                raise 

        ITF_file_name = config.get( 'global', 'ITF_file' )
        print 'Reading in jet bias and errors from file:',ITF_file_name
        self.ITF_file = R.TFile( ITF_file_name )
        assert self.ITF_file.IsOpen()

        self.q_ITFs = {}
        self.b_ITFs = {}
        for region in range(5):
            first_q_ITF = first_b_ITF = None
            for lne_ind in range(16):
                self.q_ITFs[ region, lne_ind ] = self.ITF_file.Get( 'hq{0}_{1}'.format( 1+region, 1+lne_ind ) )
                if first_q_ITF is None and self.q_ITFs[ region, lne_ind ]:
                    first_q_ITF = lne_ind
                self.b_ITFs[ region, lne_ind ] = self.ITF_file.Get( 'hb{0}_{1}'.format( 1+region, 1+lne_ind ) )
                if first_b_ITF is None and self.b_ITFs[ region, lne_ind ]:
                    first_b_ITF = lne_ind
            if first_q_ITF is None or first_b_ITF is None:
                print 'ERROR: no ITFs found for region',region
                raise Exception( 'BadInput', 'ITF file', 'region '+str(region) )

            if region == 0:
                self.ITF_example = self.q_ITFs[ 0, first_q_ITF ]
                assert self.ITF_example
                
            last_q_ITF = last_b_ITF = None
            for lne_ind in range(15, -1, -1):
                if last_q_ITF is None and self.q_ITFs[ region, lne_ind ]:
                    last_q_ITF = lne_ind
                if last_b_ITF is None and self.b_ITFs[ region, lne_ind ]:
                    last_b_ITF = lne_ind

            for lne_ind in range( first_q_ITF ):
                self.q_ITFs[ region, lne_ind ] = self.q_ITFs[ region, first_q_ITF ]
            for lne_ind in range( first_b_ITF ):
                self.b_ITFs[ region, lne_ind ] = self.b_ITFs[ region, first_b_ITF ]
            for lne_ind in range( 1 + last_q_ITF, 16 ):
                self.q_ITFs[ region, lne_ind ] = self.q_ITFs[ region, last_q_ITF ]
            for lne_ind in range( 1 + last_b_ITF, 16 ):
                self.b_ITFs[ region, lne_ind ] = self.b_ITFs[ region, last_b_ITF ]

        if dbg > 76 :
            print 'DBG77 e.g. E=100GeV, eta=1.1:',\
                  ' q(s=1)=',  self.q_ITFs[2,lne_bin( math.log( 100 ) )].GetBinContent( self.ITF_example.FindBin( 1 ) ),\
                  ' q(s=1.2)=',self.q_ITFs[2,lne_bin( math.log( 100 ) )].GetBinContent( self.ITF_example.FindBin( 1.2 ) ),\
                  ' b(s=1)=',  self.b_ITFs[2,lne_bin( math.log( 100 ) )].GetBinContent( self.ITF_example.FindBin( 1 ) ),\
                  ' b(s=1.2)=',self.b_ITFs[2,lne_bin( math.log( 100 ) )].GetBinContent( self.ITF_example.FindBin( 1.2 ) )

        eff_file_name = config.get( 'global', 'eff_file' )
        print 'Reading in parton efficiencies from file:',eff_file_name
        eff_file = R.TFile( eff_file_name )
        assert eff_file.IsOpen()

        self.q_effs = []
        self.b_effs = []
        for region in range(5):
            self.q_effs.append( eff_file.Get( 'qf{0}'.format( 1+region ) ) )
            self.b_effs.append( eff_file.Get( 'bf{0}'.format( 1+region ) ) )
            assert self.q_effs[ -1 ] and self.b_effs[ -1 ]

        met_file_name = config.get( 'global', 'MET_file' )
        print 'Reading in MET (delta nu) ITFs from file:',met_file_name,'  ( f_corr_met:',self.fcorr,"but isn't used in the fits )"
        met_file = R.TFile( met_file_name )
        assert met_file.IsOpen()

        self.dnu_PDFs = []
        for ib in range(6):
            self.dnu_PDFs.append( met_file.Get( 'fdnu{0}'.format( 1+ib ) ) )
            assert self.dnu_PDFs[ -1 ]

        if self.dbg > 9: print 'Preparing C++ fitter objects'
        self.hadFit = hadronic_fitter()
        self.hadFit._dbg = self.dbg
        self.hadFit.setup( self.prob_track_level, self.top_mass, self.widthT, self.W_mass, self.widthW )
        self.hadFit.use_simplex_first()
        self.lepFit = leptonic_fitter_algebraic()
        self.lepFit._dbg = self.dbg
        self.lepFit.setup( self.prob_track_level, self.top_mass, self.widthT, self.W_mass, self.widthW )
        if dbg > 10 and (dbg % 5) == 0 :
            self.hadFit._minuit_print_level = dbg / 10 - 5
            self.lepFit._minuit_print_level = dbg / 10 - 4
        if dbg > 10 and (dbg % 10) == 5 :
            self.lepFit._minuit_print_level += 3
       
    # ========================================================================================================

    def add_hist_vars( self, name ):

        horg = getattr( self, name )

        setattr( self, name+'_able', horg.Clone( horg.GetName()+'_able' ) )
        hable = getattr( self, name+'_able' )
        hable.SetTitle( horg.GetTitle()+' [matchable event]' )

        setattr( self, name+'_ing', horg.Clone( horg.GetName()+'_ing' ) )
        hing = getattr( self, name+'_ing' )
        hing.SetTitle( horg.GetTitle()+' [matching perm]' )
        
    # ========================================================================================================

    def fill_hist_vars( self, name, val_list ):
        horg = getattr( self, name )
        if horg:
            horg.Fill( *val_list )
            if self.matchable:
                hable = getattr( self, name+'_able' )
                hable.Fill( *val_list )
                if self.perm_is_correct:
                    hing = getattr( self, name+'_ing' )
                    hing.Fill( *val_list )
        
    # ========================================================================================================

    def write_hists( self ):
        '''Write all histograms.
           Note - caller must make ensure the same ROOT global TFile context used to create object'''
        nhist = 0
        for attr, value in self.__dict__.iteritems():
            cname = str( type( value ) )
            if cname.find( 'ROOT.TH1' ) >= 0 or cname.find( 'ROOT.TH2' ) >= 0:
                nhist += 1
                if self.dbg > 4:
                    if nhist < 5:
                        print 'Writing',attr,'of class',cname
                value.Write()
        print 'Wrote',nhist,'histograms'

    # ========================================================================================================

    def ITF( self, p4, is_b ):
        '''Returns a histogram of the Inverse Transfer Function of the parton.
        That is, the distribution of partonE/jetE for the given parton kinematics.
        '''
        region = parton_region( p4.Eta() )
        lnebin = lne_bin( math.log( p4.E() ) )
        if is_b:
            return self.b_ITFs[ region, lnebin ]
        else:
            return self.q_ITFs[ region, lnebin ]

    # ========================================================================================================

    def eff( self, p4, is_b ):
        '''Returns a TF1 of the efficiency for the parton.'''
        region = parton_region( p4.Eta() )
        if is_b:
            return self.b_effs[ region ]
        else:
            return self.q_effs[ region ]

    # ========================================================================================================

    def dnu_PDF( self, met ):
        return self.dnu_PDFs[ met_bin( met ) ]

    # ========================================================================================================
    def calc_tag_prob( self, iPQH, btag_vals ):
        '''Returns the probability of the btag_indices assuming iPQH is the correct hadronic assignment'''
        prob = 1.0
        if hasattr( self, 'b_tag_hist' ):
            iL = 6 - reduce( operator.add, iPQH )
            for iB in iL, iPQH[2]:
                prob *= self.b_tag_hist.GetBinContent( self.b_tag_hist.FindBin( btag_vals[ iB ] ) )
            cs_prob = 0.167 # a priori probability from PDG: half of W-->had
            sc_prob = 0.167 # a priori probability from PDG: half of W-->had
            qq_prob = 0.342 # a priori probability from PDG: W-->had - W-->cX = 0.676 - 0.334 = 
            for index, iQ in enumerate( iPQH[0:2] ):
                q_prob = self.q_tag_hist.GetBinContent( self.q_tag_hist.FindBin( btag_vals[ iQ ] ) )
                c_prob = self.q_tag_hist.GetBinContent( self.c_tag_hist.FindBin( btag_vals[ iQ ] ) )
                qq_prob *= q_prob
                cs_prob *= c_prob if index > 0 else q_prob
                sc_prob *= q_prob if index > 0 else c_prob
            prob *= qq_prob + cs_prob + sc_prob
        return prob

    # ========================================================================================================

    def fit( self, p4_lep, p4_jets, mex, mey, parton_matches, cheat = False, nofit = False,
             btag_indices = [], btag_vals = [], known_problem_with_lepton = False ):
        '''Does a kinematic fit given 4-vectors (the p4 inputs), lepton charge, met and btags.

           The "cheat" option means to consider only the correct assignment (if any).
           The parton_matches and known_problem_with_lepton are used to cheat and to create debug information

           If hard_btags constraints are used, only the btag_indices are needed,
           but the btag_vals may be useful to study tag_prob s.
           Otherwise, only the full btag_vals are needed.
           '''
        fit_allowed = not nofit
        Nj = len( p4_jets )
        assert Nj >= 4
        for ib in btag_indices:
            assert ib >= 0 and ib < Nj
        if not self.hard_btags:
            assert len( btag_vals ) == Nj
        if self.dbg > 50:
            print 'rocfit_algebraic.fit inputs are reasonable.'
        mex = R.TMath.Range( -300, 300, mex )
        mey = R.TMath.Range( -300, 300, mey )
        raw_nu = LV.new()
        raw_nu.SetXYZM( mex, mey, 0., 0. )
        met = raw_nu.Pt()

        lv_lep = LV.new()
        lv_lep.SetXYZM( p4_lep.X(), p4_lep.Y(), p4_lep.Z(), max( 0, p4_lep.M() ) )
        for tries in range( 14 ):
            if lv_lep.M() >= 0:
                break
            lv_lep.SetXYZM( p4_lep.X(), p4_lep.Y(), p4_lep.Z(), max( 1E-3 * 2**tries, p4_lep.M() ) )
            if tries > 9:
                print 'Warning: lepton mass insists on being negative. # of tries:',tries

        Ntag = len( btag_indices )
        self.matchable = all_matched( parton_matches ) and not known_problem_with_lepton
        true_Bs = parton_matches[ 0:2 ]
        true_Qs = parton_matches[ 2:4 ]

        j_as_q_ITF = [ self.ITF( p4, False ) for p4 in p4_jets ]
        j_as_b_ITF = [ self.ITF( p4, True ) for p4 in p4_jets ]
        j_as_q_eff = [ self.eff( p4, False ) for p4 in p4_jets ]
        j_as_b_eff = [ self.eff( p4, True ) for p4 in p4_jets ]

        recos = []
        indices = range( 4 )

        hf = self.hadFit
        lf = self.lepFit

        for iPQH in itertools.permutations(indices,3) : # loop over possible hadronic assignments
            if iPQH[0]>iPQH[1] : continue               # the two equivalent W->qq assignments
            iL = 6 - reduce( operator.add, iPQH )
            if self.hard_btags:
                if Ntag >= 2:
                    if iPQH[2] not in btag_indices or iL not in btag_indices: continue   # b-quarks can only be assigned to tagged jets
                if Ntag <= 2:
                    if iPQH[0] in btag_indices or iPQH[1] in btag_indices : continue # no spare tag to waste on W->q jets

            self.perm_is_correct = iPQH[0] in true_Qs  and  iPQH[1] in true_Qs  and \
                                   iPQH[2] == parton_matches[1]  and  iL == parton_matches[0]
            if cheat and not self.perm_is_correct:
                continue

            tag_prob = self.calc_tag_prob( iPQH, btag_vals )
            if self.dbg > 28:
                print 'DBG29 iPQH:',iPQH,'-> perm_is_correct:',self.perm_is_correct,' tag_prob:',tag_prob
            self.fill_hist_vars( 'h_tagprob', [ math.log(tag_prob) ] )

            lv_W_had = p4_jets[ iPQH[0] ] + p4_jets[ iPQH[1] ]
            lv_t_had = p4_jets[ iPQH[2] ] + lv_W_had
            m_W_had = lv_W_had.M()
            m_t_had = lv_t_had.M()
            
            self.fill_hist_vars( 'h_premass', [ m_t_had, m_W_had ] )

            this_reco = {'iPQH': iPQH, 'Ptags' : tag_prob }

            if fit_allowed:

                if self.dbg > 49:
                    print 'DBG50 will fit iPQH:',iPQH,'jITF0:',j_as_q_ITF[ iPQH[ 0 ] ]
                                                              
                perm_jets = [ p4_jets[ i ] for i in iPQH ]
                perm_ITFs = [ j_as_q_ITF[ iPQH[ 0 ] ], j_as_q_ITF[ iPQH[ 1 ] ], j_as_b_ITF[ iPQH[ 2 ] ] ]
                perm_effs = [ j_as_q_eff[ iPQH[ 0 ] ], j_as_q_eff[ iPQH[ 1 ] ], j_as_b_eff[ iPQH[ 2 ] ] ]

                corr_b_jet = LV.new( perm_jets[2] )
                self.correct_b_mass( corr_b_jet, perm_ITFs[ 2 ] )
                if self.dbg > 9:
                    print 'DBG10 TMP',lv_as_PM_4tup( perm_jets[2] ),' --> ',lv_as_PM_4tup( corr_b_jet )

                hf.fit( perm_jets[ 0 ], perm_jets[ 1 ], corr_b_jet, 
                        perm_ITFs[ 0 ], perm_ITFs[ 1 ], perm_ITFs[ 2 ],
                        perm_effs[ 0 ], perm_effs[ 1 ], perm_effs[ 2 ] )

                ratios = [ hf.ratio( ii ) for ii in range( 3 ) ]
                    
                delta_HT = reduce( operator.add, map( lambda p4,rat: p4*(self.fcorr * (rat-1)), perm_jets, ratios ))
                self.fill_hist_vars( 'h_funcalls', [hf.n_calls()] )

                if self.dbg > 12:
                    print 'DBG13',iPQH,'perm_jets:',[lv_as_3tup(lv) for lv in perm_jets],\
                          'had MLL:',hf.MLL(),'ratios:',ratios,' -->',delta_HT.Print()

                this_reco.update( { 'hadR' : ratios, 'hadMLL' : hf.MLL(), 'dHT' : delta_HT, 'hadNC' : hf.n_calls(),
                                    'hadW' : hf.fitW(), 'hadT' : hf.fitT(), 'hadConv' : hf._converged } )

                if self.prob_track_level > 0:
                    self.fill_hist_vars( 'h_Peff', [hf._effs[0]] )
                    self.fill_hist_vars( 'h_Qeff', [hf._effs[1]] )
                    self.fill_hist_vars( 'h_Heff', [hf._effs[2]] )
                    self.fill_hist_vars( 'h_Pprob', [hf._probs[0]] )
                    self.fill_hist_vars( 'h_Qprob', [hf._probs[1]] )
                    self.fill_hist_vars( 'h_Hprob', [hf._probs[2]] )
                    self.fill_hist_vars( 'h_Pcum', [hf._cums[0]] )
                    self.fill_hist_vars( 'h_Qcum', [hf._cums[1]] )
                    self.fill_hist_vars( 'h_Hcum', [hf._cums[2]] )
                    self.fill_hist_vars( 'h_WH', [hf._res[0]] )
                    self.fill_hist_vars( 'h_TH', [hf._res[1]] )
                    if self.prob_track_level > 2:
                        this_reco.update( { 'hadEffs' : hf._effs, 'hadProbs' : hf._probs, 'hadRes' : hf._res,
                                            'hadCums' : hf._cums } )

                corr_lep_b_jet = LV.new( p4_jets[iL] )
                self.correct_b_mass( corr_lep_b_jet,  j_as_b_ITF[iL] )
                if self.dbg > 9:
                    print 'DBG10 TMP leptonic ',lv_as_PM_4tup( p4_jets[ iL ] ),' --> ',lv_as_PM_4tup( corr_lep_b_jet )
                OK = lf.fit( corr_lep_b_jet, j_as_b_ITF[iL], j_as_b_eff[iL],
                             lv_lep,
                             mex, mey, self.dnu_PDF( met ) )
                if not OK:
                    print 'ERROR! leptonic_fitter_algebraic failed!?'
                    return

                fitNu = lf.fitNu()

                if self.dbg > 39:
                    print 'DBG40 lepFit has MLL:',lf.MLL(),'.nu:',lv_as_PM_4tup(fitNu),\
                          '.W:',lv_as_PM_4tup(lf.fitW()),'.T:',lv_as_PM_4tup(lf.fitT()),' --> penalty: ',lf.penalty()
                elif self.dbg>19:
                    print 'DBG20 lepFit return W & top masses of:',lf.fitW().M(),'&',lf.fitT().M(),' --> penalty: ',lf.penalty()
                assert lf.penalty() > 0 or lv_lep.P() > 2*self.top_mass or \
                       (abs( lf.fitW().M() - self.W_mass ) < 1. and abs( lf.fitT().M() - self.top_mass ) < 1.)

                MLL = kinMLL = hf.MLL() + lf.MLL()
                if not self.hard_btags:
                    MLL -= math.log( tag_prob )

                self.fill_hist_vars( 'h_lep_funcalls', [lf.n_calls()] )

                dnu = raw_nu - fitNu

                this_reco.update( {'nu'     : fitNu,               'lepW' : lf.fitW(),       'lepT' : lf.fitT(),
                                   'lepB'   : lf.param_value(),    'lepL' : 1.0,
                                   'lepX'   : dnu.X(),             'lepY' : dnu.Y(), 
                                   'lepMLL' : lf.MLL(),          'kinMLL' : kinMLL,           'MLL' : MLL,
                                   'lepNC'  : lf.n_calls(),     'penalty' : lf.penalty(), 
                                   'lepConv': lf._converged,       'swap' : lf._swapped
                                   } )
                if self.prob_track_level > 0:
                    self.fill_hist_vars( "h_Beff", [lf._eff] )
                    self.fill_hist_vars( "h_Bprob", [lf._probs[0]] )
                    self.fill_hist_vars( "h_Bcum", [lf._cums[0]] )
                    self.fill_hist_vars( "h_Xprob", [lf._probs[1]] )
                    self.fill_hist_vars( "h_Yprob", [lf._probs[2]] )
                    self.fill_hist_vars( "h_Xcum", [lf._cums[1]] )
                    self.fill_hist_vars( "h_Ycum", [lf._cums[2]] )
                    if self.prob_track_level > 2:
                        this_reco.update( { 'lepEffs' : [lf._eff], 'lepProbs' : lf._probs, 'lepCums' : lf._cums } )
            recos.append( this_reco )
            if self.dbg>89:
                print 'DBG90 this_reco:',sorted(this_reco.items())

        if fit_allowed:
            if self.hard_btags:
                recos.sort( key=lambda assignment: assignment[ 'kinMLL' ] )
            else:
                recos.sort( key=lambda assignment: assignment[ 'MLL' ] )
        
        return recos

    # ========================================================================================================

    def correct_b_mass( self, p4, ITF ):
        scale = R.TMath.Range( 0.1, 6.0, self.lepFit.likeliest_scale( ITF ) )
        p4 *= scale
        p4.SetVectM( p4.Vect(), self.b_mass )
        p4 *= 1 / scale
        
