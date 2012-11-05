#!/usr/bin/python2.7
#
# A partial reconstruction algorithm for 3-jet top pair l+jet events
#
# Working assumption: the missing jet is a light jet
#
# Working approximations:
# - though the missing jet may have merged into another jet (25-30% of the time), we always use the observed jets as is
# - we neglect the missing jet, though in ~20% of the time it has very high
#   |eta| causing heavy tails in our hadronic-y resolution
#
# Method at a glance:
# - simple reconstruction of neutrino from W mass constraint with MET rescaling
# - use the top mass constraint to choose the neutrino PZ solution
# - identify the correct permutation with the leptonic top mass and a partial
#   (two-body) invariant mass as a proxy for the hadronic top 
#
import operator
import math

import ROOT as R
import lorentzv as LV
import rocfit_util as U
import top_angle_util

if R.gROOT.GetClass('simple_reco') == None:
    R.gROOT.LoadMacro( 'simple_reco.c+' )
from ROOT import simple_reco

# ========================================================================================================
class partial_reco( object ):
    '''A partial reconstruction for top pair events in the lepton+jets channel with only 3 jets
    '''
    def __init__( self, dbg = 1, top_mass = 172.5, fQ = 0.17 ) :
        '''fQ is the fraction of 3 jet events that match bqq,
        out of all those that match the bqq or bbq hypotheses'''

        self.dbg = dbg
        if dbg > 9:
            print 'Initializing a partial_reco'

        self.top_mass = top_mass
        self.fQ = fQ

        btag_file_name = 'partial_reco_btag.root'
        print 'Reading in btag distributions from file:',btag_file_name
        self.btag_file = R.TFile( btag_file_name ) # can not use temporary variable, since the TFile manages these histograms
        assert self.btag_file.IsOpen()

        self.b_tag_hist = self.btag_file.Get( 'hbtags' )
        self.w_tag_hist = self.btag_file.Get( 'hwtags' )
        assert (not self.b_tag_hist is None) and (not self.w_tag_hist is None)

        mass_file_name = 'partial_reco_mass.root'
        print 'Reading in mass distributions from file:',mass_file_name
        self.mass_file = R.TFile( mass_file_name )
        assert self.mass_file.IsOpen()

        for name in [ 'h_mt_l', 'h_mt_h', 'h_mt_q', 'h_mp_hl', 'h_mp_hq', 'h_mp_lq', 'h_mp_qq' ]:
            setattr( self, name, self.mass_file.Get( name ) )
            hist = getattr( self, name )
            assert not hist is None
            if hist.GetMinimum() < 0 :
                raise Exception( 'BadInput', name, 'Input mass distribution histogram has negative entries' )

        self.SR = simple_reco( dbg / 2 )

    # ========================================================================================================
    def calc_tag_prob_H( self, iW, btag_vals ):
        '''Returns the 3 jet b-tagging probability assuming bbq and iW, the index of the jet from the hadronic W decay'''
        assert iW >= 0 and iW <= 2
        prob = U.hist_value_at( self.w_tag_hist, btag_vals[ iW ] )
        for iB in range( 3 ):
            if iB != iW :
                prob *= U.hist_value_at( self.b_tag_hist, btag_vals[ iB ] )
        return prob

    # ========================================================================================================
    def calc_tag_prob_Q( self, iL, btag_vals ):
        '''Returns the 3 jet b-tagging probability assuming bqq and iL, the index of the jet from the leptonic b'''
        assert iL >= 0 and iL <= 2
        prob = U.hist_value_at( self.b_tag_hist, btag_vals[ iL ] )
        for iW in range( 3 ):
            if iW != iL :
                prob *= U.hist_value_at( self.w_tag_hist, btag_vals[ iW ] )
        return prob

    # ========================================================================================================
    @staticmethod
    def averagables_array( reco, lv_lep ):
        lv_tl = reco['tl']
        lv_th = reco['proxy'] 
        acb, ach, act = top_angle_util.calc_angles( lv_tl, lv_lep )
        #print 'TMPDBG @ aa tl:',U.lv_as_PM_4tup( lv_tl ),' lv_lep:',U.lv_as_PM_4tup( lv_lep ),' (and prob:',reco['prob'],')'\
        #      '\n    -> cb ch ct:',acb, ach, act
        return [ reco['prob'], reco['cmtt'], lv_tl.Rapidity(), lv_th.Rapidity(), acb, ach, act ]

    # ========================================================================================================
    def reco( self, lv_lep, lv_jets, mex, mey, true_leptonic_b_index, cheat = False,
              btag_vals = [], known_problem_with_lepton = False ):
        '''Does a partial reconstruction of ttbar --> ( b (l nu) ) [ b q ... ].

           The "cheat" option means to consider only the correct assignment (if any).
           The parton_matches and known_problem_with_lepton are used to cheat and to create debug information
           '''
        Nj = len( lv_jets )
        assert Nj == 3
        assert (not cheat) or true_leptonic_b_index in range(3)
        dbg = self.dbg
        if dbg > 50:
            print 'partial_reco inputs are reasonable.'

        self.matchable = true_leptonic_b_index >= 0 and not known_problem_with_lepton

        # prepare the observables (one per reco-assignments)
        SR = self.SR
        OK = SR.reco_nu( mex, mey, lv_lep )
        if not OK:
            raise Exception( 'Unexpected failure of simple_reco::reco_nu' )
        if SR.two_solutions() :
            lv_nu_pzs = [ SR.nu(), SR.alt_nu() ]
        if dbg > 29:
            print 'DBG30PR SR(',mex,',',mey,',',U.lv_as_EPM_str(lv_lep),') -> '\
                  '2?',self.SR.two_solutions(),'nu:',U.lv_as_3tup(SR.nu()),'/',U.lv_as_3tup(SR.alt_nu())
        lv_nu_iLs = []
        lv_tl_iLs = []
        mtls = []
        mps = []
        lv_3j = sum( lv_jets, LV.new() )
        for iL in range( 3 ) :
            self.perm_is_correct = iL == true_leptonic_b_index
            lv_LBL = lv_jets[ iL ] + lv_lep
            if self.SR.two_solutions() :
                leptonic_tops = [ lv_LBL + lv for lv in lv_nu_pzs ] 
                dmtops = [ abs( self.top_mass - lv.M() ) for lv in leptonic_tops ] 
                iNu = 0 if dmtops[0] < dmtops[1] else 1
                if dbg > 18:
                    print 'DBG19PR iL:',iL,'-> perm_is_correct:',self.perm_is_correct,', dmtops:',dmtops,'-> iNu:',iNu
                lv_nu_iLs.append( lv_nu_pzs[ iNu ] )
                lv_tl = leptonic_tops[ iNu ]
            else:
                lv_nu = self.SR.nu()
                lv_nu_iLs.append( lv_nu )
                lv_tl = lv_LBL + lv_nu
            if dbg > 14:
                print 'DBG15PR iL:',iL,'-> nu:',U.lv_as_EPM_str(lv_nu_iLs[-1]),', tl:',U.lv_as_EPM_str(lv_tl)
            lv_tl_iLs.append( lv_tl )
            mtls.append( lv_tl.M() )
            mps.append( ( lv_3j - lv_jets[ iL ] ).M() )

        # loop over statistics-assignments
        recos = []
        for iL in range( 3 ) :
            self.perm_is_correct = iL == true_leptonic_b_index
            if cheat and not self.matchable and self.perm_is_correct:
                continue
 
            lv_nu = lv_nu_iLs[ iL ]
            lv_tl = lv_tl_iLs[ iL ]

            for iH in range( 3 ) : # loop over hadronic b assignment
                if iH == iL :
                    continue
                iW = 3 - iL - iH
                iLHW = [iL,iH,iW]
                
                # probabilities for hypothesis H: lhq   (i.e. a jet from the hadronic W is missing)
                prob_t_H = U.hist_value_at( self.h_mt_l, mtls[ iL ] ) * \
                           U.hist_value_at( self.h_mt_h, mtls[ iH ] ) * \
                           U.hist_value_at( self.h_mt_q, mtls[ iW ] )
                prob_p_H = U.hist_value_at( self.h_mp_hq, mps[ iL ] ) * \
                           U.hist_value_at( self.h_mp_lq, mps[ iH ] ) * \
                           U.hist_value_at( self.h_mp_hl, mps[ iW ] )
                prob_tag_H = self.calc_tag_prob_H( iW, btag_vals )
                prob_H = prob_t_H * prob_p_H * prob_tag_H

                # probabilities for hypothesis Q: lqq   (i.e. the hadronic b is missing)
                prob_t_Q = U.hist_value_at( self.h_mt_l, mtls[ iL ] ) * \
                           U.hist_value_at( self.h_mt_q, mtls[ iH ] ) * \
                           U.hist_value_at( self.h_mt_q, mtls[ iW ] )
                prob_p_Q = U.hist_value_at( self.h_mp_qq, mps[ iL ] ) * \
                           U.hist_value_at( self.h_mp_lq, mps[ iH ] ) * \
                           U.hist_value_at( self.h_mp_lq, mps[ iW ] )
                prob_tag_Q = self.calc_tag_prob_Q( iL, btag_vals )
                prob_Q = prob_t_Q * prob_p_Q * prob_tag_Q

                if dbg > 28:
                    print 'DBG29PR iLHW:',iLHW,'-> prob_tag_H:',prob_tag_H,'prob_tag_Q:',prob_tag_Q

                prob = self.fQ * prob_Q + ( 1-self.fQ ) * prob_H
                postriori_fQ = ( self.fQ * prob_Q / prob ) if prob > 0 else -1

                lv_th_proxy = lv_jets[ iH ] + lv_jets[ iW ]

                xx = R.TMath.Range( 0.2, 1.2, 100 / lv_th_proxy.E() )
                alpha = 0.6500 + 0.7076 * xx
                lv_ttbar = lv_tl + lv_th_proxy * alpha                

                recos.append( { 'iLHW' : iLHW, 'tagH' : prob_tag_H, 'tagQ' : prob_tag_Q,
                                'tl' : lv_tl, 'proxy' : lv_th_proxy, 'nu' : lv_nu,
                                'alpha' : alpha, 'ttbar' : lv_ttbar, 'cmtt' : ( lv_ttbar.M() - 64.44 ) / 0.7596,
                                'fQ' : postriori_fQ, 'prob' : prob } )
                if dbg > 34:
                    print 'DBG35PR alpha:',alpha,'prob_Q:',prob_Q,'prob_H',prob_H
         
        recos.sort( key=lambda reco: - reco[ 'prob' ] )

        denom = sum( [ reco['prob'] for reco in recos ] )
        if denom <= 0:
            if dbg :
                print 'Information (DBG1): partial reco probs are non-positive:',\
                      [ reco['prob'] for reco in recos ],'(legitimate as long as it is rare)'
            averages = self.averagables_array( recos[0], lv_lep ) [1:]
        else:
            cands = [ self.averagables_array( reco, lv_lep ) for reco in recos ]
            if dbg > 21 :
                print 'DBG22 partial_reco\'s cands for averaging are',cands
            weighted_cands = [ [cand[0]*obs for obs in cand[1:] ]  for cand in cands ]
            averages = [ sum(column)/denom for column in zip( *weighted_cands ) ]
            if dbg > 14 :
                print 'DBG15 partial_reco\'s averages:',averages

        return recos, { 'mtt': averages[ 0 ], 'yl': averages[ 1 ], 'yh': averages[ 2 ],
                        'cb' : averages[ 3 ], 'ch': averages[ 4 ], 'ct': averages[ 5 ] }

