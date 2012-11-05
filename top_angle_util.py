#!/usr/bin/python2.7
import ROOT as R
import lorentzv as LV
from math import cos

def calc_angles( lv_top, lv_lep, Qlep = 1 ):
    '''Calculated the cosine of the angle of the lepton in the (leptonic) top rest frame
    relative to three axes: the beam axis, the helicity axis, and the transverse axis'''
    lv_lep_tag = LV.new( lv_lep )
    lv_lep_tag.Boost( -lv_top.BoostVector() )
    lv_beam = R.TVector3( 980 * Qlep, 0, 0 ) # haven't checked the sign. Maybe other way around.
    lv_perp = lv_beam.Cross( lv_top.Vect() )
    beam_angle = lv_beam.Angle(  lv_lep_tag.Vect() )
    helicity_angle = lv_top.Vect().Angle( lv_lep_tag.Vect() )
    transverse_angle = lv_perp.Angle( lv_lep_tag.Vect() )
    return cos( beam_angle ), cos( helicity_angle ), cos( transverse_angle )

