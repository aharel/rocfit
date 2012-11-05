#!/usr/bin/python2.7
import ROOT as R

lvClass = None
def new(*args) :
    '''Create a Lorentz vector object with optimized class lookup '''
    global lvClass
    if lvClass is None : lvClass = R.TLorentzVector
    return lvClass(*args)
