
#####################################
def lv_as_3tup( lv ):
    return [lv.X(),lv.Y(),lv.Z()]

#####################################
def lv_as_PM_4tup( lv ):
    return [lv.X(),lv.Y(),lv.Z(),lv.M()]

#####################################
def lv_as_PtEtaPhiM_4tup( lv ):
    return [lv.Pt(),lv.Eta(),lv.Phi(),lv.M()]

#####################################
def lv_as_EPM_5tup( lv ):
    return [lv.E(),lv.X(),lv.Y(),lv.Z(),lv.M()]

#####################################
def lv_as_EPM_str( lv ):
    return '{:9.3f};{:8.3f},{:8.3f},{:8.3f};{:6.3f}'.format(lv.E(),lv.X(),lv.Y(),lv.Z(),lv.M())

#####################################
def hist_value_at( hist, x ):
    return hist.GetBinContent( hist.FindBin( x ) )
    


