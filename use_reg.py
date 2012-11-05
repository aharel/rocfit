#!/usr/bin/python2.7
#
# Use a TMVA-regression to combine the various reconstructed ttbar invariant masses
#
import ROOT as R
import array

# ========================================================================================================
class use_reg( object ):
    '''Use a TMVA regression to reconstruct mtt in the 4-jet channel, based on inputs
    from ttbar reconstruction algorithms.
    '''
    def __init__( self, dbg = 0, method_name = 'BDTGA', weights_file_name = 'weights/regA_BDTG.weights.xml' ) :
        self.dbg = dbg
        self.title = '{} method'.format( method_name )
        self.reader = R.TMVA.Reader('!Color:!Silent')
        input_names = [ 'r_mtt', 'mll', 'mtt', 'omtt', 'fM', 'chi2', 'aqz', 'jet1_m', 'm3f' ]
        input_expressions = [ 'r_mtt', 'rf_mll', 'sr_mtt', 'sr_omtt', 'fM', 'chi2', '(sr_qz+sr_oqz)/2', 'jet1_m', 'sr_m3f' ]
        self.inputs = {}
        for name, expression in zip( input_names, input_expressions ):
            self.inputs[ name ] = array.array( 'f', [ 0.0 ] )
            self.reader.AddVariable( expression, self.inputs[ name ] )
            if dbg > 19 : print 'DBG10 connected input "'+name+'" to the expression:',expression

        self.method = self.reader.BookMVA( self.title, weights_file_name )
        if dbg > 9 : print 'DBG10 use_reg constructor is done'

    # ----------------------------------------------------------------------------------------------------

    def reco( self, r_mtt, mll, mtt, omtt, fM, chi2, aqz, jet1_m, m3f ) : # the exact names matter...
        v = vars()
        inp = self.inputs
        for name in self.inputs:
            inp[ name ][0] = v[ name ]
        if self.dbg > 11 : print 'DBG12 inputs are:',self.inputs
        return self.reader.EvaluateRegression( self.title ).at( 0 )

    # ----------------------------------------------------------------------------------------------------

    def test( self ):
        inputs = [ \
            [ 651.05528, 40.840890, 909.46722, 912.34009, 920.44348, 24.781261, -17.41778, 41.095758, 0.9466060 ],
            [ 709.10907, 46.303872, 939.23779,  1017.849, 1041.1175, 0.9984394, -68.72002, 48.136325, 0.9271105 ],
            [ 905.40790, 40.340546, 1015.3878, 1017.6752, 1054.3913, 0.0331445, -13.61622, 117.77924, 0.9579159 ],
            [ 940.73169, 32.190947, 667.84217, 779.11416, 797.89812, 0.1227590, -29.44752, 26.629362, 0.6211921 ],
            [ 649.29452, 46.807339, 1071.6661, 1128.3834, 956.47634, 6.9037814, -79.25628, 112.28417, 0.9628700 ],
            [ 931.68012, 30.616773, 908.08721, 946.80167, 859.28443, 2.6892142, 74.806729, 26.656368, 0.9219670 ],
            [ 757.50524, 39.186823, 738.85556, 738.85556, 747.27679, 4.2093849, 71.062260,  74.30627, 0.9413877 ],
            [ 765.785  , 45.881953, 1077.6600, 1082.1603, 1070.9941, 1.3074210, 15.813060, 28.872190, 0.9197110 ],
            [ 585.52020, 54.189904, 1021.1321, 1096.7568, 999.27276, 15.895190, 61.961495, 51.574401, 0.9663761 ],
            [ 980.24503, 27.157816, 961.67987, 961.67987, 1000.8387, 0.3237953, -200.1655, 25.137892, 0.9195946 ],
            [ 419.51225, 25.003363, 428.43238, 517.77784, 551.64863, 0.1771327, 12.332707, 13.973769, 0.9414695 ],
            [ 872.58780, 30.631664, 630.49256, 937.38110, 911.89192, 1.2742829, 125.81705, 23.955900, 0.8986271 ],
            [ 655.64482,  30.19839, 739.00400, 806.61362, 715.57065, 1.9442286, -83.77151, 30.694011, 0.9138584 ],
            [ 960.08169, 26.617724, 965.41200, 973.56509, 992.80458, 1.1723325, -8.057478, 32.914131, 0.9583430 ],
            [ 788.17840, 34.242875, 724.37848,  899.5137, 766.02673, 0.1604083, 75.648040, 17.548141, 0.8919289 ] ]

        answers = [ 839.15228, 870.35113, 926.41436, 818.88604, 1006.6489, 898.00164, 740.78637, 958.86071, 808.24176, 944.02343,
                    460.61193, 821.35992, 682.85498, 941.99829, 762.40429 ]

        all_pass = True
        for input, answer in zip( inputs, answers ):
            mtt = self.reco( *input )
            print '{:7f} =? {:7f} from {:6f}, {:6f}, {:6f},   {:6f}, {:6f}, {:6f},   {:6f}, {:6f}, {:6f}'\
                  ''.format( * ([ mtt, answer ] + input ) )
            if abs( mtt - answer ) > 0.1 :
                all_pass = False
                print '**** Test failed!'

        if all_pass:
            print '\nAll',len( answers ),'tests passed'
        else:
            print '\nOOPS! failed at least one test'
        
