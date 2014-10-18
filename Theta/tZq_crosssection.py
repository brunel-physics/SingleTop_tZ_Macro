from ROOT import gROOT, TCanvas, TF1, TH1F, TFile, TGraph
from array import array

import numpy as n

dosysttable = True

dobiasscan  = False
 
options = Options()
options.set('minimizer', 'strategy', 'newton_vanilla')

 
# for model building:
def get_model(signalname):

    # Read in and build the model automatically from the histograms in the root file. 
    # This model will contain all shape uncertainties given according to the templates-
    # which also includes rate changes according to the alternate shapes.
    # For more info about this model and naming conventuion, see documentation
    # of build_model_from_rootfile.

    model = build_model_from_rootfile('../TMVA/template_theta.root',  include_mc_uncertainties=True)
    
      
    
    
    # If the prediction histogram is zero, but data is non-zero, teh negative log-likelihood
    # is infinity which causes problems for some methods. Therefore, we set all histogram
    # bin entries to a small, but positive value:
    model.fill_histogram_zerobins()

    # define what the signal processes are. All other processes are assumed to make up the 
    # 'background-only' model.
    model.set_signal_processes('tZq')
    
    
    model.add_lognormal_uncertainty('ZZ_rate', math.log(1.3), 'ZZ')
    model.add_lognormal_uncertainty('WZ_rate', math.log(1.3), 'WZ')
    model.add_lognormal_uncertainty('Zjets_rate', math.log(10), 'Zjets')
    model.add_lognormal_uncertainty('TTZ_rate', math.log(1.3), 'TTZ')
    model.add_lognormal_uncertainty('TTW_rate', math.log(1.3), 'TTW')
    
    # Add some lognormal rate uncertainties. The first parameter is the name of the
    # uncertainty (which will also be the name of the nuisance parameter), the second
    # is the 'effect' as a fraction, the third one is the process name. The fourth parameter
    # is optional and denotes the channl. The default '*' means that the uncertainty applies
    # to all channels in the same way.
    # Note that you can use the same name for a systematic here as for a shape
    # systematic. In this case, the same parameter will be used; shape and rate changes 
    # will be 100% correlated.
    
    
   
    
    return model

# -------------- TO CHANGE BY THE USER
signalname = 'tZq'
# -------------- TO CHANGE BY THE USER
model = get_model(signalname)


# first, it is a good idea to generate a summary report to make sure everything has worked
# as expected. The summary will generate quite some information which should it make easy to spot
# errors like typos in the name of uncertainties, missing shape uncertaintie, etc.
model_summary(model)


signal_shapes = {'tZq': ['tZq']}  

fit = mle(model, input = 'data', n = 1, signal_process_groups = signal_shapes, with_covariance=True, with_error=True, ks = True, chi2 = True, options = options)




one_sigma = 0.6827
two_sigma = 0.95

print ("measurement of the cross-section")
res = pl_interval(model, 'data', n=1, cls = [one_sigma], signal_process_groups = signal_shapes, options = options )
#twi keys 'tZq' and the interval "one_sigma", it returns a list of double entries : lower and upper bound
print [ "%.3f" % res['tZq'][0][0] , "%.3f" %res['tZq'][one_sigma][0][0] , "%.3f" %res['tZq'][one_sigma][0][1] ]

tZq_init_xs = 245.8


tZq_fit  = tZq_init_xs*res['tZq'][0][0]
tZq_down = tZq_init_xs*res['tZq'][one_sigma][0][0]
tZq_up   = tZq_init_xs*res['tZq'][one_sigma][0][1] 

print ["fitted cross section ", "%.1f" %tZq_fit]
print ["down variation       ", "%.1f" %tZq_down]
print ["up variation         ", "%.1f" %tZq_up]




