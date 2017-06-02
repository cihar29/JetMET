#!/usr/bin/env python
''' Analysis script for L3 residuals (Z balancing) 
'''
#
# Standard imports and bateh mode
#
import ROOT
ROOT.gROOT.SetBatch(True)
import itertools
import os
import array
import pickle

from math                                import sqrt, cos, sin, pi, atan2
from RootTools.core.standard             import *
from JetMET.tools.user                   import plot_directory as user_plot_directory
from JetMET.tools.helpers                import deltaPhi, deltaR

# Object selection
from JetMET.tools.objectSelection        import getFilterCut, getJets, jetVars

#
# Arguments
# 
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',           action='store',      default='INFO',          nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging" )
argParser.add_argument('--triggers',           action='store',      default='DiPFJetAve',    nargs='?', choices=['DiPFJetAve', 'DiPFJetAve_HFJEC', 'PFJet'], help="trigger suite" )
argParser.add_argument('--ptBinningVar',       action='store',      default='ave',           nargs='?', choices=['ave', 'tag'], help="jet pT binning variable (pT avg or pT tag)" )
argParser.add_argument('--era',                action='store',      default='Run2016',       nargs='?', choices=['Run2016', 'Run2016BCD', 'Run2016EFearly', 'Run2016FlateG', 'Run2016H'], help="era" )
argParser.add_argument('--small',                                   action='store_true',     help='Run only on a small subset of the data?')#, default = True)
argParser.add_argument('--overwrite',                               action='store_true',     help='Overwrite results.pkl?')
argParser.add_argument('--plot_directory',     action='store',      default='JEC/L2res_v3',  help="subdirectory for plots")
args = argParser.parse_args()

if args.ptBinningVar == 'tag':
    args.plot_directory += '_tagJetPtBin'
    pt_binning_variable = 'Jet_pt[tag_jet_index]'
    pt_binning_legendText = 'p_{T,tag} '
else:
    pt_binning_legendText = 'p_{T,avg} '
    pt_binning_variable = "pt_avg"
if args.small:
    args.plot_directory += '_small'

plot_directory = os.path.join( user_plot_directory, args.plot_directory, args.era, args.triggers )

# Lumi for MC
lumi = 35.9
# DrawLatex objects for plots
tex = ROOT.TLatex()
tex.SetNDC()
tex.SetTextSize(0.04)
tex.SetTextAlign(11) # align right
def drawObjects( dataMCScale, lumi ):
    lines = [
      #(0.15, 0.95, args.era), 
      (0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV) Scale %3.2f'% ( lumi, dataMCScale ) )
    ]
    return [tex.DrawLatex(*l) for l in lines] 

# Formatting for 1D plots
def draw1DPlots(plots, dataMCScale):
  for log in [ True]:
    plot_directory_ = os.path.join(plot_directory, ("log" if log else "") )
    for plot in plots:
      #if not max(l[0].GetMaximum() for l in plot.histos): continue # Empty plot
      p_drawObjects = map( lambda l:tex.DrawLatex(*l), getattr(plot, "drawObjects", [] ) )

      if hasattr( plot, "subdir"):
        plot_directory__ = os.path.join( plot_directory_, plot.subdir)
      else:
        plot_directory__ = plot_directory_

      plotting.draw(plot,
        plot_directory = plot_directory__,
        #ratio          = {'yRange':(0.6,1.4)} if len(plot.stack)>=2 else None,
        logX = False, logY = log, sorting = False,
        yRange         = (0.0003, "auto") if log else (0.001, "auto"),
        #scaling        = {0:1} if len(plot.stack)==2 else {},
        legend         = [ (0.15,0.91-0.035*5,0.95,0.91), 2 ],
        drawObjects    = drawObjects( dataMCScale , lumi ) + p_drawObjects
      )

#
# Logger
#
import JetMET.tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)


from JetMET.JEC.samples.L2res_skim import *

if args.era == 'Run2016':
    data = JetHT_Run2016
elif args.era == 'Run2016BCD':
    data = JetHT_Run2016BCD
elif args.era == 'Run2016EFearly':
    data = JetHT_Run2016EFearly
elif args.era == 'Run2016FlateG':
    data = JetHT_Run2016FlateG
elif args.era == 'Run2016H':
    data = JetHT_Run2016H


if args.triggers=='DiPFJetAve':
    triggers = [ 
        "HLT_DiPFJetAve40",
        "HLT_DiPFJetAve60",
        "HLT_DiPFJetAve80",
        "HLT_DiPFJetAve140",
        "HLT_DiPFJetAve200",
        "HLT_DiPFJetAve260",
        "HLT_DiPFJetAve320",
        "HLT_DiPFJetAve400",
        "HLT_DiPFJetAve500",
    ]
elif args.triggers == 'PFJet':
    triggers = [
        "HLT_PFJet40",
        "HLT_PFJet60",
        "HLT_PFJet80",
        "HLT_PFJet140",
        "HLT_PFJet200",
        "HLT_PFJet260",
        "HLT_PFJet320",
        "HLT_PFJet400",
        "HLT_PFJet450",
        "HLT_PFJet500",
    ]
elif args.triggers == 'DiPFJetAve_HFJEC':
    triggers = [
        "HLT_DiPFJetAve60_HFJEC",
        "HLT_DiPFJetAve80_HFJEC",
        "HLT_DiPFJetAve100_HFJEC",
        "HLT_DiPFJetAve160_HFJEC",
        "HLT_DiPFJetAve220_HFJEC",
        "HLT_DiPFJetAve300_HFJEC",
    ]
else:
    triggers = [ args.triggers ]

mc = QCD_Pt
samples = [mc, data]

selection = [
   ("tgb", "abs(Jet_eta[tag_jet_index])<1.3"),
   ("btb", "cos(Jet_phi[tag_jet_index] - Jet_phi[probe_jet_index]) < cos(2.7)"),
   ("a30", "alpha<0.3"), 
]

for s in samples:   
    s.addSelectionString( "&&".join(c[1] for c in selection))
    if args.small:
        s.reduceFiles( to = 1 )

# Add trigger selection to data
data.addSelectionString( "("+"||".join(triggers)+")")

#colors = [ ROOT.kRed, ROOT.kBlue, ROOT.kGreen, ROOT.kMagenta, ROOT.kOrange, ROOT.kViolet,  ROOT.kCyan, ROOT.kOrange - 1, ROOT.kViolet - 1, ROOT.kCyan + 3]
colors = [ j+1 for j in range(0,9) ] + [ j+31 for j in range(9,18) ]

from JetMET.JEC.L2res.thresholds import pt_avg_thresholds, pt_avg_bins, eta_thresholds, abs_eta_thresholds

thresholds = [-1.2+x*2.4/96. for x in range(97)]

weightString   = "weight"

results_file = os.path.join( plot_directory, 'results.pkl' )

h = {}
p = {}

if os.path.exists( results_file ) and not args.overwrite:
    h, p = pickle.load( file (results_file ) )
    logger.info( "Loaded %s", results_file )
else: 
    for var in [ "A", "B" ]:
        h[var] = {}
        p[var] = {}
        for s in samples:
            logger.info( "Make TH3D for sample %s and variable %s", s.name, var )
            h[var][s.name] = ROOT.TH3D( "h_%s_%s"%( var, s.name), "h_%s_%s"%(var, s.name),\
                    len(thresholds) - 1, array.array('d', thresholds), 
                    len(eta_thresholds)-1, array.array('d', eta_thresholds), 
                    len(pt_avg_thresholds) - 1, array.array('d', pt_avg_thresholds)  
                )

            weight_ = "("+s.selectionString+")*("+s.combineWithSampleWeight(weightString)+")"
            varString_ = pt_binning_variable+":Jet_eta[probe_jet_index]:%s>>h_%s_%s"%( var, var, s.name )

            logger.info("Using %s %s", varString_, weight_ ) 
            s.chain.Draw( varString_, weight_, 'goff')
            
            #logger.info( "Make TProfile2D for sample %s and variable %s", s.name, var )
            #p[var][s.name] = ROOT.TProfile2D( "p_%s_%s"%( var, s.name ), "p_%s_%s"%( var, s.name ),\
            #        len(eta_thresholds)-1, array.array('d', eta_thresholds), 
            #        len(pt_avg_thresholds) - 1, array.array('d', pt_avg_thresholds)  
            #    )
            #s.chain.Draw(pt_binning_variable+":Jet_eta[probe_jet_index]:%s>>p_%s,%s"%( var, var, s.name ),  "("+s.selectionString+")*("+s.combineWithSampleWeight(weightString)+")", 'goff')

    if not os.path.exists(os.path.dirname( results_file )): os.makedirs( os.path.dirname( results_file ) ) 
    pickle.dump( ( h, p ), file( results_file, 'w' ) )
    logger.info( "Written %s", results_file )

# Make all the projections
# x ... A,B
# y ... eta
# z ... pt_avg
projections = {}
for var in [ "A", "B" ]:
    projections[var] = {}
    for s in samples:
        projections[var][s.name] = {'neg_eta':{}, 'pos_eta':{}, 'abs_eta':{}}
        for i_aeta in range(len(abs_eta_thresholds)-1):
            eta_bin     = tuple(abs_eta_thresholds[i_aeta:i_aeta+2])
            neg_eta_bin = (-eta_bin[1], -eta_bin[0])
            bin_y       = h[var][s.name].GetYaxis().FindBin( 0.5*sum(eta_bin)  ) 
            neg_bin_y   = h[var][s.name].GetYaxis().FindBin( 0.5*sum(neg_eta_bin)  ) 

            projections[var][s.name]['neg_eta'][eta_bin] = {}
            projections[var][s.name]['pos_eta'][eta_bin] = {}
            projections[var][s.name]['abs_eta'][eta_bin] = {}
            for i_pt_avg in range(len(pt_avg_thresholds)-1):
                pt_avg_bin = tuple(pt_avg_thresholds[i_pt_avg:i_pt_avg+2])
                bin_z      = h[var][s.name].GetZaxis().FindBin( 0.5*sum(pt_avg_bin)  ) 
                projections[var][s.name]['abs_eta'][eta_bin][pt_avg_bin] = h[var][s.name].ProjectionX("abs_%s_%s_%i_%i" % ( s.name, var, bin_y, bin_z ) , bin_y, bin_y, bin_z, bin_z)
                projections[var][s.name]['pos_eta'][eta_bin][pt_avg_bin] = h[var][s.name].ProjectionX("pos_%s_%s_%i_%i" % ( s.name, var, bin_y, bin_z ) , bin_y, bin_y, bin_z, bin_z)
                projections[var][s.name]['neg_eta'][eta_bin][pt_avg_bin] = h[var][s.name].ProjectionX("neg_%s_%s_%i_%i" % ( s.name, var, neg_bin_y, bin_z ) , neg_bin_y, neg_bin_y, bin_z, bin_z)
                projections[var][s.name]['abs_eta'][eta_bin][pt_avg_bin].Add(  projections[var][s.name]['neg_eta'][eta_bin][pt_avg_bin] ) 

                # Normalize to 1
                for eta_flav in ['neg_eta', 'pos_eta', 'abs_eta']:
                    integral = projections[var][s.name][eta_flav][eta_bin][pt_avg_bin].Integral()
                    if integral>0: projections[var][s.name][eta_flav][eta_bin][pt_avg_bin].Scale(1./integral)

            for eta_flav in ['neg_eta', 'pos_eta', 'abs_eta']:
                histos = []
                for i_pt_avg_bin, pt_avg_bin in enumerate(pt_avg_bins):
                    histos.append( projections[var][s.name][eta_flav][eta_bin][pt_avg_bin])
                    histos[-1].style = styles.lineStyle( colors[ i_pt_avg_bin ] ) 
                    histos[-1].legendText = "%i #leq %s < %i" % ( pt_avg_bin[0], pt_binning_legendText, pt_avg_bin[1] )

                name = "%s_%s_%s_%i_%i" % ( s.name.replace('_'+args.era, ''), var, eta_flav, 1000*eta_bin[0], 1000*eta_bin[1] )
                plot = Plot.fromHisto( name, [ [histo] for histo in histos], texX = var, texY = "Number of Events" )    
                draw1DPlots( [plot], 1.)
