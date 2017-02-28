''' FWLiteReader example: Loop over a sample and write some data to a histogram.
'''
# Standard imports
import os
import logging
import ROOT
import array

#RootTools
from RootTools.core.standard import *

#Helper
import JetMET.tools.helpers as helpers

# argParser
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel', 
      action='store',
      nargs='?',
      choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'],
      default='INFO',
      help="Log level for logging"
)

args = argParser.parse_args()
logger = get_logger(args.logLevel, logFile = None)

max_events = -1
max_files = -1


ref    = FWLiteSample.fromDirectory("ref", "/afs/cern.ch/user/s/schoef/eos/cms/store/group/phys_jetmet/mdjordje/NEW_JEC/QCD_default/QCD_Pt-15to3000_TuneCUETP8M1_Flat_13TeV_pythia8/crab_test_QCD_default/170220_213824/0000", maxN = max_files)
sigma1 = FWLiteSample.fromDirectory("sigma1", "/afs/cern.ch/user/s/schoef/eos/cms/store/group/phys_jetmet/mdjordje/NEW_JEC/QCD_HF/QCD_Pt-15to3000_TuneCUETP8M1_Flat_13TeV_pythia8/crab_test_QCD_HF/170220_213904/0000", maxN = max_files)

preprefix = "milos_sigma1_HF_HLT_jets" 

comparisons = \
[
    [ "HFsigma1", sigma1, ref ],
]

# define TProfiles
pt_thresholds = [ 10**(x/10.) for x in range(11,36) ] 

eta_thresholds = [0, 1.3, 2.5, 3.2]

#color={0:ROOT.kBlack, 1.3:ROOT.kBlue, 2.5:ROOT.kRed, 3.2:ROOT.kMagenta, "12_21_06_38":ROOT.kBlack, "12_21_06_56":ROOT.kBlue, "12_21_07_12":ROOT.kGreen, "12_21_07_52":ROOT.kMagenta}

resp = {}
resp_eta = {}
for comp in [c[0] for c in comparisons]:
    resp_eta[comp] = ROOT.TProfile("resp_eta_%s"%comp, "resp_eta_%s"%comp, 26, 0, 5.2 )
    #resp_eta[comp].style = styles.lineStyle( color[comp] )
    resp_eta[comp].legendText = comp 

    resp[comp]={}
    for i_eta_th, eta_th in enumerate( eta_thresholds ):
        resp[comp][eta_th] = ROOT.TProfile("response_%s"%comp, "response_%s"%comp, len(pt_thresholds)-1, array.array('d', pt_thresholds) )
        #resp[comp][eta_th].style = styles.lineStyle(color[eta_th] )
        resp[comp][eta_th].legendText = "%2.1f<=#eta"%eta_th
        if eta_th!=eta_thresholds[-1]: resp[comp][eta_th].legendText += "<%2.1f"%eta_thresholds[i_eta_th+1]

products = {
    'jets':      {'type': 'vector<reco::PFJet>', 'label':"hltAK4PFJets"},
    }


for comp, legacy, ref in comparisons:
    r1 = legacy.fwliteReader( products = products )
    r2 = ref.fwliteReader( products = products )

    r1.start()
    runs_1 = set()
    position_r1 = {}
    count=0
    while r1.run( readProducts = False ):
            position_r1[r1.evt] = r1.position-1
            count+=1
            if max_events is not None and max_events>0 and count>=max_events:break

    r2.start()
    runs_2 = set()
    position_r2 = {}
    count=0
    while r2.run( readProducts = False ):
            position_r2[r2.evt] = r2.position-1
            count+=1
            if max_events is not None and max_events>0 and count>=max_events:break

    logger.info( "Have %i events in first samle and %i in second", len(position_r1), len(position_r2) )

    # Fast intersect
    intersec = set(position_r1.keys()).intersection(set(position_r2.keys()))
    positions = [(position_r1[i], position_r2[i]) for i in intersec]

    # Without sorting, there is a jump between files with almost every event -> extremly slow
    positions.sort()
    logger.info("Have %i events in common.", len(intersec))

    #Looping over common events
    for i, p in enumerate(positions):
        p1,p2 = p
        r1.goToPosition(p1)
        r2.goToPosition(p2)
        if i%10000==0: logger.info("At %i/%i of common events.", i, len(positions))
        jets1_ = [ j for j in r1.products['jets'] if helpers.jetID( j )]
        jets2_ = [ j for j in r2.products['jets'] if helpers.jetID( j )]
        jets1 = [{'pt':j.pt(), 'eta':j.eta(), 'phi':j.phi(), 'j':j} for j in jets1_]
        jets2 = [{'pt':j.pt(), 'eta':j.eta(), 'phi':j.phi(), 'j':j} for j in jets2_]
        for c in zip(jets1, jets2):
            if helpers.deltaR2(*c)<0.2**2:
                if not ( helpers.jetID(c[0]['j']) and helpers.jetID(c[1]['j']) ): continue
                resp_eta[comp].Fill( abs(c[0]['eta']), c[0]['pt']/c[1]['pt'] )
                for eta_th in reversed(eta_thresholds):
                    if abs(c[0]['eta'])>eta_th:
                        resp[comp][eta_th].Fill( c[0]['pt'], c[0]['pt']/c[1]['pt'] )
                        break

    # Make plot
    profiles = [resp[comp][t] for t in eta_thresholds]
    #profiles = [ jetResponse_NJC]
    prefix=preprefix + "_" + legacy.name + ("_max_events_%s_"%max_events if max_events is not None and max_events>0 else "" )
    histos = [ [h.ProjectionX()] for h in profiles ]
    for i, h in enumerate(histos):
        h[0].__dict__.update(profiles[i].__dict__)
    
    jetResponsePlot = Plot.fromHisto(name = prefix+"jetResponseRatio_relval", histos = histos, texX = "new Jet p_{T}" , texY = "response ratio new/old" )
    plotting.draw(jetResponsePlot, plot_directory = "/afs/hephy.at/user/r/rschoefbeck/www/etc/", ratio = None, logY = False, logX = True, yRange=(0.7,1.2))

# Make eta plot
profiles = [resp_eta[comp] for comp in [c[0] for c in comparisons] ]
prefix=preprefix + "_eta_" + ("_max_events_%s_"%max_events if max_events is not None and max_events>0 else "" )
histos = [ [h.ProjectionX()] for h in profiles ]
for i, h in enumerate(histos):
    h[0].__dict__.update(profiles[i].__dict__)

jetResponsePlot = Plot.fromHisto(name = prefix+"jetResponseRatio_eta_relval", histos = histos, texX = "new Jet #eta" , texY = "response ratio new/old" )
plotting.draw(jetResponsePlot, plot_directory = "/afs/hephy.at/user/r/rschoefbeck/www/etc/", ratio = None, logY = False, logX = False, yRange=(0.7,1.2))
