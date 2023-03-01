#!/usr/bin/python

import os, sys, time, math, datetime, pickle, itertools, getopt
import numpy as np
from array import array
from argparse import ArgumentParser

sys.path.append( os.path.dirname( "../" ) ) 

from utils import contains_category, hist_tag
import config
from xsec import xsec

parser = ArgumentParser()
parser.add_argument( "-v", "--variable", default = "HT" )
parser.add_argument( "-y", "--year", default = "17" )
parser.add_argument( "-l", "--lepton", default = "E" )
parser.add_argument( "-nh", "--nhot", default = "0p" )
parser.add_argument( "-nt", "--nt", default = "0p" )
parser.add_argument( "-nw", "--nw", default = "0p" )
parser.add_argument( "-nb", "--nb", default = "2p" )
parser.add_argument( "-nj", "--nj", default = "5p" )
parser.add_argument( "-sd", "--subDir", default = "test" )
args = parser.parse_args()

if args.year == "16APV":
  import samplesUL16APV as samples 
elif args.year == "16":
  import samplesUL16 as samples
elif args.year == "17":
  import samplesUL17 as samples
elif args.year == "18":
  import samplesUL18 as samples
else:
  quit( "[ERR] Invalid -y (--year) option used. Quitting..." )

import ROOT

ROOT.gROOT.SetBatch(1)
start_time = time.time()

category = {
  "LEPTON": [ args.lepton ],
  "NHOT": [ args.nhot ],
  "NT": [ args.nt ],
  "NW": [ args.nw ],
  "NB": [ args.nb ],
  "NJ": [ args.nj ]
}

groups = {
  "DAT": sorted( [ str( process ) for process in samples.samples[ "DAT" ] ] ),
  "SIG": sorted( [ str( process ) for process in samples.samples[ "SIG" ] ] ),
  "BKG": sorted( [ str( process ) for process in samples.samples[ "BKG" ] ] ),
  "UE": sorted( [ str( process ) for process in samples.samples[ "UE" ] ] ),
  "HD": sorted( [ str( process ) for process in samples.samples[ "HD" ] ] ),
  "TEST": [ str( process ) for process in samples.samples[ "TEST" ] ]
}

def read_tree( samplePath ):
  if not os.path.exists( samplePath ):
    print("[ERR] {} does not exist.  Exiting program...".format( samplePath ) )
    sys.exit(1)
  rootFile = ROOT.TFile.Open( samplePath, "READ" )
  rootTree = rootFile.Get( "ljmet" )
  return rootFile, rootTree

def analyze( rTree, nHist, year, process, variable, doSYST, doPDF, doABCDNN, category, verbose ):
  variableName = config.plot_params[ "VARIABLES" ][ variable ][0]
  histBins = array( "d", config.plot_params[ "VARIABLES" ][ variable ][1] )
  xLabel = config.plot_params[ "VARIABLES" ][ variable ][2]
  print( ">> Processing {} for 20{} {}".format( variable, year, process ) )

  
  # modify weights
  # scale up MC samples used in DNN training where dataset partitioned into 60/20/20 so scale isTraining==1 || isTraining==3 by 1.25
  mc_weights = { "NOMINAL": "1.25" if ( ( process.startswith( "TTTo" ) or process.startswith( "TTTW" ) or process.startswith( "TTTJ" ) ) and "DNN" in variable ) else "1" } # weights only applied to MC
  if process in xsec:
    nTrueHist = nHist[ process ]
    for splitPrefix in samples.split:
      if splitPrefix in process:
        for splitProcess in samples.split[ splitPrefix ]:
          if process != splitProcess: nTrueHist += nHist[ splitProcess ]
        print( "[INFO] {} was hadded into more than one file, consolidating numTrueHist across split files: {} --> {}".format( process, nHist[ process ], nTrueHist ) )
    mc_weights[ "PROCESS" ] = "( {:.2f} * {:.10f} / {:.1f} )".format( config.lumi[ args.year ], xsec[ process ], nTrueHist ) 
  elif process in samples.samples[ "DAT" ]:
    mc_weights[ "PROCESS" ] = "1"
  else:
    sys.exit( "[ERROR] {} is neither data nor does it have a cross-section listed. Exiting...".format( process ) )

  if doABCDNN:
    abcdnnTag      = config.params[ "ABCDNN" ][ "TAG" ]
    extABCDTF      = config.params[ "ABCDNN" ][ "TF" ][ args.year ]
    extABCDTFScale = config.params[ "ABCDNN" ][ "TF SCALE" ]
    abcdnnName     = variableName + "_{}".format( abcdnnTag )
    print( "   + Including ABCDnn Histograms with tag {}".format( abcdnnTag ) )
    mc_weights[ "ABCDNN" ] = "{} * {}".format( extABCDTF, extABCDTFScale, abcdnnTag ) 

  if process.startswith( "TTTo" ): # https://twiki.cern.ch/twiki/bin/view/Sandbox/JamesKeaveneySandbox
    mc_weights[ "NOMINAL" ] += " * topPtWeight13TeV"

  if process not in groups[ "DAT" ]:
    mc_weights[ "NOMINAL" ] += "*{}*{}".format( config.mc_weight, mc_weights[ "PROCESS" ] )

  if category["NB"] != "0p" and process not in groups[ "DAT" ]:
    mc_weights[ "NOMINAL" ] += " * btagDeepJetWeight * btagDeepJet2DWeight_HTnj" 

  if year in [ "16APV", "16", "17" ] and process not in groups ["DAT"]:
    mc_weights[ "NOMINAL" ] += " * L1NonPrefiringProb_CommonCalc"

  if year in [ "16APV" ]: # since single lepton triggers weren't produced for EOY in 2016APV, default to triggerX
    mc_weights[ "NOMINAL" ] = mc_weights[ "NOMINAL" ].replace( "triggerSF", "1" )
   
  if process not in groups[ "DAT" ] and doSYST:
    if config.systematics[ "MC" ][ "pileup" ][0]:
      mc_weights[ "PILEUP" ] = { "UP": mc_weights[ "NOMINAL" ].replace( "pileupWeight", "pileupWeightUp" ),
                                 "DN": mc_weights[ "NOMINAL" ].replace( "pileupWeight", "pileupWeightDown" ) }
    if config.systematics[ "MC" ][ "pileupJetID" ][0]:
      mc_weights[ "PILEUPJETID" ] = { "UP": mc_weights[ "NOMINAL" ].replace( "pileupJetIDWeight", "pileupJetIDWeightUp" ),
                                      "DN": mc_weights[ "NOMINAL" ].replace( "pileupJetIDWeight", "pileupJetIDWeightDown" ) }
    if year in [ "16APV", "16", "17" ] and config.systematics[ "MC" ][ "prefire" ]:
      mc_weights[ "PREFIRE" ] = { "UP": mc_weights[ "NOMINAL" ].replace("L1NonPrefiringProb_CommonCalc","L1NonPrefiringProbUp_CommonCalc"),
                                  "DN": mc_weights[ "NOMINAL" ].replace("L1NonPrefiringProb_CommonCalc","L1NonPrefiringProbDown_CommonCalc") }
    if config.systematics[ "MC" ][ "muRFcorrd" ][0]:
      mc_weights[ "MURFCORRD" ] = { "UP": "renormWeights[5] * {}".format( mc_weights[ "NOMINAL" ] ),
                                    "DN": "renormWeights[3] * {}".format( mc_weights[ "NOMINAL" ] ) }
    if config.systematics[ "MC" ][ "muR" ][0]:
      mc_weights[ "MUR" ] = { "UP": "renormWeights[4] * {}".format( mc_weights[ "NOMINAL" ] ),
                              "DN": "renormWeights[2] * {}".format( mc_weights[ "NOMINAL" ] ) }
    if config.systematics[ "MC" ][ "muF" ][0]:
      mc_weights[ "MUF" ] = { "UP": "renormWeights[1] * {}".format( mc_weights[ "NOMINAL" ] ),
                              "DN": "renormWeights[0] * {}".format( mc_weights[ "NOMINAL" ] ) }
    if config.systematics[ "MC" ][ "isr" ][0]: 
      mc_weights[ "ISR" ] = { "UP": "renormPSWeights[0] * {}".format( mc_weights[ "NOMINAL" ] ),
                              "DN": "renormPSWeights[2] * {}".format( mc_weights[ "NOMINAL" ] ) }
    if config.systematics[ "MC" ][ "fsr" ][0]:
      mc_weights[ "FSR" ] = { "UP": "renormPSWeights[1] * {}".format( mc_weights[ "NOMINAL" ] ),
                              "DN": "renormPSWeights[3] * {}".format( mc_weights[ "NOMINAL" ] ) }
    if config.systematics[ "MC" ][ "toppt" ][0]:
      mc_weights[ "TOPPT" ] = { "UP": "({}) * {}".format( "topPtWeight13TeV" if "TTTo" in process else "1", mc_weights[ "NOMINAL" ] ),
                                "DN": "(1/{}) * {}".format( "topPtWeight13TeV" if "TTTo" in process else "1", mc_weights[ "NOMINAL" ] ) }
    if config.systematics[ "MC" ][ "ABCDNNMODEL" ][0] and doABCDNN and "ABCDNNMODEL" in config.params[ "ABCDNN" ][ "SYSTEMATICS" ]:
      mc_weights[ "ABCDNN ABCDNNMODEL" ] = { "UP": mc_weights[ "ABCDNN" ].replace( abcdnnName, abcdnnName + "_MODELUP" ), 
                                             "DN": mc_weights[ "ABCDNN" ].replace( abcdnnName, abcdnnName + "_MODELDN" ) }
    if config.systematics[ "MC" ][ "ABCDNNCLOSURE" ][0] and doABCDNN and "ABCDNNCLOSURE" in config.params[ "ABCDNN" ][ "SYSTEMATICS" ]:
      mc_weights[ "ABCDNN ABCDNNCLOSURE" ] = { "UP": mc_weights[ "ABCDNN" ].replace( abcdnnName, abcdnnName + "_CLOSUREUP" ),
                                               "DN": mc_weights[ "ABCDNN" ].replace( abcdnnName, abcdnnName + "_CLOSUREDN" ) }

    # deep jet related systematics
    for syst in [ "LF", "lfstats1", "lfstats2", "HF", "hfstats1", "hfstats2", "cferr1", "cferr2" ]:
      if config.systematics[ "MC" ][ syst ]:
        mc_weights[ syst.upper() ] = {}
        for shift in [ "up", "dn" ]:
          mc_weights[ syst.upper() ][ shift.upper() ] = mc_weights[ "NOMINAL" ].replace( "btagDeepJetWeight", "btagDeepJetWeight_" + syst + shift ).replace( "btagDeepJet2DWeight_HTnj", "btagDeepJet2DWeight_HTnj_" + syst + shift )
  
  # modify cuts
  cuts = { "BASE": config.base_cut }
  if "TTToSemiLepton" in process and "HT500" in process: cuts[ "NOMINAL" ] = cuts[ "BASE" ] + " && isHTgt500Njetge9==1"
  elif "TTToSemiLepton" in process and "HT500" not in process: cuts[ "NOMINAL" ] = cuts[ "BASE" ] + " && isHTgt500Njetge9==0"
  else: cuts[ "NOMINAL" ] = cuts[ "BASE" ][:]
  if ( ( process.startswith( "TTTo" ) or process.startswith( "TTTW" ) or process.startswith( "TTTJ" ) ) and "DNN" in variable ):
    cuts[ "NOMINAL" ] += " && ( isTraining == 1 || isTraining == 3 )" # Used isTraining==1 in training and isTraining==2 in validation of training
  if year == "18" and category[ "LEPTON" ][0] == "E": # exclude electrons that fall in this HEM region which resulted in many misidentifications of jets as electrons
    cuts[ "NOMINAL" ] += " && ( leptonEta_MultiLepCalc > -1.3 || ( leptonPhi_MultiLepCalc < -1.57 || leptonPhi_MultiLepCalc > -0.87 ) )"

  cuts[ "LEPTON" ] = " && isElectron==1" if category[ "LEPTON" ][0] == "E" else " && isMuon==1"
  cuts[ "NHOT" ] = " && NresolvedTops1pFake {}= {}".format( ">" if "p" in category[ "NHOT" ][0] else "=", category[ "NHOT" ][0][:-1] if "p" in category[ "NHOT" ][0] else category[ "NHOT" ][0] )
  cuts[ "NT" ] = " && NJetsTtagged {}= {}".format( ">" if "p" in category[ "NT" ][0] else "=", category[ "NT" ][0][:-1] if "p" in category[ "NT" ][0] else category[ "NT" ][0] )
  cuts[ "NW" ] = " && NJetsWtagged {}= {}".format( ">" if "p" in category[ "NW" ][0] else "=", category[ "NW" ][0][:-1] if "p" in category[ "NW" ][0] else category[ "NW" ][0] )
  cuts[ "NB" ] = " && NJetsCSV_JetSubCalc {}= {}".format( ">" if "p" in category[ "NB" ][0] else "=", category[ "NB" ][0][:-1] if "p" in category[ "NB" ][0] else category[ "NB" ][0] )
  cuts[ "NJ" ] = " && NJets_JetSubCalc {}= {}".format( ">" if "p" in category[ "NJ" ][0] else "=", category[ "NJ" ][0][:-1] if "p" in category[ "NJ" ][0] else category[ "NJ" ][0] )
 
  cuts[ "NOMINAL" ] += cuts[ "LEPTON" ] + cuts[ "NHOT" ] + cuts[ "NT" ] + cuts[ "NW" ] + cuts[ "NB" ] + cuts[ "NJ" ]
  cuts[ "ABCDNN"] = cuts["BASE"] + cuts[ "NHOT" ] + cuts[ "NT" ] + cuts[ "NW" ] + cuts[ "NB" ] + cuts[ "NJ" ] + cuts[ "LEPTON" ] 
  # modify the cuts for shifts
  cuts[ "BTAG" ] = { "UP": cuts[ "NOMINAL" ].replace( "NJetsCSV_JetSubCalc", "NJetsCSV_JetSubCalc_bSFup" ),
                     "DN": cuts[ "NOMINAL" ].replace( "NJetsCSV_JetSubCalc", "NJetsCSV_JetSubCalc_bSFdn" ) }
  cuts[ "MISTAG" ] = { "UP": cuts[ "NOMINAL" ].replace( "NJetsCSV_JetSubCalc", "NJetsCSV_JetSubCalc_lSFup" ),
                       "DN": cuts[ "NOMINAL" ].replace( "NJetsCSV_JetSubCalc", "NJetsCSV_JetSubCalc_lSFdn" ) }
  cuts[ "TAU21" ] = { "UP": cuts[ "NOMINAL" ].replace( "NJetsWtagged", "NJetsWtagged_shifts[0]" ),
                      "DN": cuts[ "NOMINAL" ].replace( "NJetsWtagged", "NJetsWtagged_shifts[1]" ) }
  cuts[ "JMSW" ] = { "UP": cuts[ "NOMINAL" ].replace( "NJetsWtagged", "NJetsWtagged_shifts[2]" ),
                     "DN": cuts[ "NOMINAL" ].replace( "NJetsWtagged", "NJetsWtagged_shifts[3]" ) }
  cuts[ "JMRW" ] = { "UP": cuts[ "NOMINAL" ].replace( "NJetsWtagged", "NJetsWtagged_shifts[4]" ),
                     "DN": cuts[ "NOMINAL" ].replace( "NJetsWtagged", "NJetsWtagged_shifts[5]" ) }
  cuts[ "TAU21PT" ] = { "UP": cuts[ "NOMINAL" ].replace( "NJetsWtagged", "NJetsWtagged_shifts[6]" ),
                        "DN": cuts[ "NOMINAL" ].replace( "NJetsWtagged", "NJetsWtagged_shifts[7]" ) }
  cuts[ "TAU32" ] = { "UP": cuts[ "NOMINAL" ].replace( "NJetsTtagged", "NJetsTtagged_shifts[0]" ),
                      "DN": cuts[ "NOMINAL" ].replace( "NJetsTtagged", "NJetsTtagged_shifts[1]" ) }
  cuts[ "JMST" ] = { "UP": cuts[ "NOMINAL" ].replace( "NJetsTtagged", "NJetsTtagged_shifts[2]" ),
                     "DN": cuts[ "NOMINAL" ].replace( "NJetsTtagged", "NJetsTtagged_shifts[3]" ) }
  cuts[ "JMRT" ] = { "UP": cuts[ "NOMINAL" ].replace( "NJetsTtagged", "NJetsTtagged_shifts[4]" ),
                     "DN": cuts[ "NOMINAL" ].replace( "NJetsTtagged", "NJetsTtagged_shifts[5]" ) }
  cuts[ "HOTSTAT" ] = { "UP": cuts[ "NOMINAL" ].replace( "NresolvedTops1pFake", "NresolvedTops1pFake_shifts[0]" ),
                        "DN": cuts[ "NOMINAL" ].replace( "NresolvedTops1pFake", "NresolvedTops1pFake_shifts[1]" ) }
  cuts[ "HOTCSPUR" ] = { "UP": cuts[ "NOMINAL" ].replace( "NresolvedTops1pFake", "NresolvedTops1pFake_shifts[2]" ),
                         "DN": cuts[ "NOMINAL" ].replace( "NresolvedTops1pFake", "NresolvedTops1pFake_shifts[3]" ) }
  cuts[ "HOTCLOSURE" ] = { "UP": cuts[ "NOMINAL" ].replace( "NresolvedTops1pFake", "NresolvedTops1pFake_shifts[4]" ),
                           "DN": cuts[ "NOMINAL" ].replace( "NresolvedTops1pFake", "NresolvedTops1pFake_shifts[5]" ) }
    
  # declare histograms
  hists = {}
  categoryTag = "is{}nHOT{}nT{}nW{}nB{}nJ{}".format(
    category[ "LEPTON" ][0], category[ "NHOT" ][0], category[ "NT" ][0],
    category[ "NW" ][0], category[ "NB" ][0], category[ "NJ" ][0]
  )
  histTag = hist_tag( process, categoryTag )
  hists[ histTag ] = ROOT.TH1D( histTag, xLabel, len( histBins ) - 1, histBins )
  
  if doABCDNN:
    histTag = hist_tag( process, categoryTag, "ABCDNN" ) 
    hists[ histTag ] = ROOT.TH1D( histTag, xLabel, len( histBins ) - 1, histBins )
    
  if doSYST:
    for syst in config.systematics[ "MC" ]:
      if not config.systematics[ "MC" ][ syst ][0]: continue
      for shift in [ "UP", "DN" ]:
        if doABCDNN and syst.upper() in config.params[ "ABCDNN" ][ "SYSTEMATICS" ]:
          histTag = hist_tag( process, categoryTag, "ABCDNN", syst.upper() + shift )
          hists[ histTag ] = ROOT.TH1D( histTag, xLabel, len( histBins ) - 1, histBins )
        if "ABCD" not in syst.upper():
          if syst.upper() == "JEC":
            for systJEC in config.systematics[ "REDUCED JEC" ]:
              if not config.systematics[ "REDUCED JEC" ][ systJEC ]: continue
              histTag = hist_tag( process, categoryTag, "JEC" + systJEC.upper().replace( "ERA", "20" + args.year ).replace( "APV", "" ).replace( "_", "" ) + shift )
              hists[ histTag ] = ROOT.TH1D( histTag, xLabel, len( histBins ) - 1, histBins )
          else:
            histTag = hist_tag( process, categoryTag, syst.upper() + shift )
            if syst.upper() == "PREFIRE" and year not in [ "16APV", "16", "17" ]: continue
            hists[ histTag ] = ROOT.TH1D( histTag, xLabel, len( histBins ) - 1, histBins )
  if doPDF:
    for i in range( config.params[ "GENERAL" ][ "PDF RANGE" ] ):
      histTag = hist_tag( process, categoryTag, "PDF" + str(i) ) 
      hists[ histTag ] = ROOT.TH1D( histTag, xLabel, len( histBins ) - 1, histBins )
      if doABCDNN:
        histTag = hist_tag( process, categoryTag, "ABCDNN", "PDF" + str(i) )
        hists[ histTag ] = ROOT.TH1D( histTag, xLabel, len( histBins ) - 1, histBins )


  # Sumw2() tells the hist to also store the sum of squares of weights
  for histTag in hists: hists[ histTag ].Sumw2()
	
  if verbose: 
    print( ">> Applying NOMINAL weights: {}".format( mc_weights[ "NOMINAL" ] ) )
    if doABCDNN:
      print( ">> Including ABCDnn weights: {}".format( mc_weights[ "ABCDNN" ] ) )
    
    print( ">> Applying NOMINAL cuts: {}".format( cuts[ "NOMINAL" ] ) )
    if doABCDNN:
      print( ">> Applying ABCDnn cuts: {}".format( cuts[ "ABCDNN" ] ) )

  # draw histograms
  histTag = hist_tag( process, categoryTag )
  rTree[ process ].Draw( 
    "{} >> {}".format( variableName, histTag ), 
    "{} * ({})".format( mc_weights[ "NOMINAL" ], cuts[ "NOMINAL" ] ), 
    "GOFF" )

  if verbose: print( "  + NOMINAL: {} --> {}".format( rTree[ process ].GetEntries(), hists[ histTag ].Integral() ) )

  if doABCDNN:
    rTree[ process ].Draw(
      "{} >> {}".format( abcdnnName, hist_tag( process, categoryTag, "ABCDNN" ) ),
      "{} * ({})".format( mc_weights[ "ABCDNN" ], cuts[ "ABCDNN" ] ),
      "GOFF"
    )
    if verbose: print( "  + ABCDNN: {} --> {}".format( rTree[ process ].GetEntries(), hists[ hist_tag( process, categoryTag, "ABCDNN" ) ].Integral() ) )

  if process not in groups[ "DAT" ] and doSYST:
    nSyst, nSystABCDNN = 0, 0
    for syst in config.systematics[ "MC" ].keys():
      if not config.systematics[ "MC" ][ syst ][0]: continue
      for shift in [ "UP", "DN" ]:
        histTag = hist_tag( process, categoryTag, syst.upper() + shift )
        if syst.upper() in [ "PREFIRE" ] and year in [ "16APV", "16", "17" ]:
          rTree[ process ].Draw(
            "{} >> {}".format( variableName, histTag ),
            "{} * ({})".format( mc_weights[ syst.upper() ][ shift ], cuts[ "NOMINAL" ] ),
            "GOFF"
          )
          nSyst += 1
        if syst.upper() in [ "PILEUP", "PILEUPJETID", "MURFCORRD", "MUR", "MUF", "ISR", "FSR" ]:
          rTree[ process ].Draw( 
            "{} >> {}".format( variableName, histTag ), 
            "{} * ({})".format( mc_weights[ syst.upper() ][ shift ], cuts[ "NOMINAL" ] ), 
            "GOFF" 
          )
          nSyst += 1
        # hot-tagging plots
        elif ( syst.upper() in [ "HOTSTAT", "HOTCSPUR", "HOTCLOSURE" ] and category[ "NHOT" ][0] != "0p" ):
          rTree[ process ].Draw( 
            "{} >> {}".format( variableName, histTag ), 
            "{} * ({})".format( mc_weights[ "NOMINAL" ], cuts[ syst.upper() ][ shift ] ), 
            "GOFF" 
          )
          nSyst += 1
        # t-tagging plots
        elif ( syst.upper() in [ "TAU32", "JMST", "JMRT" ] and category[ "NT" ][0] != "0p" ):
          if "ttagged" in variableName.lower() or "tjet" in variableName.lower():
            shift_indx = 2*np.argwhere( np.array([ "TAU32", "JMST", "JMRT" ]) == syst.upper() )[0,0] + np.argwhere( np.array([ "UP", "DN" ]) == shift )[0,0]
            rTree[ process ].Draw( 
              "{}_shifts[{}] >> {}".format( variableName, shift_indx, histTag ), 
              "{} * ({})".format( mc_weights[ "NOMINAL" ], cuts[ syst.upper() ][ shift ] ), 
              "GOFF" 
            )
            nSyst += 1
          else: 
            rTree[ process ].Draw( 
              "{} >> {}".format( variableName, histTag ), 
              "{} * ({})".format( mc_weights[ "NOMINAL" ], cuts[ syst.upper() ][ shift ] ), 
              "GOFF" 
            )
            nSyst += 1
        # W-tagging plots
        elif ( syst in [ "TAU21", "JMSW", "JMRW", "TAU21PT" ] and category[ "NW" ][0] != "0p" ):
          if "wtagged" in variableName.lower() or "wjet" in variableName.lower():
            shift_indx = 2*np.argwhere( np.array([ "TAU21", "JMSW", "JMRW", "TAU21PT" ]) == syst.upper() )[0,0] + np.argwhere( np.array([ "UP", "DN" ]) == shift )[0,0]
            rTree[ process ].Draw( 
              "{}_shifts[{}] >> {}".format( variableName, shift_indx, histTag ), 
              "{} * ({})".format( mc_weights[ "NOMINAL" ], cuts[ syst.upper() ][ shift ] ), 
              "GOFF" 
            )
            nSyst += 1
          else: 
            rTree[ process ].Draw( 
              "{} >> {}".format( variableName, histTag ), 
              "{} * ({})".format( mc_weights[ "NOMINAL" ], cuts[ syst.upper() ][ shift ] ), 
              "GOFF" 
            )
            nSyst += 1
        # b-tagging plots
        elif ( syst.upper() in [ "LF", "LFSTATS1", "LFSTATS2", "HF", "HFSTATS1", "HFSTATS2", "CFERR1", "CFERR2" ] and category[ "NB" ][0] != "0p" ):
          rTree[ process ].Draw(
            "{} >> {}".format( variableName, histTag ), 
            "{} * ({})".format( mc_weights[ syst.upper() ][ shift ], cuts[ "NOMINAL" ] ), 
            "GOFF" 
          )
          nSyst += 1
        # process jec and jer
        elif syst.upper() in [ "JER" ]: 
          rTree[ process + syst.upper() + shift ].Draw( 
            "{} >> {}".format( variableName, histTag ), 
            "{} * ({})".format( mc_weights[ "NOMINAL" ], cuts[ "NOMINAL" ] ), 
            "GOFF" 
          )
          nSyst += 1
        elif syst.upper() in [ "JEC" ]:
          for systJEC in config.systematics[ "REDUCED JEC" ]:
            if not config.systematics[ "REDUCED JEC" ][ systJEC ]: continue
            rTree[ process + "JEC" + systJEC.upper().replace( "ERA", "20" + args.year ).replace( "APV", "" ).replace( "_", "" ) + shift ].Draw(
              "{} >> {}".format( variableName, hist_tag( process, categoryTag, "JEC" + systJEC.upper().replace( "ERA", "20" + args.year ).replace( "APV", "" ).replace( "_", "" ) + shift ) ),
              "{} * ({})".format( mc_weights[ "NOMINAL" ], cuts[ "NOMINAL" ] ),
              "GOFF" 
            )
            nSyst += 1
        elif syst.upper() in [ "TOPPT" ]:
          rTree[ process ].Draw(
            "{} >> {}".format( variableName, histTag ),
            "{} * ({})".format( mc_weights[ syst.upper() ][ shift ], cuts[ "NOMINAL" ] )
          )
          nSyst += 1
        else:
          print( "[WARN] {} turned on, but excluded for {} in traditional SF {}...".format( syst.upper() + shift, process, categoryTag ) )
        if syst.upper() in config.params[ "ABCDNN" ][ "SYSTEMATICS" ] and doABCDNN and config.systematics[ "MC" ][ syst ][0]:
          print( "[ABCDNN] Including {} for ABCDnn {} {}".format( syst.upper() + shift, process, categoryTag ) )
          if syst.upper() == "ABCDNNSAMPLE":
            rTree[ process + "JECABCDNNSAMPLE" + shift ].Draw(
              "{} >> {}".format( abcdnnName, hist_tag( process, categoryTag, "ABCDNN", syst.upper() + shift ) ),
              "{} * ({})".format( mc_weights[ "ABCDNN" ], cuts[ "ABCDNN" ] ),
              "GOFF"
            )
            nSystABCDNN += 1
          elif syst.upper() == "ABCDNNCLOSURE":
            rTree[ process ].Draw(
              "{} >> {}".format( abcdnnName + "_CLOSURE" + shift, hist_tag( process, categoryTag, "ABCDNN", syst.upper() + shift ) ),
              "{} * ({})".format( mc_weights[ "ABCDNN {}".format( syst.upper() ) ][ shift ], cuts[ "ABCDNN" ] ),
              "GOFF"
            )
            nSystABCDNN += 1
          elif syst.upper() == "ABCDNNMODEL":
            rTree[ process ].Draw(
              "{} >> {}".format( abcdnnName + "_MODEL" + shift, hist_tag( process, categoryTag, "ABCDNN", syst.upper() + shift ) ),
              "{} * ({})".format( mc_weights[ "ABCDNN {}".format( syst.upper() ) ][ shift ], cuts[ "ABCDNN" ] ),
              "GOFF"
            )
            nSystABCDNN += 1
    if verbose: print( "[DONE] Added {} systematics and {} ABCDnn systematics".format( nSyst, nSystABCDNN ) ) 
	
  if doPDF:
    for i in range( config.params[ "GENERAL" ][ "PDF RANGE" ] ):
      histTag = hist_tag( process, categoryTag, "PDF" + str(i) )
      rTree[ process ].Draw( 
        "{} >> {}".format( variableName, histTag ), 
        "pdfWeights[{}] * {} * ({})".format( i, mc_weights[ "NOMINAL" ], cuts[ "NOMINAL" ] ), 
        "GOFF" 
      )
    if verbose: print( "  + PDF" )
							
  for key in hists: hists[ key ].SetDirectory(0)
  return hists

def numTrueHist( useJES, useABCDNN ):
  def add_process( nHist, group, key, process, shift, postfix ):
    inputDir = config.inputDir[ args.year] 
    if group=="SIG":
        inputDir = config.siginputDir[ args.year ] 
    rFile = ROOT.TFile.Open( os.path.join( inputDir.replace( "step2", "step1hadds" ), shift + "/", samples.samples[ group ][ process ] + "_{}.root".format( postfix ) ) )
    #rFile = ROOT.TFile.Open( os.path.join( inputDir, shift + "/", samples.samples[ group ][ process ] + "_{}.root".format( postfix ) ) )
    print(key) 
    nHist[ key ] = rFile.Get( "NumTrueHist" ).Integral()
    rFile.Close()
    return nHist

  print( "[START] Retrieving the MC hist count" )
  nHist = {}
  for group in [ "BKG", "SIG" ]:
    for process in groups[ group ]:
      nHist = add_process( nHist, group, process, process, "nominal", "hadd" )
      for shift in [ "up", "down" ]:
        shift_ = "UP" if shift == "up" else "DN"
        if useJES:
          if config.systematics[ "MC" ][ "JEC" ][0]:
            for systJEC in config.systematics[ "REDUCED JEC" ]:
              if not config.systematics[ "REDUCED JEC" ][ systJEC ]: continue
              systJEC_ = systJEC.replace( "Era", "20" + args.year ).replace( "APV", "" )
              if systJEC.upper() == "TOTAL":
                add_process( nHist, group, "JEC" + systJEC_.upper() + shift_.upper(), process, "JEC" + shift, "hadd" )
              else:
                add_process( nHist, group, "JEC" + systJEC_.replace( "_", "" ).upper() + shift_.upper(), process, systJEC_ + shift, "hadd" )
          if config.systematics[ "MC" ][ "JER" ][0]:
            add_process( nHist, group, "JER" + shift_.upper(), process, "JER" + shift, "hadd" )
        elif useABCDNN and "ABCDNNSAMPLE" in config.params["ABCDNN"]["SYSTEMATICS"] and config.systematics[ "MC" ][ "ABCDNNSAMPLE" ][0]:
          for shift in [ "up", "down" ]:
            shift_ = "UP" if shift == "up" else "DN"
            add_process( nHist, group, "JECABCDNNSAMPLE" + shift_.upper(), process, "JEC" + shift, "ABCDnn_hadd" )
  return nHist

def make_hists( groups, group, category, nHist, useABCDNN ): 
  # only valid group arguments are DAT, SIG, BKG, TEST
  doSys = config.options[ "GENERAL" ][ "SYSTEMATICS" ] if group in [ "SIG", "BKG", "TEST" ] else False

  inputDir = config.inputDir[ args.year]
  if group=="SIG":
      inputDir = config.siginputDir[ args.year ]
  hists = {}
  for process in groups[ group ]:
    process_time = time.time()
    rFiles, rTrees = {}, {} 
    variable = args.variable
    if useABCDNN and args.variable in config.params[ "ABCDNN" ][ "TRANSFER VARIABLES" ] and process in samples.groups[ "BKG" ][ "ABCDNN" ] and contains_category( category, config.hist_bins[ "ABCDNN" ] ):
      print( "[INFO] Using ABCDnn sample for: {}".format( process ) )
      isABCDNN = True
      rFiles[ process ], rTrees[ process ] = read_tree( os.path.join( inputDir.replace( "step3", "step3_ABCDnn" ), "nominal/", samples.samples[ group ][ process ] + "_ABCDnn_hadd.root" ) ) 
    else:
      isABCDNN = False
      rFiles[ process ], rTrees[ process ] = read_tree( os.path.join( inputDir, "nominal/", samples.samples[ group ][ process ] + "_hadd.root" ) )
    if group in [ "SIG", "BKG", "TEST" ]:
      for shift in [ "up", "down" ]:
        shift_ = "UP" if shift == "up" else "DN"
        if config.systematics[ "MC" ][ "JEC" ][0]:
          if isABCDNN and process in samples.groups[ "BKG" ][ "ABCDNN" ] and config.systematics[ "MC" ][ "ABCDNNSAMPLE" ][0]:
            rFiles[ process + "JECABCDNNSAMPLE" + shift_ ], rTrees[ process + "JECABCDNNSAMPLE" + shift_ ] = read_tree( os.path.join( inputDir.replace( "step3", "step3_ABCDnn" ), "JEC" + shift, samples.samples[ group ][ process ] + "_ABCDnn_hadd.root" ) )
          for systJEC in config.systematics[ "REDUCED JEC" ]:
            systJEC_ = systJEC.replace( "Era", "20" + args.year ).replace( "APV", "" )
            if config.systematics[ "REDUCED JEC" ][ systJEC ]:
              if systJEC == "Total":
                rFiles[ process + "JEC" + systJEC_.upper() + shift_ ], rTrees[ process + "JEC" + systJEC_.upper() + shift_ ] = read_tree( os.path.join( inputDir, "JEC" + shift, samples.samples[ group ][ process ] + "_hadd.root" ) )
              else:
                rFiles[ process + "JEC" + systJEC_.upper().replace( "_", "" ) + shift_ ], rTrees[ process + "JEC" + systJEC_.upper().replace( "_", "" ) + shift_ ] = read_tree( os.path.join( inputDir, systJEC_ + shift, samples.samples[ group ][ process ] + "_hadd.root" ) )
        if config.systematics[ "MC" ][ "JER" ][0]:
          rFiles[ process + "JER" + shift_ ], rTrees[ process + "JER" + shift_ ] = read_tree( os.path.join( inputDir, "JER" + shift, samples.samples[ group ][ process ] + "_hadd.root" ) )
    hists.update( analyze( rTrees, nHist, args.year, process, variable, doSys, config.options[ "GENERAL" ][ "PDF" ], isABCDNN, category, True ) )
    print( "[OK] Added hists for {} in {:.2f} minutes".format( process, round( ( time.time() - process_time ) / 60,2 ) ) )
    del rFiles, rTrees

  if config.options[ "GENERAL" ][ "UE" ] and group in [ "UE" ]:
    for process in groups[ "UE" ]:
      process_time = time.time()
      rTree = read_tree( os.path.join( inputDir, "nominal/", samples.samples[ "BKG" ][ process ] + "_hadd.root" ) )
      hists.update( analyze( rTree, nHist, args.year, process, variable, False, config.options[ "GENERAL" ][ "PDF" ], False, category, True ) )
      print( "[OK] Added hists for {} in {:.2f} minutes".format( process, round( ( time.time() - process_time ) / 60, 2 ) ) )

  if config.options[ "GENERAL" ][ "HDAMP" ] and group in [ "HD" ]:
    for process in groups[ "HD" ]:
      process_time = time.time()
      rTree = read_tree( os.path.join( inputDir, "nominal/", samples.samples[ "BKG" ][ process ] + "_hadd.root" ) )
      hists.update( analyze( rTree, nHist, args.year, process, variable, False, config.options[ "GENERAL" ][ "PDF" ], False, category, True ) )
      print( "[OK] Added hists for {} in {:.2f} minutes".format( process, round( ( time.time() - process_time ) / 60, 2 ) ) )
  categoryDir = "is{}nHOT{}nT{}nW{}nB{}nJ{}".format( category[ "LEPTON" ][0], category[ "NHOT" ][0], category[ "NT" ][0], category[ "NW" ][0], category[ "NB" ][0], category[ "NJ" ][0] )

  if not os.path.exists( "{}/{}".format( args.subDir, categoryDir ) ): os.system( "mkdir -vp {}/{}".format( args.subDir, categoryDir ) )
  pickle.dump( hists, open( "{}/{}/{}_{}.pkl".format( args.subDir, categoryDir, group, args.variable ), "wb" ) )

def main():
  nHist = numTrueHist( config.options[ "GENERAL" ][ "SYSTEMATICS" ], config.options[ "GENERAL" ][ "ABCDNN" ] )
  if not config.options[ "GENERAL" ][ "TEST" ]:
    for group in [ "DAT", "BKG", "SIG" ]:
      group_time = time.time()
      print( "[START] Processing hists for {}".format( group ) )
      for key in category: print( "  - {}: {}".format( key, category[ key ] ) )
      make_hists( groups, group, category, nHist, config.options[ "GENERAL" ][ "ABCDNN" ] )
      print( "[DONE] Finished processing hists for {} in {} minutes".format( group, round( ( time.time() - group_time ) / 60, 2 ) ) )
  else:
    test_time = time.time() 
    print( "[START] Processing TEST hists" )
    for key in category: print( "  - {}: {}".format( key, category[ key ] ) )
    make_hists( groups, "TEST", category, nHist, config.options[ "GENERAL" ][ "ABCDNN" ] )
    print( "[DONE] Finished processing hists for TEST in {} minutes".format( round( ( time.time() - test_time ) / 60, 2 ) ) )

  print( "[DONE] Finished making hists in {}".format( round( ( time.time() - start_time ) / 60, 2 ) ) )

main()
