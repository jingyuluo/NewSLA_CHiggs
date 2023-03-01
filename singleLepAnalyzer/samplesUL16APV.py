#!/usr/bin/python

import numpy as np
import config

samples = {
  "TEST": {
    #"TTToSemiLeptonicttjj1": "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_HT0Njet0_ttjj_1",
    "Ts": "ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8",
    #"TTToSemiLeptonictt2b": "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_HT0Njet0_tt2b",
    #"TTTJ": "TTTJ_TuneCP5_13TeV-madgraph-pythia8",
    #"TTTW": "TTTW_TuneCP5_13TeV-madgraph-pythia8",
    #"DataE": "SingleElectron",
    #"TTToSemiLeptonicHT500tt2b": "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_HT500Njet9_tt2b"
  },
  "DAT": {
    "DataE": "SingleElectron",
    "DataM": "SingleMuon",
  },
  "SIG": {
    "TTTJ": "TTTJ_TuneCP5_13TeV-madgraph-pythia8",
    "TTTW": "TTTW_TuneCP5_13TeV-madgraph-pythia8",
  },
  "BKG": {
    "WW": "WW_TuneCP5_13TeV-pythia8",
    "WZ": "WZ_TuneCP5_13TeV-pythia8",
    "ZZ": "ZZ_TuneCP5_13TeV-pythia8",

    "Ts": "ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8",
    "Tt": "ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8",
    "Tbt": "ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8",
    "TtW": "ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8",
    "TbtW": "ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8",

    "TTTT": "TTTT_TuneCP5_13TeV-amcatnlo-pythia8",
    "TTHH": "TTHH_TuneCP5_13TeV-madgraph-pythia8",
    "TTWH": "TTWH_TuneCP5_13TeV-madgraph-pythia8",
    "TTWW": "TTWW_TuneCP5_13TeV-madgraph-pythia8",
    "TTWZ": "TTWZ_TuneCP5_13TeV-madgraph-pythia8",
    "TTZH": "TTZH_TuneCP5_13TeV-madgraph-pythia8",
    "TTZZ": "TTZZ_TuneCP5_13TeV-madgraph-pythia8",
    "TTWl": "TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8",
    "TTWq": "TTWJetsToQQ_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8", 
    "TTZlM10": "TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8",
    "TTZlM1to10": "TTZToLL_M-1to10_TuneCP5_13TeV-amcatnlo-pythia8",
    "TTHB": "ttHTobb_M125_TuneCP5_13TeV-powheg-pythia8",
    "TTHnoB": "ttHToNonbb_M125_TuneCP5_13TeV-powheg-pythia8"
  },
  "HD": {},
  "UE": {}
}

for DYM_HT in [ "200to400", "400to600", "600to800", "800to1200", "1200to2500", "2500toInf" ]: 
  samples[ "BKG" ][ "DYM" + DYM_HT.split( "to" )[0] ] = "DYJetsToLL_M-50_HT-{}_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8".format( DYM_HT )

for WJets_HT in [ "200To400", "400To600", "600To800", "800To1200", "1200To2500", "2500ToInf" ]: 
  samples[ "BKG" ][ "WJetsMG" + WJets_HT.split( "To" )[0] ] = "WJetsToLNu_HT-{}_TuneCP5_13TeV-madgraphMLM-pythia8".format( WJets_HT )
  
for QCD_HT in [ "200to300", "300to500", "500to700", "700to1000", "1000to1500", "1500to2000", "2000toInf" ]: 
  samples[ "BKG" ][ "QCD" + QCD_HT.split( "to" )[0] ] = "QCD_HT{}_TuneCP5_PSWeights_13TeV-madgraph-pythia8".format( QCD_HT ) 
  
shifts = {
  "TuneCP5": "",
  "TuneCP5up": "_UEUP",
  "TuneCP5down": "_UEDN",
  "hdampUP_TuneCP5": "_HDUP",
  "hdampDOWN_TuneCP5": "_HDDN"
}
  
N = {
  "TTToSemiLeptonicttjj": range(1,11)
}

for tt in [ "SemiLepton", "SemiLeptonic", "Hadronic", "2L2Nu" ]:
  for shift in shifts:
    for fs in [ "1b", "2b", "bb", "cc", "jj" ]:
      if tt == "SemiLeptonic":
        if shift == "TuneCP5":
          if fs == "jj":
            for n in N["TTToSemiLeptonicttjj"]:
              samples[ "BKG" ][ "TTTo{}{}tt{}{}".format( tt, shifts[ shift ], fs, n ) ] = "TTTo{}_{}_13TeV-powheg-pythia8_HT0Njet0_tt{}_{}".format( tt, shift, fs, n )
          else:
            samples[ "BKG" ][ "TTTo{}{}tt{}".format( tt, shifts[ shift ], fs ) ] = "TTTo{}_{}_13TeV-powheg-pythia8_HT0Njet0_tt{}".format( tt, shift, fs )
          samples[ "BKG" ][ "TTTo{}HT500{}tt{}".format( tt, shifts[ shift ], fs ) ] = "TTTo{}_{}_13TeV-powheg-pythia8_HT500Njet9_tt{}".format( tt, shift, fs )
        elif "HD" in shifts[ shift ]:
          samples[ "HD" ][ "TTTo{}{}tt{}".format( tt, shifts[ shift ], fs ) ] = "TTTo{}_{}_13TeV-powheg-pythia8_tt{}".format( tt, shift, fs )
        elif "UE" in shifts[ shift ]:
          samples[ "UE" ][ "TTTo{}{}tt{}".format( tt, shifts[ shift ], fs ) ] = "TTTo{}_{}_13TeV-powheg-pythia8_tt{}".format( tt, shift, fs )
      elif tt == "SemiLepton":
        if shift != "TuneCP5": continue
        samples[ "BKG" ][ "TTTo{}HT500tt{}".format( tt, fs ) ] = "TTTo{}_HT500Njet9_{}_13TeV-powheg-pythia8_tt{}".format( tt, shift, fs )
      else:
        if shift == "TuneCP5":
          samples[ "BKG" ][ "TTTo{}{}tt{}".format( tt, shifts[ shift ], fs ) ] = "TTTo{}_{}_13TeV-powheg-pythia8_tt{}".format( tt, shift, fs )
        elif "HD" in shifts[ shift ]:
          samples[ "HD" ][ "TTTo{}{}tt{}".format( tt, shifts[ shift ], fs ) ] = "TTTo{}_{}_13TeV-powheg-pythia8_tt{}".format( tt, shift, fs )
        elif "UE" in shifts[ shift ]:
          samples[ "UE" ][ "TTTo{}{}tt{}".format( tt, shifts[ shift ], fs ) ] = "TTTo{}_{}_13TeV-powheg-pythia8_tt{}".format( tt, shift, fs )

split = {
  "TTToSemiLeptonicttjj": [ "TTToSemiLeptonicttjj" + str(i) for i in N["TTToSemiLeptonicttjj"] ]
}

# define groups

groups = { group: {} for group in [ "BKG", "SIG", "DAT", "TEST" ] }

groups[ "TEST" ][ "PROCESS" ] = [ str( key ) for key in samples[ "TEST" ].keys() ]
 
groups[ "DAT" ][ "PROCESS" ] = [ "DataE", "DataM" ]  

groups[ "SIG" ][ "PROCESS" ] = [ "TTTW", "TTTJ" ]
  
ttbar = [ "TTToHadronic", "TTTo2L2Nu", "TTToSemiLeptonHT500", "TTToSemiLeptonicHT500", "TTToSemiLeptonic" ]
groups[ "BKG" ][ "PROCESS" ] = {
  "WJETS": [ "WJetsMG200", "WJetsMG400", "WJetsMG600", "WJetsMG800", "WJetsMG1200", "WJetsMG2500" ],
  "DYM": [ "DYM200", "DYM400", "DYM600", "DYM800", "DYM1200", "DYM2500" ],
  "QCD": [ "QCD200", "QCD300", "QCD500", "QCD700", "QCD1000", "QCD1500", "QCD2000" ],
  "VV": [ "WW", "WZ", "ZZ" ],
  "TOP": [ "Ts", "Tt", "Tbt", "TtW", "TbtW" ],
  "TTV": [ "TTWl", "TTWq", "TTZlM10", "TTZlM1to10", "TTHB", "TTHnoB" ], 
  "TTXY": [ "TTTT", "TTWW", "TTWH", "TTHH", "TTZZ", "TTWZ", "TTZH" ],
  "TTJJ": [ tt + "ttjj" for tt in ttbar if tt != "TTToSemiLeptonic" ], 
  "TTCC": [ tt + "ttcc" for tt in ttbar ],
  "TT1B": [ tt + "tt1b" for tt in ttbar ],
  "TT2B": [ tt + "tt2b" for tt in ttbar ],
  "TTBB": [ tt + "ttbb" for tt in ttbar ]
}

groups[ "BKG" ][ "PROCESS" ][ "TTJJ" ] += [ "TTToSemiLeptonicttjj" + str(num) for num in range(1,11) ]

groups[ "BKG" ][ "ALL" ] = np.concatenate( [ groups[ "BKG" ][ "PROCESS" ][ process ] for process in groups[ "BKG" ][ "PROCESS" ].keys() ] )

    
# grouped background processes
groups[ "BKG" ][ "SUPERGROUP" ] = {
  "TTNOBB": np.concatenate( [ groups[ "BKG" ][ "PROCESS" ][ process ] for process in [ "TTJJ", "TTCC", "TT1B", "TT2B" ] ] ).tolist(),
  "TTBB": groups[ "BKG" ][ "PROCESS" ][ "TTBB" ],
  "TOP": np.concatenate( [ groups[ "BKG" ][ "PROCESS" ][ process ] for process in [ "TOP", "TTV", "TTXY" ] if process not in [ "TTHB", "TTHnoB" ] ] ).tolist(),
  "EWK": np.concatenate( [ groups[ "BKG" ][ "PROCESS" ][ process ] for process in [ "WJETS", "DYM", "VV" ] ] ).tolist(),
  "QCD": [ "QCD200", "QCD300", "QCD500", "QCD700", "QCD1000", "QCD1500", "QCD2000" ],
  "TTH": [ "TTHB", "TTHnoB" ]
}

groups[ "BKG" ][ "ABCDNN" ] = []
for group in config.params[ "ABCDNN" ][ "GROUPS" ]:
  for process in groups[ "BKG" ][ "SUPERGROUP" ][ group ]:
    groups[ "BKG" ][ "ABCDNN" ].append( process )
  
groups[ "BKG" ][ "TTBAR_GROUPS" ] = {
  group: groups[ "BKG" ][ "SUPERGROUP" ][ group ] for group in [ "TTNOBB", "TTBB" ]
}
  
groups[ "BKG" ][ "TTBAR_PROCESS" ] = {
  process: groups[ "BKG" ][ "PROCESS" ][ process ] for process in [ "TTJJ", "TTCC", "TT1B", "TT2B", "TTBB" ]
}
  
groups[ "BKG" ][ "HT" ] = {
  "EWK": groups[ "BKG" ][ "SUPERGROUP" ][ "EWK" ],
  "WJETS": groups[ "BKG" ][ "PROCESS" ][ "WJETS" ],
  "QCD": groups[ "BKG" ][ "SUPERGROUP" ][ "QCD" ]
}
  
groups[ "BKG" ][ "TOPPT" ] = {
  process: np.array( [ groups[ "BKG" ][ "PROCESS" ][ process ] for process in [ "TTJJ", "TTCC", "TTBB", "TT1B", "TT2B" ] ] ).flatten().tolist()
}
