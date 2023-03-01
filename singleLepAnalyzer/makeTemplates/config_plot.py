import os, sys
sys.path.append( "../" )
import config
import ROOT

options = {
  "ALL SYSTEMATICS": True,
  "CR SYST": False,
  "REBINNED": True,
  "ABCDNN": False,
  "YIELDS": False,
  "NORM BIN WIDTH": False,
  "COMPARE SHAPES": False,
  "SCALE SIGNAL YIELD": True,
  "SCALE SIGNAL XSEC": False,
  "REAL PULL": False,
  "BLIND": False,
  "Y LOG": True,
  "SMOOTH": True,
}

params = {
  "POSTFIX TEXT": "Preliminary",
  "INCLUDE LEP": [ "E", "M", "L" ], # E,M,L
  "ERROR BAND": [ "SHAPE", "STAT", "NORM" ], # STAT, SHAPE, NORM, ALL
  "EXCLUDE SYST": [ # templates will contain some systematics that are being unused, so exclude them from the plots 
    "ABCDNNSAMPLE", "ABCDNNMODEL",
    "PDFEWK", "PDFQCD", "PDFTOP", "PDFTTBAR", "PDFTTH", "PDFSIG",
    "PSWGT", "PSWGTSIG", "PSWGTTTBAR", "PSWGTTOP", "PSWGTTTH", "PSWGTEWK", "PSWGTQCD",
    #"PILEUP", 
    #"PREFIRE",
    "MUR", "MURSIG", "MURTTBAR", "MURTOP", "MURTTH", "MUREWK", "MURQCD",
    "MUF", "MUFSIG", "MUFTTBAR", "MUFTOP", "MUFTTH", "MUFEWK", "MUFQCD",
    #"MURFSIG", "MURFTTBAR", "MURFTOP", "MURFTTH", "MURFEWK", "MURFQCD",
    "MURFCORRD", "MURFCORRDSIG", "MURFCORRDTTBAR", "MURFCORRDTOP", "MURFCORRDTTH", "MURFCORRDEWK", "MURFCORRDQCD",
    "MURF",
    "ISR",  
    "FSR",   
    #"HOTSTAT",   
    #"HOTCSPUR",
    #"HOTCLOSURE",
    #"LF", # this is fine
    #"LFSTATS2", # this one might have an issue
    #"lfstats2", # this one might have an issue
    #"HF", # this is fine
    #"hfstats1",
    #"hfstats2",
    #"cferr1",
    #"cferr2",
    #"JER", 
    #"JEC", 
  ],
  "SCALE SIGNAL YIELD": 10000,
  "DAT COLOR": ROOT.kBlack,
  "SIG COLOR": ROOT.kBlack,
  "SIG PULL COLOR": 2,
  "BKG COLORS": {
    "TT2B": ROOT.kRed + 3,
    "TT1B": ROOT.kRed - 3,
    "TTBB": ROOT.kCyan + 1, 
    "TTCC": ROOT.kRed - 5,
    "TTJJ": ROOT.kRed - 7,
    "TTNOBB": ROOT.kOrange - 2,
    "TOP": ROOT.kGreen + 1,  
    "TTH": ROOT.kBlue + 1,  
    "EWK": ROOT.kAzure + 2, 
    "QCD": ROOT.kSpring + 6,
    "TTBAR": ROOT.kOrange - 2,
    "ABCDNN": ROOT.kRed - 3,
    "ERROR": ROOT.kBlack,
  },
  "Y DIV": 0.35,
  "CANVAS": {
    "H REF": 700,
    "W REF": 800,
    "I PERIOD": 4,
    "I POSITION": 11,
  },
  "LEGEND": {
    "X1": 0.15,
    "Y1": 0.80,
    "X2": 0.50,
    "Y2": 0.88,
    "TEXT SIZE": 0.02
  }
}

for i in range( len( params[ "EXCLUDE SYST" ] ) ):
  params[ "EXCLUDE SYST" ][i] = params[ "EXCLUDE SYST" ][i].upper() 
  if options[ "SMOOTH" ]: params[ "EXCLUDE SYST" ].append( params[ "EXCLUDE SYST" ][i] + config.params[ "MODIFY BINNING" ][ "SMOOTHING ALGO" ].upper() )

params[ "CANVAS" ][ "T" ] = 0.10 * params[ "CANVAS" ][ "H REF" ] 
params[ "CANVAS" ][ "B" ] = 0.12 * params[ "CANVAS" ][ "H REF" ] if options[ "BLIND" ] else 0.35 * params[ "CANVAS" ][ "H REF" ]
params[ "CANVAS" ][ "L" ] = 0.12 * params[ "CANVAS" ][ "W REF" ]
params[ "CANVAS" ][ "R" ] = 0.04 * params[ "CANVAS" ][ "W REF" ]
params[ "CANVAS" ][ "W" ] = 1. * params[ "CANVAS" ][ "W REF" ]
params[ "CANVAS" ][ "H" ] = 1. * params[ "CANVAS" ][ "W REF" ]
params[ "CANVAS" ][ "TAG X" ] = 0.20
params[ "CANVAS" ][ "TAG Y" ] = 0.76

params[ "LATEX SIZE" ] = 0.04
