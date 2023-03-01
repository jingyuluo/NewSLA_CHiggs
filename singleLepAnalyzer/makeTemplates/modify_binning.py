#!/usr/bin/python

import os, sys, time, math, fnmatch
import numpy as np
sys.path.append( os.path.dirname( os.getcwd() ) )
from array import array
from utils import hist_tag, hist_parse
import config
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument( "-y", "--year", required = True )
parser.add_argument( "-t", "--tag", required = True )
parser.add_argument( "-r", "--region", required = True )
parser.add_argument( "-v", "--variable", required = True )
parser.add_argument( "--abcdnn", action = "store_true" )
parser.add_argument( "--verbose", action = "store_true" )
args = parser.parse_args()

import ROOT

if args.year == "16APV": 
  import samplesUL16APV as samples
elif args.year == "16":
  import samplesUL16 as samples
elif args.year == "17": 
  import samplesUL17 as samples
elif args.year == "18": 
  import samplesUL18 as samples
else: quit( "[ERR] Invalid -y (--year) argument. Quitting" )

def overflow( hist ):
  nbins = hist.GetXaxis().GetNbins()
  yields = hist.GetBinContent( nbins ) + hist.GetBinContent( nbins + 1 )
  error = math.sqrt( hist.GetBinError( nbins )**2 + hist.GetBinError( nbins + 1 )**2 )
  hist.SetBinContent( nbins, yields )
  hist.SetBinError( nbins, error )
  hist.SetBinContent( nbins + 1, 0 )
  hist.SetBinError( nbins + 1, 0 )
  
def underflow( hist ):
  nbins = hist.GetXaxis().GetNbins()
  yields = hist.GetBinContent( 0 ) + hist.GetBinContent( 1 )
  error = math.sqrt( hist.GetBinError( 0 )**2 + hist.GetBinError( 1 )**2 )
  hist.SetBinContent( 1, yields )
  hist.SetBinError( 1, error )
  hist.SetBinContent( 0, 0 )
  hist.SetBinError( 0, 0 ) 

def negative_bin_correction( hist ):
  integral = hist.Integral()
  for i in range( 0, hist.GetNbinsX() + 2 ):
    if hist.GetBinContent(i) <= 0:
      hist.SetBinContent( i, config.params[ "GENERAL" ][ "ZERO" ] )
      hist.SetBinError( i, config.params[ "GENERAL" ][ "ZERO" ] )
  if hist.Integral() != 0 and integral > 0: 
    hist.Scale( integral / hist.Integral() )


def smooth_shape( hist_n, hist_d, hist_u, syst, algo = "lowess", symmetrize = False ):
  hist_name = hist_n.GetName()
  graph_error = {
    path: {
      shift: ROOT.TGraphErrors() for shift in [ "UP", "DN" ]
    } for path in [ "IN", "OUT" ]
  }
  graph_smooth = {
    shift: ROOT.TGraphSmooth( hist_name + "_{}{}{}".format( syst, algo.upper(), shift ) ) for shift in [ "UP", "DN" ]
  }
  hist = {
    path: {
      "UP": hist_u.Clone(),
      "DN": hist_d.Clone()
    } for path in [ "IN", "OUT" ]
  }
  for shift in [ "UP", "DN" ]:
    hist[ "IN" ][ shift ].Divide( hist_n )

  for i in range( 1, hist_n.GetNbinsX() + 1 ):
    x = ( hist[ "IN" ][ "UP" ].GetBinLowEdge(i) + hist[ "IN" ][ "UP" ].GetBinLowEdge(i+1) ) / 2
    y = {
      "UP": hist[ "IN" ][ "UP" ].GetBinContent(i),
      "DN": hist[ "IN" ][ "DN" ].GetBinContent(i)
    }
    for shift in [ "UP", "DN" ]:
      graph_error[ "IN" ][ shift ].SetPoint( i - 1, x, y[ shift ] )
      if algo.upper() == "SUPER":
        graph_error[ "OUT" ][ shift ] = graph_smooth[ shift ].SmoothSuper( graph_error[ "IN" ][ shift ], "", 9, 0 )
      elif algo.upper() == "KERN":
        graph_error[ "OUT" ][ shift ] = graph_smooth[ shift ].SmoothKern( graph_error[ "IN" ][ shift ], "normal", 5.0 )
      elif algo.upper() == "LOWESS":
        graph_error[ "OUT" ][ shift ] = graph_smooth[ shift ].SmoothLowess( graph_error[ "IN" ][ shift ], "", config.params[ "MODIFY BINNING" ][ "LOWESS" ] )

  for i in range( 1, hist_n.GetNbinsX() + 1 ):
    if symmetrize:
      upShift = abs( 1 - graph_error[ "OUT" ][ "UP" ].GetY()[i-1] )
      dnShift = abs( 1 - graph_error[ "OUT" ][ "DN" ].GetY()[i-1] )
      meanShift = ( upShift + dnShift ) / 2.
    
      hist[ "OUT" ][ "UP" ].SetBinContent( i, hist_n.GetBinContent(i) * ( 1 + meanShift ) )
      hist[ "OUT" ][ "DN" ].SetBinContent( i, hist_n.GetBinContent(i) * ( 1 - meanShift ) )
    else:
      for shift in [ "UP", "DN" ]:
        hist[ "OUT" ][ shift ].SetBinContent( i, hist_n.GetBinContent(i) * graph_error[ "OUT" ][ shift ].GetY()[i-1] )
  return hist[ "OUT" ]

class ModifyTemplate():
  def __init__( self, filepath, options, params, groups, variable, doABCDNN ):
    self.filepath = filepath
    self.options = options
    self.params = params
    self.groups = groups
    self.variable = variable
    self.doABCDNN = doABCDNN
    self.rebinned = {}
    self.yields = {}
    print( "[INFO] Running ModifyTemplate() with the following options" )
    for option in self.options:
      print( "  + {}: {}".format( option, self.options[ option ] ) )
    print( "[INFO] Running ModifyTemplate() with the following parameters" )
    for param in self.params:
      print( "  + {}: {}".format( param, self.params[ param ] ) )
    self.load_histograms()
    self.get_xbins()
    self.rebin()
    
  def load_histograms( self ):
    print( "[START] Loading histograms from {}".format( self.filepath ) )
    self.rFile = { "INPUT":  ROOT.TFile( self.filepath ) }
    self.hist_names = [ hist_name.GetName() for hist_name in self.rFile[ "INPUT" ].GetListOfKeys() ]
    self.categories = list( set( hist_parse( hist_name, samples )[ "CATEGORY" ] for hist_name in self.hist_names ) )
    self.channels = list( set( category[3:] for category in self.categories ) )
    
    syst_log = { key: [] for key in [ "SIG SYST", "BKG SYST" ] }
    self.histograms = { key: {} for key in [ "BKG", "BKG SYST", "SIG", "SIG SYST", "DAT", "TOTAL BKG", "TOTAL SIG", "TOTAL DAT" ] }
    self.categories_abcdnn = [ category for category in self.categories if hist_parse( category, samples )[ "ABCDNN" ] ]
    count = { key: 0 for key in [ "BKG", "DAT", "SIG", "BKG SYST", "SIG SYST" ] } 
    for hist_name in sorted( self.hist_names ):
      parse = hist_parse( hist_name, samples ) 
      if parse[ "GROUP" ] == "DAT":
        if args.verbose: print( "   + data_obs: {}".format( hist_name ) )
        self.histograms[ "DAT" ][ hist_name ] = self.rFile[ "INPUT" ].Get( hist_name ).Clone( hist_name )
        try:
          self.histograms[ "TOTAL DAT" ][ parse[ "CATEGORY" ] ].Add( self.histograms[ "DAT" ][ hist_name ] )
        except:
          self.histograms[ "TOTAL DAT" ][ parse[ "CATEGORY" ] ] = self.histograms[ "DAT" ][ hist_name ].Clone( hist_tag( "DAT", parse[ "CATEGORY" ] ) )
        count[ "DAT" ] += 1

      elif parse[ "GROUP" ] == "BKG":
        if args.verbose and not parse[ "IS SYST" ]: print( "   + BKG: {}".format( hist_name ) )
        if parse[ "IS SYST" ]:
          self.histograms[ "BKG SYST" ][ hist_name ] = self.rFile[ "INPUT" ].Get( hist_name ).Clone( hist_name )
          if parse[ "SYST" ] not in syst_log[ "BKG SYST" ]: syst_log[ "BKG SYST" ].append( parse[ "SYST" ] )
          count[ "BKG SYST" ] += 1
        else:
          self.histograms[ "BKG" ][ hist_name ] = self.rFile[ "INPUT" ].Get( hist_name ).Clone( hist_name )
          if self.doABCDNN:
            if parse[ "CATEGORY" ] in self.categories_abcdnn:
              if "ABCDNN" in hist_name or parse[ "COMBINE" ] in config.params[ "ABCDNN" ][ "MINOR BKG" ]:
                try: self.histograms[ "TOTAL BKG" ][ parse[ "CATEGORY" ] ].Add( self.histograms[ "BKG" ][ hist_name ] )
                except: self.histograms[ "TOTAL BKG" ][ parse[ "CATEGORY" ] ] = self.histograms[ "BKG" ][ hist_name ].Clone( hist_tag( "BKG", parse[ "CATEGORY" ] ) )
            else:
              try: self.histograms[ "TOTAL BKG" ][ parse[ "CATEGORY" ] ].Add( self.histograms[ "BKG" ][ hist_name ] )
              except: self.histograms[ "TOTAL BKG" ][ parse[ "CATEGORY" ] ] = self.histograms[ "BKG" ][ hist_name ].Clone( hist_tag( "BKG", parse[ "CATEGORY" ] ) )
          else:
            try: self.histograms[ "TOTAL BKG" ][ parse[ "CATEGORY" ] ].Add( self.histograms[ "BKG" ][ hist_name ] )
            except: self.histograms[ "TOTAL BKG" ][ parse[ "CATEGORY" ] ] = self.histograms[ "BKG" ][ hist_name ].Clone( hist_tag( "BKG", parse[ "CATEGORY" ] ) )
          count[ "BKG" ] += 1

      elif parse[ "GROUP" ] == "SIG":
        if args.verbose and not parse[ "IS SYST" ]: print( "   + SIG: {}".format( hist_name ) )
        if parse[ "IS SYST" ]:
          self.histograms[ "SIG SYST" ][ hist_name ] = self.rFile[ "INPUT" ].Get( hist_name ).Clone( hist_name )
          if parse[ "SYST" ] not in syst_log[ "SIG SYST" ]: syst_log[ "SIG SYST" ].append( parse[ "SYST" ] )
          count[ "SIG SYST" ] += 1
        else:
          self.histograms[ "SIG" ][ hist_name ] = self.rFile[ "INPUT" ].Get( hist_name ).Clone( hist_name )
          try:
            self.histograms[ "TOTAL SIG" ][ parse[ "CATEGORY" ] ].Add( self.histograms[ "SIG" ][ hist_name ] )
          except:
            self.histograms[ "TOTAL SIG" ][ parse[ "CATEGORY" ] ] = self.histograms[ "SIG" ][ hist_name ].Clone( hist_tag( "SIG", parse[ "CATEGORY" ] ) )
          count[ "SIG" ] += 1
      else:
        if args.verbose: print( "  [WARN] {} does not belong to any of the groups: BKG, SIG, DAT, excluding...".format( hist_name ) )
    print( "  [DONE] Loaded {} histograms:".format( sum( [ count[ key ] for key in count ] ) ) )
    for key in count:
      print( "    + {}: {}".format( key, count[ key ] ) )
    print( "  [INFO] Total yields by category:" )
    for key in [ "TOTAL BKG", "TOTAL DAT", "SIG" ]:
      for category in sorted( self.histograms[ key ].keys() ):
        print( "    + {} > {}: {}".format( key, category, self.histograms[ key ][ category ].Integral() ) )
    print( "  [START] Creating lepton categories" )
    count = 0
    for key in self.histograms:
      print( "   o {}".format( key ) )
      hist_names = [ hist_name for hist_name in self.histograms[ key ].keys() if "isE" in hist_name ]
      for hist_name in sorted( hist_names ):
        name_lepton = hist_name.replace( "isE", "isL" )
        parse = hist_parse( name_lepton, samples )
        self.histograms[ key ][ name_lepton ] = self.histograms[ key ][ hist_name ].Clone( name_lepton )
        self.histograms[ key ][ name_lepton ].Add( self.histograms[ key ][ hist_name.replace( "isE", "isM" ) ] )
        if not parse[ "IS SYST" ]: print( "     + {}: {}".format( name_lepton, self.histograms[ key ][ name_lepton ].Integral() ) )  
        count += 1
    print( "   [DONE] Created {} lepton categories".format( count ) )

    total_count = 0
    for hist_key in self.histograms:
      for hist_name in self.histograms[ hist_key ]:
        self.histograms[ hist_key ][ hist_name ].SetDirectory(0)
        total_count += 1

    self.rFile[ "INPUT" ].Close()
    print( "[DONE] Adding {} histograms to modified Combine template".format( count, total_count ) )
  
  def get_xbins( self ): # done
    # get the new histogram bins that satisfy the requirement bin error / yield <= threshold
    print( "[START] Determining modified histogram binning" )
    self.xbins = { key: {} for key in [ "MERGED", "LIMIT", "MODIFY" ] } 
    for channel in sorted( self.channels ):
      print( "   + Channel: {}".format( channel ) )
      N_BINS = self.histograms[ "TOTAL BKG" ][ "isL" + channel ].GetNbinsX()
      self.xbins[ "MERGED" ][ channel ] = [ self.histograms[ "TOTAL BKG" ][ "isL" + channel ].GetXaxis().GetBinUpEdge( N_BINS ) ]
      bin_content = {
        key_lep: {
          key_type: {
            key_stat: 0. for key_stat in [ "YIELD", "ERROR" ]
          } for key_type in [ "TOTAL BKG", "TOTAL DAT", "TOTAL SIG" ]
        } for key_lep in [ "isE", "isM" ]
      }
      for process in config.params[ "COMBINE" ][ "SIGNALS" ]:
        for key_lep in [ "isE", "isM" ]:
          bin_content[ key_lep ][ process ] = { key_stat: 0. for key_stat in [ "YIELD", "ERROR" ] }
      N_MERGED = 0
      for i in range( 1, N_BINS + 1 ):
        N_MERGED += 1
        if self.params[ "STAT THRESHOLD" ] > 1.0:
          if N_MERGED < self.params[ "MIN MERGE" ]: 
            continue
          else:
            self.xbins[ "MERGED" ][ channel ].append( self.histograms[ "TOTAL BKG" ][ "isL" + channel ].GetXaxis().GetBinLowEdge( N_BINS + 1 - i ) )
            N_MERGED = 0
        elif self.variable in [ "NPU", "NJETS", "NBJETS", "NHOT" ]:
          if self.histograms[ "TOTAL BKG" ][ "isL" + channel ].GetBinContent( N_BINS + 1 - i ) > 0:
            self.xbins[ "MERGED" ][ channel ].append( self.histograms[ "TOTAL BKG" ][ "isL" + channel ].GetXaxis().GetBinLowEdge( N_BINS + 1 - i ) )
          N_MERGED = 0
        else:
          for key_type in [ "TOTAL BKG", "TOTAL DAT", "TOTAL SIG" ] + config.params[ "COMBINE" ][ "SIGNALS" ]:
            for key_lep in [ "isE", "isM" ]:
              if key_type in [ "TOTAL BKG", "TOTAL DAT", "TOTAL SIG" ]:
                bin_content[ key_lep ][ key_type ][ "YIELD" ] += self.histograms[ key_type ][ key_lep + channel ].GetBinContent( N_BINS + 1 - i )
                bin_content[ key_lep ][ key_type ][ "ERROR" ] += self.histograms[ key_type ][ key_lep + channel ].GetBinError( N_BINS + 1 - i )**2
              else:
                bin_content[ key_lep ][ key_type ][ "YIELD" ] += self.histograms[ "SIG" ][ hist_tag( key_type, key_lep + channel ) ].GetBinContent( N_BINS + 1 - i )
                bin_content[ key_lep ][ key_type ][ "YIELD" ] += self.histograms[ "SIG" ][ hist_tag( key_type, key_lep + channel ) ].GetBinError( N_BINS + 1 - i )**2
          if N_MERGED < self.params[ "MIN MERGE" ]: 
            continue
          else:
            bPass = False
            for key_type in [ "TOTAL BKG", "TOTAL DAT", "TOTAL SIG" ] + config.params[ "COMBINE" ][ "SIGNALS" ]:
              if not ( bin_content[ "isE" ][ key_type ][ "YIELD" ] > 0 and bin_content[ "isM" ][ key_type ][ "YIELD" ] > 0 ): bPass = True
            if bPass: continue
            ratio_e = math.sqrt( bin_content[ "isE" ][ "TOTAL BKG" ][ "ERROR" ] ) / bin_content[ "isE" ][ "TOTAL BKG" ][ "YIELD" ]
            ratio_m = math.sqrt( bin_content[ "isM" ][ "TOTAL BKG" ][ "ERROR" ] ) / bin_content[ "isM" ][ "TOTAL BKG" ][ "YIELD" ]
            ratio_sig = math.sqrt( bin_content[ "isE" ][ "TOTAL SIG" ][ "ERROR" ]**2 + bin_content[ "isM" ][ "TOTAL SIG" ][ "ERROR" ]**2  ) / ( bin_content[ "isE" ][ "TOTAL SIG" ][ "YIELD" ] + bin_content[ "isM" ][ "TOTAL SIG" ][ "YIELD" ] )
            if ratio_e <= self.params[ "STAT THRESHOLD" ] and ratio_m <= self.params[ "STAT THRESHOLD" ] and ratio_sig <= self.params[ "STAT THRESHOLD" ]:
              for key_type in [ "TOTAL BKG", "TOTAL DAT", "TOTAL SIG" ] + config.params[ "COMBINE" ][ "SIGNALS" ]:
                for key_lep in [ "isE", "isM" ]:
                  for key_stat in [ "YIELD", "ERROR" ]:
                    bin_content[ key_lep ][ key_type ][ key_stat ] = 0
                    N_MERGED = 0
              self.xbins[ "MERGED" ][ channel ].append( self.histograms[ "TOTAL BKG" ][ "isL" + channel ].GetXaxis().GetBinLowEdge( N_BINS + 1 - i ) ) 
      if self.xbins[ "MERGED" ][ channel ][-1] != self.histograms[ "TOTAL BKG" ][ "isL" + channel ].GetXaxis().GetBinLowEdge(1): 
        self.xbins[ "MERGED" ][ channel ].append( self.histograms[ "TOTAL BKG" ][ "isL" + channel ].GetXaxis().GetBinLowEdge(1) )
      if self.params[ "STAT THRESHOLD" ] <= 1.0:
        if self.variable in [ "NJETS", "NBJETS", "NPU", "NHOT" ]:
          if self.histograms[ "TOTAL BKG" ][ "isL" + channel ].GetBinContent(1) == 0.:
            del self.xbins[ "MERGED" ][ channel ][-1]
        else:
          if self.histograms[ "TOTAL BKG" ][ "isE" + channel ].GetBinContent(1) == 0. or self.histograms[ "TOTAL BKG" ][ "isM" + channel ].GetBinContent(1) == 0. or self.histograms[ "TOTAL SIG" ][ "isE" + channel ].GetBinContent(1) == 0. or self.histograms[ "TOTAL SIG" ][ "isM" + channel ].GetBinContent(1) == 0.:
            if len( self.xbins[ "MERGED" ][ channel ] ) > 2: 
              del self.xbins[ "MERGED" ][ channel ][-2]
          elif self.histograms[ "TOTAL BKG" ][ "isE" + channel ].GetBinError(1) / self.histograms[ "TOTAL BKG" ][ "isE" + channel ].GetBinContent(1) > self.params[ "STAT THRESHOLD" ] or self.histograms[ "TOTAL BKG" ][ "isM" + channel ].GetBinError(1) / self.histograms[ "TOTAL BKG" ][ "isM" + channel ].GetBinContent(1) > self.params[ "STAT THRESHOLD" ] or self.histograms[ "TOTAL SIG" ][ "isL" + channel ].GetBinError(1) / self.histograms[ "TOTAL SIG" ][ "isL" + channel ].GetBinContent(1) > self.params[ "STAT THRESHOLD" ]:
            if len( self.xbins[ "MERGED" ][ channel ] ) > 2:
              del self.xbins[ "MERGED" ][ channel ][-2]
      
      self.N_NEWBINS = len( self.xbins[ "MERGED" ][ channel ] )
      self.xbins[ "LIMIT" ][ channel ] = []
      for i in range( self.N_NEWBINS ):
        self.xbins[ "LIMIT" ][ channel ].append( self.xbins[ "MERGED" ][ channel ][ self.N_NEWBINS - 1 - i ] )
      
      self.xbins[ "LIMIT" ][ channel ][0] = max( min( config.plot_params[ "VARIABLES" ][ args.variable ][1] ), self.xbins[ "LIMIT" ][ channel ][0] )
      self.xbins[ "LIMIT" ][ channel ][-1] = min( max( config.plot_params[ "VARIABLES" ][ args.variable ][1] ), self.xbins[ "LIMIT" ][ channel ][-1] )
      for i in range( 1, len( self.xbins[ "LIMIT" ][ channel ] ) - 1 ):
        if self.xbins[ "LIMIT" ][ channel ][i] < self.xbins[ "LIMIT" ][ channel ][0] or self.xbins[ "LIMIT" ][ channel ][i] > self.xbins[ "LIMIT" ][ channel ][-1]:
          del self.xbins[ "LIMIT" ][ channel ][i]
          
      self.xbins[ "MODIFY" ][ channel ] = array( "d", self.xbins[ "LIMIT" ][ channel ] )
      print( "   >> New binning ({} bins): {}".format( 
        i, self.xbins[ "MODIFY" ][ channel ]
      ) )
        
  def rebin( self ): # done
  # merge the histogram bins using an uncertainty threshold requiring bin error / yield <= threshold
  # the merging requirements are determined in self.get_xbins()
    print( "[START] Rebinning histograms" )
    count = 0
    for hist_key in self.histograms:
      print( ">> Rebinning {}".format( hist_key ) )
      self.rebinned[ hist_key ] = {}
      for hist_name in self.histograms[ hist_key ]:
        xbins_channel = None
        for channel in self.xbins[ "MODIFY" ]:
          if channel in hist_name:
            xbins_channel = self.xbins[ "MODIFY" ][ channel ]
        self.rebinned[ hist_key ][ hist_name ] = self.histograms[ hist_key ][ hist_name ].Rebin(
          len( xbins_channel ) - 1,
          hist_name,
          xbins_channel
        )
        self.rebinned[ hist_key ][ hist_name ].SetDirectory(0)
        overflow( self.rebinned[ hist_key ][ hist_name ] )
        underflow( self.rebinned[ hist_key ][ hist_name ] )
        
        count += 1
    
    print( "[DONE] {} histograms rebinned".format( count ) )
  
  def compute_yield_stats( self ): # done
  # get the integral yield for each bin as well as the associated error
    print( "[START] Retrieving yields and errors for each histogram's bins." )
    count = 0
    for hist_key in self.rebinned:
      print( ">> Retrieving yields and errors for {}".format( hist_key ) )
      self.yields[ hist_key ] = {}
      hist_names = sorted( self.rebinned[ hist_key ].keys() )
      for hist_name in hist_names:
        parse = hist_parse( hist_name, samples )
        self.yields[ hist_key ][ hist_name ] = {
          "COUNT": self.rebinned[ hist_key ][ hist_name ].Integral(),
          "ERROR": 0
        }
        for i in range( 1, self.rebinned[ hist_key ][ hist_name ].GetXaxis().GetNbins() + 1 ):
          self.yields[ hist_key ][ hist_name ][ "ERROR" ] += self.rebinned[ hist_key ][ hist_name ].GetBinError(i)**2
        self.yields[ hist_key ][ hist_name ][ "ERROR" ] = math.sqrt( self.yields[ hist_key ][ hist_name ][ "ERROR" ] )
        if args.verbose and not parse[ "IS SYST" ]: 
          print( "   + {}: {:.2f} pm {:.2f}".format( 
            hist_name, 
            self.yields[ hist_key ][ hist_name ][ "COUNT" ],
            self.yields[ hist_key ][ hist_name ][ "ERROR" ]
          ) )
        count += 1
    print( "[DONE] Calculated yields for {} histograms".format( count ) )  
          
  def add_trigger_efficiency( self ): # done
  # specify trigger efficiencies for the single leptons
    print( "[START] Differentiating trigger efficiency histogram naming between lepton flavors" )
    count = 0
    for hist_key in [ "BKG SYST", "SIG SYST" ]:
      hist_names = self.rebinned[ hist_key ].keys()
      for hist_name in hist_names:
        parse = hist_parse( hist_name, samples )
        if parse[ "SYST" ].upper() !=  "TRIGEFF": continue
        if parse[ "CATEGORY" ].startswith( "isE" ):
          hist_name_el = self.rebinned[ hist_key ][ hist_name ].GetName().replace( "TRIGEFF", "ELTRIGGEFF" )
          self.rebinned[ hist_key ][ hist_name_el ] = self.rebinned[ hist_key ][ hist_name ].Clone( hist_name_el )
          self.rebinned[ hist_key ][ hist_name_el ].SetDirectory(0)
          count += 1
        if parse[ "CATEGORY" ].startswith( "isM" ):
          hist_name_mu = self.rebinned[ hist_key ][ hist_name ].GetName().replace( "TRIGEFF", "MUTRIGGEFF" )
          self.rebinned[ hist_key ][ hist_name_mu ] = self.rebinned[ hist_key ][ hist_name ].Clone( hist_name_mu )
          self.rebinned[ hist_key ][ hist_name_mu ].SetDirectory(0)
          count += 1
    print( "[DONE] Adjusted trigger naming for {} histograms.".format( count ) )
    
  def uncorrelate_years( self ): # done
  # differentiate the shifts by year
    print( "[START] Differentiating systematic shifts by year" )
    count = 0
    for hist_key in [ "BKG SYST", "SIG SYST" ]:
      hist_names = self.rebinned[ hist_key ].keys()
      for hist_name in hist_names:
        parse = hist_parse( hist_name, samples )
        if parse[ "IS SYST" ]:
          hist_name_new = self.rebinned[ hist_key ][ hist_name ].GetName().replace( "{}UP".format( parse[ "SYST" ] ), "{}{}UP".format( parse[ "SYST" ], args.year ) ).replace( "{}DN".format( parse[ "SYST" ] ), "{}{}DN".format( parse[ "SYST" ], args.year ) )
          self.rebinned[ hist_key ][ hist_name_new ] = self.rebinned[ hist_key ][ hist_name ].Clone( hist_name_new )
          self.rebinned[ hist_key ][ hist_name_new ].SetDirectory(0)
          count += 1
    print( "[DONE] Adjusted systematic shift names by year for {} histograms".format( count ) )
    
  def symmetrize_topPT_shift( self ): # done
  # symmetrize the up and down shifts for toppt systematic
    print( "[START] Symmetrizing the toppt systematic shifts" )
    count = 0
    for hist_key in [ "SIG SYST", "BKG SYST" ]:
      hist_names = self.rebinned[ hist_key ].keys()
      for hist_name in hist_names:
        parse = hist_parse( hist_name, samples )
        if parse[ "SYST" ] != "TOPPT" and parse[ "SHIFT" ] != "DN": continue # adjust TOPPTDN to TOPPTUP
        for i in range( 1, self.rebinned[ hist_key ][ hist_name ].GetNbinsX() + 1 ):
          self.rebinned[ hist_key ][ hist_name ].SetBinContent(
            i, 2. * self.rebinned[ hist_key.split( " " )[0] ][ hist_tag( parse[ "COMBINE" ], parse[ "CATEGORY" ] ) ].GetBinContent(i) - self.rebinned[ hist_key ][ hist_tag( parse[ "COMBINE" ], parse[ "CATEGORY" ], parse[ "SYST" ] + "UP" ) ].GetBinContent(i) 
          )
        count += 1
    print( "[DONE] Adjusted {} toppt histograms".format( count ) )

  def symmetrize_theory_shift( self ):
    print( "[START] Symmetrizing the ISR, FSR, MUR, MUF and MURFCORRD systematic shifts" )
    count = 0
    for hist_key in [ "SIG SYST", "BKG SYST" ]:
      hist_names = self.rebinned[ hist_key ].keys()
      for hist_name in hist_names:
        parse = hist_parse( hist_name, samples )
        if parse[ "SYST" ] not in [ "ISR", "FSR", "MUR", "MUF", "MURFCORRD" ] and parse[ "SHIFT" ] != "UP": continue
        for i in range( 1, self.rebinned[ hist_key ][ hist_name ].GetNbinsX() + 1 ):
          nNominal = self.rebinned[ hist_key.split( " " )[0] ][ hist_tag( parse[ "COMBINE" ], parse[ "CATEGORY" ] ) ].GetBinContent(i)
          nUp = self.rebinned[ hist_key ][ hist_tag( parse[ "COMBINE" ], parse[ "CATEGORY" ], parse[ "SYST" ] + "UP" ) ].GetBinContent(i)
          nDn = self.rebinned[ hist_key ][ hist_tag( parse[ "COMBINE" ], parse[ "CATEGORY" ], parse[ "SYST" ] + "DN" ) ].GetBinContent(i)
          scaleShift = 0.707
          shift = scaleShift * abs( nUp - nDn ) / 2.
          self.rebinned[ hist_key ][ hist_tag( parse[ "COMBINE" ], parse[ "CATEGORY" ], parse[ "SYST" ] + "UP" ) ].SetBinContent( i, nNominal + shift )
          self.rebinned[ hist_key ][ hist_tag( parse[ "COMBINE" ], parse[ "CATEGORY" ], parse[ "SYST" ] + "DN" ) ].SetBinContent( i, nNominal - shift )
        count += 1
    print( "[DONE] Adjusted {} theory histograms".format( count ) )

  def add_statistical_shapes( self ): # done
  # add shifts to the bin content for the statistical shape uncertainty
    def write_statistical_hists( category, group, i, nBB ):
      count = { key: 0 for key in [ "NEGATIVE", "ZERO" ] }
      if self.rebinned[ "TOTAL BKG" ][ category ].GetNbinsX() == 1 or group == "SIG":  
        hist_names = self.rebinned[ group ].keys()
        for hist_name in sorted( hist_names ):
          parse = hist_parse( hist_name, samples )
          if parse[ "IS SYST" ] or parse[ "CATEGORY" ] != category: continue
          yields = {
            "COUNT": self.rebinned[ group ][ hist_name ].GetBinContent(i),
            "ERROR": self.rebinned[ group ][ hist_name ].GetBinError(i)
          }
          if yields[ "COUNT" ] == 0: continue
          shift_name = { 
            shift: "{}_{}_BIN{}{}".format( parse[ "COMBINE" ], category, i, shift ) for shift in [ "UP", "DN" ] 
          }
          for shift in [ "UP", "DN" ]:
            self.rebinned[ group ][ shift_name[ shift ] ] = self.rebinned[ group ][ hist_name ].Clone( shift_name[ shift ] ) 
            if shift == "UP": self.rebinned[ group ][ shift_name[ shift ] ].SetBinContent( i, yields[ "COUNT" ] + yields[ "ERROR" ] )
            if shift == "DN": self.rebinned[ group ][ shift_name[ shift ] ].SetBinContent( i, yields[ "COUNT" ] - yields[ "ERROR" ] )
            if yields[ "COUNT" ] - yields[ "ERROR" ] < 0:
              negative_bin_correction( self.rebinned[ group ][ shift_name[ shift ] ] )
              count[ "NEGATIVE" ] += 1
            if yields[ "COUNT" ] - yields[ "ERROR" ] == 0:
              self.rebinned[ group ][ shift_name[ shift ] ].SetBinContent( i, yields[ "COUNT" ] * config.params[ "GENERAL" ][ "ZERO" ] )
              count[ "ZERO" ] += 1
            self.rebinned[ group ][ shift_name[ shift ] ].SetDirectory(0)
          nBB[ group ] += 1
      else:
        bkg_max = ""
        count_max = 0
        bkg_names = self.rebinned[ "BKG" ].keys()
        for bkg_name in bkg_names:
          parse = hist_parse( bkg_name, samples )
          if parse[ "IS SYST" ]: continue
          if count_max < self.rebinned[ "BKG" ][ bkg_name ].GetBinContent(i):
            count_max = self.rebinned[ "BKG" ][ bkg_name ].GetBinContent(i)
            bkg_max = bkg_name
        error_max = self.rebinned[ "BKG" ][ bkg_max ].GetBinError(i)
        parse = hist_parse( bkg_max, samples )
        shift_name = {
          shift: "{}_{}_BIN{}{}".format( parse[ "COMBINE" ], parse[ "CATEGORY" ], i, shift ) for shift in [ "UP", "DN" ]
        }
        count = { key: 0 for key in [ "NEGATIVE", "ZERO" ] }
        for shift in [ "UP", "DN" ]:
          self.rebinned[ "BKG" ][ shift_name[ shift ] ] = self.rebinned[ "BKG" ][ bkg_max ].Clone( shift_name[ shift ] )
          if shift == "UP": self.rebinned[ "BKG" ][ shift_name[ shift ] ].SetBinContent( i, count_max + error_max )
          if shift == "DN": self.rebinned[ "BKG" ][ shift_name[ shift ] ].SetBinContent( i, count_max - error_max )
          if count_max - error_max < 0:
            negative_bin_correction( self.rebinned[ "BKG" ][ shift_name[ shift ] ] )
            count[ "NEGATIVE" ] += 1
          if count_max - error_max == 0:
            self.rebinned[ "BKG" ][ shift_name[ shift ] ].SetBinContent( i, count_max * 0.001 )
            count[ "ZERO" ] += 1
          self.rebinned[ "BKG" ][ shift_name[ shift ] ].SetDirectory(0)
        nBB[ "BKG" ] += 1
      if args.verbose and ( count[ "NEGATIVE" ] > 0 or count[ "ZERO" ] > 0 ):
        print( "[INFO] Corrections for {}, bin {}/{}:".format( hist_name, i, self.rebinned[ group ][ hist_name ].GetNbinsX() ) )
        print( "   + Negative Correction: {}".format( count[ "NEGATIVE" ] ) )
        print( "   + Zero Correction: {}".format( count[ "ZERO" ] ) )
          
    print( "[START] Adding statistical shape systematics, excluding bins beneath {} significance:".format( self.params[ "THRESHOLD BB" ] ) )
    nBB = { "BKG": 0, "SIG": 0 }
    hist_names = self.rebinned[ "TOTAL BKG" ].keys()
    for category in self.categories:
      count = { "INCLUDE": 0, "EXCLUDE": 0 }
      for i in range( 1, self.rebinned[ "TOTAL BKG" ][ category ].GetNbinsX() + 1 ):
        if self.rebinned[ "TOTAL BKG" ][ category ].GetBinContent(i) == 0.:
          error_ratio = 0.
        else:
          error_ratio = self.rebinned[ "TOTAL BKG" ][ category ].GetBinError(i) / self.rebinned[ "TOTAL BKG" ][ category ].GetBinContent(i)
        if error_ratio <= self.params[ "THRESHOLD BB" ]: # don't include the bin shape uncertainty if it's already very small
          count[ "EXCLUDE" ] += 1
          continue
        write_statistical_hists( category, "SIG", i, nBB )
        write_statistical_hists( category, "BKG", i, nBB )
        count[ "INCLUDE" ] += 1
      if args.verbose: print( "[INFO] {}: {}/{} bins shapes included".format( category, count[ "INCLUDE" ], count[ "EXCLUDE" ] + count[ "INCLUDE" ] ) )
    print( "[DONE] {} Signal bin shapes added, {} Background bin shapes added".format( nBB[ "SIG" ], nBB[ "BKG" ] ) )
      
  def normalize_abcdnn( self ):
    print( "[START] Normalizing the ABCDnn systematic shifts" )
    count = 0
    hist_names = self.rebinned[ "BKG SYST" ].keys()
    for hist_name in hist_names:
      if not self.doABCDNN: continue
      parse = hist_parse( hist_name, samples )
      if "ABCDNN" not in parse[ "SYST" ]: continue
      count += 1
      self.rebinned[ "BKG SYST" ][ hist_name ].Scale( self.rebinned[ "BKG" ][ hist_tag( "ABCDNN", parse[ "CATEGORY" ] ) ].Integral() / self.rebinned[ "BKG SYST" ][ hist_name ].Integral() ) 

    print( "[DONE] Normalized {} ABCDnn histograms".format( count ) )


  def symmetrize_HOTclosure( self ): # done
    # make the up and down shifts of the HOTClosure systematic symmetric
    print( "[START] Symmetrizing the HOT closure systematic down shifts to match the up shifts" )
    count = 0
    for hist_key in [ "SIG SYST", "BKG SYST" ]:
      hist_names = self.rebinned[ hist_key ].keys()
      for hist_name in hist_names:
        parse = hist_parse( hist_name, samples )
        if "HOTCLOSURE" not in parse[ "SYST" ].upper() and parse[ "SHIFT" ] != "UP": continue
        count += 1
        HOT_name = {
          "NOM": hist_tag( parse[ "COMBINE" ], parse[ "CATEGORY" ] ),
          "UP": hist_tag( parse[ "COMBINE" ], parse[ "CATEGORY" ], parse[ "SYST" ] + "UP" ),
          "DN": hist_tag( parse[ "COMBINE" ], parse[ "CATEGORY" ], parse[ "SYST" ] + "DN" )
        }
        for i in range( 1, self.rebinned[ hist_key ][ hist_name ].GetNbinsX() + 1 ):
          max_shift = max(
            abs( self.rebinned[ hist_key.split( " " )[0] ][ HOT_name[ "NOM" ] ].GetBinContent(i) - self.rebinned[ hist_key ][ HOT_name[ "UP" ] ].GetBinContent(i) ),
            abs( self.rebinned[ hist_key.split( " " )[0] ][ HOT_name[ "NOM" ] ].GetBinContent(i) - self.rebinned[ hist_key ][ HOT_name[ "DN" ] ].GetBinContent(i) )
          )
          self.rebinned[ hist_key ][ HOT_name[ "UP" ] ].SetBinContent( i, self.rebinned[ hist_key.split( " " )[0] ][ HOT_name[ "NOM" ] ].GetBinContent(i) + max_shift )
          self.rebinned[ hist_key ][ HOT_name[ "DN" ] ].SetBinContent( i, self.rebinned[ hist_key.split( " " )[0] ][ HOT_name[ "NOM" ] ].GetBinContent(i) - max_shift )
    print( "[DONE] Adjusted the HOT closure systematic shift for {} histograms".format( count ) )
    
  def add_muRF_shapes( self ): 
  # adding MU R+F shape systematics
    print( "[START] Adding QCD UV Renormalization and IR Factorization Scale Factors (mu) systematic shapes" )
    count = 0
    for hist_key in [ "SIG SYST", "BKG SYST" ]:
      hist_names = self.rebinned[ hist_key ].keys()
      for hist_name in hist_names:
        parse = hist_parse( hist_name, samples )
        if "ABCDNN" in hist_name: continue
        if parse[ "SYST" ].upper() != "MUR" and parse[ "SHIFT" ] != "UP": continue 
        count += 1 
        hist_muRF = { "NOMINAL": self.rebinned[ hist_key.split( " " )[0] ][ hist_tag( parse[ "COMBINE" ], parse[ "CATEGORY" ] ) ].Clone() }
        for syst in [ "MURUP", "MURDN", "MUFUP", "MUFDN", "MURFCORRDUP", "MURFCORRDDN" ]:
          hist_muRF[ syst ] = self.rebinned[ hist_key ][ hist_tag( parse[ "COMBINE" ], parse[ "CATEGORY" ], syst ) ].Clone()
        hist_muRF[ "MURFUP" ] = self.rebinned[ hist_key.split( " " )[0] ][ hist_tag( parse[ "COMBINE" ], parse[ "CATEGORY" ] ) ].Clone()
        hist_muRF[ "MURFDN" ] = self.rebinned[ hist_key.split( " " )[0] ][ hist_tag( parse[ "COMBINE" ], parse[ "CATEGORY" ] ) ].Clone()
        for i in range( 1, hist_muRF[ "NOMINAL" ].GetNbinsX() + 1 ):
          weight_key = {
            "MAX": "NOMINAL",
            "MIN": "NOMINAL"
          }
          # only retain the largest systematic shift between muR, muF, muRFcorrd
          for key in [ "MURUP", "MURDN", "MUFUP", "MUFDN", "MURFCORRDUP", "MURFCORRDDN" ]:
            if hist_muRF[ weight_key[ "MAX" ] ].GetBinContent(i) < hist_muRF[ key ].GetBinContent(i): 
              weight_key[ "MAX" ] = key
            if hist_muRF[ weight_key[ "MIN" ] ].GetBinContent(i) > hist_muRF[ key ].GetBinContent(i): 
              weight_key[ "MIN" ] = key

          hist_muRF[ "MURFUP" ].SetBinContent( i, hist_muRF[ weight_key[ "MAX" ] ].GetBinContent(i) ) 
          hist_muRF[ "MURFDN" ].SetBinContent( i, hist_muRF[ weight_key[ "MIN" ] ].GetBinContent(i) )

        if self.options[ "NORM THEORY SIG SYST" ] and hist_key == "SIG SYST":
          hist_muRF[ "MURFUP" ].Scale( 1. / config.systematics[ "MU SF" ][ args.year ][ "UP" ] )
          hist_muRF[ "MURFDN" ].Scale( 1. / config.systematics[ "MU SF" ][ args.year ][ "DN" ] )
        if self.options[ "NORM THEORY BKG SYST" ] and hist_key == "BKG SYST":
          hist_muRF[ "MURFUP" ].Scale( hist_muRF[ "NOMINAL" ].Integral() / ( hist_muRF[ "MURFUP" ].Integral() + config.params[ "GENERAL" ][ "ZERO" ] ) )
          hist_muRF[ "MURFDN" ].Scale( hist_muRF[ "NOMINAL" ].Integral() / ( hist_muRF[ "MURFDN" ].Integral() + config.params[ "GENERAL" ][ "ZERO" ] ) )

        for shift in [ "UP", "DN" ]:
          self.rebinned[ hist_key ][ "{}_{}_MURF{}".format( parse[ "COMBINE" ], parse[ "CATEGORY" ], shift ) ] = hist_muRF[ "MURF{}".format( shift ) ].Clone( "{}_{}_MURF{}".format( parse[ "COMBINE" ], parse[ "CATEGORY" ], shift ) )
          self.rebinned[ hist_key ][ "{}_{}_MURF{}".format( parse[ "COMBINE" ], parse[ "CATEGORY" ], shift ) ].SetDirectory(0)
          if parse[ "COMBINE" ] in [ "TTNOBB", "TTBB" ]: # correlate all MURF ttbar theory systematics together
            self.rebinned[ hist_key ][ "{}_{}_MURFTTBAR{}".format( parse[ "COMBINE" ], parse[ "CATEGORY" ], shift ) ] = hist_muRF[ "MURF{}".format( shift ) ].Clone( "{}_{}_MURFTTBAR{}".format( parse[ "COMBINE" ], parse[ "CATEGORY" ], shift ) )
            self.rebinned[ hist_key ][ "{}_{}_MURFTTBAR{}".format( parse[ "COMBINE" ], parse[ "CATEGORY" ], shift ) ].SetDirectory(0)
          elif parse[ "COMBINE" ] in config.params[ "COMBINE" ][ "SIGNALS" ]: # correlate signal processes together
            self.rebinned[ hist_key ][ "{}_{}_MURFSIG{}".format( parse[ "COMBINE" ], parse[ "CATEGORY" ], shift ) ] = hist_muRF[ "MURF{}".format( shift ) ].Clone( "{}_{}_MURFSIG{}".format( parse[ "COMBINE" ], parse[ "CATEGORY" ], shift ) )
            self.rebinned[ hist_key ][ "{}_{}_MURFSIG{}".format( parse[ "COMBINE" ], parse[ "CATEGORY" ], shift ) ].SetDirectory(0)
          else: # correlate MURF theory systematics by group
            self.rebinned[ hist_key ][ "{}_{}_MURF{}{}".format( parse[ "COMBINE" ], parse[ "CATEGORY" ], parse[ "COMBINE" ], shift ) ] = hist_muRF[ "MURF{}".format( shift ) ].Clone( "{}_{}_MURF{}{}".format( parse[ "COMBINE" ], parse[ "CATEGORY" ], parse[ "COMBINE" ], shift ) )
            self.rebinned[ hist_key ][ "{}_{}_MURF{}{}".format( parse[ "COMBINE" ], parse[ "CATEGORY" ], parse[ "COMBINE" ], shift ) ].SetDirectory(0)
            
    print( "[DONE] Created {} MU Renormalization and Factorization histograms".format( count ) ) 
  
  def add_PSWeight_shapes( self ): # done
    print( "[START] Determining PS weights shape systematics" )
    count = 0
    for hist_key in [ "SIG SYST", "BKG SYST" ]:
      hist_names = self.rebinned[ hist_key ].keys()
      for hist_name in hist_names:
        parse = hist_parse( hist_name, samples )
        if parse[ "SYST" ].upper() != "ISR" and parse[ "SHIFT" ] != "UP": continue 
        if parse[ "COMBINE" ] == "ABCDNN" and ( "ISR" not in config.params[ "ABCDNN" ][ "SYSTEMATICS" ] and "FSR" not in config.params[ "ABCDNN" ][ "SYSTEMATICS" ] ): continue
        count += 1
        hist_PSWeight = { "NOMINAL": self.rebinned[ hist_key.split( " " )[0] ][ hist_tag( parse[ "COMBINE" ], parse[ "CATEGORY" ] ) ].Clone() }
        hist_PSWeight[ "PSWGTUP" ] = self.rebinned[ hist_key.split( " " )[0] ][ hist_tag( parse[ "COMBINE" ], parse[ "CATEGORY" ] ) ].Clone()
        hist_PSWeight[ "PSWGTDN" ] = self.rebinned[ hist_key.split( " " )[0] ][ hist_tag( parse[ "COMBINE" ], parse[ "CATEGORY" ] ) ].Clone()
        for syst in [ "ISR", "FSR" ]:
          for shift in [ "UP", "DN" ]:
            hist_PSWeight[ syst + shift ] = self.rebinned[ hist_key ][ hist_tag( parse[ "COMBINE" ], parse[ "CATEGORY" ], syst + shift ) ].Clone()

        for i in range( 1, hist_PSWeight[ "NOMINAL" ].GetNbinsX() + 1 ):
          weight_key = {
            "MAX": "NOMINAL",
            "MIN": "NOMINAL"
          }
          weight_limit = {
            "MAX": hist_PSWeight[ "NOMINAL" ].GetBinContent(i),
            "MIN": hist_PSWeight[ "NOMINAL" ].GetBinContent(i)
          }
          weight_error = {
            "MAX": hist_PSWeight[ "NOMINAL" ].GetBinError(i),
            "MIN": hist_PSWeight[ "NOMINAL" ].GetBinError(i)
          }

          # only retain the largest systematic shift between ISR, FSR 
          for key in hist_PSWeight: 
            if hist_PSWeight[ key ].GetBinContent(i) > weight_limit[ "MAX" ]:
              weight_limit[ "MAX" ] = hist_PSWeight[ key ].GetBinContent(i)
              weight_error[ "MAX" ] = hist_PSWeight[ key ].GetBinError(i)
              weight_key[ "MAX" ] = key
            if hist_PSWeight[ key ].GetBinContent(i) < weight_limit[ "MIN" ]:
              weight_limit[ "MIN" ] = hist_PSWeight[ key ].GetBinContent(i)
              weight_error[ "MIN" ] = hist_PSWeight[ key ].GetBinError(i)
              weight_key[ "MIN" ] = key
        
          # in-case symmetrization is needed for PSWGTUP:
          hist_PSWeight[ "PSWGTUP" ].SetBinContent( i, weight_limit[ "MAX" ] )
          hist_PSWeight[ "PSWGTUP" ].SetBinError( i, weight_error[ "MAX" ] )
          hist_PSWeight[ "PSWGTDN" ].SetBinContent( i, weight_limit[ "MIN" ] )
          hist_PSWeight[ "PSWGTDN" ].SetBinError( i, weight_error[ "MIN" ] )
        
        if self.options[ "NORM THEORY {}".format( hist_key ) ]:
          for shift in [ "UP", "DN" ]:
            for syst in [ "PSWGT", "ISR", "FSR" ]:
              hist_PSWeight[ syst + shift ].Scale( hist_PSWeight[ "NOMINAL" ].Integral() / ( hist_PSWeight[ syst + shift ].Integral() + config.params[ "GENERAL" ][ "ZERO" ] ) )
        
        for syst in [ "PSWGT", "ISR", "FSR" ]:
          for shift in [ "UP", "DN" ]:
            self.rebinned[ hist_key ][ "{}_{}_{}".format( parse[ "COMBINE" ], parse[ "CATEGORY" ], syst + shift ) ] = hist_PSWeight[ syst + shift ].Clone( "{}_{}_{}".format( parse[ "COMBINE" ], parse[ "CATEGORY" ], syst + shift ) )
            self.rebinned[ hist_key ][ "{}_{}_{}".format( parse[ "COMBINE" ], parse[ "CATEGORY" ], syst + shift ) ].SetDirectory(0)
            if parse[ "COMBINE" ] in [ "TTNOBB", "TTBB" ]:
              self.rebinned[ hist_key ][ "{}_{}_{}TTBAR{}".format( parse[ "COMBINE" ], parse[ "CATEGORY" ], syst, shift ) ] = hist_PSWeight[ syst + shift ].Clone( "{}_{}_{}TTBAR{}".format( parse[ "COMBINE" ], parse[ "CATEGORY" ], syst, shift ) )
              self.rebinned[ hist_key ][ "{}_{}_{}TTBAR{}".format( parse[ "COMBINE" ], parse[ "CATEGORY" ], syst, shift ) ].SetDirectory(0)
            elif parse[ "COMBINE" ] in config.params[ "COMBINE" ][ "SIGNALS" ]:
              self.rebinned[ hist_key ][ "{}_{}_{}SIG{}".format( parse[ "COMBINE" ], parse[ "CATEGORY" ], syst, shift ) ] = hist_PSWeight[ syst + shift ].Clone( "{}_{}_{}SIG{}".format( parse[ "COMBINE" ], parse[ "CATEGORY" ], syst, shift ) )
              self.rebinned[ hist_key ][ "{}_{}_{}SIG{}".format( parse[ "COMBINE" ], parse[ "CATEGORY" ], syst, shift ) ].SetDirectory(0)
            else:
              self.rebinned[ hist_key ][ "{}_{}_{}{}{}".format( parse[ "COMBINE" ], parse[ "CATEGORY" ], syst, parse[ "COMBINE" ].upper(), shift ) ] = hist_PSWeight[ syst + shift ].Clone( "{}_{}_{}{}{}".format( parse[ "COMBINE" ], parse[ "CATEGORY" ], syst, parse[ "COMBINE" ], shift ) )
              self.rebinned[ hist_key ][ "{}_{}_{}{}{}".format( parse[ "COMBINE" ], parse[ "CATEGORY" ], syst, parse[ "COMBINE" ].upper(), shift ) ].SetDirectory(0)
             
    print( "[DONE] Created {} PS Weight histograms".format( count ) )
          
  def add_PDF_shapes( self ): 
    print( "[START] Determining PDF shape systematics" )
    count = 0
    for hist_key in [ "SIG SYST", "BKG SYST" ]:
      hist_names = self.rebinned[ hist_key ].keys()
      donePDF = []
      for hist_name in hist_names:
        parse = hist_parse( hist_name, samples )
        if parse[ "SYST" ] != "PDF": continue
        thisPDF = hist_tag( parse[ "COMBINE" ], parse[ "CATEGORY" ] )
        if thisPDF not in donePDF:
          donePDF.append( thisPDF )
          count += 1
        else:
          continue
        hist_PDF = { 
          "NOMINAL": self.rebinned[ hist_key.split( " " )[0] ][ hist_tag( parse[ "COMBINE" ], parse[ "CATEGORY" ] ) ].Clone(),
          "PDFUP": self.rebinned[ hist_key.split( " " )[0] ][ hist_tag( parse[ "COMBINE" ], parse[ "CATEGORY" ] ) ].Clone(),
          "PDFDN": self.rebinned[ hist_key.split( " " )[0] ][ hist_tag( parse[ "COMBINE" ], parse[ "CATEGORY" ] ) ].Clone()
        }
        for i in range( config.params[ "GENERAL" ][ "PDF RANGE" ] ):
          hist_PDF[ "PDF{}".format(i) ] = self.rebinned[ hist_key ][ hist_tag( parse[ "COMBINE" ], parse[ "CATEGORY" ], "PDF" + str(i) ) ].Clone( "PDF{}".format(i) )
        for i in range( 1, hist_PDF[ "NOMINAL" ].GetNbinsX() + 1 ):
          weight_limit = {
            "MAX": hist_PDF[ "PDF0" ].GetBinContent(i),
            "MIN": hist_PDF[ "PDF0" ].GetBinContent(i)
          }
          weight_error = {
            "MAX": hist_PDF[ "PDF0" ].GetBinError(i),
            "MIN": hist_PDF[ "PDF0" ].GetBinError(i)
          }
          for j in range( config.params[ "GENERAL" ][ "PDF RANGE" ] ):
            if hist_PDF[ "PDF{}".format(j) ].GetBinContent(i) > weight_limit[ "MAX" ]:
              weight_limit[ "MAX" ] = hist_PDF[ "PDF{}".format(j) ].GetBinContent(i)
              weight_error[ "MAX" ] = hist_PDF[ "PDF{}".format(j) ].GetBinError(i)
            if hist_PDF[ "PDF{}".format(j) ].GetBinContent(i) < weight_limit[ "MIN" ]:
              weight_limit[ "MIN" ] = hist_PDF[ "PDF{}".format(j) ].GetBinContent(i)
              weight_error[ "MIN" ] = hist_PDF[ "PDF{}".format(j) ].GetBinError(i)
          
          hist_PDF[ "PDFUP" ].SetBinContent( i, weight_limit[ "MAX" ] )
          hist_PDF[ "PDFUP" ].SetBinError( i, weight_error[ "MAX" ] )
          hist_PDF[ "PDFDN" ].SetBinContent( i, weight_limit[ "MIN" ] )
          hist_PDF[ "PDFDN" ].SetBinError( i, weight_error[ "MIN" ] )

        if hist_key == "SIG SYST" and not self.options[ "NORM THEORY SIG SYST" ]:
          hist_PDF[ "PDFUP" ].Scale( 1. / config.systematics[ "PDF SF" ][ args.year ][ "UP" ] )
          hist_PDF[ "PDFDN" ].Scale( 1. / config.systematics[ "PDF SF" ][ args.year ][ "DN" ] )
        elif hist_key == "BKG SYST" and self.options[ "NORM THEORY BKG SYST" ]:
          hist_PDF[ "PDFUP" ].Scale( hist_PDF[ "NOMINAL" ].Integral() / ( hist_PDF[ "PDFUP" ].Integral() + config.params[ "GENERAL" ][ "ZERO" ] ) )
          hist_PDF[ "PDFDN" ].Scale( hist_PDF[ "NOMINAL" ].Integral() / ( hist_PDF[ "PDFDN" ].Integral() + config.params[ "GENERAL" ][ "ZERO" ] ) )
        
        for shift in [ "UP", "DN" ]:
          self.rebinned[ hist_key ][ "{}_{}_PDF{}".format( parse[ "COMBINE" ], parse[ "CATEGORY" ], shift ) ] = hist_PDF[ "PDF{}".format( shift ) ].Clone( "{}_{}_PDF{}".format( parse[ "COMBINE" ], parse[ "CATEGORY" ], shift ) )
          self.rebinned[ hist_key ][ "{}_{}_PDF{}".format( parse[ "COMBINE" ], parse[ "CATEGORY" ], shift ) ].SetDirectory(0)
          if parse[ "COMBINE" ] in [ "TTNOBB", "TTBB" ]:
            self.rebinned[ hist_key ][ "{}_{}_PDFTTBAR{}".format( parse[ "COMBINE" ], parse[ "CATEGORY" ], shift ) ] = hist_PDF[ "PDF{}".format( shift ) ].Clone( "{}_{}_PDFTTBAR{}".format( parse[ "COMBINE" ], parse[ "CATEGORY" ], shift ) )
            self.rebinned[ hist_key ][ "{}_{}_PDFTTBAR{}".format( parse[ "COMBINE" ], parse[ "CATEGORY" ], shift ) ].SetDirectory(0)
          elif parse[ "COMBINE" ] in config.params[ "COMBINE" ][ "SIGNALS" ]:
            self.rebinned[ hist_key ][ "{}_{}_PDFSIG{}".format( parse[ "COMBINE" ], parse[ "CATEGORY" ], shift ) ] = hist_PDF[ "PDF{}".format( shift ) ].Clone( "{}_{}_PDFSIG{}".format( parse[ "COMBINE" ], parse[ "CATEGORY" ], shift ) )
            self.rebinned[ hist_key ][ "{}_{}_PDFSIG{}".format( parse[ "COMBINE" ], parse[ "CATEGORY" ], shift ) ].SetDirectory(0)
          else:
            self.rebinned[ hist_key ][ "{}_{}_PDF{}{}".format( parse[ "COMBINE" ], parse[ "CATEGORY" ], parse[ "COMBINE" ], shift ) ] = hist_PDF[ "PDF{}".format( shift ) ].Clone( "{}_{}_PDF{}{}".format( parse[ "COMBINE" ], parse[ "CATEGORY" ], parse[ "COMBINE" ], shift ) )
            self.rebinned[ hist_key ][ "{}_{}_PDF{}{}".format( parse[ "COMBINE" ], parse[ "CATEGORY" ], parse[ "COMBINE" ], shift ) ].SetDirectory(0)

    print( "[DONE] Determine {} PDF shape systematic histograms".format( count ) )
          
  def add_smooth_shapes( self ): # done
    print( "[START] Smoothing systematic shapes using {} smoothing".format( self.params[ "SMOOTHING ALGO" ] ) )
    sTime = time.time()
    count = 0
    symmetrize_list = [ nSyst.upper() for nSyst in config.systematics[ "MC" ] if config.systematics[ "MC" ][ nSyst ][1] ]
    for hist_key in [ "BKG SYST", "SIG SYST" ]:
      hist_names = self.rebinned[ hist_key ].keys()
      for hist_name in hist_names:
        parse = hist_parse( hist_name, samples )
        if parse[ "SHIFT" ] != "UP": continue
        try:
          hist_syst = {
            "NOMINAL": self.rebinned[ hist_key.split( " " )[0] ][ hist_tag( parse[ "COMBINE" ], parse[ "CATEGORY" ] ) ].Clone(),
            "UP": self.rebinned[ hist_key ][ hist_tag( parse[ "COMBINE" ], parse[ "CATEGORY" ], parse[ "SYST" ] + "UP" ) ].Clone(),
            "DN": self.rebinned[ hist_key ][ hist_tag( parse[ "COMBINE" ], parse[ "CATEGORY" ], parse[ "SYST" ] + "DN" ) ].Clone()
          }
        except:
          continue
        count += 1
        bSymmetrize = True if parse[ "SYST" ].upper() in symmetrize_list else False
        smooth_hist = smooth_shape( hist_syst[ "NOMINAL" ], hist_syst[ "DN" ], hist_syst[ "UP" ], syst = parse[ "SYST" ], algo = self.params[ "SMOOTHING ALGO" ], symmetrize = bSymmetrize )
        for shift in [ "UP", "DN" ]:
          smooth_name = hist_tag( parse[ "COMBINE" ], parse[ "CATEGORY" ], parse[ "SYST" ] + self.params[ "SMOOTHING ALGO" ].upper() + shift ) 
          self.rebinned[ hist_key ][ smooth_name ] = smooth_hist[ shift ].Clone( smooth_name )
          self.rebinned[ hist_key ][ smooth_name ].SetDirectory(0)
    print( "[DONE] Added {} smoothed systematic histograms".format( count ) )
      
  def write_combine( self ):
    self.outpath = self.filepath.replace( ".root", "_rebinned_merge{}_stat{}.root".format( self.params[ "MIN MERGE" ], str( self.params[ "STAT THRESHOLD" ] ).replace( ".", "p" ) ) )
    print( "[START] Storing modified histograms in {}".format( self.outpath ) )
    self.rFile[ "OUTPUT" ] = ROOT.TFile( self.outpath, "RECREATE" )
    count = 0
    for hist_key in self.rebinned:
      if "TOTAL" in hist_key: continue
      hist_names = self.rebinned[ hist_key ].keys()
      print( ">> Loading {}".format( hist_key ) )
      for hist_name in hist_names:
        if "PDF" in hist_name and not ( hist_name.endswith( "UP" ) or hist_name.endswith( "DN" ) ): continue
        negative_bin_correction( self.rebinned[ hist_key ][ hist_name ] )
        self.rebinned[ hist_key ][ hist_name ].Write() # for plotting
        count += 1
        parse = hist_parse( hist_name, samples )
        if parse[ "SYST" ] == "ABCDNN" and not hist_name.startswith( "ABCDNN" ): continue
        shift = "Down" if hist_name.endswith( "DN" ) else "Up"
        combine_tag = "NOMINAL" if not parse[ "IS SYST" ] else parse[ "SYST" ] + shift
        combine_name = hist_tag( parse[ "COMBINE" ], parse[ "CATEGORY" ], combine_tag )
        self.rebinned[ hist_key ][ combine_name ] = self.rebinned[ hist_key ][ hist_name ].Clone( combine_name )
        negative_bin_correction( self.rebinned[ hist_key ][ combine_name ] )
        self.rebinned[ hist_key ][ combine_name ].Write()
        count += 1
    print( "[DONE] {} histograms written to Combine template.".format( count ) )
    self.rFile[ "OUTPUT" ].Close()
    
def print_tables():
  table = []
  table = yield_tables( table )
  table = systematic_tables( table )
  summary_templates()
  
def main():
  # parse options and parameters
  template_prefix = config.region_prefix[ args.region ]
  templateDir = os.path.join( os.getcwd(), "{}_UL{}_{}".format( template_prefix, args.year, args.tag ) )
  
  params = config.params[ "MODIFY BINNING" ].copy()
  options = config.options[ "MODIFY BINNING" ].copy()
  if args.region == "BASELINE":
    if "NJET" in args.variable.upper() or "NBJET" in args.variable.upper() or "NHOT" in args.variable.upper() or args.variable.upper() == "NPU":
      print( "   > MIN MERGE ({}): {} --> 1".format( args.variable, params[ "MIN MERGE" ] ) )
      params[ "MIN MERGE" ] = 1
  else:
    if not config.options[ "GENERAL" ][ "FINAL ANALYSIS" ]:
      print( "[WARN] Running {} region, turning on blinding".format( args.region ) )
      options[ "BLIND" ] = True
    
  file_name = "template_combine_{}_UL{}.root".format( args.variable, args.year ) 
  file_path = os.path.join( templateDir, file_name )

  # default rebin/merge histograms
  template = ModifyTemplate( file_path, options, params, samples.groups, args.variable, args.abcdnn )  
  
  ## handling systematics
  if options[ "SYMM TOP PT" ]:
    template.symmetrize_topPT_shift()
  if options[ "TRIGGER EFFICIENCY" ]:
    template.add_trigger_efficiency()
  if options[ "SHAPE STAT" ]:
    template.add_statistical_shapes()
  if options[ "SYMM HOTCLOSURE" ]:
    template.symmetrize_HOTclosure()
  if options[ "SYMM THEORY" ]:
    template.symmetrize_theory_shift()
  if options[ "MURF SHAPES" ]:
    template.add_muRF_shapes()
  if options[ "PS WEIGHTS" ]:
    template.add_PSWeight_shapes()
  if options[ "PDF" ]:
    template.add_PDF_shapes()
  if options[ "NORM ABCDNN" ]:
    template.normalize_abcdnn()
  if options[ "SMOOTH" ]:
    template.add_smooth_shapes()
  if options[ "UNCORRELATE YEARS" ]:
    template.uncorrelate_years()
   
  # calculate yields
  #template.compute_yield_stats()
    
  #print_tables( template.table )
        
  template.write_combine()
  
main()
