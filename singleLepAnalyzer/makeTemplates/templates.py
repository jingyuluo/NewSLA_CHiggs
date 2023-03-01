#!/usr/bin/python

import os, sys, json, time, math, datetime, pickle, itertools
sys.path.append( "../" )
sys.path.append( "../../" )
from argparse import ArgumentParser
import numpy as np
from array import array
from utils import hist_parse, hist_tag
import config

parser = ArgumentParser()
parser.add_argument( "-y", "--year",   required = True )
parser.add_argument( "-t", "--tag", required = True )
parser.add_argument( "-v", "--variables", nargs = "+", required = True )
parser.add_argument( "-r", "--region", required = True )
parser.add_argument( "--verbose", action = "store_true" )
args = parser.parse_args()

# parse options
doABCDNN = config.options[ "GENERAL" ][ "ABCDNN" ]
if doABCDNN:
  print( "[INFO] Running ABCDnn templates" )
if args.region not in list( config.region_prefix.keys() ): quit( "[ERR] Invalid region argument used. Quiting..." )
if args.year == "16APV":
  import samplesUL16APV as samples
elif args.year == "16":
  import samplesUL16 as samples
elif args.year == "17":
  import samplesUL17 as samples
elif args.year == "18":
  import samplesUL18 as samples
else:
  quit( "[ERR] Invalid -y (--year) argument used. Quitting" )

from ROOT import gROOT, TFile, TH1F, Double

gROOT.SetBatch(1)

def get_categories( directory ):
  categories = [ directory for directory in os.walk( directory ).next()[1] if directory.startswith( "isE" ) or directory.startswith( "isM" ) ]
  return categories
  
def load_histograms( variable, templateDir, categories ): 
  print( "[START] Loading histograms from {} for {}".format( templateDir, variable ) )
  sTime = time.time()
  hists =  {}
  
  for category in categories:
    if args.verbose: print( "  >> Loading category: {}".format( category ) )
    categoryDir = os.path.join( templateDir, category )
    hist_keys = [ filename.split( "_" )[0] for filename in os.listdir( categoryDir ) if filename.endswith( ".pkl" ) ]
    for hist_key in hist_keys:
      if hist_key == "TEST" and not config.options[ "GENERAL" ][ "TEST" ]: continue
      if hist_key not in hists: hists[ hist_key ] = {}
      hists[ hist_key ].update( pickle.load( open( os.path.join( categoryDir, "{}_{}.pkl".format( hist_key, variable ) ), "rb" ) ) ) 
  count = 0
  for hist_key in hists:
    count += len( hists[ hist_key ].keys() )

  print( "[DONE] Finished loading {} histograms in {:.2f} minutes".format( count, ( time.time() - sTime ) / 60 ) )
  return hists
  
def clean_histograms( hists, hist_key, scale, rebin ):
  def scale_luminosity( hists_, hist_key, scale ):
    print( "  [START] Scaling {} MC luminosity by factor: {}".format( hist_key, scale ) )
    count = 0
    for hist_name in hists_[ hist_key ]:
      parse = hist_parse( hist_name, samples )
      if parse[ "GROUP" ] in [ "BKG", "SIG" ]:
        if args.verbose: print( "  + {}".format( hist_name ) )
        hists_[ hist_name ].Scale( scale )
        count += 1
    print( "  [DONE] {} histograms scaled".format( count ) )
    return hists_
  
  def rebinning( hists_, hist_key, rebin ):
    print( "  [START] Re-binning {} histogram bins by: {}".format( hist_key, rebin ) )
    count = 0
    for hist_name in hists_[ hist_key ]:
      if args.verbose: print( "  + {}".format( hist_name ) )
      hists_[ hist_key ][ hist_name ].Rebin( rebin )
      count += 1
    print( "  [DONE] Re-binned {} histograms".format( count ) )
    return hists_
  
  def negative_correction( hists_, hist_key ):
    def function( hist_ ):
      change = False
      integral = hist_.Integral()
      for i in range( hist_.GetNbinsX() + 2 ):
        if hist_.GetBinContent( i ) < 0:
          hist_.SetBinContent( i, 0 )
          hist_.SetBinError( i, 0 )
          change = True
      if hist_.Integral() != 0 and integral > 0: hist_.Scale( integral / hist_.Integral() )
      return hist_, change
      
    print( "  [START] Correcting negative {} histogram bins".format( hist_key ) )
    count = 0
    for hist_name in hists_[ hist_key ]:
      parse = hist_parse( hist_name, samples )
      if parse[ "GROUP" ] in [ "SIG", "BKG" ]:
        hists_[ hist_key ][ hist_name ], change = function( hists_[ hist_key ][ hist_name ] )
        if change: 
          count += 1
    print( "  [DONE] Corrected {} negative bins".format( count ) )
    return hists_
  
  def bin_correction( hists_, hist_key ):
    def overflow( hist_ ):
      n = hist_.GetNbinsX()
      content_over = hist_.GetBinContent( n ) + hist_.GetBinContent( n + 1 )
      error_over = math.sqrt( hist_.GetBinError( n )**2 + hist_.GetBinError( n + 1 )**2 )
      hist_.SetBinContent( n, content_over )
      hist_.SetBinError( n, error_over )
      hist_.SetBinContent( n + 1, 0 )
      hist_.SetBinError( n + 1, 0 )
      return hist_

    def underflow( hist_ ):
      content_under = hist_.GetBinContent( 1 ) + hist_.GetBinContent( 0 )
      error_under = math.sqrt( hist_.GetBinError( 1 )**2 + hist_.GetBinError( 0 )**2 )
      hist_.SetBinContent( 1, content_under )
      hist_.SetBinError( 1, error_under )
      hist_.SetBinContent( 0, 0 )
      hist_.SetBinError( 0, 0 )
      return hist_
    
    if args.verbose: print( "  [START] Correcting {} over/under-flow bins".format( hist_key ) )
    for hist_name in hists_[ hist_key ]:
      hists_[ hist_key ][ hist_name ] = overflow( hists_[ hist_key ][ hist_name ] )
      hists_[ hist_key ][ hist_name ] = underflow( hists_[ hist_key ][ hist_name ] )
    
    print( "  [DONE]" )
    return hists_
 
  sTime = time.time()
  print( "[START] Cleaning {} histograms".format( hist_key ) )
  if scale != 1.: hists = scale_luminosity( hists, hist_key, scale )
  if rebin > 0: hists = rebinning( hists, hist_key, rebin )
  hists = negative_correction( hists, hist_key )
  hists = bin_correction( hists, hist_key )
  print( "[DONE] Finished cleaning histograms in {:.2f} minutes".format( ( time.time() - sTime ) / 60. ) )
  return hists
  
def combine_histograms( hists, variable, categories, groups, doABCDNN ):
  def scale_ttbar( hists_, scaleHF, scaleLF ):
    print( "  [START] Scaling CMB ttbb histograms by a factor of {:.3f}".format( scaleHF ) )
    count = 0
    for category in sorted( categories ):
      if hist_tag( "TTBB", category ) not in hists_[ "CMB" ].keys(): continue
      N = {
        "TTBB": hists[ "CMB" ][ hist_tag( "TTBB", category ) ].Integral(),
        "TTNOBB": hists[ "CMB" ][ hist_tag( "TTNOBB", category ) ].Integral()
      }
      if scaleLF < 0:
        try: scaleLF = max( 0, 1. + ( 1. - scaleHF ) * ( N[ "TTBB" ] / N[ "TTNOBB" ] ) )
        except ZeroDivisionError: scaleLF = 1.
      if args.verbose: print( "     + {} ({:.3f}): {} --> {}:".format( category, scaleLF, N[ "TTNOBB" ], N[ "TTNOBB" ] * scaleLF ) )

      hists_[ "CMB" ][ hist_tag( "TTBB", category ) ].Scale( scaleHF )
      hists_[ "CMB" ][ hist_tag( "TTNOBB", category ) ].Scale( scaleLF )
      count += 1

      if config.options[ "GENERAL" ][ "SYSTEMATICS" ]:
        for syst in config.systematics[ "MC" ].keys():
          if args.year == "18" and syst.upper() == "PREFIRE": continue
          if not config.systematics[ "MC" ][ syst ][0] or "ABCD" in syst: continue
          if syst.upper() == "HD" and not config.options[ "GENERAL" ][ "HDAMP" ]: continue
          if syst.upper() == "UE" and not config.options[ "GENERAL" ][ "UE" ]: continue
          for shift in [ "UP", "DN" ]:
            if syst == "JEC":
              for systJEC in config.systematics[ "REDUCED JEC" ]:
                if not config.systematics[ "REDUCED JEC" ][ systJEC ]: continue
                systJEC_ = "JEC" + systJEC.replace( "Era", "20" + args.year ).replace( "APV", "" ).replace( "_", "" )
                hists_[ "CMB" ][ hist_tag( "TTBB", category, systJEC_.upper() + shift ) ].Scale( scaleHF )
                hists_[ "CMB" ][ hist_tag( "TTNOBB", category, systJEC_.upper() + shift ) ].Scale( scaleLF )
                count += 1
            else:
              hists_[ "CMB" ][ hist_tag( "TTBB", category, syst.upper() + shift ) ].Scale( scaleHF )
              hists_[ "CMB" ][ hist_tag( "TTNOBB", category, syst.upper() + shift ) ].Scale( scaleLF )
              count += 1
      if config.options[ "GENERAL" ][ "PDF" ]:
        for i in range( config.params[ "GENERAL" ][ "PDF RANGE" ] ):
          hists_[ "CMB" ][ hist_tag( "TTBB", category, "PDF" + str(i) ) ].Scale( scaleHF )
          hists_[ "CMB" ][ hist_tag( "TTNOBB", category, "PDF" + str(i) ) ].Scale( scaleLF )
          count += 1
    print( "  [DONE] Scaled {} ttbb (ttnobb) histograms by {:.3f}".format( count, scaleHF ) )
    return hists_
    
  def set_zero( hists_, categories, zero ):
    print( "  [START] Setting 0 bins to be non-trivial ({:.3e}) in histograms".format( zero ) )
    count = 0
    for hist_key in hists_:
      for hist_name in hists_[ hist_key ]:
        parse = hist_parse( hist_name, samples )
        if parse[ "GROUP" ] == "DAT": continue
        if hists_[ hist_key ][ hist_name ].Integral() == 0: 
          hists_[ hist_key ][ hist_name ].SetBinContent( 1, zero )
          count += 1
    print( "  [DONE] Set {} 0 bins to be non-trivial".format( count ) )
    return hists_

  print( "[START] Consolidating histograms by Higgs Combine grouping" )
  sTime = time.time()
  count = {}
  hists[ "CMB" ] = {}
  for hist_key in hists:
    if hist_key == "CMB" or "TOTAL" in hist_key: continue
    count[ hist_key ] = 0
    for hist_name in sorted( hists[ hist_key ].keys() ):
      parse = hist_parse( hist_name, samples )
      if doABCDNN and "ABCDNN" in hist_name: 
        if parse[ "IS SYST" ] and parse[ "SYST" ] in config.params[ "ABCDNN" ][ "SYSTEMATICS" ]:
          try: 
            hists[ "CMB" ][ hist_tag( "ABCDNN", parse[ "CATEGORY" ], parse[ "SYST" ] + parse[ "SHIFT" ] ) ].Add( hists[ hist_key ][ hist_name ] )
          except: 
            if args.verbose: print( "  + Creating {} histogram: {}".format( hist_key, hist_tag( "ABCDNN", parse[ "CATEGORY" ], parse[ "SYST" ] + parse[ "SHIFT" ] ) ) )
            hists[ "CMB" ][ hist_tag( "ABCDNN", parse[ "CATEGORY" ], parse[ "SYST" ] + parse[ "SHIFT" ] ) ] = hists[ hist_key ][ hist_name ].Clone( hist_tag( "ABCDNN", parse[ "CATEGORY" ], parse[ "SYST" ] + parse[ "SHIFT" ] ) )
        elif not parse[ "IS SYST" ]:
          try: 
            hists[ "CMB" ][ hist_tag( "ABCDNN", parse[ "CATEGORY" ] ) ].Add( hists[ hist_key ][ hist_name ] )
          except: 
            print( "  + Creating {} histogram: {}".format( hist_key, hist_tag( "ABCDNN", parse[ "CATEGORY" ] ) ) )
            hists[ "CMB" ][ hist_tag( "ABCDNN", parse[ "CATEGORY" ] ) ] = hists[ hist_key ][ hist_name ].Clone( hist_tag( "ABCDNN", parse[ "CATEGORY" ] ) )
        else:
          continue
      if parse[ "IS SYST" ] and "ABCDNN" not in hist_name:
        try: 
          hists[ "CMB" ][ hist_tag( parse[ "COMBINE" ], parse[ "CATEGORY" ], parse[ "SYST" ] + parse[ "SHIFT" ] ) ].Add( hists[ hist_key ][ hist_name ] )
        except:
          if args.verbose: print( "  + Creating {} histogram: {}".format( hist_key, hist_tag( parse[ "COMBINE" ], parse[ "CATEGORY" ], parse[ "SYST" ] + parse[ "SHIFT" ] ) ) )
          hists[ "CMB" ][ hist_tag( parse[ "COMBINE" ], parse[ "CATEGORY" ], parse[ "SYST" ] + parse[ "SHIFT" ] ) ] = hists[ hist_key ][ hist_name ].Clone( hist_tag( parse[ "COMBINE" ], parse[ "CATEGORY" ], parse[ "SYST" ] + parse[ "SHIFT" ] ) )
      elif "ABCDNN" not in hist_name:
        try: 
          hists[ "CMB" ][ hist_tag( parse[ "COMBINE" ], parse[ "CATEGORY" ] ) ].Add( hists[ hist_key ][ hist_name ] )
        except: 
          print( "  + Creating {} histogram: {}".format( hist_key, hist_tag( parse[ "COMBINE" ], parse[ "CATEGORY" ] ) ) )
          hists[ "CMB" ][ hist_tag( parse[ "COMBINE" ], parse[ "CATEGORY" ] ) ] = hists[ hist_key ][ hist_name ].Clone( hist_tag( parse[ "COMBINE" ], parse[ "CATEGORY" ] ) )
      count[ hist_key ] += 1


  for key in hists[ "CMB" ]: hists[ "CMB" ][ key ].SetDirectory(0)
  if config.params[ "HISTS" ][ "TTHFSF" ] != 1: hists = scale_ttbar( hists, config.params[ "HISTS" ][ "TTHFSF" ], config.params[ "HISTS" ][ "TTLFSF" ] )
  hists = set_zero( hists, categories, config.params[ "GENERAL" ][ "ZERO" ] )
  print( "[DONE] Consolidated histograms into Combine groupings in {:.2f} minutes:".format( ( time.time() - sTime ) / 60. ) )
  for key in count:
    print( "   + {}: {}".format( key, count[ key ] ) )
  return hists

def write_combine( hists, variable, categories, groups, templateDir, doABCDNN ):
  print( "[START] Writing Combine templates" )
  sTime = time.time()
  combine_name = "{}/template_combine_{}_UL{}.root".format( templateDir, variable, args.year )
  combine_file = TFile( combine_name, "RECREATE" )

  for category in categories:
    print( ">> Writing category: {}".format( category ) )
    hists[ "CMB" ][ hist_tag( "data_obs", category ) ].Write()
    print( "  + DAT > {}: {}".format( hist_tag( "data_obs", category ), hists[ "CMB" ][ hist_tag( "data_obs", category ) ].Integral() ) )

    for process in groups[ "SIG" ][ "PROCESS" ]:
      hists[ "CMB" ][ hist_tag( process, category ) ].Write()
      print( "  + SIG > {}: {}".format( hist_tag( process, category ), hists[ "CMB" ][ hist_tag( process, category ) ].Integral() ) )
      if config.options[ "GENERAL" ][ "SYSTEMATICS" ]:
        if args.verbose: print( "  + SIG SYST > {}".format( hist_tag( process, category ) ) )
        for syst in config.systematics[ "MC" ].keys():
          if args.year == "18" and syst.upper() == "PREFIRE": continue
          if not config.systematics[ "MC" ][ syst ][0] or "ABCD" in syst: continue
          if syst == "HD" and not config.options[ "GENERAL" ][ "HDAMP" ]: continue
          if syst == "UE" and not config.options[ "GENERAL" ][ "UE" ]: continue
          if args.verbose: print( "[INFO] Including {} to Combine template".format( syst ) )
          for shift in [ "UP", "DN" ]:
            if syst == "JEC":
              for systJEC in config.systematics[ "REDUCED JEC" ]:
                if not config.systematics[ "REDUCED JEC" ][ systJEC ]: continue
                systJEC_ = "JEC" + systJEC.replace( "Era", "20" + args.year ).replace( "APV", "" ).replace( "_", "" )
                hists[ "CMB" ][ hist_tag( process, category, systJEC_.upper() + shift ) ].Write()
            else:
              hists[ "CMB" ][ hist_tag( process, category, syst.upper() + shift ) ].Write()
      if config.options[ "GENERAL" ][ "PDF" ]:
        if args.verbose: print( "  + SIG PDF > {}".format( hist_tag( process, category ) ) )
        for i in range( config.params[ "GENERAL" ][ "PDF RANGE" ] ):
          hists[ "CMB" ][ hist_tag( process, category, "PDF" + str(i) ) ].Write()

    yield_total = sum( [ hists[ "CMB" ][ hist_tag( group, category ) ].Integral() for group in groups[ "BKG" ][ "SUPERGROUP" ] if hist_tag( group, category ) in hists[ "CMB" ].keys() ] )
    min_bkg_yield = 0 if args.region in [ "BASELINE", "ABCDNN" ] else config.params[ "HISTS" ][ "MIN BKG YIELD" ]
    max_bkg_error = 1.0 if args.region in [ "BASELINE", "ABCDNN" ] else config.params[ "HISTS" ][ "MAX BKG ERROR" ]
    for group in groups[ "BKG" ][ "SUPERGROUP" ]:
      scale_group = False
      error_group = Double(0)
      if hist_tag( group, category ) not in hists[ "CMB" ].keys(): continue
      yield_group = hists[ "CMB" ][ hist_tag( group, category ) ].IntegralAndError( 1, hists[ "CMB" ][ hist_tag( group, category ) ].GetXaxis().GetNbins(), error_group )
      if ( yield_group / yield_total <= min_bkg_yield or error_group / yield_group >= max_bkg_error ):
        scale_group = True
        if args.verbose and yield_group / yield_total <= min_bkg_yield: print( "  [WARN] {} beneath yield threshold, scaling by {:.1e} in Combine template".format( group, config.params[ "GENERAL" ][ "ZERO" ] ) )
        if args.verbose and error_group / yield_group >= max_bkg_error: print( "  [WARN] {} above error threshold, scaling by {:.1e} in Combine template".format( group, config.params[ "GENERAL" ][ "ZERO" ] ) )
        hists[ "CMB" ][ hist_tag( group, category ) ].Scale( config.params[ "GENERAL" ][ "ZERO" ] )
      hists[ "CMB" ][ hist_tag( group, category ) ].Write()
      print( "  + BKG > {}: {}".format( hist_tag( group, category ), hists[ "CMB" ][ hist_tag( group, category ) ].Integral() ) ) 
      if config.options[ "GENERAL" ][ "SYSTEMATICS" ]:
        for syst in config.systematics[ "MC" ].keys():
          if args.year == "18" and syst.upper() == "PREFIRE": continue
          if syst.upper() == "HD" and not config.options[ "GENERAL" ][ "HDAMP" ]: continue
          if syst.upper() == "UE" and not config.options[ "GENERAL" ][ "UE" ]: continue
          if "ABCD" in syst or not config.systematics[ "MC" ][ syst ][0]: continue
          for shift in [ "UP", "DN" ]:
            if syst == "JEC":
              for systJEC in config.systematics[ "REDUCED JEC" ]:
                if not config.systematics[ "REDUCED JEC" ][ systJEC ]: continue
                sysTag = "JEC" + systJEC.replace( "Era", "20" + args.year ).replace( "APV", "" ).replace( "_", "" ).upper() + shift
                if scale_group:
                  hists[ "CMB" ][ hist_tag( group, category, sysTag ) ].Scale( config.params[ "GENERAL" ][ "ZERO" ] )
                hists[ "CMB" ][ hist_tag( group, category, sysTag ) ].Write()
            else:
              sysTag = syst.upper() + shift
              if scale_group:
                hists[ "CMB" ][ hist_tag( group, category, sysTag ) ].Scale( config.params[ "GENERAL" ][ "ZERO" ] )
              hists[ "CMB" ][ hist_tag( group, category, sysTag ) ].Write()
        if args.verbose: print( "  + BKG SYST > {}".format( group ) )
            
      if config.options[ "GENERAL" ][ "PDF" ]:
        for i in range( config.params[ "GENERAL" ][ "PDF RANGE" ] ):
          if scale_group:
            hists[ "CMB" ][ hist_tag( group, category, "PDF" + str(i) ) ].Scale( config.params[ "GENERAL" ][ "ZERO" ] )
          hists[ "CMB" ][ hist_tag( group, category, "PDF" + str(i) ) ].Write()
        if args.verbose: print( "  + BKG PDF > {}".format( group ) )
  
    if doABCDNN and hist_parse( category, samples )[ "ABCDNN" ]:
      hists[ "CMB" ][ hist_tag( "ABCDNN", category ) ].Write()
      print( "  + BKG > {}: {}".format( hist_tag( "ABCDNN", category ), hists[ "CMB" ][ hist_tag( "ABCDNN", category ) ].Integral() ) )
      if config.options[ "GENERAL" ][ "SYSTEMATICS" ]:
        for syst in config.systematics[ "MC" ].keys():
          if syst.upper() not in config.params[ "ABCDNN" ][ "SYSTEMATICS" ]: continue
          if args.year == "18" and syst.upper() == "PREFIRE": continue
          for shift in [ "UP", "DN" ]:
            if config.systematics[ "MC" ][ syst ][0]:
              hists[ "CMB" ][ hist_tag( "ABCDNN", category, syst.upper() + shift ) ].Write()
        if config.options[ "GENERAL" ][ "PDF" ] and "PDF" in config.params[ "ABCDNN" ][ "SYSTEMATICS" ]:
          for i in range( config.params[ "GENERAL" ][ "PDF RANGE" ] ):
            hists[ "CMB" ][ hist_tag( "ABCDNN", category, "PDF" + str(i) ) ].Write()
        if args.verbose: print( "  + BKG SYST > ABCDNN" )

  combine_file.Close()
  print( "[DONE] Finished writing Combine templates in {:.2f} minutes".format( ( time.time() - sTime ) / 60. ) ) 
   
def make_tables( hists, categories, groups, variable, templateDir, lumiStr, doABCDNN ):
  def initialize():
    yield_table = { "YIELD": {}, "ERROR": {} }
    for category in categories:
      for stat in yield_table: yield_table[ stat ][ category ] = {}
      if config.options[ "GENERAL" ][ "SYSTEMATICS" ]:
        for syst in config.systematics[ "MC" ].keys():
          if not config.systematics[ "MC" ][ syst ][0]: continue
          if args.year == "18" and syst.upper() == "PREFIRE": continue
          if syst.upper() == "HD" and not config.options[ "GENERAL" ][ "HDAMP" ]: continue
          if syst.upper() == "UE" and not config.options[ "GENERAL" ][ "UE" ]: continue
          for shift in [ "UP", "DN" ]:
            if syst == "JEC":
              for systJEC in config.systematics[ "REDUCED JEC" ]:
                if not config.systematics[ "REDUCED JEC" ][ systJEC ]: continue
                systJEC_ = "JEC" + systJEC.replace( "Era", "20" + args.year ).replace( "APV", "" ).replace( "_", "" )
                for stat in yield_table: yield_table[ stat ][ hist_tag( category, systJEC_.upper() + shift ) ] = {}
            else:
              for stat in yield_table: yield_table[ stat ][ hist_tag( category, syst.upper() + shift ) ] = {}
    return yield_table
  
  def fill_yield( tables ):
    def get_yield( table, hist_key, group, category ):
      for process in group: 
        if args.verbose: ( "   + {}".format( process ) )
        table[ "YIELD" ][ category ][ process ] = hists[ hist_key ][ hist_tag( process, category ) ].Integral()
        if config.options[ "GENERAL" ][ "SYSTEMATICS" ] and "DAT" not in process.upper():
          for syst in config.systematics[ "MC" ].keys():
            if not config.systematics[ "MC" ][ syst ][0]: continue
            if args.year == "18" and syst.upper() == "PREFIRE": continue
            if process == "ABCDNN" and syst.upper() not in config.params[ "ABCDNN" ][ "SYSTEMATICS" ]: continue
            if "ABCD" in syst.upper() and group not in [ "BKG", "BKG SYST" ]: continue
            if syst.upper() == "HD" and not config.options[ "GENERAL" ][ "HDAMP" ]: continue
            if syst.upper() == "UE" and not config.options[ "GENERAL" ][ "UE" ]: continue
            for shift in [ "UP", "DN" ]:
              if syst == "JEC":
                for systJEC in config.systematics[ "REDUCED JEC" ]:
                  if not config.systematics[ "REDUCED JEC" ][ systJEC ]: continue
                  systJEC_ = "JEC" + systJEC.replace( "Era", "20" + args.year ).replace( "APV", "" ).replace( "_", "" )
                  table[ "YIELD" ][ hist_tag( category, systJEC_.upper() + shift ) ][ process ] = hists[ hist_key ][ hist_tag( process, category, systJEC_.upper() + shift ) ].Integral()
              else:
                table[ "YIELD" ][ hist_tag( category, syst.upper() + shift ) ][ process ] = hists[ hist_key ][ hist_tag( process, category, syst.upper() + shift ) ].Integral()
      return table

    for category in categories:
      if args.verbose: print( ">> Computing yields for DAT {}".format( category ) )
      tables = get_yield( tables, "CMB", [ "data_obs" ], category )
      if args.verbose: print( ">> Computing yields for SIGNAL {}".format( category ) )
      tables = get_yield( tables, "CMB", groups[ "SIG" ][ "PROCESS" ], category )
      #if args.verbose: print( ">> Computing yields for BACKGROUND physics groups" )
      #tables = get_yield( tables, "BKG", groups[ "BKG" ][ "PROCESS" ].keys(), category )
      if doABCDNN and hist_parse( category, samples )[ "ABCDNN" ]:
        if args.verbose: print( ">> Computing yields for ABCDNN BACKGROUND {}".format( category ) )
        tables = get_yield( tables, "CMB", [ "ABCDNN" ] + config.params[ "ABCDNN" ][ "MINOR BKG" ], category )
      else:
        if args.verbose: print( ">> Computing yields for BACKGROUND Combine {}".format( category ) )
        tables = get_yield( tables, "CMB", groups[ "BKG" ][ "SUPERGROUP" ].keys(), category )
      
      tables[ "YIELD" ][ category ][ "TOTAL BKG" ] = 0
      for group in groups[ "BKG" ][ "SUPERGROUP" ]:
        if hist_tag( group, category ) not in hists[ "CMB" ].keys() or ( doABCDNN and hist_parse( category, samples )[ "ABCDNN" ] ): continue
        tables[ "YIELD" ][ category ][ "TOTAL BKG" ] += hists[ "CMB" ][ hist_tag( group, category ) ].Integral()
      if doABCDNN and hist_parse( category, samples )[ "ABCDNN" ]:
        tables[ "YIELD" ][ category ][ "TOTAL BKG" ] += hists[ "CMB" ][ hist_tag( "ABCDNN", category ) ].Integral()
        for group in config.params[ "ABCDNN" ][ "MINOR BKG" ]:
          tables[ "YIELD" ][ category ][ "TOTAL BKG" ] += hists[ "CMB" ][ hist_tag( group, category ) ].Integral()

      tables[ "YIELD" ][ category ][ "DATA:BKG" ] = tables[ "YIELD" ][ category ][ "data_obs" ] / ( tables[ "YIELD" ][ category ][ "TOTAL BKG" ] + config.params[ "GENERAL" ][ "ZERO" ] )
   
    return tables
  
  def fill_error( tables ):
    def get_error( table, hist_key, group, category ):
      for process in group:
        if args.verbose: print( "   + {}".format( process ) )
        table[ "ERROR" ][ category ][ process ] = 0
        for i in range( 1, hists[ hist_key ][ hist_tag( process, category ) ].GetXaxis().GetNbins() + 1 ):
          table[ "ERROR" ][ category ][ process ] += hists[ hist_key ][ hist_tag( process, category ) ].GetBinError(i)**2
        table[ "ERROR" ][ category ][ process ] = math.sqrt( table[ "ERROR" ][ category ][ process ] )
        if config.options[ "GENERAL" ][ "SYSTEMATICS" ] and "DAT" not in process.upper():
          for syst in config.systematics[ "MC" ].keys():
            if not config.systematics[ "MC" ][ syst ][0]: continue
            if "ABCD" in process and syst.upper() not in config.params[ "ABCDNN" ][ "SYSTEMATICS" ]: continue
            if "ABCD" in syst.upper() and group not in [ "BKG", "BKG SYST" ]: continue
            if syst.upper() == "HD" and not config.options[ "GENERAL" ][ "HDAMP" ]: continue
            if syst.upper() == "UE" and not config.options[ "GENERAL" ][ "UE" ]: continue
            if syst.upper() == "PREFIRE" and args.year not in [ "16APV", "16", "17" ]: continue
            for shift in [ "UP", "DN" ]:
              if syst == "JEC":
                for systJEC in config.systematics[ "REDUCED JEC" ]:
                  if not config.systematics[ "REDUCED JEC" ][ systJEC ]: continue
                  systJEC_ = "JEC" + systJEC.replace( "Era", "20" + args.year ).replace( "APV", "" ).replace( "_", "" )
                  table[ "ERROR" ][ hist_tag( category, systJEC_.upper() + shift ) ][ process ] = 0
                  for i in range( 1, hists[ hist_key ][ hist_tag( process, category, systJEC_.upper() + shift ) ].GetXaxis().GetNbins() + 1 ):
                    table[ "ERROR" ][ hist_tag( category, systJEC_.upper() + shift ) ][ process ] += hists[ hist_key ][ hist_tag( process, category, systJEC_.upper() + shift ) ].GetBinError(i)**2
                  table[ "ERROR" ][ hist_tag( category, systJEC_.upper() + shift ) ][ process ] = hists[ hist_key ][ hist_tag( process, category, systJEC_.upper() + shift ) ].Integral()
              else:
                table[ "ERROR" ][ hist_tag( category, syst.upper() + shift ) ][ process ] = 0
                for i in range( 1, hists[ hist_key ][ hist_tag( process, category, syst.upper() + shift ) ].GetXaxis().GetNbins() + 1 ):
                  table[ "ERROR" ][ hist_tag( category, syst.upper() + shift ) ][ process ] += hists[ hist_key ][ hist_tag( process, category, syst.upper() + shift ) ].GetBinError(i)**2
                table[ "ERROR" ][ hist_tag( category, syst.upper() + shift ) ][ process ] = hists[ hist_key ][ hist_tag( process, category, syst.upper() + shift ) ].Integral()
      return table

    for category in categories:
      if args.verbose: print( ">> Computing errors for DATA histograms" )
      tables = get_error( tables, "CMB", [ "data_obs" ], category )
      if args.verbose: print( ">> Computing errors for SIGNAL histograms" )
      tables = get_error( tables, "CMB", groups[ "SIG" ][ "PROCESS" ], category )
      #if args.verbose: print( ">> Computing errors for BACKGROUND physics groups" )
      #tables = get_error( tables, "BKG", groups[ "BKG" ][ "PROCESS" ].keys(), category )
      if doABCDNN and hist_parse( category, samples )[ "ABCDNN" ]:
        if args.verbose: print( ">> Computing errors for ABCDNN BACKGROUND Combine groups" )
        tables = get_error( tables, "CMB", [ "ABCDNN" ] + config.params[ "ABCDNN" ][ "MINOR BKG" ], category )
      else:
        if args.verbose: print( ">> Computing errors for BACKGROUND Combine groups" )
        tables = get_error( tables, "CMB", groups[ "BKG" ][ "SUPERGROUP" ].keys(), category )
    

      tables[ "ERROR" ][ category ][ "TOTAL BKG" ] = 0
      tables[ "ERROR" ][ category ][ "DATA:BKG" ] = 0
      for group in groups[ "BKG" ][ "SUPERGROUP" ]:
        if hist_tag( group, category ) not in hists[ "CMB" ].keys(): continue
        if doABCDNN and hist_parse( category, samples )[ "ABCDNN" ] and group not in config.params[ "ABCDNN" ][ "MINOR BKG" ]: continue
        for i in range( 1, hists[ "CMB" ][ hist_tag( group, category ) ].GetXaxis().GetNbins() + 1 ):
          tables[ "ERROR" ][ category ][ "TOTAL BKG" ] += hists[ "CMB" ][ hist_tag( group, category ) ].GetBinError(i)**2
      if doABCDNN and hist_parse( category, samples )[ "ABCDNN" ]:
        for i in range( 1, hists[ "CMB" ][ hist_tag( "ABCDNN", category ) ].GetXaxis().GetNbins() + 1 ):
          tables[ "ERROR" ][ category ][ "TOTAL BKG" ] += hists[ "CMB" ][ hist_tag( "ABCDNN", category ) ].GetBinError(i)**2

      tables[ "ERROR" ][ category ][ "TOTAL BKG" ] = math.sqrt( tables[ "ERROR" ][ category ][ "TOTAL BKG" ] ) 
      yield_d = tables[ "YIELD" ][ category ][ "data_obs" ] 
      yield_b = tables[ "YIELD" ][ category ][ "TOTAL BKG" ]
      error_d = math.sqrt( yield_d )
      error_b = math.sqrt( yield_b )
      tables[ "ERROR" ][ category ][ "DATA:BKG" ] = math.sqrt( ( error_d / yield_b )**2 + ( yield_d * error_d / yield_b**2 )**2 )      
      tables[ "ERROR" ][ category ][ "DATA:BKG" ] = math.sqrt( tables[ "ERROR" ][ category ][ "DATA:BKG" ] )

    return tables

  tables = initialize()
  tables = fill_yield( tables )
  tables = fill_error( tables )

  return tables

def print_tables( tables, categories, groups, variable, templateDir ):
  def yield_table():
    def format_section( table, title, columns, systematic, precision ):
      pm = "{:." + precision + "f} pm {:." + precision + "f}"
      table.append( [ "YIELD:", title ] )
      table.append( [ "CATEGORY" ] + columns )
      sum_row = [ "TOTAL" ]
      process_stat = {
        key: {
          process: 0 for process in columns
        } for key in [ "YIELD", "ERROR" ]
      }
      for category in sorted( categories ):
        category_mc = category if systematic == "NOMINAL" else hist_tag( category, systematic )
        row = [ category_mc ]
        for process in columns:
          if process == "data_obs": category_mc = category
          if process not in tables[ "YIELD" ][ category_mc ].keys():
            row.append( pm.format(
              0,0
            ) )
          else:
            row.append( pm.format(
              tables[ "YIELD" ][ category_mc ][ process ],
              tables[ "ERROR" ][ category_mc ][ process ]
            ) )
            process_stat[ "YIELD" ][ process ] += tables[ "YIELD" ][ category_mc ][ process ]
            process_stat[ "ERROR" ][ process ] += tables[ "ERROR" ][ category_mc ][ process ]**2
        table.append( row )
      for process in columns:
        process_stat[ "ERROR" ][ process ] = -1 if ":" in process else math.sqrt( process_stat[ "ERROR" ][ process ] )
        process_stat[ "YIELD" ][ process ] = -1 if ":" in process else process_stat[ "YIELD" ][ process ]
        sum_row.append( pm.format( 
          process_stat[ "YIELD" ][ process ],
          process_stat[ "ERROR" ][ process ]
        ) )
      table.append( sum_row )
      table.append( [ "" ] )
      table.append( [ "" ] )
      return table


    print( ">> Printing out the nominal table" )
    table = []

    #table = format_section( table, "PHYSICS PROCESS", [ "DAT" ] + list( groups[ "BKG" ][ "PROCESS" ].keys() ), "2" )
    table = format_section( table, "COMBINE ANALYSIS", [ "data_obs" ] + list( groups[ "BKG" ][ "SUPERGROUP" ].keys() ), "NOMINAL", "2" )
    table = format_section( table, "COMBINE ABCDNN", [ "data_obs", "ABCDNN" ] + config.params[ "ABCDNN" ][ "GROUPS" ], "NOMINAL", "2" )
    table = format_section( table, "SIGNAL", groups[ "SIG" ][ "PROCESS" ], "NOMINAL", "3" )
    table = format_section( table, "SUMMARY", [ "TOTAL BKG", "DATA:BKG" ], "NOMINAL", "4" )

    return table

  def an_table():
    pass
  def pas_table( tables ):
    pass
  def systematic_table( tables ):
    pass
  def print_table( table, section_header, file_name ):
    def get_max_width( table, index ):
      max_width = 0
      for row in table:
        try:
          n = len( format( row[ index ] ) )
          if n > max_width: max_width = n
        except:
          pass
      return max_width

    def print_section( section ):
      column_padding = []
      max_columns = 0
      for row in section: 
        if len( row ) > max_columns: max_columns = len( row )
      for i in range( max_columns ): column_padding.append( get_max_width( section, i ) )
      for row in section:
        print >> file_name, format( row[0] ).ljust( column_padding[0] + 1 ),
        for i in range( 1, len( row ) ):
          column = format( row[i] ).ljust( column_padding[i] + 2 )
          print >> file_name, column,
        print >> file_name

    def get_section( table ):
      sections = []
      section = []
      for i, row in enumerate( table ):
        if section_header in row or i == len( table ) - 1: 
          if i != 0: sections.append( section )
          section = [ row ]
        else: section.append( row )
      return sections

    print( ">> Writing out {} table to: {}".format( section_header, file_name ) )
    for section in get_section( table ): print_section( section )    


  if not os.path.exists( os.path.join( templateDir, "tables/" ) ): os.system( "mkdir -vp {}".format( os.path.join( templateDir, "tables/" ) ) )
  print_table( yield_table(), "YIELD:", open( os.path.join( templateDir, "tables/", "yield_table.txt" ), "w" ) )


def main():
  start_time = time.time()

  template_prefix = config.region_prefix[ args.region ]
  templateDir = os.path.join( os.getcwd(), "{}_UL{}_{}".format( template_prefix, args.year, args.tag ) )
  categories = get_categories( templateDir )
  groups = samples.groups 

  for variable in args.variables:
    hists = load_histograms( variable, templateDir, categories )
    for hist_key in hists:
      if len( hists[ hist_key ].keys() ) <= 0: continue
      hists = clean_histograms( hists, hist_key, config.params[ "HISTS" ][ "LUMISCALE" ], config.params[ "HISTS" ][ "REBIN" ] )
    hists = combine_histograms( hists, variable, categories, groups, config.options[ "GENERAL" ][ "ABCDNN" ] )
    write_combine( hists, variable, categories, groups, templateDir, config.options[ "GENERAL" ][ "ABCDNN" ] )
    tables = make_tables( hists, categories, groups, variable, templateDir, config.lumiStr[ args.year ], config.options[ "GENERAL" ][ "ABCDNN" ] )
    print_tables( tables, categories, groups, variable, templateDir )
    del hists
    #del tables

main()
