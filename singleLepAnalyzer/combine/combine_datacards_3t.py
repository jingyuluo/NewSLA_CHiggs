import CombineHarvester.CombineTools.ch as ch
from argparse import ArgumentParser
import os
import datetime
import yaml
import ROOT

parser = ArgumentParser()
parser.add_argument( "-c", "--config", help = "YAML configuration file" )
parser.add_argument( "--validation", action = "store_true", help = "Combine custom datacard combinations" )
parser.add_argument( "--correlations", action = "store_true", help = "Adjust bin naming to account for correlations" )
parser.add_argument( "--single", action = "store_true", help = "Create datacards for single systematics" )
parser.add_argument( "--positive", action = "store_true", help = "Turn off zero and/or negative bins" )
parser.add_argument( "--verbose", action = "store_true", help = "Verbosity of messages" )
args = parser.parse_args()

def hist_key( *args ):
  key = args[0]
  for arg in args[1:]: key += "_{}".format( arg )
  return key
  
class file(object):
  def __init__( self, tag, workingDir ):
    self.tag = tag
    self.workingDir = workingDir
    self.count = 0
    pass
  def write_shell( self, cards ):
    self.file = open( "combined_{}/start_T2W.sh".format( self.tag ), "w" )
    self.file.write( "#!/bin/sh\n" )
    self.file.write( "set -e\n" )
    self.file.write( "ulimit -s unlimited\n" )
    self.file.write( "\n" )
    self.file.write( "cd {}".format( self.workingDir ) )
    self.file.write( "source /cvmfs/cms.cern.ch/cmsset_default.sh\n" )
    self.file.write( "eval `scramv1 runtime -sh`\n" )
    self.file.write( "cd {}/combined_{}/\n\n".format( self.workingDir, self.tag ) )
    self.file.write( "DATACARD=\n" )
    self.file.write( "case $1 in\n" )
    for card in cards:
      self.add_datacard( card )
    self.file.write( "esac\n" )
    self.file.write( "text2workspace.py $DATACARD -o ${DATACARD:0:-4}.root -m 125 --channel-masks\n" )
    self.file.close()
    pass
  def write_condor( self ):
    self.file = open( "combined_{}/condor_T2W.sub".format( self.tag ), "w" )
    self.file.write( "executable = start_T2W.sh\n" )
    self.file.write( "arguments = ${PROCID}\n" )
    self.file.write( "output = start_T2W.out\n" )
    self.file.write( "error = start_T2W.err\n" )
    self.file.write( "log = start_T2W.log\n" )
    self.file.write( "queue\n" )
    self.file.close()
    pass
  def write_to( self, text ):
      self.file.write( text )
  def add_datacard( self, datacard ):
    self.file.write( "  {})\n".format( self.count ) )
    self.file.write( "    DATACARD={}\n".format( datacard ) )
    self.file.write( "    ;;\n" )
    self.count += 1
  def close( self ):
    self.file.close()
  
class harvester(object):
  def __init__( self, datacards, tag, postfix, combinations, correlations, removeUnct, verbose, analysis = "TTTX", era = "Run2UL" ):
    print( "[COMBINE HARVESTER] Instantiating harvester class for {} {} analysis".format( era, analysis ) )
    print( "  [INFO] Verbosity {}".format( "ON" if verbose else "OFF" ) )
    self.verbose = verbose            # verbosity of messages
    self.analysis = analysis          # analysis tag, arbitrary
    self.era = era                    # era tag, arbitrary
    self.datacards = datacards        # dictionary of datacard paths associated to each analysis channel
    self.tag = tag                    # tag used for combination campaign
    self.postfix = postfix            # dictionary of postfixes associated to each analysis channel
    self.combinations = combinations  # dictionary of desired combinations of datacards
    self.correlations = correlations  # dictionary of correlated uncertainties across different datacards
    self.removeUnct = removeUnct      # dictionary of uncertainties to remove from each datacard
    self.cards_to_add = []
    self.ch = ch.CombineHarvester()
    pass
  
  def load_datacards( self ):
    print( "[COMBINE HARVESTER] Parsing through datacards, and extracting and organizing bin information" )
    if self.verbose: print( "  [INFO] Considering {} datacards".format( len( self.datacards ) ) )
    self.bins_sorted = { "ALL": [] }
    for datacard in self.datacards:
      self.bins_sorted[ datacard ] = []
      self.ch.ParseDatacard( self.datacards[ datacard ], self.analysis, self.era, datacard, 0, "125" )
      self.ch.cp().channel( [ datacard ] ).ForEachObj( lambda x: self.bins_sorted[ "ALL" ].append( ( x.bin(), hist_key( datacard, x.bin() ) ) ) )
      self.ch.cp().channel( [ datacard ] ).ForEachObj( lambda x: self.bins_dict[ datacard ].append( hist_key( datacard, x.bin() ) ) )
    pass
  
  def remove_uncertainties( self ):
    for key in self.removeUnct:
      for source in self.removeUnct[ key ]:
        self.ch.FilterSysts( lambda x: x.channel() == channel and x.name() == source )
    pass
  
  def format_bin_names( self ):
    print( "[COMBINE HARVESTER] Renaming bins to distinguish channels in final datacard..." )
    for datacard in self.datacards:
      self.ch.cp().channel( [ datacard ] ).ForEachObj( lambda x: x.set_bin( hist_key( datacard, x.bin() ) ) )
    print( "[COMBINE HARVESTER] Renaming autoMCstats bins..." )
    for old, new in list( set( self.bins_sorted[ "ALL" ] ) ):
      self.ch.RenameAutoMCStatsBin( old, new )
    pass
  
  def rename_uncertainties( self ):
    def rename_uncertainty( postfix ):
      def rename_function( x ):
        if x.type() != "rateParam": x.set_name( hist_key( x.name(), suffix ) )
      return rename_function
     
    for analysis in self.postfix:
      postfix = self.config[ "POSTFIX" ][ analysis ]
      self.ch.cp().channel( [ analysis ] ).ForEachSyst( rename_uncertainty( postfix ) )
    pass
  
  def modify_uncertainties( self, bPositive_ ):
    if bPositive_:
      print( "[COMBINE HARVESTER] Turning off zero and/or negative bins" )
      self.ch.FilterProces( lambda x: x.rate() <= 0 )
      self.ch.FilterSysts( lambda x: x.type() == "shape" and x.value_u() <= 0 )
      self.ch.FilterSysts( lambda x: x.type() == "shape" and x.value_d() <= 0 )
    pass
  
  def write_datacard( self, ch_, tag_ ):
    print( "[COMBINE HARVESTER] Preparing main datacard" )
    rOut = ROOT.TFile( "combined_{}/{}_combined_{}.root".format( self.tag, self.analysis, tag_ ), "RECREATE" )
    ch_.WriteDatacard( "combined_{}/{}_combined_{}.txt".format( self.tag, self.analysis, tag_ ), rOut )
    rOut.Close()
    self.cards_to_add.append( "{}_combined_{}.txt".format( self.analysis, tag_ ) )
  
  def combine_correlations( self, bIndividual_ ):
    print( "[COMBINE HARVESTER] Writing datacards for combined correlations" )
    os.system( "mkdir combined_{}/correlated_cards/".format( self.tag ) )
    self.ch_all = self.ch.deep()
    for category in self.correlations.keys():
      print( "  + {}".format( category ) ) 
      if bIndividual_:
        ch_syst = self.ch.deep()
      for source in self.correlations[ category ].keys():
        for content in self.correlations[ category ][ source ]:
          self.ch_all.cp().channel( [ content[ "TAG" ] ] ).syst_name( [
            hist_key( content[ "NAME" ], self.postfix[ content[ "TAG" ] ] )
          ] ).ForEachSyst(
            lambda x: x.set_name( source ) 
          )
          if bIndividual_:
            ch_syst.cp().channel( [ content[ "TAG" ] ] ).syst_name( [
              hist_key( content[ "NAME" ], self.postfix[ content[ "TAG" ] ] )
            ] ).ForEachSyst(
              lambda x: x.set_name( source )
            )
      if bIndividual_:
        self.write_datacard( ch_syst, category )
        del ch_syst
    self.write_datacard( self.ch_all, "ALLSYST" )
    pass
 
  def combine_analysis( self ):
    print( "[COMBINE HARVESTER] Writing datacards for analysis combinations" )
    os.system( "mkdir combined_{}/correlated_analysis/".format( self.tag ) )
    for combination in self.combinations:
      ch_combination = self.ch_all.deep()
      for datacard in self.datacards:
        if datacard not in combination[ "CONTENT" ]:
          ch_combination.channel( [ datacard ], False )
      self.write_datacard( ch_combination, combination[ "NAME" ] )
      
  def autoMCstats( self, bIndividual_ ):
    print( "[COMBINE HARVESTER] Preparing autoMCStats configs from main datacard to individual cards" )
    os.system( "grep autoMCStats combined_{}/{}_combined_{}.txt > combined_{}/autoMCStats_list.txt".format( self.tag, self.analysis, self.tag ) )
    if bIndividual_:
      print( "  + Propagating autoMCStats to single correlated systematic datacards" )
      os.system( "for file in combined_{0}/correlated_cards/*.txt; do cat combined_{0}/autoMCStats_list.txt >> $file; done".format( self.tag ) )
    for combination in self.combinations:
      for datacard in self.datacards:
        if datacard in combination[ "CONTENT" ]:
          os.system( "grep '^{2}' combined_{0}/autoMCStats_lists.txt >> combined_{0}/correlated_cards/{1}_combined_{3}.txt".format(
            self.tag, self.analysis, datacard, combination[ "NAME" ]
          ) )
  
def main():
  tStart = datetime.time.now()
  print( "[START] Running combine_datacards_3t.py with the following settings:" )
  print( "  + Config file: {}".format( args.config ) )
  print( "  + Negative/Zero Bin Correction: {}".format( args.positive ) )
  print( "  + Verbose: {}".format( args.verbose ) )
  os.system( "mkdir -vp combined_{}".format( tag ) )
  with open( args.config ) as f:
    config = yaml.safe_load( f )

  out = file( config[ "TAG" ] )
  condor = file( config[ "TAG" ] )
 
  ch_ttt = harvester()
  
  ch_ttt.load_datacards(
    datacards = config[ "DATACARDS" ]
  )
  ch_ttt.remove_uncertainties()
  ch_ttt.modify_uncertainties( args.positive )
  ch_ttt.format_bin_names()
  ch_ttt.rename_uncertainties()
  ch_ttt.combine_correlations( args.single )
  ch_ttt.combine_analysis()
  ch_ttt.autoMCstats( args.single )
  
  out.write_shell( ch_ttt.cards_to_add )
  condor.write_condor()
  
  tEnd = datetime.time.now()
  print( "[DONE] Finished running in {:.2f} minutes".format( ( tEnd - tStart ) / 60. ) )
