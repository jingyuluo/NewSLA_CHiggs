import os,sys, time
import config
from argparse import ArgumentParser

cmsswbase = os.path.join( os.getcwd(), ".." )

parser = ArgumentParser()
parser.add_argument( "-s", "--step", required = True, help = "Options: 1-5,impact,fit,correlation,breakdown" )
parser.add_argument( "-y", "--years", nargs = "+",  help = "Options: 16APV, 16, 17, 18, Run2" )
parser.add_argument( "-t", "--tags", nargs = "+", required = True )
parser.add_argument( "-v", "--variables", nargs = "+", required = True )
parser.add_argument( "-r", "--region", default = "SR" )
parser.add_argument( "--html", help = "Path to html page on BRUX" )
parser.add_argument( "--verbose", action = "store_true" )
args = parser.parse_args()

if args.region not in list( config.region_prefix.keys() ): quit( "[ERR] Invalid option used for -r (--region). Quitting." )

verbose = "--verbose" if args.verbose else ""

def get_trainings( tags, years, variables ):
  trainings = []
  for tag in tags:
    for year in years:
      trainings.append( {
        "year": year,
        "variable": variables,
        "tag": tag,
        "path": config.inputDir[ year ]
      } )
  return trainings
    
def condor_template( logName, nameCondor ):
  jobName = os.path.join( logName, nameCondor + ".job" )
  jobFile = open( jobName, "w" )
  jobFile.write(
"""universe = vanilla \n\
Executable = {0}.sh \n\
Should_Transfer_Files = YES \n\
WhenToTransferOutput = ON_EXIT \n\
JobBatchName = {1} \n\
request_memory = 3072 \n\
Output = {0}.out \n\
Error = {0}.err \n\
Log = {0}.log \n\
Notification = Error \n\
Arguments = \n\
Queue 1""".format(
  os.path.join( os.getcwd(), logName, nameCondor ),
  nameCondor
)
  )
  jobFile.close()
  os.system( "condor_submit {}".format( jobName ) )
  
def make_templates():
  trainings = get_trainings( args.tags, args.years, args.variables )
  os.chdir( "makeTemplates" )
  for training in trainings:
    for variable in training[ "variable" ]:
      command = "python condor_templates.py -y {} -v {} -p {} -r {}".format(
        training[ "year" ],
        variable,
        training[ "tag" ],
        args.region
      )
      os.system( command ) 
      time.sleep( 1 )
  os.chdir( ".." )
                
def format_templates():
  trainings = get_trainings( args.tags, args.years, args.variables )
  doABCDnn = "--abcdnn" if config.options[ "GENERAL" ][ "ABCDNN" ] else ""
  for training in trainings:
    for variable in training[ "variable" ]:
      os.chdir( "makeTemplates" )
      nameLog = "log_UL{}_{}_{}_{}".format( training[ "year" ], variable, args.region, training[ "tag" ] )
      nameCondor = "SLA_step2_{}_{}_{}_{}".format( training[ "year" ], variable, args.region, training[ "tag" ] )
      if not os.path.exists( nameLog ): os.system( "mkdir -vp {}".format( nameLog ) )
      shell = open( "{}/{}.sh".format( nameLog, nameCondor ), "w" )
      shell.write(
"#!/bin/bash \n\
source /cvmfs/cms.cern.ch/cmsset_default.sh \n\
cd {0} \n\
eval `scramv1 runtime -sh` \n\
cd {1} \n\
python templates.py -y {2} -t {3} -v {4} -r {5} \n\
python modify_binning.py -y {2} -t {3} -v {4} -r {5} {6}".format( 
  cmsswbase, os.getcwd(), training[ "year" ], training[ "tag" ], variable, args.region, doABCDnn )
      )
      shell.close()
      condor_template( nameLog, nameCondor )
      os.chdir( "../" )
  
def plot_templates():
  trainings = get_trainings( args.tags, args.years, args.variables )
  try:
    argHTML = "--html {}".format( args.html )
  except:
    argHTML = ""
  for training in trainings:
    for variable in training[ "variable" ]:
      os.chdir( "makeTemplates" )
      nameLog = "log_UL{}_{}_{}_{}".format( training[ "year" ], variable, args.region, training[ "tag" ] )
      nameCondor = "SLA_step3_{}_{}_{}_{}".format( training[ "year" ], variable, args.region, training[ "tag" ] )
      shell = open( "{}/{}.sh".format( nameLog, nameCondor ), "w" )
      shell.write(
"#!/bin/bash \n\
source /cvmfs/cms.cern.ch/cmsset_default.sh \n\
cd {0} \n\
eval `scramv1 runtime -sh` \n\
cd {1} \n\
python plot_templates.py -y {2} -v {3} -t {4} -r {5} --ratios --shifts --systematics {6} \n\
python plot_templates.py -y {2} -v {3} -t {4} -r {5} --templates {6}".format( 
  cmsswbase, os.getcwd(), 
  training[ "year" ], variable, training[ "tag" ], args.region, argHTML
)
      )
      shell.close()
      condor_template( nameLog, nameCondor )
      os.chdir( ".." )
 
def run_combine():
  systCombo = {
    "": "--normSyst --shapeSyst --theorySyst",
    "noShapenoNorm": "--theorySyst",
    "noShapenoTheory": "--normSyst",
    "noNormnoTheory": "--shapeSyst",
    "noTheory": "--normSyst --shapeSyst",
    "noNorm": "--shapeSyst --theorySyst",
    "noShape": "--normSyst --theorySyst",
    "noShapenoNormnoTheory": ""
  }
  trainings = get_trainings( args.tags, args.years, args.variables )
  tagSmooth = "" if not config.options[ "COMBINE" ][ "SMOOTH" ] else config.params[ "MODIFY BINNING" ][ "SMOOTHING ALGO" ].upper()
  tagABCDnn = "" if not config.options[ "COMBINE" ][ "ABCDNN" ] else "ABCDNN"
  for training in trainings:
    for variable in training[ "variable" ]:
      for systTag in systCombo:
        postfix = tagABCDnn + systTag + tagSmooth
        if ( systTag != "" and systTag != "noShapenoNormnoTheory" ) and not config.options[ "COMBINE" ][ "GROUPS" ]: continue
        os.chdir( "combine" )
        nameLog = "log_UL{}_{}_{}_{}{}".format( training[ "year" ], variable, args.region, training[ "tag" ], systTag )
        nameCondor = "SLA_step4_{}_{}_{}_{}".format( training[ "year" ], variable, args.region, training[ "tag" ] )
        if not os.path.exists( nameLog ): os.system( "mkdir {}".format( nameLog ) )
        shell = open( os.path.join( nameLog, nameCondor + ".sh" ), "w" )
        shell.write(
"#!/bin/bash\n\
source /cvmfs/cms.cern.ch/cmsset_default.sh\n\
cd {0} \n\
eval `scramv1 runtime -sh`\n\
cd {1} \n\
python create_datacard.py -y {2} -v {3} -r {4} -t {5} {6} \n\
cd limits_UL{2}_{3}_{4}_{5}_{7}\n\
combineTool.py -M T2W -i cmb/ -o workspace.root --parallel 8\n\
ValidateDatacards.py cmb/combined.txt.cmb --printLevel 2\n\
mv validation.json validation_UL{2}_{3}_{4}_{5}_{7}.json\n\
combine -M FitDiagnostics cmb/combined.txt.cmb {8}\n\
mkdir -vp cmb/FitDiagnostics_UL{2}_{3}_{4}_{5}_{7}\n\
mv *.png cmb/FitDiagnostics_UL{2}_{3}_{4}_{5}_{7}/ \n\
combine -M Significance cmb/workspace.root {9} > significance_merge{11}_stat{12}.txt\n\
combine -M AsymptoticLimits cmb/workspace.root {10} > limits_merge{11}_stat{12}.txt\n\
cd ..\n\
python systematicsAnalyzer.py limits_UL{2}_{3}_{4}_{5}_{7}/cmb/combined.txt.cmb -a > limits_UL{2}_{3}_{4}_{5}_{7}/cmb/datacard_UL{2}_{3}_{4}_{5}_{7}.html".format(
            cmsswbase, os.getcwd(), training[ "year" ], variable, args.region, training[ "tag" ], 
            systCombo[ systTag ], postfix, 
            " ".join( config.params[ "COMBINE" ][ "FITS" ][ "ARGS" ] ),
            " ".join( config.params[ "COMBINE" ][ "SIGNIFICANCE" ][ "ARGS" ] ),
            " ".join( config.params[ "COMBINE" ][ "LIMITS" ][ "ARGS" ] ),
            config.params[ "MODIFY BINNING" ][ "MIN MERGE" ],
            str( config.params[ "MODIFY BINNING" ][ "STAT THRESHOLD" ] ).replace( ".", "p" )
          )
        )
        shell.close()
        condor_template( nameLog, nameCondor )
        os.chdir( ".." )
  
def combine_years():
  tagSmooth = "" if not config.options[ "COMBINE" ][ "SMOOTH" ] else config.params[ "MODIFY BINNING" ][ "SMOOTHING ALGO" ].upper()
  tagABCDnn = "" if not config.options[ "COMBINE" ][ "ABCDNN" ] else "ABCDNN" 
  for variable in args.variables:
    for tag in args.tags:
      os.chdir( "combine" )
      nameLog = "log_{}_{}_{}".format( variable, args.region, tag )
      nameCondor = "SLA_step5_{}_{}_{}_{}".format( variable, args.region, tag, tagABCDnn + tagSmooth )
      if not os.path.exists( nameLog ): os.system( "mkdir -vp {}".format( nameLog ) )
      tagAllSyst = "{}_{}_{}_{}{}".format( variable, args.region, tag, tagABCDnn, tagSmooth )
      tagNoSyst = "{}_{}_{}_{}noShapenoNormnoTheory{}".format( variable, args.region, tag, tagABCDnn, tagSmooth )
      shell = open( "{}/{}.sh".format( nameLog, nameCondor ), "w" )
      shell.write(
"#!/bin/bash\n\
source /cvmfs/cms.cern.ch/cmsset_default.sh\n\
cd {0}\n\
eval `scramv1 runtime -sh`\n\
cd {1}\n\
mkdir -vp Results/{2}/\n\
combineCards.py UL16APV=limits_UL16APV_{2}/cmb/combined.txt.cmb UL16=limits_UL16_{2}/cmb/combined.txt.cmb UL17=limits_UL17_{2}/cmb/combined.txt.cmb  UL18=limits_UL18_{2}/cmb/combined.txt.cmb > Results/{2}/workspace.txt\n\
text2workspace.py Results/{2}/workspace.txt -o Results/{2}/workspace.root --channel-masks\n\
combine -M Significance Results/{2}/workspace.root {5} > Results/{2}/significance_merge{7}_stat{8}.txt\n\
combine -M AsymptoticLimits Results/{2}/workspace.root {6} > Results/{2}/limits_merge{7}_stat{8}.txt\n\
cd Results/{2}/\n\
ValidateDatacards.py workspace.root --printLevel 2\n\
combine -M FitDiagnostics workspace.root {4}\n\
mkdir -vp FitDiagnostics_{2}\n\
mv *.png FitDiagnostics_{2}/ \n\
cd ../../\n\
python systematicsAnalyzer.py Results/{2}/workspace.txt -a > Results/{2}/datacard_{2}.html\n\
mkdir -vp Results/{3}/\n\
combineCards.py UL16APV=limits_UL16APV_{3}/cmb/combined.txt.cmb UL16=limits_UL16_{3}/cmb/combined.txt.cmb UL17=limits_UL17_{3}/cmb/combined.txt.cmb  UL18=limits_UL18_{3}/cmb/combined.txt.cmb > Results/{3}/workspace.txt \n\
text2workspace.py Results/{3}/workspace.txt -o Results/{3}/workspace.root --channel-masks\n\
combine -M Significance Results/{3}/workspace.root {5} > Results/{3}/significance_merge{7}_stat{8}.txt \n\
combine -M AsymptoticLimits Results/{3}/workspace.root {6} > Results/{3}/limits_merge{7}_stat{8}.txt\n".format(
          cmsswbase, os.getcwd(), tagAllSyst, tagNoSyst, 
          " ".join( config.params[ "COMBINE" ][ "FITS" ][ "ARGS" ] ), 
          " ".join( config.params[ "COMBINE" ][ "SIGNIFICANCE" ][ "ARGS" ] ), 
          " ".join( config.params[ "COMBINE" ][ "LIMITS" ][ "ARGS" ] ), 
          config.params[ "MODIFY BINNING" ][ "MIN MERGE" ], 
          str( config.params[ "MODIFY BINNING" ][ "STAT THRESHOLD" ] ).replace( ".", "p" )
        )
      )
      shell.close()
      condor_template( nameLog, nameCondor )
      os.chdir( ".." )

def systematics_correlation_combined():
  tagSmooth = "" if not config.options[ "COMBINE" ][ "SMOOTH" ] else config.params[ "MODIFY BINNING" ][ "SMOOTHING ALGO" ].upper()
  tagABCDnn = "" if not config.options[ "COMBINE" ][ "ABCDNN" ] else "ABCDNN" 
  systematics = [
    "ABCDNNCLOSURE",
    "CFERR1",
    "CFERR2",
    "FSREWK",
    "FSRQCD",
    "FSRSIG",
    "FSRTOP",
    "FSRTTBAR",
    "FSRTTH",
    "HF",
    "HFSTATS1",
    "HFSTATS2",
    "HOTCLOSURE",
    "HOTCSPUR",
    "HOTSTAT",
    "ISREWK",
    "ISRQCD",
    "ISRSIG",
    "ISRTOP",
    "ISRTTBAR",
    "ISRTTH",
    "JECABSOLUTE",
    "JECABSOLUTE20$ERA",
    "JECBBEC1",
    "JECBBEC120$ERA",
    "JECEC2",
    "JECEC220$ERA",
    "JECFLAVORQCD",
    "JECHF",
    "JECHF20$ERA",
    "JECRELATIVEBAL",
    "JECRELATIVESAMPLE20$ERA",
    "JER",
    "LF",
    "LFSTATS1",
    "LFSTATS2",
    "MURFEWK",
    "MURFQCD",
    "MURFSIG",
    "MURFTOP",
    "MURFTTBAR",
    "MURFTTH",
    "PDF",
    "PILEUP",
    "PILEUPJETID",
    "PREFIRE",
    "TOPPT"
  ]
  def write_shell( systName, variable, tag, postfix ):
    os.chdir( "combine" )
    nameLog = "log_{}_{}_{}".format( variable, args.region, tag )
    nameCondor = "SLA_correlation_{}_{}_{}_{}".format( variable, args.region, tag, systName )
    if not os.path.exists( nameLog ): os.system( "mkdir {}".format( nameLog ) )
    shell = open( os.path.join( nameLog, nameCondor + ".sh" ), "w" )
    shell.write(
"#!/bin/bash\n\
source /cvmfs/cms.cern.ch/cmsset_default.sh\n\
cd {0}\n\
eval `scramv1 runtime -sh`\n\
cd {1}\n\
mkdir Results/{3}_{4}_{5}_{6}/correlations\n\
python diffNuisances.py -f html -p {2} Results/{3}_{4}_{5}_{6}/fitDiagnosticsTest.root > correlation_{2}_{3}_{4}_{5}_{6}.html\n\
mv correlation_{2}_{3}_{4}_{5}_{6}.html Results/{3}_{4}_{5}_{6}/correlations/".format( 
      cmsswbase, os.getcwd(), systName, variable, args.region, tag, postfix
      )
    )
    shell.close()
    condor_template( nameLog, nameCondor )
    os.chdir( ".." )

  for tag in args.tags:
    for variable in args.variables:
      for systematic in systematics:
        systName = systematic + tagSmooth
        if systematic in [ "HFSTATS1", "HFSTATS2", "HOTCLOSURE", "HOTCSPUR", "HOTSTAT", "JECABSOLUTE20$ERA", "JECBBEC120$ERA", "JECEC220$ERA", "JECHF20$ERA", "JECRELATIVESAMPLE20$ERA", "LFSTATS1", "LFSTATS2", "PREFIRE" ]:
          for era in [ "16APV", "16", "17", "18" ]:
            systNameEra = systName.replace( "$ERA", era ).replace( "2016APV", "2016" ) # replace 2016APV with 2016 for JEC reduced source splitting
            systNameEra += era
            write_shell( systNameEra, variable, tag, tagABCDnn + tagSmooth )
        else:
          write_shell( systName, variable, tag, tagABCDnn + tagSmooth )

def systematics_correlation():
  trainings = get_trainings( args.tags, args.years, args.variables )
  tagSmooth = "" if not config.options[ "COMBINE" ][ "SMOOTH" ] else config.params[ "MODIFY BINNING" ][ "SMOOTHING ALGO" ].upper()
  tagABCDnn = "" if not config.options[ "COMBINE" ][ "ABCDNN" ] else "ABCDNN" 
  systematics = [
    "ABCDNNCLOSURE",
    "CFERR1",
    "CFERR2",
    "FSREWK",
    "FSRQCD",
    "FSRSIG",
    "FSRTOP",
    "FSRTTBAR",
    "FSRTTH",
    "HF",
    "HFSTATS1",
    "HFSTATS2",
    "HOTCLOSURE",
    "HOTCSPUR",
    "HOTSTAT",
    "ISREWK",
    "ISRQCD",
    "ISRSIG",
    "ISRTOP",
    "ISRTTBAR",
    "ISRTTH",
    "JECABSOLUTE",
    "JECABSOLUTE20$ERA",
    "JECBBEC1",
    "JECBBEC120$ERA",
    "JECEC2",
    "JECEC220$ERA",
    "JECFLAVORQCD",
    "JECHF",
    "JECHF20$ERA",
    "JECRELATIVEBAL",
    "JECRELATIVESAMPLE20$ERA",
    "JER",
    "LF",
    "LFSTATS1",
    "LFSTATS2",
    "MURFEWK",
    "MURFQCD",
    "MURFSIG",
    "MURFTOP",
    "MURFTTBAR",
    "MURFTTH",
    "PDF",
    "PILEUP",
    "PILEUPJETID",
    "PREFIRE",
    "TOPPT"
  ]
  for training in trainings:
    for variable in training[ "variable" ]:
      for systematic in systematics:
        systName = systematic + tagSmooth
        if "$ERA" in systematic:
          systName = systName.replace( "$ERA", training[ "year" ] ).replace( "2016APV", "2016" ) # replace 2016APV with 2016 for JEC reduced source splitting
        if systematic in [ "HFSTATS1", "HFSTATS2", "HOTCLOSURE", "HOTCSPUR", "HOTSTAT", "JECABSOLUTE20$ERA", "JECBBEC120$ERA", "JECEC220$ERA", "JECHF20$ERA", "JECRELATIVESAMPLE20$ERA", "LFSTATS1", "LFSTATS2", "PREFIRE" ]:
          systName += training[ "year" ]
        os.chdir( "combine" )
        nameLog = "log_UL{}_{}_{}_{}".format( training[ "year" ], variable, args.region, training[ "tag" ] )
        nameCondor = "SLA_correlation_{}_{}_{}_{}_{}".format( training[ "year" ], variable, args.region, training[ "tag" ], systName )
        if not os.path.exists( nameLog ): os.system( "mkdir {}".format( nameLog ) )
        shell = open( os.path.join( nameLog, nameCondor + ".sh" ), "w" )
        shell.write(
"#!/bin/bash\n\
source /cvmfs/cms.cern.ch/cmsset_default.sh\n\
cd {0}\n\
eval `scramv1 runtime -sh`\n\
cd {1}\n\
python diffNuisances.py -f html -a -p {2} limits_UL{3}_{4}_{5}_{6}_{7}/fitDiagnosticsTest.root > correlation_UL{3}_{2}.html\n\
mv correlation_UL{3}_{2}.html limits_UL{3}_{4}_{5}_{6}_{7}/cmb/".format( 
  cmsswbase, os.getcwd(), systName, training[ "year" ], variable, args.region, training[ "tag" ], tagABCDnn + tagSmooth
)
        )
        shell.close()
        condor_template( nameLog, nameCondor )
        os.chdir( ".." )

def impact_plots_era():
  trainings = get_trainings( args.tags, args.years, args.variables )
  freezeParams = {
    "NOMINAL": "",
    #"CONSTRAINED": "JECFLAVORQCDLOWESS,MURFTTBARLOWESS,ISRTTBARLOWESS,FSRTTBARLOWESS,ISRTOPLOWESS,TOPPTLOWESS"
    #"FSRTTBAR": "FSRTTBARLOWESS",
    #"THEORYTTH": "ISRTTHLOWESS,FSRTTHLLOWESS,MURFTTHLOWESS,XSEC_TTH",
    #"ABCDNN": "ABCDNNCLOSURE",
    #"THEORYTTBAR": "ISRTTBARHLOWESS,FSRTTBARLOWESS,MURFTTBARLOWESS,XSEC_TTBAR",
    #"THEORYSIG": "ISRSIG,FSRSIGLOWESS,MURFSIGLOWESS"
  }
  systCombo = [
    "", 
    "noShape", 
    "noNorm",
    "noTheory",
    "noShapenoNorm",
    "noShapenoTheory",
    "noNormnoTheory"
  ]
  tagSmooth = "" if not config.options[ "COMBINE" ][ "SMOOTH" ] else config.params[ "MODIFY BINNING" ][ "SMOOTHING ALGO" ].upper()
  tagABCDnn = "" if not config.options[ "COMBINE" ][ "ABCDNN" ] else "ABCDNN"
 
  for training in trainings:
    for variable in training[ "variable" ]:
      for freezeTag in freezeParams:
        if freezeTag != "NOMINAL" and not config.options[ "COMBINE" ][ "IMPACTS" ][ "FREEZE" ]: continue
        for systTag in systCombo:
          if systTag != "" and not config.options[ "COMBINE" ][ "GROUPS" ]: continue
          postfix = tagABCDnn + systTag + tagSmooth + "_merge" + str( config.params[ "MODIFY BINNING" ][ "MIN MERGE" ] ) + "_stat" + str( config.params[ "MODIFY BINNING" ][ "STAT THRESHOLD" ] ).replace( ".", "p" )
          freezeParam = "" if freezeTag == "NOMINAL" else "--freezeParameters {}".format( freezeParams[ freezeTag ] )
          tagFreeze = "noFreeze" if freezeTag == "NOMINAL" else "freeze" + freezeTag
          os.chdir( "combine" )
          nameCondor = "SLA_impacts_{}_{}_{}_{}_{}".format( training[ "year" ], training[ "tag" ], variable, postfix, tagFreeze )
          nameLog = "log_UL{}_{}_{}_{}".format( training[ "year" ], variable, args.region, training[ "tag" ] )
          if not os.path.exists( nameLog ): os.system( "mkdir -vp {}".format( nameLog ) )
          html_line = ""
          if len( args.html ) > 0:
            if not os.path.exists( os.path.join( args.html, "impacts_UL{}_{}_{}_{}".format( training[ "year" ], variable, training[ "tag" ], args.region ) ) ): os.mkdir( os.path.join( args.html, "impacts_UL{}_{}_{}_{}".format( training[ "year" ], variable, training[ "tag" ], args.region ) ) )
            html_line = "cp impacts*.pdf {}".format( os.path.join( args.html, "impacts_UL{}_{}_{}_{}".format( training[ "year" ], variable, training[ "tag" ], args.region ) ) ) 
          shell = open( "{}/{}.sh".format( nameLog, nameCondor ), "w" )
          shell.write(
"#!/bin/bash\n\
source /cvmfs/cms.cern.ch/cmsset_default.sh\n\
cd {0} \n\
eval `scramv1 runtime -sh` \n\
cd {1} \n\
cd limits_UL{2}_{3}_{4}_{5}_{6}/cmb/ \n\
mkdir {7} \n\
cd {7} \n\
combineTool.py -M Impacts -d ../workspace.root -m 125 --doInitialFit --parallel 40 {10}\n\
combineTool.py -M Impacts -d ../workspace.root -m 125 --doFits --parallel 40 --exclude rgx{8} {9} {10} \n\
combineTool.py -M Impacts -d ../workspace.root -m 125 -o impacts_UL{2}_{4}_{3}_{5}_{12}.json --exclude rgx{8} \n\
plotImpacts.py -i impacts_UL{2}_{4}_{3}_{5}_{12}.json -o impacts_UL{2}_{4}_{3}_{5}_{12} \n\
{11} \n".format(
  cmsswbase, os.getcwd(), training[ "year" ], variable, args.region, training[ "tag" ], tagABCDnn + systTag + tagSmooth, tagFreeze, "\{prop_bin.*\}", freezeParam, " ".join( config.params[ "COMBINE" ][ "FITS" ][ "ARGS" ] ), html_line, postfix 
)
          )
          shell.close()
          condor_template( nameLog, nameCondor )
          os.chdir( ".." )
  
def impact_plots_era_masked():
  cPrefix = "TTTX"
  print( "[INFO] Looking for datacards that begin with the prefix {}".format( cPrefix ) )
  trainings = get_trainings( args.tags, args.years, args.variables )
  freezeTags = {
    "NOMINAL": "",
    "FSRTTBAR": "FSRTTBARLOWESS",
    "THEORYTTH": "ISRTTHLOWESS,FSRTTHLLOWESS,MURFTTHLOWESS,XSEC_TTH",
    "ABCDNN": "ABCDNNCLOSURE",
    "THEORYTTBAR": "ISRTTBARHLOWESS,FSRTTBARLOWESS,MURFTTBARLOWESS,XSEC_TTBAR",
    "THEORYSIG": "ISRSIG,FSRSIGLOWESS,MURFSIGLOWESS"
  }
  systCombo = [
    "", 
    "noShape", 
    "noNorm",
    "noTheory",
    "noShapenoNorm",
    "noShapenoTheory",
    "noNormnoTheory"
  ]
  tagSmooth = "" if not config.options[ "COMBINE" ][ "SMOOTH" ] else config.params[ "MODIFY BINNING" ][ "SMOOTHING ALGO" ].upper()
  tagABCDnn = "" if not config.options[ "COMBINE" ][ "ABCDNN" ] else "ABCDNN"
  for training in trainings:
    for variable in training[ "variable" ]:
      for freezeTag in freezeTags:
        if freezeTag != "NOMINAL" and not config.options[ "COMBINE" ][ "IMPACTS" ][ "FREEZE" ]: continue
        freezeParam = "" if freezeTag == "NOMINAL" else "--freezeParameters {}".format( freezeTags[ freezeTag ] )
        tagFreeze = "noFreeze" if freezeTag == "NOMINAL" else "freeze" + freezeTag
        for systTag in systCombo:
          if systTag != "" and not config.options[ "COMBINE" ][ "GROUPS" ]: continue
          postfix = tagABCDnn + systTag + tagSmooth
          dirContent = os.listdir( "combine/limits_UL{}_{}_{}_{}_{}/cmb/".format( training[ "year" ], variable, args.region, training[ "tag" ], postfix ) )
          cardNames = [ item for item in dirContent if ( item.startswith( cPrefix ) and item.endswith( ".txt" ) ) ] 
          print( "[INFO] Found {} datacards in limits_UL{}_{}_{}_{}_{}/cmb/".format( len( cardNames ), training[ "year" ], variable, args.region, training[ "tag" ], postfix ) )
          for cardName in cardNames:
            if "isM" in cardName: continue # mask electron and muon channels together
            tagCard = cardName.split( "." )[0].replace( "isE", "" ).split( "_" )[1]
            os.chdir( "combine" )
            nameLog = "log_UL{}_{}_{}_{}".format( training[ "year" ], variable, args.region, training[ "tag" ] )
            if not os.path.exists( nameLog ): os.system( "mkdir -vp {}".format( nameLog ) )
            nameCondor = "SLA_impacts_{}_{}_{}_{}_{}_mask{}".format( training[ "year" ], training[ "tag" ], variable, postfix, tagFreeze, tagCard )
            shell = open( "{}/{}.sh".format( nameLog, nameCondor ), "w" )
            shell.write(
"#!/bin/bash\n\
source /cvmfs/cms.cern.ch/cmsset_default.sh\n\
cd {0} \n\
eval `scramv1 runtime -sh` \n\
cd {1} \n\
cd limits_UL{2}_{3}_{4}_{5}_{6}/cmb/ \n\
text2workspace.py combined.txt.cmb --channel-masks -o workspace_mask.root \n\
mkdir {7}/mask{8}/ \n\
cd {7}/mask{8}/ \n\
combineTool.py -M Impacts -d ../../workspace_mask.root -m 125 {9} --doInitialFit --parallel 40 --setParameters mask_{10}=1,mask_{11}=1 {13}\n\
combineTool.py -M Impacts -d ../../workspace_mask.root -m 125 {9} --doFits --exclude rgx{12} --parallel 40 --setParameters mask_{10}=1,mask_{11}=1 {13}\n\
combineTool.py -M Impacts -d ../../workspace_mask.root -m 125 -o impacts_UL{2}_{4}_{3}_{5}_{6}_{7}_mask{8}.json --exclude rgx{12}\n\
plotImpacts.py -i impacts_UL{2}_{4}_{3}_{5}_{6}_{7}_mask{8}.json -o ../impacts_UL{2}_{4}_{3}_{5}_{6}_{7}_mask{8} \n".format(
  cmsswbase, os.getcwd(), training[ "year" ], variable, args.region, training[ "tag" ], postfix, tagFreeze, tagCard, " ".join( config.params[ "COMBINE" ][ "FITS" ][ "ARGS" ] ), cardName.split( "." )[0], cardName.replace( "isE", "isM" ).split( "." )[0], "\{prop_bin.*\}", freezeParam 
)
            )
            shell.close()
            condor_template( nameLog, nameCondor )
            os.chdir( ".." )
  
def impact_plots_combined():
  freezeParams = {
    "NOMINAL": "",
    #"JECFLAVORQCD": "JECFLAVORQCDLOWESS",
    #"ISRTTBAR": "ISRTTBARLOWESS",
    #"HF": "HFLOWESS",
    #"ISRTOP": "ISRTOPLOWESS",
    #"MURFTTBAR": "MURFTTBARLOWESS",
    #"ISRQCD": "ISRQCDLOWESS",
    #"FSRTTBAR": "FSRTTBARLOWESS",
    #"THEORYTTH": "ISRTTHLOWESS,FSRTTHLLOWESS,MURFTTHLOWESS,XSEC_TTH"
  }
  tagSmooth = "" if not config.options[ "COMBINE" ][ "SMOOTH" ] else config.params[ "MODIFY BINNING" ][ "SMOOTHING ALGO" ].upper()
  tagABCDnn = "" if not config.options[ "COMBINE" ][ "ABCDNN" ] else "ABCDNN"
  for tag in args.tags:
    for variable in args.variables:
      for freezeTag in freezeParams:
        freezeParam = "" if freezeTag == "NOMINAL" else "--freezeParameters {}".format( freezeParams[ freezeTag ] )
        tagFreeze = "noFreeze" if freezeTag == "NOMINAL" else "freeze" + freezeTag
        os.chdir( "combine" )
        postfix = tagABCDnn + tagSmooth
        nameCondor = "SLA_impacts_{}_{}_{}_{}".format( variable, args.region, tag, postfix + tagFreeze )  
        nameLog = "log_{}_{}_{}".format( variable, args.region, tag )
        if not os.path.exists( nameLog ): os.system( "mkdir -vp {}".format( nameLog ) )
        tagAllSyst = "{}_{}_{}_{}".format( variable, args.region, tag, postfix )
        if not os.path.exists( os.path.join( "Results/", tagAllSyst, tagFreeze ) ): os.system( "mkdir -vp Results/{}".format( os.path.join( tagAllSyst, tagFreeze ) ) )
        shell = open( "{}/{}.sh".format( nameLog, nameCondor ), "w" )
        shell.write(
"#!/bin/bash\n\
source /cvmfs/cms.cern.ch/cmsset_default.sh\n\
cd {0} \n\
eval `scramv1 runtime -sh` \n\
cd {1} \n\
cd Results/{2}/{3}/ \n\
combineTool.py -M Impacts -d ../workspace.root -m 125 --doInitialFit --parallel 40 {4} {5} \n\
combineTool.py -M Impacts -d ../workspace.root -m 125 --doFits --exclude rgx{6} --parallel 40 {4} {5} \n\
combineTool.py -M Impacts -d ../workspace.root -m 125 -o impacts_{2}.json --exclude rgx{6} \n\
plotImpacts.py -i impacts_{2}.json -o impacts_{2} \n\
cd .. \n".format(
            cmsswbase, os.getcwd(), tagAllSyst, tagFreeze, freezeParam,
            " ".join( config.params[ "COMBINE" ][ "FITS" ][ "ARGS" ] ),
            "\{prop_bin.*\}"
          )
        )
        shell.close()
        condor_template( nameLog, nameCondor )
        os.chdir( ".." )

def uncertainty_breakdown_era():
  trainings = get_trainings( args.tags, args.years, args.variables )
  tagSmooth = "" if not config.options[ "COMBINE" ][ "SMOOTH" ] else config.params[ "MODIFY BINNING" ][ "SMOOTHING ALGO" ].upper()
  tagABCDnn = "" if not config.options[ "COMBINE" ][ "ABCDNN" ] else "ABCDNN"
  postfix = tagABCDnn + tagSmooth
  for training in trainings:
    for variable in training[ "variable" ]:
      os.chdir( "combine" )
      nameCondor = "SLA_step8_{}_{}_{}_{}".format( training[ "year" ], training[ "tag" ], variable, postfix )
      nameLog = "log_UL{}_{}_{}_{}".format( training[ "year" ], variable, args.region, training[ "tag" ] )
      if not os.path.exists( nameLog ): os.system( "mkdir -vp {}".format( nameLog ) )
      shell = open( "{}/{}.sh".format( nameLog, nameCondor ), "w" )
      shell.write(
"#!/bin/bash\n\
source /cvmfs/cms.cern.ch/cmsset_default.sh\n\
cd {0} \n\
eval `scramv1 runtime -sh` \n\
cd {1} \n\
cd limits_UL{2}_{3}_{4}_{5}_{6}/cmb/ \n\
combine workspace.root -M MultiDimFit -t -1 --saveWorkspace --X-rtd MINIMIZER_analytic --cminDefaultMinimizerStrategy 0 --cminDefaultMinimizerType Minuit -n postfit \n\
combine workspace.root -M MultiDimFit -t -1 --algo grid --snapshotName MultiDimFit --setParameterRanges r=-100,100 -n workspace.total \n\
combine workspace.root -M MultiDimFit -t -1 --algo grid --snapShotName MultiDimFit --setParameterRanges r=-100,100 --freezeNuisanceGroups SHAPE,NORM -n workspace.freeze_shape_norm \n\
combine workspace.root -M MultiDimFit -t -1 --algo grid --snapShotName MultiDimFit --setParameterRanges r=-100,100 --freezeNuisanceGroups SHAPE,NORM,THEORY -n workspace.freeze_shape_norm_theory \n\
combine workspace.root -M MultiDimFit -t -1 --algo grid --snapShotName MultiDimFit --setParameterRanges r=-100,100 --freezeParameters allConstrainedNuisances -n workspace.freeze_all \n\
plot1DScan.py workspace.root --main-label \"Total Uncertainty\" --others workspace.freeze_shape.MultiDimFit.mH120.root:\"Freeze Shape\":4 workspace.freeze_shape_norm.MultiDimFit.mH120.root:\"Freeze Shape+Norm\":7 workspace.freeze_shape_norm_theory.MultiDimFit.mH120.root:\"Stat Only\":6 --output breakdown_UL{2}_{3}_{4}_{5}_{6} --y-max 10 --y-cut 40 --breakdown \"Shape,Norm,Theory,Rest,Stat\" \n\
cd .. \n".format(
  cmsswbase, os.getcwd(), training[ "year" ], variable, args.region, training[ "tag" ], postfix, "\{prop_bin.*\}"
)
      )
      shell.close()
      condor_template( nameLog, nameCondor )
      os.chdir( ".." )
  

  return

def likelihood_fit_era():
  trainings = get_trainings( args.tags, args.years, args.variables )
  tagSmooth = "" if not config.options[ "COMBINE" ][ "SMOOTH" ] else config.params[ "MODIFY BINNING" ][ "SMOOTHING ALGO" ].upper()
  tagABCDnn = "" if not config.options[ "COMBINE" ][ "ABCDNN" ] else "ABCDNN"
  postfix = tagABCDnn + tagSmooth
  fit_setting = (-20,20,2,2,True)
  fits = { # r-min, r-max, y-cut, y-max, shape
    "ABCDNNCLOSURE":        fit_setting,
    #"CFERR1":               fit_setting,
    #"CFERR2":               fit_setting,
    #"ISREWK":               fit_setting,
    #"FSREWK":               fit_setting,
    "MURFEWK":              fit_setting,
    #"ISRQCD":               fit_setting,
    "FSRQCD":               fit_setting,
    #"MURFQCD":              fit_setting,
    #"ISRTOP":               fit_setting,
    "FSRTOP":               fit_setting,
    #"MURFTOP":              fit_setting,
    "ISRTTBAR":             fit_setting,
    "FSRTTBAR":             fit_setting,
    #"MURFTTBAR":            fit_setting,
    #"ISRTTH":               fit_setting,
    #"FSRTTH":               fit_setting,
    #"MURFTTH":              fit_setting,
    #"HF":                   fit_setting,
    #"HFSTATS1ERA":          fit_setting,
    "HFSTATS2ERA":          fit_setting,
    #"HOTCLOSUREERA":        fit_setting,
    #"HOTCSPURERA":          fit_setting,
    #"HOTSTATERA":           fit_setting,
    #"JECABSOLUTE":          fit_setting,
    #"JECABSOLUTEERA":       fit_setting,
    #"JECBBEC1":             fit_setting,
    #"JECBBEC1ERA":          fit_setting,
    #"JECEC2":               fit_setting,
    #"JECEC2ERA":            fit_setting,
    #"JECFLAVORQCD":         fit_setting,
    #"JECHF":                fit_setting,
    #"JECHFERA":             fit_setting,
    #"JECRELATIVEBAL":       fit_setting,
    #"JECRELATIVESAMPLEERA": fit_setting,
    #"JERERA":               fit_setting,
    #"LF":                   fit_setting,
    #"LFSTATS1ERA":          fit_setting,
    #"LFSTATS2ERA":          fit_setting,
    #"PDF":                  fit_setting,
    #"PILEUP":               fit_setting,
    "PREFIREERA":           fit_setting,
    #"PILEUPJETID":          fit_setting,
    #"TOPPT":                fit_setting
  }
  for training in trainings:
    for variable in training[ "variable" ]:
      for fitVar in fits:
        fitName = fitVar + config.params[ "MODIFY BINNING" ][ "SMOOTHING ALGO" ].upper() if ( config.options[ "COMBINE" ][ "SMOOTH" ] and fits[ fitVar ][-1] ) else fitVar
        if "PREFIRE" in fitName and training[ "year" ] not in [ "16APV", "16", "17" ]: continue
        if "ERA" in fitVar:
          if "JEC" in fitVar:
            fitName = fitName.replace( "ERA", "20" + training[ "year" ] ).replace( "APV", "" ) + training[ "year" ]
          else:
            fitName = fitName.replace( "ERA", "" ) + training[ "year" ]
        nameCondor = "SLA_MDFit_{}_{}_{}_{}_{}".format( training[ "year" ], training[ "tag" ], variable, postfix, fitName )
        nameLog = "log_UL{}_{}_{}_{}".format( training[ "year" ], variable, args.region, training[ "tag" ] )
        os.chdir( "combine" ) 
        if not os.path.exists( nameLog ): os.system( "mkdir -vp {}".format( nameLog ) )
        shell = open( "{}/{}.sh".format( nameLog, nameCondor ), "w" )
        shell.write(
"#!/bin/bash\n\
source /cvmfs/cms.cern.ch/cmsset_default.sh\n\
cd {0} \n\
eval `scramv1 runtime -sh` \n\
cd {1}/limits_UL{2}_{3}_{4}_{5}_{6}/cmb/ \n\
combine workspace.root -M MultiDimFit --algo grid --floatOtherPOIs 1 --saveInactivePOI 1 --autoRange 1.2 --squareDistPoiStep --points 100 {7} -n _scan_{8} -P {8} --setParameterRange r={9},{10}:{8}=-3.5,4\n\
plot1DScan.py higgsCombine_scan_{8}.MultiDimFit.mH125.root --main-label {8} -o scan_UL{2}_{3}_{4}_{5}_{6}_{8} --y-cut {11} --y-max {12}\n\
mkdir scan_UL{2}_{3}_{4}_{5}_{6}/ \n\
mv scan_UL{2}_{3}_{4}_{5}_{6}_{8}.* scan_UL{2}_{3}_{4}_{5}_{6}/ \n\
rm *{8}*".format(
          cmsswbase, os.getcwd(), training[ "year" ], variable, args.region, training[ "tag" ], postfix, " ".join( config.params[ "COMBINE" ][ "FITS" ][ "ARGS" ] ), fitName, fits[ fitVar ][0], fits[ fitVar ][1], fits[ fitVar ][2], fits[ fitVar ][3] 
          )
        ) 
        shell.close()
        condor_template( nameLog, nameCondor )
        os.chdir( "../" )

  return

def likelihood_fit_combined():
  tagSmooth = "" if not config.options[ "COMBINE" ][ "SMOOTH" ] else config.params[ "MODIFY BINNING" ][ "SMOOTHING ALGO" ].upper()
  tagABCDnn = "" if not config.options[ "COMBINE" ][ "ABCDNN" ] else "ABCDNN"
  postfix = tagABCDnn + tagSmooth
  fit_setting = (-7,7,2,2,True)
  fits = { # r-min, r-max, y-cut, y-max, shape
    "ABCDNNCLOSURE":         fit_setting,
    "CFERR1":                fit_setting,
    "CFERR2":                fit_setting,
    "ISREWK":                fit_setting,
    "FSREWK":                fit_setting,
    "MURFEWK":               fit_setting,
    "ISRQCD":                fit_setting,
    "FSRQCD":                fit_setting,
    "MURFQCD":               fit_setting,
    "ISRTOP":                fit_setting,
    "FSRTOP":                fit_setting,
    "MURFTOP":               fit_setting,
    "ISRTTBAR":              fit_setting,
    "FSRTTBAR":              fit_setting,
    "MURFTTBAR":             fit_setting,
    "ISRTTH":                fit_setting,
    "FSRTTH":                fit_setting,
    "MURFTTH":               fit_setting,
    "HF":                    fit_setting,
    "HFSTATS1":              fit_setting,
    "HFSTATS2":              fit_setting,
    "HOTCLOSURE":            fit_setting,
    "HOTCSPUR":              fit_setting,
    "HOTSTAT":               fit_setting,
    "JECABSOLUTE":           fit_setting,
    "JECABSOLUTE":           fit_setting,
    "JECBBEC1":              fit_setting,
    "JECBBEC1":              fit_setting,
    "JECEC2":                fit_setting,
    "JECEC2$ERA":            fit_setting,
    "JECFLAVORQCD":          fit_setting,
    "JECHF":                 fit_setting,
    "JECHF$ERA":             fit_setting,
    "JECRELATIVEBAL":        fit_setting,
    "JECRELATIVESAMPLE$ERA": fit_setting,
    "JER":                   fit_setting,
    "LF":                    fit_setting,
    "LFSTATS1":              fit_setting,
    "LFSTATS2":              fit_setting,
    "PDF":                   fit_setting,
    "PILEUP":                fit_setting,
    "PREFIRE":               fit_setting,
    "PILEUPJETID":           fit_setting,
    "TOPPT":                 fit_setting
  }
  def write_shell( fitName, fitVar, fits, tag, variable, postfix ):
    nameCondor = "SLA_MDFit_{}_{}_{}_{}".format( variable, tag, postfix, fitName )
    nameLog = "log_{}_{}_{}".format( variable, args.region, tag )
    os.chdir( "combine" ) 
    if not os.path.exists( nameLog ): os.system( "mkdir -vp {}".format( nameLog ) )
    shell = open( "{}/{}.sh".format( nameLog, nameCondor ), "w" )
    shell.write(
"#!/bin/bash\n\
source /cvmfs/cms.cern.ch/cmsset_default.sh\n\
cd {0} \n\
eval `scramv1 runtime -sh` \n\
cd {1}/Results/{2}_{3}_{4}_{5}/ \n\
combine workspace.root -M MultiDimFit --algo grid --floatOtherPOIs 1 --saveInactivePOI 1 --autoRange 1.5 --squareDistPoiStep --points 200 {6} -n _scan_{7} -P {7} --setParameterRange r={8},{9}:{10}=-6,6\n\
plot1DScan.py higgsCombine_scan_{7}.MultiDimFit.mH125.root --main-label {7} -o scan_{2}_{3}_{4}_{5}_{7} --y-cut {10} --y-max {11}\n\
mkdir scan_{2}_{3}_{4}_{5}/ \n\
mv scan_{2}_{3}_{4}_{5}_{7}.* scan_{2}_{3}_{4}_{5}/ \n\
rm *{8}*".format(
      cmsswbase, os.getcwd(), variable, args.region, tag, postfix, 
      " ".join( config.params[ "COMBINE" ][ "FITS" ][ "ARGS" ] ), 
      fitName, fits[ fitVar ][0], fits[ fitVar ][1], fits[ fitVar ][2], fits[ fitVar ][3] 
      )
    ) 
    shell.close()
    condor_template( nameLog, nameCondor )
    os.chdir( "../" )

  for tag in args.tags:
    for variable in args.variables:
      for fitVar in fits:
        fitName = fitVar + config.params[ "MODIFY BINNING" ][ "SMOOTHING ALGO" ].upper() if ( config.options[ "COMBINE" ][ "SMOOTH" ] and fits[ fitVar ][-1] ) else fitVar
        if fitVar in [ "HFSTATS1", "HFSTATS2", "HOTCLOSURE", "HOTCSPUR", "HOTSTAT", "JECABSOLUTE20$ERA", "JECBBEC120$ERA", "JECEC220$ERA", "JECHF20$ERA", "JECRELATIVESAMPLE20$ERA", "LFSTATS1", "LFSTATS2", "PREFIRE" ]:
          for era in [ "16APV", "16", "17", "18" ]:
            if "PREFIRE" in fitName and era not in [ "16APV", "16", "17" ]: continue
            fitNameEra = fitName.replace( "$ERA", era ).replace( "2016APV", "2016" )
            fitNameEra += era
            write_shell( fitNameEra, fitVar, fits, tag, variable, postfix )
        else:
          write_shell( fitName, fitVar, fits, tag, variable, postfix )


if args.step == "1": make_templates()
elif args.step == "2": format_templates()
elif args.step == "3": plot_templates()
elif args.step == "4": run_combine()
elif args.step == "5": combine_years()
elif args.step == "impact":
  if "Run2" in args.years:
    impact_plots_combined()
  if "16APV" in args.years or "16" in args.years or "17" in args.years or "18" in args.years:
    impact_plots_era()
    if config.options[ "COMBINE" ][ "IMPACTS" ][ "MASKED" ]:
      impact_plots_era_masked()
elif args.step == "correlation":
  if "Run2" in args.years:
    systematics_correlation_combined()
  if "16APV" in args.years or "16" in args.years or "17" in args.years or "18" in args.years:
    systematics_correlation()
elif args.step == "fit":
  if "Run2" in args.years:
    likelihood_fit_combined()
  if "16APV" in args.years or "16" in args.years or "17" in args.years or "18" in args.years:
    likelihood_fit_era()
elif args.step == "breakdown":
  if "Run2" in args.years:
    uncertainty_breakdown_combined()
  if "16APV" in args.years or "16" in args.years or "17" in args.years or "18" in args.years:
    uncertainty_breakdown_era()
  
else:
  print( "[ERR] Invalid step option used" ) 
