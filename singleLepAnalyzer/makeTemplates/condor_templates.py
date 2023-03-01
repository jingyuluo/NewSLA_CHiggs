#!/usr/bin/python

import os,sys,datetime,itertools
from argparse import ArgumentParser
sys.path.append( "../" )
sys.path.append( "../singleLepAnalyzer/" )
import config
import utils

parser = ArgumentParser()
parser.add_argument( "-y", "--year", required = True, help = "[16,17,18]" )
parser.add_argument( "-v", "--variables", nargs = "+", required = True )
parser.add_argument( "-p", "--postfix", default = "test" )
parser.add_argument( "-r", "--region", required = True, help = "[SR,PS,TTCR,WJCR]" )
args = parser.parse_args()

thisDir = os.getcwd()

if args.region == "TTCR": 
  bins = config.hist_bins[ "VR" ]
elif args.region == "WJCR":
  bins = config.hist_bins[ "VR" ]
elif args.region == "SR":
  bins = config.hist_bins[ "SR" ]
elif args.region == "VR":
  bins = config.hist_bins[ "VR" ]
elif args.region == "BASELINE":
  bins = config.hist_bins[ "BASELINE" ]
elif args.region == "ABCDNN":
  bins = config.hist_bins[ "ABCDNN" ]
else:
  quit( "[ERR] Invalid region argument used. Quitting." )

categories = list(
  itertools.product(
    bins[ "LEPTON" ],
    bins[ "NHOT" ],
    bins[ "NT" ],
    bins[ "NW" ],
    bins[ "NB" ],
    bins[ "NJ" ]
  )
)
	
subDir = "{}_UL{}_{}".format( config.region_prefix[ args.region ], args.year, args.postfix )
outputPath = os.path.join( os.getcwd(), subDir )
if not os.path.exists( outputPath ): os.system( "mkdir -vp {}".format( outputPath ) )

nJobs = 0
for variable in args.variables:
  print( ">> Generating templates for {}".format( variable ) )
  for category in categories:
    categoryTag = "is{}nHOT{}nT{}nW{}nB{}nJ{}".format( 
      category[0],
      category[1],
      category[2],
      category[3],
      category[4],
      category[5]
    )
    if ( int(category[1][0]) + int(category[2][0]) + int(category[3][0]) + int(category[4][0]) ) > int(category[5][0]):
      print( "[WARN] {} is not topologically possible, skipping...".format( categoryTag ) )
      continue
      
    if not os.path.exists( os.path.join( outputPath, categoryTag ) ): os.system( "mkdir -vp {}".format( os.path.join( outputPath, categoryTag ) ) )
    os.chdir( os.path.join( outputPath, categoryTag ) )

    jobParams = {
      "VARIABLE": variable,
      "YEAR": args.year,
      "LEPTON": category[0],
      "NHOT":  category[1],
      "NT": category[2],
      "NW": category[3],
      "NB": category[4],
      "NJ": category[5],
      "EXEDIR": thisDir,
      "SUBDIR": subDir,
      "POSTFIX": args.postfix
    }

    jdf = open( "condor_step1_{}.job".format( variable ), "w" )
    jdf.write(
"""universe = vanilla
Executable = %(EXEDIR)s/condor_templates.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
request_memory = 5000
Output = condor_step1_%(VARIABLE)s.out
Error = condor_step1_%(VARIABLE)s.err
Log = condor_step1_%(VARIABLE)s.log
JobBatchName = SLA_step1_%(YEAR)s_%(VARIABLE)s_%(POSTFIX)s
Notification = Error
Arguments = %(VARIABLE)s %(YEAR)s %(LEPTON)s %(NHOT)s %(NT)s %(NW)s %(NB)s %(NJ)s %(EXEDIR)s %(SUBDIR)s
Queue 1"""%jobParams
    )
    jdf.close()
    os.system( "condor_submit condor_step1_{}.job".format( variable ) )
    os.chdir( ".." )
    nJobs += 1
    if config.options[ "GENERAL" ][ "TEST" ]: 
      print( "[OPT] Testing one job." )
      break

print( "[DONE] Total jobs submitted: {}".format( nJobs ) )
