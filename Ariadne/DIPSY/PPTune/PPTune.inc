# HOW TO TUNE
# 
# We need to generate a number of yoda files for different energies
# different observables and different parameter settings. This can be
# done in an semi-automated way.
#
# For the total and elastic cross sections we do not need to generate
# final states and these runs can in principle be done once and for
# all. In fact running PP16XSECYM??.run (generated in the SAVERUNYM
# just after the line "# 16 inclusive observables" below) will produce
# all necessary data.
#
# This input file is (to make a gross understatement) a terrible
# mess. It needs to beprocessed by the script ./expand.pl before going
# to setupThePEG and it will produce an awfull amount of .run
# files. However doing eg.
# > make blaha
# in this directory yo will only produce the .run files matching "blaha".
# 
# For the final state you have to run the files generated after the
# line "# 16 final state observables" below. You will have one run per
# energy in PP16XSECYM@??.run (where @ is some arbitrary letter used
# to differentiate between settings). The 7TeV will take quite a while
# so it is best to divide it up in smaller runs. Using the scripts
# "sub" and "runtag" in ~leif/bin you can start the programs as
# follows:
#
# > sub runtag ./runThePEG, PP16XSECYM@02, PP16XSECYM@09,
#                           PP16XSECYM@71_61-80, PP16XSECYM@71_41-60,
#                           PP16XSECYM@71_21-40, PP16XSECYM@71_1-20
#                           PP16XSECYM@71_81-100
#
# After these runs are done you can start the actual tuning with
# professor. First you need to copy the relevant files from the
# inclusive runs to match the names of the final states runs as
# follows:
# > cp PP16XSECYM01.log PP16XSECYM@01.log
# > cp PP16XSECYM05.log PP16XSECYM@05.log
# > cp PP16XSECYM18.log PP16XSECYM@18.log
# > cp PP16XSECYM70.log PP16XSECYM@70.log
# then run the mkprof.pl script as follows to get a "global" tune
#
# > ./mkprof.pl -t -w PPTuneFS00.weights -o PP16XSECYM@00.prof PP16XSECYM@
# preferably redirecting both stdout and stderr to some file (there
# will be a lot of output). Similarly
# > ./mkprof.pl -t -w PPTuneFS70.weights -o PP16XSECYM@70.prof PP16XSECYM@

# will give the tune restricting the final state data to 7TeV, and
# changing all "70" to "09" or "02" will give 900 GeV or 200 GeV to
# look at the result check the file
# PP16XSECYM@00.prof/tunes/results.pkl or use the following script
# > ./profsum.pl PP16XSECYM@00.prof/tunes/results.pkl
# To get just a summary.
#
# After the line "# 16 pipeing the tunes" below you can then pipe the
# different tunes, just set the parameters according to the fit, do
# > make PPT16XSECYM@
# and run eg.
# > sub runtag ./runThePEG, PPT16XSECYM@02, PPT16XSECYM@09, PPT16XSECYM@71,
#                           PPT16XSECYM@0202, PPT16XSECYM@0909, PPT16XSECYM@7170
#
#
#
#
cd /DIPSY

## First we setup some previously tuned tuned parameters
#read CurrentTune.in
read ../Tune31.in

## Now we set up an event generator.
cp EventHandler PPEventHandler
set stdProton:R0 0.0
set stdAntiProton:R0 0.0
set Proton:R0 0.0
set AntiProton:R0 0.0
set PPEventHandler:WFL stdProton
set PPEventHandler:WFR stdProton
set PPEventHandler:ConsistencyLevel 0
set PPEventHandler:XSecFn:CheckOffShell false
set PPEventHandler:CascadeHandler NULL
set PPEventHandler:HadronizationHandler NULL
set PPEventHandler:DecayHandler NULL
create ThePEG::FixedCMSLuminosity PPLumi
set PPEventHandler:LuminosityFunction PPLumi
cp Generator PPGenerator
erase PPGenerator:AnalysisHandlers[0]
set PPGenerator:HistogramFactory NULL
set PPGenerator:EventHandler PPEventHandler
set PPGenerator:NumberOfEvents 0
set PPGenerator:EventHandler:EventFiller:PTCut 0.6
set PPEventHandler:BGen:Width 5
set PPEventHandler:EffectivePartonMode Colours

## These are the analysess we will run
## Some semi-inclusive cross section for DIPSY which need at least
## four combinations of left- and right-moving cascades.
erase PPEventHandler:AnalysisHandlers[0]
erase PPEventHandler:AnalysisHandlers[0]
create DIPSY::SemiInclusiveXSecAnalysis SemiIncl SemiInclusiveXSecAnalysis.so
insert PPEventHandler:AnalysisHandlers[0] SemiIncl
set PPEventHandler:PreSampleL 1
set PPEventHandler:PreSampleR 1
## This is just to keep track of the progress of a run
create DIPSY::AnalysisProgress AnaLog AnalysisProgress.so
set AnaLog:Interval 600
insert PPEventHandler:AnalysisHandlers[0] AnaLog

## The sample rates need to be adjusted so that we get a reasonable
## statistics in a reasonable time. It is typically efficient to
## sample a number of impact parameter values for each pair of DIPSY
## cascades.
set PPEventHandler:PreSampleB 1
## But we need a good sample of cascades
set PPEventHandler:PreSamples 1000

## We need the same set of parameters for all different energies, so
## we use a separate random generator.
cp  /Defaults/Random RandomArg
set PPGenerator:SeparateRandom RandomArg

## These are the parameters we want to tune
do PPGenerator:AddRndInterface /DIPSY/PPEventHandler:RMax 100 2.0 4.0 2.9 0.3
do PPGenerator:AddRndInterface /DIPSY/PPEventHandler:LambdaQCD 1 0.10 0.3 0.22 0.04
do PPGenerator:AddRndInterface /DIPSY/PPEventHandler:BaryonSize 1 2.0 4.0 2.9 0.3
do PPGenerator:AddRndInterface /DIPSY/PPEventHandler:Emitter:PMinusOrdering 1 0.5 2.0 1.0 0.2
do PPGenerator:AddRndInterface /DIPSY/PPEventHandler:Emitter:PSInflation 1 0.5 2.0 1.0 0.2

## Now we want to run for all energies and two values for YFrametest

SAVERUN PP01XSEC PPGenerator

## Set the pp energy and interaction frame we want to run with

do PPGenerator:RemoveInterface /DIPSY/PPEventHandler:BaryonSize
set /DIPSY/PPEventHandler:BaryonSize 0.0
do PPGenerator:AddRndInterface /DIPSY/PPEventHandler:Emitter:PSInflation 1 0.5 2.0 1.0 0.2

SAVERUN PP02XSEC PPGenerator

do PPGenerator:RemoveInterface /DIPSY/PPEventHandler:Emitter:PSInflation
do PPGenerator:AddRndInterface /DIPSY/PPEventHandler:Swinger:Lambda 1 0.5 4.0 1.0 0.5

SAVERUN PP05XSEC PPGenerator

do PPGenerator:RemoveInterface /DIPSY/PPEventHandler:Swinger:Lambda
do PPGenerator:AddRndInterface /DIPSY/PPEventHandler:Swinger:Lambda 2 0.0 6.0

SAVERUN PP06XSEC PPGenerator

do PPGenerator:RemoveInterface /DIPSY/PPEventHandler:Swinger:Lambda
set PPEventHandler:Swinger:Lambda 5.0
do PPGenerator:RemoveInterface /DIPSY/PPEventHandler:RMax
do PPGenerator:AddRndInterface /DIPSY/PPEventHandler:RMax 100 1.0 5.0 2.9 0.6
SAVERUN PP07XSEC PPGenerator
do PPGenerator:AddRndInterface /DIPSY/PPEventHandler:BaryonSize 1 2.0 4.0 2.9 0.5
SAVERUN PP08XSEC PPGenerator
set PPEventHandler:Swinger:Lambda 2.0

# 09

SAVERUN PP09XSEC PPGenerator
SAVERUNYR PP09XSECYR PPGenerator

set PPGenerator:EventHandler:DecayHandler /Defaults/Handlers/StandardDecayHandler
set PPGenerator:EventHandler:HadronizationHandler Frag8
set PPGenerator:EventHandler:CascadeHandler AriadneCascade
set PPGenerator:EventHandler:PreSamples 0
set PPGenerator:NumberOfEvents 10000
create ThePEG::ProgressLog Logger ProgressLog.so
set Logger:Interval 600
insert PPGenerator:AnalysisHandlers[0] Logger
create ThePEG::RivetAnalysis RivetTune RivetAnalysis.so
insert RivetTune:Paths[0] .
insert RivetTune:Analyses[0] ATLAS_2010_S8918562
insert RivetTune:Analyses[0] STAR_2008_S7869363
insert RivetTune:Analyses[0] CMS_2011_S8978280
insert RivetTune:Analyses[0] PYTHIA_TUNING
insert PPGenerator:AnalysisHandlers[0] RivetTune
set PPGenerator:EventHandler:LuminosityFunction:Energy 900
saverun PP09XSECYR09 PPGenerator
set PPGenerator:EventHandler:LuminosityFunction:Energy 7000
saverun PP09XSECYR71 PPGenerator
set PPGenerator:EventHandler:YFrametest 0.5
set PPGenerator:EventHandler:LuminosityFunction:Energy 900
saverun PP09XSEC09 PPGenerator
set PPGenerator:EventHandler:LuminosityFunction:Energy 7000
saverun PP09XSEC71 PPGenerator
saverun dPP09XSEC71 PPGenerator

# 10

set Frag8:StringZ_aLund 0.42
set Frag8:StringZ_bLund 0.40
set Frag8:StringPT_sigma 0.32
set Frag8:StringFlav_probQQtoQ 0.084
set Frag8:StringFlav_probStoUD 0.22
set Frag8:Collapser /Ariadne5/Defaults/Collapser
cp Frag8 Frag8Sw0
set Frag8Sw0:StringZ_aLund 0.30
set Frag8Sw0:StringZ_bLund 0.36
set Frag8Sw0:StringPT_sigma 0.32
set Frag8Sw0:StringFlav_probQQtoQ 0.082
set Frag8Sw0:StringFlav_probStoUD 0.22
cp Frag8 Frag8Sw2
set Frag8Sw2:StringZ_aLund 0.45
set Frag8Sw2:StringZ_bLund 0.42
set Frag8Sw2:StringPT_sigma 0.32
set Frag8Sw2:StringFlav_probQQtoQ 0.085
set Frag8Sw2:StringFlav_probStoUD 0.22
cp Frag8 Frag8Rope
set Frag8Rope:FragmentationScheme dipole
set Frag8Rope:StringR0 1.0
set Frag8Rope:Stringm0 0.2
set Frag8Rope:Average false
set Frag8Rope:ThrowAway true
set Frag8Rope:BaryonSuppression 0.25
set Frag8Rope:StringZ_aLund 0.41
set Frag8Rope:StringZ_bLund 0.37
set Frag8Rope:StringPT_sigma 0.31
set Frag8Rope:StringFlav_probQQtoQ 0.073
set Frag8Rope:StringFlav_probStoUD 0.21

set PPEventHandler:WFL Proton
set PPEventHandler:WFR Proton
set PPEventHandler:EventFiller:SoftRemove NoValence

set PPGenerator:EventHandler:YFrametest 0.5
set PPGenerator:EventHandler:LuminosityFunction:Energy 900
saverun PP10XSEC09 PPGenerator
set PPGenerator:EventHandler:LuminosityFunction:Energy 7000
saverun PP10XSEC71 PPGenerator

set PPGenerator:EventHandler:YFrametest -1.2
set PPGenerator:EventHandler:LuminosityFunction:Energy 900
saverun PP10XSECYR09 PPGenerator
set PPGenerator:EventHandler:LuminosityFunction:Energy 7000
saverun PP10XSECYR71 PPGenerator

set Frag8:FragmentationScheme pipe
set PPGenerator:EventHandler:LuminosityFunction:Energy 900
saverun PP10XSECYRR09 PPGenerator
set PPGenerator:EventHandler:LuminosityFunction:Energy 7000
saverun PP10XSECYRR71 PPGenerator
set Frag8:FragmentationScheme none

# 11

set PPEventHandler:EventFiller:PTCut 1.0
set PPEventHandler:Swinger:Lambda 1.0
set PPEventHandler:EffectivePartonMode Colours
do PPGenerator:RemoveInterface /DIPSY/PPEventHandler:BaryonSize
set /DIPSY/PPEventHandler:BaryonSize 0.0
do PPGenerator:AddRndInterface /DIPSY/PPEventHandler:CoherenceRange 1 0.2 5.0 3.0 2.0
erase PPGenerator:AnalysisHandlers[0]
erase PPGenerator:AnalysisHandlers[0]
set PPGenerator:EventHandler:PreSamples 1000
set PPGenerator:NumberOfEvents 0
SAVERUNYR PP11XSECYR PPGenerator
SAVERUNYM PP11XSECYM PPGenerator
set PPGenerator:EventHandler:FudgeME true
SAVERUNYM PP11XSECYm PPGenerator
set PPGenerator:EventHandler:Emitter:PTScale 2.0
SAVERUNYM PP11XSECYP PPGenerator
set PPGenerator:EventHandler:Emitter:PTScale 0.5
SAVERUNYM PP11XSECYp PPGenerator
set PPGenerator:EventHandler:Emitter:PTScale 1.4
SAVERUNYM PP11XSECYs PPGenerator
set PPGenerator:EventHandler:Emitter:PTScale 1.0
set PPGenerator:EventHandler:FudgeME false
insert PPGenerator:AnalysisHandlers[0] Logger
insert PPGenerator:AnalysisHandlers[0] RivetTune
set PPGenerator:EventHandler:PreSamples 0
set PPGenerator:NumberOfEvents 10000


set PPGenerator:EventHandler:YFrametest -1.2
set PPGenerator:EventHandler:LuminosityFunction:Energy 900
saverun PP11XSECYR09 PPGenerator
set PPGenerator:EventHandler:LuminosityFunction:Energy 7000
saverun PP11XSECYR71 PPGenerator

set Frag8:FragmentationScheme pipe
set PPGenerator:EventHandler:LuminosityFunction:Energy 900
saverun PP11XSECYRR09 PPGenerator
set PPGenerator:EventHandler:LuminosityFunction:Energy 7000
saverun PP11XSECYRR71 PPGenerator
set Frag8:FragmentationScheme none

set PPGenerator:EventHandler:YFrametest 0.5
set PPGenerator:EventHandler:LuminosityFunction:Energy 900
saverun PP11XSECYM09 PPGenerator
set PPGenerator:EventHandler:LuminosityFunction:Energy 7000
saverun PP11XSECYM71 PPGenerator
set PPGenerator:EventHandler:FudgeME true
set PPGenerator:EventHandler:LuminosityFunction:Energy 900
saverun PP11XSECYm09 PPGenerator
set PPGenerator:EventHandler:LuminosityFunction:Energy 7000
saverun PP11XSECYm71 PPGenerator
set PPGenerator:EventHandler:Emitter:PTScale 2.0
set PPGenerator:EventHandler:LuminosityFunction:Energy 900
saverun PP11XSECYP09 PPGenerator
set PPGenerator:EventHandler:LuminosityFunction:Energy 7000
saverun PP11XSECYP71 PPGenerator
set PPGenerator:EventHandler:Emitter:PTScale 1.0
set PPGenerator:EventHandler:FudgeME false

set Frag8:FragmentationScheme pipe
set PPGenerator:EventHandler:LuminosityFunction:Energy 900
saverun PP11XSECYMR09 PPGenerator
set PPGenerator:EventHandler:LuminosityFunction:Energy 7000
saverun PP11XSECYMR71 PPGenerator
set PPGenerator:EventHandler:FudgeME true
set PPGenerator:EventHandler:LuminosityFunction:Energy 900
saverun PP11XSECYmR09 PPGenerator
set PPGenerator:EventHandler:LuminosityFunction:Energy 7000
saverun PP11XSECYmR71 PPGenerator
set PPGenerator:EventHandler:Emitter:PTScale 2.0
set PPGenerator:EventHandler:LuminosityFunction:Energy 900
saverun PP11XSECYPR09 PPGenerator
set PPGenerator:EventHandler:LuminosityFunction:Energy 7000
saverun PP11XSECYPR71 PPGenerator
set PPGenerator:EventHandler:Emitter:PTScale 0.5
set PPGenerator:EventHandler:LuminosityFunction:Energy 900
saverun PP11XSECYpR09 PPGenerator
set PPGenerator:EventHandler:LuminosityFunction:Energy 7000
saverun PP11XSECYpR71 PPGenerator
set PPGenerator:EventHandler:Emitter:PTScale 1.4
set PPGenerator:EventHandler:LuminosityFunction:Energy 900
saverun PP11XSECYsR09 PPGenerator
set PPGenerator:EventHandler:LuminosityFunction:Energy 7000
saverun PP11XSECYsR71 PPGenerator
set PPGenerator:EventHandler:Emitter:PTScale 1.0
set PPGenerator:EventHandler:FudgeME false

set FSSwinger:SetRmax -2.0
set PPGenerator:EventHandler:LuminosityFunction:Energy 900
saverun PP11XSECYMr09 PPGenerator
set PPGenerator:EventHandler:LuminosityFunction:Energy 7000
saverun PP11XSECYMr71 PPGenerator
set FSSwinger:SetRmax -2.7
set Frag8:FragmentationScheme none
erase AriadneCascade:Emitters[0]
set PPGenerator:EventHandler:LuminosityFunction:Energy 900
saverun PP11XSECYN09 PPGenerator
set PPGenerator:EventHandler:LuminosityFunction:Energy 7000
saverun PP11XSECYN71 PPGenerator
insert AriadneCascade:Emitters[0] FSSwinger
set PPGenerator:EventHandler:CascadeHandler NULL
set PPGenerator:EventHandler:LuminosityFunction:Energy 900
saverun PP11XSECYn09 PPGenerator
set PPGenerator:EventHandler:LuminosityFunction:Energy 7000
saverun PP11XSECYn71 PPGenerator
set PPGenerator:EventHandler:CascadeHandler AriadneCascade

cp PPGenerator PPTestTune
do PPTestTune:RemoveInterface /DIPSY/PPEventHandler:CoherenceRange
do PPTestTune:RemoveInterface /DIPSY/PPEventHandler:LambdaQCD
do PPTestTune:RemoveInterface /DIPSY/PPEventHandler:RMax
do PPTestTune:RemoveInterface /DIPSY/stdEmitter:PMinusOrdering
set PPTestTune:NumberOfEvents 100000

set PPEventHandler:CoherenceRange 2.7
set PPEventHandler:LambdaQCD 0.26
set PPEventHandler:RMax 2.7 
set PPEventHandler:Emitter:PMinusOrdering 0.73
set PPTestTune:EventHandler:LuminosityFunction:Energy 900
saverun PPT11XSECYM09 PPTestTune
set PPTestTune:EventHandler:LuminosityFunction:Energy 7000
saverun PPT11XSECYM71 PPTestTune

set PPTestTune:NumberOfEvents 10000
saverun PPT11XSECYM710 PPTestTune
saverun PPT11XSECYM711 PPTestTune
set PPTestTune:EventHandler:FudgeME true
saverun PPT11XSECYM712 PPTestTune
saverun PPT11XSECYM713 PPTestTune
set PPEventHandler:CoherenceRange 2.2
set PPEventHandler:LambdaQCD 0.27
set PPEventHandler:RMax 3.0
set PPEventHandler:Emitter:PMinusOrdering 0.77
set PPTestTune:NumberOfEvents 100000
saverun PPT11XSECYm71 PPTestTune

set PPTestTune:EventHandler:FudgeME false
set PPEventHandler:CoherenceRange 2.7
set PPEventHandler:LambdaQCD 0.26
set PPEventHandler:RMax 2.7 
set PPEventHandler:Emitter:PMinusOrdering 0.73

set PPEventHandler:CoherenceRange 3.0
set PPEventHandler:LambdaQCD 0.27
set PPEventHandler:RMax 2.7 
set PPEventHandler:Emitter:PMinusOrdering 0.76
set Frag8:FragmentationScheme pipe
set PPTestTune:EventHandler:LuminosityFunction:Energy 900
saverun PPT11XSECYMR09 PPTestTune
set PPTestTune:EventHandler:LuminosityFunction:Energy 7000
saverun PPT11XSECYMR71 PPTestTune

set PPTestTune:EventHandler:FudgeME true
set PPEventHandler:CoherenceRange 2.5
set PPEventHandler:LambdaQCD 0.27
set PPEventHandler:RMax 3.0
set PPEventHandler:Emitter:PMinusOrdering 0.81
saverun PPT11XSECYmR71 PPTestTune
set PPEventHandler:CoherenceRange 2.4
set PPEventHandler:LambdaQCD 0.12
set PPEventHandler:RMax 3.3
set PPEventHandler:Emitter:PMinusOrdering 0.93
set PPGenerator:EventHandler:Emitter:PTScale 2.0
saverun PPT11XSECYPR71 PPTestTune
set PPEventHandler:CoherenceRange 2.5
set PPEventHandler:LambdaQCD 0.18
set PPEventHandler:RMax 3.2
set PPEventHandler:Emitter:PMinusOrdering 0.71
set PPGenerator:EventHandler:Emitter:PTScale 1.4
saverun PPT11XSECYsR71 PPTestTune

set PPTestTune:EventHandler:FudgeME false
set PPEventHandler:CoherenceRange 3.0
set PPEventHandler:LambdaQCD 0.27
set PPEventHandler:RMax 2.7 
set PPEventHandler:Emitter:PMinusOrdering 0.76
set Frag8:FragmentationScheme none

erase AriadneCascade:Emitters[0]
set PPEventHandler:CoherenceRange 3.9
set PPEventHandler:LambdaQCD 0.26
set PPEventHandler:RMax 2.78
set PPEventHandler:Emitter:PMinusOrdering 1.1
set PPTestTune:EventHandler:LuminosityFunction:Energy 900
saverun PPT11XSECYN09 PPTestTune
set PPTestTune:EventHandler:LuminosityFunction:Energy 7000
saverun PPT11XSECYN71 PPTestTune
insert AriadneCascade:Emitters[0] FSSwinger
set PPGenerator:EventHandler:CascadeHandler NULL
saverun PPT11XSECYn71 PPTestTune
set PPGenerator:EventHandler:CascadeHandler AriadneCascade



set PPEventHandler:EventFiller:PTCut 0.6



# 12

set PPEventHandler:Swinger:Lambda 3.0
set PPEventHandler:CoherenceRange 1.0
set PPEventHandler:EffectivePartonMode Colours
do PPGenerator:RemoveInterface /DIPSY/PPEventHandler:CoherenceRange
set PPEventHandler:CoherenceRange 1.0
do PPGenerator:AddRndInterface /DIPSY/PPEventHandler:Emitter:PSInflation 1 0.5 2.0 1.0 0.2
erase PPGenerator:AnalysisHandlers[0]
erase PPGenerator:AnalysisHandlers[0]
set PPGenerator:EventHandler:PreSamples 1000
set PPGenerator:NumberOfEvents 0
SAVERUNYR PP12XSECYR PPGenerator
insert PPGenerator:AnalysisHandlers[0] Logger
insert PPGenerator:AnalysisHandlers[0] RivetTune
set PPGenerator:EventHandler:PreSamples 0
set PPGenerator:NumberOfEvents 10000


set PPGenerator:EventHandler:YFrametest -1.2
set PPGenerator:EventHandler:LuminosityFunction:Energy 900
saverun PP12XSECYR09 PPGenerator
set PPGenerator:EventHandler:LuminosityFunction:Energy 7000
saverun PP12XSECYR71 PPGenerator

set Frag8:FragmentationScheme pipe
set PPGenerator:EventHandler:LuminosityFunction:Energy 900
saverun PP12XSECYRR09 PPGenerator
set PPGenerator:EventHandler:LuminosityFunction:Energy 7000
saverun PP12XSECYRR71 PPGenerator
set Frag8:FragmentationScheme none

# 13

set PPEventHandler:EventFiller:PTCut 0.6
set PPEventHandler:Swinger:Lambda 3.0
set PPEventHandler:CoherenceRange 1.0
set PPEventHandler:EffectivePartonMode Colours
#do PPGenerator:RemoveInterface /DIPSY/PPEventHandler:CoherenceRange
set PPEventHandler:CoherenceRange 1.0
#do PPGenerator:AddRndInterface /DIPSY/PPEventHandler:Emitter:PSInflation 1 0.5 2.0 1.0 0.2
erase PPGenerator:AnalysisHandlers[0]
erase PPGenerator:AnalysisHandlers[0]
set PPGenerator:EventHandler:PreSamples 1000
set PPGenerator:NumberOfEvents 0
SAVERUNYR PP13XSECYR PPGenerator
insert PPGenerator:AnalysisHandlers[0] Logger
insert PPGenerator:AnalysisHandlers[0] RivetTune
set PPGenerator:EventHandler:PreSamples 0
set PPGenerator:NumberOfEvents 10000


set PPGenerator:EventHandler:YFrametest -1.2
set PPGenerator:EventHandler:LuminosityFunction:Energy 900
saverun PP13XSECYR09 PPGenerator
set PPGenerator:EventHandler:LuminosityFunction:Energy 7000
saverun PP13XSECYR71 PPGenerator

set Frag8:FragmentationScheme pipe
set PPGenerator:EventHandler:LuminosityFunction:Energy 900
saverun PP13XSECYRR09 PPGenerator
set PPGenerator:EventHandler:LuminosityFunction:Energy 7000
saverun PP13XSECYRR71 PPGenerator
set Frag8:FragmentationScheme none

# 14

set PPEventHandler:EventFiller:PTCut 1.0
set PPEventHandler:Swinger:Lambda 3.0
set PPEventHandler:CoherenceRange 1.0
set PPEventHandler:EffectivePartonMode Colours
#do PPGenerator:RemoveInterface /DIPSY/PPEventHandler:CoherenceRange
set PPEventHandler:CoherenceRange 1.0
#do PPGenerator:AddRndInterface /DIPSY/PPEventHandler:Emitter:PSInflation 1 0.5 2.0 1.0 0.2
erase PPGenerator:AnalysisHandlers[0]
erase PPGenerator:AnalysisHandlers[0]
set PPGenerator:EventHandler:PreSamples 1000
set PPGenerator:NumberOfEvents 0
SAVERUNYR PP14XSECYR PPGenerator
insert PPGenerator:AnalysisHandlers[0] Logger
insert PPGenerator:AnalysisHandlers[0] RivetTune
set PPGenerator:EventHandler:PreSamples 0
set PPGenerator:NumberOfEvents 10000


set PPGenerator:EventHandler:YFrametest -1.2
set PPGenerator:EventHandler:LuminosityFunction:Energy 900
saverun PP14XSECYR09 PPGenerator
set PPGenerator:EventHandler:LuminosityFunction:Energy 7000
saverun PP14XSECYR71 PPGenerator

set Frag8:FragmentationScheme pipe
set PPGenerator:EventHandler:LuminosityFunction:Energy 900
saverun PP14XSECYRR09 PPGenerator
set PPGenerator:EventHandler:LuminosityFunction:Energy 7000
saverun PP14XSECYRR71 PPGenerator
set Frag8:FragmentationScheme none

# 15

#set PPEventHandler:EventFiller:PTCut 1.0
set PPEventHandler:Emitter:PSInflation 1.0
set PPEventHandler:Swinger:Lambda 3.0
set PPEventHandler:CoherenceRange 1.0
set PPEventHandler:EffectivePartonMode Colours
do PPGenerator:RemoveInterface /DIPSY/PPEventHandler:Emitter:PSInflation
set PPEventHandler:CoherenceRange 1.0
do PPGenerator:AddRndInterface /DIPSY/PPEventHandler:EventFiller:PTCut 1 0.3 1.7 1.0 0.5
erase PPGenerator:AnalysisHandlers[0]
erase PPGenerator:AnalysisHandlers[0]
set PPGenerator:EventHandler:PreSamples 1000
set PPGenerator:NumberOfEvents 0
SAVERUNYR PP15XSECYR PPGenerator
insert PPGenerator:AnalysisHandlers[0] Logger
insert PPGenerator:AnalysisHandlers[0] RivetTune
set PPGenerator:EventHandler:PreSamples 0
set PPGenerator:NumberOfEvents 10000


set PPGenerator:EventHandler:YFrametest -1.2
set PPGenerator:EventHandler:LuminosityFunction:Energy 900
saverun PP15XSECYR09 PPGenerator
set PPGenerator:EventHandler:LuminosityFunction:Energy 7000
saverun PP15XSECYR71 PPGenerator

set Frag8:FragmentationScheme pipe
set PPGenerator:EventHandler:LuminosityFunction:Energy 900
saverun PP15XSECYRR09 PPGenerator
set PPGenerator:EventHandler:LuminosityFunction:Energy 7000
saverun PP15XSECYRR71 PPGenerator
set Frag8:FragmentationScheme none

# 16 inclusive observables

set PPGenerator:EventHandler:FudgeME true
set PPEventHandler:EventFiller:PTCut 1.0
set PPEventHandler:Emitter:PSInflation 1.0
set PPEventHandler:Swinger:Lambda 1.0
set PPEventHandler:CoherenceRange 2.5
set PPEventHandler:EffectivePartonMode Colours
do PPGenerator:RemoveInterface /DIPSY/PPEventHandler:EventFiller:PTCut
set PPEventHandler:CoherenceRange 2.5
do PPGenerator:AddRndInterface /DIPSY/PPEventHandler:Emitter:PTScale 1 0.5 2.0 1.4 0.5
erase PPGenerator:AnalysisHandlers[0]
erase PPGenerator:AnalysisHandlers[0]
set PPGenerator:EventHandler:PreSamples 1000
set PPGenerator:NumberOfEvents 0
SAVERUNYM PP16XSECYM PPGenerator
insert PPGenerator:AnalysisHandlers[0] Logger
insert PPGenerator:AnalysisHandlers[0] RivetTune
set PPGenerator:EventHandler:PreSamples 0
set PPGenerator:NumberOfEvents 10000


# 16 final state observables

set PPGenerator:EventHandler:YFrametest 0.5
SAVERUNFS PP16XSECYM PPGenerator

set Frag8:FragmentationScheme dep2
set PPGenerator:EventHandler:LuminosityFunction:Energy 200
SAVERUNFS PP16XSECYMQ PPGenerator
SAVERUNFS PP16XSECYMq PPGenerator
SAVERUNFS PP16XSECYMo PPGenerator

set Frag8:FragmentationScheme dipole
set Frag8:ThrowAway true
set PPGenerator:DumpPeriod 0
set Frag8:StringR0 0.8
set Frag8:Stringm0 0.2
set Frag8:BaryonSuppression 0.22 
SAVERUNFS PP16XSECYMF PPGenerator

set Frag8:BaryonSuppression 0.33
SAVERUNFS PP16XSECYMf PPGenerator

set Frag8:StringR0 1.2
set Frag8:Stringm0 0.2
set Frag8:BaryonSuppression 0.167
SAVERUNFS PP16XSECYMK PPGenerator

set Frag8:StringR0 1.2
set Frag8:Stringm0 0.2
set Frag8:BaryonSuppression 0.22 
SAVERUNFS PP16XSECYMJ PPGenerator
SAVERUNFS PP16XSECYMW PPGenerator

set Frag8:StringR0 1.0
set Frag8:Stringm0 0.2
set Frag8:BaryonSuppression 0.25 
SAVERUNFS PP16XSECYMV PPGenerator

set Frag8:StringR0 1.0
set Frag8:Stringm0 0.2
set Frag8:BaryonSuppression 0.25
set Frag8:Average false
SAVERUNFS PP16XSECYMX PPGenerator
SAVERUNFS PP16XSECYXX PPGenerator
set Frag8:BaryonSuppression 0.5
SAVERUNFS PP16XSECYXX PPGenerator
set Frag8:StringR0 0.8
set Frag8:Stringm0 0.2
set Frag8:BaryonSuppression 0.2
set Frag8:Average false
set PPGenerator:EventHandler:HadronizationHandler Frag8Rope
SAVERUNFS PP16XSECYXR PPGenerator
set PPGenerator:EventHandler:HadronizationHandler Frag8
SAVERUNFS PP16XSECYMx PPGenerator
set Frag8:Average true

set Frag8:StringR0 1.2
set Frag8:Stringm0 0.2
set Frag8:BaryonSuppression 0.5
SAVERUNFS PP16XSECYMj PPGenerator
SAVERUNFS PP16XSECYMw PPGenerator

set Frag8:StringR0 1.0
set Frag8:Stringm0 0.2
set Frag8:BaryonSuppression 0.5
SAVERUNFS PP16XSECYMv PPGenerator

set Frag8:FragmentationScheme average
set Frag8:ThrowAway true
set PPGenerator:DumpPeriod 0
set Frag8:StringR0 1.0
set Frag8:Stringm0 1.0
SAVERUNFS PP16XSECYMD PPGenerator

set Frag8:StringR0 0.5
set Frag8:Stringm0 0.5
SAVERUNFS PP16XSECYMd PPGenerator

set Frag8:Stringm0 1.0
set PPGenerator:EventHandler:LuminosityFunction:Energy 200
SAVERUNFS PP16XSECYMt PPGenerator

set Frag8:StringR0 1.2
set Frag8:Stringm0 0.6
set Frag8:BaryonSuppression 0.33
SAVERUNFS PP16XSECYMT PPGenerator

set Frag8:ThrowAway false
SAVERUNFS PP16XSECYMu PPGenerator

set Frag8:StringR0 0.5
set Frag8:Stringm0 1.0
set Frag8:BaryonSuppression 0.5

set Frag8:FragmentationScheme pipe
set Frag8:ThrowAway false
SAVERUNFS PP16XSECYMR PPGenerator
SAVERUNFS PP16XSECYMr PPGenerator
SAVERUNFS PP16XSECYOr PPGenerator

set FSSwinger:SizeOpt 2
set FSOrdering:HardSuppression 2
set PPGenerator:EventHandler:LuminosityFunction:Energy 200
SAVERUNFS PP16XSECYoR PPGenerator

set FSSwinger:MaxRho 1.0
set PPGenerator:EventHandler:HadronizationHandler Frag8Sw2
SAVERUNFS PP16XSECY2R PPGenerator

set PPGenerator:EventHandler:HadronizationHandler Frag8
set FSSwinger:MaxRho -2.0
erase AriadneCascade:Emitters[0]
SAVERUNFS PP16XSECYsR PPGenerator

set PPGenerator:EventHandler:HadronizationHandler Frag8Sw0
SAVERUNFS PP16XSECYMS PPGenerator
SAVERUNFS PP16XSECYOS PPGenerator

SAVERUNFS PP16XSECYoS PPGenerator


set PPGenerator:EventHandler:HadronizationHandler Frag8
insert AriadneCascade:Emitters[0] FSSwinger
set FSSwinger:SizeOpt 1
set FSOrdering:HardSuppression 0
set Frag8:FragmentationScheme none
set Frag8:ThrowAway false

# 16 testing the tunes

set PPEventHandler:LambdaQCD 0.225
set PPEventHandler:RMax 3.0
set PPEventHandler:Emitter:PMinusOrdering 0.67
set PPEventHandler:Emitter:PTScale 1.25
SAVERUNTUNE PPT16XSECYMp PPTestTune PP16XSECYM

set Frag8:FragmentationScheme pipe
set PPEventHandler:CoherenceRange 2.5
set PPEventHandler:LambdaQCD 0.25
set PPEventHandler:RMax 2.9
set PPEventHandler:Emitter:PMinusOrdering 0.69
set PPEventHandler:Emitter:PTScale 1.22
saverun PPT16XSECYMR71 PPTestTune
set PPTestTune:EventHandler:LuminosityFunction:Energy 200
saverun PPT16XSECYMR02 PPTestTune
set PPTestTune:EventHandler:LuminosityFunction:Energy 900
saverun PPT16XSECYMR09 PPTestTune
set PPTestTune:EventHandler:LuminosityFunction:Energy 7000
set PPEventHandler:LambdaQCD 0.25
set PPEventHandler:RMax 2.9
set PPEventHandler:Emitter:PMinusOrdering 0.71
set PPEventHandler:Emitter:PTScale 1.22
saverun PPT16XSECYMr71 PPTestTune
set PPTestTune:EventHandler:LuminosityFunction:Energy 200
saverun PPT16XSECYMr02 PPTestTune
set PPTestTune:EventHandler:LuminosityFunction:Energy 900
saverun PPT16XSECYMr09 PPTestTune
set PPTestTune:EventHandler:LuminosityFunction:Energy 7000

set Frag8:FragmentationScheme dep2
SAVERUNTUNE PPT16XSECYMo PPTestTune PP16XSECYMo

set Frag8:FragmentationScheme dipole
set Frag8:ThrowAway true
set PPGenerator:DumpPeriod 0
set Frag8:StringR0 0.8
set Frag8:Stringm0 0.2
set Frag8:BaryonSuppression 0.22 
SAVERUNTUNE PPT16XSECYMF PPTestTune PP16XSECYMF

set Frag8:FragmentationScheme dipole
set Frag8:ThrowAway true
set PPGenerator:DumpPeriod 0
set Frag8:StringR0 0.8
set Frag8:Stringm0 0.2
set Frag8:BaryonSuppression 0.33 
SAVERUNTUNE PPT16XSECYMf PPTestTune PP16XSECYMf
set PPTestTune:EventHandler:LuminosityFunction:Energy 900
set Frag8:FragmentationScheme dep1
set PPTestTune:NumberOfEvents 1000
saverun quicktest PPTestTune
set Frag8:FragmentationScheme dipole
set PPTestTune:EventHandler:LuminosityFunction:Energy 7000
set PPTestTune:NumberOfEvents 100000

set Frag8:StringR0 1.2
set Frag8:Stringm0 0.2
set Frag8:BaryonSuppression 0.167 
SAVERUNTUNE PPT16XSECYMk PPTestTune PP16XSECYMK

SAVERUNTUNE PPT16XSECYMK PPTestTune PP16XSECYMK

set Frag8:StringR0 1.2
set Frag8:Stringm0 0.2
set Frag8:BaryonSuppression 0.22 
SAVERUNTUNE PPT16XSECYMJ PPTestTune PP16XSECYMJ 
SAVERUNTUNE PPT16XSECYMI PPTestTune PP16XSECYMJ 
SAVERUNTUNE PPT16XSECYMi PPTestTune PP16XSECYMF 
SAVERUNTUNE PPT16XSECYMW PPTestTune PP16XSECYMW 

set Frag8:StringR0 1.0
set Frag8:Stringm0 0.2
set Frag8:BaryonSuppression 0.25 
SAVERUNTUNE PPT16XSECYMV PPTestTune PP16XSECYMV 

set Frag8:StringR0 1.0
set Frag8:Stringm0 0.2
set Frag8:BaryonSuppression 0.25
set Frag8:Average false
#set PPTestTune:LogNonDefault Yes
SAVERUNTUNE PPT16XSECYMX PPTestTune PP16XSECYMX 
SAVERUNTUNE PPT16XSECYMXX PPTestTune PP16XSECYXX
SAVERUNTUNE PPT16XSECYMXP PPTestTune PP16XSECYXX
set Frag8:StringFlav_popcornRate 0.8
SAVERUNTUNE PPT16XSECYMXx PPTestTune PP16XSECYXX 
set Frag8:StringFlav_popcornRate 0.5
set Frag8:BaryonSuppression 0.5
SAVERUNTUNE PPT16XSECYMXZ PPTestTune PP16XSECYXZ 
set Frag8:BaryonSuppression 0.25
set PPGenerator:EventHandler:HadronizationHandler Frag8Rope
SAVERUNTUNE PPT16XSECYMXR PPTestTune PP16XSECYXR
set PPGenerator:EventHandler:HadronizationHandler Frag8
set Frag8:StringR0 0.8
set Frag8:Stringm0 0.2
set Frag8:BaryonSuppression 0.2
set Frag8:Average false
#set PPTestTune:LogNonDefault Yes
SAVERUNTUNE PPT16XSECYMx PPTestTune PP16XSECYMx 
set Frag8:Average true

set Frag8:StringR0 1.2
set Frag8:Stringm0 0.2
set Frag8:BaryonSuppression 0.5
SAVERUNTUNE PPT16XSECYMj PPTestTune PP16XSECYMj 
SAVERUNTUNE PPT16XSECYMw PPTestTune PP16XSECYMw 

set Frag8:StringR0 1.0
set Frag8:Stringm0 0.2
set Frag8:BaryonSuppression 0.5
SAVERUNTUNE PPT16XSECYMv PPTestTune PP16XSECYMv

set Frag8:FragmentationScheme average
set Frag8:ThrowAway true
set PPGenerator:DumpPeriod 0
set Frag8:StringR0 1.0
set Frag8:Stringm0 1.0
set Frag8:BaryonSuppression 0.5
set PPEventHandler:LambdaQCD 0.245
set PPEventHandler:RMax 3.0
set PPEventHandler:Emitter:PMinusOrdering 0.69
set PPEventHandler:Emitter:PTScale 1.2
set PPTestTune:EventHandler:LuminosityFunction:Energy 200
saverun PPT16XSECYMD02 PPTestTune
set PPTestTune:EventHandler:LuminosityFunction:Energy 900
saverun PPT16XSECYMD09 PPTestTune
set PPTestTune:EventHandler:LuminosityFunction:Energy 7000
saverun PPT16XSECYMD71 PPTestTune
set PPEventHandler:LambdaQCD 0.245
set PPEventHandler:RMax 3.0
set PPEventHandler:Emitter:PMinusOrdering 0.88
set PPEventHandler:Emitter:PTScale 1.3
set PPTestTune:EventHandler:LuminosityFunction:Energy 200
saverun PPT16XSECYMD02z PPTestTune
set PPTestTune:EventHandler:LuminosityFunction:Energy 900
saverun PPT16XSECYMD09z PPTestTune
set PPTestTune:EventHandler:LuminosityFunction:Energy 7000
saverun PPT16XSECYMD71z PPTestTune

set PPEventHandler:LambdaQCD 0.21
set PPEventHandler:RMax 3.1
set PPEventHandler:Emitter:PMinusOrdering 0.75
set PPEventHandler:Emitter:PTScale 1.25
set PPTestTune:EventHandler:LuminosityFunction:Energy 200
saverun PPT16XSECYMD0202 PPTestTune
set PPEventHandler:LambdaQCD 0.23
set PPEventHandler:RMax 3.05
set PPEventHandler:Emitter:PMinusOrdering 0.78
set PPEventHandler:Emitter:PTScale 1.2
set PPTestTune:EventHandler:LuminosityFunction:Energy 200
saverun PPT16XSECYMD0202z PPTestTune

set PPEventHandler:LambdaQCD 0.255
set PPEventHandler:RMax 3.05
set PPEventHandler:Emitter:PMinusOrdering 0.85
set PPEventHandler:Emitter:PTScale 1.1
set PPTestTune:EventHandler:LuminosityFunction:Energy 900
saverun PPT16XSECYMD0909 PPTestTune
set PPEventHandler:LambdaQCD 0.23
set PPEventHandler:RMax 3.15
set PPEventHandler:Emitter:PMinusOrdering 0.97
set PPEventHandler:Emitter:PTScale 1.2
set PPTestTune:EventHandler:LuminosityFunction:Energy 900
saverun PPT16XSECYMD0909z PPTestTune

set PPEventHandler:LambdaQCD 0.245
set PPEventHandler:RMax 2.95
set PPEventHandler:Emitter:PMinusOrdering 0.66
set PPEventHandler:Emitter:PTScale 1.25
set PPTestTune:EventHandler:LuminosityFunction:Energy 7000
saverun PPT16XSECYMD7170 PPTestTune
set PPEventHandler:LambdaQCD 0.245
set PPEventHandler:RMax 3.0
set PPEventHandler:Emitter:PMinusOrdering 0.95
set PPEventHandler:Emitter:PTScale 1.4
set PPTestTune:EventHandler:LuminosityFunction:Energy 7000
saverun PPT16XSECYMD7170z PPTestTune


set Frag8:StringR0 0.5
set Frag8:Stringm0 0.5
set PPEventHandler:LambdaQCD 0.26
set PPEventHandler:RMax 2.9
set PPEventHandler:Emitter:PMinusOrdering 0.73
set PPEventHandler:Emitter:PTScale 1.2
set PPTestTune:EventHandler:LuminosityFunction:Energy 200
saverun PPT16XSECYMd02 PPTestTune
set PPTestTune:EventHandler:LuminosityFunction:Energy 900
saverun PPT16XSECYMd09 PPTestTune
set PPTestTune:EventHandler:LuminosityFunction:Energy 7000
saverun PPT16XSECYMd71 PPTestTune
set PPEventHandler:LambdaQCD 0.215
set PPEventHandler:RMax 3.05
set PPEventHandler:Emitter:PMinusOrdering 0.64
set PPEventHandler:Emitter:PTScale 1.3
set PPTestTune:EventHandler:LuminosityFunction:Energy 200
saverun PPT16XSECYMd02z PPTestTune
set PPTestTune:EventHandler:LuminosityFunction:Energy 900
saverun PPT16XSECYMd09z PPTestTune
set PPTestTune:EventHandler:LuminosityFunction:Energy 7000
saverun PPT16XSECYMd71z PPTestTune

set PPEventHandler:LambdaQCD 0.22
set PPEventHandler:RMax 3.1
set PPEventHandler:Emitter:PMinusOrdering 0.74
set PPEventHandler:Emitter:PTScale 1.2
set PPTestTune:EventHandler:LuminosityFunction:Energy 200
saverun PPT16XSECYMd0202 PPTestTune
set PPEventHandler:LambdaQCD 0.22
set PPEventHandler:RMax 3.05
set PPEventHandler:Emitter:PMinusOrdering 0.795
set PPEventHandler:Emitter:PTScale 1.25
set PPTestTune:EventHandler:LuminosityFunction:Energy 200
saverun PPT16XSECYMd0202z PPTestTune

set PPEventHandler:LambdaQCD 0.205
set PPEventHandler:RMax 3.25
set PPEventHandler:Emitter:PMinusOrdering 0.89
set PPEventHandler:Emitter:PTScale 1.2
set PPTestTune:EventHandler:LuminosityFunction:Energy 900
saverun PPT16XSECYMd0909 PPTestTune
set PPEventHandler:LambdaQCD 0.19
set PPEventHandler:RMax 3.25
set PPEventHandler:Emitter:PMinusOrdering 0.795
set PPEventHandler:Emitter:PTScale 1.25
set PPTestTune:EventHandler:LuminosityFunction:Energy 900
saverun PPT16XSECYMd0909z PPTestTune

set PPEventHandler:LambdaQCD 0.255
set PPEventHandler:RMax 2.9
set PPEventHandler:Emitter:PMinusOrdering 0.69
set PPEventHandler:Emitter:PTScale 1.26
set PPTestTune:EventHandler:LuminosityFunction:Energy 7000
saverun PPT16XSECYMd7170 PPTestTune
set PPEventHandler:LambdaQCD 0.25
set PPEventHandler:RMax 2.9
set PPEventHandler:Emitter:PMinusOrdering 0.65
set PPEventHandler:Emitter:PTScale 1.3
set PPTestTune:EventHandler:LuminosityFunction:Energy 7000
saverun PPT16XSECYMd7170z PPTestTune

set Frag8:StringR0 1.2
set Frag8:Stringm0 0.6
set Frag8:BaryonSuppression 0.33
set PPEventHandler:LambdaQCD 0.25
set PPEventHandler:RMax 2.95
set PPEventHandler:Emitter:PMinusOrdering 0.72
set PPEventHandler:Emitter:PTScale 1.2
set PPTestTune:EventHandler:LuminosityFunction:Energy 200
saverun PPT16XSECYMT02 PPTestTune
set PPTestTune:EventHandler:LuminosityFunction:Energy 900
saverun PPT16XSECYMT09 PPTestTune
set PPTestTune:EventHandler:LuminosityFunction:Energy 7000
saverun PPT16XSECYMT71 PPTestTune
set PPEventHandler:LambdaQCD 0.255
set PPEventHandler:RMax 2.9
set PPEventHandler:Emitter:PMinusOrdering 0.67
set PPEventHandler:Emitter:PTScale 1.2
set PPTestTune:EventHandler:LuminosityFunction:Energy 200
saverun PPT16XSECYMT02z PPTestTune
set PPTestTune:EventHandler:LuminosityFunction:Energy 900
saverun PPT16XSECYMT09z PPTestTune
set PPTestTune:EventHandler:LuminosityFunction:Energy 7000
saverun PPT16XSECYMT71z PPTestTune

set PPEventHandler:LambdaQCD 0.225
set PPEventHandler:RMax 3.05
set PPEventHandler:Emitter:PMinusOrdering 0.78
set PPEventHandler:Emitter:PTScale 1.2
set PPTestTune:EventHandler:LuminosityFunction:Energy 200
saverun PPT16XSECYMT0202 PPTestTune
set PPEventHandler:LambdaQCD 0.215
set PPEventHandler:RMax 3.15
set PPEventHandler:Emitter:PMinusOrdering 0.805
set PPEventHandler:Emitter:PTScale 1.2
set PPTestTune:EventHandler:LuminosityFunction:Energy 200
saverun PPT16XSECYMT0202z PPTestTune

set PPEventHandler:LambdaQCD 0.245
set PPEventHandler:RMax 3.1
set PPEventHandler:Emitter:PMinusOrdering 0.78
set PPEventHandler:Emitter:PTScale 1.15
set PPTestTune:EventHandler:LuminosityFunction:Energy 900
saverun PPT16XSECYMT0909 PPTestTune
set PPEventHandler:LambdaQCD 0.22
set PPEventHandler:RMax 3.15
set PPEventHandler:Emitter:PMinusOrdering 0.73
set PPEventHandler:Emitter:PTScale 1.2
set PPTestTune:EventHandler:LuminosityFunction:Energy 900
saverun PPT16XSECYMT0909z PPTestTune

set PPEventHandler:LambdaQCD 0.25
set PPEventHandler:RMax 2.95
set PPEventHandler:Emitter:PMinusOrdering 0.68
set PPEventHandler:Emitter:PTScale 1.3
set PPTestTune:EventHandler:LuminosityFunction:Energy 7000
saverun PPT16XSECYMT7170 PPTestTune
set PPEventHandler:LambdaQCD 0.265
set PPEventHandler:RMax 2.85
set PPEventHandler:Emitter:PMinusOrdering 0.685
set PPEventHandler:Emitter:PTScale 1.25
set PPTestTune:EventHandler:LuminosityFunction:Energy 7000
saverun PPT16XSECYMT7170z PPTestTune

set Frag8:ThrowAway false
set PPEventHandler:LambdaQCD 0.23
set PPEventHandler:RMax 3.0
set PPEventHandler:Emitter:PMinusOrdering 0.705
set PPEventHandler:Emitter:PTScale 1.25
set PPTestTune:EventHandler:LuminosityFunction:Energy 200
saverun PPT16XSECYMu02 PPTestTune
set PPTestTune:EventHandler:LuminosityFunction:Energy 900
saverun PPT16XSECYMu09 PPTestTune
set PPTestTune:EventHandler:LuminosityFunction:Energy 7000
saverun PPT16XSECYMu71 PPTestTune

set PPEventHandler:LambdaQCD 0.215
set PPEventHandler:RMax 3.15
set PPEventHandler:Emitter:PMinusOrdering 0.82
set PPEventHandler:Emitter:PTScale 1.2
set PPTestTune:EventHandler:LuminosityFunction:Energy 200
saverun PPT16XSECYMu0202 PPTestTune

set PPEventHandler:LambdaQCD 0.195
set PPEventHandler:RMax 3.3
set PPEventHandler:Emitter:PMinusOrdering 0.83
set PPEventHandler:Emitter:PTScale 1.2
set PPTestTune:EventHandler:LuminosityFunction:Energy 900
saverun PPT16XSECYMu0909 PPTestTune

set PPEventHandler:LambdaQCD 0.245
set PPEventHandler:RMax 2.9
set PPEventHandler:Emitter:PMinusOrdering 0.675
set PPEventHandler:Emitter:PTScale 1.3
set PPTestTune:EventHandler:LuminosityFunction:Energy 7000
saverun PPT16XSECYMu7170 PPTestTune

# set PPTestTune:NumberOfEvents 10000
set Frag8:StringyCut 0.0
saverun QuickRun0 PPTestTune
set Frag8:StringyCut 1.0
saverun QuickRun1 PPTestTune
set Frag8:StringyCut 0.0
set Frag8:BaryonSuppression 1.0
saverun QuickRun2 PPTestTune
set Frag8:BaryonSuppression 0.5
saverun QuickRun3 PPTestTune
set PPTestTune:NumberOfEvents 100000

set Frag8:StringR0 0.5
set Frag8:Stringm0 1.0
set PPEventHandler:LambdaQCD 0.23
set PPEventHandler:RMax 3.0
set PPEventHandler:Emitter:PMinusOrdering 0.67
set PPEventHandler:Emitter:PTScale 1.5
set PPTestTune:EventHandler:LuminosityFunction:Energy 200
saverun PPT16XSECYMt02 PPTestTune
set PPTestTune:EventHandler:LuminosityFunction:Energy 900
saverun PPT16XSECYMt09 PPTestTune
set PPTestTune:EventHandler:LuminosityFunction:Energy 7000
saverun PPT16XSECYMt71 PPTestTune

set PPEventHandler:LambdaQCD 0.22
set PPEventHandler:RMax 3.1
set PPEventHandler:Emitter:PMinusOrdering 0.80
set PPEventHandler:Emitter:PTScale 1.2
set PPTestTune:EventHandler:LuminosityFunction:Energy 200
saverun PPT16XSECYMt0202 PPTestTune

set PPEventHandler:LambdaQCD 0.21
set PPEventHandler:RMax 3.2
set PPEventHandler:Emitter:PMinusOrdering 0.82
set PPEventHandler:Emitter:PTScale 1.2
set PPTestTune:EventHandler:LuminosityFunction:Energy 900
saverun PPT16XSECYMt0909 PPTestTune

set PPEventHandler:LambdaQCD 0.235
set PPEventHandler:RMax 2.95
set PPEventHandler:Emitter:PMinusOrdering 0.63
set PPEventHandler:Emitter:PTScale 1.3
set PPTestTune:EventHandler:LuminosityFunction:Energy 7000
saverun PPT16XSECYMt7170 PPTestTune



set Frag8:FragmentationScheme pipe
set Frag8:ThrowAway false
set Frag8:StringR0 0.5

set PPEventHandler:LambdaQCD 0.24
set PPEventHandler:RMax 3.0
set PPEventHandler:Emitter:PMinusOrdering 0.68
set PPEventHandler:Emitter:PTScale 1.24
saverun PPT16XSECYMR7102 PPTestTune
set PPTestTune:EventHandler:LuminosityFunction:Energy 900
saverun PPT16XSECYMR0902 PPTestTune
set PPTestTune:EventHandler:LuminosityFunction:Energy 200
saverun PPT16XSECYMR0202 PPTestTune
set PPTestTune:EventHandler:LuminosityFunction:Energy 7000
set PPEventHandler:LambdaQCD 0.22
set PPEventHandler:RMax 3.1
set PPEventHandler:Emitter:PMinusOrdering 0.78
set PPEventHandler:Emitter:PTScale 1.21
saverun PPT16XSECYMr7102 PPTestTune
set PPTestTune:EventHandler:LuminosityFunction:Energy 900
saverun PPT16XSECYMr0902 PPTestTune
set PPTestTune:EventHandler:LuminosityFunction:Energy 200
saverun PPT16XSECYMr0202 PPTestTune
set PPTestTune:EventHandler:LuminosityFunction:Energy 7000

set PPEventHandler:LambdaQCD 0.19
set PPEventHandler:RMax 3.3
set PPEventHandler:Emitter:PMinusOrdering 0.87
set PPEventHandler:Emitter:PTScale 1.22
saverun PPT16XSECYMR7109 PPTestTune
set PPTestTune:EventHandler:LuminosityFunction:Energy 900
saverun PPT16XSECYMR0909 PPTestTune
set PPTestTune:EventHandler:LuminosityFunction:Energy 7000
set PPEventHandler:LambdaQCD 0.17
set PPEventHandler:RMax 3.4
set PPEventHandler:Emitter:PMinusOrdering 0.93
set PPEventHandler:Emitter:PTScale 1.30
saverun PPT16XSECYMr7109 PPTestTune
set PPTestTune:EventHandler:LuminosityFunction:Energy 900
saverun PPT16XSECYMr0909 PPTestTune
set PPTestTune:EventHandler:LuminosityFunction:Energy 7000

set PPEventHandler:LambdaQCD 0.25
set PPEventHandler:RMax 2.9
set PPEventHandler:Emitter:PMinusOrdering 0.67
set PPEventHandler:Emitter:PTScale 1.28
saverun PPT16XSECYMR7170 PPTestTune
set PPTestTune:EventHandler:LuminosityFunction:Energy 900
saverun PPT16XSECYMR0970 PPTestTune
set PPTestTune:EventHandler:LuminosityFunction:Energy 7000
set PPEventHandler:LambdaQCD 0.24
set PPEventHandler:RMax 2.9
set PPEventHandler:Emitter:PMinusOrdering 0.65
set PPEventHandler:Emitter:PTScale 1.28
saverun PPT16XSECYMr7170 PPTestTune
set PPTestTune:EventHandler:LuminosityFunction:Energy 900
saverun PPT16XSECYMr0970 PPTestTune
set PPTestTune:EventHandler:LuminosityFunction:Energy 7000

set PPEventHandler:LambdaQCD 0.215
set PPEventHandler:RMax 3.0
set PPEventHandler:Emitter:PMinusOrdering 0.64
set PPEventHandler:Emitter:PTScale 1.27
set PPTestTune:EventHandler:LuminosityFunction:Energy 200
saverun PPT16XSECYOr02 PPTestTune
set PPTestTune:EventHandler:LuminosityFunction:Energy 900
saverun PPT16XSECYOr09 PPTestTune
set PPTestTune:EventHandler:LuminosityFunction:Energy 7000
saverun PPT16XSECYOr71 PPTestTune
set PPEventHandler:LambdaQCD 0.22
set PPEventHandler:RMax 3.15
set PPEventHandler:Emitter:PMinusOrdering 0.8
set PPEventHandler:Emitter:PTScale 1.2
set PPTestTune:EventHandler:LuminosityFunction:Energy 200
saverun PPT16XSECYOr0202 PPTestTune
set PPEventHandler:LambdaQCD 0.18
set PPEventHandler:RMax 3.35
set PPEventHandler:Emitter:PMinusOrdering 0.91
set PPEventHandler:Emitter:PTScale 1.25
set PPTestTune:EventHandler:LuminosityFunction:Energy 900
saverun PPT16XSECYOr0909 PPTestTune
set PPEventHandler:LambdaQCD 0.23
set PPEventHandler:RMax 2.95
set PPEventHandler:Emitter:PMinusOrdering 0.65
set PPEventHandler:Emitter:PTScale 1.3
set PPTestTune:EventHandler:LuminosityFunction:Energy 7000
saverun PPT16XSECYOr7170 PPTestTune



set PPEventHandler:LambdaQCD 0.23
set PPEventHandler:RMax 3.3
set PPEventHandler:Emitter:PMinusOrdering 0.65
set PPEventHandler:Emitter:PTScale 1.25
saverun PPT16XSECYMR71FS PPTestTune
set PPTestTune:EventHandler:LuminosityFunction:Energy 900
saverun PPT16XSECYMR09FS PPTestTune
set PPTestTune:EventHandler:LuminosityFunction:Energy 7000
set FSSwinger:SizeOpt 2
set FSOrdering:HardSuppression 2
set PPEventHandler:LambdaQCD 0.22
set PPEventHandler:RMax 3.0
set PPEventHandler:Emitter:PMinusOrdering 0.63
set PPEventHandler:Emitter:PTScale 1.27
saverun PPT16XSECYoR71 PPTestTune
set FSSwinger:MaxRho 1.0
set PPEventHandler:LambdaQCD 0.20
set PPEventHandler:RMax 3.1
set PPEventHandler:Emitter:PMinusOrdering 0.67
set PPEventHandler:Emitter:PTScale 1.31
set PPGenerator:EventHandler:HadronizationHandler Frag8Sw2
saverun PPT16XSECY2R71 PPTestTune
set FSSwinger:MaxRho -2.0
set PPGenerator:EventHandler:HadronizationHandler Frag8
erase AriadneCascade:Emitters[0]
set PPEventHandler:LambdaQCD 0.17
set PPEventHandler:RMax 3.3
set PPEventHandler:Emitter:PMinusOrdering 1.0
set PPEventHandler:Emitter:PTScale 1.42
saverun PPT16XSECYsR71 PPTestTune
set PPGenerator:EventHandler:HadronizationHandler Frag8Sw0
set PPEventHandler:LambdaQCD 0.16
set PPEventHandler:RMax 3.3
set PPEventHandler:Emitter:PMinusOrdering 1.0
set PPEventHandler:Emitter:PTScale 1.44
saverun PPT16XSECYSR71 PPTestTune

SAVERUNTUNE PPT16XSECYoS PPTestTune PP16XSECYoS

set PPEventHandler:LambdaQCD 0.16
set PPEventHandler:RMax 3.3
set PPEventHandler:Emitter:PMinusOrdering 1.05
set PPEventHandler:Emitter:PTScale 1.45
set PPTestTune:EventHandler:LuminosityFunction:Energy 200
saverun PPT16XSECYOS02 PPTestTune
set PPTestTune:EventHandler:LuminosityFunction:Energy 900
saverun PPT16XSECYOS09 PPTestTune
set PPTestTune:EventHandler:LuminosityFunction:Energy 7000
saverun PPT16XSECYOS71 PPTestTune
set PPEventHandler:LambdaQCD 0.23
set PPEventHandler:RMax 3.1
set PPEventHandler:Emitter:PMinusOrdering 0.8
set PPEventHandler:Emitter:PTScale 1.15
set PPTestTune:EventHandler:LuminosityFunction:Energy 200
saverun PPT16XSECYOS0202 PPTestTune
set PPEventHandler:LambdaQCD 0.17
set PPEventHandler:RMax 3.5
set PPEventHandler:Emitter:PMinusOrdering 1.15
set PPEventHandler:Emitter:PTScale 1.25
set PPTestTune:EventHandler:LuminosityFunction:Energy 900
saverun PPT16XSECYOS0909 PPTestTune
set PPEventHandler:LambdaQCD 0.17
set PPEventHandler:RMax 3.2
set PPEventHandler:Emitter:PMinusOrdering 1.25
set PPEventHandler:Emitter:PTScale 1.6
set PPTestTune:EventHandler:LuminosityFunction:Energy 7000
saverun PPT16XSECYOS7170 PPTestTune




set Frag8:FragmentationScheme none
saverun PPT16XSECYOD7170 PPTestTune

#

set PPGenerator:EventHandler:YFrametest 0.5
erase PPGenerator:AnalysisHandlers[0]
erase PPGenerator:AnalysisHandlers[0]
set PPGenerator:EventHandler:PreSamples 1000
set PPGenerator:NumberOfEvents 0
set PPEventHandler:WFL stdProton
set PPEventHandler:WFR stdProton




do PPGenerator:RemoveInterface /DIPSY/PPEventHandler:Swinger:Lambda
do PPGenerator:RemoveInterface /DIPSY/PPEventHandler:Emitter:PMinusOrdering
do PPGenerator:RemoveInterface /DIPSY/PPEventHandler:RMax
do PPGenerator:RemoveInterface /DIPSY/PPEventHandler:LambdaQCD
do PPGenerator:RemoveInterface /DIPSY/PPEventHandler:BaryonSize
set PPEventHandler:RMax 2.95
set PPEventHandler:Swinger:Lambda 1.3
set PPEventHandler:Emitter:PMinusOrdering 0.52
set PPEventHandler:LambdaQCD 0.20
SAVERUN PP05XSECt PPGenerator

set PPEventHandler:RMax 2.24
set PPEventHandler:BaryonSize 3.45
set PPEventHandler:Emitter:PMinusOrdering 0.55
set PPEventHandler:LambdaQCD 0.25
set PPEventHandler:Swinger:Lambda 2.0
SAVERUN PP09XSECt PPGenerator

set PPEventHandler:RMax 2.04
set PPEventHandler:BaryonSize 3.64
set PPEventHandler:Emitter:PMinusOrdering 0.69
set PPEventHandler:LambdaQCD 0.28
set PPEventHandler:Swinger:Lambda 2.0
SAVERUN PP09XSECT PPGenerator


set PPEventHandler:RMax 3.0
set PPEventHandler:Swinger:Lambda 2.0
set PPEventHandler:Emitter:PMinusOrdering 0.48
set PPEventHandler:LambdaQCD 0.20
SAVERUN PP03XSECt PPGenerator

set Frag8:StringZ_aLund 0.42
set Frag8:StringZ_bLund 0.40
set Frag8:StringPT_sigma 0.32
set Frag8:StringFlav_probQQtoQ 0.084
set Frag8:StringFlav_probStoUD 0.22

set PPGenerator:EventHandler:YFrametest 0.5
set PPGenerator:EventHandler:DecayHandler /Defaults/Handlers/StandardDecayHandler
set PPGenerator:EventHandler:HadronizationHandler Frag8
set PPGenerator:EventHandler:CascadeHandler AriadneCascade
set PPGenerator:EventHandler:PreSamples 0
set PPGenerator:NumberOfEvents 100000
create ThePEG::ProgressLog Logger ProgressLog.so
set Logger:Interval 600
insert PPGenerator:AnalysisHandlers[0] Logger
create ThePEG::RivetAnalysis RivetLHC RivetAnalysis.so
insert RivetLHC:Analyses[0] ATLAS_2010_S8918562
insert PPGenerator:AnalysisHandlers[0] RivetLHC

saverun PP01Full70 PPGenerator
set PPEventHandler:WFL:RemnantDecayer /Defaults/Partons/BaryonRemnantDecayer
set PPEventHandler:EventFiller:SoftRemove NoValence
saverun PP03Full70 PPGenerator

set PPEventHandler:WFL Proton
set PPEventHandler:WFR Proton
set PPEventHandler:EventFiller:SoftRemove NoValence
# set PPEventHandler:EventFiller:ValenceChargeNormalisation None
insert RivetLHC:Paths[0] .
insert RivetLHC:Analyses[0] MC_MYFB
set PPGenerator:NumberOfEvents 1000000
saverun PP02Full70 PPGenerator


