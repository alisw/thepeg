# -*- eval: (ThePEG-repository-mode) -*-
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
set PPEventHandler:ConsistencyLevel 1
set PPEventHandler:XSecFn:CheckOffShell false
set PPEventHandler:CascadeHandler AriadneCascade
set PPEventHandler:HadronizationHandler Frag8
set PPEventHandler:DecayHandler /Defaults/Handlers/StandardDecayHandler
create ThePEG::FixedCMSLuminosity PPLumi
set PPEventHandler:LuminosityFunction PPLumi
cp Generator PPGenerator
erase PPGenerator:AnalysisHandlers[0]
set PPGenerator:HistogramFactory NULL
set PPGenerator:EventHandler PPEventHandler
set PPEventHandler:BGen:Width 5
set PPEventHandler:WFL Proton
set PPEventHandler:WFR Proton
set PPEventHandler:EventFiller:SoftRemove NoValence
set PPEventHandler:FudgeME 1
set PPEventHandler:EventFiller:PTCut 1.0
set PPEventHandler:Emitter:PSInflation 1.0
set PPEventHandler:Swinger:Lambda 1.0
set PPEventHandler:EffectivePartonMode Colours
set PPEventHandler:CoherenceRange 2.5

## These are the analysess we will run
## Some semi-inclusive cross section for DIPSY which need at least
## four combinations of left- and right-moving cascades.
erase PPEventHandler:AnalysisHandlers[0]
erase PPEventHandler:AnalysisHandlers[0]
create DIPSY::SemiInclusiveXSecAnalysis SemiIncl SemiInclusiveXSecAnalysis.so
insert PPEventHandler:AnalysisHandlers[0] SemiIncl
set PPEventHandler:PreSampleL 2
set PPEventHandler:PreSampleR 2
## This is just to keep track of the progress of a run
create DIPSY::AnalysisProgress AnaLog AnalysisProgress.so
set AnaLog:Interval 600
insert PPEventHandler:AnalysisHandlers[0] AnaLog

## The sample rates need to be adjusted so that we get a reasonable
## statistics in a reasonable time. It is typically efficient to
## sample a number of impact parameter values for each pair of DIPSY
## cascades.
set PPEventHandler:PreSampleB 1

## We need the same set of parameters for all different energies, so
## we use a separate random generator.
cp  /Defaults/Random RandomArg
set PPGenerator:SeparateRandom RandomArg

## These are the parameters we want to tune

## Now we want to run for all energies and two values for YFrametest

set Frag8:Collapser /Ariadne5/Defaults/Collapser

# This is the fragmentation parameters tuned without rope and with swing at LEP
# ../../src/FTuneSwLEP01.prin
set Frag8:FragmentationScheme none
set Frag8:StringZ_aLund 0.42
set Frag8:StringZ_bLund 0.40
set Frag8:StringPT_sigma 0.32
set Frag8:StringFlav_probQQtoQ 0.084
set Frag8:StringFlav_probStoUD 0.22

# This is the fragmentation parameters tuned without rope and swing at LEP
# ../../src/FTuneLEP01.prin
cp Frag8 Frag8Sw0
set Frag8Sw0:StringZ_aLund 0.30
set Frag8Sw0:StringZ_bLund 0.36
set Frag8Sw0:StringPT_sigma 0.32
set Frag8Sw0:StringFlav_probQQtoQ 0.082
set Frag8Sw0:StringFlav_probStoUD 0.22

# TEST fragmentation parameters tuned at LEP without rope and swing with MaxRho = 1
# ../../src/FTuneswLEP03.prin
cp Frag8 Frag8Sw2
set Frag8Sw2:StringZ_aLund 0.45
set Frag8Sw2:StringZ_bLund 0.42
set Frag8Sw2:StringPT_sigma 0.32
set Frag8Sw2:StringFlav_probQQtoQ 0.085
set Frag8Sw2:StringFlav_probStoUD 0.22

# This is the fragmentation parameters tuned to rope and default swing at LEP
#  ../../src/FTuneSwLEP07.prin
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

cp AriadneCascade AriadneNoSwing
erase AriadneNoSwing:Emitters[0]
set PPGenerator:EventHandler:PreSamples 0
set PPGenerator:EventHandler:YFrametest 0.5
set /DIPSY/PPEventHandler:BaryonSize 0.0
set FSSwinger:SetRmax -2.7
set FSSwinger:MaxRho -2.0




create ThePEG::ProgressLog Logger ProgressLog.so
set Logger:Interval 600
create ThePEG::RivetAnalysis RivetTune RivetAnalysis.so
insert RivetTune:Paths[0] .
insert RivetTune:Analyses[0] ATLAS_2010_S8918562
insert RivetTune:Analyses[0] STAR_2008_S7869363
insert RivetTune:Analyses[0] CMS_2011_S8978280
insert RivetTune:Analyses[0] PYTHIA_TUNING

insert PPGenerator:AnalysisHandlers[0] Logger
insert PPGenerator:AnalysisHandlers[0] RivetTune





set PPGenerator:EventHandler:LuminosityFunction:Energy 7000
set PPGenerator:DumpPeriod 0
cp PPGenerator PPTestTune
set PPGenerator:NumberOfEvents 10000

set PPTestTune:NumberOfEvents 100000
create ThePEG::LWHFactory LWHFactory LWHFactory.so
set LWHFactory:Suffix dat
set LWHFactory:StoreType flat
set PPTestTune:HistogramFactory LWHFactory
create DIPSY::PTAnalysis PTAnalysis PTAnalysis.so
insert PPTestTune:AnalysisHandlers[0] PTAnalysis

do PPGenerator:AddRndInterface /DIPSY/PPEventHandler:RMax 100 1.0 5.0 2.9 0.6
do PPGenerator:AddRndInterface /DIPSY/PPEventHandler:LambdaQCD 1 0.10 0.3 0.22 0.04
do PPGenerator:AddRndInterface /DIPSY/PPEventHandler:Emitter:PTScale 1 0.5 2.0 1.4 0.5
do PPGenerator:AddRndInterface /DIPSY/PPEventHandler:Emitter:PMinusOrdering 1 0.5 2.0 1.0 0.2

cp PPGenerator PPXSecGenerator
erase PPXSecGenerator:AnalysisHandlers[0]
erase PPXSecGenerator:AnalysisHandlers[0]
create ThePEG::Settings OnlyXSec
do OnlyXSec:set PPXSecGenerator:EventHandler:PreSamples 1000
insert PPXSecGenerator:DefaultObjects[0] OnlyXSec
set PPXSecGenerator:NumberOfEvents 0
cp OnlyXSec AddXSec
insert PPTestTune:DefaultObjects[0] AddXSec


create ThePEG::Settings Default
cp Default NoSwing
cp Default DefRope

do Default:set PPEventHandler:HadronizationHandler Frag8
do Default:set PPEventHandler:CascadeHandler AriadneCascade
do NoSwing:set PPEventHandler:HadronizationHandler Frag8Sw0
do NoSwing:set PPEventHandler:CascadeHandler AriadneNoSwing
do DefRope:set PPEventHandler:HadronizationHandler Frag8Rope
do DefRope:set PPEventHandler:CascadeHandler AriadneCascade

insert PPGenerator:DefaultObjects[0] Default
insert PPXSecGenerator:DefaultObjects[0] Default
insert PPTestTune:DefaultObjects[0] Default

cp Default NewSing
do NewSing:set PPEventHandler:EventFiller:Mode NewSingle
do NewSing:set PPTestTune:EventHandler:EventFiller:OnlyOnce true
do NewSing:set PPEventHandler:FudgeME 2
do NewSing:set PPEventHandler:XSecFn:SinFunction Average

cp Default NewSinF
do NewSinF:set PPEventHandler:XSecFn:SinFunction Average

cp Default EffPRel
do EffPRel:set PPEventHandler:EffectivePartonMode Relatives

cp Default DefaY7F
do DefaY7F:set PPGenerator:EventHandler:YFrametest 0.7

cp DefaY7F Def01Y7
do Def01Y7:set PPEventHandler:Swinger:Lambda 0.1

cp Default DefEven
do DefEven:set PPEventHandler:EventFiller:EffectiveWeights PlusEvenWeighted

cp Default DefSwL4
do DefSwL4:set PPEventHandler:Swinger:Lambda 4.0

cp Default DefSw01
do DefSw01:set PPEventHandler:Swinger:Lambda 0.1

cp Default DefFSGO
do DefFSGO:set FSOrdering:Generous Generous
do DefFSGO:set FSOrdering:OnlyOriginal OnlyOriginalEmissions

cp DefFSGO DefFSG2
do DefFSG2:set FSOrdering:Fudge 1.6

cp Default Shadows
do Shadows:set PPEventHandler:EffectivePartonMode Shadows

cp Default DefSing
do DefSing:set PPEventHandler:FudgeME 2
do DefSing:set PPTestTune:EventHandler:EventFiller:OnlyOnce true
do DefSing:set PPTestTune:EventHandler:EventFiller:SingleMother 1

cp DefSing DefComp
do DefComp:set PPTestTune:EventHandler:EventFiller:CompatMode 1

cp Default Shadow8
do Shadow8:set PPEventHandler:EffectivePartonMode Shadows
do Shadow8:set PPEventHandler:FudgeME 2
do Shadow8:set PPTestTune:EventHandler:EventFiller:OnlyOnce true
do Shadow8:set PPTestTune:EventHandler:EventFiller:SingleMother 1
do Shadow8:set stdEmitter:MinusOrderingMode OrderedShadow
do Shadow8:set stdXSec:IntOrdering ShadowOpen

cp Default Shadow9
do Shadow9:set PPEventHandler:EffectivePartonMode Shadows
do Shadow9:set PPEventHandler:FudgeME 2
do Shadow9:set PPTestTune:EventHandler:EventFiller:OnlyOnce true
do Shadow9:set PPTestTune:EventHandler:EventFiller:SingleMother 1
do Shadow9:set stdEmitter:MinusOrderingMode OrderedShadow
do Shadow9:set stdXSec:IntOrdering ShadowColour

# Run for tuning

SAVERUNFSXD Default PPGenerator Default PPXSecGenerator
SAVERUNFSXD DefRope PPGenerator Default PPXSecGenerator
SAVERUNFSXD NoSwing PPGenerator Default PPXSecGenerator
SAVERUNFSXD NewSinF PPGenerator NewSinF PPXSecGenerator
SAVERUNFSXD NewSing PPGenerator NewSing PPXSecGenerator
SAVERUNFSXD EffPRel PPGenerator EffPRel PPXSecGenerator
SAVERUNFSXD DefaY7F PPGenerator DefaY7F PPXSecGenerator
SAVERUNFSXD Def01Y7 PPGenerator Def01Y7 PPXSecGenerator
SAVERUNFSXD DefEven PPGenerator DefEven PPXSecGenerator
SAVERUNFSXD DefSwL4 PPGenerator DefSwL4 PPXSecGenerator
SAVERUNFSXD DefSw01 PPGenerator DefSw01 PPXSecGenerator
SAVERUNFSXD DefFSGO PPGenerator Default PPXSecGenerator
SAVERUNFSXD DefFSG2 PPGenerator Default PPXSecGenerator
SAVERUNFSXD Shadows PPGenerator Shadows PPXSecGenerator
SAVERUNFSXD DefSing PPGenerator DefSing PPXSecGenerator
SAVERUNFSXD DefComp PPGenerator DefComp PPXSecGenerator
SAVERUNFSXD Shadow8 PPGenerator Shadow8 PPXSecGenerator
SAVERUNFSXD Shadow9 PPGenerator Shadow9 PPXSecGenerator

# testing the tunes

SAVERUNTUNED Default PPTestTune
SAVERUNTUNED DefRope PPTestTune
SAVERUNTUNED NewSinF PPTestTune
SAVERUNTUNED NoSwing PPTestTune
SAVERUNTUNED NewSing PPTestTune
SAVERUNTUNED EffPRel PPTestTune
SAVERUNTUNED DefaY7F PPTestTune
SAVERUNTUNED Def01Y7 PPTestTune
SAVERUNTUNED DefEven PPTestTune
SAVERUNTUNED DefSwL4 PPTestTune
SAVERUNTUNED DefSw01 PPTestTune
SAVERUNTUNED DefFSGO PPTestTune
SAVERUNTUNED DefFSG2 PPTestTune
SAVERUNTUNED Shadows PPTestTune
SAVERUNTUNED DefSing PPTestTune
SAVERUNTUNED DefComp PPTestTune
SAVERUNTUNED Shadow8 PPTestTune
SAVERUNTUNED Shadow9 PPTestTune

##### playground #####
erase PPTestTune:DefaultObjects[0]
set PPTestTune:DefaultObjects[0] Default

set PPEventHandler:LambdaQCD 0.22
set PPEventHandler:RMax 3.0
set stdEmitter:PMinusOrdering 1.0
set stdEmitter:PTScale 1.0

set PPGenerator:EventHandler:PreSamples 0
# set PPGenerator:EventHandler:PreSamples 0
set PPGenerator:EventHandler:EventFiller:DebugHist 1
set PPTestTune:NumberOfEvents 10000

saverun compare PPTestTune
set PPEventHandler:FudgeME 2
set PPTestTune:EventHandler:EventFiller:OnlyOnce true
set PPTestTune:EventHandler:EventFiller:SingleMother 1
saverun compar2 PPTestTune
set PPTestTune:EventHandler:EventFiller:CompatMode 1
saverun compar3 PPTestTune
set PPTestTune:EventHandler:EventFiller:CompatMode 0
set PPTestTune:EventHandler:EffectivePartonMode Relatives
saverun compar4 PPTestTune
set PPTestTune:EventHandler:EventFiller:OnlyOnce false
set PPEventHandler:FudgeME 1
set PPTestTune:EventHandler:EventFiller:SingleMother 0
set PPTestTune:EventHandler:EventFiller:CompatMode 0


set PPTestTune:EventHandler:EffectivePartonMode Shadows
saverun shadows PPTestTune
set PPEventHandler:FudgeME 2
set PPTestTune:EventHandler:EventFiller:OnlyOnce true
set PPTestTune:EventHandler:EventFiller:SingleMother 1
saverun shadow2 PPTestTune
saverun shadow3 PPTestTune
set stdEmitter:MinusOrderingMode EffectiveParton
saverun shadow4 PPTestTune
set stdEmitter:MinusOrderingMode TrueShadow
saverun shadow6 PPTestTune
set stdEmitter:MinusOrderingMode OrderedShadow
saverun shadow7 PPTestTune
set stdXSec:IntOrdering ShadowOpen
saverun shadow8 PPTestTune
set stdXSec:IntOrdering ShadowColour
saverun shadow9 PPTestTune
set PPTestTune:NumberOfEvents 2
set PPGenerator:EventHandler:PreSamples 1000
saverun shadow10 PPTestTune
set PPGenerator:EventHandler:PreSamples 0
set PPTestTune:NumberOfEvents 10000
set stdXSec:IntOrdering VeryOpen
set stdEmitter:MinusOrderingMode EffectivePT
saverun shadow5 PPTestTune
set PPEventHandler:FudgeME 1
set PPTestTune:EventHandler:EventFiller:OnlyOnce false
set PPTestTune:EventHandler:EventFiller:SingleMother 0
set PPTestTune:EventHandler:EffectivePartonMode Colours

set PPTestTune:EventHandler:EventFiller:Mode NewSingle
saverun compsng PPTestTune
set PPEventHandler:FudgeME 2
set PPTestTune:EventHandler:EventFiller:OnlyOnce true
set PPTestTune:EventHandler:EventFiller:SingleMother 1
saverun compsn2 PPTestTune
set PPTestTune:EventHandler:EventFiller:CompatMode 1
saverun compsn3 PPTestTune
set PPEventHandler:FudgeME 1
set PPTestTune:EventHandler:EventFiller:OnlyOnce false
set PPTestTune:EventHandler:EventFiller:SingleMother 0
set PPTestTune:EventHandler:EventFiller:CompatMode 0

