#!/bin/bash

#python makeDiphoPlots_Jan.py -k ggH_PU0 -w
#python makeDiphoPlots_Jan.py -k ggH_PU200 -w
#python makeDiphoPlots_Jan.py -k VBF_PU0 -w
#python makeDiphoPlots_Jan.py -k VBF_PU200 -w
#python makeDiphoPlots_Jan.py -k DiPhoJetsBox_PU200 -w
#python makeDiphoPlots_Jan.py -k GJet_PU200 -w

#python makeDiphoPlots_Jan.py -k ggH_PU0 -w --doLoose
#python makeDiphoPlots_Jan.py -k ggH_PU200 -w --doLoose
#python makeDiphoPlots_Jan.py -k VBF_PU0 -w --doLoose
#python makeDiphoPlots_Jan.py -k VBF_PU200 -w --doLoose
#python makeDiphoPlots_Jan.py -k GJet_PU200 -w --doLoose
#python makeDiphoPlots_Jan.py -k DiPhoJetsBox_PU200 -w --doLoose

#python makeFlatTree.py -k ggH_PU200 -w -t 10000
#python makeFlatTree.py -k VBF_PU200 -w -t 10000
#python makeFlatTree.py -k DiPhoJetsBox_PU200 -w -t 10000
#python makeFlatTree.py -k GJet_PU200 -w -t 10000

#python makeFlatTree.py -k ggH_PU200 -w --forTesting
#python makeFlatTree.py -k VBF_PU200 -w --forTesting
#python makeFlatTree.py -k DiPhoJetsBox_PU200 -w --forTesting
#python makeFlatTree.py -k GJet_PU200 -w --forTesting

python getEfficiencies.py -k ggH_PU200
python getEfficiencies.py -k VBF_PU200
#python getEfficiencies.py -k DiPhoJetsBox_PU200
#python getEfficiencies.py -k GJet_PU200
