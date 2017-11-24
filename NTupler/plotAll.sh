#!/bin/bash

python makeDiphoPlots_Jan.py -k ggH_PU0 -w
python makeDiphoPlots_Jan.py -k ggH_PU200 -w
python makeDiphoPlots_Jan.py -k VBF_PU0 -w
python makeDiphoPlots_Jan.py -k VBF_PU200 -w

python makeDiphoPlots_Jan.py -k ggH_PU0 -w --doLoose
python makeDiphoPlots_Jan.py -k ggH_PU200 -w --doLoose
python makeDiphoPlots_Jan.py -k VBF_PU0 -w --doLoose
python makeDiphoPlots_Jan.py -k VBF_PU200 -w --doLoose
