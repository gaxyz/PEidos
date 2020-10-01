#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 16:25:06 2020

@author: gaston
"""

import argparse
# Load peidos

from peidos import trees 
from peidos import utils as utils

# Argparse first
parser = argparse.ArgumentParser( description = 'Compute theoretical covariance matricex from nwk tree. Outgroup is excluded and its label assumed to be "p1".')
parser.add_argument("popsize", type = str, help = "population size" )
parser.add_argument("firstSplitGen", type = str , help ="generation where the first split occurs (outgroup)"  )
parser.add_argument("outMatrix", type =str, help = "Output matrix file compatible with hapFLK.", default = "covariance.tab")

args = parser.parse_args()
N=args.popsize
firstSplitGen=args.firstSplitGen
outfile=args.outMatrix

tree = trees.read_tree( "serial_splits.nwk", firstSplitGen )
cov_matrix = utils.theoretical_covariance(tree, N)
cov_matrix.to_csv(outfile, header = False, sep  = " ")

