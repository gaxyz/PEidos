#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 19:43:41 2020

@author: gaston
"""
import os
import sys

sys.path.insert( 0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from peidos import trees 
from peidos import simulations as sim

tree = trees.read_tree("test.nwk", 1800)


config_dict ={
    "mutation_rate":1e-7,
    "recombination_rate":1e-8,
    "genome_size":999999,
    "burnin_time":1800,
    "pop_size_filename":"pop_sizes.txt",
    "simulation_popsize":200,
    "output_sample_size":50
    }


script_file = "treelike_neutral.slim"

sim.treelike_neutral(tree,config_dict,script_file)

sim.slimsim(script_file, 2020,".","testing")


