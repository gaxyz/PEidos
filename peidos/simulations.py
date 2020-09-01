#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 18:18:41 2020

@author: gaston
"""

from . import events as ee
import os
from pathlib import Path
import subprocess


def treelike_neutral(tree, d, script_file ): 
    """
    Write a very specific simulation scheme (tree-like, migration).
    
    Tree is properly labeled DendroPy tree (use read_tree()).
    d is a configuration dictionary with simulation parameters.
    Script file is the output name of the file.
    """    
    slim_script = ee.SlimScript()

    # Frequently used parameters
    popsize = d["simulation_popsize"]

    # Initialize
    slim_script.initialize(d["mutation_rate"],
                           d["recombination_rate"],
                           d["genome_size"])  
    # Setup
    slim_script.setup(d["pop_size_filename"],
                      tree.seed_node.slabel,
                      popsize)    
    
    # Traverse tree and write splits
    for node in tree.preorder_node_iter():
        
        # If node has child nodes (omit leaves)
        if node.child_nodes():
            # List child nodes
            childs = node.child_nodes()
            # Tag source population
            source_pop = node.slabel
            # Tag generation of split
            generation = node.generation
            # For every child that is not the source population, write a split
            for child in childs:
                if source_pop != child.slabel:
                    slim_script.split( child.slabel,
                                       source_pop,
                                       generation,
                                       popsize )
    # End simulation
    total_tree_length = d["burnin_time"] + int(tree.length()) # This is total number of generations
    slim_script.end_simulation(total_tree_length,
                               tree.leaves,
                               d["output_sample_size"]) 
    
    slim_script.write_script( script_file )
    




def migration_pulse_neutral(tree, d, script_file ): 
    """
    Write a very specific simulation scheme (tree-like, migration).
    
    Tree is properly labeled DendroPy tree (use read_tree()).
        d is a configuration dictionary with simulation parameters.
    Script file is the output name of the file.
    """    
    slim_script = ee.SlimScript()

    # Frequently used parameters
    popsize = d["simulation_popsize"]

    # Initialize
    slim_script.initialize(d["mutation_rate"],
                           d["recombination_rate"],
                           d["genome_size"])  
    # Setup
    slim_script.setup(d["pop_size_filename"],
                      tree.seed_node.slabel,
                      popsize)    
    
    # Traverse tree and write splits
    for node in tree.preorder_node_iter():
        
        # If node has child nodes (omit leaves)
        if node.child_nodes():
            # List child nodes
            childs = node.child_nodes()
            # Tag source population
            source_pop = node.slabel
            # Tag generation of split
            generation = node.generation
            # For every child that is not the source population, write a split
            for child in childs:
                if source_pop != child.slabel:
                    slim_script.split( child.slabel,
                                       source_pop,
                                       generation,
                                       popsize )
    
    # Write migration pulse
    
    slim_script.admixture_pulse(d["pulse_destination"],
                                d["pulse_source"],
                                d["pulse_rate"],
                                d["pulse_generation"])
    
    
    
    # End simulation
    total_sim_generations = d["pulse_generation"] + d["post_admixture_time"]  # This is total number of generations
    total_sim_generations = int(total_sim_generations)
    slim_script.end_simulation(total_sim_generations,
                               tree.leaves,
                               d["output_sample_size"]) 
    
    slim_script.write_script( script_file )
    
    
    
    

def slimsim(script_file, seed, outdir,tag):
    
    """
    Function for running a single instance of a given SLiM simulation.
    Written in a way that can be easily paralellizable, each function call a process.
    """
    
    
    # Create simulation - specific output folder
    out_folder = Path(outdir) / "{0}_{1}".format(tag,seed)
    try:
        os.mkdir(out_folder)
    except FileExistsError:
        print( "Directory {0} already exists.".format(str(out_folder)))
    except OSError:
        print( "Creation of directory {0} failed.".format(str(out_folder)))
    
    
    
    # Run command
    command = "slim -d seed={0} -d outdir='{1}' {2}".format(seed,
                                                          str(out_folder)+"/",
                                                          script_file).split()
    subprocess.run(command)
    


