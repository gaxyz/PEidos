#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 13:44:15 2020

@author: gaston
"""

def initialize(mutation_rate,recombination_rate,genome_size):
    s = '''initialize(){{
    setSeed(seed);
    initializeMutationRate({0});
    initializeMutationType("m1",0.5,"f",0.0);
    initializeMutationType("m2",0.5,"f",0.0);
    initializeGenomicElementType("g1",m1,1.0);
    
    initializeGenomicElement(g1,0,{2});
    
    initializeRecombinationRate({1});
    
    }}'''.format(mutation_rate,
    recombination_rate,
    genome_size)
    
    
    return s


def setup(pop_size_filename, first_pop_name, first_pop_size):
    s = '''1 {{
    writeFile(paste0(c(outdir,"{0}")), "population size generation", append = F);
    defineConstant("simID", getSeed());
    sim.addSubpop("{1}",{2});
    }}'''.format(pop_size_filename,
                first_pop_name,
                first_pop_size)
    return s

def split(newpop,source,generation,popsize):
    s = '{0} {{ sim.addSubpopSplit("{1}",{2}, {3} );}}'.format(generation,
                                                            newpop,
                                                            popsize,
                                                            source )
    return s   
    
def migration(destination,source,rate,generation):
    """ Migration is unidirectional """
    s ='{0} {{{1}.setMigrationRates( {2}, {3} );}}'.format(generation,
                                                         destination,
                                                         source,
                                                         rate)
    return s
    
def add_mutation(mutation_name, population, generation, mutation_site):
    s= ''''{0} late(){{sim.outputFull("/tmp/slim_"+simID+".txt");
                      target = sample({1}.genomes,1);
                      target.addNewDrawnMutation({2},{3});
                      }}'''.format(generation,
                                    population,
                                    mutation_name,
                                    mutation_site
                                    )
    return s

def supress_mutations(generation, mutation_name):    
    s = '{0}: mutation({1}){{return F;}}'.format(generation, mutation_name)
    return s
    
    
def modify_fitness( mutation_name,population,fitness,generation):
    s = '{0}: fitness({1},{2}){{ return {3} ; }} '.format(generation,
                                           mutation_name,
                                           population,
                                           fitness)
    return s
    
    
def end_simulation(generation,populations,outputSampleSize):
    
    writeFile = '{0}.outputVCFSample(sampleSize={1},replace=F,outputMultiallelics=F,filePath=paste0(c(outdir,"{0}.vcf")));'
    blockList = [] 
    for i in populations:
        blockList.append(writeFile.format(i,outputSampleSize))
    blockList.append("sim.simulationFinished();")
    
    s = '{0} late(){{\n'.format(generation) + "\n".join(blockList) + "}"
    
    return s
    

def check_for_establishment(generation,population, mutation_name,mutation_site):
    """Check for establishment of a given mutation. Also called burn-in."""
    s = '''{0}: late(){{
                    mut = sim.mutationsOfType( {1} ) ;
                    
                    if ( size(mut) == 1 )
                    {{
                    
                        if (sim.mutationFrequencies(NULL,mut) > 0.2 )
                        {{
                            cat( simID + ": ESTABLISHED -- PROCEEDING WITH SPLIT\n" ) ;
                            sim.deregisterScriptBlock(self);
                        }}
                    }}
                    else
                    {{
                            cat(simID + " GENERATION: "+ sim.generation +": LOST BEFORE ESTABLISHMENT -- RESTARTING \n");
                            // back to gen 1800
                            sim.readFromPopulationFile( "/tmp/slim_" + simID + ".txt" );
                            // start newly seeded run
                            setSeed( rdunif(1,0,asInteger(2^62) - 1 ));
                    
                            // re introduce the sweep mutation
                            target = sample( {2}.genomes, 1 );
                            target.addNewDrawnMutation( {1}, {3} );
                    }}
                    }}'''.format(generation,mutation_name,population,mutation_site)

    return s


