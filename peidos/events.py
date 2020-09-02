#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 13:44:15 2020

@author: gaston
"""

class SlimScript():
    
    def __init__(self):
        self.initialization = ""
        self.setup_config = ""
        self.splits = []
        self.migrations = []
        self.mutation = ""
        self.supressed_mutation = ""
        self.fitness = ""
        self.end = ""
        self.check_establishment = ""
        

    def initialize(self,mutation_rate,recombination_rate,genome_size):
        s = '''initialize(){{
        setSeed(seed);
        initializeMutationRate({0});
        initializeMutationType("m1",0.5,"f",0.0);
        initializeMutationType("m2",0.5,"f",0.0);
        initializeGenomicElementType("g1",m1,1.0);
        
        initializeGenomicElement(g1,0,{2});
        
        initializeRecombinationRate({1});
        
        }}\n'''.format(mutation_rate,
        recombination_rate,
        genome_size)
        
        self.initialization = s



    def setup(self,pop_size_filename, first_pop_name, first_pop_size):
        s = '''1 {{
        writeFile(paste0(c(outdir,"{0}")), "population size generation",append = F);
        defineConstant("simID", getSeed());
        sim.addSubpop("{1}",{2});
        }}\n'''.format(pop_size_filename,
                    first_pop_name,
                    first_pop_size)
        self.setup_config = s


    def split(self,newpop,source,generation,popsize):
        s = '{0} {{ sim.addSubpopSplit("{1}",{2}, {3} );}}\n'.format(generation,
                                                                newpop,
                                                                popsize,
                                                                source )
        self.splits.append(s)   
    

    def admixture_pulse(self,destination,source,rate,generation):
        
        """This function turns a population leaf into an admixed leaf.
           Source population will still be present as and unadmixed population,
           while the destination population will be the admixed one.
           
           As it is a pulse, the migration rate should be relatively
           high and shortlasting."""
        s ="""{0} {{{1}.setMigrationRates( {2}, {3} );}}
{4} {{{1}.setMigrationRates( {2}, {5} );}}\n""".format(generation,
                                                             destination,
                                                             source,
                                                             rate,
                                                             generation + 1,
                                                             0)

    def migration(self,destination,source,rate,generation):
        """ Migration is unidirectional """
        s ='{0} {{{1}.setMigrationRates( {2}, {3} );}}\n'.format(generation,
                                                             destination,
                                                             source,
                                                             rate)
        self.migrations.append(s)
    
    def add_mutation(self,mutation_name, population, generation, mutation_site):
        s= ''''{0} late(){{sim.outputFull("/tmp/slim_"+simID+".txt");
                          target = sample({1}.genomes,1);
                          target.addNewDrawnMutation({2},{3});
                          }}\n'''.format(generation,
                                        population,
                                        mutation_name,
                                        mutation_site
                                        )
        self.mutation = s

    def supress_mutations(self,generation, mutation_name):    
        s = '{0}: mutation({1}){{return F;}}\n'.format(generation, mutation_name)
        self.supressed_mutation = s
    

    
    
    def modify_fitness( self,mutation_name,population,fitness,generation):
        s = '{0}: fitness({1},{2}){{ return {3} ; }}\n'.format(generation,
                                               mutation_name,
                                               population,
                                               fitness)
        self.fitness = s
    
    
    def end_simulation(self,generation,populations,outputSampleSize):
        
        writeFile = '{0}.outputVCFSample(sampleSize={1},replace=F,outputMultiallelics=F,filePath=paste0(c(outdir,"{0}.vcf")));'
        blockList = [] 
        for i in populations:
            blockList.append(writeFile.format(i,outputSampleSize))
        blockList.append("sim.simulationFinished();")
        
        s = '{0} late(){{\n'.format(generation) + "\n".join(blockList) + "}\n"
        
        self.end = s
    

    def check_for_establishment(self,generation,population, mutation_name,mutation_site):
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
                        }}\n'''.format(generation,mutation_name,population,mutation_site)
    
        self.check_establishment = s

    def write_script(self,script_file_name):
        with open( script_file_name , 'w' ) as handle:
            
            separator = "////////////////////////\n"
            # Write initialization
            handle.write(self.initialization )
            handle.write(separator)
            # Write setup
            handle.write(self.setup_config )
            handle.write(separator)
            
            # Write splits
            if self.splits:
                for split in self.splits:
                    handle.write(split)
                    handle.write(separator)
            # Write migrations
            if self.migrations:
                for mig in self.migrations:
                    handle.write(mig)
                    handle.write(separator)
            # Write mutation
            if self.mutation:
                handle.write(self.mutation )
                handle.write(separator)
            # Supress mutation
            if self.supressed_mutation:
                handle.write(self.supressed_mutation)
            
            # Modify fitness
            if self.fitness:
                handle.write(self.fitness)
            # Check establishment
            if self.check_establishment:
                handle.write(self.check_establishment)
                
            # End simulation
            if self.end:
                handle.write(self.end)
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            