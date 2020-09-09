#!/usr/bin/env python3
import pandas as pd
import os
import sys
from datetime import date
import subprocess
import numpy as np


class VCF():
    
    def __init__(self):
        self.breakpoints = []
        self.original_files = []
        self.header = ""
        self.pandasVCF = None
        
        
    def read_breakpoints(self,rmap_file):
        """
        This function chromosome breakpoints from a recombination map file (custom from SLiM).
        """
        with open(rmap_file, 'r' ) as r_file:
            next(r_file)
            for line in r_file:
                end, rate = line.rstrip().split()
                if rate == "0.5":
                    self.breakpoints.append( int(end) + 1 )

    def create_header(self):
        """
        This function creates an appropiate VCF header.
        """
        fileFormat = "##fileformat=VCFv4.2"
        fileDate = "##fileDate=" + date.today().strftime("%Y%m%d")
        source = "##source=" + sys.argv[0].split("/")[-1]
        f = '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
        self.header = "\n".join([fileFormat,fileDate,source,f]) + "\n"
    

    def modify_chromosomes(self,input_vcf):
        """
        This function takes a -position sorted- VCF pandas df and a breakpoints list.
        It converts 1-chromosome coordinates (SLiM) to multiple chromosomes, defined by breaking points of a recombination rate of 0.5.
        Breakpoints shold be a list of 1-chromosome coordinates that delimits different chromosomes.
        Function assumes vcf is position-sorted.
        """ 
        new_df = input_vcf
        # Preallocate position and chromosome lists
        positions = list(new_df["POS"])
        new_positions = [None]*len(positions)
        new_chromosomes = [None]*len(positions)
        brk_counter = 0
        newchrom = 1
        end = False
        for j in range(0,len(positions)):
            oldpos = positions[j]
            if oldpos < self.breakpoints[brk_counter]:
                if brk_counter == 0:
                    newpos = oldpos
                else:
                    newpos = oldpos - self.breakpoints[brk_counter - 1]
                new_positions[j] = int(newpos)
                new_chromosomes[j] = newchrom
            else:
                brk_counter += 1
                # Set maximum break counter number
                if brk_counter > len(self.breakpoints) - 1:
                    brk_counter = len(self.breakpoints) - 1
                    newpos = oldpos - self.breakpoints[brk_counter]
                    newchrom = len(self.breakpoints) + 1 
                    end = True
                if not end:
                    newchrom += 1
                    newpos = oldpos - self.breakpoints[brk_counter -1 ]
                new_positions[j] = int(newpos)
                new_chromosomes[j] = newchrom
        # Update POS,  CHROM, and rename properly (vcf standards)
        new_df["POS"] = new_positions
        new_df["#CHROM"] = new_chromosomes


    

    def read_merge_vcf(self,vcf_files, rmap=None ):
        """
        This function merges multiple VCF files into one. This is not a general purpose function. The VCF files that this function expects are VCF outputted by SLiM with the outputVCFSampel() function.
    
        This function assumes that VCF headers and 0-8 columns specifications are the same in every file (specific rows may not be the same, i.e. different loci present)
    
        Input: 
            
          vcf_files  list of vcf files with population name as their basename (${population}.vcf)
          outfile output file name
        """
        prettyprinter = "#"*20 + "\n" +  "{0}" + "\n"
        if rmap:
            sys.stdout.write(  prettyprinter.format("READING RECOMBINATION MAP")  )
            ## First read recombination map
            self.read_breakpoints( rmap )
    
        ## Now parse vcf
        sys.stdout.write(  prettyprinter.format("PARSING VCF FILES...")  )
             
        vcf_keys = "#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT".split()
        
        # Create header
        self.create_header() 
        
        for i in range(len(vcf_files)):
            pop = os.path.splitext( os.path.basename(vcf_files[i]) )[0]
            vcf = vcf_files[i]
            sys.stdout.write(  "--> PARSING FILE {0}\n".format(vcf)  )
            if i == 0:
                # Parse first file into table
                final_vcf =  pd.read_table( vcf, skiprows = 12 ,sep = "\t")
                # Remove info fields
                final_vcf["INFO"] = "."
                # Set indices
                final_vcf = final_vcf.set_index( vcf_keys )
                # Rename appending population of origin
                final_vcf = final_vcf.rename( mapper =  lambda x: pop + ":" + pop+"-"+x if (x.count("i")) == 1 else x , axis = 1)
            else:
                # Parse file into table
                tmp =  pd.read_table( vcf, skiprows = 12 ,sep = "\t" )
                # Remove info fields
                tmp["INFO"] = "." 
                 # Set indices
                tmp = tmp.set_index( vcf_keys )
                # Rename appending population of origin
                tmp = tmp.rename( mapper =  lambda x: pop + ":" + pop+"-"+x if (x.count("i")) ==     1 else x , axis = 1)
                # Merge outer, keep all genomic regions 
                final_vcf = pd.merge( tmp, final_vcf, how = "outer" , on  = tmp.index.names )
    
        # Fill nas with 0|0 (Not recorded by slim. Phased b/c homozygozyty)
        sys.stdout.write(  "-> ADDING MISSING (ANCESTRAL) SNPS\n"  )
        final_vcf = final_vcf.fillna( value = "0|0" )
        # reset indices and sort
        sys.stdout.write(  "-> RESETTING INDICES AND SORTING\n"  )
        
        final_vcf = final_vcf.reset_index()
        final_vcf = final_vcf.astype({"POS":"int64"})
        final_vcf = final_vcf.sort_values(by="POS", axis = 0 )
    
        # Modify chromosomes    
        if rmap:
            sys.stdout.write(  "-> ADJUSTING CHROMOSOME NAMES AND POSITIONS\n"  )
            self.modify_chromosomes( final_vcf )
    
        # add SNP ids
        sys.stdout.write(  "-> MODIFYING SNP IDs\n"  )
        final_vcf["ID"] = "snp_" + final_vcf["#CHROM"].astype(str)  + "_" + final_vcf["POS"].astype(str)
           
        
        self.pandasVCF = final_vcf
        
        
    def write_vcf(self, outfile):
        """
        Parameters
        ----------
        outfile : string
            Output filename for the VCF output.

        Returns
        -------
        None.
        """
        self.pandasVCF.to_csv( outfile, sep = "\t", header = True, index = False )
        # Write header
        sys.stdout.write(  "WRITING HEADER\n"  )
        with open( outfile, 'r+') as f:
            content = f.read()
            f.seek(0, 0)
            f.write( "".join( self.header ) + content)
    
        sys.stdout.write(  "#"*20 + "\n" + "DONE!\n" + "#"*20 + "\n"  )
    

def theoretical_covariance(tree, popsize):
    """
    Parameters
    ----------
    tree : dendropy.datamomdel.treemodel.Tree (dendropy tree)
        Dendropy tree that has been read by trees.read_tree(). Which means nodes are properly labeled.
    
    popsize: int.
        (Constant) effectuve population size throughout simulations.
    
    Returns
    -------
    Theoretical covariance matrix for a given tree.
    """
    
    init_gen = tree.nodes()[0].generation # Get first generation (actually burnin time)
    pops = tree.leaves # List populations
    m = np.zeros( ( len(pops),len(pops) ) )   # Matrix for storing distances
    pdm = tree.phylogenetic_distance_matrix() # Compute patristic distances 
    for idx, itx in enumerate(tree.taxon_namespace): # Loop through every leaf twice
        for jdx , jtx in enumerate(tree.taxon_namespace):
            if idx >= jdx:          # As matrix is symmetrical, only compute lower half
                mrca_time = pdm.mrca( itx , jtx ).generation - init_gen
                m[idx,jdx] = mrca_time
                m[jdx,idx] = mrca_time
                
            else: 
                pass
    
    m = pd.DataFrame( m )
    m.columns = pops
    m.index = pops
    
    c = lambda t,Ne: t/(2*Ne)
    m = c(m, popsize )
    m = m.drop(labels="p1", axis = 0)
    m = m.drop(labels="p1", axis = 1)
    
    
    return m 







def merge_vcf( vcf_files , outfile):
    """
    Merge several vcf files into a single VCF file using pandas.
    
    Parameters
    ----------
    vcf_files : list
        List of vcf filenames for merging and writing.
    outfile : str
        Output filename.

    Returns
    -------
    None.
    """
    vcf = VCF()
    vcf.read_merge_vcf(vcf_files)
    vcf.write_vcf(outfile)
    
def vcf_to_bed(vcf_file, output_prefix):
    """
    Use plink to convert a vcf file into a bed file with accessory files included.
    Parameters
    ----------
    vcf_file : str
        Input VCF file.
    output_prefix : str
        output prefix for plink file generation.

    Returns
    -------
    None.
    """
    
    s = ["plink", "--vcf", vcf_file,
         "--id-delim", ":", "--make-bed","--out",output_prefix]
    sys.stdout.write(  "CONVERTING {0} to bed format...\n".format(vcf_file)  )
    subprocess.run(s)
    sys.stdout.write(  "DONE!\n"  )


def ld_prune(bed_prefix, list_prefix,
             window_size,
             step_size,
             r_squared):
    """
    Wrapper for LD pruning of SNPs using PLINK.
    Generates lists of excluded and included SNPs.
    """
    
    sys.stdout.write("PERFORMING LD PRUNING...\n")    

    s = ["plink", "--bfile", bed_prefix, "--out", list_prefix,
         "--indep-pairwise", str(window_size), 
         str(step_size), str(r_squared) ]

    sys.stdout.write("--> COMPUTING CORRELATIONS...\n")
    subprocess.run(s)

    sys.stdout.write("DONE!\n")
        


def plink_extract( bed_prefix , list_file , bed_outfile_prefix  ):
    """
    Wrapper for extracting SNPs from a bed file specified in a list.
    """
    sys.stdout.write("--> EXTRACTING INDEPENDENT SNPs...\n")
    
    list_file = list_file.with_suffix('.prune.in')
    
    
    s = ["plink", "--bfile", bed_prefix, 
         "--extract", list_file, "--make-bed",
         "--out", bed_outfile_prefix ]
    subprocess.run(s)

def run_hapFLK(file_prefix,outfile_prefix,d):
    """
    Wrapper for running hapFLK.
    
    Parameters
    ----------
    file_prefix: str
        File prefix of input file(s).
    d : dict
        Dictionary that contains parameters for running hapflk.
    
    Returns
    -------
    None.
    """
    assert type(d) is dict
    
    ncpu = str(d["ncpu"])
    reynold_snps = str(d["reynold_snps"])
    K = str(d["K"])
    nfit = str(d["nfit"])
    covariance_file = d["covariance_file"]
    
    if covariance_file:
        s = ["hapflk", "--ncpu", ncpu,
             "--reynolds-snps", reynold_snps,
             "--bfile", file_prefix,
             "--prefix",outfile_prefix,
             "--outgroup", "p1",
             "-K", K , "--nfit", nfit,
             "--kinship", covariance_file,
             "--keep-outgroup"]
    else:
        s = ["hapflk", "--ncpu", ncpu,
             "--reynolds-snps", reynold_snps,
             "--bfile", file_prefix,
             "--prefix",outfile_prefix,
             "--outgroup", "p1",
             "-K", K , "--nfit", nfit,
             "--keep-outgroup"]
        
    
            
    sys.stdout.write(  "RUNNING hapFLK on {0}...\n".format(file_prefix)  )
    subprocess.run(s)
    sys.stdout.write(  "DONE!\n"  )
    
    
    
    
    
    
    
    
    