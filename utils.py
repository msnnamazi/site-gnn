import os
import glob
import time
import subprocess
import multiprocessing

from AutoDockTools.MoleculePreparation import AD4ReceptorPreparation
from MolKit import Read

from exceptions import *
from config import *


def run_pdbqt(io_pair, quiet=False):
    ''' 
        generates and saves pdbqt file from input file.
        params:
            io_pair: A tuple consists of (input_file, output_file)
    '''
    try:
        mol = Read(io_pair[0])
        AD4ReceptorPreparation(mol[0], 
                               charges_to_add=None, 
                               cleanup="", 
                               outputfilename=io_pair[1])
    except:
        if not(quiet):
            raise PDBQTGenerationError

def run_all_pdbqt(io_pairs):
    """
        generates pdbqt files for a list of tuples(io_pairs)
        params:
            io_pairs: list of tuples. each tuple is a (input_file, output_file) pair     
    """
    for io_pair in io_pairs:
        run_pdbqt(io_pair)

def run_all_pdbqt_p(io_pairs):
    """
        generates pdbqt files for a list of tuples(io_pairs)
        params:
            io_pairs: list of tuples. each tuple is a (input_file, output_file) pair     
    """
    with multiprocessing.Pool() as pool:
        pool.map(run_pdbqt, io_pairs)

def run_dssp(io_pair):  
    """
        runs dssp for computing secondary structure and sasa
        params:
            io_pair: a tuple consists of input_file address and output_file address
    """  
    try:
        subprocess.run([DSSP_EXECUTABLE, '-i', io_pair[0], '-o', io_pair[1]])
    except:
        pass

def run_all_dssp(io_pairs):
    """
        run dssp for a list of io_pairs in a single cpu core
        params:
            io_pairs: list of tuples. each tuple is a (input_file, output_file) pair            
    """
    for io_pair in io_pairs:
        if not(os.path.isfile(io_pair[1])):
            run_dssp(io_pair)

def run_all_dssp_p(io_pairs):
    """
        runs dssp on a list of io_pairs in parallel
        params:
            io_pairs: list of tuples. each tuple is a (input_file, output_file) pair     
    """
    with multiprocessing.Pool() as pool:
        pool.map(run_dssp, io_pairs)

def run_msms(io_pair):
    """
        runs msms for computing molecular surface
        params:
            io_pair: a tuple (input_pdb_file, output_name(without extension))
    """
    try:
        with open(io_pair[1] + '.xyzr', 'w') as f:
            subprocess.run([PDB_TO_XYZR_EXECUTABLE, io_pair[0]], stdout=f)
    except:
        return
    try:
        subprocess.run([MSMS_EXECUTABLE, "-if", io_pair[1] + '.xyzr', "-of", io_pair[1]], stdout=open(os.devnull, 'wb'))
    except:
        pass

def run_all_msms(io_pairs):
    """
        runs msms for computing molecular surfce of a list of io_pairs on a single cpu core
        params:
            io_pairs: list of (input_pdb_file, output_name(without extension)) tuples
    """
    for io_pair in io_pairs:
        run_msms(io_pair)

def run_all_msms_p(io_pairs):
    """
        runs msms for computing molecular surfce of a list of io_pairs in parallel
        params:
            io_pairs: list of (input_pdb_file, output_name(without extension)) tuples
    """
    with multiprocessing.Pool() as pool:
        pool.map(run_msms, io_pairs)






