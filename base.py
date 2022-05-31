from exceptions import PDBQTReadingError
import os.path as osp
import itertools
from itertools import groupby


class Atom(object):
    '''
    This class is basic Atom class in package.
    params:
        serial_number   : atom serial number
        name            : atom name
        altloc          : atom alternative location
        res_name        : atom residue name
        chain_id        : atom chain identifier
        res_seq_num     : atom residue sequence number
        res_insertion   : code for insertion of residues
        x               : atom X coordinate
        y               : atom Y coordinate
        z               : atom Z coordinate
        occupancy       : atom occupancy value
        temp_factor     : atom temperature factor
        partial_charge  : atom partial charge
        atom_type       : atom type
    '''

    def __init__(self, *args):
        self.serial_number = args[0]
        self.name = args[1]
        self.altloc = args[2]
        self.res_name = args[3]
        self.chain_id = args[4]
        self.res_seq_num = args[5]
        self.res_insertion = args[6]
        self.x = args[7]
        self.y = args[8]
        self.z = args[9]
        self.occupancy = args[10]
        self.temp_factor = args[11]
        self.partial_charge = args[12]
        self.atom_type = args[13]

    def __str__(self):
        out = "=" * 50 + "\n" + \
            "Serial Number: " + str(self.serial_number) + "\n" + \
            "Name: " + self.name + "\n" + \
            "Alternative Location: " + self.altloc + "\n" + \
            "Residue Name: " + self.res_name + "\n" + \
            "Chain Identifier: " + self.chain_id + "\n" + \
            "Residue Sequence Number: " + str(self.res_seq_num) + \
            "\n" + \
            "Code for Insertion of Residues: " + self.res_insertion + \
            "\n" + \
            "X Coordinate: " + str(self.x) + "\n" + \
            "Y Coordinate: " + str(self.y) + "\n" + \
            "Z Coordinate: " + str(self.z) + "\n" + \
            "Occupancy: " + str(self.occupancy) + "\n" + \
            "Temperature Factor: " + str(self.temp_factor) + "\n" + \
            "Partial Charge: " + str(self.partial_charge) + "\n" + \
            "Atom Type: " + self.atom_type + "\n" + \
            "=" * 50 + "\n"
        return out

def int1(instring):
    try:
        return int(instring)
    except:
        return -1

def float1(instring):
    try:
        return float(instring)
    except:
        return 0.

class Structure(object):
    '''
    This is the basic Structure class in this package
    params:
        pdbqt_file  : The address of pdbqt file on disk
        site_file   : mol2 file contains site atoms (applicable just for training phase)
        mol2_file   : mol2 file contains site
    '''
    def __init__(self, pdbqt_file, site_file=None):
        self.binding_residues = []
        self.pdbqt_file = pdbqt_file
        self.site_file = site_file
        self.n_atoms = 0
        self.n_chains = 0
        self.residues = None

    def read_pdbqt(self):
        '''
        this method reads pdbqt file and returns list of atom information
        '''
        try:
            with open(self.pdbqt_file, "r") as f:
                lines = f.readlines()
            self.atoms = [Atom(
                int1(line[6:11]),  # 0  atom serial number
                line[12:16].strip(),  # 1  atom name
                line[16],  # 2  alternate location indicator
                line[17:21].strip(),  # 3  residue name
                line[21],  # 4  chain identifier
                int1(line[22:26]),  # 5  residue sequence number
                line[26],  # 6  code for insertion of residues
                float1(line[30:38]),  # 7  x
                float1(line[38:46]),  # 8  y
                float1(line[46:54]),  # 9  z
                float1(line[54:60]),  # 10 occupancy
                float1(line[60:66]),  # 11 temperature factor
                float1(line[66:76]),  # 12 partial charge
                line[77:].strip())  # 13 autodock atom type
                for line in lines if line.startswith("ATOM")]
            self.n_atoms = len(self.atoms)
            self.n_chains = len(set([atom.chain_id for atom in self.atoms]))
        except:
            raise PDBQTReadingError

    def _correct_residues(self, original_mol2_file):
        """autodock_tools_py3 has some problems with residue id and chain assignment
            this method solves it. Just applicable on training phase"""
        with open(original_mol2_file, "r") as f:
            lines = f.readlines()
        start = 0
        end = -1
        for i, line in enumerate(lines):
            if line.startswith('@<TRIPOS>SUBSTRUCTURE'):
                start = i + 1
            elif line.startswith('@<TRIPOS>SET'):
                end = i
        lines = lines[start:end]
        lines = [line.strip().split() for line in lines]
        
        for i in range(len(lines) - 1):
            start = int(lines[i][2]) - 1
            end = int(lines[i+1][2]) - 1
            res_seq_num = int(lines[i][-1])
            chain_id = lines[i][-2]
            for atom in self.atoms[start:end]:
                atom.res_seq_num = res_seq_num
                atom.chain_id = chain_id
        res_seq_num = int(lines[i + 1][-1])
        chain_id = lines[i + 1][-2]
        for atom in self.atoms[end:]:
            atom.res_seq_num = res_seq_num
            atom.chain_id = chain_id  
 
    def clean_altloc(self):
        '''
        This method handles altlocs
        '''

        atoms = [[atom.chain_id, atom.res_seq_num, atom.name, i]
                 for i, atom in enumerate(self.atoms)]
        atoms = sorted(atoms, key=lambda el: (el[0], el[1], el[2]))
        indices = []
        groups = itertools.groupby(atoms, lambda el: (el[0], el[1], el[2]))

        for k, g in groups:
            g = list(g)
            index = sorted(g, key=lambda el: el[-1])[0][-1]
            indices.append(index)

        self.atoms = [self.atoms[i] for i in indices]
        self.n_atoms = len(self.atoms)

    def read_site(self):
        '''
        This method reads site.mol2 file
        '''

        with open(self.site_file, "r") as f:
            lines = f.readlines()
        start = 0
        end = -1
        for i, line in enumerate(lines):
            if line.startswith('@<TRIPOS>SUBSTRUCTURE'):
                start = i + 1
            elif line.startswith('@<TRIPOS>SET'):
                end = i
        lines = lines[start:end]
        for line in lines:
            line = line.strip()
            self.binding_residues.append(
                (line[-5], int(line[-4:])))

    def compute_residues(self):
        '''
        This method fills self.residues attribute
        '''
        if not self.residues:
            self.residues = list(
                set([(atom.chain_id, atom.res_seq_num) for atom in self.atoms]))

    def get_res_name_groups(self, res_name):
        '''
        This method groups atoms based on (chain_id, res_seq_num) for a 
        specific residue (e.g. ALA, PHE, ...)
        params:
            res_name:   residue name (e.e. ALA, PHE, ...)
        '''

        res_name_atoms = list(
            filter(lambda atom: atom.res_name == res_name, self.atoms))
        res_name_atoms = sorted(res_name_atoms, key=lambda atom: (
            atom.chain_id, atom.res_seq_num, atom.name))
        groups = groupby(res_name_atoms, key=lambda atom: (
            atom.chain_id, atom.res_seq_num))
        res_name_groups = [(k, tuple(g)) for k, g in groups]
        return res_name_groups

    def __len__(self):
        '''
        returns number of atoms
        '''
        return len(self.atoms)

    def __str__(self):
        out = "=" * 50 + "\n" \
            "structure name: " + osp.basename(self.pdbqt_file)[:-6] + "\n" + \
            "number of atoms: " + str(self.n_atoms) + "\n" + \
            "number of chains: " + str(self.n_chains) + "\n"
        return out
