#!/usr/bin/env python

import sys
from subprocess import Popen, PIPE


def convert_to_pdbqt(receptor_pdb):
    pdbqt_outf_path = receptor_pdb.replace('.pdb', '.pdbqt')
    cmd = 'babel -ipdb %s -opdbqt %s --partialcharge gasteiger --AddPolarH' % (receptor_pdb, pdbqt_outf_path)
    p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = p.communicate()
    print '\nCompleted! Output PDBQT file has been written.\n'


if __name__ == '__main__':
    try:
        receptor_pdb = sys.argv[1]
        convert_to_pdbqt(receptor_pdb)
    except ValueError:
        print '\nUsage: prepare_receptor.py <receptor_pdb_file>\n'