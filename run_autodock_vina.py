#!/usr/bin/env python

##########################################################################
# author: joemartaganna (joemar.ct@gmail.com)
#
# Script for running AutoDock Vina directly against gzipped mol2 files
# to be more memory- and disk space-efficient. Individual mol2 files are
# extracted from a gzipped file downloaded from the ZINC database. Each
# file is then converted to pdbqt, then docked against a receptor protein.
# An affinity cut-off is set so that the intermediate files produced
# during docking with results below this cut-off are discarded. Only the
# outfiles files from high affinity conformers are stored.
#
# Requires: vina, babel, joblib, redis, redispy, prepare_ligands
##########################################################################

import os
import redis
import shutil
from joblib import Parallel, delayed
from subprocess import Popen, PIPE
from prepare_ligands import split_gzipped_mol2, convert_to_pdbqt


def dock_ligand(vina_conf, mol2_string, output_dir, affinity_cutoff):
    global db
    mol_id = mol2_string.split('\n')[1]
    try:
        # Check if this has been docked
        db[mol_id]
    except KeyError:        
        ligand_paths = convert_to_pdbqt(mol2_string, output_dir)
        ligand = ligand_paths['pdbqt']
        cmd = 'vina --config %s --ligand %s' % (vina_conf, ligand)
        p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
        stdout, stderr = p.communicate()
        if 'Refining results ... done' in stdout:
            lines = stdout.splitlines()
            best = None
            for i, line in enumerate(lines):
                if '-+-' in line:
                    best = lines[i+1]
            if best:
                affinity = float(best.split()[1]) * -1
            else:
                affinity = 0
        # Put this into the record
        db[mol_id] = affinity
        # Remove the input and output files if affinity is below cut-off
        if affinity < affinity_cutoff:
            os.remove(ligand)
            os.remove(ligand.replace('.pdbqt', '_out.pdbqt'))
        else:
            binders_dir = os.path.join(output_dir, 'binders')
            if not os.path.exists(binders_dir):
                os.mkdir(binders_dir)
            shutil.move(ligand, binders_dir)
            shutil.move(ligand.replace('.pdbqt', '_out.pdbqt'), binders_dir)
        
        
def execute_virtual_screening(vina_conf, gzipped_mol2, output_dir=None, affinity_cutoff=0):
    mols = split_gzipped_mol2(gzipped_mol2, output_dir=output_dir)
    print 'There are %s ligands in this file.' % len(mols)
    Parallel(n_jobs=-1, backend='multiprocessing', verbose=55)(delayed(dock_ligand)(vina_conf, mol, output_dir, affinity_cutoff) for mol in mols)
    
        
if __name__ == '__main__':
    import sys
    conf = sys.argv[1]
    infiles = sys.argv[2:]
    # Set the output directory
    output_dir = os.path.join(os.getcwd(), 'docking_results')
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    # Start the redis db
    redis_conf = os.path.join(output_dir, 'redis.conf')
    if os.path.exists(redis_conf):
        db = redis.Redis(host='localhost', port=6379, db=7)
        for i, f in enumerate(infiles):
            print '\n\n%s/%s - dealing with %s' % (i+1, len(infiles), f)
            execute_virtual_screening(conf, f, output_dir=output_dir, affinity_cutoff=7.5)
    else:
        rc = open(redis_conf, 'w')
        rc.write('dbfilename records.rdb\ndir %s\nsave 900 1' % output_dir)
        rc.close()
        print 'Redis configuration file has been created.'
        print 'Run redis-server using the redis.conf file and rerun this script to start the virtual screening.'