#!/usr/bin/env python

##########################################################################
# Author: joemartaganna (joemar.ct@gmail.com)
#
# Description:
#
# Script for running AutoDock Vina directly against gzipped mol2 files
# to be more memory- and disk space-efficient. Individual mol2 files are
# extracted from a gzipped file downloaded from the ZINC database. Each
# file is then converted to pdbqt, then docked against a receptor protein.
# An affinity cut-off is set so that the intermediate files produced
# during docking with results below this cut-off are discarded. Only the
# output files from high affinity conformers (very few) are stored. 
#
# The ZINC id of the ligands docked are saved into a Redis database. Each
# docking run checks this database to avoid redocking the same ligands 
# if the script is rerun because of errors or accidental shutdown.
#
# There is an optional support for delivery of binders into Dropbox. This
# is useful for easy monitoring of results in cases when the virtual
# screening is run remotely on a server.
#
# Requires: vina, babel, joblib, redis, redispy, prepare_ligands
# Optional: dropbox
##########################################################################

import os
import sys
import json
import redis
import shutil
from joblib import Parallel, delayed
from subprocess import Popen, PIPE
from prepare_ligands import split_gzipped_mol2, convert_to_pdbqt


try:
    import dropbox
    dropbox_installed = True
except ImportError:
    dropbox_installed = False
    
    
def setup_dropbox_integration():
    script_path = os.path.realpath(__file__)
    conf_path = os.path.join(os.path.split(script_path)[0], '.DockingWatch.conf')
    app_key = '0raajjtnfqaq6ag'
    app_secret = 'ruc4zlu5gyk0rn3'
    flow = dropbox.client.DropboxOAuth2FlowNoRedirect(app_key, app_secret)
    if os.path.exists(conf_path):
        dbox_args = json.load(open(conf_path, 'r'))
        access_token = dbox_args['access_token']
        folder = dbox_args['folder']
    else:
        authorize_url = flow.start()
        # Authorize the Dropbox integration
        msg = '\nFor Dropbox integration, visit the following URL in your browser, '
        msg += 'grant access, and copy the resulting authorization code.'
        msg += ' Just hit enter without the code to skip.'
        print msg, '\n'
        print authorize_url
        code = raw_input("\nEnter the authorization code here: ").strip()
        folder = raw_input("Enter custom folder under DockingWatch app: ").strip()
        access_token, user_id = flow.finish(code)
        conf = {'access_token': access_token, 'folder': folder}
        json.dump(conf, open(conf_path, 'w'))
    return access_token, folder
    
    
def upload_to_dropbox(f, client, folder):
    fname = os.path.split(f)[-1]
    f = open(f, 'r')
    client.put_file(folder+'/binders/'+fname, f)
    

def dock_ligand(vina_conf, mol2_string, output_dir, affinity_cutoff):
    global db
    global dropbox_installed
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
            docking_out = ligand.replace('.pdbqt', '_out.pdbqt')
            if affinity < affinity_cutoff:
                os.remove(ligand)
                os.remove(docking_out)
            else:
                binders_dir = os.path.join(output_dir, 'binders')
                if not os.path.exists(binders_dir):
                    os.mkdir(binders_dir)
                shutil.move(ligand, binders_dir)
                shutil.move(ligand.replace('.pdbqt', '_out.pdbqt'), binders_dir)
                if dropbox_installed:
                    global client
                    global dropbox_folder
                    upload_to_dropbox(os.path.join(binders_dir, mol_id+'.pdbqt'), client, dropbox_folder)
                    upload_to_dropbox(os.path.join(binders_dir, mol_id+'_out.pdbqt'), client, dropbox_folder)
        else:
            print 'Error on %s: The conversion from mol2 to pdbqt may not have succeeded.' % mol_id
        
        
def execute_virtual_screening(vina_conf, gzipped_mol2, output_dir=None, affinity_cutoff=0):
    mols = split_gzipped_mol2(gzipped_mol2, output_dir=output_dir)
    print 'There are %s ligands in this file.' % len(mols)
    Parallel(n_jobs=-1, backend='multiprocessing', verbose=55)(delayed(dock_ligand)(vina_conf, mol, output_dir, affinity_cutoff) for mol in mols)
    
    
def connect_to_redisdb():
    db = redis.Redis(host='localhost', port=6379, db=7)
    return db

        
if __name__ == '__main__':

    conf = sys.argv[1]
    infiles = sys.argv[2:]
    # Set the output directory
    output_dir = os.path.join(os.getcwd(), 'docking_results')
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    # Start the redis db
    redis_conf = os.path.join(output_dir, 'redis.conf')
    # Setup Dropbox integration
    if dropbox_installed:
        access_token, dropbox_folder = setup_dropbox_integration()
        if access_token:
            client = dropbox.client.DropboxClient(access_token)
    # Execute virtual screening
    if os.path.exists(redis_conf):
        db = connect_to_redisdb()
        for i, f in enumerate(infiles):
            print '\n\n%s/%s - dealing with %s' % (i+1, len(infiles), f)
            execute_virtual_screening(conf, f, output_dir=output_dir, affinity_cutoff=7.5)
    else:
        rc = open(redis_conf, 'w')
        rc.write('dbfilename records.rdb\ndir %s\nsave 900 1' % output_dir)
        rc.close()
        print '\nRedis configuration file has been created.'
        print 'Run redis-server using the redis.conf file and rerun this script to start the virtual screening.\n'
            