import os
import unittest
from prepare_ligands import convert_to_pdbqt


class LigandPrepTests(unittest.TestCase):
    
    def test_babel(self):
        mol2 = open('test_files/ZINC67842136.mol2', 'r').read()
        ligand_paths = convert_to_pdbqt(mol2, 'test_files', remove_input=False)
        pdbqt = ligand_paths['pdbqt']
        pdbqt_contents = open(pdbqt, 'r').read()
        os.remove(pdbqt)
        self.assertTrue(len(pdbqt_contents) > 0)
        

if __name__ == '__main__':
    unittest.main()