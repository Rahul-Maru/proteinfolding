import unittest
import sys
import os

# Add the project root to the Python path to allow importing lib and consts
# Adjust path to go up two levels (from lib/tests to project root)
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, project_root)

# Now imports relative to project root should work
from lib.Protein import Protein
# Attempt to import AA_MAP, provide a mock if it fails or is unavailable
try:
    from consts import AA_MAP
except ImportError:
    print("Warning: consts.AA_MAP not found. Using mock AA_MAP for tests.")
    AA_MAP = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
        'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
        'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
        'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
    }

# Test subclass that skips file reading and allows setting residues directly
class MockProtein(Protein):
    def __init__(self, test_residues):
        """Initializes a mock protein with predefined residues, skipping file I/O."""
        # Skip the original Protein.__init__ entirely to avoid file operations
        self.id = "TEST"
        self.path = "mock/path.pdb" # Dummy path
        self.display_mode = [False]*3
        self.lines = [] # Dummy value
        self.atoms = [] # Dummy value
        self.hetatms = [] # Dummy value

        # Set the residues directly. This works by setting the instance dictionary
        # entry that @cached_property would normally use to store the result.
        self.__dict__['residues'] = test_residues

    # Override methods that might be called by __init__ or others if they error
    def get_record(self, record: str) -> list[str]:
        return [] # Return empty list to avoid errors if called

    def get_xyzlist(self, query=None, qtype="", triplet = False):
         # Return dummy value to avoid errors if called unexpectedly
         return ([], [], []) if not triplet else []


class TestProteinSeq(unittest.TestCase):
    """Test suite for the Protein.seq() method."""

    def test_single_chain_no_gaps(self):
        """Test sequence with consecutive residues in one chain."""
        test_residues = [
            {'code': 'ALA', 'id': 'A_1', 'type': 'atom', 'atoms': []},
            {'code': 'CYS', 'id': 'A_2', 'type': 'atom', 'atoms': []},
            {'code': 'ASP', 'id': 'A_3', 'type': 'atom', 'atoms': []},
        ]
        protein = MockProtein(test_residues)
        expected_seq = "ACD"
        self.assertEqual(protein.seq(), expected_seq)

    def test_multiple_chains(self):
        """Test sequence spanning multiple chains with chain break marker."""
        test_residues = [
            {'code': 'ALA', 'id': 'A_1', 'type': 'atom', 'atoms': []},
            {'code': 'CYS', 'id': 'A_2', 'type': 'atom', 'atoms': []}, # End chain A
            {'code': 'ASP', 'id': 'B_5', 'type': 'atom', 'atoms': []}, # Start chain B
            {'code': 'GLU', 'id': 'B_6', 'type': 'atom', 'atoms': []},
        ]
        protein = MockProtein(test_residues)
        expected_seq = "AC|----DE" # Chain break + leading dashes for chain B
        self.assertEqual(protein.seq(), expected_seq)

    def test_gaps_within_chain(self):
        """Test sequence with gaps within a single chain."""
        test_residues = [
            {'code': 'ALA', 'id': 'A_1', 'type': 'atom', 'atoms': []},
            {'code': 'CYS', 'id': 'A_2', 'type': 'atom', 'atoms': []},
            {'code': 'GLU', 'id': 'A_5', 'type': 'atom', 'atoms': []}, # Gap of 2 (3, 4 missing)
            {'code': 'PHE', 'id': 'A_6', 'type': 'atom', 'atoms': []},
        ]
        protein = MockProtein(test_residues)
        expected_seq = "AC--EF" # Two dashes for missing 3 and 4
        self.assertEqual(protein.seq(), expected_seq)

    def test_single_residue(self):
        """Test sequence with only one residue and leading dashes."""
        test_residues = [
            {'code': 'GLY', 'id': 'A_10', 'type': 'atom', 'atoms': []},
        ]
        protein = MockProtein(test_residues)
        expected_seq = "---------G" # Nine dashes for positions 1-9
        self.assertEqual(protein.seq(), expected_seq)

    def test_empty_residues(self):
        """Test sequence with an empty residue list."""
        test_residues = []
        protein = MockProtein(test_residues)
        expected_seq = ""
        self.assertEqual(protein.seq(), expected_seq)

    def test_custom_residue_list(self):
        """Test calling seq() with a specific list of residues."""
        base_residues = [
            {'code': 'ALA', 'id': 'A_1', 'type': 'atom', 'atoms': []}, #0
            {'code': 'CYS', 'id': 'A_2', 'type': 'atom', 'atoms': []}, #1
            {'code': 'ASP', 'id': 'A_3', 'type': 'atom', 'atoms': []}, #2
            {'code': 'GLU', 'id': 'A_5', 'type': 'atom', 'atoms': []}, #3
        ]
        protein = MockProtein(base_residues) # Mock protein still needed

        # Pass only a subset to seq()
        custom_list = [
             base_residues[1], # CYS (A_2)
             base_residues[3], # GLU (A_5) - Gap relative to CYS (A_2)
        ]
        # Gaps are calculated based *only* on the list provided to seq()
        expected_seq = "-C--E" # Leading dash + gap between 2 and 5
        self.assertEqual(protein.seq(residues=custom_list), expected_seq)

    def test_no_gap_if_consecutive(self):
        """Test sequence with consecutive residues; ensure no gap is added."""
        test_residues = [
           {'code': 'ALA', 'id': 'A_10', 'type': 'atom', 'atoms': []},
           {'code': 'CYS', 'id': 'A_11', 'type': 'atom', 'atoms': []},
        ]
        protein = MockProtein(test_residues)
        expected_seq = "---------AC" # Nine dashes for positions 1-9
        self.assertEqual(protein.seq(), expected_seq)

    def test_gap_at_start_of_custom_list(self):
        """Test gap calculation when the gap is between the first two elements."""
        protein = MockProtein([]) # Base residues don't matter here
        custom_list = [
           {'code': 'ALA', 'id': 'A_3', 'type': 'atom', 'atoms': []},
           {'code': 'CYS', 'id': 'A_5', 'type': 'atom', 'atoms': []}, # Gap of 1 (residue 4)
        ]
        expected_seq = "--A-C" # Two leading dashes + gap between 3 and 5
        self.assertEqual(protein.seq(residues=custom_list), expected_seq)

    def test_multiple_chain_breaks(self):
        """Test sequence with multiple chain breaks and leading dashes."""
        test_residues = [
            {'code': 'ALA', 'id': 'A_2', 'type': 'atom', 'atoms': []},
            {'code': 'CYS', 'id': 'A_3', 'type': 'atom', 'atoms': []},
            {'code': 'ASP', 'id': 'B_3', 'type': 'atom', 'atoms': []},
            {'code': 'GLU', 'id': 'B_4', 'type': 'atom', 'atoms': []},
            {'code': 'PHE', 'id': 'C_2', 'type': 'atom', 'atoms': []},
        ]
        protein = MockProtein(test_residues)
        expected_seq = "-AC|--DE|-F" # Leading dashes for each chain + chain breaks
        self.assertEqual(protein.seq(), expected_seq)

    # Example of testing for expected errors (e.g., non-standard residue codes)
    def test_non_amino_acid_code_raises_error(self):
        """Test that a KeyError occurs if a code is not in AA_MAP."""
        test_residues = [
            {'code': 'ALA', 'id': 'A_1', 'type': 'atom', 'atoms': []},
            {'code': 'HOH', 'id': 'A_2', 'type': 'hetatm', 'atoms': []}, # Water
            {'code': 'CYS', 'id': 'A_3', 'type': 'atom', 'atoms': []},
        ]
        protein = MockProtein(test_residues)
        # Assumes AA_MAP does not contain 'HOH' and seq() doesn't handle it
        with self.assertRaises(KeyError):
            protein.seq()

if __name__ == '__main__':
    unittest.main() 