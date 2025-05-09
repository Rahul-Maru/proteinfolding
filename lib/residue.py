"""Residue type definitions for protein structures."""

from typing import TypedDict, List

#TODO do this

class Residue(TypedDict):
    """
    A TypedDict representing a residue in a protein structure.
    
    Attributes:
        code (str): 3-letter residue code (e.g., 'ALA', 'GLY')
        id (str): Chain ID + sequence number (e.g., 'A-42')
        type (str): Record type, typically 'atom' or 'hetatm'
        atoms (List[str]): List of ATOM or HETATM records for this residue
    """
    code: str
    id: str
    type: str
    atoms: List[str]