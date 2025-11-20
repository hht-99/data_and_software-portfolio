#!/usr/bin/env python
# coding: utf-8

# dogs_growth.py (modular fragment growth engine)
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit.Chem.rdMolAlign import GetBestRMS
from copy import deepcopy
import heapq

# --- Data Classes ---
class FragmentRecord:
    def __init__(self, mol, smiles, origin, reactions):
        self.mol = mol
        self.smiles = smiles
        self.origin = origin  # 'original' or FGA/FGI name
        self.reactions = reactions  # list of applicable coupling reactions

class ReactionStep:
    def __init__(self, reaction_name, partner_fragment_smiles, resulting_smiles):
        self.reaction_name = reaction_name
        self.partner_fragment_smiles = partner_fragment_smiles
        self.resulting_smiles = resulting_smiles

class ProductTrail:
    def __init__(self, mol, history, mass, score):
        self.mol = mol
        self.history = history  # list of ReactionStep
        self.mass = mass
        self.score = score

    def __lt__(self, other):
        return self.score > other.score  # max-heap by score

# --- Shape scoring function ---
def compute_shape_similarity(ref_mol, query_mol):
    ref = deepcopy(ref_mol)
    query = deepcopy(query_mol)
    try:
        AllChem.EmbedMolecule(ref, AllChem.ETKDG())
        AllChem.EmbedMolecule(query, AllChem.ETKDG())
        AllChem.UFFOptimizeMolecule(ref)
        AllChem.UFFOptimizeMolecule(query)
        rms = GetBestRMS(ref, query)
        return 1 / (1 + rms)  # Convert RMSD to similarity-like score
    except:
        return 0.0

# --- Growth controller ---
def grow_from_fragment(fragment, library, ref_mol, max_steps=4, mass_bounds=(0.7, 1.3)):
    trails = []
    queue = [ProductTrail(fragment.mol, [], Descriptors.MolWt(fragment.mol), compute_shape_similarity(ref_mol, fragment.mol))]

    ref_mass = Descriptors.MolWt(ref_mol)
    min_mass = mass_bounds[0] * ref_mass
    max_mass = mass_bounds[1] * ref_mass

    while queue:
        current = queue.pop(0)
        current_mass = current.mass
        current_score = current.score

        if len(current.history) >= max_steps:
            trails.append(current)
            continue

        # Re-annotate current mol with new reactive sites (to implement externally)
        new_reactions = annotate_mol_reactions(current.mol)  # Placeholder
        if not new_reactions:
            trails.append(current)
            continue

        for reaction_name in new_reactions:
            partners = [f for f in library if reaction_name in f.reactions]
            for partner in partners:
                rxn = REACTION_DICT[reaction_name]  # Provided externally
                try:
                    products = rxn.RunReactants((current.mol, partner.mol))
                    for prod_tuple in products:
                        prod = prod_tuple[0]
                        Chem.SanitizeMol(prod)
                        new_mass = Descriptors.MolWt(prod)
                        new_score = compute_shape_similarity(ref_mol, prod)

                        if new_mass > max_mass:
                            trails.append(current)
                            continue

                        new_history = current.history + [ReactionStep(reaction_name, partner.smiles, Chem.MolToSmiles(prod))]

                        if new_mass < min_mass or new_score > current_score:
                            queue.append(ProductTrail(prod, new_history, new_mass, new_score))
                        else:
                            trails.append(current)

                except Exception:
                    continue

    return trails

# --- Batch wrapper for simplified use ---
def grow_from_library(fragments, ref_mol, max_steps=4, mass_bounds=(0.7, 1.3)):
    all_trails = []
    fragment_records = [
        FragmentRecord(mol, Chem.MolToSmiles(mol), "original", reactions)
        for mol, reactions in fragments
    ]

    for fragment in fragment_records:
        trails = grow_from_fragment(
            fragment=fragment,
            library=fragment_records,
            ref_mol=ref_mol,
            max_steps=max_steps,
            mass_bounds=mass_bounds
        )
        all_trails.extend(trails)

    return all_trails
