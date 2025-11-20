#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# dogs_core.py
import csv
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import rdChemReactions

# Allowed atoms
ALLOWED_ELEMENTS = {"C", "N", "O", "S", "P", "F", "Cl", "Br", "I", "B", "Si", "Se"}

# Hann filter SMARTS (problematic groups)
HANN_SMARTS = [
    "[N;R0;$(N=*)]",         # nitroso
    "[N;R0;$(N#N)]",         # azide
    "[N;R0;$(N-[N]=N)]",     # diazo
    "[N;R0;$(N-[C]=[N])]",   # isocyanide
    "[N;R0;$(N=O)]",         # nitro
    "[N;R0;$(N-[O])]",       # hydroxylamine
    "[Cl,Br,I][Cl,Br,I]",    # dihalides
    "[O;R0;$(O-O)]",         # peroxide
]
hann_patterns = [Chem.MolFromSmarts(smarts) for smarts in HANN_SMARTS]

# SMARTS-based FGA/FGI reactions
REACTION_SMARTS = [
    "[C:1](=O)[OX2H1:2]>>[C:1](=O)Cl",
    "[C:1](=O)[O-]>>[C:1](=O)Cl",
    "[C:1][OX2H:2]>>[C:1]Br",
    "[C:1][OX2H:2]>>[C:1]Cl",
    "[$([#16X4](=[OX1])(=[OX1])([#6])[OX2H,OX1H0-]),$([#16X4+2]([OX1-])([OX1-])([#6])[OX2H,OX1H0-])]>>S(=O)(=O)Cl",
    "[$([#16X4](=[OX1])(=[OX1])([#6])[OX2H0]),$([#16X4+2]([OX1-])([OX1-])([#6])[OX2H0])]>>S(=O)(=O)Cl",
    "[c:1][Cl,Br]>>[c:1]C#N",
    "[C:1][OX2H:2]>>[C:1]C#N",
    "[C:1][NX3;H2,H1;!$(NC=O):2]>>[C:1]C#N",
    "[C:1][NX3+1;H3;!$(NC=O):2]>>[C:1]C#N",
    "[C:1][$([CX2]#C):2]>>[C:1]C#N",
]
REACTIONS = [rdChemReactions.ReactionFromSmarts(s) for s in REACTION_SMARTS]

# ---------- Core Functions ----------

def load_fragments_from_sdf(path):
    suppl = Chem.SDMolSupplier(path)
    mols = [mol for mol in suppl if mol is not None]
    smiles = [Chem.MolToSmiles(mol, isomericSmiles=True) for mol in mols]
    return mols, smiles

def is_element_allowed(mol):
    return all(atom.GetSymbol() in ALLOWED_ELEMENTS for atom in mol.GetAtoms())

def has_incorrect_valence(mol):
    try:
        Chem.SanitizeMol(mol)
        return False
    except:
        return True

def has_unwanted_substructure(mol):
    return any(mol.HasSubstructMatch(p) for p in hann_patterns)

def standardize_protonation(mol):
    uncharger = rdMolStandardize.Uncharger()
    mol = uncharger.uncharge(mol)
    Chem.SanitizeMol(mol)
    return mol

def filter_fragments(mols):
    filtered = []
    seen = set()

    for mol in mols:
        if mol is None:
            continue
        try:
            mw = Descriptors.ExactMolWt(mol)
            if not (30 <= mw <= 300):
                continue
            if rdMolDescriptors.CalcNumRings(mol) > 4:
                continue
            if not is_element_allowed(mol):
                continue
            if sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == "F") > 3:
                continue
            if has_incorrect_valence(mol):
                continue
            if has_unwanted_substructure(mol):
                continue
            mol = standardize_protonation(mol)
            smi = Chem.MolToSmiles(mol, isomericSmiles=True)
            if smi in seen:
                continue
            seen.add(smi)
            filtered.append((mol, "original"))
        except:
            continue
    return filtered

def expand_fragments_with_reactions(mol_entries):
    all_smiles = set()
    result_entries = []

    for mol, origin in mol_entries:
        try:
            Chem.SanitizeMol(mol)
            smi = Chem.MolToSmiles(mol)
            if smi not in all_smiles:
                all_smiles.add(smi)
                result_entries.append((mol, origin))
        except:
            continue

        for i, rxn in enumerate(REACTIONS):
            try:
                products = rxn.RunReactants((mol,))
                for prod_tuple in products:
                    prod = prod_tuple[0]
                    Chem.SanitizeMol(prod)
                    smi = Chem.MolToSmiles(prod)
                    label = f"FGA/FGI_{i+1}"
                    if smi not in all_smiles:
                        all_smiles.add(smi)
                        result_entries.append((prod, label))
            except:
                continue

    return result_entries

def annotate_reactive_sites(mol_entries):
    annotated = []
    for mol, label in mol_entries:
        matches = []
        for idx, rxn in enumerate(REACTIONS):
            try:
                reactant_pattern = rxn.GetReactants()[0]
                if mol.HasSubstructMatch(reactant_pattern):
                    if len(mol.GetSubstructMatches(reactant_pattern)) == 1:
                        matches.append(idx)
            except:
                continue
        annotated.append((mol, label, matches))
    return annotated

def write_output_csv(mol_entries, path="processed_fragments.csv"):
    with open(path, mode="w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["SMILES", "Origin"])
        for mol, label in mol_entries:
            try:
                smi = Chem.MolToSmiles(mol)
                writer.writerow([smi, label])
            except:
                continue

# ---------- Pipeline Entry Point ----------

def prepare_fragment_library(sdf_path, output_csv=True):
    raw_mols, _ = load_fragments_from_sdf(sdf_path)
    filtered = filter_fragments(raw_mols)
    expanded = expand_fragments_with_reactions(filtered)
    annotated = annotate_reactive_sites(expanded)
    if output_csv:
        write_output_csv(expanded)
    return annotated

