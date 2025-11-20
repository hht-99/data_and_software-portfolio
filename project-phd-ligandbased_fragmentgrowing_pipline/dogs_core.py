#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# dogs_core.py
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors, Draw
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import rdChemReactions
import pandas as pd
import base64
from io import BytesIO
import logging
import os

# Optional logging control
logger = logging.getLogger("dogs_core")
logger.setLevel(logging.INFO)

def configure_logging(enable_rdkit_logs=False):
    if not enable_rdkit_logs:
        RDLogger.DisableLog('rdApp.*')
    else:
        RDLogger.EnableLog('rdApp.*')

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

# SMARTS-based FGA/FGI reactions with proper atom mapping
REACTION_SMARTS = [
    ("Acyl chloride (neutral)", "[C:1](=O)[OX2H1:2]>>[C:1](=O)[Cl:2]"),
    ("Acyl chloride (charged)", "[C:1](=O)[O-:2]>>[C:1](=O)[Cl:2]"),
    ("Alcohol to bromide", "[C:1][OX2H:2]>>[C:1][Br:2]"),
    ("Alcohol to chloride", "[C:1][OX2H:2]>>[C:1][Cl:2]"),
    ("Sulfonic acid to sulfonyl chloride (neutral)", "[$([#16X4:1](=[OX1:2])(=[OX1:3])([#6:4])[OX2H,OX1H0-:5]),$([#16X4+2:1]([OX1-:2])([OX1-:3])([#6:4])[OX2H,OX1H0-:5])]>>[S:1](=[OX1:2])(=[OX1:3])([#6:4])[Cl:5]"),
    ("Sulfonic acid to sulfonyl chloride (charged)", "[$([#16X4:1](=[OX1:2])(=[OX1:3])([#6:4])[OX2H0:5]),$([#16X4+2:1]([OX1-:2])([OX1-:3])([#6:4])[OX2H0:5])]>>[S:1](=[OX1:2])(=[OX1:3])([#6:4])[Cl:5]"),
    ("Aryl halide to nitrile", "[c:1][Cl,Br:2]>>[c:1][C:2]#N"),
    ("Alcohol to nitrile", "[C:1][OX2H:2]>>[C:1][C:2]#N"),
    ("Primary amine to nitrile (neutral)", "[C:1][C:2][NX3;H2;!$(NC=[!#6]);!$(NC#[!#6])]>>[C:1][C:2]#N"),
    ("Primary amine to nitrile (charged)", "[C:1][C:2][N+1H3;!$(NC=[!#6]);!$(NC#[!#6])]>>[C:1][C:2]#N"),
    ("Alkyne to nitrile", "[C:1][C:2]#C>>[C:1][C:2]#N"),
]
REACTIONS = [(name, rdChemReactions.ReactionFromSmarts(smarts)) for name, smarts in REACTION_SMARTS]

# Coupling reactions from DOGS Table S1 (minimal SMARTS for reactant matching)
# Each entry corresponds to a coupling reaction with its required input substructure patterns
# These are simplified representations of required substructures for each coupling reaction partner
COUPLING_REACTIONS = [
    # Format: (reaction_name, [reactant1_smiles, reactant2_smiles, ...])
    
    ("Bischler-Napieralski", ["c1cc(CCNC(=O)C)ccc1"]),
    ("Pictet-Gams", ["c1cc(C(O)CNC(=O)C)ccc1"]),
    ("Pictet-Spengler (neutral amine)", ["c1cc(CCN)ccc1", "CC(=O)"]),
    ("Pictet-Spengler (charged amine)", ["c1cc(CC[NH3+])ccc1", "CC(=O)"]),
    ("Bischler Indole", ["c1c(N)cccc1", "c1c(C(=O)CBr)cccc1"]),
    ("Benzimidazol (charge 1)", ["c1cccc(N)c1NC", "CC(=O)O"]),
    ("Benzimidazol (charge 2)", ["c1cccc(N)c1NC", "CC(=O)[O-]"]),
    ("Aminothiazol", ["CC(=O)C(Br)C"]),
    ("Benzoxazol (charge 1)", ["c1ccc(N)c(O)c1", "CC(=O)O"]),
    ("Benzoxazol (charge 2)", ["c1ccc(N)c(O)c1", "CC(=O)[O-]"]),
    ("Benzothiazol", ["c1ccc(N)c(S)c1", "c1cc(C(=O))ccc1"]),
    ("Rap-Stoermer", ["c1ccc(O)c(C(=O))c1", "CC(=O)CCl"]),
    ("Niementowski (charge 1)", ["c1ccc(N)c(C(=O)O)c1"]),
    ("Niementowski (charge 2)", ["c1ccc(N)c(C(=O)[O-])c1"]),
    ("Quinazolinone (1)", ["c1ccc(N)c(C(=O)O)c1", "CN"]),
    ("Quinazolinone (2)", ["c1ccc(N)c(C(=O)[O-])c1", "CN"]),
    ("Quinazolinone (3)", ["c1ccc(N)c(C(=O)O)c1", "C[NH3+]"]),
    ("Quinazolinone (4)", ["c1ccc(N)c(C(=O)[O-])c1", "C[NH3+]"]),
    ("Chinolin-2-one intramol.", ["c1cc(NC(=O)CC(=O)C)ccc1"]),
    ("Tetrazol", ["CC#N"]),
    ("Tetrahydro-Indole (neutral amine)", ["CN", "CC(=O)C(O)C"]),
    ("Tetrahydro-Indole (charged amine)", ["C[NH3+]", "CC(=O)C(O)C"]),
    ("3-nitrile pyridine (symmetry 1)", ["CC(=O)CC(=O)C"]),
    ("3-nitrile pyridine (symmetry 2)", ["CC(=O)CC(=O)C"]),
    ("Triazole", ["c1ccccc1C#N", "NNC(=O)c1ccccc1"]),
    ("Huisgen 1-3 dipolar (azid in_situ)", ["CCl", "CC#C"]),
    ("Diels-Alder (symmetry 1)", ["C=CC=C", "C=C"]),
    ("Diels-Alder (symmetry 2)", ["C=CC=C", "C=C"]),
    ("Diels-Alder Alkyne (1)", ["C=CC=C", "C#C"]),
    ("Diels-Alder Alkyne (2)", ["C=CC=C", "C#C"]),
    ("Spiro-piperidine", ["c1cc(C(=O)C)c(O)cc1", "C1C(=O)CCNC1"]),
    ("Pyrazole (1)", ["CC(=O)CC(=O)C", "NNC"]),
    ("Pyrazole (2)", ["CC(=O)CC(=O)C", "NNC"]),
    ("Phthalazinone (1)", ["c1c(C(=O)O)c(C(=O)C)ccc1", "NNC"]),
    ("Phthalazinone (2)", ["c1c(C(=O)[O-])c(C(=O)C)ccc1", "NNC"]),
    ("Paal-Knorr pyrrole (neutral amine)", ["CC(=O)CCC(=O)C", "CN"]),
    ("Paal-Knorr pyrrole (charged amine)", ["CC(=O)CCC(=O)C", "C[NH3+]"]),
    ("Triaryl-imidazol (1,2 diketone)", ["c1ccccc1C(=O)C(=O)c1ccccc1", "c1ccc(C(=O))cc1"]),
    ("Triarylimidazol (alpha hydroxy-ketone)", ["c1ccccc1C(O)C(=O)c1ccccc1", "c1ccc(C(=O))cc1"]),
    ("Fischer indole", ["c1ccc(NN)cc1", "CC(=O)CC"]),
    ("Friedlaender chinoline", ["c1cc(C=O)c(N)cc1", "CC(=O)CC"]),
    ("Pechmann coumarine", ["c1cc(O)ccc1", "CC(=O)CC(=O)OCC"]),
    ("Benzofuran", ["c1cc(O)c(I)cc1", "CC#C"]),
    ("Imidazol-Acetamid", ["CC(=O)C(Br)"]),
    ("Dieckmann 5-ring (1)", ["CCOC(=O)CCCCC(=O)OCC"]),
    ("Dieckmann 5-ring (2)", ["CCOC(=O)CCCCC(=O)OCC"]),
    ("Dieckmann 6-ring (1)", ["CCOC(=O)CCCCCC(=O)OCC"]),
    ("Dieckmann 6-ring (2)", ["CCOC(=O)CCCCCC(=O)OCC"]),
    ("Flavone", ["c1cc(O)c(C(=O)C)cc1", "c1ccc(C(=O)Cl)cc1"]),
    ("Oxadiazole (charge 1)", ["c1cc(C#N)ccc1", "CC(=O)O"]),
    ("Oxadiazole (charge 2)", ["c1cc(C#N)ccc1", "CC(=O)[O-]"]),
    ("Michael addition", ["CC(=O)CC(=O)C", "C=CC(=O)C"]),
    ("crossed Claissen", ["c1ccccc1C(=O)OC", "CCC(=O)OC"]),
    ("Williamson ether", ["c1cc(O)ccc1", "CCBr"]),
    ("reductive amination ketone (neutral amine)", ["CC(=O)C", "CN"]),
    ("reductive amination ketone (charged amine)", ["CC(=O)C", "C[NH3+]"]),
    ("Suzuki", ["c1cc(B(O)(O))ccc1", "c1cc(Cl)ccc1"]),
    ("Negishi", ["CI", "CCBr"]),
    ("Mitsunobu (imide)", ["CC(O)C", "C(=O)NC(=O)"]),
    ("Mitsunobu (carboxylic acid, neutral)", ["CC(O)C", "CC(=O)O"]),
    ("Mitsunobu (carboxylic acid, charged)", ["CC(O)C", "CC(=O)[O-]"]),
    ("Mitsunobu (sulfonic amide)", ["CC(O)C", "CNS(=O)(=O)C"]),
    ("Heck", ["CBr", "CC(=CC)C"]),
    ("Amide (neutral amine)", ["CC(=O)Cl", "CN"]),
    ("Amide (charged amine)", ["CC(=O)Cl", "C[NH3+]"]),
    ("Ester", ["CC(=O)Cl", "CO"]),
    ("Thioether", ["c1ccccc1C=C", "CS"]),
    ("Ketone", ["CC(=O)Cl", "CI"]),
    ("Sulfonamid (neutral)", ["CS(=O)(=O)Cl", "CN"]),
    ("Sulfonamid (charged)", ["CS(=O)(=O)Cl", "C[NH3+]"]),
    ("Ar-Pyrazole", ["c1cc(B(O)(O))ccc1", "C1=CC=NN1"]),
    ("Ar-Imidazole", ["c1cc(B(O)(O))ccc1", "C1=CN=CN1"]),
    ("Alkyne alkylation", ["CCl", "C#C"]),
    ("Alkyne acylation", ["CC(=O)Cl", "C#C"])
]

# Converts each SMILES in the list to RDKit Mol objects for matching
COUPLING_PATTERNS = {
    name: [Chem.MolFromSmiles(sm) for sm in smiles]
    for name, smiles in COUPLING_REACTIONS
}

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

def mol_to_image_tag(mol):
    try:
        if mol.GetNumAtoms() == 0:
            raise ValueError("Empty molecule")

        AllChem.Compute2DCoords(mol)
        img = Draw.MolToImage(mol, size=(200, 150))
        buffered = BytesIO()
        img.save(buffered, format="PNG")
        img_str = base64.b64encode(buffered.getvalue()).decode("utf-8")
        return f"<img src='data:image/png;base64,{img_str}'/>"
    except Exception as e:
        logger.warning(f"Failed to render molecule: {e}")
        return "Rendering failed"

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
            filtered.append((mol, smi, "original", "none"))
        except:
            continue
    return filtered

def expand_fragments_with_reactions(filtered_tuples):
    all_smiles = set(smi for _, smi, _, _ in filtered_tuples)
    results = filtered_tuples[:]

    for mol, _, _, _ in filtered_tuples:
        for rxn_name, rxn in REACTIONS:
            try:
                products = rxn.RunReactants((mol,))
                for prod_tuple in products:
                    prod = prod_tuple[0]
                    Chem.SanitizeMol(prod)
                    smi = Chem.MolToSmiles(prod)
                    if smi not in all_smiles:
                        all_smiles.add(smi)
                        results.append((prod, smi, "reaction", rxn_name))
            except:
                continue
    return results

def annotate_reactive_sites(mols):
    annotated = []
    for mol in mols:
        matches = []
        for idx, (name, rxn) in enumerate(REACTIONS):
            try:
                reactant_pattern = rxn.GetReactants()[0]
                if mol.HasSubstructMatch(reactant_pattern):
                    if len(mol.GetSubstructMatches(reactant_pattern)) == 1:
                        matches.append(name)
            except:
                continue
        if matches:
            annotated.append((mol, matches))
    return annotated

def annotate_coupling_reactivity(mol):
    matches = []
    for rxn_name, patterns in COUPLING_PATTERNS.items():
        for patt in patterns:
            hits = mol.GetSubstructMatches(patt)
            if len(hits) == 1:
                matches.append(rxn_name)
                break
    return matches

# ---------- Pipeline Entry Point ----------

def prepare_fragment_library(sdf_path, enable_logging=False, output_csv=True):
    if not enable_logging:
        RDLogger.DisableLog("rdApp.*")

    raw_mols, _ = load_fragments_from_sdf(sdf_path)
    filtered = filter_fragments(raw_mols)
    expanded = expand_fragments_with_reactions(filtered)

    if output_csv:
        sdf_name = os.path.splitext(os.path.basename(sdf_path))[0]
        output_filename = f"fragment_library_processed_{sdf_name}.csv"
        
        records = []
        for mol, smi, origin, reaction in expanded:
            image_tag = mol_to_image_tag(mol)
            coupling_matches = annotate_coupling_reactivity(mol)
            records.append({
                "SMILES": smi,
                "Origin": origin,
                "Reaction": reaction,
                "Coupling_Reactivity": ";".join(coupling_matches),
                "Source SDF": sdf_path,
                "Structure": image_tag
            })
        df = pd.DataFrame(records)
        df.to_csv(output_filename, index=False)

    annotated = annotate_reactive_sites([mol for mol, _, _, _ in expanded])
    return annotated

