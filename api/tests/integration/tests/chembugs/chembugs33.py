import os
import sys

##YQ: REMOVE ME BEFORE TESTING. This is just for IDE highlighting.
#from api.tests.integration.common.env_indigo import *

sys.path.append(
    os.path.normpath(
        os.path.join(os.path.abspath(__file__), "..", "..", "..", "common")
    )
)
from env_indigo import *
from sapiopycommons.chem.IndigoMolecules import indigo
indigo.setOption("rpe-bypass-saturation-check", True)
def test_chembugs33():
    print("*** Test CHEMBUGS-33 ***")
    with open(joinPathPy("33/reaction_alt.rxn", "r")) as f:
        query_reaction_rxn: str = f.read()
        query_reaction = indigo.loadQueryReaction(query_reaction_rxn)
        query_reaction.dearomatize()
        query_reaction.aromatize()
        print("Reaction SMARTS: " + query_reaction.smarts())
    with open(joinPathPy("33/r1.mol", "r")) as f:
        mol_file: str = f.read()
        mol1 = indigo.loadMolecule(mol_file)
        mol1.dearomatize()
        mol1.aromatize()
        print("Reactant 1: " + mol1.smiles())
    with open(joinPathPy("33/r2.mol", "r")) as f:
        mol_file: str = f.read()
        mol2 = indigo.loadMolecule(mol_file)
        mol2.dearomatize()
        mol2.aromatize()
        print("Reactant 2: " + mol2.smiles())

    monomers_table = indigo.createArray()
    monomers_table.arrayAdd(indigo.createArray())
    monomers_table.at(0).arrayAdd(mol1)
    monomers_table.arrayAdd(indigo.createArray())
    monomers_table.at(1).arrayAdd(mol2)

    output_reactions = indigo.reactionProductEnumerate(query_reaction, monomers_table)
    output_list = []

    for enum_reaction in output_reactions.iterateArray():
        output_list.append(enum_reaction)

    if len(output_list) == 1:
        print("CHEMBUGS-33 OK")
    else:
        raise ValueError("CHEMBUGS-33 Failed.")