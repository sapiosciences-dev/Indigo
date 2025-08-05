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

def test_chembugs42():
    print("*** Test CHEMBUGS-42 ***")
    with open(joinPathPy("42/reaction.rxn", "r")) as f:
        query_reaction_rxn: str = f.read()
        query_reaction = indigo.loadQueryReaction(query_reaction_rxn)
        query_reaction.aromatize()
        query_reaction.automap("keep")
    with open(joinPathPy("42/r1.mol", "r")) as f:
        mol_file: str = f.read()
        mol1 = indigo.loadMolecule(mol_file)
        mol1.aromatize()
    with open(joinPathPy("42/r2.mol", "r")) as f:
        mol_file: str = f.read()
        mol2 = indigo.loadMolecule(mol_file)
        mol2.aromatize()

    print("Reactant 1: " + mol1.smiles())
    print("Reactant 2: " + mol2.smiles())

    monomers_table = indigo.createArray()
    monomers_table.arrayAdd(indigo.createArray())
    monomers_table.at(0).arrayAdd(mol1)
    monomers_table.arrayAdd(indigo.createArray())
    monomers_table.at(1).arrayAdd(mol2)

    output_reactions = indigo.reactionProductEnumerate(query_reaction, monomers_table)
    output_list = []

    for enum_reaction in output_reactions.iterateArray():
        enum_reaction.aromatize()
        output_list.append(enum_reaction)

    if len(output_list) == 1:
        print("CHEMBUGS-42 OK")
    else:
        raise ValueError("CHEMBUGS-42 Failed.")