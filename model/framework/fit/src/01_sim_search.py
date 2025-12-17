import os
import pandas as pd
from rdkit import Chem


act = pd.read_csv("../data/actives.csv")
all = pd.read_csv("../data/train.csv")

inact = all[all["bin"]==0]

print(len(all), len(act), len(inact))

active_smi = act["smiles"].tolist()
inactive_smi = inact["st_smiles"].tolist()

from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs

def smiles_to_fp(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)

active_fps = [(s, smiles_to_fp(s)) for s in active_smi]
inactive_fps = [(s, smiles_to_fp(s)) for s in inactive_smi]

matches = []

for s_act, fp_act in active_fps:
    sims = DataStructs.BulkTanimotoSimilarity(
        fp_act, [fp for _, fp in inactive_fps]
    )
    for (s_inact, _), sim in zip(inactive_fps, sims):
        if sim > 0.90:
            matches.append((s_act, s_inact, sim))

print(f"Found {len(matches)} pairs with Tanimoto > 0.90")

for a, i, sim in matches:
    print(a, i, sim)

inactives_drop = {s_inact for _, s_inact, _ in matches}
all_filtered = all[~all["st_smiles"].isin(inactives_drop)].copy()

print("Rows before:", len(all))
print("Rows removed:", len(all) - len(all_filtered))
print("Rows after:", len(all_filtered))
print("Actives final:", len(all_filtered[all_filtered["bin"]==1]))

all_filtered.to_csv("../data/train_filtered.csv", index=False)