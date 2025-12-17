import os
import pandas as pd
from rdkit import Chem
from standardiser import standardise
import numpy as np

def standardise_smiles(smiles):
    st_smiles = []
    for smi in smiles:
        if smi is None:
            st_smi = np.nan
            st_smiles += [st_smi]
            continue
        smi = str(smi)
        smi = smi.strip()
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            st_smi = np.nan
            st_smiles += [st_smi]
            continue
        try:
            std_mol = standardise.run(mol)
            st_smi = Chem.MolToSmiles(std_mol, canonical=True)
            st_smiles += [st_smi]
        except:
            st_smi = Chem.MolToSmiles(mol, canonical=True)
            st_smiles += [st_smi]
        if std_mol is None:
            st_smi = Chem.MolToSmiles(mol, canonical=True)
            st_smiles += [st_smi]
            continue
    return st_smiles

root = os.path.dirname(os.path.abspath(__file__))

# Specs compounds obtained from: https://academic.oup.com/database/article/doi/10.1093/database/baaf081/8374761?login=false

df =pd.read_excel(os.path.join(root, "..", "data", "R4A_input.xlsx"), sheet_name="Sheet1")
df["st_smiles"] = standardise_smiles(df["smiles"].tolist())
df = df[~df["st_smiles"].isna()] 
print("Specs:", len(set(df["st_smiles"].tolist())), len(df))

# Prestwick compounds obtained from: https://www.lsi.umich.edu/science/centers-technologies/center-chemical-genomics/services/sample-libraries

df2 = pd.read_csv(os.path.join(root, "..", "data", "Prestwick.csv"))
df2["st_smiles"] = standardise_smiles(df2["Structure"].tolist())
df2 = df2[~df2["st_smiles"].isna()] 
print("Prestwick",len(set(df2["st_smiles"].tolist())), len(df2))

all_smi = list(set(df["st_smiles"])) + list(set(df2["st_smiles"]))
print("Specs + Prestwick by smiles", len(all_smi), len(set(all_smi)))

all_smi = list(set(all_smi))

inchis_all = []
for smi in all_smi:
    mol = Chem.MolFromSmiles(smi)
    inchis_all += [Chem.inchi.MolToInchiKey(mol)]

df_all = pd.DataFrame({"st_smiles": all_smi, "inchikey": inchis_all})
print("Specs and Prestwick libraries together:",df_all.shape)
df_all = df_all.drop_duplicates(subset=["inchikey"], keep="first")
print("Specs and Prestwick fter deduplicating by InchiKey:",df_all.shape)

df3 = pd.read_csv(os.path.join(root, "..", "data", "actives.csv"))
df3["st_smiles"] = standardise_smiles(df3["smiles"].tolist())
df3 = df3[~df3["st_smiles"].isna()] 
active_smi = df3["st_smiles"].tolist()

inchis_active = []
for smi in active_smi:
    mol = Chem.MolFromSmiles(smi)
    inchis_active += [Chem.inchi.MolToInchiKey(mol)]

df3["inchikey"]=inchis_active
print("Actives:", df3.shape)
df3 = df3.drop_duplicates(subset=["inchikey"], keep="first")
print("After deduplicating by InchiKey:",df3.shape)
df3 = df3.drop(columns=["name", "smiles"])

for inchi in inchis_active:
    if inchi not in df_all["inchikey"].tolist():
        print(inchi)


df_joint = pd.concat([df_all, df3])
print("Joint:",df_joint.shape)
df_joint = df_joint.drop_duplicates(subset=["inchikey"],keep="first")
print("Joint dedup:",df_joint.shape)


activity = []
for s in df_joint["inchikey"].tolist():
    if s in inchis_active:
        activity += [1]
    else:
        activity += [0]

df_joint["bin"]=activity

print("Final data")
print(len(df_joint[df_joint["bin"]==1]))
print(df_joint.shape)
df_joint.to_csv(os.path.join(root, "..", "data", "train.csv"), index=False)
