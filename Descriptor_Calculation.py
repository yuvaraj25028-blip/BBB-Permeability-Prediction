import pandas as pd
import numpy as np
from rdkit import Chem
from tqdm import tqdm
import warnings
import urllib.request


warnings.filterwarnings('ignore')


print(" Downloading and loading the BBBP dataset...")
url = "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/BBBP.csv"
file_path = "BBBP.csv"
urllib.request.urlretrieve(url, file_path)

df = pd.read_csv(file_path)

smiles_col = 'smiles' if 'smiles' in df.columns else df.columns[2]

print(f"✅ Loaded {len(df)} drugs. Starting calculations...\n")


def calculate_14_descriptors(smiles):
    try:

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        n = mol.GetNumAtoms()
        m = mol.GetNumBonds()


        if n < 2:
            return None


        A = Chem.GetAdjacencyMatrix(mol, useBO=False).astype(float)


        dist_mat = Chem.GetDistanceMatrix(mol)


        su = np.sum(dist_mat, axis=1)


        degrees = np.sum(A, axis=1)
        D = np.diag(degrees)

        L = D - A
        Q = D + A

        ES = np.zeros((n, n))
        RS = np.zeros((n, n))
        SS = np.zeros((n, n))
        GAS = np.zeros((n, n))

        for bond in mol.GetBonds():
            i = bond.GetBeginAtomIdx()
            j = bond.GetEndAtomIdx()

            sui = su[i]
            suj = su[j]

            if sui == 0 or suj == 0:
                continue

            es_val = 0.5 * ((sui / suj) + (suj / sui))
            ES[i, j] = es_val
            ES[j, i] = es_val

            rs_val = 1.0 / np.sqrt(sui * suj)
            RS[i, j] = rs_val
            RS[j, i] = rs_val

            ss_val = 1.0 / np.sqrt(sui + suj)
            SS[i, j] = ss_val
            SS[j, i] = ss_val

            gas_val = (2.0 * np.sqrt(sui * suj)) / (sui + suj)
            GAS[i, j] = gas_val
            GAS[j, i] = gas_val

        eig_A = np.linalg.eigvalsh(A)
        eig_L = np.linalg.eigvalsh(L)
        eig_Q = np.linalg.eigvalsh(Q)
        eig_ES = np.linalg.eigvalsh(ES)
        eig_RS = np.linalg.eigvalsh(RS)
        eig_SS = np.linalg.eigvalsh(SS)
        eig_GAS = np.linalg.eigvalsh(GAS)

        def get_energy(eig):
            return np.sum(np.abs(eig))

        def get_L_energy(eig):
            avg_degree = (2 * m) / n
            return np.sum(np.abs(eig - avg_degree))

        def get_radius(eig):
            return np.max(np.abs(eig))

        E_A = get_energy(eig_A)
        Rho_A = get_radius(eig_A)

        E_L = get_L_energy(eig_L)
        Rho_L = get_radius(eig_L)

        E_Q = get_L_energy(eig_Q)
        Rho_Q = get_radius(eig_Q)

        E_ES = get_energy(eig_ES)
        Rho_ES = get_radius(eig_ES)

        E_RS = get_energy(eig_RS)
        Rho_RS = get_radius(eig_RS)

        E_SS = get_energy(eig_SS)
        Rho_SS = get_radius(eig_SS)

        E_GAS = get_energy(eig_GAS)
        Rho_GAS = get_radius(eig_GAS)

        return {
            'E_A': E_A, 'Rho_A': Rho_A,
            'E_L': E_L, 'Rho_L': Rho_L,
            'E_Q': E_Q, 'Rho_Q': Rho_Q,
            'E_ES': E_ES, 'Rho_ES': Rho_ES,
            'E_RS': E_RS, 'Rho_RS': Rho_RS,
            'E_SS': E_SS, 'Rho_SS': Rho_SS,
            'E_GAS': E_GAS, 'Rho_GAS': Rho_GAS
        }

    except Exception as e:
        return None

results = []

for idx, row in tqdm(df.iterrows(), total=len(df), desc="Processing Molecules", ncols=100):
    desc = calculate_14_descriptors(row[smiles_col])
    results.append(desc)

results = [r for r in results if r is not None]

df_desc = pd.DataFrame(results)

df_final = pd.concat([df, df_desc], axis=1)

initial_count = len(df_final)
df_final = df_final.dropna()
final_count = len(df_final)

output_file = 'BBBP_14_Descriptors_Fresh.csv'
df_final.to_csv(output_file, index=False)
