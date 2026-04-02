import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import NearestNeighbors

print("Loading data...")
# Load the dataset
df = pd.read_csv("unified_variant_feature_table.txt", sep="\t")

# Coerce columns used for calculation to numeric types; invalid parsing will be set as NaN, then drop missing values
cols_to_convert = ['maf', 'ld_proxy', 'dist_tss']
for col in cols_to_convert:
    df[col] = pd.to_numeric(df[col], errors='coerce')
df = df.dropna(subset=cols_to_convert)

# 1. Log-transform data with long-tailed distributions
df['ld_proxy_log'] = np.log1p(df['ld_proxy'])
df['dist_tss_log'] = np.log1p(np.abs(df['dist_tss']))

# 2. Extract features for matching and ensure they are float type in the DataFrame
features = ['maf', 'ld_proxy_log', 'dist_tss_log']
df[features] = df[features].astype(float)

# 3. Z-score normalization (standardize scales to ensure equal weighting for all features)
scaler = StandardScaler()
df[features] = scaler.fit_transform(df[features])

# 4. Separate datasets by variant type (reset index for direct index-based retrieval later)
sv_df = df[df['variant_type'] == 'SV'].reset_index(drop=True)
indel_df = df[df['variant_type'] == 'INDEL'].reset_index(drop=True)
snp_df = df[df['variant_type'] == 'SNP'].reset_index(drop=True)

print(f"Total SVs: {len(sv_df)}, INDELs: {len(indel_df)}, SNPs: {len(snp_df)}")
print("Matching 5 INDELs and 40 SNPs per SV (With Replacement)...")

# ================= Core matching section =================
# Use KD-Tree for efficient KNN search. This is "with replacement" by default 
# because it merely queries the tree and does not remove items from the original pool.

# Build INDEL tree and find nearest neighbors
if len(indel_df) >= 5:
    nn_indel = NearestNeighbors(n_neighbors=5, metric='euclidean', algorithm='kd_tree')
    nn_indel.fit(indel_df[features].values)
    # Batch calculate indices of the nearest 5 INDELs for all SVs
    _, indices_indel = nn_indel.kneighbors(sv_df[features].values)
else:
    indices_indel = []

# Build SNP tree and find nearest neighbors (set to 40)
if len(snp_df) >= 40:
    nn_snp = NearestNeighbors(n_neighbors=40, metric='euclidean', algorithm='kd_tree')
    nn_snp.fit(snp_df[features].values)
    # Batch calculate indices of the nearest 40 SNPs for all SVs
    _, indices_snp = nn_snp.kneighbors(sv_df[features].values)
else:
    indices_snp = []

# ================= Result aggregation section =================
matched_records = []

# Extract ID arrays to speed up the iteration process
sv_ids = sv_df['variant_id'].values
indel_ids = indel_df['variant_id'].values
snp_ids = snp_df['variant_id'].values

for i, sv_id in enumerate(sv_ids):
    # Append the 5 matched INDELs
    if len(indices_indel) > 0:
        for idx in indices_indel[i]:
            matched_records.append({
                'sv_id': sv_id, 
                'matched_variant': indel_ids[idx], 
                'type': 'INDEL'
            })
            
    # Append the 40 matched SNPs
    if len(indices_snp) > 0:
        for idx in indices_snp[i]:
            matched_records.append({
                'sv_id': sv_id, 
                'matched_variant': snp_ids[idx], 
                'type': 'SNP'
            })

# 6. Save output results
out_df = pd.DataFrame(matched_records)
out_df.to_csv("matched_variants_bootstrap.txt", sep="\t", index=False)
print("Matching complete. Output saved to 'matched_variants_bootstrap.txt'.")