from pathlib import Path

import h5py
import pandas as pd
from pyscenic.rss import regulon_specificity_scores

# define directories
root = Path("/home/sejj745/Scratch/pySCENIC")
output = root / "output"

# load cell metadata
metadata = pd.read_csv(root / "input/pyscenic_scoring_metadata.csv", index_col="cell")

# extract AUCell activity matrix
with h5py.File(output / "pyscenic_scoring_auc.loom", "r") as f:
    cell_ids = [x.decode() if isinstance(x, bytes) else x 
                for x in f["col_attrs"]["CellID"][:]]

    auc = pd.DataFrame.from_records(f["col_attrs"]["RegulonsAUC"][:], index=cell_ids)

auc.index.name = "cell"

print("AUCell dimensions:", auc.shape)

# validate cell identifiers
assert auc.index.is_unique
assert metadata.index.is_unique
assert auc.index.difference(metadata.index).empty
assert metadata.index.difference(auc.index).empty

# align metadata to AUCell matrix
metadata = metadata.loc[auc.index]

assert metadata["scenic_group"].notna().all()

# export AUCell matrix
auc.to_csv(output / "pyscenic_auc_matrix.csv.gz", compression="gzip")

print("AUCell matrix exported.")

# calculate regulon specificity scores
rss = regulon_specificity_scores(auc, metadata["scenic_group"])

print("RSS dimensions:", rss.shape)

# export RSS matrix
rss.T.rename_axis("regulon").to_csv(output / "pyscenic_regulon_RSS.csv")

print("RSS matrix exported.")

# calculate RSS across mature neuronal populations
mature_groups = ["PVN: non-NEN", 
                 "ARC: non-NEN", 
                 "PVN: OXT",
                 "PVN: AVP-OXT", 
                 "PVN: AVP-CRH", 
                 "PVN: TRH",
                 "ARC: GHRH"]

mature_metadata = metadata[metadata["scenic_group"].isin(mature_groups)]

mature_auc = auc.loc[mature_metadata.index]

# define NEN and non-NEN populations
nen_groups = {"PVN: OXT": "PVN NEN",
              "PVN: AVP-OXT": "PVN NEN",
              "PVN: AVP-CRH": "PVN NEN",
              "PVN: TRH": "PVN NEN",
              "PVN: non-NEN": "PVN non-NEN",
              "ARC: GHRH": "ARC NEN",
              "ARC: non-NEN": "ARC non-NEN"}

# prepare NEN metadata and AUCell matrix
nen_metadata = metadata[metadata["scenic_group"].isin(nen_groups)].copy()
assert set(nen_metadata["scenic_group"]) == set(mature_groups)
nen_metadata["NEN_group"] = nen_metadata["scenic_group"].map(nen_groups)

nen_auc = auc.loc[nen_metadata.index]

# check population sizes
print(nen_metadata["NEN_group"].value_counts())

# calculate RSS using all cells
nen_rss = regulon_specificity_scores(nen_auc, nen_metadata["NEN_group"])

nen_rss.T.rename_axis("regulon").to_csv(output / "pyscenic_NEN_RSS.csv")

# calculate balanced RSS
n_cells = nen_metadata["NEN_group"].value_counts().min()

balanced_rss = []

for seed in range(20):

    sampled_cells = (nen_metadata
                     .groupby("NEN_group")
                     .sample(n=n_cells, random_state=seed)
                     .index)

    rss_sample = regulon_specificity_scores(nen_auc.loc[sampled_cells], nen_metadata.loc[sampled_cells, "NEN_group"])

    balanced_rss.append(rss_sample)

# average RSS across iterations
balanced_nen_rss = sum(balanced_rss) / len(balanced_rss)

balanced_nen_rss.T.rename_axis("regulon").to_csv(output / "pyscenic_NEN_RSS_balanced.csv")

print("Balanced NEN RSS dimensions:", balanced_nen_rss.shape)

# check population sizes
print(mature_metadata["scenic_group"].value_counts())

# calculate RSS using all mature neurons
mature_rss = regulon_specificity_scores(mature_auc, mature_metadata["scenic_group"])

mature_rss.T.rename_axis("regulon").to_csv(output / "pyscenic_mature_RSS.csv")

print("Mature RSS dimensions:", mature_rss.shape)

# calculate balanced mature-neuron RSS
n_cells = mature_metadata["scenic_group"].value_counts().min()

balanced_rss = []

for seed in range(20):

    sampled_cells = (mature_metadata
        .groupby("scenic_group")
        .sample(n=n_cells, random_state=seed)
        .index)

    rss_sample = regulon_specificity_scores(mature_auc.loc[sampled_cells],
                                            mature_metadata.loc[sampled_cells, "scenic_group"])

    balanced_rss.append(rss_sample)

# average RSS across iterations
balanced_mature_rss = sum(balanced_rss) / len(balanced_rss)

balanced_mature_rss.T.rename_axis("regulon").to_csv(output / "pyscenic_mature_RSS_balanced.csv")

print("Balanced mature RSS dimensions:", balanced_mature_rss.shape)