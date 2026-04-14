import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from typing import Optional


###############################################
# Directories
###############################################

SS_DIR = "/maps/projects/prohaska/people/vpt968/Data_analysis_revised_data/Results_of_three_data_types/Updated_mapping/lca_files_with_age/_metadmg_SS"
DS_DIR = "/maps/projects/prohaska/people/vpt968/Data_analysis_revised_data/Results_of_three_data_types/Updated_mapping/lca_files_with_age/metadmg_DS"
LIB_and_ENC_DIR = "/maps/projects/prohaska/people/vpt968/Data_analysis_revised_data/Original_files/Data_SS_lower_threshold/results/metadmg/lca"
WHITELIST = "/maps/projects/prohaska/people/vpt968/Data_analysis_revised_data/Genera_passing_damage_QC.xlsx"
OUTDIR = "/maps/projects/prohaska/people/vpt968/Scripts_for_GitHub/Results/Passing_QC"
os.makedirs(OUTDIR, exist_ok=True)

###############################################
# This it the first filter applied. The other ones come from the LIB_and_ENC_DIR and the WHITELIST, and
# I have defined it here because it makes it easier to change the min reads limit later. This value is for 
# reads of a speficic genus across all samples. 
###############################################

MIN_READS_TOTAL = 50

###############################################
# Here I read the LCA-files and pass them to the pd.DataFrame. This is basically the same function
# as from the Structure_finalized script except that I pass it to the pd.DataFrame 
###############################################

def read_lca_file(lca_files: str) -> pd.DataFrame:
    return pd.read_csv(lca_files, sep="\t", compression="gzip" if lca_files.endswith(".gz") else None, dtype=str, comment="#")

###############################################
# Here I reuse the normalization function. 
# I strip the names for aný trailing whitespaces and tabs and make it all lowercase.
###############################################

def normalization_of_names(species_name):
    if pd.isna(species_name):
        return None
    return " ".join(str(species_name).strip().split()).lower()

###############################################
# This function is also reused from the Structure_finalized script:
# Now it's time to get the species extracted from the taxa_path. The species are the highest 
# taxonomix level possibel for the taxa_string, so if present they will always be first. Therefore
# I simply extract the first entry for every line in the LCA's and later filter those with the 
# Structure Score List. So, each taxonomic level in the LCA's are divided by a ; and that's why I
# split it there. Further, each taxonomic level is split (":") and only the center part [1] is used
# as that is the organisms' name at the taxonomic level (species, genus etc.) that cal be called out
# when using this function later in the script. 
###############################################

def extract_rank(taxa_path, rank):
    if pd.isna(taxa_path):
        return None
    for entry in str(taxa_path).split(";"):
        if f':"{rank}"' in entry:
            return entry.split(":")[1].strip('"')
    return None

###############################################
# Here I load the WHITELIST excel sheet and normalize the names by stripping them of blank spaces
# and colons. The column names are cleaned as well 
###############################################

def load_genus_whitelist(excel_path: str) -> set:
    df = pd.read_excel(excel_path)
    df.columns = [str(c).strip().rstrip(":") for c in df.columns]
    return set(df["Genus"].dropna().map(normalization_of_names))


###############################################
# Here I define a function to load the LCA-files and pass them to pandas dataframe. 
###############################################

def load_samples(LCA_files: str, data_type: str) -> pd.DataFrame:
    rows = []
    for file in os.listdir(LCA_files):
        if not file.endswith(".lca.gz"):
            continue
        df = read_lca_file(os.path.join(LCA_files, file))
        df.columns = [str(c).strip() for c in df.columns]
        if "Age (ka) Brandon et al. 2015" not in df.columns or "taxa_path" not in df.columns:
            continue
        age = pd.to_numeric(df["Age (ka) Brandon et al. 2015"], errors="coerce").dropna()
        if age.empty:
            continue
        age_val = round(float(age.iloc[0]), 2)
        genus = df["taxa_path"].apply(lambda x: extract_rank(x, "genus")).dropna().map(normalization_of_names).dropna()
        if genus.empty:
            continue
        rows.append(pd.DataFrame({
            "type": data_type,
            "age": age_val,
            "genus_norm": genus.values
        }))
    return pd.concat(rows, ignore_index=True)

###############################################
# Here I define a function to load the genera from the controls (ENC and LIB-files) and put in an
# empty set. All ENC and LibPTC files are examined and all the sample lca files are skipped. I extract
# all the occurences of "genus" and use the extract rank function from above. That way I get a
# nice and cleaned set that I can use to remove from the LCA-files later.   
###############################################

def load_control_genera(LIB_and_ENC_files: str) -> set:
    present = set()
    for file in os.listdir(LIB_and_ENC_files):
        if not ((file.startswith("ENC") or file.startswith("LibPTC")) and file.endswith("collapsed.lca.gz")):
            continue
        dataframe = read_lca_file(os.path.join(LIB_and_ENC_files, file))
        dataframe.columns = [str(c).strip() for c in dataframe.columns]
        if "taxa_path" not in dataframe.columns:
            continue
        genera = dataframe["taxa_path"].apply(lambda x: extract_rank(x, "genus")).dropna().map(normalization_of_names).dropna().unique()
        present.update(genera.tolist())
    present.discard("")
    present.discard(None)
    return present


###############################################
# Here I load the SS and DS LCA files, and merge the two dataframes ss and ds. The genera from the
# ENCs and LIBs are loaded as well as the whitelist genera.  
###############################################
ss = load_samples(SS_DIR, "SS")
ds = load_samples(DS_DIR, "DS")
samples = pd.concat([ss, ds], ignore_index=True)

CONTROL_GENERA = load_control_genera(LIB_and_ENC_DIR)
WHITELIST = load_genus_whitelist(WHITELIST)

###############################################
# The fist line here group samples by type (SS or DS) and age and finds all the unique 
# normalized genus names per age and type and count them. 
# Next the genus_totals counts how many rows each genus appeas in because I then filter away all
# genera with >= 50 reads, genera present in ENC and LIB files and whitelist genera. After each step
# the number of remaining unique genera per age and datatype is counted. Lastly with the merged
# I create a table with all four steps: before filtering, 50 n_reads, controls removed and whitelist
# aka Passing damage QC   
###############################################
before_filtering = samples.groupby(["type", "age"])["genus_norm"].nunique().reset_index(name="Before filtering")

genus_totals = samples.groupby(["type", "genus_norm"]).size().reset_index(name="total_rows")
keep = genus_totals.loc[genus_totals["total_rows"] >= MIN_READS_TOTAL, ["type", "genus_norm"]]
genera_after_50_nreads_removal = samples.merge(keep, on=["type", "genus_norm"], how="inner")
counted_genera_after_50_nreads_removal = genera_after_50_nreads_removal.groupby(["type", "age"])["genus_norm"].nunique().reset_index(name="≥50 reads")

genera_after_CONTROL_removal = genera_after_50_nreads_removal[~genera_after_50_nreads_removal["genus_norm"].isin(CONTROL_GENERA)]
counted_genera_after_CONTROL_removal = genera_after_CONTROL_removal.groupby(["type", "age"])["genus_norm"].nunique().reset_index(name="Controls removed")

genera_after_whitelist_removal = genera_after_CONTROL_removal[genera_after_CONTROL_removal["genus_norm"].isin(WHITELIST)]
counted_genera_after_whitelist_removal = genera_after_whitelist_removal.groupby(["type", "age"])["genus_norm"].nunique().reset_index(name="Passing damage QC")

merged = (
    before_filtering.merge(counted_genera_after_50_nreads_removal, on=["type", "age"], how="left")
          .merge(counted_genera_after_CONTROL_removal, on=["type", "age"], how="left")
          .merge(counted_genera_after_whitelist_removal, on=["type", "age"], how="left")
          .fillna(0)
)

step_cols = ["Before filtering", "≥50 reads", "Controls removed", "Passing damage QC"]
for c in step_cols:
    merged[c] = merged[c].astype(int)

###############################################
# Here I create an excel file with the genera remaining after the different QC steps.  
###############################################
qc_out = merged.copy()
qc_out["Removed at ≥50 reads"] = qc_out["Before filtering"] - qc_out["≥50 reads"]
qc_out["Removed at Controls removed"] = qc_out["≥50 reads"] - qc_out["Controls removed"]
qc_out["Removed at Passing damage QC"] = qc_out["Controls removed"] - qc_out["Passing damage QC"]

qc_out["% removed at ≥50 reads"] = np.where(qc_out["Before filtering"] > 0,
                                            100 * qc_out["Removed at ≥50 reads"] / qc_out["Before filtering"], np.nan)
qc_out["% removed at Controls removed"] = np.where(qc_out["≥50 reads"] > 0,
                                                   100 * qc_out["Removed at Controls removed"] / qc_out["≥50 reads"], np.nan)
qc_out["% removed at Passing damage QC"] = np.where(qc_out["Controls removed"] > 0,
                                                    100 * qc_out["Removed at Passing damage QC"] / qc_out["Controls removed"], np.nan)

qc_out = qc_out[[
    "type", "age",
    "Before filtering", "≥50 reads", "Controls removed", "Passing damage QC",
    "Removed at ≥50 reads", "Removed at Controls removed", "Removed at Passing damage QC",
    "% removed at ≥50 reads", "% removed at Controls removed", "% removed at Passing damage QC"
]].sort_values(["type", "age"])

loss_cols = ["Removed at ≥50 reads", "Removed at Controls removed", "Removed at Passing damage QC"]
summary = qc_out.groupby("type", as_index=False)[step_cols + loss_cols].sum()

summary["% removed at ≥50 reads"] = np.where(summary["Before filtering"] > 0,
                                            100 * summary["Removed at ≥50 reads"] / summary["Before filtering"], np.nan)
summary["% removed at Controls removed"] = np.where(summary["≥50 reads"] > 0,
                                                    100 * summary["Removed at Controls removed"] / summary["≥50 reads"], np.nan)
summary["% removed at Passing damage QC"] = np.where(summary["Controls removed"] > 0,
                                                     100 * summary["Removed at Passing damage QC"] / summary["Controls removed"], np.nan)

excel_out = os.path.join(OUTDIR, "SS_DS_genus_QC_removed_per_step.xlsx")
with pd.ExcelWriter(excel_out, engine="openpyxl") as xw:
    qc_out.to_excel(xw, sheet_name="Removed_per_age", index=False)
    summary.to_excel(xw, sheet_name="Totals_by_type", index=False)

###############################################
# Here I plot the figure. I get alle the ages for the Y axis and reverse them so that the oldest
# at the top of the plot. I define the colours and how the bars should look. Next I create the
# figure in which the subplots go into: one row, two columns and make them share the y-axis. 
# I define the x-axis and that SS data should be used for SS plot and DS for DS. Age us used to
# group the n_genera for each QC-step for so there now is a dataset for DS and SS. By using reindex
# the same order is sequred and any occurence of zeroes or N/A gets a 0. I define the positioning on
# the y-axis by arranging all ages. Now the bars representing the different QC-steps are defined
# and they are drawn to the plots. Then titles, labels and legends are included.          
###############################################
all_ages = sorted(merged["age"].unique(), reverse=True)

colors = ["#bbb6b6", "#6f8fba", "#366c92", "#042639"]
bar_h = 0.25
offsets = [(-1.5)*bar_h, (-0.5)*bar_h, (0.5)*bar_h, (1.5)*bar_h]

fig, axes = plt.subplots(1, 2, figsize=(14, 10), sharey=True)

for ax, type in zip(axes, ["SS", "DS"]):
    merged_subset = merged[merged["type"] == type].set_index("age").reindex(all_ages).fillna(0)
    y = np.arange(len(all_ages))

    for lab, col, off in zip(step_cols, colors, offsets):
        ax.barh(y + off, merged_subset[lab].values, height=bar_h, color=col, label=lab)

    ax.set_title(type, fontsize=18)
    ax.tick_params(axis="x", labelsize=18)
    ax.set_xlabel("Number of genera", fontsize=18)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{a:.2f}" for a in all_ages], fontsize=18)
    ax.invert_yaxis()
    ax.grid(axis="x", alpha=0.3)

axes[0].set_ylabel("Age (ka)", fontsize=18)

handles, leg_labels = axes[0].get_legend_handles_labels()
fig.legend(handles, leg_labels, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, -0.03), fontsize=18)
fig.suptitle("Genus retention across QC steps", fontsize=20)

plt.tight_layout()
plt.subplots_adjust(bottom=0.12)

png_out = os.path.join(OUTDIR, "SS_DS_genus_QC_barplot.png")
plt.savefig(png_out, dpi=300, bbox_inches="tight")
plt.close()
print("✅ QC bar plot created:", png_out)

