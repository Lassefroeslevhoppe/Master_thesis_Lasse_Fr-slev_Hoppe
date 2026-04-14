import os
import sys
import gzip
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np
from statsmodels.nonparametric.smoothers_lowess import lowess

###############################################
# Function for calling out taxon and directories that I will use in the panel plot:
###############################################
if len(sys.argv) != 2:
    sys.exit("Usage: python panel_plot_2.py <TAXON>")

TAXON = sys.argv[1]
###############################################
# Directories
###############################################

LCA_SS = "/maps/projects/prohaska/people/vpt968/Data_analysis_revised_data/Results_of_three_data_types/Updated_mapping/lca_files_with_age/_metadmg_SS"
LCA_DS = "/maps/projects/prohaska/people/vpt968/Data_analysis_revised_data/Results_of_three_data_types/Updated_mapping/lca_files_with_age/metadmg_DS"

DFIT_SS = "/maps/projects/prohaska/people/vpt968/Data_analysis_revised_data/Original_files/Data_SS_lower_threshold/results/metadmg/dfit"
DFIT_DS = "/maps/projects/prohaska/people/vpt968/Data_analysis_revised_data/Original_files/Data_DS_lower_treshold/results/metadmg/dfit"

FASTQ_FILE = "/maps/projects/prohaska/people/vpt968/Data_analysis_revised_data/fastq_read_counts.xlsx"

###############################################
# And output
###############################################
OUTFIG = f"panel_{TAXON}_LOWESS.png"

###############################################
# Here the functions used in my script are defined
# First a function to extract the CGG ID - this is used to link dfit, fastq n_reads and lca files  .  
###############################################

def extract_cgg_id(filename):
    import re
    m = re.search(r"(CGG-\d-\d{6})", filename)
    return m.group(1) if m else None

###############################################
# Here I define the function to read the lca files even though they are compressed 
###############################################

def load_lca_file(path):
    return pd.read_csv(path, sep="\t", compression="gzip", low_memory=False)

###############################################
# A function to find any occurence of the genus called out in the start of the script 
# inside the column taxa_path of the given file .
###############################################

def taxon_in_path(taxa_path, taxon_name):
    if not isinstance(taxa_path, str):
        return False
    return taxon_name.lower() in taxa_path.lower()

###############################################
#  A function to count the taxa_path column where the called taxon is present
# The taxa_path is divided by ";", so I split the string per ";" and look for the occurence
# "genus" and the called taxon.
###############################################
def get_taxid_counts(df_lca, taxon_name):
    taxid_counts = defaultdict(int)
    taxon_lower = taxon_name.lower()

    for _, row in df_lca.iterrows():
        path = row["taxa_path"]
        if not isinstance(path, str):
            continue
        if taxon_lower not in path.lower():
            continue

        for seg in path.split(";"):
            if taxon_lower in seg.lower() and '"genus"' in seg.lower():
                try:
                    tid = seg.split(":")[0].strip()
                    taxid_counts[tid] += 1
                except:
                    pass
    return taxid_counts if taxid_counts else None

###############################################
# Function for loading the dfit files
###############################################

def load_dfit(path):
    try:
        return pd.read_csv(path, sep="\t", compression="gzip")
    except:
        return None

###############################################
# Function to add a confidence interval around the points
###############################################


def lowess_ci(x, y, frac=0.5):
    sm = lowess(y, x, frac=frac, return_sorted=True)
    fit = np.interp(x, sm[:, 0], sm[:, 1])
    resid = np.abs(y - fit)
    sigma = np.nanstd(resid)
    lower = sm[:, 1] - 1.96 * sigma
    upper = sm[:, 1] + 1.96 * sigma
    return sm, lower, upper

###############################################
# Here I read the excel file containing the reads for fastq files 
# and group them according to the data type (SS or DS) and sample ID and
# extract average n_reads per sample pair (there was an occurence of an R1 and R2 
# not having the same value) 
###############################################
fastq = pd.read_excel(FASTQ_FILE)

fastq["Data_type"] = fastq["Data_type"].astype(str).str.strip().str.upper()

fastq_sum = (
    fastq
    .groupby(["sample_id", "Data_type"])["n_reads"]
    .mean()
    .reset_index()
)

fastq_sum.columns = ["CGG", "DataType", "FASTQreads"]


###############################################
# Looping through all the .lca.gz files for the LCA_SS and LCA_DS.
# As a backup, just to get the script running, I added a part where the script would 
# continue even though no CGG ID was found. Further I extract the age per sample 
# (That is found in the column I added to each LCA file earlier) and round it down to two 
# decimals. I further count all the lines (mapped DNA-strings) for each LCA-file for the
# normalization later. Lastly the "taxids =" counts the number of rows per LCA file that 
# the called taxon occurs in. 
###############################################
rows = []

for dtype, folder in [("SS", LCA_SS), ("DS", LCA_DS)]:
    for f in sorted(os.listdir(folder)):
        if not f.endswith(".lca.gz"):
            continue

        cgg = extract_cgg_id(f)
        if not cgg:
            continue

        df = load_lca_file(os.path.join(folder, f))
        age = round(float(df["Age (ka) Brandon et al. 2015"].iloc[0]), 2)
        total_lca = len(df)
        tax_reads = df["taxa_path"].apply(lambda x: taxon_in_path(x, TAXON)).sum()
        taxids = get_taxid_counts(df, TAXON)

        rows.append({
            "CGG": cgg,
            "DataType": dtype,
            "Age": age,
            "TotalLCA": total_lca,
            "TaxReads": tax_reads,
            "TaxIDs": taxids
        })

###############################################
# And then I add the fastq reads 
###############################################

df = pd.DataFrame(rows)
df = df.merge(fastq_sum, on=["CGG", "DataType"], how="left").fillna(0)


###############################################
# Now, for the damage parameter A. As I am gathering reads for a genus that might
# hold many species, each with their individual damage parameter A-value, I take 
# the median of A across all occurences.
###############################################
A_vals = []

for _, r in df.iterrows():
    if not isinstance(r["TaxIDs"], dict):
        A_vals.append(None)
        continue

    dfit_dir = DFIT_SS if r["DataType"] == "SS" else DFIT_DS

    dfit_file = next(
        (os.path.join(dfit_dir, f)
         for f in os.listdir(dfit_dir)
         if r["CGG"] in f and f.endswith(".dfit.gz")),
        None
    )

    if dfit_file is None:
        A_vals.append(None)
        continue

    dfit = load_dfit(dfit_file)
    if dfit is None:
        A_vals.append(None)
        continue

    A_list = []
    for tid in r["TaxIDs"].keys():
        match = dfit[dfit["taxid"] == int(tid)]
        if not match.empty:
            A_list.append(float(match["A"].iloc[0]))

    if len(A_list) == 0:
        A_vals.append(None)
    else:
        A_vals.append(np.median(A_list))

df["A"] = A_vals


###############################################
# Here I normalize the n_reads of the called genus to both n_lca_reads and n_fastq_reads
###############################################
df["RPM_LCA"] = (df["TaxReads"] / df["TotalLCA"]) * 1e6
df["RPM_FASTQ"] = df.apply(
    lambda r: (r["TaxReads"] / r["FASTQreads"] * 1e6) if r["FASTQreads"] > 0 else 0,
    axis=1
)

###############################################
# To plot the samples in the correct order and categorical (so same distance between)
# the samples in the final plot, I build a list that treat the ages as categories, not
# continious values.  
###############################################
ages = sorted(df["Age"].unique())
df["AgeCat"] = pd.Categorical(df["Age"], categories=ages, ordered=True)
ypos = df["AgeCat"].cat.codes

###############################################
# Here I define that I wish to create a 1x4 grid in which each of my four sub-figures 
# go into. Also, the shared title is set
###############################################
fig, axs = plt.subplots(1, 4, figsize=(26, 10))
fig.suptitle(f"Reads and damage patterns of {TAXON}", fontsize=20)

bar_h = 0.35
offset = bar_h / 2
###############################################
# Now, the first three sub-plots are quite similar, so they go together here:
###############################################

for i, (col, title, xlabel) in enumerate([
    ("TaxReads", "Raw reads", "Reads"),
    ("RPM_LCA", "RPM (LCA)", "RPM"),
    ("RPM_FASTQ", "RPM (FASTQ)", "RPM")
]):
    for dt, color, sign in [("SS", "blue", -1), ("DS", "orange", 1)]:
        sub = df[df["DataType"] == dt]
        axs[i].barh(
            sub["AgeCat"].cat.codes + sign * offset,
            sub[col],
            height=bar_h,
            color=color,
            label=dt if i == 0 else None
        )

    axs[i].set_title(title)
    axs[i].set_xlabel(xlabel)
    axs[i].set_yticks(range(len(ages)))
    axs[i].set_yticklabels([f"{a:.2f}" for a in ages])
    axs[i].invert_yaxis()
    axs[i].set_ylabel("Age (ka)")

axs[0].legend()

###############################################
# Panel four has more nuance to it, so I designated a seperate part for that.
# I added bubbles instead of dots in the scatter plot. The largest dot is defined as
# the sample/age with most reads and the others represent the proportionate difference 
# between samples - but each sample gets a minimum bubble size of 30 (if the genus is present)  
###############################################

ax = axs[3]

for dt, color in [("SS", "blue"), ("DS", "orange")]:
    sub = df[(df["DataType"] == dt) & df["A"].notna()]
    if len(sub) < 3:
        continue

    x = sub["A"].values
    y = sub["AgeCat"].cat.codes
    sizes = sub["TaxReads"] / sub["TaxReads"].max() * 500 + 30

    ax.scatter(x, y, s=sizes, color=color, alpha=0.7, label=dt)

    sm, lo, hi = lowess_ci(x, y)
    ax.plot(sm[:, 0], sm[:, 1], color=color, lw=2)
    ax.fill_between(sm[:, 0], lo, hi, color=color, alpha=0.2)

ax.set_title("Damage (median A)")
ax.set_xlabel("Damage parameter A")
ax.set_xlim(0, 0.8)
ax.set_yticks(range(len(ages)))
ax.set_yticklabels([f"{a:.2f}" for a in ages])
ax.invert_yaxis()
ax.legend()



plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig(OUTFIG, dpi=300)
print(f"\nFigure saved as: {OUTFIG}\n")
plt.show()
