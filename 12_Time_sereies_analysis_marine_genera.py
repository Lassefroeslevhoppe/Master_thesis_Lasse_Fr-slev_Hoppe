import os
import re
import gzip
import shutil
import subprocess
from collections import defaultdict
from typing import Optional, List, Tuple
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from sklearn.decomposition import PCA

###############################################
# Directories
###############################################
LCA_SS = "/maps/projects/prohaska/people/vpt968/Data_analysis_revised_data/Results_of_three_data_types/Updated_mapping/lca_files_with_age/_metadmg_SS"
LCA_DS = "/maps/projects/prohaska/people/vpt968/Data_analysis_revised_data/Results_of_three_data_types/Updated_mapping/lca_files_with_age/metadmg_DS"
GENUS_LIST_EXCEL = "/maps/projects/prohaska/people/vpt968/Scripts_for_GitHub/Marine_genera_damage_QC.xlsx"
OUTDIR = "/maps/projects/prohaska/people/vpt968/Scripts_for_GitHub/Results/Time_series_Marine"
os.makedirs(OUTDIR, exist_ok=True)  

###############################################
# Setting for STARS and figure colour bands.
###############################################
MAXIXIMUM_AGE = 19.0
STARS_L = 2
STARS_P = 0.05

FIND_GENUS = re.compile(r':"([^"]+)":"([^"]+)"')

CLR_PSEUDOCOUNT = 1e-6

RSCRIPT_BIN = "Rscript"


AGE_BANDS = [
    ("Late Pleistocene", 11.7, 19.0, "#6A3D9A"),
    ("Early Holocene",    8.2, 11.7, "#FFD92F"),
    ("Mid Holocene",      4.2,  8.2, "#4DAF4A"),
    ("Late Holocene",     1.0,  4.2, "#E5C494"),
    ("Historical",        0.1,  1.0, "#377EB8"),
    ("Recent",            0.0,  0.1, "#A6CEE3"),
]

###############################################
# Here I define a function to load the genera from the whitelist. This function is different from
# the one in the heatmap-script because I have already here secured the order of the genera. What
# happens here is that I create a string List[str] containing all the genera while trimming trailing
# whitespaces, make it all lowercase. It turned out it was easier to do it this way, instead of 
# trying to match the "Genus:" columns name perfectly.       
###############################################

def load_genus_whitelist(excel_path: str) -> List[str]:
    df = pd.read_excel(excel_path)
    genus_col = df.columns[0]
    for c in df.columns:
        if str(c).strip().lower().startswith("genus"):
            genus_col = c
            break
    genera = df[genus_col].dropna().astype(str).str.strip().tolist()

    seen = set()
    out = []
    for g in genera:
        if g and g not in seen:
            seen.add(g)
            out.append(g)
    return out

###############################################
# Now i run the function above and put it in WHITELIST_SET   
###############################################

WHITELIST_GENERA = load_genus_whitelist(GENUS_LIST_EXCEL)
WHITELIST_SET = set(WHITELIST_GENERA)

###############################################
# Here I define a function to extract the CGG ID (and only the CGG ID)   
###############################################

def extract_sample_id(filename: str) -> Optional[str]:
    m = re.search(r"CGG-\d-\d{6}", filename)
    return m.group(0) if m else None

###############################################
#  Here I defnine a function to extract "Genus" from the taxa_path of the LCA-files.
#  Th LCA-files taxa-path is divided by ";" and each taxonomic level has a double " 
#  around it and my script did not like that. So I remove them here and strip them of 
#  trailing whitespaces (just as above. Next I split the now uninterupted string by 
#  the ";"s and find the segments that contain "genus". 
###############################################

def extract_genus_from_taxa_path(taxa_path: str) -> Optional[str]:
    if not isinstance(taxa_path, str) or taxa_path.strip() == "":
        return None
    s = taxa_path.replace('""', '"').strip()
    if s.startswith('"') and s.endswith('"'):
        s = s[1:-1]
    genus = None
    for seg in s.split(";"):
        m = FIND_GENUS.search(seg)
        if m and m.group(2).lower() == "genus":
            genus = m.group(1).strip()
    return genus

###############################################
#  Here I defnine a function to make the CLR transformation that goes into pd.DataFrame.
#  Here I build on the X step by step, first ensuring that everything is numerical, i add
#  pseudocount deined at the top of the script, normalize the rows (so they add up to 1).
#  I make sure that all 0-values are replaced with the pseudocount and ensure that the rows
#  sum to 1. Now that I know there are no zero values (I had som issues with that), i apply
#  the log function and center the data. 
###############################################

def clr_transform(rel_abund_df: pd.DataFrame, pseudocount: float = 1e-6) -> pd.DataFrame:
    X = rel_abund_df.astype(float).copy()
    X = X + pseudocount
    X = X.div(X.sum(axis=1).replace(0, np.nan), axis=0).fillna(0.0)  
    X = X.replace(0, pseudocount)
    X = X.div(X.sum(axis=1), axis=0)

    logX = np.log(X)
    return logX.sub(logX.mean(axis=1), axis=0)

###############################################
#  I define a function to make the shaded bands from the setting at the top of the script.
#  The alpha was set to 0,25 to ensure enough transparancy and zorder to put it behind all 
#  other objects of the figure.  
###############################################

def shade_age_bands(ax):
    for _, lo, hi, color in AGE_BANDS:
        ax.axvspan(lo, hi, color=color, alpha=0.25, zorder=0)

###############################################
# Here I define a function to extract both age, count the total rows in the LCA files,
# the ages and counts of genera from the whitelist and prepare containers for counts, ages
# and totals. The age values are stripped for trailing spaces and mulitplied with 1000 to get
# years. I extract all rows with genus it it and create a string with all genera listed. This
# string is now used to find the occurence of g (a genus from the whitelist) and count them with
# the +=1. The string of counted genera now looks like {"GenusA":100, "GenusB":30,...}.and finally
# returns three outputs as defined for the counts, ages and totals    
###############################################

def load_lca_folder(folder: str, dtype: str, whitelist_set: set):
    counts, ages, totals = {}, {}, {}
    for fn in os.listdir(folder):
        if not fn.endswith(".lca.gz"):
            continue
        sample_id = extract_sample_id(fn)
        if not sample_id:
            continue

        key = (sample_id, dtype)
        path = os.path.join(folder, fn)

        with gzip.open(path, "rt") as f:
            df = pd.read_csv(f, sep="\t", low_memory=False)

        ages[key] = float(df["Age (ka) Brandon et al. 2015"].dropna().iloc[0]) * 1000.0
        totals[key] = int(len(df))

        df["genus"] = df["taxa_path"].apply(extract_genus_from_taxa_path)
        df["genus"] = df["genus"].astype(str).str.strip()
        df.loc[df["genus"].isin(["None", "nan", "NaN"]), "genus"] = np.nan

        cdict = defaultdict(int)
        for g in df["genus"].dropna():
            if g in whitelist_set:
                cdict[g] += 1
        counts[key] = cdict

    return counts, ages, totals


###############################################
# Here I define a function to run PCA on SS or DS seperately, make a 
# CSV for R and return the CSV path and sample ages.   
###############################################

def make_pca_scores_csv_for_dtype(
    clr_df: pd.DataFrame,
    meta: pd.DataFrame,
    dtype: str,
    prefix: str
) -> Tuple[str, np.ndarray]:
    mask = meta["Data_type"].astype(str) == dtype
    ages_ka = meta.loc[mask, "age_bp"].astype(float) / 1000.0
    ages_ka = ages_ka[ages_ka <= MAXIXIMUM_AGE]

    idx = ages_ka.index
    X = clr_df.loc[idx].copy()

    min_n = 2 * STARS_L + 1
    if X.shape[0] < min_n:
        raise RuntimeError(f"[ERROR] Not enough {dtype} samples for STARS: n={X.shape[0]}, need >= {min_n}")

    pca = PCA(n_components=2)
    scores = pca.fit_transform(X.values)

    df_scores = (
        pd.DataFrame({"age_ka": ages_ka.values, "PC1": scores[:, 0], "PC2": scores[:, 1]})
        .sort_values("age_ka", ascending=False)
        .reset_index(drop=True)
    )

    out_csv = os.path.join(OUTDIR, f"{prefix}_PCA_scores_{dtype}.csv")
    df_scores.to_csv(out_csv, index=False)
    return out_csv, df_scores["age_ka"].to_numpy()

###############################################
#  I define a function to help me make sure that R exsist in the environment that
#  I have created to run this script.  It was quite a hassle, but the environement to 
#  use can be found along with the R-script in the supplementary material. Instead of
#  opening R, I have found out that it can be run directly through this script. This is
#  the most complicated part - here I define a function to run the R script directly. 
#  The secure_file_path ensures that the paths are readable by in the R-script.
#  The R-code itself, starts by loading libraries - we need the rshift to preform
#  the STARS analysis, the readr library is used for reading and creating csv files
#  and finally dplyr provides many of the functions used in the rscript.
#  First set the paths to the pca_scores and the out_global containing the ages where
#  significant shifts are found. 
#  Next I define the L (lenght) and P (significance level) as an integer number and as 
#  a float as the defined regime lenght (STARS_L) is a "full" number (e.g. 2) and my
#  significance level is 0.05 (a float). Because I do this in python, the values defined
#  are passed to the R_code before it is run. The "pca" gets the values of the csv file and
#  only keeps rows with a value (is.finite = with a value) for both age, PC1 and 2 and finally
#  I arrange them in a descending age - because I wish to apply the STARS on the oldest to the
#  youngest ages. Just before I run the Rodionov(), I tell the R-script to examine the axis columns
#  (PC1 or PC2) and run them one at a time. I rename them to y (both PC1 and PC2) so they can be
#  compared later in the script. Age is kept as it is, mostly because there is only one column of
#  age and it apllies to both PC1 and PC2. Again, just to make sure that the order I defined before
#  is kept, I arrange the rows after descending age. Now I apply the Rodionov() function (the basis 
#  of the entire STARS analysis) with the constants L (STARS_L) and P (STARS_P) defined at the top 
#  of the python script here. I have done it this way, because I am running the analysis with varying
#  values of STARS_L. I then turn the output of the Rodionov() into a tibble and build the output containing
#  the shifts and at what ages they occur as well as the RSI (how strong they are). Now the "run_axis" is
#  preformed for for both axes (PC1 and PC2) and written to the CSV.    
#  Back in Python, I then define the placement of the R-script (that I just composed) and it is run fromm the
#  terminal with an "error" if I did something wrong. I finally managed to get it to go though, so that is 
#  nice.             
###############################################

def run_rshift(pca_scores_csv: str, out_global_csv: str):

    def secure_file_path(p: str) -> str:
        return p.replace("\\", "/").replace('"', '\\"')

    r_code = f"""
suppressPackageStartupMessages({{
  library(rshift)
  library(readr)
  library(dplyr)
}})

PCA_CSV <- "{secure_file_path(pca_scores_csv)}"
OUT_GLOBAL <- "{secure_file_path(out_global_csv)}"

L <- {int(STARS_L)}
P <- {float(STARS_P)}

pca <- read_csv(PCA_CSV) %>%
  filter(is.finite(age_ka), is.finite(PC1), is.finite(PC2)) %>%
  arrange(desc(age_ka))

run_axis <- function(axis_col) {{
  d <- pca %>% select(age_ka, y = all_of(axis_col)) %>% arrange(desc(age_ka))
  STARS_OUTPUT <- Rodionov(d, col="y", time="age_ka", l=L, prob=P, merge=FALSE)
  STARS_OUTPUT <- as_tibble(STARS_OUTPUT)
  tibble(axis = axis_col,
         shift_age_ka = as.numeric(STARS_OUTPUT[[1]]),
         RSI = as.numeric(STARS_OUTPUT[[2]]))
}}

global <- bind_rows(run_axis("PC1"), run_axis("PC2")) %>%
  arrange(axis, desc(shift_age_ka))

write_csv(global, OUT_GLOBAL)
"""
    r_path = os.path.join(OUTDIR, "run_rshift.R")
    with open(r_path, "w") as f:
        f.write(r_code)

    cp = subprocess.run([RSCRIPT_BIN, r_path], capture_output=True, text=True)
    if cp.returncode != 0:
        raise RuntimeError("Failed")

###############################################
# Here I define a function for plotting. I filter the ages and the shifts so that they
# align when I am plotting them. I create a figure with room for two subplots (2 rows, 1 column)
# Sharex makes them share the same X-axis. I make the plots seperately SS and DS, add the shade bands
# as defined in the function somewhere above and then I plot the sample ages with the fixed Y-value
# of 0.5, just to show where they are. I plot the STARS shifts as vertical lines and some other visual
# editing. I invert the X-axis so it goes from old to younger samples and add a legend and a 
# figure title and save the figures.      
###############################################

def plot_global_shift(shifts_ss: List[float], shifts_ds: List[float], ages_ss: np.ndarray, ages_ds: np.ndarray, out_prefix: str):
    shifts_ss = sorted([a for a in shifts_ss if np.isfinite(a) and a <= MAXIXIMUM_AGE])
    shifts_ds = sorted([a for a in shifts_ds if np.isfinite(a) and a <= MAXIXIMUM_AGE])
    ages_ss = np.array([a for a in ages_ss if np.isfinite(a) and a <= MAXIXIMUM_AGE], dtype=float)
    ages_ds = np.array([a for a in ages_ds if np.isfinite(a) and a <= MAXIXIMUM_AGE], dtype=float)

    fig, axes = plt.subplots(2, 1, figsize=(10.5, 5.2), sharex=True)

    for ax, label, shifts, ages in [(axes[0], "SS", shifts_ss, ages_ss), (axes[1], "DS", shifts_ds, ages_ds)]:
        shade_age_bands(ax)
        if len(ages) > 0:
            ax.scatter(ages, np.full_like(ages, 0.5), s=22, color="black", alpha=0.75, zorder=4)
        for a in shifts:
            ax.axvline(a, color="black", lw=1.4, alpha=0.8, zorder=3)
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.set_ylabel(label, rotation=0, ha="right", va="center", labelpad=18, fontsize=16)
        ax.grid(True, axis="x", alpha=0.15, zorder=1)
        ax.text(0.01, 0.90,
                f"{len(ages)} samples - {len(shifts)} shift(s)",
                transform=ax.transAxes, ha="left", va="top", fontsize=16)

    axes[0].set_xlim(MAXIXIMUM_AGE, 0.0)
    axes[-1].set_xlabel("Age (ka BP)", fontsize=16)
    band_handles = [Patch(facecolor=c, edgecolor="none", alpha=0.25, label=name)
                for name, _, _, c in AGE_BANDS]

    fig.legend(handles=band_handles, loc="upper center", ncol=3, frameon=True, fontsize=16, bbox_to_anchor=(0.5, 0.98))
    fig.tight_layout(rect=[0, 0, 1, 0.8])

    png = os.path.join(OUTDIR, f"{out_prefix}_global_STARS_shifts_SS_vs_DS_0to{MAXIXIMUM_AGE:g}ka.png")
    pdf = os.path.join(OUTDIR, f"{out_prefix}_global_STARS_shifts_SS_vs_DS_0to{MAXIXIMUM_AGE:g}ka.pdf")
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    plt.close(fig)


###############################################
# Here I wrap up the functions I have defined above to create the resulting figures. 
###############################################

ss_counts, ss_ages, ss_totals = load_lca_folder(LCA_SS, "SS", WHITELIST_SET)
ds_counts, ds_ages, ds_totals = load_lca_folder(LCA_DS, "DS", WHITELIST_SET)

all_keys = set(ss_counts) | set(ds_counts)

rows = []
for sample_id, dtype in sorted(all_keys):
    key = (sample_id, dtype)
    row = {
        "sample_id": sample_id,
        "Data_type": dtype,
        "age_bp": ss_ages.get(key, ds_ages.get(key)),
        "total_reads": ss_totals.get(key, ds_totals.get(key)),
    }
    csrc = ds_counts.get(key, {}) if dtype == "DS" else ss_counts.get(key, {})
    for g in WHITELIST_GENERA:
        row[g] = int(csrc.get(g, 0))
    rows.append(row)

raw = pd.DataFrame(rows)

taxa_cols = WHITELIST_GENERA[:]
denom = raw["total_reads"].replace(0, np.nan)

rel = raw[["sample_id", "Data_type", "age_bp"] + taxa_cols].copy()
for g in taxa_cols:
    rel[g] = raw[g] / denom
rel[taxa_cols] = rel[taxa_cols].fillna(0.0)

meta = rel[["sample_id", "Data_type", "age_bp"]].copy()
clr_df = clr_transform(rel[taxa_cols], pseudocount=CLR_PSEUDOCOUNT)
clr_df.index = rel.index
meta.index = rel.index

prefix = "rshift_stars_global_CLR_all_genera"

pca_ss_csv, ages_ss = make_pca_scores_csv_for_dtype(clr_df, meta, "SS", prefix)
pca_ds_csv, ages_ds = make_pca_scores_csv_for_dtype(clr_df, meta, "DS", prefix)

out_ss = os.path.join(OUTDIR, f"{prefix}_stars_global_SS.csv")
out_ds = os.path.join(OUTDIR, f"{prefix}_stars_global_DS.csv")

run_rshift(pca_ss_csv, out_ss)
run_rshift(pca_ds_csv, out_ds)

global_ss_df = pd.read_csv(out_ss)
global_ds_df = pd.read_csv(out_ds)

shifts_ss = sorted(set(global_ss_df["shift_age_ka"].dropna().astype(float).tolist()))
shifts_ds = sorted(set(global_ds_df["shift_age_ka"].dropna().astype(float).tolist()))

plot_global_shift(shifts_ss, shifts_ds, ages_ss, ages_ds, prefix)

