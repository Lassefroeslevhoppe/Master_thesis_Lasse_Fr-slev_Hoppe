import pandas as pd
import numpy as np
import gzip
import os
import re
from collections import defaultdict
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from typing import Optional
import matplotlib.patches as mpatches

try:
    from adjustText import adjust_text
    HAS_ADJUSTTEXT = True
except ImportError:
    HAS_ADJUSTTEXT = False

###############################################
# Directories
###############################################
LCA_SS = "/maps/projects/prohaska/people/vpt968/Data_analysis_revised_data/Results_of_three_data_types/Updated_mapping/lca_files_with_age/_metadmg_SS"
LCA_DS = "/maps/projects/prohaska/people/vpt968/Data_analysis_revised_data/Results_of_three_data_types/Updated_mapping/lca_files_with_age/metadmg_DS"
OUTDIR = "/maps/projects/prohaska/people/vpt968/Scripts_for_GitHub/Results"
GENUS_LIST_EXCEL = "/maps/projects/prohaska/people/vpt968/Scripts_for_GitHub/Genera_passing_damage_QC.xlsx"
os.makedirs(OUTDIR, exist_ok=True)

###############################################
# As a double check, I only allow genera from the whitelist to be included
# if they have more than 50 reads, the only allowed rank from the LCA-files
# are "genus" and I define that I wish to extract the top 25 genera from the 
# PCA's here so that I can change it easily if needing to rerun-      
###############################################
ALLOWED_RANKS = {"genus"}
MIN_TOTAL_READS = 50
TOP_N = 25

###############################################
# For the Biplot, I define that I wish to plot the top ten drivers for each 
# PCA biplot as arrows, showing the eigenvector directions.       
###############################################

TOP_N_ARROWS = 10
ARROW_SCALE_FALLBACK = 3.0
FIGSIZE_BIPLOT = (11, 9)

AGE_COL = "Age (ka) Brandon et al. 2015"

###############################################
# Here I define groups that I wish to colour on the final figure (From the excel sheet):
###############################################
Cnidaria = ["Porites", "Acropora", "Turbinaria", "Montipora", "Palythoa", "Leptoseris", "Millepora",
            "Echinopora", "Isopora", "Pocillopora", "Stylophora", "Pachyseris", "Madracis",
            "Duncanopsammia", "Micromussa"]
Green_algae = ["Chloropicon", "Bathycoccus", "Micromonas", "Halimeda", "Ostreobium", "Pycnococcus"]
Sharks = ["Cetorhinus", "Isurus", "Carcharodon", "Sphyrna", "Mustelus"]
Ascidian = ["Diplosoma", "Trididemnum"]
Fish = ["Benthosema", "Caesio", "Engraulis", "Plotosus"]
Sponges = ["Cliona", "Spheciospongia", "Haliclona"]
Mangroves = ["Rhizophora", "Ceriops", "Bruguiera", "Avicennia", "Excoecaria", "Heritiera"]
Terrestrial_plant = ["Solanum", "Ficus", "Eucalyptus", "Themeda", "Calamus", "Corymbia", "Hibiscus",
                     "Acacia", "Urtica", "Casuarina", "Oryza"]
Coastal_plant = ["Phragmites", "Carex"]
Marine_plants = ["Halophila"]
Nematode = ["Cylicocyclus"]
Mollusca = ["Octopus", "Pinctada", "Eledone", "Acanthosepion", "Mactra", "Sepioteuthis"]
Diatomes = ["Pseudo-nitzschia", "Chaetoceros", "Skeletonema"]
Crustaceans = ["Caprella"]
Haptophyta = ["Emiliania"]
Micro_algae = ["Bigelowiella"]
Fungi = ["Sporisorium", "Annulohypoxylon", "Paraphaeosphaeria", "Hypomontagnella", "Periconia", "Xylaria"]

GROUPS_TO_COLOUR = {
    "Cnidaria":        (Cnidaria,         "#C4685D"),
    "Green algae":     (Green_algae,      "#07FE38"),
    "Sharks":          (Sharks,           "#DE33D3"),
    "Ascidians":       (Ascidian,         "#C4E50B"),
    "Fish":            (Fish,             "#206CDF"),
    "Sponges":         (Sponges,          "#E0AB0D"),
    "Mangroves":       (Mangroves,        "#39888A"),
    "Land plants":     (Terrestrial_plant,"#5AA35C"),
    "Coastal plants":  (Coastal_plant,    "#0CF304"),
    "Nematodes":       (Nematode,         "#824545"),
    "Mollusca":        (Mollusca,         "#D10BF0"),
    "Diatoms":         (Diatomes,         "#08CAE3"),
    "Marine plants":   (Marine_plants,    "#4EF60C"),
    "Haptophyta":      (Haptophyta,       "#0124E8"),
    "Crustaceans":     (Crustaceans,      "#EE9105"),
    "Micro algae":     (Micro_algae,      "#8A2BE2"),
    "Fungi":           (Fungi,            "#F70C5A"),
}


GENUS_TO_GROUP = {}
for group_name, (genus_list, color) in GROUPS_TO_COLOUR.items():
    for g in genus_list:
        GENUS_TO_GROUP[g] = (group_name, color)

###############################################
# Here I load the genera from the whitelist. I had some trouble with it because
# I created the sheet manually, and apparently some of the genera differntiating
# formatting or something, making it the script malfunction then looking for them
# in the LCA-files. Also, the genus column was difficult to find, so I had to work
# around that obstacle somehow. This is how I did    
###############################################
def load_genus_whitelist(excel_path: str) -> set:
    df = pd.read_excel(excel_path)

    if df.shape[1] == 0:
        raise ValueError(f"No columns found in {excel_path}")

    norm_map = {}
    for c in df.columns:
        c_str = str(c)
        norm = c_str.strip().lower()
        if norm.endswith(":"):
            norm = norm[:-1].strip()
        norm_map[c] = norm

        genus_col = None
    for original, norm in norm_map.items():
        if norm == "genus":
            genus_col = original
            break

    whitelist = set(
        df[genus_col]
        .dropna()
        .astype(str)
        .str.strip()
    )
    whitelist = {g for g in whitelist if g}

    return whitelist

###############################################
# Here define a function to extract the CGG IDs
###############################################
def extract_sample_id(filename: str) -> Optional[str]:
    match = re.search(r"(CGG[-_]?\d{1,2}[-_]?\d{6})", filename)
    if not match:
        return None
    return match.group(1).replace("-", "_")

###############################################
# Here deifine a function to extract second taxonomic rank from the taxa_path string.
# Later I filter to make sure that it is only genus-level that is counted. The function
# looped through for both the SS and DS LCA-files   
###############################################

def extract_taxon_and_rank(taxa_path):
    if pd.isna(taxa_path):
        return None, None

    for e in str(taxa_path).split(";"):
        parts = e.split(":")
        if len(parts) < 3:
            continue

        name = parts[1].strip('"')
        rank = parts[2].strip('"')

        if rank == "genus":
            return name, rank

    return None, None

###############################################
# Here I define a function to loop through the LCA-files and get the age, 
# sample ID grouped for data types (either SS or DS libraries). I multiply
# the age with 1000 to get it in years instead of ka (some guy on Reddit did that
# because it apparently had some influence on how the PCA worked). I extract genera
# from the taxa_path as defined in the ALLOWED RANK above and return for each sample
# a file with the number of counts for each genus, the sample ages and the number of reads
# in that sample. 
###############################################

def load_lca_folder(folder, data_type):
    taxon_counts = {}
    sample_ages = {}
    sample_total_reads = {}
    for fname in os.listdir(folder):
        if not fname.endswith(".lca.gz"):
            continue
        sample_id = extract_sample_id(fname)
        if sample_id is None:
            continue
        key = (sample_id, data_type)
        with gzip.open(os.path.join(folder, fname), "rt") as f:
            df = pd.read_csv(f, sep="\t", low_memory=False)
        if AGE_COL not in df.columns:
            continue
        age_series = pd.to_numeric(df[AGE_COL], errors="coerce").dropna()
        if age_series.empty:
            continue
        sample_ages[key] = float(age_series.iloc[0]) * 1000.0
        sample_total_reads[key] = int(len(df))
        if "taxa_path" not in df.columns:
            continue
        df[["taxon", "rank"]] = df["taxa_path"].apply(lambda x: pd.Series(extract_taxon_and_rank(x)))
        taxon_counts.setdefault(key, defaultdict(int))

        for _, row in df.iterrows():
            t = row["taxon"]
            r = row["rank"]
            if t and (r in ALLOWED_RANKS):
                taxon_counts[key][t] += 1

    return taxon_counts, sample_ages, sample_total_reads


ss_counts, ss_ages, ss_totals = load_lca_folder(LCA_SS, "SS")
ds_counts, ds_ages, ds_totals = load_lca_folder(LCA_DS, "DS")
all_keys = set(ss_counts) | set(ds_counts)

###############################################
# Here I build the data matrix that is required to run PCA
###############################################

all_taxa = set()
for k in all_keys:
    all_taxa |= set(ss_counts.get(k, {}))
    all_taxa |= set(ds_counts.get(k, {}))

rows = []
for sample_id, dtype in sorted(all_keys):
    key = (sample_id, dtype)
    row = {
        "sample_id": sample_id,
        "Data_type": dtype,
        "age_bp": ss_ages.get(key, ds_ages.get(key)),
        "total_reads": ss_totals.get(key, ds_totals.get(key)),
    }
    for t in all_taxa:
        row[t] = ss_counts.get(key, {}).get(t, 0) + ds_counts.get(key, {}).get(t, 0)
    rows.append(row)

raw = pd.DataFrame(rows)

###############################################
# All genera that are not in the whitelist are excluded by the following part, including 
# the 50 reads bar - if there by any chance were a mixup in one of my previous scripts, 
# there is a safety mechanisms to ensure that only genera with 50 reads or more are 
# included. For the Marine only version as here, I wish to exclude the Octopus reads so
# I do that in the last part here. These raw counts are saved as a .csv for a manual check  
###############################################

meta_cols = {"sample_id", "Data_type", "age_bp", "total_reads"}
taxa_cols = [c for c in raw.columns if c not in meta_cols]
keep_taxa = [t for t in taxa_cols if t in GENUS_WHITELIST]
raw = raw[list(meta_cols) + keep_taxa]

taxa_cols = [c for c in raw.columns if c not in meta_cols]
low_abundance = raw[taxa_cols].sum()[raw[taxa_cols].sum() < MIN_TOTAL_READS].index
raw = raw.drop(columns=low_abundance)
taxa_cols = [c for c in raw.columns if c not in meta_cols]

raw.to_csv(os.path.join(OUTDIR, "filtered_raw_counts_whitelist.csv"), index=False)

###############################################
# Here I normalize the data to RPM by dividing the n_reads of a genus by the
# total number of reads in the corresponding LCA-file and multiplying with one
# million. This is done on a sample based level for each entry in the matrix "raw"
# The output is further saved as a .cvs for manual inspection. 
###############################################

raw["total_reads"] = pd.to_numeric(raw["total_reads"], errors="coerce")
for t in taxa_cols:
    raw[t] = np.where(raw["total_reads"] > 0, raw[t] / raw["total_reads"] * 1e6, 0.0)

rpm = raw[["sample_id", "Data_type", "age_bp"] + taxa_cols]
rpm.to_csv(os.path.join(OUTDIR, "rpm_normalized_to_LCA_rowcounts.csv"), index=False)

###############################################
# Here I use a hellinger transformation on the rpm from the step above.
###############################################
prop = rpm[taxa_cols].div(rpm[taxa_cols].sum(axis=1), axis=0)
prop = prop.replace([np.inf, -np.inf], np.nan).fillna(0.0)
hell = np.sqrt(prop)
hell_df = pd.DataFrame(hell, columns=taxa_cols, index=rpm.index)

###############################################
# Here create a functiomn to apply the PCA function from the sklearn.decomposition imported at the
# top of the script. And a function to run it.
###############################################
def do_pca_with_loadings(X_df: pd.DataFrame, n_components: int = 2):
    pca = PCA(n_components=n_components)
    scores = pca.fit_transform(X_df.values)
    pc_cols = [f"PC{i+1}" for i in range(n_components)]
    scores_df = pd.DataFrame(scores, columns=pc_cols, index=X_df.index)
    loadings_df = pd.DataFrame(pca.components_.T, index=X_df.columns, columns=pc_cols)
    return scores_df, loadings_df, pca

def run_pca(X_df: pd.DataFrame, meta_df: pd.DataFrame):
    scores_df, loadings_df, pca_model = do_pca_with_loadings(X_df, n_components=2)
    return meta_df.join(scores_df), loadings_df, pca_model

###############################################
# Here I run the PCA's for both data types
###############################################

meta = rpm[["sample_id", "Data_type", "age_bp"]].copy()
meta.index = rpm.index

df_hell_pooled, load_hell_pooled, pca_hell_pooled = run_pca(hell_df, meta)
df_hell_ss, load_hell_ss, pca_hell_ss = run_pca(hell_df.loc[meta["Data_type"] == "SS"], meta.loc[meta["Data_type"] == "SS"])
df_hell_ds, load_hell_ds, pca_hell_ds = run_pca(hell_df.loc[meta["Data_type"] == "DS"], meta.loc[meta["Data_type"] == "DS"])

###############################################
# For adding the arrows to the PCA plots and making them bi-plots, I need to add
# arrows - here for the top 10 genera. The following function scales the lenght of
# the arrows automatically so that they represent to what degree they influence the
# positioning in the PCA plots. To keep the arrows within the figure, I defined that
# they should only be represented by 35% of the computed lenght.      
###############################################

def autoscale_arrow_scale(scores_df: pd.DataFrame,
                          loadings_df: pd.DataFrame,
                          top_n: int = TOP_N_ARROWS,
                          fallback: float = ARROW_SCALE_FALLBACK) -> float:
    score_radius = max(
        np.nanmax(np.abs(scores_df["PC1"].values)),
        np.nanmax(np.abs(scores_df["PC2"].values))
    )
    mag = np.sqrt(loadings_df["PC1"]**2 + loadings_df["PC2"]**2)
    top_mag = float(mag.sort_values(ascending=False).head(top_n).max())
    if top_mag <= 0 or not np.isfinite(top_mag) or score_radius <= 0 or not np.isfinite(score_radius):
        return fallback
    return 0.35 * score_radius / top_mag

###############################################
#  The DS biplot and the SS biplot should use only data corresponding to the data type
#  The explained variance is multiplied by 100 to get the value in percent, not a
#  decimal number. And the functon I defined above is used in this function for the 
#  scaling of the arrows.Further it is specified that only the top 10 drivers 
#  should be included and they are coloured according to the groups and colour definitions
#  set at the top of the script. I further add a legend.          
###############################################

def plot_pca_biplot_separate(df_scores_meta: pd.DataFrame,
                             pca_model: PCA,
                             loadings_df: pd.DataFrame,
                             title: str,
                             outpath: str,
                             top_n_arrows: int = TOP_N_ARROWS,
                             arrow_scale: Optional[float] = None):

    pc1_pct = pca_model.explained_variance_ratio_[0] * 100.0
    pc2_pct = pca_model.explained_variance_ratio_[1] * 100.0

    mags = np.sqrt(loadings_df["PC1"]**2 + loadings_df["PC2"]**2)
    top_taxa = mags.sort_values(ascending=False).head(top_n_arrows).index
    top_load = loadings_df.loc[top_taxa]

    sub = df_scores_meta.sort_values("age_bp")
    if sub.empty:
        return

    arrow_scale_use = (
        autoscale_arrow_scale(sub[["PC1", "PC2"]], loadings_df, top_n=top_n_arrows, fallback=ARROW_SCALE_FALLBACK)
        if arrow_scale is None else arrow_scale
    )

    fig, ax = plt.subplots(figsize=FIGSIZE_BIPLOT)
    ax.plot(sub["PC1"], sub["PC2"], color="black", linewidth=1.2, alpha=0.6, zorder=1)

    sc = ax.scatter(
        sub["PC1"], sub["PC2"],
        c=sub["age_bp"], cmap="viridis",
        edgecolor="black", linewidth=0.4, s=130, zorder=3
    )
    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label("Age (BP)", fontsize=16)
    cbar.ax.tick_params(labelsize=16)

    head_w = 0.03 * arrow_scale_use
    head_l = 0.03 * arrow_scale_use
    texts = []

    groups_seen = {} 

    for taxon, row in top_load.iterrows():
        x = row["PC1"] * arrow_scale_use
        y = row["PC2"] * arrow_scale_use

        if taxon in GENUS_TO_GROUP:
            group_name, color = GENUS_TO_GROUP[taxon]
            groups_seen[group_name] = color
        else:
            group_name, color = "Other", "black"

        ax.arrow(
            0, 0, x, y,
            head_width=head_w, head_length=head_l,
            fc=color, ec=color,
            lw=1.7, alpha=0.80,
            length_includes_head=True,
            zorder=2
        )
        texts.append(
            ax.text(x * 1.08, y * 1.08, str(taxon),
                    fontsize=16, color=color, fontweight="bold", zorder=4)
        )

    if HAS_ADJUSTTEXT and len(texts) > 1:
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="gray", lw=0.4, alpha=0.7))

    handles = [mpatches.Patch(color=c, label=g) for g, c in sorted(groups_seen.items())]
    if handles:
        ax.legend(handles=handles, loc="upper right", frameon=False, fontsize=12)

    ax.axhline(0, color="gray", linestyle="--", alpha=0.5)
    ax.axvline(0, color="gray", linestyle="--", alpha=0.5)

    ax.set_xlabel(f"PC1 ({pc1_pct:.1f}%)", fontsize=16)
    ax.set_ylabel(f"PC2 ({pc2_pct:.1f}%)", fontsize=16)
    ax.set_title(f"{title} (Top {top_n_arrows} loadings)", fontsize=16)
    ax.tick_params(axis="both", which="major", labelsize=16)

    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    plt.close(fig)

###############################################
# Here the biplots are finished by using the function made above, one for SS libraries and 
# one for DS libraries 
###############################################

plot_pca_biplot_separate(
    df_hell_ss, pca_hell_ss, load_hell_ss,
    title="PCA biplot – Hellinger – SS",
    outpath=os.path.join(OUTDIR, f"PCA_biplot_Hellinger_SS_top{TOP_N_ARROWS}.png"),
    top_n_arrows=TOP_N_ARROWS
)

plot_pca_biplot_separate(
    df_hell_ds, pca_hell_ds, load_hell_ds,
    title="PCA biplot – Hellinger – DS",
    outpath=os.path.join(OUTDIR, f"PCA_biplot_Hellinger_DS_top{TOP_N_ARROWS}.png"),
    top_n_arrows=TOP_N_ARROWS
)

