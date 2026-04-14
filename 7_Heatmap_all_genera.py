import os
import re
import gzip
from pathlib import Path
from typing import Optional, Tuple, List, Dict
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


###############################################
# Directories
###############################################

LCA_SS = "/maps/projects/prohaska/people/vpt968/Data_analysis_revised_data/Results_of_three_data_types/Updated_mapping/lca_files_with_age/_metadmg_SS"
LCA_DS = "/maps/projects/prohaska/people/vpt968/Data_analysis_revised_data/Results_of_three_data_types/Updated_mapping/lca_files_with_age/metadmg_DS"
GENUS_LIST_EXCEL = "/maps/projects/prohaska/people/vpt968/Scripts_for_GitHub/Genera_passing_damage_QC.xlsx"

OUTDIR = "/maps/projects/prohaska/people/vpt968/Scripts_for_GitHub/Results"
os.makedirs(OUTDIR, exist_ok=True)

###############################################
# Settings, finding the correct column in the whitelist, age column and taxa_path. 
###############################################

WHITELIST_COL = "Genus:"
TAXA_PATH_COL = "taxa_path"
AGE_COL = "Age (ka) Brandon et al. 2015"

FIND_GENUS = re.compile(r':"([^"]+)":"([^"]+)"')
EPS_PERCENT = 1e-3


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
Fungi = ["Sporisorium", "Annulohypoxylon", "Paraphaeosphaeria", "Hypomontagnella", "Periconia", "Xylaria", "Neoroussoella"]

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

###############################################
# Function to read the excel whitelist column "genus", removing any blank spaces,
# trailing whitespaces (I made the excel manually, and the script just would not run
# without taking all these measures into account. Thanks to a guy somewhere on Reddit
# I found a way to work around it) and form a standard python list (tolist)  
###############################################

def load_whitelist(excel_path: str, genus_col: str) -> List[str]:
    df = pd.read_excel(excel_path)
    return (
        df[genus_col]
        .dropna()
        .astype(str)
        .str.strip()
        .replace("", np.nan)
        .dropna()
        .unique()
        .tolist()
    )

###############################################
#  A function to extract the sample ID from the LCA-files
###############################################

def extract_sample_id(fname: str) -> str:
    for suf in [".lca.gz"]:
        if fname.endswith(suf):
            return fname[: -len(suf)]
    return fname

###############################################
#  Here I defnine a function to extract "Genus" from the taxa_path of the LCA-files.
#  Th LCA-files taxa-path is divided by ";" and each taxonomic level has a double " 
#  around it and my script did not like that. So I remove them here and strip them of 
#  trailing whitespaces (just as above. Next I split the now uninterupted string by 
#  the ";"s and find the segments that contain "genus"
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
#  Here I define a function to open the compressed lca files and put in a 
#  Pandas dataframe.
###############################################

def read_lca_readlevel_table(path: str) -> pd.DataFrame:
    with gzip.open(path, "rt") as f:
        return pd.read_csv(f, sep="\t", low_memory=False)

###############################################
#  Here I extract ages of the LCA-files and put them in a dictionary string
#  with two empty slots.
#  Next they are ordered to fit the Sample ID's (from the first function in 
#  this script) ann return fill the sample ID's and age labels into the dictionary
###############################################

def make_age_labels(sample_ids: List[str], ages: pd.Series) -> Dict[str, str]:
    ages = ages.reindex(sample_ids)
    labels = [
        f"{a:.2f}" if pd.notna(a) else "NA"
        for a in ages.values
    ]
    return dict(zip(sample_ids, labels))

###############################################
#  Now I define a function, using some of the functions above to
#  create a matrix by gathering all the strings made from each 
#  LCA file and put it into  Tuple strings as they much better hold different
#  data types. I use the functions I defined above to extract sample ID's,
#  to put the lca files into pandas dataframe into df (dataframe). Further I make the
#  age column numeric values and simply use the first valid value (anything other than 0) 
#  for age from each LCA-file. Then I use the previously defined function to extract genera
#  from the taxa_path, creating a string of these. Next I count the number of times each
#  genus occur in said string structured like "Sample_ID, Genus_name, Number_of_occurences"
###############################################

def build_full_genus_matrix(lca_dir: str) -> Tuple[pd.DataFrame, pd.Series]:
    rows, ages = [], {}
    for lp in sorted(Path(lca_dir).glob("*.lca.gz")):
        sample_id = extract_sample_id(lp.name)
        df = read_lca_readlevel_table(str(lp))
        age_vals = pd.to_numeric(df[AGE_COL], errors="coerce").dropna()
        ages[sample_id] = float(age_vals.iloc[0]) if len(age_vals) else np.nan
        genera = (
            df[TAXA_PATH_COL]
            .astype(str)
            .map(extract_genus_from_taxa_path)
            .dropna()
        )
        rows.append({"sample": sample_id, **genera.value_counts().to_dict()})
    mat = pd.DataFrame(rows).set_index("sample").fillna(0).astype(float)
    age_s = pd.Series(ages).reindex(mat.index)
    return mat, age_s

###############################################
#  This function takes the pd.Dataframe (mat) from above and turn the reads per genus
#  into a % of the entire LCA-file. So if there are e.g. 500 reads and one genus has 5,
#  the new value is 1 (%). whatsoever, each row sum to 100 (%). This is simply a way of
#  normalizing the data so taxa that appear across samples can be compared.    
###############################################

def percent_per_sample(mat: pd.DataFrame) -> pd.DataFrame:
    X = mat.to_numpy(float)
    rs = X.sum(axis=1, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        P = np.divide(X, rs, out=np.zeros_like(X), where=rs > 0) * 100.0
    return pd.DataFrame(P, index=mat.index, columns=mat.columns)

###############################################
#  Before plotting, I add a function to reciognize the groups that gets their respective colour.
# The second part here adds the colour.  
###############################################


def colour_to_each_group(ax, groups_to_colour, alpha=0.5, pad=1.5):
    genus_to_colour = {}
    for group_name, (genera, colour) in groups_to_colour.items():
        for g in genera:
            genus_to_colour.setdefault(g, colour)

    for lable in ax.get_xticklabels():
        g = lable.get_text()
        if g in genus_to_colour:
            lable.set_bbox(dict(
                facecolor=genus_to_colour[g],
                edgecolor="none",
                alpha=alpha,
                pad=pad
            ))

###############################################
#  To add a legend, I define the following function to find the correct groups defined in
#  groups to colour and add the colour of them to the legend along the   
###############################################

def legends_for_x_axis_colour(groups_to_colour, shown_genera):
    handles = []
    for group_name, (genera, colour) in groups_to_colour.items():
        if any(g in shown_genera for g in genera):
            handles.append(Patch(facecolor=colour, edgecolor="none", label=group_name, alpha=0.5))
    return handles



###############################################
#  Here I define a function to create the heatmap. Here gahter the functions
#  defined above. I make sure to preserve the order of the whitelist, as they are structured
#  after phyla. The sample order is set and the percentage matix is restructured to match.
#  The Z=-part makes the logtransformation adding the ESP_PERCENT to ensure no "0"-values.
#  This might inflate the very rare genera a little, but not to a degree where it has any real
#  influence. Otherwise, the rest of the function is just to make sure the heatmap itself fits
#  the figure, colouring and so forth.    
###############################################

def plot_heatmap(
    counts: pd.DataFrame,
    ages: pd.Series,
    whitelist: List[str],
    title: str,
    outpath: str
):
    perc = percent_per_sample(counts)
    genus_columns = [g for g in whitelist if g in perc.columns]
    if not genus_columns:
        return

    perc = perc[genus_columns]
    order = ages.sort_values().index.tolist()
    perc = perc.loc[order]
    ages_ord = ages.reindex(order)

    Z = np.log10(perc.to_numpy() + EPS_PERCENT)
    finite = Z[np.isfinite(Z)]
    vmin, vmax = np.percentile(finite, [2, 98])

    fig, ax = plt.subplots(
        figsize=(max(10, 0.16 * Z.shape[1] + 6),
                 max(8, 0.20 * Z.shape[0] + 4))
    )

    x = np.arange(Z.shape[1] + 1)
    y = np.arange(Z.shape[0] + 1)
    mesh = ax.pcolormesh(
        x, y, Z,
        shading="flat",
        edgecolors="lightgrey",
        linewidth=0.35,
        cmap="viridis",
        vmin=vmin,
        vmax=vmax
    )

    ax.invert_yaxis()
    ax.set_title(title)

 
    ax.set_xticks(np.arange(Z.shape[1]) + 0.5)
    ax.set_xticklabels(perc.columns.tolist(), rotation=90, fontsize=16, fontstyle="italic")

  
    colour_to_each_group(ax, GROUPS_TO_COLOUR)


    cbar = fig.colorbar(mesh, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("log10(percentage)")


    legend_handles = legends_for_x_axis_colour(GROUPS_TO_COLOUR, set(perc.columns))
    if legend_handles:
        ax.legend(handles=legend_handles, title="Ecological group", loc="upper left", bbox_to_anchor=(1.1, 1), frameon=False, fontsize=16, title_fontsize=16)

    label_map = make_age_labels(perc.index.tolist(), ages_ord)
    ax.set_yticks(np.arange(Z.shape[0]) + 0.5)
    ax.set_yticklabels([label_map[s] for s in perc.index], fontsize=16)

    ax.set_ylabel("Age (ka)", fontsize=16)

    for spine in ax.spines.values():
        spine.set_visible(False)

    ax.tick_params(length=0)

    fig.subplots_adjust(bottom=0.30)
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)
###############################################
#  Finally i define the main that is called in the bottom to create the figures
###############################################

def main():
    whitelist = load_whitelist(GENUS_LIST_EXCEL, WHITELIST_COL)
    ss_counts, ss_ages = build_full_genus_matrix(LCA_SS)
    ds_counts, ds_ages = build_full_genus_matrix(LCA_DS)

    plot_heatmap(
        ss_counts, ss_ages, whitelist,
        "SS library - genera",
        os.path.join(OUTDIR, "SS_whitelist_log10percent_kapkstyle.png")
    )

    plot_heatmap(
        ds_counts, ds_ages, whitelist,
        "DS library - genera",
        os.path.join(OUTDIR, "DS_whitelist_log10percent_kapkstyle.png")
    )


if __name__ == "__main__":
    main()
