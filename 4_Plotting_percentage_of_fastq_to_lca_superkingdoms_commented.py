import os
import gzip
import pandas as pd
import re
import matplotlib.pyplot as plt

###############################################
# Directories 
###############################################

metadmg_ss_dir = "/maps/projects/prohaska/people/vpt968/Data_analysis_revised_data/Results_of_three_data_types/Updated_mapping/lca_files_with_age/_metadmg_SS"
metadmg_ds_dir = "/maps/projects/prohaska/people/vpt968/Data_analysis_revised_data/Results_of_three_data_types/Updated_mapping/lca_files_with_age/metadmg_DS"

prefilter_ss_dir = "/maps/projects/prohaska/people/vpt968/Data_analysis_revised_data/Original_files/Data_SS_lower_threshold/results/prefilter_metadmg/lca/lca_with_age"
prefilter_ds_dir = "/maps/projects/prohaska/people/vpt968/Data_analysis_revised_data/Original_files/Data_DS_lower_treshold/results/prefilter_metadmg/lca/lca_with_age"

fastq_xlsx = "/maps/projects/prohaska/people/vpt968/Scripts_for_GitHub/fastq_read_counts.xlsx"

output_dir = "/maps/projects/prohaska/people/vpt968/Data_analysis_revised_data/Results"
os.makedirs(output_dir, exist_ok=True)

###############################################
# Superkingdoms and the colors they have in the final plot.
###############################################
CATEGORIES = ["Bacteria", "Archaea", "Eukaryota", "Viruses", "unclassified", "misclassified"]

COLORS = {
    "Bacteria":      "#1f77b4",
    "Archaea":       "#ff7f0e",
    "Eukaryota":     "#2ca02c",
    "Viruses":       "#d62728",
    "unclassified":  "#7f7f7f",
    "misclassified": "#bcbd22",
}

###############################################
# Here I load the FASTQ counts from the fastq_read_counts.xlsx and extract the sample ID's, Data type
# (either SS or DS libraries) and the number of reads after bioinformatics QC. In my excel file, there are
# both R1 and R2 values, so I specifiy that any duplicates should be dropped, as I otherwise overestimate 
# the number of DNA-strings in the FASTQ files. The R1 and R2 are there simply because the strings have been 
# run in both directions - these are used in mapping, but usually either R1 or R2 return a taxon. Lastly the 
# data type is stripped of trailing spaces and renamed to "group"
###############################################

fastq_df = pd.read_excel(fastq_xlsx)
fastq_df = fastq_df[["sample_id", "Data_type", "n_reads_after_QC"]].drop_duplicates()
fastq_df["group"] = fastq_df["Data_type"].str.strip()

###############################################
# First I define a function to open the zipped lca-files. 
###############################################

def open_lca_file(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

###############################################
# Next I define a function to extract the age from each file - located in the column with "Age (ka)"
###############################################

def find_age_column(columns):
    for col in columns:
        if "Age (ka" in col:
            return col
    return None

###############################################
# Here I define a function to extract the superkingdom/domain from the taxa_path string. Depending 
# on the LCA file (either prefilter or metadmg), the superkingdom is either called "superkingdom" 
# or "domain". Also, the LCA files containing prokaryote deads, have a prefix in front of the domain
# abd tgus is removed. So bascically, I first look for superkingdom, but if that is not present, I
# use domain instead. It is not present for the eukaryote LCA files. 
###############################################

def get_superkingdom(taxa_path):
    raw = None
    for level in taxa_path.split(";"):
        parts = level.split(":")
        if len(parts) < 3:
            continue

        name = parts[1].strip('"')
        rank = parts[2].strip('"').lower()

        if rank == "superkingdom":
            if name.startswith("d__"):
                name = name[3:]
            raw = name
        elif rank == "domain" and raw is None:
            raw = name

    return raw

###############################################
# So here I defne a function to deal with entries in the prefilter LCA-files that have no 
# domain are returned as unclassified and any eukaryotes are returned as misclassified.
###############################################

def classify_prefilter(raw):
    if raw is None:
        return "unclassified"
    if raw == "Eukaryota":
        return "misclassified"
    if raw in ("Bacteria", "Archaea", "Viruses"):
        return raw
    return "unclassified"

###############################################
# And here for the eukaryote LCA-files (metadmg), I classify all eukaryotes as what they are. If there is no data
# for superkingdom in the files, I return unclassified and if there is a superkingdom, but it is not Eukaryota, I 
# return misclassified. 
###############################################

def classify_metadmg(raw):
    if raw is None:
        return "unclassified"
    if raw == "Eukaryota":
        return "Eukaryota"
    return "misclassified"

###############################################
# Now, I define a function to extract the file names from the LCA files - so that I can match them to the 
# fastq-reads
###############################################

def extract_sample_id(filename):
    m = re.search(r"(CGG-[0-9]-[0-9]{6})", filename)
    return m.group(1) if m else None

###############################################
# Here I defiene a function to load the lca files (with the renamed data_type as group) and extract age, 
# superkingdom/domain as well as the sample ID. I do so by using the functions defined above. 
###############################################

def load_lca_dir(input_dir, group, prefiler_or_metadmg):

    files = sorted(
        os.path.join(input_dir, f)
        for f in os.listdir(input_dir)
        if f.endswith(".lca") or f.endswith(".lca.gz")
    )

    if not files:
        return pd.DataFrame(columns=["group", "age", "superkingdom", "sample_id"])

    records = []

    for path in files:
        sample_id = extract_sample_id(os.path.basename(path))

        with open_lca_file(path) as fh:
            df = pd.read_csv(fh, sep="\t", low_memory=False)

        age_col = find_age_column(df.columns)

        df["age"] = pd.to_numeric(df[age_col], errors="coerce")
        df = df.dropna(subset=["age"])

        superkingdom_both_types = df["taxa_path"].astype(str).apply(get_superkingdom)

        if prefiler_or_metadmg == "prefilter":
            superkingdom = superkingdom_both_types.apply(classify_prefilter)
        else:
            superkingdom = superkingdom_both_types.apply(classify_metadmg)

        records.append(pd.DataFrame({"group": group, "age": df["age"], "superkingdom": superkingdom, "sample_id": sample_id
        }))

    return pd.concat(records, ignore_index=True)

###############################################
# Now I define a function to turn the data into pivot tables for plotting. The data is divided
# into SS and DS libraries (data type/group) and grouped by age and superkingdom and any values  
# with no values gets a "0"
###############################################

def make_normalized_pivots(df_all, pct_col):
    df = df_all.copy()
    df["value"] = df[pct_col]

    ss = df[df["group"] == "SS"].groupby(["age", "superkingdom"])["value"].sum().unstack(fill_value=0)
    ds = df[df["group"] == "DS"].groupby(["age", "superkingdom"])["value"].sum().unstack(fill_value=0)

    ages = sorted(ss.index.unique())
    ss = ss[CATEGORIES].reindex(ages, fill_value=0.0)
    ds = ds[CATEGORIES].reindex(ages, fill_value=0.0)
    return ss, ds, ages

###############################################
# Now I define a function for plotting, first by creating the figure with room for two subplots (SS and DS)
# Nect, I simply define the labels, sizes and structure of how I wish the plot to look, making sure that 
# font size is always size 16, so that it is easy to read. Also, I add a legend.+
###############################################

def plot_two_panel(ss_pivot, ds_pivot, ages, title_suffix, filename, x_label):
    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(14, 8))

    ss_pivot.plot(kind="barh", stacked=True,
                  color=[COLORS[k] for k in CATEGORIES],
                  ax=ax1, width=0.8, legend=False)

    ax1.set_title(f"SS {title_suffix}", fontsize=16)
    ax1.set_xlabel(x_label, fontsize=16)
    ax1.set_ylabel("Age (ka)", fontsize=16)
    ax1.set_yticks(range(len(ages)))
    ax1.set_yticklabels([f"{a:.2f}" for a in ages], fontsize=16)
    ax1.tick_params(axis="x", labelsize=16)
    ax1.tick_params(axis="y", labelsize=16)
    ax1.invert_yaxis()

    ds_pivot.plot(kind="barh", stacked=True,
                  color=[COLORS[k] for k in CATEGORIES],
                  ax=ax2, width=0.8, legend=True)

    ax2.set_title(f"DS {title_suffix}", fontsize=16)
    ax2.set_xlabel(x_label, fontsize=16)
    ax2.set_yticks(range(len(ages)))
    ax2.set_yticklabels([f"{a:.2f}" for a in ages], fontsize=16)
    ax2.tick_params(axis="x", labelsize=16)
    ax2.tick_params(axis="y", labelsize=16)
    ax2.invert_yaxis()

    handles, labels = ax2.get_legend_handles_labels()
    fig.legend(handles, labels, title="Category",
               fontsize=16, title_fontsize=16,
               bbox_to_anchor=(1.02, 0.5), loc="center left")

    plt.tight_layout()
    plt.subplots_adjust(right=0.85)
    plt.savefig(os.path.join(output_dir, filename), dpi=300)
    plt.close(fig)

###############################################
# Finally I get to wrap all the functions above into the main where i load all the datasets, merge with fastq reads
# and normalize the number of lca reads with the number of fastq reads and generate the final figure. 
###############################################

if __name__ == "__main__":

    ss_meta = load_lca_dir(metadmg_ss_dir, "SS", "metadmg")
    ds_meta = load_lca_dir(metadmg_ds_dir, "DS", "metadmg")
    ss_pref = load_lca_dir(prefilter_ss_dir, "SS", "prefilter")
    ds_pref = load_lca_dir(prefilter_ds_dir, "DS", "prefilter")

    dfs = [df for df in [ss_meta, ds_meta, ss_pref, ds_pref] if not df.empty]
    if not dfs:
        raise SystemExit("No valid LCA data found.")

    df_all = pd.concat(dfs, ignore_index=True)

    df_all = df_all.merge(
        fastq_df[["sample_id", "group", "n_reads_after_QC"]],
        on=["sample_id", "group"],
        how="left"
    )

    # Each LCA row represents one read → normalize by total FASTQ reads
    df_all["pct_reads"] = 100 / df_all["n_reads_after_QC"]

    ss_pct_reads, ds_pct_reads, ages = make_normalized_pivots(df_all, "pct_reads")

    plot_two_panel(
        ss_pct_reads,
        ds_pct_reads,
        ages,
        "% FASTQ reads (unique)",
        "SS_DS_pct_fastq_unique_reads.png",
        "% of FASTQ reads"
    )
