import os
import re
import pandas as pd
from glob import glob

# ---------------------------------------------------------
# Directories
# ---------------------------------------------------------

metadata_file = "/maps/projects/prohaska/people/vpt968/Data_analysis_revised_data/Data/Age_files/GBR_Hole58A_samples_list_APB.xlsx"
lca_folder = "/maps/projects/prohaska/people/vpt968/Data_analysis_revised_data/Original_files/Data_SS_lower_threshold/results/metadmg/lca"
output_folder = "/maps/projects/prohaska/people/vpt968/Data_analysis_revised_data/Results_of_three_data_types/Updated_mapping/lca_files_with_age/prefilter_metadmg_SS"

os.makedirs(output_folder, exist_ok=True)

# ---------------------------------------------------------
# Here I load and read the excel file, skipping the first row. 
# ---------------------------------------------------------

metadata_df = pd.read_excel(metadata_file, skiprows=1) 
metadata_df.columns = metadata_df.columns.str.strip()

print("Columns loaded:", metadata_df.columns.tolist())

# ---------------------------------------------------------
# There were some issues loading the CGG-IDs so I normalize them as follows
# ---------------------------------------------------------
metadata_df["Normalized_CGG_ID"] = metadata_df["CGG ID"].astype(str).str.replace("-", "_", regex=False)

# ---------------------------------------------------------
# Here I find the LCA-files
# ---------------------------------------------------------

lca_files = glob(os.path.join(lca_folder, "*lca.gz"))
print(f"Found {len(lca_files)} LCA files in {lca_folder}")

# ---------------------------------------------------------
# And extract the CGG IDs that are in the LCA-filenames
# ---------------------------------------------------------
lca_ids = set()
for f in lca_files:
    match = re.search(r"(CGG[-_]?\d{1,2}[-_]?\d{6})", os.path.basename(f))
    if match:
        lca_ids.add(match.group(1).replace("-", "_"))

# ---------------------------------------------------------
# Here I create a dictionary based on the headers of the lca-files
# ---------------------------------------------------------

metadata_filtered = metadata_df[metadata_df["Normalized_CGG_ID"].isin(lca_ids)].copy()
metadata_dict = metadata_filtered.set_index("Normalized_CGG_ID").to_dict(orient="index")

# ---------------------------------------------------------
# Next is the processing, where I first find the lca files based on their name
# Then I normalize the CGG IDs from the excel sheet so that they can be matched to
# the corresponding LCA-files And I add the meta data columns (age and depth) to each row in the corresponding LCA-file
# ---------------------------------------------------------

for filepath in lca_files:
    filename = os.path.basename(filepath)
    match = re.search(r"(CGG[-_]?\d{1,2}[-_]?\d{6})", filename)
    if not match:
        continue

    cgg_id = match.group(1).replace("-", "_")

    if cgg_id not in metadata_dict:
        continue

    try:
        df = pd.read_csv(filepath, sep="\t", comment="#")
        df.columns = df.columns.str.strip()

        meta = metadata_dict[cgg_id]
        for col in ["Top depth (cm)", "Bottom depth (cm)", "Age (ka) Brandon et al. 2015"]:
            df[col] = meta[col]

        new_filename = filename.replace("collapsed", "with_depth")
        output_path = os.path.join(output_folder, new_filename)
        df.to_csv(output_path, sep="\t", index=False)

