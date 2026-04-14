import os
import re
import gzip
import pandas as pd
from pathlib import Path

###############################################
# Directories
###############################################

LCA_DIRS = {"SS": "//maps/projects/prohaska/people/vpt968/Data_analysis_revised_data/Results_of_three_data_types/Updated_mapping/lca_files_with_age/_metadmg_SS",
    "DS": "/maps/projects/prohaska/people/vpt968/Data_analysis_revised_data/Results_of_three_data_types/Updated_mapping/lca_files_with_age/metadmg_DS"
}
OUTPUT_EXCEL = "all_phyla_genus_summary.xlsx"

###############################################
# Here I define the different ranks found in the taxa_path 
###############################################

TAXON_RANKS = [
    "species", "genus", "family", "order", "class", "subclass",
    "superorder", "superfamily", "subfamily", "infraorder",
    "cohort", "infraclass", "phylum", "kingdom", "domain"
]

###############################################
# Here I define a function to clean the taxa_path, because as it is now, there are a lot of ""s
# and they make it really difficult to cleanly extract the ranks. Also trailing whitespaces and
# anything with no value is handled. 
###############################################


def clean_taxpath(taxa_path):
    if taxa_path is None:
        return ""
    clean_string = str(taxa_path).strip()
    if clean_string.startswith('"') and clean_string.endswith('"'):
        clean_string = clean_string[1:-1]
    while '""' in clean_string:
        clean_string = clean_string.replace('""', '"')
    return clean_string

###############################################
# Next is a function to parse the taxa_path - in other words to split the taxa_path string (s)
# into segments (seg) containing the taxonomic level and the name of it. There were some issues 
# with the ""'s again, so I remove them once more. The two "mached" parts were quite the hussle. But
# I noticed that some of the taxa were not passed to the resulting excel. And that was because some of 
# the taxa taxonomic levels in the taxa path did not have ""'s. So the second part is to handle them.
###############################################

def parse_taxapath(taxa_path):
    s = clean_taxpath(taxa_path)
    taxa = {}
    for seg in s.split(";"):
        seg = seg.strip()
        while '""' in seg:
            seg = seg.replace('""', '"')

        matched = re.match(r'(\d+):"([^"]+)"\s*:\s*"([^"]+)"', seg)
        if matched:
            taxid, name, rank = matched.groups()
            taxa[rank] = name
            continue

        matched2 = re.match(r'(\d+):"([^"]+)"\s*:\s*([A-Za-z0-9 _-]+)', seg)
        if matched2:
            taxid, name, rank = matched2.groups()
            taxa[rank] = name
            continue

    return taxa

###############################################
# Now I define a function to extract the ranks from the taxa_path string that I parsed above
# Into out, I now extract the higest possible rank and the name for that. 
###############################################

def extract_taxonomy_all_ranks(taxa_path):
    raw = parse_taxapath(taxa_path)
    out = {r: raw.get(r) for r in TAXON_RANKS}

    for r in [
        "species","genus","family","order","class","subclass",
        "superorder","superfamily","subfamily","infraorder",
        "cohort","infraclass"
    ]:
        if out.get(r) is not None:
            out["highest_taxonomic_resolution"] = out[r]
            out["highest_taxonomic_rank"] = r
            break
    else:
        out["highest_taxonomic_resolution"] = None
        out["highest_taxonomic_rank"] = None

    return out

###############################################
# Here is a function to extract the file names from the LCA files. 
###############################################

SAMPLE_ID_REGEX = re.compile(r"(CGG-1-\d+)")

def extract_sample_id(filename):
    m = SAMPLE_ID_REGEX.search(filename)
    return m.group(1) if m else None

###############################################
# Now I load all the LCA-files. This is done for all files ending with lca.gz. All lines in all files
# are read and taxa_path and the age of the samples are extracted. all lines without "phyla" in it are discarded.
# All Suidae and Homo sapiens reads are removed, because i am certain that they are contaminations. They took up a LOT
# of space in the LCA-files when I inspected them manually. Lastly, the the information is stored in a list and return
# it as a pandas dataframe.  
###############################################

def load_all_lca_records():
    records = []
    for datatype, folder in LCA_DIRS.items():
        folder = Path(folder)
        if not folder.exists():
            continue
        for fn in os.listdir(folder):
            if not fn.endswith(".lca.gz"):
                continue
            file_path = folder / fn
            sample_id = extract_sample_id(fn)
            with gzip.open(file_path, "rt") as f:
                header = True
                for line in f:
                    if header:
                        header = False
                        continue

                    parts = line.strip().split("\t")
                    if len(parts) < 10:
                        continue

                    taxa_path = parts[6]
                    age = float(parts[9])

                    tax = extract_taxonomy_all_ranks(taxa_path)
                    phylum = tax.get("phylum")

                    if phylum is None:
                        continue

                    # Filter out Suidae / Homo sapiens contamination
                    cleaned = clean_taxpath(taxa_path)
                    if (phylum == "Chordata" and
                        ('9821:"Suidae"' in cleaned or '9606:"Homo sapiens"' in cleaned)):
                        continue

                    row = {
                        "sample_id": sample_id,
                        "data_type": datatype,
                        "phylum": phylum,
                        "age": age,
                        "taxa_path": taxa_path,
                    }
                    # add all ranks + highest
                    row.update(tax)
                    records.append(row)

    return pd.DataFrame(records)

###############################################
# Now I deine a function to build the excel sheet. I define that the each phylum gets
# its own sheet. I get the phyla and build the excel sheets after their names. 
# Next, I create a subdataframe (sub) to contain the reads belonging to the respecive phyla. 
# and I create an empty list to contain the summary rows. I did this quite some time ago, 
# and there are definetly some thing I would have done in another way if I knew what I know now. 
# BUT, I then basically group the highest taxonomic rank (both the unique name and the taxonomic level) 
# with the dataframe for the specific phyla. Further the "g" sums all the mappings with the same highest 
# taxonomic level. Then i create a string to show in what samples (what ages) the reads of a specific 
# taxonimic lecel is present at. Also I make the output specify if the reads are only present in my 
# SS libraries, bu DS libraries or both. I then append all the taxonomic levels along with the stuff above for the specific 
# unique entries  I should just have done that to begin with, because I later figured out that species level
# was a no-go in aeDNA.  
###############################################

def build_phylum_genus_excel(df, output_excel):

    writer = pd.ExcelWriter(output_excel, engine="openpyxl")

    phyla = sorted(df["phylum"].unique())

    for phylum in phyla:
        sub = df[df["phylum"] == phylum]

        summary_rows = []

   
        group_cols = ["highest_taxonomic_resolution", "highest_taxonomic_rank"]

        for (res, rank), g in sub.groupby(group_cols):
            if res is None:
                continue

            ages = sorted(g["age"].unique())
            types = sorted(g["data_type"].unique())


            def uniq(col):
                vals = sorted(x for x in g[col].dropna().unique())
                return ", ".join(vals) if vals else None

            summary_rows.append({
                "phylum": phylum,
                "highest_taxonomic_resolution": res,
                "highest_taxonomic_rank": rank,
                "n_reads": len(g),
                "ages_present": ", ".join(f"{a:.2f}" for a in ages),
                "data_types_present": "+".join(types),
                "species": uniq("species"),
                "genus": uniq("genus"),
                "family": uniq("family"),
                "order": uniq("order"),
                "class": uniq("class"),
                "subclass": uniq("subclass"),
                "superorder": uniq("superorder"),
                "superfamily": uniq("superfamily"),
                "subfamily": uniq("subfamily"),
                "infraorder": uniq("infraorder"),
                "cohort": uniq("cohort"),
                "infraclass": uniq("infraclass"),
            })

        summary_df = pd.DataFrame(summary_rows)
        if summary_df.empty:
            continue

        summary_df = summary_df.sort_values(
            by=["highest_taxonomic_rank","highest_taxonomic_resolution"]
        )

        sheet_name = phylum[:31]

        summary_df.to_excel(writer, index=False, sheet_name=sheet_name)

    writer.save()


###############################################
# Here I define the main to run the functions above and create the excel sheet. 
###############################################

def main():
    df = load_all_lca_records()

    print("\nPreview of loaded records:")
    print(df.head())

    build_phylum_genus_excel(df, OUTPUT_EXCEL)


if __name__ == "__main__":
    main()
