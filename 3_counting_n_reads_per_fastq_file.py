import gzip
import glob
import os
import pandas as pd

###############################################
# Directives
###############################################

FASTQ_DIR = "/maps/projects/prohaska/people/vpt968/Data_analysis_revised_data/Original_files/Fastq"

###############################################
# Each fastq file has four lines per. Therefore, the function for counting reads must divide the number of rows with 4
# Therefore, the function for counting reads must divide the number of rows with 4.
# the += adds the count of one line to the following, starting at 0, adding 1 to the value of line_count.
###############################################

def count_reads_fastq(file_path):
    line_count = 0
    with gzip.open(file_path, "rt", errors="ignore") as f:
        for _ in f:
            line_count += 1
    return line_count // 4  


# --------------------------------------------------
# In order to match the fastq files to the LCA files i need to extract something they have in common.
# Here I define the function to extract the CGG-ID and the number of rows/4. 
# --------------------------------------------------

def extract_sample_id(filename):
    parts = filename.split("-")
    for i in range(len(parts)):
        if parts[i].startswith("CGG"):
            return "-".join(parts[i:i+3])
    return None


# --------------------------------------------------
# Here is the function find the fastq files in the folder (ending with ".fastq.gz") 
# --------------------------------------------------

def main():
    fastq_files = []
    for ext in *.fastq.gz:
        fastq_files.extend(glob.glob(os.path.join(FASTQ_DIR, "**", ext), recursive=True))

    if not fastq_files:
        print("No FASTQ files found.")
        return

    records = []

# --------------------------------------------------
# Extracting the sample ID and counting using the counting function defined above. 
# --------------------------------------------------
    for fq in fastq_files:
        fname = os.path.basename(fq)
        sample_id = extract_sample_id(fname)

        if sample_id is None:
            print(f"WARNING: Could not extract sample ID from {fname}")
            continue

        print(f"Counting reads in {fname} ...")
        reads = count_reads_fastq(fq)

        records.append({
            "file": fname,
            "full_path": fq,
            "sample_id": sample_id,
            "reads": reads
        })

    df = pd.DataFrame(records)

# --------------------------------------------------
# Here I sum the r1 and r2 files each sample
# --------------------------------------------------
    summary = (
        df.groupby("sample_id")["reads"]
        .sum()
        .reset_index()
        .rename(columns={"reads": "total_reads"})
    )

# --------------------------------------------------
# And write it to an excel sheet with pd-ExcelWriter + defining the name of the excel file 
# --------------------------------------------------
    out_path = "fastq_read_counts.xlsx"
    with pd.ExcelWriter(out_path) as writer:
        df.to_excel(writer, sheet_name="files", index=False)
        summary.to_excel(writer, sheet_name="summary_per_sample", index=False)

    print(f"\n✔ DONE! Output written to: {out_path}")


# --------------------------------------------------
if __name__ == "__main__":
    main()
