import pandas as pd
import sys

def filter_BiAllelic_from_Hapmap(input_file, output_file):
    # read the tab-delimited file
    # Hapmap delimiter is '\t'
    df = pd.read_csv(input_file, sep="\t")

    # Keep default columns
    df['assembly#']= "NA"
    df["center"]="NA"
    df["protLSID"]="NA"
    df["assayLSID"]="NA"
    df["panelLSID"]="NA"
    df["QCcode"]="NA"

    # ensures that only bi-allelic SNPs are kept.
    # removes those with indels
    df = df[df["alleles"].astype(str).str.match(r"^[^-]/[^-]$")] 
    
    # write to output file
    df.to_csv(output_file, sep="\t", index=False)
    print(f"Filtered HapMap file saved as {output_file}")

def change_to_iupac(input_file, output_file):
    # Read Hapmap file
    df = pd.read_csv(input_file, sep = "\t")

    iupac_dict = {
        "A/T": "W", "T/A": "W",
        "C/G": "S", "G/C": "S",
        "A/C": "M", "C/A": "M",
        "G/T": "K", "T/G": "K",
        "A/G": "R", "G/A": "R",
        "C/T": "Y", "T/C": "Y",
        "A/A": "A", "C/C": "C", 
        "G/G" : "G", "T/T": "T",
        "-/-": "N", "FAIL": "N", "N/N": "N"  # Handle missing data
    }

    #genotype_data = df.iloc[:, [0] + list(range(11, df.shape[1]))]  # readfter Qcode onwards
    # Replace alleles with IUPAC codes
    #genotype_data = genotype_data.replace(iupac_dict)

    # Replace alleles directly in df, keeping the first column and modifying from 12th column onwards
    df.iloc[:, 11:] = df.iloc[:, 11:].replace(iupac_dict)

    # write to output file
    df.to_csv(output_file, sep="\t", index=False)
    print(f"Filtered HapMap file saved as {output_file}")



if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 filterBiAllelicHapmap.py <hapmap_file> <filtered_file>")
    else:
        #filter_BiAllelic_from_Hapmap(sys.argv[1], sys.argv[2])
        change_to_iupac(sys.argv[1], sys.argv[2])