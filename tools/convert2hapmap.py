import pandas as pd
import sys

def convert_agriplex_to_hapmap(input_file, output_file):
    # Load the Excel file
    xls = pd.ExcelFile(input_file)
    df = pd.read_excel(xls, sheet_name=xls.sheet_names[0], header=None)

    # IUPAC encoding dictionary
    iupac_dict = {
        "A / T": "W", "T / A": "W",
        "C / ": "S", "G / C": "S",
        "A / C": "M", "C / A": "M",
        "G / T": "K", "T / G": "K",
        "A / G": "R", "G / A": "R",
        "C / T": "Y", "T / C": "Y",
        "A": "A", "C": "C", "G" : "G", "T": "T",
        "--": "-", "FAIL": "N", "NN": "N"  # Handle missing data
    }
    
    # Extract data based on provided structure
    snp_markers = df.iloc[3, 3:].values  # D4 onwards
    alt_snp_names = df.iloc[4, 3:].values  # D5 onwards
    ref_alleles = df.iloc[5, 3:].values  # D6 onwards
    alt_alleles = df.iloc[6, 3:].values  # D7 onwards
    sample_names = df.iloc[7:, 2].values  # C7 onwards
    genotype_data = df.iloc[7:, 3:].values  # D7 onwards
    
    # Remove invalid header row if present
    if snp_markers[0] == "Customer Marker ID":
        snp_markers = snp_markers[1:]
        alt_snp_names = alt_snp_names[1:]
        ref_alleles = ref_alleles[1:]
        alt_alleles = alt_alleles[1:]
        genotype_data = genotype_data[:, 1:]
    
    # Construct HapMap format headers
    hapmap_headers = [
        "rs#", "alleles", "chrom", "pos", "strand", "assembly#", "center",
        "protLSID", "assayLSID", "panelLSID", "QCcode"
    ] + list(sample_names)

    # Prepare the HapMap data
    hapmap_rows = []
    dropped_markers = []
    for i, marker in enumerate(snp_markers):
        try:
            alt_snp_names_tokenized = alt_snp_names[i].split("_")
            # Formats:
            #   starts with chr: chrom#_pos#
            #   starts with IRGSP1: IRRI_SNP#_IRGSP1_C#_pos#
            #   starts with MSU7: IRRI_SNP#_MSU7_Chrom#_pos#_ref-alt
            #   only QTL names: IRRI_SNP#_QTLname_meta

            if marker.startswith('IRGSP1'):
                chrom = alt_snp_names_tokenized[3]
                if chrom.startswith('C'):
                    chrom = chrom.replace('C','')
                pos = alt_snp_names_tokenized[4]
            elif marker.startswith('MSU7'):
                chrom = alt_snp_names_tokenized[3]
                pos = alt_snp_names_tokenized[4]
            elif marker.startswith('chr'):
                chrom, pos = marker.split("_")[0][3:], marker.split("_")[1]  # Extract chromosome and position
            else:
                chrom, pos = "NA", "NA"  # if chromosome and position not specified
        except IndexError:
            chrom, pos = "NA", "NA"  # if chromosome and position not specified
        
        alleles = f"{ref_alleles[i]}/{alt_alleles[i]}"

        # Drop multi-allelic and indel sites (only keep single A/T/C/G alleles)
        if len(ref_alleles[i]) > 1 or len(alt_alleles[i]) > 1:
            dropped_markers.append(marker)
            continue  # Skip this marker

        # Convert genotype data to IUPAC format
        iupac_genotypes = []
        for allele in genotype_data[:, i]:
            iupac_allele = iupac_dict.get(allele, "N")  # Default to 'N' if not found
            iupac_genotypes.append(iupac_allele)

        row = [
            alt_snp_names[i], alleles, chrom.lstrip("0"), pos, "+", "NA", "NA", "NA", "NA", "NA", "NA"
        #] + list(genotype_data[:, i])
        ] + iupac_genotypes
        hapmap_rows.append(row)
    
    # Create DataFrames
    hapmap_df = pd.DataFrame(hapmap_rows, columns=hapmap_headers)
    filtered_df = hapmap_df[~(hapmap_df['chrom'] == 'NA') & ~(hapmap_df['pos'] == 'NA')]
    
    # Save as TSV file
    filtered_file = output_file + '.filtered'
    hapmap_df.to_csv(output_file, sep="\t", index=False)      # all bi-allelic
    filtered_df.to_csv(filtered_file, sep="\t", index=False)    # all non-NA chromosome
    print(f"HapMap file saved to: {output_file}")

    # Output Dropped Markers
    #dropped_markers_df = pd.DataFrame(dropped_markers)
    #dropped_markers_df.to_csv("dropped_markers.tsv", sep="\t", index=False)
    #print(f"Dropped markers file saved to: dropped_markers.tsv")


def convert_dart_to_hapmap(input_file, output_file):
    # Load spreadsheet file/csv file
    data = pd.read_csv(input_file, sep = ",", header = 1)
    data_transposed = data.transpose()

    data_transposed.to_csv(output_file, sep="\t", index=False)
    print(f"HapMap file saved to: {output_file}")
    #print("DaRT Input")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python3 convert2hapmap.py <input_type: 1- 1kRiCA, 2- DaRT> <input_excel_file> <output_tsv_file>")
    else:
        if sys.argv[1] == "1":
            #print("1KRiCA Input")
            convert_agriplex_to_hapmap(sys.argv[2], sys.argv[3])
        elif sys.argv[1] == "2":
            convert_dart_to_hapmap(sys.argv[2], sys.argv[3])
        else:
            print("Usage: python3 convert2hapmap.py <input_type: 1- 1kRiCA, 2- DaRT> <input_excel_file> <output_tsv_file>")
