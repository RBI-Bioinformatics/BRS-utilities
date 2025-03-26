import pandas as pd
import sys

def convert_agriplex_to_hapmap(input_file, output_file):
    # Load the Excel file
    xls = pd.ExcelFile(input_file)
    df = pd.read_excel(xls, sheet_name=xls.sheet_names[0], header=None)

    # IUPAC encoding dictionary
    global iupac_dict 
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
    # add UN as chromosome (unknown+counter, example: Unknown1, Unknown2,...)
    # add Pos1 as pos(Pos+ counter, example: Pos1, Pos2, Pos3)
    count = 1
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
                chrom, pos = 999, int(count)  # if chromosome and position not specified
                count += 1
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
            alt_snp_names[i], alleles, int(chrom), int(pos), "+", "NA", "NA", "NA", "NA", "NA", "NA"
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


def create_snptable_dictionary(snptable_file):
    global snptable
    # Read the file and create dictionary
    # keys are made lowercase
    with open(snptable_file, 'r') as file:
        snptable = {key.lower(): value for key, value in (line.strip().split(',', 1) for line in file)}


def convert_dart_to_hapmap(input_file,  output_file):
    # Read the CSV file
    df = pd.read_csv(input_file)

    # Extract metadata and SNP data
    # Columns will not change
    plate_name = df.iloc[:, 0].astype(str)  
    well = df.iloc[:, 1].astype(str)  
    subject_id = df.iloc[:, 2].astype(str)  
    snp_data = df.iloc[:, 3:] 

    # Concatenate sample names and metadata for unique identifiers
    #sample_ids = plate_name + "_" + well + "_" + subject_id
    sample_ids = subject_id

    # Extract SNP marker names from column headers
    marker_names = snp_data.columns.tolist()

    # Prepare HapMap format columns with the usual Hapmap headers
    # If data has strand and assembly/reference genome info, specify.
    hapmap_headers = ["rs#", "alleles", "chrom", "pos","strand", "assembly#", "center",
        "protLSID", "assayLSID", "panelLSID", "QCcode"] + sample_ids.tolist()
    
    # Initialize HapMap formatted data
    hapmap_rows = []
    for i, marker in enumerate(marker_names):
        try:
            # get the marker from the snp table dictionary
            # marker names are changed to lowercase to compare to snptable[key]
            val = snptable.get(marker.lower())
            if val is not None:
                chrom, pos = val.split(',')[2], val.split(',')[3]
                alleles = val.split(',')[4]+ "/"+ (val.split(',')[5])
            else:
                alleles, chrom, pos = "N/N", "999", "NA"
        except IndexError:
            alleles, chrom, pos = "N/N", "999", "NA"  # if chromosome and position not specified

        # Get genotype calls for this marker
        genotypes = snp_data[marker].replace({".:.": "N/N", "-:-": "N/N"}).tolist()
        # Format the row
        hapmap_rows.append([marker, alleles, chrom, pos, "+", "NA", "NA", "NA", "NA", "NA", "NA"] + genotypes)

    # Create a DataFrame for HapMap output
    hapmap_df = pd.DataFrame(hapmap_rows, columns=hapmap_headers)

    # Save to file (tab-delimited)
    hapmap_df.to_csv(output_file, sep="\t", index=False)
    print(f"HapMap file saved as {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python3 convert2hapmap.py <input_type: 1- 1kRiCA, 2- DaRT> <input_excel_file> <output_tsv_file>")
    else:
        if sys.argv[1] == "1":
            #print("1KRiCA Input")
            convert_agriplex_to_hapmap(sys.argv[2], sys.argv[3])
        elif sys.argv[1] == "2":
            create_snptable_dictionary("DArTSNPs_v4.2.csv")
            convert_dart_to_hapmap(sys.argv[2], sys.argv[3])
        else:
            print("Usage: python3 convert2hapmap.py <input_type: 1- 1kRiCA, 2- DaRT> <input_excel_file> <output_tsv_file>")