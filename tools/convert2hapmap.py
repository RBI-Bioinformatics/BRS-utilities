import pandas as pd
import sys

def check_marker_snptable(marker,none_count):
    marker_info = snptable.get(marker.lower())
    if marker_info is not None:
        chrom,pos = marker_info.split(",")[2], marker_info.split(",")[3]
        alleles = marker_info.split(",")[4] + "/" + marker_info.split(",")[5]
    else:
        alleles, chrom, pos = "N/N", 999, none_count  # if chromosome and position not specified
    
    return(alleles, chrom, pos)

def convert_agriplex_to_hapmap(input_file, output_file):
    # Load the Excel file
    xls = pd.ExcelFile(input_file)
    df = pd.read_excel(xls, sheet_name=xls.sheet_names[0], header=None)

    # IUPAC encoding dictionary
    global iupac_dict
    iupac_dict = {
        "A/T": "W", "T/A": "W",
        "C/G": "S", "G/C": "S",
        "A/C": "M", "C/A": "M",
        "G/T": "K", "T/G": "K",
        "A/G": "R", "G/A": "R",
        "C/T": "Y", "T/C": "Y",
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

    # add 999 if chromosome unknown(
    # add counter if pos unknown(Pos+ counter, example: 1, 2, 3,...,n)
    count = 1
    for i, marker in enumerate(alt_snp_names):
        try:
            alleles, chrom, pos = check_marker_snptable(marker, count)
            count += 1
        except IndexError:
            # if chromosome and position not specified
            alleles, chrom, pos = "N/N", 999, int(count)  

        # Convert genotype data to IUPAC format
        iupac_genotypes = []
        for allele in genotype_data[:, i]:
            iupac_allele = iupac_dict.get(allele.replace(" ", ""), allele)  # Default to allele if not found
            iupac_genotypes.append(iupac_allele)

        row = [
            alt_snp_names[i], alleles, chrom, pos, "+", "NA", "NA", "NA", "NA", "NA", "NA"
        ] + iupac_genotypes
        hapmap_rows.append(row)
        
    # output dataframe
    hapmap_df = pd.DataFrame(hapmap_rows, columns=hapmap_headers)
    hapmap_df.to_csv(output_file, sep="\t", index=False)
    print(f"HapMap file saved to: {output_file}")


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
    # Columns as is
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
    count = 1
    for i, marker in enumerate(marker_names):
        try:
            alleles, chrom, pos = check_marker_snptable(marker, count)
            count +=1
        except IndexError:
            alleles, chrom, pos = "N/N", "999", int(count)  # if chromosome and position not specified

        # Get genotype calls for this marker
        genotypes = snp_data[marker].replace({".:.": "N/N", "-:-": "N/N"}).tolist()
    
        # Format the row
        hapmap_rows.append([marker, alleles, chrom, pos, "+", "NA", "NA", "NA", "NA", "NA", "NA"] + genotypes)

    # Create a DataFrame for HapMap output
    # Save to file (tab-delimited)
    hapmap_df = pd.DataFrame(hapmap_rows, columns=hapmap_headers)
    hapmap_df.to_csv(output_file, sep="\t", index=False)
    print(f"HapMap file saved as {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python3 convert2hapmap.py <input_type: 1- 1kRiCA, 2- DaRT> <snptable_file> <input_excel_file> <output_tsv_file>")
    else:
        # create snp table
        create_snptable_dictionary(sys.argv[2])
        if sys.argv[1] == "1":
            convert_agriplex_to_hapmap(sys.argv[3], sys.argv[4])
        elif sys.argv[1] == "2":
            convert_dart_to_hapmap( sys.argv[3], sys.argv[4])
        else:
            print("Usage: python3 convert2hapmap.py <input_type: 1- 1kRiCA, 2- DaRT> <snptable_file> <input_excel_file> <output_tsv_file>")