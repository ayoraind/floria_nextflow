#!/usr/bin/env python3

import gzip
import sys

def main(vcf_file, new_vcf):
    refs = set()

    # Read the VCF file to extract contig names
    with gzip.open(vcf_file, 'rt') as f:
        refs = {line.split()[0] for line in f if line != "" and line[0] != '#'}

    # Sort the reference names
    refs = sorted(refs)

    # Write the new VCF file with the contig information
    with gzip.open(vcf_file, 'rt') as f:
        with open(new_vcf, 'w') as fo:
            for idx, line in enumerate(f):
                if idx == 2:  # Write contig information after the header lines
                    for ref in refs:
                        fo.write(f"##contig=<ID={ref}>\n")
                
                fo.write(line)

if __name__ == "__main__":
    # Check if the correct number of arguments is passed
    if len(sys.argv) != 3:
        print("Usage: python floria_vcf_header.py <input_vcf.gz> <output_vcf>")
        sys.exit(1)

    vcf_file = sys.argv[1]  # First argument is the input VCF file
    new_vcf = sys.argv[2]   # Second argument is the output VCF file

    main(vcf_file, new_vcf)
