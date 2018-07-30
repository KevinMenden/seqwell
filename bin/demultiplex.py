####################################################
## Demultiplex FastQ file using barcodes ###########
####################################################

# Import libraries
import argparse
import os

# Parsoe cl arguments
parser = argparse.ArgumentParser()
parser.add_argument("-bc", help="The barcode file")
parser.add_argument("-f", help="The fastq file")
parser.add_argument("-p", help="Number of cores to use. Default 1", default=1)

args = parser.parse_args()
fq_file = args.f
barcodes = open(args.bc, "r")

# Call process for every barcode
for i, line in enumerate(barcodes):
    cores = args.p
    sp = line.split()
    sample = sp[0]
    barc = sp[1]
    print("Demultiplexing " + sample)

    if i == 0:
        cmd = "cutadapt -g %(sample)s=^%(barc)s -e 0 --no-indels %(fq_file)s -o 'demult_{name}.fastq.gz'" %locals()
    else:
        cmd = "cutadapt -g %(sample)s=^%(barc)s -e 0 --no-indels %(fq_file)s -o 'demult_{name}.fastq.gz'" %locals()

    print(cmd)
    os.system(cmd)


