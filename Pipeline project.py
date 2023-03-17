import os
from Bio import SeqIO
import subprocess
from pathlib import Path
from Bio import Entrez
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

### to generate a pipeline that isn't hard coded ###

# home: would contain something like "/Users/james"
home = str(Path.home())
# path: the path to where the output files will be stored
path = home + "/Desktop/PipelineProject/Results/"

# using the sample data (Donor 1) for testing
    # SRA accession numbers for the two samples
sample = "Donor 1"
accessions = ["SRX2896360", "SRX2896363"]

    # loop over the accessions to download the paired-end fastq files
for item in accessions:
    # define the output directory for the fastq files
    output_dir = f"{item}_fastq"
    os.makedirs(output_dir, exist_ok=True)

    # use fasterq-dump to convert the SRA files to paired-end fastq files
    os.system(f"fasterq-dump --split-3 --outdir {output_dir} {item}")


### Track 2 Pipeline ####

#1. Using Bowtie2, create an index for HCMV (NCBI accession NC_006273.2).
# Next, save only the reads that map to the HCMV index for use in assembly.
# Write to your log file the number of reads in each transcriptome before and after the Bowtie2 mapping

# to create the index (assumes a user already has bowtie2 installed):
def question1(file):
    accession_number = "NC_006273.2"

    # create a log file to store the output
    logfile = path + "bowtie2.log.txt"

    # open the log file in append mode
    with open(logfile, "a") as f:
        # use f"bowtie2-build command to run bowtie 2 (code adapted from Github)
        # redirect the stdout and stderr to the log file
        subprocess.run(f"bowtie2-build {path} {accession_number}", shell=True, stdout=f, stderr=f)

    # open the Bowtie2 log file
    log_file = open("logfile", 'r')

    total_reads = 0 # total number of reads
    aligned_pairs = 0 # number of aligned reads
    filtered_pairs = 0 # number of filtered reads

    # parse through log file
    for line in log_file:
        # the line that reports the number of aligned read pairs will say "aligned concordantly exactly 1 time"
        # to count aligned reads, count when in the file that phrase appears
        if 'aligned concordantly exactly 1 time' in line:
            aligned_pairs = int(line.split()[0])
        # same as above for phrase "pairs aligned concordantly 0 times"
        elif 'pairs aligned concordantly 0 times' in line:
            filtered_pairs = int(line.split()[0])

    # calculate the total number of read pairs
    total_reads = (aligned_pairs + filtered_pairs) * 2

    # to write a log file that defines the number of reads before and after bowtie2
    # used the w+ command to create new file and read or write to that file (if not already existing)
    # this stores the log file in the path extablished above

    logFile = open(path + "PipeLineProject.log", "w+")
    logFile.write(sample + " has " + total_reads + " before BowTie2 filtering and " + filtered_pairs + " read pairs after.")


# 2. Using the Bowtie2 output reads, assemble all four transcriptomes together to produce 1 assembly via SPAdes.
# Write the SPAdes command you used to the log file.
    getAssembly = "/Users/zhasan/Downloads/SPAdes-3.15.4-Darwin/bin/spades.py	-k	55,77,99,127	-t	2	--only-assembler	-s/Users/zhasan/Desktop/PipeLineProject/SRAdata/SRR8185310_pass.fastq.gz	-o	/Users/zhasan/Desktop/PipeLineProject/SRAdata/SRA_assembly"
    subprocess.call(getAssembly, shell=True)

    # write SPAdes command to a log file
    logFile = open(path + "PipeLineProject.log", "w+")
    logFile.write("The spades command is " + getAssembly)

# 3a. Write Python code to calculate the number of contigs with a length > 1000
# and write the # to the log file as follows (replace # with the correct integer):

# initilize an empty list to store contigs > 1000
    contigList1000 = []

    # retrive fasta files from the folder designated in the path
    with open("path + contigs.fasta/") as handle:
        # using seqIO, parse through the fasta file
        # from seq/io documentation, use record.seq to find the contig number
        # add the record of the contigs > 100 to the list
        for record in SeqIO.parse(handle, "fasta"):
            if len(record.seq) > 1000:
                contigList1000.append(record)

    # write number of records that have contigs > 1000 in log file
    logFile = open(path + "PipeLineProject.log", "w+")
    logFile.write("There are %i contigs > 1000 in the assembly" % getAssembly)

    # store the contents of contigList1000 to a fasta file to be used below using seqIO write (per documentation)
    SeqIO.write(contigList1000, path + "1000contigs.fasta", "fasta")

# 3b. Write Python code to calculate the length of the assembly (the total number of bp in all of the contigs > 1000 bp in length)
# and write this # to the log file as follows (replace # with the correct integer):

    num = 0
    file1000contigs = SeqIO.parse("/Users/zhasan/Desktop/PipeLineProject/Results/1000contigs.fasta", "fasta")
    for x in file1000contigs:
        num = num + len(x)
    countlength = ("There are " + str(num) + " bp in the assembly.")

    # write total number of bp in all the contigs > 1000 bp in length in log file
    logFile = open(path + "PipeLineProject.log", "w+")
    logFile.write(countlength)

# 4. To determine if the assembly aligns with other virus strains
    # a. Python code to retrieve the longest contig from your SPAdes assembly
# Open the assembly file, and parse through the assembly file using SeqIO
assembly_file = "1000contigs.fasta"
assembly = SeqIO.parse(assembly_file, "fasta")

# Initialize variables
max_length = 0
longest_contig = ""

# Iterate over the contigs in the assembly
for contig in assembly:
    # Check if the length of the contig is longer than the current max length
    if len(contig.seq) > max_length:
        # If so, update the max length and the longest contig
        max_length = len(contig.seq)
        longest_contig = contig

    # b. Python code to use the longest contig as blast+ input to query the nr nucleotide database limited to members of the Betaherpesvirinae subfamily:

# creating a local database of just sequences from the Betaherpesvirinae subfamily
Entrez.email = "zhasan@luc.edu" # Set email address for Entrez

search_terms = "Betaherpesvirinae[Subtree] AND srcdb_genbank[PROP]" # Define search terms for Betaherpesvirinae subfamily in GenBank, code taken from documentation

# Use Entrez to search and fetch IDs of matching sequences
# code taken from documentation
handle = Entrez.esearch(db="nucleotide", term=search_terms, retmax=1000)
record = Entrez.read(handle)
handle.close()

# Using the IDs obtained from the search, fetch the specific family's sequences
handle = Entrez.efetch(db="nucleotide", id=record["IdList"], rettype="fasta", retmode="text")
sequences = list(SeqIO.parse(handle, "fasta"))
handle.close()

# Write the sequences to a local FASTA file to be used to parse through
with open("betaherpesvirinae.fasta", "w") as output_file:
    SeqIO.write(sequences, output_file, "fasta")

# using a blastn to query the database
# making sure to keep the best alignment for single query-subject pair of sequences

# run BLAST+ on the query sequence against the nr nucleotide database
# using NcbiblastnCommandline from the package Blast Applications
# code adapted from documentation
blastn_cline = NcbiblastnCommandline(query="1000contigs.fasta", db="nr", evalue=0.001, outfmt=5, out="blast_results.xml")
stdout, stderr = blastn_cline()

# parse the BLAST results and filter for matches to Betaherpesvirinae subfamily
# using NCBIXML from the package BLAST
# code adapted from documentation
blast_results = NCBIXML.parse(open("blast_results.xml"))
hits = []
for result in blast_results:
    for alignment in result.alignments:
        if "Betaherpesvirinae" in alignment.hit_def:
            hits.append(alignment.hit_def)

# add top 10 hits to log file
header = "sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle\n" # header for log file
with open(path + ".log", "w+") as f:
    f.write(header)
    # top 10 items and their counts to the log file
    for i, hit in enumerate(sorted(hits, key=lambda x: -hits.count(x))[:10]):
        f.write(f"{hit}\t{hits.count(hit)}\n")






