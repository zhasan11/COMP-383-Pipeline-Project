# COMP-383-Pipeline-Project

COMP 383 Pipeline Project 

Before you begin, have the following installed: 
1. Python 3.6 
2. SPAdes 3.15.4
3. SRA Tool Kit 3.0.0
4. Biopython 1.81 

This is a Github Repo for py script. It functions as python wrapper to automate execution of software tools for genome assembly. 

The py script will automate and product output file to "PipeLineProject.log" and files in a folder named "PipelineProject_Zarah_Hasan". All results generated will be written to this folder. This folder will be created via an os.system call.

** This Pipeline was developed for Track 2, and answers the following questions: **

1. Determine which strains are most similar to patient samples. 
2. Create 1 assembly using SPAdes.
3. Calculate the number of contigs > 1000 and the total length of the assembly.
4. Determine if the assembly developed aligns with other virus strains. 
