# What is Gsub ?

Gusb :

Gsub is a Graphical User Interface (GUI) tools written in python which allow to annotate and submit a large number of viral sequences.
Gsub uses python package *ORFfinder* to identify and annotate openreading frames (ORFs) in 
the sequences. Then, python package *pyHMMER* is used to detect potential polymerase-encoding ORFs 
(detection threshold can be modify in parameters). 
[Table2asn](https://www.ncbi.nlm.nih.gov/genbank/table2asn/) , a GenBank tools, is then used to collect the different 
information and to generate the corresponding sqn file for each sequence. 
The python packages *Gooey* and [Pyinstaller 5.2](https://pyinstaller.org/) allow respectively for 
a graphical interface and to run Gsub in Windows or Linux without installing a Python interpreter.

# Install

## Download exec

<a href="./Gsub/exec/Gsub_linux" download>Gsub for Linux</a>

<a href="./Gsub/exec/Gsub_windows.exe" download>Gsub for Windows 10</a>


## With python 

Command line

``
pip install gsub
``


# RUNNING GSUB

Image with explanation of all field.

## Input file

* FASTA
* Source 
* Template
* Output

## Output file

* GBF file
* sqn file 

And a error summary file for all sequence with XX error.

## Option

* Length of contigs
* Length of ORF
* Remove overlaps
* Keep only one strand
* PFAM score for identification of polymerase
* Assembly information, becarful Genbank format assembler v. x.x.x.  