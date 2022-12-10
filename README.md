# KMerFineRabbit
KmerFineRabbit is a fast python script for filtering NGS datasets.  Remove or keep reads with a k-mer match to another dataset or a reference genome.

*********
KMerFineRabbit ('KFR')-- Reference-indepdendent manipulation of NGS Datasets

->Fuction: KFR Filters an input NGS dataset by keeping/removing reads with k-mer that match a "Filter" dataset
 - KFR is used to discover reads within a metagenomic dataset that either share homology with a second
 - dataset or lack such homology.  Datasets can be reference genomes or relatively large fastq/fasta read files
 - from high throughput sequencing runs.  Matches are defined by setting a k-mer length (default=32, which ensures
 - that few spurious matches are identified.

->Syntax
 - python KMerFineRabbit<ver>.py SustrateFiles=<file,file...> FilterFiles=<file,file...> Mode=<positive/negative>

->Required Parameters
 - SubstrateFiles=<file[s]> : Starting dataset(s) (reads from these files are filtered to generate your output). 
 - FilterFiles=<file[s]> : Datasets to gather kmers for filtering Substrate
 -   File lists can be comma separated (no spaces) and/or specified with a wildcard ('*')
 -   Pairs of "Substrate" files (e.g.,R1,R2) with paired reads should be separated by a semicolon (not comma)
 -   FilterFile entries that are also in the SubstrateFile list will be ignored (not used for filtering)
 -   Input files can be gzipped or not (or a mixture). Please no spaces, commas, semicolons in filenames.
 -   Sequences can also be inputted directly in command line (e.g. Filter=AAGAAGAGGG,GAGGAAG,..)
 - Mode = <Positive/Negative>
 -   Positive: Keep only reads in StartFiles with at least one kmer from FilterFiles.
 -   Negative: Keep only reads in StartFiles with no kmers from FilterFiles.

->Optional Parameters
 - KLen=<int> : is the Kmer length (default 32)
 - OutFile=<file> : Allows user to specify filenames for data output (Otherwise defaults filename is assigned) 
 - ReportGranularity=<int> : How often to report progress (default every 100000 reads)

->Other Notes
 - For SRR/ERR/DRR datasets not present locally, KFC will try to download the files from NCBI (linux/mac only)
 - All features should work on Linux/Mac/Windows with Python 3.7+.
 - The pypy interpreter (www.pypy.org) may improve speed up to several fold.
