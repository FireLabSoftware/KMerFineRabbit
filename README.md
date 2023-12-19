## KMerFineRabbit ('KFR')
+ Reference-indepdendent manipulation of NGS Datasets
+ Version AN01 10_25_2023
+ ->Fuction: KFR Filters an input NGS dataset by keeping/removing reads with k-mer that match a "Filter" dataset
+  - KFR is used to discover reads within a metagenomic dataset that either share homology with a second
+  - dataset or lack such homology.  Datasets can be reference genomes or relatively large fastq/fasta read files
+  - from high throughput sequencing runs.  Matches are defined by setting a k-mer length (default=32, which ensures
+  - that few spurious matches are identified.
+
+ ->Syntax
+  - python KMerFineRabbit<ver>.py  Substrate=<file,file..>  Filter=<file,file..>  Mode=<positive/negative>
+
+ ->Required Parameters
+  - Substrate=<file[s]> : Starting dataset(s) (reads from these files are filtered to generate your output). 
+  - Filter=<file[s]> : Datasets to gather kmers for filtering Substrate
+  -   File lists can be comma separated (no spaces) and/or specified with a wildcard (*) [wildcards require quotes around name]
+  -   Pairs of "Substrate" files (e.g.,R1,R2) with paired reads should be separated by a semicolon (not comma)
+  -   FilterFile entries that are also in the SubstrateFile list will be ignored (not used for filtering)
+  -   Input files can be gzipped or not (or a mixture). Please no spaces, commas, semicolons in filenames.
+  -   Filter sequences can also be inputted directly in command line (e.g. Filter=AAGAAGAGGG,GAGGAAG,..)
+  - Mode = <Positive/Negative>
+  -   Positive: Keep only reads in StartFiles with at least one kmer from FilterFiles.
+  -   Negative: Keep only reads in StartFiles with no kmers from FilterFiles.
+
+ ->Optional Parameters
+  - KLen=<int> : is the Kmer length (default 32, realistically can be between 16 and read-length)
+  - OutFile=<file> : Allows user to specify filenames for data output (Otherwise defaults filename is assigned) 
+  - ReportGranularity=<int> : How often to report progress (default every 100000 reads)
+  - CircularFilter=<true/false> : Setting this to true will instruct KFR to treat filter files as circles (default is false)
+  - MultiplicityFilterFile=<mff_file> : KFR will restrict filtering Kmers by multiplicity in this file
+  - MultiplicityFilterMin=<int> : Restrict filtering Kmers to those whose multiplicity in <mff_file> is >= this number
+  - MultiplicityFilterMax=<int> : Restrict filtering Kmers to those whose multiplicity in <mff_file> is <= this number
+  - NoNs = <True/False> : Filter out any read with an N
+
+ ->Other Notes
+  - For SRR/ERR/DRR datasets not present locally, KFC will try to download the files from NCBI (linux/mac only)
+  - All other features should work on Linux/Mac/Windows.
+  - An (optional) parameter dealing with filter architectory PrimaryIndex:
+  -   This parameter determines whether the main Kmer index is from 'filter' or 'substrate' file
+  -   KFR will guess which will be most effective based on the smaller of the two files (the choice will be indicated)
+  -   But the "guess" is not perfect and in some cases of memory limitation you will be able to avoid an "out of memory" error by setting PrimaryIndex='substrate' or PrimaryIndex='filter'
+  - Requires Python 3.7+.  The pypy interpreter (www.pypy.org) may improve speed up to several fold.
