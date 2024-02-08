## KMerFineRabbit ('KFR')-- Reference-indepdendent manipulation of NGS Datasets
+ Version ar11 02_07_2024
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
+  - KLen=<32> : is the Kmer length (default 32, realistically can be between 16 and read-length)
+  - SearchCadence=<1> : Setting this to 1 searches every k-mer when filtering the substrate file.  Setting to a higher value n searches only every nth kmer
+  - NoNs=<False> : Setting this to true ignores any read with an N
+  - OutFile=<file name assigned by program if not specified> : Allows user to specify base filenames for data output 
+  - ReportGranularity=<100000> : How often to report progress (default every 100000 reads)
+  - MaxReadsPerFile1=<0> : Max number to check for each substrate file (default is zero (all reads)
+  - GzipOut = <True/False>: Gzips the output (default is false, can also be set by providing an output file with extension '.gz')
+  - CircularFilter=<False> : Setting this to true will instruct KFR to treat filter files as circles (default is false)
+  -
+ ->Additional Feature Options
+  - ***Trimming Function-- Option to trim Reads as they are saved.
+  -   Two types of trimming are used
+  -    (a) To avoid reading into linkers, any read that contains the complement of the first k-mer (of len TrimFk1) from its paired read will be trimmed to remove all sequences after then end of that k-mer
+  -    (b) Any read that contains a k-mer of len TrimTk1 matching a linker will be trimmed to remove all sequences from the start of that k-mer
+  -   TrimOnTheFly=<False>: setting this to True turns on the trimming on the fly option
+  -   TrimFk=<16> : sets a value for k-mer length to avoid paired reads that extend into linker
+  -   TrimTk=<13>: sets a value for k-mer to remove tetritis
+  -   TrimTf=<illuminatetritis1223.fa> provides a FastA file name to extract k-mers for tetritis removal. (a file  that can be used for this is available on the GITHub site)

+  - ***Multiplicity Restriction Option (restricts k-mers for filtering based on their incidence in a separate genome or sequencing dataset file)
+  -   MultiplicityFilterFile=<mff_file> : KFR will restrict filtering Kmers by multiplicity in this file
+  -   MultiplicityFilterMin=<int> : Restrict filtering Kmers to those whose multiplicity in <mff_file> is >= this number
+  -   MultiplicityFilterMax=<int> : Restrict filtering Kmers to those whose multiplicity in <mff_file> is <= this number
+  - 
+  - ***Aggregation Option (default is not to aggregate, so the following can be ignored unless you are using this experimental feature)
+  -   AggregateOn = <False>:  Experimental tool that aggregates k-mers as the program runs. 0 [no aggregation is the default].  Only designed to run in 'positive' mode.
+  -   In aggregation mode, the program collects k-mers in reads that already match the query and includes them in subsequent searching, thereby broadening the search as it goes
+  -   A set of filters is generally required to remove linkers and other junk and these are configured by the parameters below
+  -   Aggregation runs are determinative but can be very messy.  At best they extract sequences related to a query not only through shared K-mers
+  -    but also through bridging reads that share some k-mers with a query and others with the substrate
+  -    Issues with this approach arise if a few k-mer bridges cause KMerFineRabbit to "walk" into a space of sequences that are highly represented in the dataset
+  -    These 'spurious' bridges can lead to a substatial background-- they can be artefactual chimeric reads or true shared k-mers between a major and a minor species in the underlying biological sample
+  -    Either way, such a bridge can result in a rapid "explosion" of aggregated sequences as one or more major components of the initial
+  -    species mixture quickly becomes a target for aggregation.  A periodic readout of captured sequence count and current k-mers being sought will show any such explosion
+  -    The key parameters for designing an Aggregation run are AggregateIt and AggregateMinBridge
+  -    >Iteration Number: Aggregation is generally a multi-pass process.  AggregateIt values of >1 instruct KFRto repeatedly go the the dataset aggregating k-mers through each iteration
+  -    >Minimun Bridge Number: Default is to require two independent reads (bridges) that link a new k-mer to a "positive read" to include the new k-mer as part of the query set. This can be too sensitive, resulting in frequent explosions of the k-mer search set. Setting AggregateMinBridge to a value of n insists that any k-mer to be included in the query pool must have been identified in at least n independent contexts as associated with query k-mers to itself be considered a new query k-mer.
+  -    Full set of options for Aggregation Utility are as follows:
+  -    >AggregateTb  (default is 0)         ## Provides a trim length away from nearest TetritisKmer in adding k-mers to the list
+  -    >AggregateRt (default 20000)         ## A rarefication threshold.  Any k-mer above 1/this-frequency in a sample dataset will be ignored  -- so that for value of 20000 any k-mer with frequency <1/20000 kmers will be removed from aggregation)
+  -    >AggregateRk (default 16)         ## kmer length for rarefication threshold during aggregation
+  -    >AggregateMaxK default 999999999)    ## Maximum number of Kmers used for aggregating reads.  Default setting 999999999 will keep going until the program runs out of memory
+  -    >AggregateIt (default 1)             ## Number of pre-production iterations for aggregation (1 or 2 should be fine, but can add more)
+  -    >AggregateMultiplier (default 4)     ## Dgree of sampling to obtain a list of k-mers that are insufficiently rare.  Can provide a larger number but would result in only minor advantages
+  -    >AggregateMinBridge (default 2)      ## Minimum bridging number (number of different read combinations) to include a k-mer as part of aggregation 
+  -    >AggregateThroughProductionRound (default True)   ##  Keeps aggregating even during the production round 
+  -    >AggregateMaxK (default 999999999)   ## Maximum number of Kmers used curing aggregation of reads.  Default setting essentially will keep going until the program runs out of memory
+  -    >AggregateCadence (default 1)        ## Setting this to a higher number allows aggregation that checks only every nth k-mer
+  -    >PreScreen (default False)           ## Setting this to true applies a text-based terminal (start of R1 and R2) prescreen in an attempt to speed up aggregation
+  -    >PreScreenK (default 16)             ## K mer size used to prescreen reads for speedier aggregation
+  -    >PreScreenProduction (default False) ## Setting this to true includes the "production" run for prescreening.  This may lose some damaged reads but will speed up the production run
+  - *** Other sundry features
+  -    Memory management
+  -    >If Substrate complexity is extensive (>500M different kmers), while Filter is a smaller dataset, then
+  -    >setting PrimaryIndex1 = 'Filter' can avoid memory issues (avoids preassembling a preliminary kmer Set from Substrates)
+  -    >PrimaryIndex = 'filter' turns off Turbo options
+  -    > PrimaryIndex=<assigned by program by default>: This relates to memory management and is assigned by the program and shouldn't normally need to be set manually.  It can be set to "Substrate" or "Filter" (whichever has fewer k-mers) if there are out-of-memory errors which in some cases can overcome insufficient memory situations
+  -    Turbo Option
+  -    >Turbo (default is False)             ## Turbo=True speeds (~40%) filtering of a modest substrate dataset with huge filter datasets, otherwise not useful
+  -    >FilterScale1 (default = 28)          ## This shouldn't need to be changed but is an option to potentially optimize speed for large datasets.  FilterScale controls memory buffer used for Turbo mode. Total bytes used is 2**FilterScale/8 [so 32MB for FilterScale=28 but may be larger due to python overhead)
+
+ ->Other Notes
+  - For SRR/ERR/DRR datasets not present locally, KFR will try to download the files from NCBI (linux/mac only).
+  -   You can provide the location of fastq-dump with FastQDumpProgram= (otherwise, KFR will look for it, which can be slow)
+  - All other features should work on Linux/Mac/Windows.
+  - An (optional) parameter dealing with filter architectory PrimaryIndex:
+  -   This parameter determines whether the main Kmer index is from 'filter' or 'substrate' file
+  -   KFR will guess which will be most effective based on the smaller of the two files (the choice will be indicated)
+  -   But the "guess" is not perfect and in some cases of memory limitation you will be able to avoid an "out of memory" error
+  -    by setting PrimaryIndex='substrate' or PrimaryIndex='filter'
+  - Requires Python 3.7+.
+  - The pypy interpreter (www.pypy.org) may improve speed up to several fold.
