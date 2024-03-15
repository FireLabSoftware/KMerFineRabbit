#!/usr/bin/env python
## 
## ######################
## KMerFineRabbit ('KFR')-- Reference-indepdendent manipulation of NGS Datasets
## Version as01 03_14_2024
## ->Fuction: KFR Filters an input NGS dataset by keeping/removing reads with k-mer that match a "Filter" dataset
##  - KFR is used to discover reads within a metagenomic dataset that either share homology with a second
##  - dataset or lack such homology.  Datasets can be reference genomes or relatively large fastq/fasta read files
##  - from high throughput sequencing runs.  Matches are defined by setting a k-mer length (default=32, which ensures
##  - that few spurious matches are identified.
##
## ->Syntax
##  - python KMerFineRabbit<ver>.py  Substrate=<file,file..>  Filter=<file,file..>  Mode=<positive/negative>
##
## ->Required Parameters
##  - Substrate=<file[s]> : Starting dataset(s) (reads from these files are filtered to generate your output). 
##  - Filter=<file[s]> : Datasets to gather kmers for filtering Substrate
##  -   File lists can be comma separated (no spaces) and/or specified with a wildcard (*) [wildcards require quotes around name]
##  -   Pairs of "Substrate" files (e.g.,R1,R2) with paired reads should be separated by a semicolon (not comma)
##  -   FilterFile entries that are also in the SubstrateFile list will be ignored (not used for filtering)
##  -   Input files can be gzipped or not (or a mixture). Please no spaces, commas, semicolons in filenames.
##  -   Filter sequences can also be inputted directly in command line (e.g. Filter=AAGAAGAGGG,GAGGAAG,..)
##  - Mode = <Positive/Negative>
##  -   Positive: Keep only reads in StartFiles with at least one kmer from FilterFiles.
##  -   Negative: Keep only reads in StartFiles with no kmers from FilterFiles.
##
## ->Optional Parameters
##  - klen=<32> : is the Kmer length (default 32, realistically can be between 16 and read-length)
##  - SearchCadence=<1> : Setting this to 1 searches every k-mer when filtering the substrate file.  Setting to a higher value n searches only every nth kmer
##  - NoNs=<False> : Setting this to true ignores any read with an N
##  - MinReadLen=<default = klen1> : Save only this-length-or-longer reads in output. For paired reads, both must be at least this long. default=klen1  
##  - OutFile=<file name assigned by program if not specified> : Allows user to specify base filenames for data output 
##  - ReportGranularity=<100000> : How often to report progress (default every 100000 reads)
##  - MaxReadsPerFile1=<0> : Max number to check for each substrate file (default is zero (all reads)
##  - GzipOut = <True/False>: Gzips the output (default is false, can also be set by providing an output file with extension '.gz')
##  - CircularFilter=<False> : Setting this to true will instruct KFR to treat filter files as circles (default is false)
##  -
## ->Additional Feature Options
##  - ***Trimming Function-- Option to trim Reads as they are saved.
##  -   Two types of trimming are used
##  -    (a) To avoid reading into linkers, any read that contains the complement of the first k-mer (of len TrimFk1) from its paired read will be trimmed to remove all sequences after then end of that k-mer
##  -    (b) Any read that contains a k-mer of len TrimTk1 matching a linker will be trimmed to remove all sequences from the start of that k-mer
##  -   TrimOnTheFly=<False>: setting this to True turns on the trimming on the fly option
##  -   TrimFk=<16> : sets a value for k-mer length to avoid paired reads that extend into linker
##  -   TrimTk=<13>: sets a value for k-mer to remove tetritis
##  -   TrimTf=<illuminatetritis1223.fa> provides a FastA file name to extract k-mers for tetritis removal. (a file  that can be used for this is available on the GITHub site)

##  - ***Multiplicity Restriction Option (restricts k-mers for filtering based on their incidence in a separate genome or sequencing dataset file)
##  -   MultiplicityFilterFile=<mff_file> : KFR will restrict filtering Kmers by multiplicity in this file
##  -   MultiplicityFilterMin=<int> : Restrict filtering Kmers to those whose multiplicity in <mff_file> is >= this number
##  -   MultiplicityFilterMax=<int> : Restrict filtering Kmers to those whose multiplicity in <mff_file> is <= this number
##  - 
##  - ***Aggregation Option (default is not to aggregate, so the following can be ignored unless you are using this experimental feature)
##  -   AggregateOn = <False>:  Experimental tool that aggregates k-mers as the program runs. 0 [no aggregation is the default].  Only designed to run in 'positive' mode.
##  -   In aggregation mode, the program collects k-mers in reads that already match the query and includes them in subsequent searching, thereby broadening the search as it goes
##  -   A set of filters is generally required to remove linkers and other junk and these are configured by the parameters below
##  -   Aggregation runs are determinative but can be very messy.  At best they extract sequences related to a query not only through shared K-mers
##  -    but also through bridging reads that share some k-mers with a query and others with the substrate
##  -    Issues with this approach arise if a few k-mer bridges cause KMerFineRabbit to "walk" into a space of sequences that are highly represented in the dataset
##  -    These 'spurious' bridges can lead to a substatial background-- they can be artefactual chimeric reads or true shared k-mers between a major and a minor species in the underlying biological sample
##  -    Either way, such a bridge can result in a rapid "explosion" of aggregated sequences as one or more major components of the initial
##  -    species mixture quickly becomes a target for aggregation.  A periodic readout of captured sequence count and current k-mers being sought will show any such explosion
##  -    The key parameters for designing an Aggregation run are AggregateIt and AggregateMinBridge
##  -     Iteration Number: Aggregation is generally a multi-pass process.  AggregateIt values of >1 instruct KFRto repeatedly go the the dataset aggregating k-mers through each iteration
##  -     Minimun Bridge Number: Default is to require two independent reads (bridges) that link a new k-mer to a "positive read" to include the new k-mer as part of the query set. This can be too sensitive, resulting in frequent explosions of the k-mer search set. Setting AggregateMinBridge to a value of n insists that any k-mer to be included in the query pool must have been identified in at least n independent contexts as associated with query k-mers to itself be considered a new query k-mer.
##  -    Full set of options for Aggregation Utility are as follows:
##  -      AggregateTb  (default is 0)         ## Provides a trim length away from nearest TetritisKmer in adding k-mers to the list
##  -      AggregateRt (default 20000)         ## A rarefication threshold.  Any k-mer above 1/this-frequency in a sample dataset will be ignored  -- so that for value of 20000 any k-mer with frequency <1/20000 kmers will be removed from aggregation)
##  -      AggregateRk (default 16)         ## kmer length for rarefication threshold during aggregation
##  -      AggregateMaxK default 999999999)    ## Maximum number of Kmers used for aggregating reads.  Default setting 999999999 will keep going until the program runs out of memory
##  -      AggregateIt (default 1)             ## Number of pre-production iterations for aggregation (1 or 2 should be fine, but can add more)
##  -      AggregateMultiplier (default 4)     ## Dgree of sampling to obtain a list of k-mers that are insufficiently rare.  Can provide a larger number but would result in only minor advantages
##  -      AggregateMinBridge (default 2)      ## Minimum bridging number (number of different read combinations) to include a k-mer as part of aggregation 
##  -      AggregateThroughProductionRound (default True)   ##  Keeps aggregating even during the production round 
##  -      AggregateMaxK (default 999999999)   ## Maximum number of Kmers used curing aggregation of reads.  Default setting essentially will keep going until the program runs out of memory
##  -      AggregateCadence (default 1)        ## Setting this to a higher number allows aggregation that checks only every nth k-mer
##  -      PreScreen (default False)           ## Setting this to true applies a text-based terminal (start of R1 and R2) prescreen in an attempt to speed up aggregation
##  -      PreScreenK (default 16)             ## K mer size used to prescreen reads for speedier aggregation
##  -      PreScreenProduction (default False) ## Setting this to true includes the "production" run for prescreening.  This may lose some damaged reads but will speed up the production run
##  - *** Other sundry features
##  -    Memory management
##  -     If Substrate complexity is extensive (>500M different kmers), while Filter is a smaller dataset, then
##  -     setting PrimaryIndex1 = 'Filter' can avoid memory issues (avoids preassembling a preliminary kmer Set from Substrates)
##  -     PrimaryIndex = 'filter' turns off Turbo options
##  -      PrimaryIndex=<assigned by program by default>: This relates to memory management and is assigned by the program and shouldn't normally need to be set manually.  It can be set to "Substrate" or "Filter" (whichever has fewer k-mers) if there are out-of-memory errors which in some cases can overcome insufficient memory situations
##  -    Turbo Option
##  -      Turbo (default is False)             ## Turbo=True speeds (~40%) filtering of a modest substrate dataset with huge filter datasets, otherwise not useful
##  -      FilterScale1 (default = 28)          ## This shouldn't need to be changed but is an option to potentially optimize speed for large datasets.  FilterScale controls memory buffer used for Turbo mode. Total bytes used is 2**FilterScale/8 [so 32MB for FilterScale=28 but may be larger due to python overhead)
##
## ->Other Notes
##  - For SRR/ERR/DRR datasets not present locally, KFR will try to download the files from NCBI (linux/mac only).
##  -   You can provide the location of fastq-dump with FastQDumpProgram= (otherwise, KFR will look for it, which can be slow)
##  - All other features should work on Linux/Mac/Windows.
##  - An (optional) parameter dealing with filter architectory PrimaryIndex:
##  -   This parameter determines whether the main Kmer index is from 'filter' or 'substrate' file
##  -   KFR will guess which will be most effective based on the smaller of the two files (the choice will be indicated)
##  -   But the "guess" is not perfect and in some cases of memory limitation you will be able to avoid an "out of memory" error
##  -     by setting PrimaryIndex='substrate' or PrimaryIndex='filter'
##  - Requires Python 3.7+.
##  - The pypy interpreter (www.pypy.org) may improve speed up to several fold.
## ###############
## End Help

klen1 = 32
SubstrateFiles1 = ''
FilterFiles1 = ''
PositiveMode1 = 'default'
OutFileName1 = 'default'
ReportGranularity1 = 100000
CircularFilter1 = False
NoNs1 = False ## True filters out all reads with an 'N' anywhere in the read
MinReadLen1 = -1 ## Save only this-length-or-longer reads in final output. For paired reads, both must be at least this long.  Defaults to klen1 if not set.    
MaxReadsPerFile1 = 0 ## Set to a positive number to limit the number of reads per file

## Some additional user-specifiable values
PairedInput1 = 'default'  ## Are input files paired as R1/R2 (default will autodetect)
FastQDumpProgram1 = 'fasterq-dump' ## Full path of fastq-dump/fasterq-dump program to download from SRA if needed.
GzipOut1 = 'False' ## Setting this to true uses a 'gzip' format for output
## Some routines for multithreading== My suggestion here is not to use this unless you are desparate for additional speed
##  - Threads=<int> : Max number of processor threads to be dedicated to the activity (default = 1)
##  -   If FilterSets have already been assembled, they can be used as filters by specifying '.pck' files
Threads1 = 1
## Memory issues will need to be carefully managed here or your system can hang due to too many processes fighting over memory
WriteSet1 = False  ## Instructs this instance of KFR to write a FilterSet and stop
MyThread1 = 0             ## For use if KFR is running in multithreaded mode: not implemented in current KFR
Turbo1 = False            ## Turbo=True speeds (~40%) filtering of a modest substrate dataset with huge filter datasets
FilterScale1 = 28         ## FilterScale controls memory buffer used for Turbo mode (30 uses 125MB)
PrimaryIndex1 = 'default' ## Memory management for larger files.  PrimaryIndex can be from "Substrate" or "Filter"
MultiplicityFilterFile1 = False ## Specifies file for (optional) multiplicity restriction of filtering k-mers
MultiplicityFilterMin1 = 0 ## Specifies lower cutoff for multiplicity restriction of filtering k-mers
MultiplicityFilterMax1 = 999999999999 ## Specifies upper cutoff for multiplicity restriction of filtering k-mers
SearchCadence1 = 1        ##Setting this to 1 searches every k-mer when filtering the substrate file.  Setting to a higher value n searches only every nth kmer


TrimOnTheFly1 = False ##  Setting this to True causes KFR to trim reads when saving
TrimFk1 = 16     ##  sets a value for k-mer length to avoid paired reads that extend into linker
TrimTk1 = 13     ##  sets a value for k-mer to avoid tetritis
TrimTf1 = ''     ##  'illuminatetritis1223.fa'     ##  provides a FastA file to extract k-mers for tetritis removal

AggregateOn1 = False  ##  Experimental tool that aggregates k-mers as the program runs. 0 [no aggregation] is the default.
AggregateTb1 = 0      ##  provides a trim length away from nearest TetritisKmer in adding k-mers to the list
AggregateRt1 = 20000  ##  a rarefication threshold.  Any k-mer above this frequency in a sample dataset (10*AggregateRt kmers) will be ignored
AggregateRk1 = 16     ##  kmer length for rarefication threshold during aggregation
AggregateIt1 = 1      ##  number of pre-production iterations for aggregation (1 or 2 should be fine, but can add more)
AggregateMultiplier1 = 4  ##  Degree of sampling to obtain a list of k-mers that are insufficiently rare.  Can provide a larger number but would result in only minor advantages
AggregateMinBridge1 = 2   ##  Minimum bridging number (number of different read combinations) to include a k-mer as part of aggregation 
AggregateThroughProductionRound1 = True   ##  Keeps aggregating even during the production round 
AggregateMaxK1 = 999999999 ## Maximum number of Kmers used for aggregating reads.  Default setting essentially will keep going until the program runs out of memory
AggregateCadence1 = 1 ## Setting this to a higher number allows aggregation that checks only every nth k-mer
PreScreen1 =False     ## Setting this to true applies a text-based terminal (start of R1 and R2) prescreen in an attempt to speed up aggregation
PreScreenK1 =16     ## K mer size used to prescreen reads for speedier aggregation
PreScreenProduction1 = False ## Setting this to true includes the "production" run for prescreening.  This may lose some damaged reads but will speed up the production run
from array import array
import gzip
from itertools import chain
from glob import glob
import os,sys
from time import sleep, time, strftime, localtime
import pickle
import subprocess
from collections import Counter
Thread1 = 0
MaxReadLen1 =   200   ##the maximal read length that has been seen
def MakeCheckRange1(Start, Cadence, Max, Offset):
    MyList = [0]*MaxReadLen1
    for i in range(klen1-1+Offset,MaxReadLen1,SearchCadence1):
        MyList[i] = 1
    return MyList
        

def multiglob(x):
    l = []
    for n in x.split(','):
        n = n.strip("'").strip('"')
        l.extend(glob(n))
    return l
        
def LineEnding1(L):
    c = 0
    if len(L)>=1 and L[-1] in '\r\n':
        c+=1
    if len(L)>=2 and L[-2] in '\r\n':
        c+=1
    if c==0:
        return ''
    else:
        return L[-c:]
    
        
for a0 in sys.argv[1:]:
    if a0.startswith('=') or a0.endswith('='):
        print('Encountered equals sign in unexpected position in command line; syntax should be free of spaces')
        exit()
    if a0.startswith('h') or a0.startswith('-h'):
        for L1 in open(sys.argv[0],mode='rt').readlines():
            if L1.startswith('## End Help'): break
            if L1.startswith('##'):
                print(L1.strip('#').strip())
        exit()
    s1 = a0.split('=')[0].strip().strip("'").strip('"')
    v1 = a0.split('=')[1].strip().strip("'").strip('"')
    if s1.lower().startswith('substrate'):
        SubstrateFiles1 =  v1
    elif s1.lower().startswith('filter'):
        FilterFiles1 =  v1
    elif s1.lower().startswith('fastq'):
        FastQDumpProgram1 =  v1
    elif s1.lower().startswith('TrimTf'):
        TrimTf1 =  v1   
    elif s1.lower().startswith('mod'):
        if v1[0].lower() in 'p+ikt':
            PositiveMode1 = True
        elif v1[0].lower() in 'n-erf':
            PositiveMode1 = False
    elif s1.lower().startswith('aggregateon'):
        if v1[0].lower() in 'p+ikt':
            AggregateOn1 = True
        elif v1[0].lower() in 'n-erf':
            AggregateOn1 = False
    elif s1.lower().startswith('prescreenk'):
        PreScreenK1 =  int(v1)
    elif s1.lower().startswith('prescreenproduction'):
        if v1[0].lower() in 'p+ikt':
            PreScreenProduction1 = True
        elif v1[0].lower() in 'n-erf':
            PreScreenProduction1 = False
    elif s1.lower().startswith('prescreen'):
        if v1[0].lower() in 'p+ikt':
            PreScreen1 = True
        elif v1[0].lower() in 'n-erf':
            PreScreen1 = False
    elif s1.lower().startswith('aggregatetwostep'):
        if v1[0].lower() in 'p+ikt':
            AggregateTwoStep1 = True
        elif v1[0].lower() in 'n-erf':
            AggregateTwoStep1 = False
    elif s1.lower().startswith('gzipout'):
        if v1[0].lower() in 'p+ikt':
            GzipOut1 = True
        elif v1[0].lower() in 'n-erf':
            GzipOut1 = False
    elif s1.lower().startswith('trimonthefly'):
        if v1[0].lower() in 'p+ikt':
            TrimOnTheFly1 = True
        elif v1[0].lower() in 'n-erf':
            TrimOnTheFly1 = False
    elif s1.lower().startswith('aggregatethroughproductionround'):
        if v1[0].lower() in 'p+ikt':
            AggregateThroughProductionRound1 = True
        elif v1[0].lower() in 'n-erf':
            AggregateThroughProductionRound1 = False
    elif s1.lower().startswith('report'):
        ReportGranularity1 =  int(v1)
    elif s1.lower().startswith('searchcadence'):
        SearchCadence1 =  int(v1)
    elif s1.lower().startswith('aggregatecadence'):
        AggregateCadence1 =  int(v1)
    elif s1.lower().startswith('maxreadsperfile'):
        MaxReadsPerFile1 =  int(v1)
    elif s1.lower().startswith('minreadlen'):
        MinReadLen1 =  int(v1)
    elif s1.lower().startswith('trimfk'):
        TrimFk1 =  int(v1)
    elif s1.lower().startswith('aggregatemaxk'):
        AggregateMaxK1 =  int(v1)
    elif s1.lower().startswith('aggregateit'):
        AggregateIt1 =  int(v1)
    elif s1.lower().startswith('aggregatemultiplier'):
        AggregateMultiplier1 =  int(v1)
    elif s1.lower().startswith('aggregateminbridge'):
        AggregateMinBridge1 =  int(v1)
    elif s1.lower().startswith('trimtk1'):
        TrimTk1 =  int(v1)
    elif s1.lower().startswith('aggregatetb'):
        AggregateTb1 =  int(v1)
    elif s1.lower().startswith('aggregatert'):
        AggregateRt1 =  int(v1)
    elif s1.lower().startswith('aggregaterk'):
        AggregateRk1 =  int(v1)
    elif s1.lower().startswith('klen'):
        klen1 =  int(v1)        
    elif s1.lower().startswith('threads'):
        Threads1 =  int(v1)        
    elif s1.lower().startswith('out'):
        OutFileName1 =  str(v1)        
    elif s1.lower().startswith('paired'):
        if v1.lower()[0] in ('fn0'):
            PairedInput1 = False
        if v1.lower()[0] in ('ty1'):
            PairedInput1 = True
    elif s1.lower().startswith('turbo'):
        if v1.lower()[0] in ('fn0'):
            Turbo1 = False
        if v1.lower()[0] in ('ty1'):
            Turbo1 = True
    elif s1.lower().startswith('circ'):
        if v1.lower()[0] in ('fn0'):
            CircularFilter1 = False
        if v1.lower()[0] in ('ty1'):
            CircularFilter1 = True
    elif s1.lower().startswith('nons'):
        if v1.lower()[0] in ('fn0'):
            NoNs1 = False
        if v1.lower()[0] in ('ty1'):
            NoNs1 = True
    elif s1.lower().startswith('prim'):
        if v1.lower()[0] in ('f'):
            PrimaryIndex1 = 'filter'
        if v1.lower()[0] in ('s'):
            PrimaryIndex1 = 'substrate'
    elif s1.lower().startswith('write'):
        WriteSet1 = str(v1)
    elif s1.lower().startswith('multiplicityfilterfile'):
        MultiplicityFilterFile1 = str(v1)
    elif s1.lower().startswith('multiplicityfiltermax'):
        MultiplicityFilterMax1 = int(v1)
    elif s1.lower().startswith('multiplicityfiltermin'):
        MultiplicityFilterMin1 = int(v1)
    else:
        print()
        print('--------------------------------------------------------------------------------------------')
        print('  Your key ')
        print('  ->   '+s1)
        print('  in the command line key=value pair')
        print('  ->   '+s1+'='+v1)
        print('  was not recognized.')
        print('  Please try again, or run KMerFineRabbit with keyword "-help"')
        print('  The error is likely severe but will try to continue nonetheless')
        print('--------------------------------------------------------------------------------------------')
        print()
        sleep(6)

t0 = time()
vnow = strftime("%m%d%y_%H%M%S",localtime())
mask1 = 4**(klen1-1)-1
ksam1 = 4**(klen1-1)
Rmask1 = 4**(AggregateRk1-1)-1
Rksam1 = 4**(AggregateRk1-1)
BaseD1 = Counter({'G':0, 'A':1, 'T':2, 'C':3, 'N':0, 'g':0, 'a':1, 't':2, 'c':3, 'n':0, 'U':2, 'u':2})
BaseL1 = [0]*256
for b1 in BaseD1:
    BaseL1[ord(b1)] = BaseD1[b1]

if TrimTf1 == '' and (AggregateOn1 or TrimOnTheFly1) and os.path.isfile('illuminatetritis1223.fa'):
    TrimTf1 = 'illuminatetritis1223.fa'
if (AggregateOn1 or TrimOnTheFly1) and TrimTf1:   
    TrimTf1=multiglob(TrimTf1)
def antisense(s):
    return s.upper().replace('G','c').replace('C','g').replace('A','t').replace('T','a').upper()[::-1]
AggregateD1 = {} ## Keys are k-mers, values are sets of found positions in read (negative for r2)

LogFileName1="LogSummary_"+os.path.basename(sys.argv[0]).split('.')[0]+'.tdf'
if os.path.isfile(LogFileName1):
    open(LogFileName1,mode='a').write('\n')

def LogNote1(*notes):
    note = ' '.join(map(str,notes)).replace(': ',':')
    if Threads1>1:
        note+='; Thread'+str(myThread1)
    LogFile=open(LogFileName1,mode='a')
    LogFile.write(note+'\t'+'; t='+"{0:.3f}".format(time()-t0)+'\t'+strftime("D_%m_%d_%y_T_%H_%M_%S",localtime())+' \n')
    LogFile.close()
    print(note.split('#')[0].replace('\t',' ').replace('\r','\n').strip(',') + '; t='+"{0:.3f}".format(time()-t0))

LogNote1('Running ',os.path.basename(sys.argv[0]),'with parameters:','\r  '.join(sys.argv),'\r  #Python Version',sys.version)
if MinReadLen1 == -1:
    MinReadLen1 = klen1

S2 = set()
def ListToFileSpec1(myFiles):
    if type(myFiles)==str:
        myFiles = [myFiles,]
    mySpec = ''
    for m in myFiles:
        newSpec = ''
        for c in m:
            if c in '._': break
            if c=='*':
                newSpec+='x'
            if c.isalnum():
                newSpec+=c
        newSpec+='_'
        if not newSpec in mySpec:
            mySpec += newSpec        
    return mySpec[:100]
def findfastqdump(candidate):
    ver = 0
    try:
        ver = subprocess.check_output([candidate,'-V'])
        return candidate
    except:
        pass
    if os.path.isfile('~/FastQDumpLocation.txt'):
        candidate = open('~/FastQDumpLocation.txt', mode='rt').read()
        try:
            ver = subprocess.check_output([candidate,'-V'])
            return candidate
        except:
            pass
    LogNote1('Looking for fasterq-dump-- if this fails, provide a location in the command line (fastqdump=<path>)')
    LogNote1('or reinstall and allow execution (chmod +X <path>)')
    LogNote1('Note finding fast(er)q-dump can take some real time, get a cup of coffee or tea')
    NewCands = find1('fasterq-dump')
    NewItems = []
    for candidate in NewCands:
        try:
            ver = subprocess.check_output([candidate,'-V'])
            NewItems.append([versioner1(ver),candidate])
        except:
            pass
    if NewItems:
        candidate = sorted(NewItems, reverse=True)[0][1]
        try:
            open('~/FastQDumpLocation.txt', mode='w').write(candidate)
        except:
            pass
        return candidate
    LogNote1('Unable to find fast-q dump.  Recommend that you reinstall this from ncbi or download fastq/fasta files directly')
    return '' 
def FileListMnemonic1(FL):
    mn  = ''
    FLx = [os.path.basename(x) for x in FL]
    for i,c in enumerate(FLx[0]):
        for n in FLx:
            if len(n)<=i or n[i]!=c or n[i]=='.':
                if len(FL)>1:
                    return mn+'x'
                else:
                    return mn
        mn+=c
    return mn

PreScreenS1 = set()
PreScreenK2 = 2*PreScreenK1
def PreScreenAdd1(s):
    a = antisense(s)
    for i in range(len(s)-PreScreenK1+1):
        PreScreenS1.add(s[i:i+PreScreenK1])
        PreScreenS1.add(a[i:i+PreScreenK1])
        

SubstrateFileList1 = []
for s11 in SubstrateFiles1.split(','):
    if not(s11): continue
    if os.path.isfile(s11):
        SubstrateFileList1.append(s11)
    elif list(glob(s11)):
        SubstrateFileList1.extend(sorted(list(glob(s11))))
    elif ';' in s11:
        (sf1,sf2) = s11.split(';',1)
        if os.path.isfile(sf1):
            SubstrateFileList1.extend([sf1,sf2])
        else:
            SubstrateFileList1.extend(sorted(list(glob(sf1)+list(glob(sf2)))))
    if not(SubstrateFileList1):
        SubstrateFileList1.append(s11)
if not(SubstrateFileList1):
    try:
        History1 = open(os.path.basename(sys.argv[0]).split('.')[0]+'_recents.txt',mode='rt').readlines()
        SubstrateFileList1 = History1[-2].strip().split(',')
        if PositiveMode1 == 'default':
            PositiveMode1 = bool(History1[-1].strip().split('=')[-1])
        LogNote1('No NGS dataset files found as substrate, so trying default from history')
        LogNote1('SubstrateFile(s)=',SubstrateFileList1)
    except:
        pass
if not(SubstrateFileList1):   
    LogNote1('No NGS dataset files found as substrate [you must specify Substrate=<fasta or fastq files>]')
    exit()
SubstrateFileD1 = dict.fromkeys(SubstrateFileList1)

for Fn1 in list(SubstrateFileD1.keys()):
    if Fn1[1:3].lower()=='rr' and not(os.path.isfile(Fn1)) and not('.' in Fn1) and Fn1[:3].isalpha() and Fn1[3:].isdigit():
        if not(os.path.isfile(Fn1+'_1.fastq')):
            LogNote1(Fn1+" looks like a non-fasta, non-fastq filename; will assume it's an NCBI SRA link and try to download")
            FastQDumpProgram1 = findfastqdump(FastQDumpProgram1)
            TryFastQDump1 = subprocess.check_output([FastQDumpProgram1,Fn1])
            LogNote1("Result of "+Fn1+" NCBI Download " +str(TryFastQDump1))
        del(SubstrateFileD1[Fn1])
        if os.path.isfile(Fn1+'_1.fastq'):
            SubstrateFileD1[Fn1+'_1.fastq'] = 0
        else:
            LogNote1('Failed to download',Fn1,'- Will continue for now though')
        if os.path.isfile(Fn1+'_2.fastq'):
            SubstrateFileD1[Fn1+'_2.fastq'] = 0
SubstrateFileList1 = list(SubstrateFileD1.keys())

if PairedInput1=='default':
    if len(SubstrateFileList1)%2==0 and all((SubstrateFileList1[i].replace('_1.','_2.').replace('_R1.','_R2.').replace('_R1_','_R2_')==SubstrateFileList1[i+1].replace('_1.','_2.').replace('_R1.','_R2.').replace('_R1_','_R2_') for i in range(0,len(SubstrateFileList1),2))):
        PairedInput1==True
    elif ';' in ''.join(SubstrateFiles1):
        SubstrateFileList1 = sum([x.split(';') for x in SubstrateFileList1],[])
        PairedInput1==True
    else:
        PairedInput1 = False

FilterFileList1 = []
for s11 in FilterFiles1.split(','):
    if os.path.isfile(s11):
        FilterFileList1.append(s11)
    elif list(glob(s11)):
        FilterFileList1.extend(sorted(list(glob(s11))))
    elif ';' in s11:
        (sf1,sf2) = s11.split(';',1)
        if os.path.isfile(sf1):
            FilterFileList1.extend([sf1,sf2])
        else:
            FilterFileList1.extend(sorted(list(glob(sf1)+list(glob(sf2)))))
    if not(FilterFileList1):
        FilterFileList1.append(s11)
FilterFileList1 = [x for x in FilterFileList1 if not(x in SubstrateFileList1)] ## Substrate files are removed from filter list
if not(FilterFileList1):
    LogNote1('No NGS dataset files found as filter [you must specify Filter=<fasta or fastq files>]')
    exit()

FilterFileD1 = dict.fromkeys(FilterFileList1)
FilterFileList1 = list(FilterFileD1.keys())
for Fn1 in list(FilterFileD1.keys()):
    if Fn1[1:3].lower()=='rr' and not(os.path.isfile(Fn1)) and not('.' in Fn1) and Fn1[:3].isalpha() and Fn1[3:].isdigit():
        if not(os.path.isfile(Fn1+'_1.fastq')):
            LogNote1(Fn1+" looks like a non-fasta, non-fastq filename; will assume it's an NCBI SRA link and try to download")
            FastQDumpProgram1 = findfastqdump(FastQDumpProgram1)
            TryFastQDump1 = subprocess.check_output([FastQDumpProgram1,Fn1])
            LogNote1("Result of "+Fn1+" NCBI Download " +str(TryFastQDump1))
        del(FilterFileD1[Fn1])
        if os.path.isfile(Fn1+'_1.fastq'):
            FilterFileD1[Fn1+'_1.fastq'] = 0
        else:
            LogNote1('Failed to download',Fn1,'- Will continue for now though')
        if os.path.isfile(Fn1+'_2.fastq'):
            FilterFileD1[Fn1+'_2.fastq'] = 0
FilterFileList1 = list(FilterFileD1.keys())

if PositiveMode1 == 'default':
    PositiveMode1 = False
if PrimaryIndex1=='filter' and Turbo1:
        LogNote1("Turbo mode incompatible with PrimaryIndex='filter', turning off Turbo")
        Turbo1 = False
if Turbo1:
    A0 = array('Q',[0]*((2**FilterScale1)//64))
    mask3 = 63
    bshift0 = klen1*2-FilterScale1+6
    PrimaryIndex1='substrate'
def myGetSize(Fn):
    if os.path.isfile(Fn):
        return os.path.getsize(Fn)
    else:
        return len(Fn)
TotalFilterSize1 = sum(map(myGetSize,FilterFileList1))
LogNote1('Total Filter File Size =',TotalFilterSize1)
TotalSubstrateSize1 = sum(map(os.path.getsize,SubstrateFileList1))
LogNote1('Total Substrate File Size =',TotalSubstrateSize1)
if PrimaryIndex1=='default':
    if TotalFilterSize1<TotalSubstrateSize1:
        PrimaryIndex1 = 'filter'
        LogNote1("Primary Index mode set to 'filter' based on file sizes. Rerunning with PrimaryIndex='substrate' might help in some (rare) cases if you get a memory error")
    else:
        PrimaryIndex1 = 'substrate'
        LogNote1("Primary Index mode set to 'substrate' based on file sizes. Rerunning with PrimaryIndex='filter' might help in some (rare) cases if you get a memory error")
if OutFileName1 == 'default':
    Chaser1 = '.'+SubstrateFileList1[0].strip('.gz')[::-1].split('.',1)[0][::-1]
    if PositiveMode1:
        pm1 = '-keepKfrom-'
    else:
        pm1 = '-loseKfrom-'
    OutFileBase1 = FileListMnemonic1(SubstrateFileList1)
    FilterMnemonic1 = FileListMnemonic1(FilterFileList1)
    if PairedInput1:
        OutFileName1 = FileListMnemonic1(SubstrateFileList1[::2])+pm1+FileListMnemonic1(FilterFileList1)+'_'+vnow+'_1'+Chaser1
        OutFileName2 = FileListMnemonic1(SubstrateFileList1[1::2])+pm1+FileListMnemonic1(FilterFileList1)+'_'+vnow+'_2'+Chaser1
    else:
        OutFileName1 = FileListMnemonic1(SubstrateFileList1)+pm1+FileListMnemonic1(FilterFileList1)+'_'+vnow+Chaser1
        OutFileName2 = ''
    if GzipOut1:
        OutFileName1 = OutFileName1+'.gz'
        if OutFileName2:
            OutFileName2 = OutFileName2+'.gz'
            
elif ',' in OutFileName1:
    OutFileName1 = OutFileName1.split(',')[0]
    OutFileName2 = OutFileName1.split(',')[1]
else:
    OutFileName2 = ''
FilterSetList1 = []
FilterDataList1 = []
oldLen1 = 0
if (PrimaryIndex1=='substrate'):
    for F1 in SubstrateFileList1:
        DataCycle1 = 2 ## 2 for fasta, 4 for fastq
        if 'fastq' in os.path.basename(F1.lower()).split('.',1)[-1]:DataCycle1 = 4
        ReportLineGranularity1 = DataCycle1*ReportGranularity1
        if F1.lower().endswith('.gz'):
            F1o = gzip.open(F1, mode='rt')
        else:
            F1o = open(F1, mode='rt')
        for i1,L1 in enumerate(F1o):
            if i1%DataCycle1==1:
                if i1%ReportLineGranularity1==1:
                    LogNote1('Substrate File:',os.path.basename(F1),' Read:',i1//DataCycle1, ' Kmers In Primary:', len(S2))
                v0 = 0
                a0 = 0
                PreScreenAdd1(L1.strip())
                for j1,c1 in enumerate(L1.strip()):
                    cd1 = BaseL1[ord(c1)]
                    v0 = ((v0&mask1)<<2)+cd1
                    a0 = (a0>>2)+(3-cd1)*ksam1
                    if j1>=klen1-1:
                        vamin0 = min(v0,a0)
                        S2.add(vamin0)
                        if Turbo1:
                            APos1 = vamin0>>bshift0
                            APos2 = vamin0&mask3
                            A0[APos1] = A0[APos1] | (1<<APos2)
        F1o.close()
        LogNote1('Finishing Substrate File:',os.path.basename(F1),' Reads:',i1//DataCycle1, ' Kmers In Primary:', len(S2))

            
for f1 in FilterFileList1:
    if f1.endswith('.pck'):
        FilterSetList1.append(f1)
    else:
        FilterDataList1.append(f1)
if Threads1>1:
    FilterDataList1 = FilterDataList1[MyThread::Threads1]
SkipSubstratePreIndex1 = False
if PrimaryIndex1 == 'filter':
    SkipSubstratePreIndex1 = True
    
NotTooCommon1 = set()
if FilterDataList1:    
    for f1,F1 in enumerate(FilterFileList1):
        NB1 = 0
        NR1 = 0
        if 'fastq' in os.path.basename(F1.lower()).split('.',1)[-1]:
            fastq1 = True
        else:
            fastq1 = False
        if F1.lower().endswith('.gz'):
            Fo1 = gzip.open(F1, mode='rt')
        elif os.path.isfile(F1):
            Fo1 = open(F1, mode='rt')
        else:
            Fo1 = ['>CommandLine_'+str(f1)+'_'+F1,F1]
        if fastq1:
            ReportLineGranularity1 = ReportGranularity1*4
            for j1,L1 in enumerate(Fo1):
                if j1%4==1:
                    if j1%ReportLineGranularity1==1:
                        LogNote1('Processing',os.path.basename(F1),'at read',j1//4)
                    v1 = 0
                    a1 = 0
                    s11 = L1.strip()
                    for i1,c1 in enumerate(s11):
                        cd1 = BaseL1[ord(c1)]
                        v1 = ((v1&mask1)<<2)+cd1
                        a1 = (a1>>2)+(3-cd1)*ksam1
                        if PreScreen1: PreScreenAdd1(s11)
                        if (i1>=klen1-1):
                            vamin1 = min(a1,v1)
                            if Turbo1:
                                APos1 = vamin1>>bshift0
                                APos2 = vamin1&mask3
                                if (A0[APos1] & (1<<APos2)):
                                    S2.discard(vamin1)
                            elif SkipSubstratePreIndex1:
                                S2.add(vamin1)
                            else:
                                S2.discard(vamin1)
                    if AggregateOn1 or TrimOnTheFly1:
                        v = 0; a = 0
                        for j,c in enumerate(s11):
                            b = BaseL1[ord(c)]
                            v = ((v&Rmask1)<<2)+b
                            a = (a>>2)+(3-b)*Rksam1
                            if j>AggregateRk1-1:
                                NotTooCommon1.add(min(a,v))
                    NB1 += i1
                    NR1 += 1
        else:
            ReportLineGranularity1 = ReportGranularity1*2
            s1 = []
            for j1,L1 in enumerate(chain(Fo1,['>'])):
                if L1[0]=='>':
                    if s1:
                        if len(s1)>1:
                            MultiLineFastA1 = True
                        else:
                            MultiLineFastA1 = True
                        if MultiLineFastA1 and j1<100:
                            LogNote1('Started processing ',os.path.basename(F1),'segment',n1)
                        v1 = 0
                        a1 = 0
                        s11 = ''.join(s1)
                        LS11 = len(s11)
                        if CircularFilter1:
                            s11+=s11[:klen1-1]
                        for i1,c1 in enumerate(s11):
                            cd1 = BaseL1[ord(c1)]
                            v1 = ((v1&mask1)<<2)+cd1
                            a1 = (a1>>2)+(3-cd1)*ksam1
                            if PreScreen1: PreScreenAdd1(s11)
                            if i1>=klen1-1:
                                vamin1 = min(a1,v1)
                                if Turbo1:
                                    APos1 = vamin1>>bshift0
                                    APos2 = vamin1&mask3
                                    if (A0[APos1] & (1<<APos2)):
                                        S2.discsard(vamin1)
                                elif SkipSubstratePreIndex1:
                                    S2.add(vamin1)
                                else:
                                    S2.discard(vamin1)
                        if AggregateOn1 or TrimOnTheFly1:
                            v = 0; a = 0
                            for j,c in enumerate(s11):
                                b = BaseL1[ord(c)]
                                v = ((v&Rmask1)<<2)+b
                                a = (a>>2)+(3-b)*Rksam1
                                if (j>AggregateRk1-1):
                                    NotTooCommon1.add(min(a,v))
                        if MultiLineFastA1 and j1<100:
                            LogNote1('Finished Processing ',os.path.basename(F1),'segment',n1,'bases',LS11)
                        elif i1%ReportLineGranularity1==0:
                            LogNote1('Processing ',os.path.basename(F1),'segment',n1,'read',i1//ReportLineGranularity1)
                        NB1 += i1
                        NR1 += 1
                    s1 = []
                    n1 = L1.strip()[1:]
                else:
                    s1.append(L1.strip())
        Fo1.close()
    if SkipSubstratePreIndex1:
        LogNote1('Finished Preindex File:',os.path.basename(F1),
              ' file',f1+1,'of',len(FilterDataList1),
              ' File Lines:',j1,
              ' Segments:',NR1,
              ' Bases:',NB1+1,
              ' Kmer Addition:',oldLen1,'to',len(S2))
        oldLen1 = len(S2)
    else:
        LogNote1('Finished File:',os.path.basename(F1),
              ' file',f1+1,'of',len(FilterDataList1),
              ' Lines:',j1+1,
              ' Reads:',NR1+1,
              ' Bases:',NB1+1,
              ' Kmer Reduction:',oldLen1,'to',len(S2))
    
TooCommon1 = set() ## Digital values of k-mers that are too common to use for aggregation
BreakHere1 = set() ## Alphabetical names of k-mers that are in various linkers and thus not used for aggregation (program will stop when these are encountered)
if AggregateOn1 or TrimOnTheFly1:
    LogNote1('Starting Rarefication')
    for F1 in SubstrateFileList1:
        AggregateRCounter1 = Counter()
        kCount1 = 0
        MaxK1 = AggregateRt1*AggregateMultiplier1
        DataCycle1 = 2 ## 2 for fasta, 4 for fastq
        if 'fastq' in os.path.basename(F1.lower()).split('.',1)[-1]:DataCycle1 = 4
        ReportLineGranularity1 = DataCycle1*ReportGranularity1
        if F1.lower().endswith('.gz'):
            F1o = gzip.open(F1, mode='rt')
        else:
            F1o = open(F1, mode='rt')
        for i1,L1 in enumerate(F1o):
            if i1%DataCycle1==1:
                v = 0; a = 0
                for j,c in enumerate(L1.strip()):
                    b = BaseL1[ord(c)]
                    v = ((v&Rmask1)<<2)+b
                    a = (a>>2)+(3-b)*Rksam1
                    if (j>=AggregateRk1-1):
                        AggregateRCounter1[min(v,a)]+=1
                        kCount1 += 1
                if kCount1>=MaxK1:
                    break
        for v1 in AggregateRCounter1:
            if (AggregateRCounter1[v1]>=AggregateMultiplier1) and not(v1 in NotTooCommon1):
                TooCommon1.add(v1)
    LogNote1('Finished Rarefication')
    for FnT1 in TrimTf1:
        if FnT1.upper().endswith('.gz'):
            FT1 = gzip.open(FnT1, mode='rt')
        else:
            FT1 = open(FnT1, mode='rt')
        s1 = ''
        for L1 in chain(FT1,['>@']):
            if L1[0]=='>':
                if s1:
                    s1 = L1.strip()
                    a1 = antisense(s1)
                    for i1 in range(len(s1)-TrimTk1+1):
                        BreakHere1.add(s1[i1:i1+TrimTk1])
                        BreakHere1.add(a1[i1:i1+TrimTk1])
                s1 = ''
            else:
                s1+=L1.strip()            
        FT1.close()

if WriteSet1:
    pf1 = open(WriteSet1, mode='wb')
    pp1 = pickle.Pickler(pf1)
    pp1.dump(S2)
    pf1.close()
    exit()
for pfN0 in FilterSetList1:
    pf0 = open(pfN0, mode='rb')
    pp0 = pickle.Unpickler(pf0)
    MyD1 = pp0.load()
    pf0.close()
    OldLen1 = len(S2)
    S2 = S2.union(MyD1)
    LogNote1('Applied FilterSet Filter:',os.path.basename(pfN0),' StartKCount:',OldLen1,' FilteredKCount',len(S2))

if MultiplicityFilterFile1:
    KMersServed1 = 0
    MultiplicityCounter1 = Counter()
    s1 =[]
    LogNote1('Restricting filter k-mers based on copy number in ',MultiplicityFilterFile1,'.  MultiplicityRange:',MultiplicityFilterMin1,MultiplicityFilterMax1)
    if MultiplicityFilterFile1.lower().endswith('.gz'):
        MFF11 = gzip.open(MultiplicityFilterFile1,mode='rt')
    else:
        MFF11 = open(MultiplicityFilterFile1,mode='rt')
    for L1 in chain(MFF1,['>']):
        if L1[0]=='>':
            if s1:
                v1 = 0
                a1 = 0
                s11 = ''.join(s1)
                LS11 = len(s11)
                if CircularFilter1:
                    s11+=s11[:klen1-1]
                for i1,c1 in enumerate(s11):
                    cd1 = BaseL1[ord(c1)]
                    v1 = ((v1&mask1)<<2)+cd1
                    a1 = (a1>>2)+(3-cd1)*ksam1
                    if i1>=klen1-1:
                        vamin1 = min(a1,v1)
                        if vamin1 in S2:
                            MultiplicityCounter1[vamin1]+=1
                    KMersServed1 += 1
                    if KMersServed1%10000000==0:
                        LogNote1('Checked',str(KMersServed1//1000000)+'M','kmers','in multiplicity check')
            s1 = []
            n1 = L1.strip()[1:]
        else:
            s1.append(L1.strip())
    for s2 in set(S2):
        if not(MultiplicityFilterMin1<=MultiplicityCounter1[s2]<=MultiplicityFilterMax1):
            S2.remove(s2)
    LogNote1('Completed Multiplicity Restriction')
    
if OutFileName1.lower().endswith('.gz'):
    OutFile1 = gzip.open(OutFileName1, mode='wt')
else:
    OutFile1 = open(OutFileName1, mode='w')
if OutFileName2:
    if OutFileName2.lower().endswith('.gz'):
        OutFile2 = gzip.open(OutFileName2, mode='wt')    
    else:
        OutFile2 = open(OutFileName2, mode='w')
else:
    OutFile2 = ''


def AmIAProblemReadPair1(s1,s2):
    if not(AggregateOn1) and not(TrimOnTheFly1):
        return False
    for s in s1,s2:
        v = 0; a = 0
        for j,c in enumerate(s):
            b = BaseL1[ord(c)]
            v = ((v&Rmask1)<<2)+b
            a = (a>>2)+(3-b)*Rksam1
            if (j>AggregateRk1-1) and min(v,a) in TooCommon1:
                return True
    return False
LastIter1 = 0
if AggregateOn1:
    LastIter1 += AggregateIt1
iter1 = 0
CurrentCadence1 = AggregateCadence1
CheckRange1 = MakeCheckRange1(klen1,CurrentCadence1,MaxReadLen1,0)
SeekingK1 = len(S2)
DoneAggregating1 = False
if OutFile2:
    while iter1<=LastIter1:
        iter1+=1
        NR1 = 0
        TR1 = 0
        NK1 = 0 ## final number of written reads
        if iter1 > LastIter1:
            CurrentCadence1 = SearchCadence1
            CheckRange1 = MakeCheckRange1(klen1,CurrentCadence1,MaxReadLen1,0)
            LogNote1('Running production cycle to save reads')
        else:
            LogNote1('Running with aggregate iteration '+str(iter1))
        for F1,F2 in zip(SubstrateFileList1[::2],SubstrateFileList1[1::2]):
            DataCycle1 = 2
            if 'fastq' in os.path.basename(F1.lower()).split('.',1)[-1]:DataCycle1 = 4
            ReportLineGranularity1 = DataCycle1*ReportGranularity1
            if F1.lower().endswith('.gz'):
                F1o = gzip.open(F1, mode='rt')
                F2o = gzip.open(F2, mode='rt')
            else:
                F1o = open(F1, mode='rt')
                F2o = open(F2, mode='rt')
            myFileMnemonic1 = F1o.name.split('_')[0].split('.')[0]            
            for i1,(L1,L2) in enumerate(zip(F1o,F2o)):
                if i1==0:
                    myLineEnding1 = LineEnding1(L1)
                if i1%DataCycle1 == 0:
                    if MaxReadsPerFile1 and DataCycle1*MaxReadsPerFile1==i1: break
                    accu1 = []
                    accu2 = []
                    TR1 += 1
                    KeepRead1 = not(PositiveMode1) ## default is to keep reads for Negative Mode, Lose for Positive Mode
                    if i1%ReportLineGranularity1==0:
                        if iter1<=LastIter1:
                            LogNote1('Substrate Aggregation:',myFileMnemonic1,' Read:',i1//DataCycle1, ' Found:', NR1, ' SeekingKmers:',len(S2))
                        else:
                            LogNote1('Substrate Filtering:',myFileMnemonic1,' Read:',i1//DataCycle1, ' Kept:', NK1, ' SeekingKmers:',len(S2))
                            
                accu1.append(L1)
                accu2.append(L2)
                if i1%DataCycle1 == 1:
                    s1 = L1.strip()
                    s2 = L2.strip()
                    ls1 = len(s1)
                    ls2 = len(s2)
                    if ls1>MaxReadLen1 or ls2>MaxReadLen1:
                        MaxReadLen1 = max(ls1,ls2)
                        CheckRange1 = MakeCheckRange1(klen1-1,CurrentCadence1,MaxReadLen1,0)                    
                    if NoNs1 and (('N' in L1) or ('N' in L2)):
                        KeepRead1 = False
                        continue
                    v0 = 0
                    a0 = 0
                    if PreScreen1 and (iter1<=LastIter1 or PreScreenProduction1):
                        if ls1<PreScreenK2 or ls2<PreScreenK2: continue
                        if not(s1[:PreScreenK1] in PreScreenS1) and not(s1[PreScreenK1:PreScreenK2] in PreScreenS1) and not(s2[:PreScreenK1] in PreScreenS1) and not(s2[PreScreenK1:PreScreenK2] in PreScreenS1): continue
                    for j1,c1 in enumerate(s1):
                        cd1 = BaseL1[ord(c1)]
                        v0 = ((v0&mask1)<<2)+cd1
                        a0 = (a0>>2)+(3-cd1)*ksam1
                        if CheckRange1[j1] and SkipSubstratePreIndex1==(min(v0,a0) in S2):  ## This statement reads out true if exeptions leading to inclusion/exclusion are caught
                            KeepRead1 = PositiveMode1
                            break
                    if KeepRead1!=PositiveMode1:
                        v2 = 0
                        a2 = 0
                        for j2,c2 in enumerate(s2):
                            cd2 = BaseL1[ord(c2)]
                            v2 = ((v2&mask1)<<2)+cd2
                            a2 = (a2>>2)+(3-cd2)*ksam1
                            if CheckRange1[j2] and SkipSubstratePreIndex1==(min(v2,a2) in S2):                            
                                KeepRead1 = PositiveMode1
                                break
                if (i1+1)%DataCycle1==0 and KeepRead1:
                    NR1 += 1
                    if AggregateOn1 or TrimOnTheFly1:
                        readA1 = accu1[1].strip().upper()
                        readA2 = accu2[1].strip().upper()
                        for iA1 in range(len(readA1)-TrimTk1):
                            if readA1[iA1:iA1+TrimTk1] in BreakHere1:
                                readA1 = readA1[:iA1-AggregateTb1]
                                break
                        for iA1 in range(len(readA2)-TrimTk1):
                            if readA2[iA1:iA1+TrimTk1] in BreakHere1:
                                readA2 = readA2[:iA1-AggregateTb1]
                                break
                        p1 = readA1.rfind(antisense(readA2[:TrimFk1]))
                        p2 = readA2.rfind(antisense(readA1[:TrimFk1]))
                        if p1>0 and p2>0:
                            newLen1 = min(p1+TrimFk1,p2+TrimFk1)
                            readA1 = readA1[:newLen1]
                            readA2 = readA2[:newLen1]
                        elif p1>0:
                            readA1 = readA1[:p1+TrimFk1]
                            readA2 = readA2[:p1+TrimFk1]
                        elif p2>0:
                            readA1 = readA1[:p2+TrimFk1]
                            readA2 = readA2[:p2+TrimFk1]
                        if AggregateOn1 and (AggregateThroughProductionRound1 or iter1<=LastIter1) and len(S2)<AggregateMaxK1 and not(AmIAProblemReadPair1(readA1,readA2)):
                            idString1 = accu1[1][:5]+accu2[1][:5]
                            if not('N' in idString1):
                                for rA0 in readA1,readA2:
                                    vA0 = 0; aA0 = 0
                                    for jA0,cA0 in enumerate(rA0):
                                        bD0 = BaseL1[ord(cA0)]
                                        vA0 = ((vA0&mask1)<<2)+bD0
                                        aA0 = (aA0>>2)+(3-bD0)*ksam1
                                        vA1 = min(vA0,aA0)
                                        if (jA0>=klen1-1):
                                            if AggregateMinBridge1>1:
                                                if not(vA1 in AggregateD1):
                                                    AggregateD1[vA1] = [idString1]
                                                elif not(idString1 in AggregateD1[vA1]):
                                                    AggregateD1[vA1].append(idString1)
                                                    if len(AggregateD1[vA1])>=AggregateMinBridge1:
                                                        if PreScreen1: PreScreenAdd1(rA0[jA0-klen1+1:jA0+1])
                                                        S2.add(vA1)
                                                        if len(S2)>=AggregateMaxK1: DoneAggregating1 = True
                                            else:
                                                S2.add(vA1)
                                                if PreScreen1: PreScreenAdd1(rA0[jA0-klen1+1:jA0+1])
                                                if len(S2)>=AggregateMaxK1: DoneAggregating1 = True
                    if iter1>LastIter1:
                        if TrimOnTheFly1 and readA1!=accu1[1].strip():
                            accu1[0] = accu1[0].strip()+'_trim[:'+str(len(readA1))+']'+myLineEnding1
                            accu1[1] = readA1 + myLineEnding1
                            if len(accu1)==4:
                                accu1[3] = accu1[3][:len(readA1)] + myLineEnding1
                        if TrimOnTheFly1 and readA2!=accu2[1].strip():
                            accu2[0] = accu2[0].strip()+'_trim[:'+str(len(readA2))+']'+myLineEnding1
                            accu2[1] = readA2 + myLineEnding1
                            if len(accu2)==4:
                                accu2[3] = accu2[3][:len(readA2)] + myLineEnding1
                        if len(accu1[1].strip())>=MinReadLen1 and len(accu2[1].strip())>=MinReadLen1:
                            NK1 += 1
                            OutFile1.write(''.join(accu1))
                            OutFile2.write(''.join(accu2))
                if iter1 <= LastIter1 and DoneAggregating1: break
            F1o.close(); F2o.close()
            if iter1 <= LastIter1 and DoneAggregating1: break
        if iter1 <= LastIter1 and (DoneAggregating1 or len(S2)==SeekingK1):
            LogNote1('Aggregation stopped at iteration ',iter1,' Kmers @',len(S2),' previous iteration @',SeekingK1)
            iter1 = LastIter1
        SeekingK1 = len(S2)
        if iter1<=LastIter1:
            LogNote1('Finished iteration',str(iter1)+' :'+myFileMnemonic1,' Reads:',TR1, ' Found:', NR1, ' SoughtKmers:',len(S2))
        else:
            LogNote1('Finished Filtering:',myFileMnemonic1,' Reads:',TR1, ' Kept:', NK1, ' SoughtKmers:',len(S2))

else:
    while iter1<=LastIter1:
        iter1+=1
        NR1 = 0
        TR1 = 0
        if iter1 > LastIter1:
            CurrentCadence1 = SearchCadence1
            CheckRange1 = MakeCheckRange1(klen1,CurrentCadence1,MaxReadLen1,0)
            LogNote1('Running production cycle to save reads')
        else:
            LogNote1('Running with aggregate iteration '+str(iter1))
        for F1 in SubstrateFileList1:
            DataCycle1 = 2
            if 'fastq' in os.path.basename(F1.lower()).split('.',1)[-1]:DataCycle1 = 4
            ReportLineGranularity1 = DataCycle1*ReportGranularity1
            if F1.lower().endswith('.gz'):
                F1o = gzip.open(F1, mode='rt')
            else:
                F1o = open(F1, mode='rt')
            myFileMnemonic1 = F1o.name.split('_')[0].split('.')[0] 
            for i1,L1 in enumerate(F1o):
                if i1==0:
                    myLineEnding1 = LineEnding1(L1)
                if i1%DataCycle1 == 0:
                    if MaxReadsPerFile1 and DataCycle1*MaxReadsPerFile1==i1: break
                    accu1 = []
                    TR1 += 1
                    KeepRead1 = not(PositiveMode1)
                    if i1%ReportLineGranularity1==0:                    
                        if iter1<=LastIter1:
                            LogNote1('Substrate Aggregation:',myFileMnemonic1,'(R1only). Read:',i1//DataCycle1, ' Found:', NR1, ' SeekingKmers:',len(S2))
                        else:
                            LogNote1('Substrate Filtering:',myFileMnemonic1,'(R1only). Read:',i1//DataCycle1, ' Kept:', NK1, ' SeekingKmers:',len(S2))              
                accu1.append(L1)
                if i1%DataCycle1 == 1:
                    if NoNs1 and ('N' in L1):
                        KeepRead1 = False
                        continue
                    v0 = 0
                    a0 = 0
                    s1 = L1.strip()
                    ls1 = len(L1)
                    if ls1>MaxReadLen1:
                        MaxReadLen1 = ls1
                        CheckRange1 = MakeCheckRange1(klen1-1,CurrentCadence1,MaxReadLen1,0)
                    for j1,c1 in enumerate(s1):
                        cd1 = BaseL1[ord(c1)]
                        v0 = ((v0&mask1)<<2)+cd1
                        a0 = (a0>>2)+(3-cd1)*ksam1
                        if CheckRange1[j1] and SkipSubstratePreIndex1==(min(v0,a0) in S2):
                            KeepRead1 = PositiveMode1
                            break
                if (i1+1)%DataCycle1==0 and KeepRead1:
                    NR1 += 1
                    if AggregateOn1 or TrimOnTheFly1:
                        readA1 = accu1[1].strip().upper()
                        for iA1 in range(len(readA1)-TrimTk1):
                            if readA1[iA1:iA1+TrimTk1] in BreakHere1:
                                readA1 = readA1[:iA1-AggregateTb1]
                                break
                        if AggregateOn1 and (AggregateThroughProductionRound1 or iter1<=LastIter1) and (len(S2)<AggregateMaxK1) and not(AmIAProblemReadPair1(readA1,'')):
                            idString1 = accu1[1][:10]
                            if not('N' in idString1):
                                vA0 = 0; aA0 = 0
                                for jA0,cA0 in enumerate(readA1):
                                    bD0 = BaseL1[ord(cA0)]
                                    vA0 = ((vA0&mask1)<<2)+bD0
                                    aA0 = (aA0>>2)+(3-bD0)*ksam1
                                    vA1 = min(vA0,aA0)
                                    if (jA0>=klen1-1):
                                        if AggregateMinBridge1>1:
                                            if not(vA1 in AggregateD1):
                                                AggregateD1[vA1] = [idString1]
                                            elif not(idString1 in AggregateD1[vA1]):
                                                AggregateD1[vA1].append(idString1)
                                                if len(AggregateD1[vA1])>=AggregateMinBridge1:
                                                    S2.add(vA1)
                                                    if PreScreen1: PreScreenAdd1(readA1[jA0-klen1+1:jA0+1])
                                                    if len(S2)>=AggregateMaxK1: DoneAggregating1 = True
                                        else:
                                            S2.add(vA1)
                                            if PreScreen1: PreScreenAdd1(readA1[jA0-klen1+1:jA0+1])
                                            if len(S2)>=AggregateMaxK1: DoneAggregating1 = True
                    if iter1>LastIter1:
                        if TrimOnTheFly1 and readA1!=accu1[1].strip():
                            accu1[0] = accu1[0].strip()+'_trim[:'+str(len(readA1))+']'+myLineEnding1
                            accu1[1] = readA1 + myLineEnding1
                            if len(accu1)==4:
                                accu1[3] = accu1[3][:len(readA1)] + myLineEnding1
                        if len(accu1[1].strip())>=MinReadLen1:
                            NK1 += 1
                            OutFile1.write(''.join(accu1))
                if iter1<=LastIter1 and DoneAggregating1: break
            F1o.close()
            if iter1<=LastIter1 and DoneAggregating1: break
        if iter1 <= LastIter1 and (DoneAggregating1 or len(S2)==SeekingK1 or iter1==AggregateIt1):
            LogNote1('Aggregation stopped at iteration ',iter1,' Kmers @',len(S2),' previous iteration @',SeekingK1)
            iter1 = LastIter1
        SeekingK1 = len(S2)
        if iter1<=LastIter1:
            LogNote1('Finished iteration',str(iter1)+' :'+myFileMnemonic1,' Reads:',TR1, ' Founds:', NR1, ' SoughtKmers:',len(S2))
        else:
            LogNote1('Finished Filtering:',myFileMnemonic1,' Reads:',TR1, ' Kepts:', NK1, ' SoughtKmers:',len(S2))
LogNote1('KFR Filter Finished.  Wrote ',NK1,' filtered reads from ',TR1,' total reads, to file(s) ',OutFileName1, OutFileName2)
OutFile1.close()
if OutFile2:
    OutFile2.close()
HistoryFile1 = open(os.path.basename(sys.argv[0]).split('.')[0]+'_recents.txt',mode='a')
HistoryFile1.write('\r'+' '.join(sys.argv)+vnow+'\r')
HistoryFile1.write(OutFileName1+','+OutFileName2+'\r')
HistoryFile1.write('PositiveMode='+str(PositiveMode1))
HistoryFile1.close()

## Copywrite 2022-2024 Andrew Fire and Stanford University, All Rights Reserved
## Version Adjustments starting 12-10-22
## 12-10-22 Bug in importing subprocess fixed (program was crashing if fasterq-dump was present)
##   No effect on input/output of this change
## 12-17-22 Version aj02 Many bugs including errors spotted by Karen Artiles and Drew Galls have been fixed
##   Please report other bugs
## 12-17-22 Additional functionality added
##   Sets are now used instead of Dictionaries, increasing memory efficiency
##   A new option Primary=Filter skips the production of a k-mer set from the Substrate.  This can conserve memory
##       if the Substrate set is large (e.g. large metagenomic dataset or complex DNA) and the filter is relatively simple (<1G k-mers)
##   A new setting option "circular=True" will treat all filter fasta files as circles (default is Circular=False)
## 01-05-23 version aj05.  Fixed a bug that prevented inline SRA downloads
## 05-27-23 version al00.  Added automated choice of PrimaryFilter based on smaller file size (substrate or filter)
##        Also added the option to restrict k-mers based on copy number in an additional file
##        As an example, setting MultiplicityFilterFile='hs1.fa' MultiplicityFilterMax1=1 will use only k-mers unique in hs1.fa for filtering

## 01-27-2024 version Adds an experimental aggregation feature that will aggregate kmers in "Guilt-by-association" manner.  S
##        Also added some adjustments to speed up the underlying searches
##        And squshed a few bugs (Thanks to Drew Galls, Karen Artiles)
## 03-14-24 version adds a filter that discards reads shorter than klen1 in final output.  This can be reset to any value (e.g. 2*klen1)

    
    
