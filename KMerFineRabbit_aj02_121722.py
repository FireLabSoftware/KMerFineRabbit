#!/usr/bin/env python -i
## 
## ######################
## KMerFineRabbit ('KFR')-- Reference-indepdendent manipulation of NGS Datasets
##
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
##  - KLen=<int> : is the Kmer length (default 32, realistically can be between 16 and read-length)
##  - OutFile=<file> : Allows user to specify filenames for data output (Otherwise defaults filename is assigned) 
##  - ReportGranularity=<int> : How often to report progress (default every 100000 reads)
##  - CircularFilter=<true/false> : Setting this to true will instruct KFR to treat filter files as circles (default is false)
##
## ->Other Notes
##  - For SRR/ERR/DRR datasets not present locally, KFC will try to download the files from NCBI (linux/mac only)
##  - All features should work on Linux/Mac/Windows with Python 3.7+.
##  - The pypy interpreter (www.pypy.org) may improve speed up to several fold.
## ###############
## End Help

klen1 = 32
SubstrateFiles1 = ''
FilterFiles1 = ''
PositiveMode1 = 'default'
OutFileName1 = 'default'
ReportGranularity1 = 100000
CircularFilter1 = 100000

## Some additional user-specifiable values
PairedInput1 = 'default'  ## Are input files paired as R1/R2 (default will autodetect)
FastQDumpProgram1 = 'fasterq-dump' ## Full path of fastq-dump/fasterq-dump program to download from SRA if needed.

## Some routines for multithreading== My suggestion here is not to use this unless you are desparate for additional speed
##  - Threads=<int> : Max number of processor threads to be dedicated to the activity (default = 1)
##  -   If FilterSets have already been assembled, they can be used as filters by specifying '.pck' files
Threads1 = 1
## Memory issues will need to be carefully managed here or your system can hang due to too many processes fighting over memory
WriteSet1 = False  ## Instructs this instance of KFR to write a FilterSet and stop
MyThread1 = 0             ## For use if KFR is running in multithreaded mode
Turbo1 = False            ## Turbo=True speeds (~40%) filtering of a modest substrate dataset with huge filter datasets
FilterScale1 = 28         ## FilterScale controls memory buffer used for Turbo mode (30 uses 125MB)
PrimaryIndex1 = 'substrate' ## Memory management for larger files.  PrimaryIndex can be from "Substrate" or "Filter" 
##  If Substrate complexity is extensive (>500M different kmers), while Filter is a smaller dataset, then
##  setting PrimaryIndex1 = 'Filter' can avoid memory issues (avoids preassembling a preliminary kmer Set from Substrates)
##  PrimaryIndex1 = 'filter' turns off Turbo options

from array import array
import gzip
from itertools import chain
from glob import glob
import os,sys
from time import sleep, time, strftime, localtime
import pickle
import subprocess
Thread1 = 0

            
        
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
    if s1.lower().startswith('filter'):
        FilterFiles1 =  v1
    if s1.lower().startswith('fastq'):
        FastQDumpProgram1 =  v1
    elif s1.lower().startswith('mod'):
        if v1[0].lower() in 'p+ik':
            PositiveMode1 = True
        elif v1[0] in 'n-er':
            PositiveMode1 = False
    elif s1.lower().startswith('report'):
        ReportGranularity1 =  int(v1)
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
    elif s1.lower().startswith('prim'):
        if v1.lower()[0] in ('f'):
            PrimaryIndex1 = 'filter'
        if v1.lower()[0] in ('s'):
            PrimaryIndex1 = 'substrate'
    elif s1.lower().startswith('write'):
        WriteSet1 = str(v1)

CircularFilter1
t0 = time()
vnow = strftime("%m%d%y_%H%M%S",localtime())
mask1 = 4**(klen1-1)-1
ksam1 = 4**(klen1-1)
BaseD1 = {'G':0, 'A':1, 'T':2, 'C':3, 'N':0, 'g':0, 'a':1, 't':2, 'c':3, 'n':0}

LogFileName1="LogSummary_"+sys.argv[0].split('.')[0]+'.tdf'
if os.path.isfile(LogFileName1):
    open(LogFileName1,mode='a').write('\n')

def LogNote1(*notes):
    note = ' '.join(map(str,notes))
    LogFile=open(LogFileName1,mode='a')
    LogFile.write(note+'\t'+'; t='+"{0:.3f}".format(time()-t0)+'\t'+strftime("D_%m_%d_%y_T_%H_%M_%S",localtime())+' \n')
    LogFile.close()
    print(note.split('#')[0].replace('\t',' ').replace('\r','\n').strip(',') + '; t='+"{0:.3f}".format(time()-t0))

LogNote1('Running ',os.path.basename(sys.argv[0]),'with parameters:','\r  '.join(sys.argv),'\r  #Python Version',sys.version)

if PrimaryIndex1=='filter':
    if Turbo1:
        LogNote1("Turbo mode incompatible with PrimaryIndex='filter', turning off Turbo")
        Turbo1 = False

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
    vLog('Looking for fasterq-dump-- if this fails, provide a location in the command line (fastqdump=<path>)')
    vLog('or reinstall and allow execution (chmod +X <path>)')
    vLog('Note finding fast(er)q-dump can take some real time, get a cup of coffee or tea')
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
    vLog('Unable to find fast-q dump.  Recommend that you reinstall this from ncbi or download fastq/fasta files directly')
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

if Turbo1:
    A0 = array('Q',[0]*((2**FilterScale1)//64))
    mask3 = 63
    bshift0 = klen1*2-FilterScale1+6


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
    if Fn1[1:3].lower()=='rr' and not(os.path.isfile(Fn1)):
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
        SubstrateFileList1 = sum([x.split(';') for x in SubstrateFileList1])
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
    if Fn1[1:3].lower()=='rr' and not(os.path.isfile(Fn1)):
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

if OutFileName1 == 'default':
    Chaser1 = '.'+SubstrateFileList1[0].strip('.gz')[::-1].split('.',1)[0][::-1]
    if PositiveMode1:
        pm1 = '-keepKfrom-'
    else:
        pm1 = '-loseKfrom-'
    OutFileBase1 = FileListMnemonic1(SubstrateFileList1)
    FilterMnemonic1 = FileListMnemonic1(FilterFileList1)
    if PairedInput1:
        OutFileName1 = FileListMnemonic1(SubstrateFileList1[::2])+pm1+FileListMnemonic1(FilterFileList1)+'_'+vnow+Chaser1
        OutFileName2 = FileListMnemonic1(SubstrateFileList1[1::2])+pm1+FileListMnemonic1(FilterFileList1)+'_'+vnow+Chaser1
    else:
        OutFileName1 = FileListMnemonic1(SubstrateFileList1)+pm1+FileListMnemonic1(FilterFileList1)+'_'+vnow+Chaser1
        OutFileName2 = ''
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
                for j1,c1 in enumerate(L1.strip()):
                    cd1 = BaseD1[c1]
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
                        if Threads1>1:
                            LogNote1('Thread',str(MyThread1+1)+'/'+str(Threads1),'processing',os.path.basename(F1),'at read',j1//4)
                        else:
                            LogNote1('Processing',os.path.basename(F1),'at read',j1//4)
                    v1 = 0
                    a1 = 0
                    for i1,c1 in enumerate(L1.strip()):
                        cd1 = BaseD1[c1]
                        v1 = ((v1&mask1)<<2)+cd1
                        a1 = (a1>>2)+(3-cd1)*ksam1
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
                        if MultiLineFastA1:
                            if Threads1>1:
                                LogNote1('Thread',str(MyThread1+1)+'/'+str(Threads1),'Started processing ',os.path.basename(F1),
                              'chromosome',n1)
                            else:
                                LogNote1('Started processing ',os.path.basename(F1),
                              'chromosome',n1)
                        v1 = 0
                        a1 = 0
                        s11 = ''.join(s1)
                        LS11 = len(s11)
                        if CircularFilter1:
                            s11+=s11[:klen1-1]
                        for i1,c1 in enumerate(s11):
                            cd1 = BaseD1[c1]
                            v1 = ((v1&mask1)<<2)+cd1
                            a1 = (a1>>2)+(3-cd1)*ksam1
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
                        if MultiLineFastA1:
                            if Threads1>1:
                                LogNote1('Thread',str(MyThread1+1)+'/'+str(Threads1),'Finished Processing ',os.path.basename(F1),
                              'chromosome',n1,
                              'bases',LS11)
                            else:
                                LogNote1('Finished Processing ',os.path.basename(F1),
                              'chromosome',n1,
                              'bases',LS11)
                        elif i1%ReportLineGranularity1==0:
                            if Threads1>1:
                                LogNote1('Thread',str(MyThread1+1)+'/'+str(Threads1),'Processing ',os.path.basename(F1),
                              'chromosome',n1,
                              'read',i1//ReportLineGranularity1)
                            else:
                                LogNote1('Processing ',os.path.basename(F1),
                              'chromosome',n1,
                              'read',i1//ReportLineGranularity1)
                        NB1 += i1
                        NR1 += 1
                    s1 = []
                    n1 = L1.strip()[1:]
                else:
                    s1.append(L1.strip())
        Fo1.close()
    if SkipSubstratePreIndex1:
        if Threads1>1:
            LogNote1('Thread',str(MyThread1+1)+'/'+str(Threads1), 'Finished File:',os.path.basename(F1),
              ' file',f1+1,'of',len(FilterDataList1),
              ' Lines:',j1,
              ' Reads:',NR1,
              ' Bases:',NB1,
              ' Kmer Addition:',oldLen1,'to',len(S2))
        else:
            LogNote1('Finished File:',os.path.basename(F1),
              ' file',f1+1,'of',len(FilterDataList1),
              ' Lines:',j1,
              ' Reads:',NR1,
              ' Bases:',NB1,
              ' Kmer Addition:',oldLen1,'to',len(S2))
        oldLen1 = len(S2)
    else:
        oldLen1 = len(S2)
        if Threads1>1:
            LogNote1('Thread',str(MyThread1+1)+'/'+str(Threads1), 'Finished File:',os.path.basename(F1),
              ' file',f1+1,'of',len(FilterDataList1),
              ' Lines:',j1,
              ' Reads:',NR1,
              ' Bases:',NB1,
              ' Kmer Reduction:',oldLen1,'to',len(S2))
        else:
            LogNote1('Finished File:',os.path.basename(F1),
              ' file',f1+1,'of',len(FilterDataList1),
              ' Lines:',j1,
              ' Reads:',NR1,
              ' Bases:',NB1,
              ' Kmer Reduction:',oldLen1,'to',len(S2))
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
    S2 = S2 & MyD1
    LogNote1('Applied FilterSet Filter:',os.path.basename(pfN0),' StartKCount:',OldLen1,' FilteredKCount',len(S2))

if OutFileName2:
    OutFile1 = open(OutFileName1, mode='w')
    OutFile2 = open(OutFileName2, mode='w')
else:
    OutFile1 = open(OutFileName1, mode='w')
    OutFile2 = ''
NR1 = 0
TR1 = 0
if OutFile2:
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
        for i1,(L1,L2) in enumerate(zip(F1o,F2o)):
            if i1%DataCycle1 == 0:
                accu1 = []
                accu2 = []
                TR1 += 1
                KeepRead1 = not(PositiveMode1) ## default is to keep reads for Negative Mode, Lose for Positive Mode
                if i1%ReportLineGranularity1==0:
                    LogNote1('Interim Report, Substrate File:',os.path.basename(F1),' Read:',i1//DataCycle1, ' KeptReads:', NR1)
            accu1.append(L1)
            accu2.append(L2)
            if i1%DataCycle1 == 1:
                v0 = 0
                a0 = 0
                for j1,c1 in enumerate(L1.strip()):
                    cd1 = BaseD1[c1]
                    v0 = ((v0&mask1)<<2)+cd1
                    a0 = (a0>>2)+(3-cd1)*ksam1
                    if (j1>=klen1-1) and SkipSubstratePreIndex1==(min(v0,a0) in S2):  ## This statement reads out true if exeptions leading to inclusion/exclusion are caught
                        KeepRead1 = PositiveMode1
                        break
                v1 = 0
                if KeepRead1!=PositiveMode1:
                    v2 = 0
                    a2 = 0
                    for j2,c2 in enumerate(L2.strip()):
                        cd2 = BaseD1[c2]
                        v2 = ((v2&mask1)<<2)+cd2
                        a2 = (a2>>2)+(3-cd2)*ksam1
                        if (j2>=klen1-1) and SkipSubstratePreIndex1==(min(v2,a2) in S2):                            
                            KeepRead1 = PositiveMode1
                            break
            if (i1+1)%DataCycle1==0 and KeepRead1:
                NR1 += 1
                OutFile1.write(''.join(accu1))
                OutFile2.write(''.join(accu2))
        F1o.close(); F2o.close()
else:
    for F1 in SubstrateFileList1:
        DataCycle1 = 2
        if 'fastq' in os.path.basename(F1.lower()).split('.',1)[-1]:DataCycle1 = 4
        ReportLineGranularity1 = DataCycle1*ReportGranularity1
        if F1.lower().endswith('.gz'):
            F1o = gzip.open(F1, mode='rt')
        else:
            F1o = open(F1, mode='rt')
        for i1,L1 in enumerate(F1o):
            if i1%DataCycle1 == 0:
                accu1 = []
                TR1 += 1
                KeepRead1 = not(PositiveMode1)
                if i1%ReportLineGranularity1==0:
                    LogNote1('KFR Filter Progress:',os.path.basename(F1),' Read:',i1//DataCycle1, ' KeptReads:', NR1)
            accu1.append(L1)
            if i1%DataCycle1 == 1:
                v0 = 0
                a0 = 0
                for j1,c1 in enumerate(L1.strip()):
                    cd1 = BaseD1[c1]
                    v0 = ((v0&mask1)<<2)+cd1
                    a0 = (a0>>2)+(3-cd1)*ksam1
                    if (j1>=klen1-1) and SkipSubstratePreIndex1==(min(v0,a0) in S2):
                        KeepRead1 = PositiveMode1
                        break
            if (i1+1)%DataCycle1==0 and KeepRead1:
                NR1 += 1
                OutFile1.write(''.join(accu1))
        F1o.close()
LogNote1('KFR Filter Finished.  Wrote ',NR1,' filtered reads from ',TR1,' total reads, to file ',OutFileName1)
OutFile1.close()
if OutFile2:
    LogNote1('Finishing KFR.  Wrote ',NR1,' filtered reads from ',TR1,' total reads, to file ',OutFileName2)
    OutFile2.close()
HistoryFile1 = open(os.path.basename(sys.argv[0]).split('.')[0]+'_recents.txt',mode='a')
HistoryFile1.write('\r'+' '.join(sys.argv)+vnow+'\r')
HistoryFile1.write(OutFileName1+','+OutFileName2+'\r')
HistoryFile1.write('PositiveMode='+str(PositiveMode1))
HistoryFile1.close()

## Copywrite 2022 Andrew Fire and Stanford University, All Rights Reserved
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


    
    
