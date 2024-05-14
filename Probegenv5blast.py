import sys
sys.path.insert(1, './Parameters')
    
## importing parameters    
from default_parameters import *
from parameters import *
from Bio.SeqUtils import MeltingTemp as mt
import os 
import shutil
import Extractseq
import ProbeGenfuncfasta
import RunBlast
import Finalizingblast

### the script will extract transcript sequences based on the ensembl transcript ID and will find probes with parameters defined in the Parameters folder. THe final probes will be generated with primer and readout sequences at either end (as defined in the Primersseq.csv and Readoutseq.csv)
folder = './MLBlast2/' #folder path with input data
inputFile = folder+'/Input/Panel2.csv' #pahtgene list
pathTxShort="./GRCm26/gencode_vM26_transcripts_ShortHeader.fa" ## list of transcript sequences with ensembl IDs - this can be downloaded from gencode website
blastdb = "./mouse_GencodevM26/TxShortHeader2/gencode_vM26_transcripts_ShortHeader.fa" #blast database created from the list of transcripts
blasttrnadb = "./mouse_GencodevM26/tRNA/mm10tRNAs.fa" # blast database for tRNAs
blastrrnadb = "./mouse_GencodevM26/rRNA/M28_rRNAs.fa" # blast database for rRNAs 
readoutprobedoc = "./Parameters/Readoutseq.csv" # File with readout (secondary probes) sequences
primerdoc = "./Parameters/Primerseq.csv" # File with primer sequences
transcriptlist = "./Parameters/Transcriptlist.txt" #File containing transcripts with ENSEMBL transcript IDs and the corresponding ENSEMBL gene IDs - required for later analysis with blast to remove probes matching to other genes

primernum=1 #

exec ('nn_table = mt.%s' % nntable)

fold = [x[0] for x in os.walk(folder)][2:]

for t in range(0,len(fold)):
    shutil.rmtree(fold[t])

#### Extracting the sequences for individual transcripts
Extractseq.seqextract(inputFile,pathTxShort,folder)

### Generating encoding probe sequences for each transcript matching the parameters in the default_parameters and parameters files
ProbeGenfuncfasta.probegen(l, L, gcPercent, GCPercent, nn_table, tm, TM,
                      X, sal, form, sp, concA, concB, headerVal, bedVal,
                      OverlapModeVal, verbocity, reportVal, debugVal, metaVal,
                      outNameVal,folder)

## run blast against transcriptome, tRNA and rRNA
RunBlast.blastrun(blastdb,blasttrnadb,blastrrnadb,folder)

## finalising - includes adding removing probes matching to multiple genes, adding readout and primer sequences
Finalizingblast.finalprobe(readoutprobedoc,primerdoc,nprobes,folder,transcriptlist,primernum)

print("Done")    