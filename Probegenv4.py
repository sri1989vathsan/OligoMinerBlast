import sys
sys.path.insert(1, './Parameters')

from default_parameters import *
from parameters import *
from Bio.SeqUtils import MeltingTemp as mt
import os 
import shutil

# variables = vars().keys()
# vars =[]

# for name in variables:
    # vars.append(name)

import Extractseq
import ProbeGenfunc
import RunBowtie
import Finalizing

folder = './MLBowtie/'
inputFile = folder+'/Input/final26probes.csv'
pathTxShort="./GRCm26/gencode_vM26_transcripts_ShortHeader.fa"
bowtiedb = "./GRCm38/GRCm38"
readoutprobedoc = "./Parameters/Readoutseq.csv"
primerdoc = "./Parameters/Primerseq.csv"

exec ('nn_table = mt.%s' % nntable)

fold = [x[0] for x in os.walk(folder)][2:]

for t in range(0,len(fold)):
    shutil.rmtree(fold[t])


#if __name__=='___main__':
Extractseq.seqextract(inputFile,pathTxShort,folder)

ProbeGenfunc.probegen(l, L, gcPercent, GCPercent, nn_table, tm, TM,
                      X, sal, form, sp, concA, concB, headerVal, bedVal,
                      OverlapModeVal, verbocity, reportVal, debugVal, metaVal,
                      outNameVal,folder)


RunBowtie.bowtierun(bowtiedb,folder)

Finalizing.finalprobe(readoutprobedoc,primerdoc,nprobes,folder)

print("Done")    

#if __name__=='___main__':
#    ProbeGenfunc.probegen(l, L, gcPercent, GCPercent, nntable, tm, TM, 
#                          X, sal, form, sp, concA, concB, headerVal, bedVal,
#                          OverlapModeVal, verbocity, reportVal, debugVal, metaVal,
#                          outNameVal) 