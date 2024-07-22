GCPercent = 75  ## GC Percent high
gcPercent = 35  ## GC Percent low
TM=80   ## Melting temperature high
tm=45   ## Melting temperature low
form=30  ## Formamide concentration
l=30   ## Length high
L=30   ## Length low
output="3.fastq"
sp=3
bedVal=False
nprobes=50  ## Number of probes
sp2 = 7   ##spacing between probes in case of low number of probes detected - this is between start of previous probe to start of next probe - so overlapping probes
blastrun = 0

X = 'AAAAA,TTTTT,CCCCC,GGGGG'   ## Excluded sequences