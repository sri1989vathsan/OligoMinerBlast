def seqextract(inputFile,pathTxShort,folder):   
    
    #from Bio import SeqIO
    #import blockParse
    # Import Biopython modules.
    #from Bio.SeqUtils import MeltingTemp as mt
    #from Bio.Seq import Seq
    #from Bio.Alphabet import IUPAC
    #from Bio.SeqUtils import GC
    from Bio import SeqIO
    import numpy as np
    import os
    import openpyxl
    import pandas as pd

    # Import regex module.
    import re

    paths = os.path.join(folder,"Sequences")
    
    if not os.path.isdir(paths):
        os.mkdir(paths)

    name =[]
    seq = []
    
    for x in SeqIO.parse(pathTxShort, 'fasta'):
        seq.append(str(x.seq).upper())
        name.append(str(x.name).upper())
    names = np.array(name)
    pos2 =[]
    list1 = pd.read_csv(inputFile)
    list2 = list1.to_numpy()
    for i in range(0,len(list2[:,0])):
        if(isinstance(list2[i,3],str)==True):
            pos2.append("abc")
        else:
            pos = np.argwhere(list2[i,0]==names)
            pos2.append((pos[0][0]))
            
       
    for j in range(0,len(pos2)):
        rna_file = open(paths+"/%s_%s.fa"%(list2[j,2],list2[j,0]), 'w+')
        if pos2[j] !="abc":
            out = '\n'.join(['>%s'%names[int(pos2[j])]+" 100:200" +"\n" + seq[int(pos2[j])]])
        else:
            out = '\n'.join(['>%s'%list2[j,0]+" 100:200" +"\n" + list2[j,3]])
        rna_file.write(out)
        rna_file.close()
            
    print('Finished extracting sequences')
        
if __name__ == '__main__':
    # test1.py executed as script
    # do something
    seqextract()
