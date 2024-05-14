def finalprobe(readoutprobedoc,primerdoc,nprobes,folder):
    import os
    import glob
    import pandas as pd
    import numpy as np
    from Bio.Seq import Seq
    from Bio import SeqIO

    #from Bio.Alphabet import IUPAC
    import subprocess
    import re

    #readoutprobedoc = "/Users/asrivath/Github/OligoMiner/Readoutseq.csv"
    #primerdoc = "/Users/asrivath/Github/OligoMiner/Primerseq.csv"

    primelist1 = pd.read_csv(primerdoc,header=None,sep=",")
    primelist2 = primelist1.to_numpy()

    readlist1 = pd.read_csv(readoutprobedoc,header=None,sep=",")
    readlist2 = readlist1.to_numpy()
    
    seqlistpath = os.path.join(folder,"Sequences/*.fa")
    seqfiles = glob.glob(seqlistpath)
    
    #print(seqfiles)

    path = os.path.join(folder,"InitialProbes/*.sam")
    files = glob.glob(path)
    files.sort()
    
    path2 = os.path.join(folder,"Sequences/*.fa")
    files2 = glob.glob(path2)
    files2.sort()
    
    if not os.path.isdir(folder+"Output"):
        os.mkdir(folder+"Output")

    revcom = []
    rev = []
    rej = []
    lists = []
    num=0

    for fil in range(0,len(files)):
        print(fil)
        list1 = pd.read_csv(files[fil],header=None,sep="\t",usecols=range(0,10))
        list2 = list1.to_numpy()
        seqs = list2[:,9]
        a,b = np.unique(seqs,return_counts=True)
        #print(files2[fil])
        with open(files2[fil]) as handle:
            for seqs in SeqIO.parse(handle,"fasta"):
                lengen = len(seqs.seq)
                
        fina = a[b < 2]
        final2 = a[b>1]
        randlist = []
        if(len(fina)>=nprobes):
            for i in range(10000):
                randnum = np.random.randint(0,len(fina))
                if randnum not in randlist:
                    randlist.append(randnum)
            randlist.sort()
            randlist2 = randlist[0:35]
            final = fina[randlist2]

        else:
            final = fina
        gene = files[fil].split("/")[-1].split("_")[0]
        #print(len(final))
        
        temp1 = readlist2[fil,0].split("{")[1]
        temp2 = temp1.split("}")[0]

        #ll = os.system("grep "+final[0]+" "+files2[fil]+">/dev/null 2>&1")
        try:
            ll = subprocess.check_output("grep "+final[0]+" "+files2[fil],shell=True)
        except Exception:
            ll = None
        for t in range(0,len(final)):
            if(ll!=None):
                temp = str(Seq(final[t]).reverse_complement())
                #print("Yay")
            else:
                temp = str(final[t])
                #print("Nay")
            rev.append([gene+"_RO_"+temp2+"_"+str(t+1), temp])
            revcom.append([gene+"_RO_"+temp2+"_"+str(t+1), primelist2[0,1]+readlist2[fil,1]+temp+readlist2[fil,1]+primelist2[0,2]])
            num=t
            
        for m in range(0,len(final2)):
            temp2 = str(final2[m])
            rej.append([gene+"_"+str(m),temp2])
            
        
        lists.append([gene,str(num+1),lengen])
            
    revcoms = np.array(revcom)
    revs = np.array(rev)
    rej = np.array(rej)
    lists = np.array(lists)


    pd.DataFrame(revcoms).to_csv(folder+"/Output/MultiFISHProbes_Bowtie.csv")
    pd.DataFrame(revs).to_csv(folder+"/Output/BindingProbes_Bowtie.csv")
    pd.DataFrame(rej).to_csv(folder+"/Output/RejectedProbes_Bowtie.csv")
    pd.DataFrame(lists).to_csv(folder+"/Output/NumberProbelist_Bowtie.csv")
    print('Completed Probe generation')

#print(type(revcom))

#print(np.array(revcom).shape)

#print(len(dup))
#list3 = pd.DataFrame(list2[:,9].T).T.drop_duplicates(keep=False).as_matrix()
#print(list3)
#print(len(list3[0]))