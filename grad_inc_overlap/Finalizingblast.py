def finalprobe(readoutprobedoc,primerdoc,nprobes,folder,transcriptlist,primernum):
    import os
    import glob
    import pandas as pd
    import numpy as np
    from Bio.Seq import Seq
    from Bio import SeqIO
    from natsort import natsorted

    #from Bio.Alphabet import IUPAC
    import subprocess
    import re

    #readoutprobedoc = "/Users/asrivath/Github/OligoMiner/Readoutseq.csv"
    #primerdoc = "/Users/asrivath/Github/OligoMiner/Primerseq.csv"

    primelist1 = pd.read_csv(primerdoc,header=None,sep=",")
    primelist2 = primelist1.to_numpy()

    readlist1 = pd.read_csv(readoutprobedoc,header=None,sep=",")
    readlist2 = readlist1.to_numpy()
    
    seqlistpath = os.path.join(folder,"1_Sequences/*.fa")
    seqfiles = glob.glob(seqlistpath)
    seqfiles.sort()
    
    #print(seqfiles)
    path = os.path.join(folder,"2_InitialProbes/*.fa")
    files = glob.glob(path)
    files.sort()
    
    print(folder)
    
    # path2 = os.path.join(folder,"2_InitialProbes/*.txt")
    # files2 = glob.glob(path2)
    # files2.sort()
    
    outputfolder = os.path.join(folder,"3_Output")
    
    if not os.path.isdir(outputfolder):
        os.mkdir(outputfolder)

    revcom = []
    rev = []
    rej = []
    lists = []
    num=0
    transcriptlist1 = pd.read_csv(transcriptlist,sep=",")
    transcriptlist2 = transcriptlist1.to_numpy()

    for fil in range(0,len(files)):
        seq=[]
        name=[]
        for x in SeqIO.parse(files[fil], 'fasta'):
            seq.append(str(x.seq).upper())
            name.append(str(x.name).upper())
        name2 = np.array(name)
        seq2 = np.array(seq)
        
        with open(seqfiles[fil]) as handle:
            for seqs in SeqIO.parse(handle,"fasta"):
                lengen = len(seqs.seq)
        
        #print(lengen)
        
        out=files[fil].split("/")[-1].split(".")[0].split("_")[1]
        pos = np.argwhere(out==transcriptlist2[:,0])
        if(len(pos)>0):
            pos2 = np.argwhere(transcriptlist2[pos,1]==transcriptlist2)
            transcriptlist3 = transcriptlist2[pos2[:,0],0]
        else:
            transcriptlist3 = out
            
        blastfilename = "."+files[fil].split(".")[1]+".txt"
        
        if os.path.exists(blastfilename):
            blastres = pd.read_csv(blastfilename,sep="\t",header=None)
            blastres2 = blastres.to_numpy()
            blastres3 = blastres2[np.argwhere(blastres2[:,11]>40).flatten(),:]
            nongen = np.setdiff1d(blastres3[:,1], transcriptlist3)
            #uniqueprobes = blastres3[coor[:,0],0]
            #blastres3 = blastres2[list(np.argwhere(blastres2[:,11]>50)x,:)]
            coor =[]
            if len(nongen)>0:
                for tt in range(0,len(nongen)):
                    coor.append(np.argwhere(nongen[tt]==blastres3[:,1]))
                    #transcriptlist3[tt(dists >= r) & (dists <= r+dr)]
                    #print(coor)
                coor = np.vstack(coor)
                delprobes = blastres3[coor,0]
                blastname = np.unique(blastres3[:,0])
                coor1 = np.unique(np.argwhere(delprobes==blastname)[:,1])
                probename = np.delete(blastname,coor1,0)
            else:
                probename = np.unique(blastres3[:,0])
        
        else:
            probename=name2

        blasttrnares = pd.read_csv("."+files[fil].split(".")[1]+"tRNA.txt",sep="\t",header=None)
        blasttrnares2 = blasttrnares.to_numpy()
        blasttrnares3 = blasttrnares2[np.argwhere(blasttrnares2[:,11]>30).flatten(),:]

        
        blastrrnares = pd.read_csv("."+files[fil].split(".")[1]+"rRNA.txt",sep="\t",header=None)
        blastrrnares2 = blastrrnares.to_numpy()
        blastrrnares3 = blastrrnares2[np.argwhere(blastrrnares2[:,11]>30).flatten(),:]
        
        coors1 = []
        coors2 = []
        if len(blasttrnares3)>0:
            coors1 = blasttrnares3[:,0] 
        if len(blastrrnares3)>0:
            coors2 = blastrrnares3[:,0] 
        
        probes = np.unique(np.append(coors1, coors2))
        coor3 = []
        if len(probes)>0:
            for rr in range(0,len(probes)):
                coor3.append(np.argwhere(probename==probes[rr]))
            finprobes = np.delete(probename,coor3)
        else:
            finprobes = probename
        
        finprobes = natsorted(finprobes)
        
        fina = []
        randlist=[]
        for r in range(0,len(finprobes)):
            coor4 = np.argwhere(name2==finprobes[r].upper())
            fina.append(seq2[coor4][0][0])
            
        if(len(fina)>nprobes):
            ### if the total number of probes = nprobes should be generated randomly distributed throughout the sequence of the mrna
            if(int(np.vstack((np.char.split(name2,"_")))[0,1]) == 0):
                 randlist2 = list(range(0,nprobes))  ###if the first nprobes are needed 
            if(int(np.vstack((np.char.split(name2,"_")))[0,1]) == 1): 
                for i in range(10000):
                    randnum = np.random.randint(0,len(fina))
                    if randnum not in randlist:
                        randlist.append(randnum)
                randlist2 = randlist[0:nprobes]
                randlist2 = natsorted(randlist2)
            final = [fina[i] for i in randlist2]
        else:
            final = fina
                
        gene = files[fil].split("/")[-1].split("_")[0]
        
        temp1 = readlist2[fil,0].split("{")[1]
        temp2 = temp1.split("}")[0]

        for t in range(0,len(final)):
            temp = str(Seq(final[t]).reverse_complement())
            rev.append([gene+"_RO"+temp2+"_"+(finprobes[t].split("#")[1]).split("_")[0], temp])
            revcom.append([gene+"_RO"+temp2+"_"+(finprobes[t].split("#")[1]).split("_")[0], primelist2[primernum,1]+readlist2[fil,1]+temp+readlist2[fil,1]+primelist2[primernum,2]])
        
        num = len(final)    
            
        # for m in range(0,len(final2)):
        #     temp2 = str(final2[m])
        #     rej.append([gene+"_"+str(m),temp2])
            
        
        lists.append([gene,str(num),lengen, bool(int(np.vstack((np.char.split(name2,"_")))[0,2])),len(name2),int(np.vstack((np.char.split(name2,"_")))[0,1])])
            
    revcoms = pd.DataFrame(revcom)
    revcoms.columns = ['ProbeName','Sequence']
    revs = pd.DataFrame(rev)
    revs.columns = ['ProbeName','Sequence']
    #rej = np.array(rej)
    lists = pd.DataFrame(lists)
    lists.columns = ['Gene','Final#probes','Length(bp)','Overlap','Initial#Probes',"Starting#probes_beforeoverlapfn"]
    
    # Apply the function row-wise
    lists.loc[lists['Starting#probes_beforeoverlapfn']==0,'Starting#probes_beforeoverlapfn'] = lists['Initial#Probes']

    revcoms.to_csv(outputfolder+"/MultiFISHProbes.csv")
    revs.to_csv(outputfolder+"/BindingProbes.csv")
    #pd.DataFrame(rej).to_csv(folder+"/Output/RejectedProbes.csv")
    lists.to_csv(outputfolder+"/NumberProbelist.csv")
    print('Completed Probe generation')

#print(type(revcom))

#print(np.array(revcom).shape)

#print(len(dup))
#list3 = pd.DataFrame(list2[:,9].T).T.drop_duplicates(keep=False).as_matrix()
#print(list3)
#print(len(list3[0]))