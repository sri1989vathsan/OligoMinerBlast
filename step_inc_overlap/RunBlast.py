def blastrun(blastrun, blastdb,blasttrnadb,blastrrnadb,folder):
    import os
    import glob
    path = os.path.join(folder,"2_InitialProbes/*.fa")
    files = glob.glob(path)
    files.sort()
    
    print(os.getcwd())

    #bowtiedb = "/Users/asrivath/Github/OligoMiner/GRCm38/GRCm38"
    #print(bowtiedb)

    for fil in files:
        out="."+fil.split(".")[1]+".txt"
        outtrna = "."+fil.split(".")[1]+"tRNA.txt"
        outrrna = "."+fil.split(".")[1]+"rRNA.txt"
        #print(fil)
        print("Blasting "+ fil)
        #os.system("bowtie2 -x" +bowtiedb+" -U"+ fil+" --no-hd -t -k 100 --very-sensitive-local -S"+ out)
        if (blastrun == 1):
            os.system("blastn -query "+ fil +" -db "+ blastdb +" -outfmt 6 -out "+ out +" -task blastn-short")
        os.system("blastn -query "+ fil +" -db "+ blasttrnadb +" -outfmt 6 -out "+ outtrna +" -task blastn-short")
        os.system("blastn -query "+ fil +" -db "+ blastrrnadb +" -outfmt 6 -out "+ outrrna +" -task blastn-short")
        
    print('Completed Blast alignment')
    
