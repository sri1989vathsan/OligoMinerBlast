def bowtierun(bowtiedb,folder):
    import os
    import glob
    path = os.path.join(folder,"InitialProbes/*.fastq")
    files = glob.glob(path)
    files.sort()
    
    print(os.getcwd())

    #bowtiedb = "/Users/asrivath/Github/OligoMiner/GRCm38/GRCm38"
    #print(bowtiedb)

    for fil in files:
        out="./"+fil.split(".")[1]+".sam"
        logout = "./"+fil.split(".")[1]+".log"
        print("Aligning File "+fil)
        #print(out)
        os.system("bowtie2 -x" +bowtiedb+" -U"+ fil+" --no-hd -t -k 100 --local -D 20 -R 3 -N 1 -L 15 -i S,1,0.5 -S"+ out + " 2>"+logout)
        #os.system("bowtie -x" +bowtiedb+" -i"+ fil+" -v 5 -t -k 100 -S"+ out)
    print('Completed Bowtie alignment')
    
# if __name__ == '__main__':
#     # test1.py executed as script
#     # do something
#     bowtierun(bowtiedb)
    
#bowtie2 -x /path_to_hg38_index/hg38 -U 3.fastq --no-hd -t -k 100 --very-sensitive-local -S 3_u.sam
