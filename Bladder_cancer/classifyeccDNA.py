import glob
import re
## first detect circular numbers
def getOrignalAMP(inputPath):
    shuju={}
    st=""
    for file in sorted(glob.glob(f"{inputPath}/*_annotated_cycles.txt")):
        tmp=re.findall(r"\d+_amplicon\d+",file)[0]
        sample,amp=tmp.split("_")
        if st!=sample:
            shuju[sample]=[]
            shuju[sample].append(amp)
            st=sample
        else:
            shuju[sample].append(amp)
    return shuju

### second define class

weight={"BFB":1,"ecDNA":2,"Complex cyclic":3,"Complex non-cyclic":4,"Linear amplification":5,"unknown":6}

def defineCycle(inputPath,shuju):
    resut=[]
    fin=[]
    for file in sorted(glob.glob(f"{inputPath}/*_intervals.bed")):
        sample,amp,type,*ot=file.replace(f"{inputPath}/","").split("_")
        resut.append({"sample":sample,"amp":amp,"type":type})
    for (k,v) in shuju.items():
        for woco in v:
            rst=[x["type"] for x in resut if x["sample"]==k and x["amp"]==woco]
            if len(rst)>0:
                heihei= sorted(rst,key=lambda x:weight[x])[0]
                if heihei=="unknown":
                    fin.append([k,woco,"Invalid"])
                else:
                    fin.append([k,woco,heihei])
            else:
                fin.append([k,woco,"Invalid"]) 
    return fin
if __name__ == "__main__":
    fin=defineCycle("all_detective_classification_bed_files",getOrignalAMP("all_detective_annotated_cycles_files"))
    for wori in fin:
        print(wori[0],"\t",wori[1],"\t",wori[2],sep="")





