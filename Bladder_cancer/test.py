import sys
import re
import operator
from itertools import chain


#Declare lists
results_TRANS_FINAL = []
insertion_candidates = []
TRANSLOCATION_same = []

#Declare arguments
input_file_SVs = "test9.3.txt"
patient_id = "1"


########## RETRIEVE TRANSLOCATIONS ##########

#Open and read the  variant callers results files
file = open(input_file_SVs,'r')
SVs_results = re.split('\n',file.read())
for res in range(1,len(SVs_results)-1):
    res_line = re.split('\t',SVs_results[res])
    results_TRANS_FINAL.append(res_line)

#Change chrX and chrY for 23 and 24 to order them as integers
for T in range(0,len(results_TRANS_FINAL)-1):
    if results_TRANS_FINAL[T][1]=='X' and results_TRANS_FINAL[T][3]=='Y':
        results_TRANS_FINAL[T][1] = '23'
        results_TRANS_FINAL[T][3] = '24'
    elif results_TRANS_FINAL[T][1]=='Y' and results_TRANS_FINAL[T][3]=='X':
        results_TRANS_FINAL[T][1] = '24'
        results_TRANS_FINAL[T][3] = '23'
    elif results_TRANS_FINAL[T][1]=='X':
        results_TRANS_FINAL[T][1] = '23'
    elif results_TRANS_FINAL[T][1]=='Y':
        results_TRANS_FINAL[T][1] = '24'
    elif results_TRANS_FINAL[T][3]=='X':
        results_TRANS_FINAL[T][3] = '23'
    elif results_TRANS_FINAL[T][3]=='Y':
        results_TRANS_FINAL[T][3] = '24'
    else:
        pass


########## FIRST ROUND TO RETRIEVE INTEGRATIONS - ORDER CHR1 CHR2 POS1 POS2 ##########

#Re-sort list
results_TRANS_FINAL = sorted(results_TRANS_FINAL, key=operator.itemgetter(1,3,2,4),  reverse=False)

#Get translocations that can be putative insertions
for b in range(len(results_TRANS_FINAL)):
    sv = b
    num = sv + 1
    #Save translocations of same insertion in TRANSLOCATION_same
    if num < len(results_TRANS_FINAL) and sv < len(results_TRANS_FINAL):
        #same chromosomes involved - inserted region bigger than acceptor region - circle-genome x2 - inserted region < 2411260bp (bigger circle)
        if results_TRANS_FINAL[sv][3] == results_TRANS_FINAL[num][3] and results_TRANS_FINAL[sv][1] == results_TRANS_FINAL[num][1] and abs(int(results_TRANS_FINAL[sv][2])-int(results_TRANS_FINAL[num][2]))>=abs(int(results_TRANS_FINAL[sv][4])-int(results_TRANS_FINAL[num][4])) and results_TRANS_FINAL[sv][7]=="circle-genome" and results_TRANS_FINAL[num][7]=="circle-genome" and abs(int(results_TRANS_FINAL[sv][2])-int(results_TRANS_FINAL[num][2]))<2411260:
            if num < len(results_TRANS_FINAL):
                insertion_candidate1 = list(results_TRANS_FINAL[sv])
                insertion_candidate2 = list(results_TRANS_FINAL[num])
                insertion_candidate1.append("insertion_cluster_" + str(results_TRANS_FINAL[sv][1]) + "_" + str(results_TRANS_FINAL[sv][2]) + "_" + str(results_TRANS_FINAL[sv][3]) + "_" + str(results_TRANS_FINAL[sv][4]) + ".." + str(results_TRANS_FINAL[num][1]) + "_" + str(results_TRANS_FINAL[num][2]) + "_" + str(results_TRANS_FINAL[num][3]) + "_" + str(results_TRANS_FINAL[num][4]))
                insertion_candidate2.append("insertion_cluster_" + str(results_TRANS_FINAL[sv][1]) + "_" + str(results_TRANS_FINAL[sv][2]) + "_" + str(results_TRANS_FINAL[sv][3]) + "_" + str(results_TRANS_FINAL[sv][4]) + ".." + str(results_TRANS_FINAL[num][1]) + "_" + str(results_TRANS_FINAL[num][2]) + "_" + str(results_TRANS_FINAL[num][3]) + "_" + str(results_TRANS_FINAL[num][4]))
                insertion_candidate1.append("insertion_" + str(results_TRANS_FINAL[sv][1]) + "_" + str(results_TRANS_FINAL[sv][2]) + "_" + str(results_TRANS_FINAL[sv][3]) + "_" + str(results_TRANS_FINAL[sv][4]))
                insertion_candidate2.append("insertion_" + str(results_TRANS_FINAL[num][1]) + "_" + str(results_TRANS_FINAL[num][2]) + "_" + str(results_TRANS_FINAL[num][3]) + "_" + str(results_TRANS_FINAL[num][4]))    
                TRANSLOCATION_same.append(insertion_candidate1)
                TRANSLOCATION_same.append(insertion_candidate2)
                num = num+1
            if num == len(results_TRANS_FINAL):
                break
        #same chromosomes involved - inserted region bigger than acceptor region - circle-genome + genome-genome - inserted region < 2411260bp (bigger circle)
        elif results_TRANS_FINAL[sv][3] == results_TRANS_FINAL[num][3] and results_TRANS_FINAL[sv][1] == results_TRANS_FINAL[num][1] and abs(int(results_TRANS_FINAL[sv][2])-int(results_TRANS_FINAL[num][2]))>=abs(int(results_TRANS_FINAL[sv][4])-int(results_TRANS_FINAL[num][4])) and results_TRANS_FINAL[sv][7]=="circle-genome" and results_TRANS_FINAL[num][7]=="genome-genome" and abs(int(results_TRANS_FINAL[sv][2])-int(results_TRANS_FINAL[num][2]))<2411260:
            if num < len(results_TRANS_FINAL):
                insertion_candidate1 = list(results_TRANS_FINAL[sv])
                insertion_candidate2 = list(results_TRANS_FINAL[num])
                insertion_candidate1.append("insertion_cluster_" + str(results_TRANS_FINAL[sv][1]) + "_" + str(results_TRANS_FINAL[sv][2]) + "_" + str(results_TRANS_FINAL[sv][3]) + "_" + str(results_TRANS_FINAL[sv][4]) + ".." + str(results_TRANS_FINAL[num][1]) + "_" + str(results_TRANS_FINAL[num][2]) + "_" + str(results_TRANS_FINAL[num][3]) + "_" + str(results_TRANS_FINAL[num][4]))
                insertion_candidate2.append("insertion_cluster_" + str(results_TRANS_FINAL[sv][1]) + "_" + str(results_TRANS_FINAL[sv][2]) + "_" + str(results_TRANS_FINAL[sv][3]) + "_" + str(results_TRANS_FINAL[sv][4]) + ".." + str(results_TRANS_FINAL[num][1]) + "_" + str(results_TRANS_FINAL[num][2]) + "_" + str(results_TRANS_FINAL[num][3]) + "_" + str(results_TRANS_FINAL[num][4]))
                insertion_candidate1.append("insertion_" + str(results_TRANS_FINAL[sv][1]) + "_" + str(results_TRANS_FINAL[sv][2]) + "_" + str(results_TRANS_FINAL[sv][3]) + "_" + str(results_TRANS_FINAL[sv][4]))
                insertion_candidate2.append("insertion_" + str(results_TRANS_FINAL[num][1]) + "_" + str(results_TRANS_FINAL[num][2]) + "_" + str(results_TRANS_FINAL[num][3]) + "_" + str(results_TRANS_FINAL[num][4]))    
                TRANSLOCATION_same.append(insertion_candidate1)
                TRANSLOCATION_same.append(insertion_candidate2)
                num = num+1
            if num == len(results_TRANS_FINAL):
                break
        #same chromosomes involved - inserted region bigger than acceptor region - genome-genome + circle-genome - inserted region < 2411260bp (bigger circle)
        elif results_TRANS_FINAL[sv][3] == results_TRANS_FINAL[num][3] and results_TRANS_FINAL[sv][1] == results_TRANS_FINAL[num][1] and abs(int(results_TRANS_FINAL[sv][2])-int(results_TRANS_FINAL[num][2]))>=abs(int(results_TRANS_FINAL[sv][4])-int(results_TRANS_FINAL[num][4])) and results_TRANS_FINAL[sv][7]=="genome-genome" and results_TRANS_FINAL[num][7]=="circle-genome" and abs(int(results_TRANS_FINAL[sv][2])-int(results_TRANS_FINAL[num][2]))<2411260:
            if num < len(results_TRANS_FINAL):
                insertion_candidate1 = list(results_TRANS_FINAL[sv])
                insertion_candidate2 = list(results_TRANS_FINAL[num])
                insertion_candidate1.append("insertion_cluster_" + str(results_TRANS_FINAL[sv][1]) + "_" + str(results_TRANS_FINAL[sv][2]) + "_" + str(results_TRANS_FINAL[sv][3]) + "_" + str(results_TRANS_FINAL[sv][4]) + ".." + str(results_TRANS_FINAL[num][1]) + "_" + str(results_TRANS_FINAL[num][2]) + "_" + str(results_TRANS_FINAL[num][3]) + "_" + str(results_TRANS_FINAL[num][4]))
                insertion_candidate2.append("insertion_cluster_" + str(results_TRANS_FINAL[sv][1]) + "_" + str(results_TRANS_FINAL[sv][2]) + "_" + str(results_TRANS_FINAL[sv][3]) + "_" + str(results_TRANS_FINAL[sv][4]) + ".." + str(results_TRANS_FINAL[num][1]) + "_" + str(results_TRANS_FINAL[num][2]) + "_" + str(results_TRANS_FINAL[num][3]) + "_" + str(results_TRANS_FINAL[num][4]))
                insertion_candidate1.append("insertion_" + str(results_TRANS_FINAL[sv][1]) + "_" + str(results_TRANS_FINAL[sv][2]) + "_" + str(results_TRANS_FINAL[sv][3]) + "_" + str(results_TRANS_FINAL[sv][4]))
                insertion_candidate2.append("insertion_" + str(results_TRANS_FINAL[num][1]) + "_" + str(results_TRANS_FINAL[num][2]) + "_" + str(results_TRANS_FINAL[num][3]) + "_" + str(results_TRANS_FINAL[num][4]))    
                TRANSLOCATION_same.append(insertion_candidate1)
                TRANSLOCATION_same.append(insertion_candidate2)
                num = num+1
            if num == len(results_TRANS_FINAL):
                break
        #same chromosomes involved - inserted region bigger than acceptor region - genome-genome + genome-genome - inserted region < 2411260bp (bigger circle)
        elif results_TRANS_FINAL[sv][3] == results_TRANS_FINAL[num][3] and results_TRANS_FINAL[sv][1] == results_TRANS_FINAL[num][1] and (abs(int(results_TRANS_FINAL[sv][2])-int(results_TRANS_FINAL[num][2]))>=abs(int(results_TRANS_FINAL[sv][4])-int(results_TRANS_FINAL[num][4])) or abs(int(results_TRANS_FINAL[sv][2])-int(results_TRANS_FINAL[num][2]))<=abs(int(results_TRANS_FINAL[sv][4])-int(results_TRANS_FINAL[num][4]))) and results_TRANS_FINAL[sv][7]=="genome-genome" and results_TRANS_FINAL[num][7]=="genome-genome" and (abs(int(results_TRANS_FINAL[sv][2])-int(results_TRANS_FINAL[num][2]))<2411260 and abs(int(results_TRANS_FINAL[sv][4])-int(results_TRANS_FINAL[num][4]))<2411260):
            if num < len(results_TRANS_FINAL):
                insertion_candidate1 = list(results_TRANS_FINAL[sv])
                insertion_candidate2 = list(results_TRANS_FINAL[num])
                insertion_candidate1.append("insertion_cluster_" + str(results_TRANS_FINAL[sv][1]) + "_" + str(results_TRANS_FINAL[sv][2]) + "_" + str(results_TRANS_FINAL[sv][3]) + "_" + str(results_TRANS_FINAL[sv][4]) + ".." + str(results_TRANS_FINAL[num][1]) + "_" + str(results_TRANS_FINAL[num][2]) + "_" + str(results_TRANS_FINAL[num][3]) + "_" + str(results_TRANS_FINAL[num][4]))
                insertion_candidate2.append("insertion_cluster_" + str(results_TRANS_FINAL[sv][1]) + "_" + str(results_TRANS_FINAL[sv][2]) + "_" + str(results_TRANS_FINAL[sv][3]) + "_" + str(results_TRANS_FINAL[sv][4]) + ".." + str(results_TRANS_FINAL[num][1]) + "_" + str(results_TRANS_FINAL[num][2]) + "_" + str(results_TRANS_FINAL[num][3]) + "_" + str(results_TRANS_FINAL[num][4]))
                insertion_candidate1.append("insertion_" + str(results_TRANS_FINAL[sv][1]) + "_" + str(results_TRANS_FINAL[sv][2]) + "_" + str(results_TRANS_FINAL[sv][3]) + "_" + str(results_TRANS_FINAL[sv][4]))
                insertion_candidate2.append("insertion_" + str(results_TRANS_FINAL[num][1]) + "_" + str(results_TRANS_FINAL[num][2]) + "_" + str(results_TRANS_FINAL[num][3]) + "_" + str(results_TRANS_FINAL[num][4]))    
                TRANSLOCATION_same.append(insertion_candidate1)
                TRANSLOCATION_same.append(insertion_candidate2)
                num = num+1
            if num == len(results_TRANS_FINAL):
                break    
        else:
            pass
    else:
        break


########## SECOND ROUND TO RETRIEVE INTEGRATIONS - ORDER CHR1 CHR2 POS2 POS1 ##########

#Same as previously but with re-sorted list to get all the possible integrations
#Re-sort list
results_TRANS_FINAL = sorted(results_TRANS_FINAL, key=operator.itemgetter(1,3,4,2),  reverse=False)

for b in range(len(results_TRANS_FINAL)):
    sv = b
    num = sv + 1
    if num < len(results_TRANS_FINAL) and sv < len(results_TRANS_FINAL):
        if results_TRANS_FINAL[sv][3] == results_TRANS_FINAL[num][3] and results_TRANS_FINAL[sv][1] == results_TRANS_FINAL[num][1] and abs(int(results_TRANS_FINAL[sv][2])-int(results_TRANS_FINAL[num][2]))>=abs(int(results_TRANS_FINAL[sv][4])-int(results_TRANS_FINAL[num][4])) and results_TRANS_FINAL[sv][7]=="circle-genome" and results_TRANS_FINAL[num][7]=="circle-genome" and abs(int(results_TRANS_FINAL[sv][2])-int(results_TRANS_FINAL[num][2]))<2411260:
            if num < len(results_TRANS_FINAL):
                insertion_candidate1 = list(results_TRANS_FINAL[sv])
                insertion_candidate2 = list(results_TRANS_FINAL[num])
                insertion_candidate1.append("insertion_cluster_" + str(results_TRANS_FINAL[sv][1]) + "_" + str(results_TRANS_FINAL[sv][2]) + "_" + str(results_TRANS_FINAL[sv][3]) + "_" + str(results_TRANS_FINAL[sv][4]) + ".." + str(results_TRANS_FINAL[num][1]) + "_" + str(results_TRANS_FINAL[num][2]) + "_" + str(results_TRANS_FINAL[num][3]) + "_" + str(results_TRANS_FINAL[num][4]))
                insertion_candidate2.append("insertion_cluster_" + str(results_TRANS_FINAL[sv][1]) + "_" + str(results_TRANS_FINAL[sv][2]) + "_" + str(results_TRANS_FINAL[sv][3]) + "_" + str(results_TRANS_FINAL[sv][4]) + ".." + str(results_TRANS_FINAL[num][1]) + "_" + str(results_TRANS_FINAL[num][2]) + "_" + str(results_TRANS_FINAL[num][3]) + "_" + str(results_TRANS_FINAL[num][4]))
                insertion_candidate1.append("insertion_" + str(results_TRANS_FINAL[sv][1]) + "_" + str(results_TRANS_FINAL[sv][2]) + "_" + str(results_TRANS_FINAL[sv][3]) + "_" + str(results_TRANS_FINAL[sv][4]))
                insertion_candidate2.append("insertion_" + str(results_TRANS_FINAL[num][1]) + "_" + str(results_TRANS_FINAL[num][2]) + "_" + str(results_TRANS_FINAL[num][3]) + "_" + str(results_TRANS_FINAL[num][4]))    
                TRANSLOCATION_same.append(insertion_candidate1)
                TRANSLOCATION_same.append(insertion_candidate2)
                num = num+1
            if num == len(results_TRANS_FINAL):
                break
        elif results_TRANS_FINAL[sv][3] == results_TRANS_FINAL[num][3] and results_TRANS_FINAL[sv][1] == results_TRANS_FINAL[num][1] and abs(int(results_TRANS_FINAL[sv][2])-int(results_TRANS_FINAL[num][2]))>=abs(int(results_TRANS_FINAL[sv][4])-int(results_TRANS_FINAL[num][4])) and results_TRANS_FINAL[sv][7]=="circle-genome" and results_TRANS_FINAL[num][7]=="genome-genome" and abs(int(results_TRANS_FINAL[sv][2])-int(results_TRANS_FINAL[num][2]))<2411260:
            if num < len(results_TRANS_FINAL):
                insertion_candidate1 = list(results_TRANS_FINAL[sv])
                insertion_candidate2 = list(results_TRANS_FINAL[num])
                insertion_candidate1.append("insertion_cluster_" + str(results_TRANS_FINAL[sv][1]) + "_" + str(results_TRANS_FINAL[sv][2]) + "_" + str(results_TRANS_FINAL[sv][3]) + "_" + str(results_TRANS_FINAL[sv][4]) + ".." + str(results_TRANS_FINAL[num][1]) + "_" + str(results_TRANS_FINAL[num][2]) + "_" + str(results_TRANS_FINAL[num][3]) + "_" + str(results_TRANS_FINAL[num][4]))
                insertion_candidate2.append("insertion_cluster_" + str(results_TRANS_FINAL[sv][1]) + "_" + str(results_TRANS_FINAL[sv][2]) + "_" + str(results_TRANS_FINAL[sv][3]) + "_" + str(results_TRANS_FINAL[sv][4]) + ".." + str(results_TRANS_FINAL[num][1]) + "_" + str(results_TRANS_FINAL[num][2]) + "_" + str(results_TRANS_FINAL[num][3]) + "_" + str(results_TRANS_FINAL[num][4]))
                insertion_candidate1.append("insertion_" + str(results_TRANS_FINAL[sv][1]) + "_" + str(results_TRANS_FINAL[sv][2]) + "_" + str(results_TRANS_FINAL[sv][3]) + "_" + str(results_TRANS_FINAL[sv][4]))
                insertion_candidate2.append("insertion_" + str(results_TRANS_FINAL[num][1]) + "_" + str(results_TRANS_FINAL[num][2]) + "_" + str(results_TRANS_FINAL[num][3]) + "_" + str(results_TRANS_FINAL[num][4]))    
                TRANSLOCATION_same.append(insertion_candidate1)
                TRANSLOCATION_same.append(insertion_candidate2)
                num = num+1
            if num == len(results_TRANS_FINAL):
                break
        elif results_TRANS_FINAL[sv][3] == results_TRANS_FINAL[num][3] and results_TRANS_FINAL[sv][1] == results_TRANS_FINAL[num][1] and abs(int(results_TRANS_FINAL[sv][2])-int(results_TRANS_FINAL[num][2]))>=abs(int(results_TRANS_FINAL[sv][4])-int(results_TRANS_FINAL[num][4])) and results_TRANS_FINAL[sv][7]=="genome-genome" and results_TRANS_FINAL[num][7]=="circle-genome" and abs(int(results_TRANS_FINAL[sv][2])-int(results_TRANS_FINAL[num][2]))<2411260:
            if num < len(results_TRANS_FINAL):
                insertion_candidate1 = list(results_TRANS_FINAL[sv])
                insertion_candidate2 = list(results_TRANS_FINAL[num])
                insertion_candidate1.append("insertion_cluster_" + str(results_TRANS_FINAL[sv][1]) + "_" + str(results_TRANS_FINAL[sv][2]) + "_" + str(results_TRANS_FINAL[sv][3]) + "_" + str(results_TRANS_FINAL[sv][4]) + ".." + str(results_TRANS_FINAL[num][1]) + "_" + str(results_TRANS_FINAL[num][2]) + "_" + str(results_TRANS_FINAL[num][3]) + "_" + str(results_TRANS_FINAL[num][4]))
                insertion_candidate2.append("insertion_cluster_" + str(results_TRANS_FINAL[sv][1]) + "_" + str(results_TRANS_FINAL[sv][2]) + "_" + str(results_TRANS_FINAL[sv][3]) + "_" + str(results_TRANS_FINAL[sv][4]) + ".." + str(results_TRANS_FINAL[num][1]) + "_" + str(results_TRANS_FINAL[num][2]) + "_" + str(results_TRANS_FINAL[num][3]) + "_" + str(results_TRANS_FINAL[num][4]))
                insertion_candidate1.append("insertion_" + str(results_TRANS_FINAL[sv][1]) + "_" + str(results_TRANS_FINAL[sv][2]) + "_" + str(results_TRANS_FINAL[sv][3]) + "_" + str(results_TRANS_FINAL[sv][4]))
                insertion_candidate2.append("insertion_" + str(results_TRANS_FINAL[num][1]) + "_" + str(results_TRANS_FINAL[num][2]) + "_" + str(results_TRANS_FINAL[num][3]) + "_" + str(results_TRANS_FINAL[num][4]))    
                TRANSLOCATION_same.append(insertion_candidate1)
                TRANSLOCATION_same.append(insertion_candidate2)
                num = num+1
            if num == len(results_TRANS_FINAL):
                break
        elif results_TRANS_FINAL[sv][3] == results_TRANS_FINAL[num][3] and results_TRANS_FINAL[sv][1] == results_TRANS_FINAL[num][1] and (abs(int(results_TRANS_FINAL[sv][2])-int(results_TRANS_FINAL[num][2]))>=abs(int(results_TRANS_FINAL[sv][4])-int(results_TRANS_FINAL[num][4])) or abs(int(results_TRANS_FINAL[sv][2])-int(results_TRANS_FINAL[num][2]))<=abs(int(results_TRANS_FINAL[sv][4])-int(results_TRANS_FINAL[num][4]))) and results_TRANS_FINAL[sv][7]=="genome-genome" and results_TRANS_FINAL[num][7]=="genome-genome" and (abs(int(results_TRANS_FINAL[sv][2])-int(results_TRANS_FINAL[num][2]))<2411260 and abs(int(results_TRANS_FINAL[sv][4])-int(results_TRANS_FINAL[num][4]))<2411260):
            if num < len(results_TRANS_FINAL):
                insertion_candidate1 = list(results_TRANS_FINAL[sv])
                insertion_candidate2 = list(results_TRANS_FINAL[num])
                insertion_candidate1.append("insertion_cluster_" + str(results_TRANS_FINAL[sv][1]) + "_" + str(results_TRANS_FINAL[sv][2]) + "_" + str(results_TRANS_FINAL[sv][3]) + "_" + str(results_TRANS_FINAL[sv][4]) + ".." + str(results_TRANS_FINAL[num][1]) + "_" + str(results_TRANS_FINAL[num][2]) + "_" + str(results_TRANS_FINAL[num][3]) + "_" + str(results_TRANS_FINAL[num][4]))
                insertion_candidate2.append("insertion_cluster_" + str(results_TRANS_FINAL[sv][1]) + "_" + str(results_TRANS_FINAL[sv][2]) + "_" + str(results_TRANS_FINAL[sv][3]) + "_" + str(results_TRANS_FINAL[sv][4]) + ".." + str(results_TRANS_FINAL[num][1]) + "_" + str(results_TRANS_FINAL[num][2]) + "_" + str(results_TRANS_FINAL[num][3]) + "_" + str(results_TRANS_FINAL[num][4]))
                insertion_candidate1.append("insertion_" + str(results_TRANS_FINAL[sv][1]) + "_" + str(results_TRANS_FINAL[sv][2]) + "_" + str(results_TRANS_FINAL[sv][3]) + "_" + str(results_TRANS_FINAL[sv][4]))
                insertion_candidate2.append("insertion_" + str(results_TRANS_FINAL[num][1]) + "_" + str(results_TRANS_FINAL[num][2]) + "_" + str(results_TRANS_FINAL[num][3]) + "_" + str(results_TRANS_FINAL[num][4]))
                TRANSLOCATION_same.append(insertion_candidate1)
                TRANSLOCATION_same.append(insertion_candidate2)
                num = num+1
            if num == len(results_TRANS_FINAL):
                break    
        else:
            pass
    else:
        break


#Remove VCF duplicates from lists (function from https://stackoverflow.com/questions/34630772/removing-duplicate-elements-by-their-attributes-in-python)
#Duplicates are selected by the 5th element, the info which includes the SV id

def deduplicate(items):
    seen = set()
    for item in items:
        if not item[9] in seen:
            seen.add(item[9])
            yield item


insertion_candidates.append(list(deduplicate(TRANSLOCATION_same)))
results_TRANS_FINAL = list(chain.from_iterable(insertion_candidates))


########## SAVE RESULTS IN FILE ##########

#Open the file to write in
results_patient = open('Insertion_candidates_from_SVs_circleseqcircles_patient_' + patient_id + '.txt','w')

#Change chr 23 and 24 for X and Y again
for T in range(0,len(results_TRANS_FINAL)):
    if results_TRANS_FINAL[T][1]=='23' and results_TRANS_FINAL[T][3]=='24':
        results_TRANS_FINAL[T][1] = 'X'
        results_TRANS_FINAL[T][3] = 'Y'
    elif results_TRANS_FINAL[T][1]=='24' and results_TRANS_FINAL[T][3]=='23':
        results_TRANS_FINAL[T][1] = 'Y'
        results_TRANS_FINAL[T][3] = 'X'
    elif results_TRANS_FINAL[T][1]=='23':
        results_TRANS_FINAL[T][1] = 'X'
    elif results_TRANS_FINAL[T][1]=='24':
        results_TRANS_FINAL[T][1] = 'Y'
    elif results_TRANS_FINAL[T][3]=='23':
        results_TRANS_FINAL[T][3] = 'X'
    elif results_TRANS_FINAL[T][3]=='24':
        results_TRANS_FINAL[T][3] = 'Y'
    else:
        pass

#Print results

for item in results_TRANS_FINAL:
    print(item)