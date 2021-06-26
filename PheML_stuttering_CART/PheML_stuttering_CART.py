import sys, getopt
import math
import numpy as np
import argparse
import os
from scipy import stats


parser = argparse.ArgumentParser()
parser.add_argument('-i','--phecodes', help='Input phecodes file', required=True)
parser.add_argument('-o','--output', help='Output file name of positive ID samples (list of predicted stuttering subjects)', required=True)

args = parser.parse_args()
print("Phecode file: ", args.phecodes)
print("Output file name: ", args.output)
input_phecode_file = str(args.phecodes)
output_file = str(args.output)


phe_enrich = open("./model/stuttering_extended_phewas_enrichment.csv","r")
phe_list = []
p_value_threshold = 0.000001
count_threshold = 19
phecode_index = {}
phe_ind_dic = {}
xin = 0
for line in phe_enrich:
    if(xin>1):
        spline = line.split(",")
        if(float(spline[3])<p_value_threshold and int(spline[2])>count_threshold):
            phe_list.append(spline[0])
            phe_ind_dic[xin-1] = spline[0]
    xin+=1
phe_list.sort()
ind_x = 0
for phe in phe_list:
    phecode_index[ind_x] = phe
    ind_x+=1


from sklearn import tree
from joblib import dump, load

model = load('./model/PheML_stutteringmodel.joblib')


BioVU_phe_file = open(input_phecode_file,"r")

#Get the phecodes for all the GRIDs so we can construct out firthtalbe
GRID_phe_dic = {}
x = 0
print("Getting phecodes for BioVU patients... this will take awhile")
BioVU_phe_file.readline()
for line in BioVU_phe_file:
    spline = (line.rstrip()).split(",")
    if(spline[0] not in GRID_phe_dic):
        GRID_phe_dic[spline[0]] = [spline[1]]
    else:
        GRID_phe_dic[spline[0]].append(spline[1])
    x+=1
####Making the table to feed into the prediction

total_BV_table = []
GRID_index_reference_list = []
for GRID in GRID_phe_dic:
    GRID_index_reference_list.append(GRID)
    ind_list_tab = []
    for inx in phe_ind_dic:
        if(phe_ind_dic[inx] in GRID_phe_dic[GRID]):
            ind_list_tab.append(1)
        else:
            ind_list_tab.append(0)
    total_BV_table.append(ind_list_tab)
BioVU_X = np.array(total_BV_table)
BioVU_y_predict = model.predict(BioVU_X)

pred_ind = 0
positive_hits = []
for person in BioVU_y_predict:
    if(person==1):
        positive_hits.append(GRID_index_reference_list[pred_ind])
    pred_ind+=1



print("Total identified patients: ",len(positive_hits))

f = open("%s.txt" % (output_file),"w+")

for goober in positive_hits:
    f.write(goober)
    f.write("\n")


