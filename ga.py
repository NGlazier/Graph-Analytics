import numpy as npa
import matplotlib.pyplot as plt
import seaborn as sns
from random import randrange
import csv
import plotly.express as px
import pandas as pd


#open file removing the np names, chromsome number, start and stop points
def read_data(fname):
    with open(fname) as textFile:
        data = [line.split() for line in textFile]
    data.pop(0)
    for x in range(len(data)):
        data[x].pop(0)
        data[x].pop(0)
        data[x].pop(0)
    return data

#function that calculates dmax using D, fA, and fB
def getdmax(det, a, b):
    if det < 0:
        return min(a*b, (1-a)*(1-b))
    elif det > 0:
        return min(b*(1-a), a*(1-b))
    else: 
        return 0

#open file containing hist1 region data
hist1 = read_data("Hist1.txt")
#print("Number of Genomic Windows in Hist1 Region:", len(hist1))
columns = [0] * (len(hist1[0]))


#remove unused columns and first 3 columns
for x in range(len(hist1)):
    for np in range(len(hist1[x])):
        if hist1[x][np] == '1':
            columns[np] = columns[np] + 1

difference = 0
#remove NP's not shown in Hist1 Region
for x in range(len(columns)):
    if columns[x] == 0:
        for y in range(len(hist1)):
            hist1[y].pop(x - difference)
        
        difference += 1
    
    
#calculate fA / fB for all windows
det_freq = [0] * (len(hist1))
norm_link = npa.zeros([len(hist1), len(hist1)])
for x in range(len(hist1)):
    for y in range(len(hist1[x])):
        if hist1[x][y] == '1':
            det_freq[x] = det_freq[x] + 1
for x in range(len(det_freq)):
    det_freq[x] = det_freq[x] / len(hist1[0])

#print(det_freq)

#calculate fAB
#calculate normalized linkage matrix and fill in the table for every combination of windows
for i in range(len(hist1)):
    for j in range(i, 81):
        fab = 0
        for k in range(len(hist1[i])):
            if hist1[i][k] == '1' and hist1[j][k] == '1':
                fab = fab + 1
        fab = fab / len(hist1[0])
        d = fab - (det_freq[i] * det_freq[j])
        dmax = getdmax(d, det_freq[i], det_freq[j])
        if dmax == 0:
            norm_link[i][j] = 0
            norm_link[j][i] = 0
        else:
            norm_link[i][j] = d / dmax
            norm_link[j][i] = d / dmax

#print the normalized linkage table
print(norm_link)
#for x in range(len(norm_link)):
#    for y in range(len(norm_link[x])):
#        print(round(norm_link[x][y], 2), end=" ", flush=True)
#    print()

#display normalized linkage table as a heatmap
plt.imshow(norm_link, cmap='hot', interpolation='nearest')
plt.title("Normalized Linkage Table")
plt.show()

