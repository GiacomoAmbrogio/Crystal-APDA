###############################################################################
#                                                                             #
#                                                                             #
#                   ALL PROCESS DATA ANALYSIS  -  APDA                        #
#                                                                             #
#                                                                             #
#    Author: GIACOMO AMBROGIO - 2022                                          #
#    Usage:  Basic script for reading .meminfo files from MEMOPRT             #
#            'ALLPINT' and 'ALLPEXT' Crystal calculations                     #
#                                                                             #
#---------------------------------------------------                          #
#                 BEFORE USE!!!                                               #
#---------------------------------------------------                          #
#   Change in Variables section the Working Directories (wdir)                #
#   wdir is the main directory, in_dir and out_dir are the subdirectories     #
#   useed for INPUT and OUTPUT respectively                                   #
#                                                                             #
#                                                                             #
#   INSTRUCTIONS:                                                             #
#                                                                             #
#   Drop all *.meminfo files in the in_dir directory                          #
#                                                                             #
#   *.meminfo files from a single crystal calculation should have the same    #
#   name and extension ".pX.meminfo" where X is the number of the processors  #
#                                                                             #
#   This program also works fine with multiple files at the same time         #
#   (different names)                                                         #
#                                                                             #
#                                                                             #
#   There must always be a *.p0.meminfo file and *.p1.meminfo file            #
#   (at least 2 procs)                                                        #
#                                                                             #
###############################################################################



import os
import numpy as np
from matplotlib import pyplot as plt
from datetime import datetime
import re
import pandas as pd

#---------------------------------------------------
#                 VARIABLES
#---------------------------------------------------
#Colors for plot
iam0      = 'red'
altri     = 'k'

#Working directories
wdir      = "C:\\Folder\python\APDA"+"\\"
in_dir    = wdir + "DROP_IN"+"\\"
out_dir   = wdir + "OUTPUTS"+"\\"




#---------------------------------------------------
#                 FUNCTIONS
#---------------------------------------------------
#function to exract mem info from a .meminfo file

#return TIMVRS calls and DVMEM  
def DVMEM(infile):
    path    = in_dir+infile
    DVMEM   = []
    calls   = []
    file    = open(path, "r")
    for line in file:
        if "DVMEM" in line:
            element = line[8:19]
            calls.append(element.strip())
            element = line.split()
            ind     = element.index('DVMEM') + 1
            DVMEM.append(float(element[ind]))
    file.close()
    return calls, DVMEM

#return System VmPreak
def SYSMEM(infile):
    path = in_dir+infile
    file = open(path, "r")
    vm   = []
    for line in file:
        if "SYS VMEM(MB)" in line:
            element = line.split()
            ind     = element.index('VMEM(MB)') + 1
            vm.append(float(element[ind]))
    file.close()
    return vm

#return Node istant Memory
def NODE(infile):
    path = in_dir+infile
    nod  = []
    file = open(path, "r")
    for line in file:
        if "NOD TMEM(MB)" in line:
            element = line.split()
            ind     = element.index('TMEM(MB)') + 1
            nod.append(float(element[ind]))
    file.close()
    return nod

#return PID and chek if all pids are the same
def PID(infile):
    path = in_dir+infile
    pids = []
    file = open(path, "r")
    for line in file:
        if "PROCESS PID" in line:
            element = line.split()
            ind = element.index('PID') + 1
            pids.append(element[ind])
    for n in pids:
        if not pids[0]==n:
            print ('###############################')
            print ('WARNING: PID ERROR IN: '+infile)
    file.close()
    return pids[0]



#--------------------------------------------------
#                 START
#--------------------------------------------------
#Check directories
if not os.path.exists(wdir):
    os.mkdir(wdir)
if not os.path.exists(out_dir):
    os.mkdir(out_dir)
if not os.path.exists(in_dir):
    os.mkdir(in_dir)

#Time definition
now      = datetime.now()
day      = now.strftime("%b-%d-%Y")
#searching for input names
allfiles = os.listdir(in_dir)
r        = re.compile('.*\.p0\.')
p0files  = list(filter(r.match, allfiles))
names    = []
for file in p0files:
    a=file.replace('.p0.meminfo','')
    names.append(a)
print('Spotted files: ',names)

#---------------------------------------------------
#                 LOOP
#---------------------------------------------------
i=0
for name in names:
    i += 1
#counter for nproc
    nproc = sum(name+'.' in s for s in allfiles)
    proc = np.arange(1, nproc, 1, dtype=int)
    allproc = np.arange(0, nproc, 1, dtype=int)
    
#---------------------------------------------------
#                 ARRAY - DATAFRAME
#---------------------------------------------------
#opening file IAMEQ0
    file     = name+'.p0.meminfo'
    sub, mem = DVMEM(file)
    ARRAY    = np.array(sub)
    MEMORY   = np.array(mem)
    MEMORY   = np.vstack(MEMORY)

#loop for processors
    for p in proc:
        file     = name+'.p'+str(p)+'.meminfo'
        sub, mem = DVMEM(file)
        MEMORY2  = np.array(mem)
        MEMORY2  = np.vstack(MEMORY2)
        MEMORY   = np.append(MEMORY, MEMORY2, axis=1)
#identification of maximum memory/core used        
    maxmem = np.amax(MEMORY)
    DF     = pd.DataFrame(MEMORY)

 #---------------------------------------------------
 #                 PLOT
 #---------------------------------------------------
 #defining the 'call timvrs' axis
    x = np.array(range(len(ARRAY)))
 #ratio for figure sides
    ratio = len(ARRAY)*0.4
 #creating plot figure
    plt.figure   (num=i, figsize=(ratio, 4), edgecolor='k', dpi=300)  
    plt.plot     (x, DF[0], linewidth=1, marker=',', color=iam0, zorder=10)
    for p in proc:
        plt.plot (x, DF[p], linewidth=1, marker=',', color=altri, alpha=0.4)
#axis parameters
    plt.xlim   (left=1, right=len(ARRAY)-1)
    plt.ylim   (bottom=0)
    plt.ylabel ('Max MEM (MB)', fontsize=8)
    plt.yticks (fontsize=7)
    plt.xticks (x, sub, fontsize=4, rotation = 60)
    plt.vlines (x, 0, maxmem, alpha=0.2, color='k', linewidth=0.3)
    caption  = day+' - '+name
    plt.title  (caption, fontsize=8)
#saving figure in output directory
    plt.savefig(out_dir+name+'.png', dpi=300)

#---------------------------------------------------
#                 TABULAR.txt
#---------------------------------------------------
#creating tabular output file
    tabpath = out_dir+name+'.txt'
    if os.path.exists(tabpath):
        os.remove(tabpath)
    print('Creating output file in: ', tabpath)
    tab      = open(tabpath, 'a')
    tab.write('This file was generated with Python')
    tab.write('\nAuthor:        GIACOMO AMBROGIO')
    tab.write('\nCreation date: '+day)
    tab.write('\n\nThis file combines all the informations obtained from')
    tab.write('\n\"Chiara Ribaldone\'s subroutines\" in Crystal')
    tab.write('\n\n')
    tab.write('\n============================================================')
    tab.write('\nCrystal calculation:       '+name)
    tab.write('\n============================================================')
    tab.write('\n\n')
    tab.write('\n============================================================')
    tab.write('\nGeneral Informations:')
    tab.write('\n============================================================')
    tab.write('\n')
    tab.write('\nNumber of Processors:           '+str(nproc))
    tab.write('\nMaximum Memory per core:        '+str(maxmem)+' MB')
    tab.write('\n   Reached by processor:       ')
    for p in allproc:
        pp  = str(p).zfill(len(str(nproc-1)))
        hhh = DF[p].max()
        if hhh == maxmem:
            tab.write(' '+pp) 
    tab.write('\n')
    tab.write('\n============================================================')
    tab.write('\nProcessors Informations:')
    tab.write('\n============================================================')
    tab.write('\n\n')
    procindex = []   # Processors index
    indcalls  = []   # Number and name of 'max Mem' call to timvrs 
    pmaxmems  = []   # Max Memory usage
    BaseLine  = []   # VmPeak Baseline
    Pids      = []   # Pids
    for p in allproc:
        pp  = str(p).zfill(len(str(nproc-1)))
        procindex.append('p'+str(pp))
        indcall = DF[p].idxmax()
        pmaxmem = DF[p].max()
        howmanycalls = len(str(len(ARRAY)))
        CALLS = '['+str(indcall).zfill(howmanycalls)+'] '+str(ARRAY[indcall])
        indcalls.append(CALLS)
        pmaxmems.append(str(pmaxmem))
#for other info in .meminfo
        file   = name+'.p'+str(p)+'.meminfo'
        VmPeak = SYSMEM(file)
        BaseLine.append(VmPeak[0])
        Pid = PID(file)
        Pids.append(Pid)
#Dataframe creation: 
    dati = pd.DataFrame(list(zip(procindex,pmaxmems, indcalls,Pids, BaseLine)))
    Indici = ['Proc', 'MaxMEM', 'TIMVRS call', 'PID', 'BaseLine']
    tab.write(pd.DataFrame.to_string(dati, header=Indici, index=False, justify='justify'))
    tab.write('\n')
    DF.set_index(ARRAY, inplace=True, append=True)
    tab.write('\n============================================================')
    tab.write('\nAll Memory data acquired:')
    tab.write('\n============================================================')
    tab.write('\n\n')
    tab.write(pd.DataFrame.to_string(DF, header=procindex, index=True))
    tab.write('\n\n')
    tab.write('\n============================================================')
    tab.write('\nTermination')
    tab.write('\n============================================================')
    tab.close()
print('...')
print('...')
print('TERMINATION')

#---------------------------------------------------
#                 END
#---------------------------------------------------















    
