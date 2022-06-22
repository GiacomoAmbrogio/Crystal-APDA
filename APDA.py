###############################################################################
#                                                                             #
#                                                                             #
#                   ALL PROCESS DATA ANALYSIS  -  APDA                        #
#                                                                             #
#                                                                             #
#    Author: GIACOMO AMBROGIO - 2022                                          #
#            University of Turin                                              #
#                                                                             #
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
from tabulate import tabulate
import configparser as cnfg
import matplotlib.font_manager as font_manager



#---------------------------------------------------
#                 CONFIGURATION
#---------------------------------------------------
path_py = os.path.abspath(__file__)
wdir    = os.path.dirname(path_py)+"\\"
#reading config.ini file
cp = cnfg.ConfigParser()
#check config.ini file
if not os.path.exists(wdir+'config.ini'):
    config = open(wdir+'config.ini', 'a')
    config.write('[APDA]\n')
    config.write('DEBUG = False\n')
    config.write('PLOT= True\n')
    config.write('TABULAR = True\n')
    config.write('INTERNAL = True\n')
    config.write('EXTERNAL = True\n')
    config.write('#    Only accept True or False\n\n\n')
    config.write('[DIRECTORIES]\n')
    config.write('INPUT = DROP_IN\n')
    config.write('OUTPUT = OUTPUTS\n\n\n')
    config.write('[PLOT]\n')
    config.write('SEPARATEPLOTS = True\n')
    config.write('#    Only accept True or False, ')
    config.write('create two different plot figures\n\n')
    config.write('IAMEQ0_COLOR = default\n')
    config.write('OTHER_PROCS_COLOR = default\n')
    config.write('#    Only accept two Hex colors ')
    config.write('(example IAMEQ0_COLOR = #00FF00 #FF0000)\n')
    config.write('#    First one is for internal, second is for external\n\n\n')
    config.close()
cp.read(wdir+'config.ini')
#setting variables
debug    = cp['APDA'].getboolean('DEBUG')
makeplot = cp['APDA'].getboolean('PLOT')
maketab  = cp['APDA'].getboolean('TABULAR')
intoverride = cp['APDA'].getboolean('INTERNAL')
extoverride = cp['APDA'].getboolean('EXTERNAL')
in_dir   = wdir+cp['DIRECTORIES']['INPUT']+"\\"
out_dir  = wdir+cp['DIRECTORIES']['OUTPUT']+"\\"
sepplots = cp['PLOT'].getboolean('SEPARATEPLOTS')
#plot lines colors settings
if cp['PLOT']['IAMEQ0_COLOR'] == 'default':
    iam0=['green', 'red']
else:
    iam0=cp['PLOT']['IAMEQ0_COLOR'].split()
if cp['PLOT']['OTHER_PROCS_COLOR'] == 'default':
    altri=['#004000', '#400000']
else:
    altri=cp['PLOT']['OTHER_PROCS_COLOR'].split()
#font for labels in plots
font = font_manager.FontProperties(family='Consolas',style='normal', size=7)
#Console output
print('--------------------------')
print('SETTINGS:')
print('--------------------------')
print('Debug        : ', debug)
print('Make PLOT    : ', makeplot)
print('Make TABULAR : ', maketab)
print('Internal mem : ', intoverride)
print('External mem : ', extoverride)
print('--------------------------')
print()
if debug: print('config.ini   : ', wdir+'config.ini')
if debug: print('Working dir  : ', wdir)
if debug: print('Input dir    : ', in_dir)
if debug: print('Output dir   : ', out_dir)
if debug: print('Make 2 plots : ', sepplots)
if debug: print('IAMEQ0 colors: ', iam0)
if debug: print('OTHER colors : ', altri)

#---------------------------------------------------
#                 FUNCTIONS
#---------------------------------------------------
#function to exract mem info from a .meminfo file

#define if ALLPINT and/or ALLPEXT was used
def MEMOPRT(infile):
    path = in_dir+infile
    file    = open(path, "r")
    internal = False
    external = False
    for line in file:
        if "MM_INT" in line: internal = True
        if "MM_EXT" in line: external = True
    if not (internal or external):
        print ('##################################################')
        print ('WARNING: ERROR SEARCHING LINES IN: '+infile)
    file.close()
    return internal, external

#return TIMVRS calls and DVMEM  
def DVMEM(infile):
    path    = in_dir+infile
    dvmem   = []
    calls   = []
    file    = open(path, "r")
    for line in file:
        if "DVMEM" in line:
            element = line[8:19]
            calls.append(element.strip())
            element = line.split()
            ind     = element.index('DVMEM') + 1
            dvmem.append(float(element[ind]))
    file.close()
    return calls, dvmem

#return TIMVRS calls and MAX MEM  
def USED(infile):
    path    = in_dir+infile
    used    = []
    calls   = []
    file    = open(path, "r")
    for line in file:
        if "MAX MEM" in line:
            element = line[8:19]
            calls.append(element.strip())
            element = line.split()
            ind     = element.index('MEM') + 1
            used.append(float(element[ind]))
    file.close()
    return calls, used

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

#return Node istant Memory  - NOT IMPLEMENTED
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
            print ('#######################################')
            print ('WARNING: PID ERROR IN: '+infile)
    file.close()
    return pids[0]





#--------------------------------------------------
#                 START
#--------------------------------------------------
#Check directories
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
if debug: print('Spotted files: ',names)

#---------------------------------------------------
#                 LOOP
#---------------------------------------------------
if debug: print ()
if debug: print ('LOOP start')
i=0
for name in names:
    i += 1
    if debug: print ('--------------------------')
    if debug: print ('File         : ', name)
#counter for nproc
    nproc = sum(name+'.' in s for s in allfiles)
    if debug: print ('nProc        : ', nproc)
    proc    = np.arange(1, nproc, 1, dtype=int)
    allproc = np.arange(0, nproc, 1, dtype=int)
#---------------------------------------------------
#                 ARRAY - DATAFRAME
#---------------------------------------------------
#opening file IAMEQ0
    file = name+'.p0.meminfo'
#looking for internal and/or external
    inter, exter = MEMOPRT(file)
#check if config block some kind of memory analysis
    internal = inter and intoverride
    external = exter and extoverride
    if debug: print ('Internal mem :', internal)
    if debug: print ('External mem :', external)
#☻extracting .p0.meminfo MEM values
    if internal:
        sub, mem = USED(file)
        TIMVRS   = np.array(sub)
        MEMINT   = np.array(mem)
        MEMINT   = np.vstack(MEMINT)
    if external:
        sub, mem = DVMEM(file)
        if not internal:
            TIMVRS    = np.array(sub)
        MEMEXT = np.array(mem)
        MEMEXT = np.vstack(MEMEXT)
#loop for processors
    for p in proc:
        file     = name+'.p'+str(p)+'.meminfo'
        if internal:
            sub, mem = USED(file)
            MEMINT2  = np.array(mem)
            MEMINT2  = np.vstack(MEMINT2)
            MEMINT   = np.append(MEMINT, MEMINT2, axis=1)
        if external:
            sub, mem = DVMEM(file)
            MEMEXT2  = np.array(mem)
            MEMEXT2  = np.vstack(MEMEXT2)
            MEMEXT   = np.append(MEMEXT, MEMEXT2, axis=1)
#identification of maximum memory/core used        
    maxmemext = 0
    maxmemint = 0
    if internal: 
        maxmemint = np.amax(MEMINT)
        DFint = pd.DataFrame(MEMINT)
    if external:
        maxmemext = np.amax(MEMEXT)
        DFext = pd.DataFrame(MEMEXT)
    
#---------------------------------------------------
#                 PLOT
#---------------------------------------------------
    if makeplot:
#defining the 'call timvrs' axis
        x = np.array(range(len(TIMVRS)))
#ratio for figure sides
        ratio = len(TIMVRS)*0.4
#creating first plot figure
        plt.figure (num=i, figsize=(ratio, 4), edgecolor='k', dpi=300)
#setup number of plot figures
        mode = internal and external
        twoplots = mode and sepplots
#plotting internal values
        if internal:
            plt.plot (x, DFint[0], linewidth=1, marker=',',
                      color=iam0[0], zorder=10)
            # hack for label
            plt.plot (0,0, label='IAMEQ0int', marker=',', color=iam0[0])
            plt.plot (0,0, label='otherint', marker=',', color=altri[0])
            for p in proc:
                plt.plot (x, DFint[p], linewidth=1, marker=',',
                          color=altri[0], alpha=0.4)
#saving only internal if SEPARATEPLOTS is true
        if twoplots:
            #axis parameters
            plt.xlim   (left=1, right=len(TIMVRS)-1)
            plt.ylim   (bottom=0)
            plt.ylabel ('Max MEM (MB)', fontsize=8)
            plt.yticks (fontsize=7)
            plt.xticks (x, sub, fontsize=4, rotation = 60)
            plt.vlines (x, 0, maxmemint, alpha=0.2, color='k', linewidth=0.3)
            #red_patch = mpatches.Patch(color='black', label='Other')
            plt.legend(loc=2, labelcolor='mfc', frameon=False, prop=font)      
            caption  = day+' - '+name+' INTERNAL'
            plt.title  (caption, fontsize=8)
            #saving figure in output directory
            plt.savefig(out_dir+name+'_int'+'.png', dpi=300)
            if debug: print('Creating figure in : ', out_dir+name+'_int.png')
#creating second plot
            i+=1
            plt.figure (num=i, figsize=(ratio, 4), edgecolor='k', dpi=300)
#poltting external values
        if external:
            plt.plot (x, DFext[0], linewidth=1, marker=',',
                      color=iam0[1], zorder=10)
            # hack for label
            plt.plot (0,0, label='IAMEQ0ext', marker=',', color=iam0[1])
            plt.plot (0,0, label='otherext', marker=',', color=altri[1])
            for p in proc:
                plt.plot (x, DFext[p], linewidth=1, marker=',',
                          color=altri[1], alpha=0.4)
#saving second plot, external
        if twoplots:
            #axis parameters
            plt.xlim   (left=1, right=len(TIMVRS)-1)
            plt.ylim   (bottom=0)
            plt.ylabel ('Max MEM (MB)', fontsize=8)
            plt.yticks (fontsize=7)
            plt.xticks (x, sub, fontsize=4, rotation = 60)
            plt.vlines (x, 0, maxmemext, alpha=0.2, color='k', linewidth=0.3)
            plt.legend(loc=2, labelcolor='mfc', frameon=False, prop=font)      
            caption  = day+' - '+name+' EXTERNAL'
            plt.title  (caption, fontsize=8)
            #saving figure in output directory
            plt.savefig(out_dir+name+'_ext'+'.png', dpi=300)
            if debug: print('Creating figure in : ', out_dir+name+'_ext.png')
#saving yhe only figure if SEPARATEPLOTS is false
#or if there is only one memory analysis type
        if not twoplots:
            #axis parameters
            plt.xlim   (left=1, right=len(TIMVRS)-1)
            plt.ylim   (bottom=0)
            plt.ylabel ('Max MEM (MB)', fontsize=8)
            plt.yticks (fontsize=7)
            plt.xticks (x, sub, fontsize=4, rotation = 60)
            plt.vlines (x, 0, max(maxmemext,maxmemint), alpha=0.2,
                        color='k', linewidth=0.3)
            plt.legend(loc=2, labelcolor='mfc', frameon=False, prop=font)      
            caption  = day+' - '+name
            plt.title  (caption, fontsize=8)
            #saving figure in output directory
            plt.savefig(out_dir+name+'.png', dpi=300)
            if debug: print('Creating figure in : ', out_dir+name+'.png')






















#---------------------------------------------------
#                 TABULAR.txt
#---------------------------------------------------
    if maketab:
#creating tabular output file
            tabpath = out_dir+name+'.txt'
            if os.path.exists(tabpath):
                os.remove(tabpath)
            if debug: print('Creating tabular in: ', tabpath)
            tab = open(tabpath, 'a')
            tab.write('This file was generated with Python')
            tab.write('\nAuthor:        GIACOMO AMBROGIO')
            tab.write('\n               Univerity of Turin')
            tab.write('\nCreation date: '+day)
            tab.write('\n\nThis file combines all the informations obtained')
            tab.write('\nwith MEMOPRT keyword in Crystal calculations')
            tab.write('\n\n')
            tab.write('\n'+60*'=')
            tab.write('\n'+60*'=')
            tab.write('\nCRYSTAL JOB NAME:   '+name)
            tab.write('\n'+60*'=')
            tab.write('\n'+60*'=')
            tab.write('\n\n')
            tab.write('\n'+60*'=')
            tab.write('\nGeneral Informations:')
            tab.write('\n'+60*'=')
            tab.write('\n')
            tab.write('\nNumber of Processors   :   '+str(nproc))
            tab.write('\n')
            if internal: 
                tab.write('\nINTERNAL memory analysys detected')
                tab.write('\nMaximum Memory per core:   '+str(maxmemint)+' MB')
                tab.write('\n   Reached by processor:  ')
                for p in allproc:
                    pp  = str(p).zfill(len(str(nproc-1)))
                    hhh = DFint[p].max()
                    if hhh == maxmemint:
                        tab.write(' '+pp)
            tab.write('\n\n')
            if external: 
                tab.write('\nEXTERNAL memory analysys detected')
                tab.write('\nMaximum Memory per core:   '+str(maxmemext)+' MB')
                tab.write('\n   Reached by processor:  ')
                for p in allproc:
                    pp  = str(p).zfill(len(str(nproc-1)))
                    hhh = DFext[p].max()
                    if hhh == maxmemext:
                        tab.write(' '+pp)
            tab.write('\n\n')
            tab.write('\n'+60*'=')
            tab.write('\nProcessors Informations:')
            tab.write('\n'+60*'=')
            tab.write('\n')
            
            
            if internal:
                tab.write('\n-----------------------------------------')
                tab.write('\nINTERNAL Memory Analisys:')
                tab.write('\n-----------------------------------------')
                tab.write('\n--All memory values are expressed in MB--')
                tab.write('\n-----------------------------------------')
                tab.write('\n')
                procindex = []   # Processors index
                indcalls  = []   # Number and name of 'max Mem' call to timvrs 
                pmaxmems  = []   # Max Memory usage
                for p in allproc:
                    pp  = str(p).zfill(len(str(nproc-1)))
                    procindex.append('p'+str(pp))
                    indcall = DFint[p].idxmax()
                    pmaxmem = DFint[p].max()
                    howmanycalls = len(str(len(TIMVRS)))
                    CALLS = '['+str(indcall).zfill(howmanycalls)+'] '+str(TIMVRS[indcall])
                    indcalls.append(CALLS)
                    pmaxmems.append(str(pmaxmem))
#Dataframe creation: 
                dati = pd.DataFrame(list(zip(procindex,indcalls,pmaxmems)))
                Indici = ['Proc', 'TIMVRS call', 'MaxMEM']
                tab.write(tabulate(dati, headers=Indici, tablefmt='github'))
                tab.write('\n')
            if external:
                tab.write('\n-----------------------------------------')
                tab.write('\nEXTERNAL Memory Analisys:')
                tab.write('\n-----------------------------------------')
                tab.write('\n--All memory values are expressed in MB--')
                tab.write('\n-----------------------------------------')
                tab.write('\n')
                procindex = []   # Processors index
                indcalls  = []   # Number and name of 'max Mem' call to timvrs 
                pmaxmems  = []   # Max Memory usage
                BaseLine  = []   # VmPeak Baseline
                Pids      = []   # Pids
                for p in allproc:
                    pp  = str(p).zfill(len(str(nproc-1)))
                    procindex.append('p'+str(pp))
                    indcall = DFext[p].idxmax()
                    pmaxmem = DFext[p].max()
                    howmanycalls = len(str(len(TIMVRS)))
                    CALLS = '['+str(indcall).zfill(howmanycalls)+'] '+str(TIMVRS[indcall])
                    indcalls.append(CALLS)
                    pmaxmems.append(str(pmaxmem))
#for other info in .meminfo
                    file   = name+'.p'+str(p)+'.meminfo'
                    VmPeak = SYSMEM(file)
                    BaseLine.append(VmPeak[0])
                    Pid = PID(file)
                    Pids.append(Pid)
#Dataframe creation: 
                dati = pd.DataFrame(list(zip(procindex,Pids,indcalls,pmaxmems,BaseLine)))
                Indici = ['Proc', 'PID', 'TIMVRS call', 'MaxMEM', 'BaseLine']
                tab.write(tabulate(dati, headers=Indici, tablefmt='github'))
                tab.write('\n\n\n')
            tab.write('\n'+60*'=')
            tab.write('\nAll Memory data acquired:')
            tab.write('\n'+60*'=')
            tab.write('\n\n')
            if internal:
                tab.write('\n--------')
                tab.write('\nINTERNAL')
                tab.write('\n--------')
                tab.write('\n')
                DFint.set_index(TIMVRS, inplace=True, append=True)
                tab.write(pd.DataFrame.to_string(DFint, header=procindex, index=True))
                tab.write('\n\n')
            if external:
                tab.write('\n--------')
                tab.write('\nEXTERNAL')
                tab.write('\n--------')
                tab.write('\n')
                DFext.set_index(TIMVRS, inplace=True, append=True)
                tab.write(pd.DataFrame.to_string(DFext, header=procindex, index=True))
            tab.write('\n\n')
            tab.write('\n'+60*'=')
            tab.write('\nTermination')
            tab.write('\n'+60*'=')
            tab.close()
# Special thanks to Chiara Ribaldone who wrote all the MEMOPRT-related
# subroutines that made theese extended memory analysis possible
#---------------------------------------------------
#                 END
#---------------------------------------------------


print('...')
print('...')
print('TERMINATION')

















    