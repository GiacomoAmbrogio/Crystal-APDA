# Crystal-APDA
Crystal- All Process Data Analysis - Python
------------------------------------------------------------------------
Author: Giacomo AMBROGIO
        @University of Turin
        2022
________________________________________________________________________
This python tool is used to manage all output files generated by crystal
calculations using the keyword 'MEMOPRT' with 'ALLPINT' and/or 'ALLPEXT'
options.


#   SETUP
Simply put APDA.py in a folder and run it.
It will automatically generate the confing.ini file (default settings)
and input-output folders. Please note that these folder must stay in the
same directories of the APDA.py tool.


#   USE
Drop in the input folder all the .pX.meminfo files from crystal and run
APDA.py to generate plots and .txt file with all the datas combined.

Note that an unlimited number of different crystal calculations can be
processed (as long as they have different names), but you must always
have at least one .p0.meminfo and .p1.meminfo file (at least 2 cores)
for any calculations.

The program will automatically detect if 'ALLPINT' and/or 'ALLPEXT'
oprions have been used, but it is possible to disable the analysis of
these memory results via the config.ini file (see next section).


#   CONFIGURATION
Using the config.ini file you can modify the behavor of the program:

APDA: turn on/off features of the APDA tool
    PLOT     - making of flot images
    TABULAR  - making of .txt outputs file
    INTERNAL - analyse internal memory
    EXTERNAL - analyse external memory
    
DIRECTORIES: select names for directories

PLOT: customize plot outputs
    SEPARATEPLOTS     - if True APDA will generate two different plots
                        (internal and external memories)
    IAMEQ0_COLOR      - select colors for the first processors in plots
    OTHER_PROCS_COLOR - select colors for other processors




#   CONTACTS
For any problems, bug reports or feature requests please contact:
   GIACOMO AMBROGIO
   University of Turin
   giacomo.ambrogio@edu.unito.it
