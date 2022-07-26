#!/usr/bin/python
print "-----------------------------------------------"
print "Hello! I am FluroStronger"
print "-----------------------------------------------"
print "I take Bruker fids, process them using nmrPipe"
print "Alternatively, save the spectrum as a text file"
print "and ignore the processing parts of the code."
print "Then take the data and fit to a series of singlets"
print "and quartets."
print "I hope this works for you. Email us if not."
print "-----------------------------------------------"
print "A.Baldwin July 2021"
print "(c) University of Oxford"
print 

import sys,os   #I need these standard modules.
sys.path.append('./py')   #adding py to the path (this contains the source code)
from fit import fitter    #bring in a fitting module

PROCESS=False #True  #set to true of fitting
############################################################
if(PROCESS):
    from anal import process  #bring in a process module 

    #O1 takes 'temp', 'guess' or a value for centre of spectrum
    wf=15       #line broadening. (lb)
    o1 ='guess' #figure out a sensible o1 from temperature using acqus. 
    p0='auto'   #autophase p0. pretty slow: designed to handle nasty fluorine baselines.
    #p1='auto'   #autophase p1

    pars={}


    #parameters for the zoom in.
    pars={}  #zoom in and ignore big peak.
    pars['p1']=120 #starting guess for phasing.
    pars['base']=((-65,-70,-85,-90,-111,-108),50,False)

    expts=[]
    mode=[]

    #if you omit the '707' it will try and find all 1D data in bruker folders for processing.
    expts.append('./CF2Lys(Me3)_GdHCl/707') ;mode.append('e') 


    import os,subprocess 
    for i,expt in enumerate(expts): #run over array of experiments. look for Bruker files.

        process(expt,'NMR','find',700,(1.92,3.53),p0,pars['p1'],wf,ws=8,O1=o1,base='lin',basePars=pars['base'],amx=True)
        print expt
        bef=os.getcwd()
        os.chdir(expt)
        subprocess.call("pipe2txt.tcl test.ft2 -index ppm > test.txt", shell=True, executable='/bin/csh') #turn processed data into a text file.
        os.chdir(bef)
        
    sys.exit(100)   #quit.


    
#run fits to make strong coupling figure fit.



filo='./CF2Lys(Me3)_GdHCl/707/test.txt'  #output file.

#the lists are initial conditions for the parameters.
#if initial conditions are reasonably good, the fit will have no problem.
spectrumRange=(-105,-70)  #low and high ppms to fit

#singlets:  #keep adding entries to the lists as required
singletChemShifts=(-75.5,-101.2,-98.9)   #initial chemical shifts
singletIntensities=(0.1,0.2,0.1)         #initial relative intensities
singletLinewidths=(300,400,600)          #initial linewidths (deliberately set these to be broad to help fitting)

#quartets:  #keep adding entries to the lists as required
quartetChemShifts=((-99.6,-98.0),(-97.9,-98.6))  #locations of each pair of spins
quartetLinewidths=(200,200)              #initial linewidths
quartetJvalues=(250,200)                 #initial scalar coupling constants

frq=564.680          #operating frequency of nucleus of interest (MHz)

mainQuartetIntensity=3E9  #intensity of the first quartet (Sets overall scale for all - relative intensity of this guy will be 1)
otherQuartetIntensities=(1,) #relative intensity of remaining quartets

jiggles=False  #set to true to try and make the fit get out of a local minimum (slow)
fit=True       #set to false to just simulate the data, to help you figure out what the initial conditions should be

inst=fitter(filo,'CF2lysMe3Gd',spectrumRange,mainQuartetIntensity,otherQuartetIntensities,singletIntensities,quartetChemShifts,singletChemShifts,quartetLinewidths,singletLinewidths,quartetJvalues,frq,fit=fit,jiggles=jiggles)



