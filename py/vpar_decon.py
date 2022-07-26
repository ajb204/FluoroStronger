#!/usr/bin/python
#####################################################
# Functions to work out an fid.com conversion script
# Use at your peril.

import os,sys
from scipy.optimize import leastsq
import nmrglue as ng
import numpy

#A.Baldwin July 2021
#(c) University of Oxford


class vpar():
    def __init__(self,outdir,dim,labb,rk,phase,p1,wf,nuslist='',o1p='auto',quant=False,ws=4,base=False,amx=False,basePars=((9.5,10,-3.,-2),5,0.1)):
        self.outdir=outdir
        self.quant=quant
        self.phase=phase
        self.p1=p1
        self.wf=wf


        self.ws=ws
        self.base=base
        self.amx=amx

        
        
        #self.p1='auto'
        #self.phase='auto'

        #unpack basepars
        #if self.phase is auto, this will be used.
        self.posNoise=numpy.min(basePars[0][:2]),numpy.max(basePars[0][:2])
        if(len(basePars[0])!=4):
            self.midNoise=[]
            print len(basePars[0])/2-2
            for i in range(len(basePars[0])/2-2):
                self.midNoise.append((numpy.min(basePars[0][2+2*i:2*i+4]),numpy.max(basePars[0][2+2*i:2*i+4])))
            self.negNoise=numpy.min(basePars[0][-2:]),numpy.max(basePars[0][-2:])
        else:
            self.midNoise=False
            self.negNoise=numpy.min(basePars[0][2:4]),numpy.max(basePars[0][2:4])
        self.LB=basePars[1] #used for autophasing
        self.centreFilt=basePars[2]



        
        #self.posNoise=9.5,10  #range of positive noise
        #self.negNoise=-3,-2   #range of negative noise
        #self.centreFilt=0.1   #region to screen out.

        
        #print outdir
        #print os.path.exists(os.path.join(outdir,'ser'))

        #possible names of raw fid file
        testfiles=[]
        testfiles.append('fid')
        testfiles.append('ser')
        testfiles.append('fid.gz')
        testfiles.append('ser.gz')
        print testfiles
        tick=0
        for test in testfiles:
            print test,os.path.join(outdir,test),os.path.exists(os.path.join(outdir,test))==True
            if(os.path.exists(os.path.join(outdir,test))==True):
               tick=1
               infid=os.path.join(outdir,test)
               break

        if('ser.gz' in os.listdir(self.outdir)):
            self.GZ='ser'
        elif('fid.gz' in os.listdir(self.outdir)):
            self.GZ='fid'
        else:
            self.GZ=False

               
        if(tick==0):
            print 'Cannot find fid. Aborting.'
            
            return 

        if(infid[-3:]=='.gz'): #fid is zipped up...
            #unzipping.
            os.system('gunzip '+infid)
               

               
               
        if(os.path.exists(os.path.join(outdir,'procpar'))==0 and os.path.exists(os.path.join(outdir,'acqus'))==0):
            print 'Cannot find procpar (varian) or acqus (bruker)'
            print 'in folder',outdir
            return
        self.dim=dim
        self.labb=labb
        self.rk=rk
        self.nuslist=nuslist
        self.o1p=o1p


        self.abort=0
        self.GetSpectrometerType(path=outdir)
        self.GetSequence()
        self.Convert()
        if(self.abort==1):
            print 'Conversion aborted - errors encountered'
            return

        self.PipeParse()
        os.system('csh '+self.outdir+'/fid.test.com')
        
        if(self.phase=='auto'): #run autophase
            self.autoPhase()



        self.nmrPipe()
        os.system('csh '+self.outdir+'/nmrproc.test.com')

    def TransPhase(self):

        Fft=numpy.fft.fftshift(numpy.fft.fft(self.sog))  #Fourier Transform

        p0 = self.p0_0 * numpy.pi / 180.  #convert to radians
        p1 = self.p1_0 * numpy.pi / 180.  #convert to radians
        apod = numpy.exp(1.0j * (p0 + ((p1*1.) * self.phaseArr )))
        self.Fft=numpy.flip(apod)*Fft #apply phase correction

        #now baseline.
        #self.Fft=ng.proc_bl.baseline_corrector(self.Fft,wd=300)
        #base=numpy.median(self.Fft[self.baseMask].real)
        #self.Fft-=base
        
        if(self.midNoise==False):
            x2=numpy.average(self.posNoise)
            x1=numpy.average(self.negNoise)

            y2= numpy.average(self.Fft[(self.Frq>numpy.min(self.posNoise) * (self.Frq<numpy.max(self.posNoise)))].real)
            y1= numpy.average(self.Fft[(self.Frq>numpy.min(self.negNoise)) * (self.Frq<numpy.max(self.negNoise))].real)

            m=(y2-y1)/(x2-x1)
            c=y2-m*x2
            yvals=self.Frq*m+c #get baseline
        else:
            yvals=numpy.zeros_like(self.Frq)
            for i in range(len(self.xrefs)): #do each region, one at a time.
                ya= numpy.average(self.Fft[self.maskref[i][0]].real) #most positive
                yb= numpy.average(self.Fft[self.maskref[i][1]].real) #least positive
                xa=self.xrefs[i][0]  #most positive
                xb=self.xrefs[i][1]  #least positive
                m1=(ya-yb)/(xa-xb)   #gradient
                c1=yb-m1*xb

                #print 'xa',xa,'xb',xb,'ya',ya,'yb',yb
                #print self.Frq[self.maskref[i][0]]
                #print self.Frq[self.maskref[i][1]]

                #c2=ya-m1*xa
                #print c1,c2
                yvals[self.regiref[i]]=self.freqref[i]*m1+c1 #calc baseline for region
            
        #outy=open('tmp.out','w')
        #for i in range(len(self.Fft)):
        #    outy.write('%f\t%f\t%e\t%e\n' % (self.Frq[i],yvals[i],self.Fft[i],(self.Fft-yvals)[i]))
        #outy.close()
        #sys.exit(100)
        self.Fft-=yvals    #subtract baseline

        
    def ChiPhase(self,x):
        self.unpackPhase(x)
        #print 'chi:',self.p0_0,self.p1_0
        self.TransPhase()
        pos=numpy.sum(self.Fft[(self.Fft.real>0)*self.PhaseMask].real)
        neg=(numpy.sum(self.Fft[(self.Fft.real<0)*self.PhaseMask].real*-1)) #find the minimum in this.

        #return neg,neg
        return neg,1/pos #,1/pos
        
    def packPhase(self):
        x=[]
        x.append(self.p0_0)
        if(self.p1=='auto'):
            x.append(self.p1_0)
        return x
            
    def unpackPhase(self,x):
        self.p0_0=x[0]
        if(self.p1=='auto'):
            self.p1_0=x[1]


    def DoWalkP0(self,phbest0,phmax0,grids):
        print 'phase:',self.p0_0,self.p1_0 #'%.4e' % self.ChiPhase((self.p0_0,self.p1_0))[0]
        for k in range(3): #repeat thrice.
            print 'walking p0',phbest0,'+/-',phmax0,self.p1_0
            phaseTest0=numpy.linspace(phbest0-phmax0,phbest0+phmax0,grids)
            #print phaseTest0
            grid=[]
            phs=[]
            for i,p in enumerate(phaseTest0):
                q=self.p1_0
                opt=self.ChiPhase((p,q))
                #outy.write('%e\t%e\t%e\t%e\n' % (p,q,opt[0],opt[1]))
                grid.append(opt[0])
                phs.append((p,q))
                #print x0
                
            grid=numpy.array(grid)
            phs=numpy.array(phs)

            self.p0_0=phs[numpy.argmin(grid)][0]
            self.p1_0=phs[numpy.argmin(grid)][1]

            phbest0=self.p0_0
            #phbest1=self.p1_0
            
            phmax0=phmax0/10
            #phmax1=phmax1/10

        print 'best p0 found:',self.p0_0,'%.4e' % numpy.min(grid)
        return self.p0_0,numpy.min(grid)
            
    def DoWalkP1(self,phbest1,phmax1,grids):
        print 'phase:',self.p0_0,self.p1_0 #,'%.4e' % self.ChiPhase((self.p0_0,self.p1_0))[0]
        for k in range(3): #repeat thrice.
            print 'walking p1',phbest1,'+/-',phmax1
            phaseTest1=numpy.linspace(phbest1-phmax1,phbest1+phmax1,grids)
            grid=[]
            phs=[]
            #for i,p in enumerate(phaseTest0):
            p=self.p0_0
            for j,q in enumerate(phaseTest1):
                opt=self.ChiPhase((p,q))
                #outy.write('%e\t%e\t%e\t%e\n' % (p,q,opt[0],opt[1]))
                #print 'go:',p,q,'%.6e' % opt[0],opt[0]
                grid.append(opt[0])
                phs.append((p,q))
                #print x0
                #outy.write('\n')

            grid=numpy.array(grid)
            phs=numpy.array(phs)

            self.p0_0=phs[numpy.argmin(grid)][0]
            self.p1_0=phs[numpy.argmin(grid)][1]

            phbest0=self.p0_0
            phbest1=self.p1_0
            
            #phmax0=phmax0/10
            phmax1=phmax1/10
            #outy.write('\n\n')
        print 'best p1 found:',self.p1_0,'%.4e' % numpy.min(grid)
        return self.p1_0,numpy.min(grid)
            
    def DoWalkGrid(self):
        for k in range(3): #repeat thrice.
            phaseTest0=numpy.linspace(phbest0-phmax0,phbest0+phmax0,20)
            print phaseTest0
            if(self.p1=='auto'):
                phaseTest1=numpy.linspace(phbest1-phmax1,phbest1+phmax1,20)
            grid=[]
            phs=[]

            for i,p in enumerate(phaseTest0):
                if(self.p1=='auto'):
                    for j,q in enumerate(phaseTest1):
                        opt=self.ChiPhase((p,q))
                        outy.write('%e\t%e\t%e\t%e\n' % (p,q,opt[0],opt[1]))
                        grid.append(opt[0])
                        phs.append((p,q))
                        #print x0
                    outy.write('\n')
                else:
                    q=self.p1_0
                    opt=self.ChiPhase((p,q))
                    outy.write('%e\t%e\t%e\t%e\n' % (p,q,opt[0],opt[1]))
                    grid.append(opt[0])
                    phs.append((p,q))
                    #print x0
                

            grid=numpy.array(grid)
            phs=numpy.array(phs)

            self.p0_0=phs[numpy.argmin(grid)][0]
            self.p1_0=phs[numpy.argmin(grid)][1]

            phbest0=self.p0_0
            phbest1=self.p1_0
            
            phmax0=phmax0/10
            phmax1=phmax1/10
            outy.write('\n\n')

            
    def SetBaselineBoundary(self):
        print 'setting baseline boundaries..'
        x2=numpy.average(self.posNoise)  #positive extrema
        xm=numpy.average(self.midNoise[0]) #first boundary from positive
            
        self.xrefs=[]
        self.regiref=[]
        self.freqref=[]
        
        self.xrefs.append((x2,xm))       #between first boundary and positive extrema
        self.regiref.append(self.Frq>=xm) #frequencies from boundary average to extrema
        self.freqref.append(self.Frq[self.regiref[-1]])
        print self.xrefs[-1],self.regiref[-1],self.freqref[-1]
        
        for i in range(len(self.midNoise)-1): #loop over interior boundaries
            xb=numpy.average(self.midNoise[i])   #positive
            xa=numpy.average(self.midNoise[i+1]) #less positive
            self.xrefs.append((xb,xa))          #most positive, then less positive
            self.regiref.append((self.Frq<xb)*(self.Frq>=xa)) #region is between these two.
            self.freqref.append(self.Frq[self.regiref[-1]])
            print self.xrefs[-1],self.regiref[-1],self.freqref[-1]
                
        x1=numpy.average(self.negNoise)       #negative extrema
        xm=numpy.average(self.midNoise[-1])   #least positive boundary
        self.xrefs.append((xm,x1))            #most positive, then less positive
        self.regiref.append((self.Frq<xm ))
        self.freqref.append(self.Frq[self.regiref[-1]])
        print self.xrefs[-1],self.regiref[-1],self.freqref[-1]

            
        self.maskref=[]
        self.maskref.append(( (self.Frq>numpy.min(self.posNoise)) * (self.Frq<numpy.max(self.posNoise)),(self.Frq>numpy.min(self.midNoise[0])) * (self.Frq<numpy.max(self.midNoise[0])))) #most positive, then less positive.
        print numpy.min(self.posNoise),numpy.max(self.posNoise),numpy.min(self.midNoise[0]),numpy.max(self.midNoise[0])
        for i in range(len(self.midNoise)-1):
            self.maskref.append(((self.Frq>numpy.min(self.midNoise[i])) * (self.Frq<numpy.max(self.midNoise[i])),(self.Frq>numpy.min(self.midNoise[i+1])) * (self.Frq<numpy.max(self.midNoise[i+1]))))
            print numpy.min(self.midNoise[i]),numpy.max(self.midNoise[i]),numpy.min(self.midNoise[i+1]),numpy.max(self.midNoise[i+1])
        self.maskref.append(((self.Frq>numpy.min(self.midNoise[-1])) * (self.Frq<numpy.max(self.midNoise[-1])),(self.Frq>numpy.min(self.negNoise)) * (self.Frq<numpy.max(self.negNoise))))
        print numpy.min(self.midNoise[-1]),numpy.max(self.midNoise[-1]),numpy.min(self.negNoise),numpy.max(self.negNoise)
        print 'done'

        #for i in range(len(self.maskref)):
        #    print self.Frq[self.maskref[i][0]],self.Frq[self.maskref[i][1]]
        #print numpy.min(self.posNoise),numpy.max(self.posNoise)
        #print 'ff',self.Frq[(self.Frq>numpy.min(self.posNoise)) * (self.Frq<numpy.max(self.posNoise))]

            
    def DoAutoPhase(self):

        if(self.midNoise!=False):
            self.SetBaselineBoundary()

        
        
        self.p0_0=0 #self.phase
        #self.p0_0=0#330
        #self.p1=1028
        #self.p1=1024


        self.p1_0=self.p1
        self.p1='auto' #need to set to auto for p1 optimisation.

        #330 1236
        
        phbest0=self.p0_0  #if gridding on p0
        phmax0=180 #if gridding on p0
        phbest1=self.p1_0  #if gridding on p1
        phmax1=2000  #if gridding on p1



        phbest0,curr=self.DoWalkP0(phbest0,phmax0,21)
        phbest1,curr=self.DoWalkP1(phbest1,phmax1,21)

        #phmax0=phmax0/10
        #phmax1=phmax1/10
        
        phbest0,curr=self.DoWalkP0(phbest0,phmax0,31)

        phmax1=phmax1/50
        phbest1,curr=self.DoWalkP1(phbest1,phmax1,31)




        
        #first gridsearch over p0 to get minimum.

        """
        print self.phase,self.p1

        outy=open('min.opt','w')
        for k in range(3): #repeat thrice.

            phaseTest0=numpy.linspace(phbest0-phmax0,phbest0+phmax0,20)
            print phaseTest0
            if(self.p1=='auto'):
                phaseTest1=numpy.linspace(phbest1-phmax1,phbest1+phmax1,20)
            grid=[]
            phs=[]

            for i,p in enumerate(phaseTest0):
                if(self.p1=='auto'):
                    for j,q in enumerate(phaseTest1):
                        opt=self.ChiPhase((p,q))
                        outy.write('%e\t%e\t%e\t%e\n' % (p,q,opt[0],opt[1]))
                        grid.append(opt[0])
                        phs.append((p,q))
                        #print x0
                    outy.write('\n')
                else:
                    q=self.p1_0
                    opt=self.ChiPhase((p,q))
                    outy.write('%e\t%e\t%e\t%e\n' % (p,q,opt[0],opt[1]))
                    grid.append(opt[0])
                    phs.append((p,q))
                    #print x0
                

            grid=numpy.array(grid)
            phs=numpy.array(phs)

            self.p0_0=phs[numpy.argmin(grid)][0]
            self.p1_0=phs[numpy.argmin(grid)][1]

            phbest0=self.p0_0
            phbest1=self.p1_0
            
            phmax0=phmax0/10
            phmax1=phmax1/10
            outy.write('\n\n')
        outy.close()
        """

        print self.p0_0,self.p1_0

        self.ChiPhase(self.packPhase())

        #self.p1='auto'
        
        print 'grid:',self.p0_0,self.p1_0
        #now optimise, including p1 if desired.

        outy=open('test.out','w')
        print self.o1p
        for i in range(len(self.Frq)):
            outy.write('%e\t%e\n' % (self.Frq[i],self.Fft[i]))
        outy.close()
        
        x0=leastsq(self.ChiPhase,self.packPhase())
        #x0=(self.p0_0,self.p1_0),0
        print x0[0]
        print 'after optimisation:'
        print 'p0:',self.p0_0
        print 'p1:',self.p1_0

        #update globals for nmrproc.com.
        self.p1=self.p1_0
        self.phase=self.p0_0


        print len(self.PhaseData)
        outy=open('test.out','w')
        print self.o1p
        for i in range(len(self.Frq)):
            outy.write('%e\t%e\n' % (self.Frq[i],self.Fft[i]))
        outy.close()

    def autoPhase(self):
        #read in fid and figure out phase.

        #self.LB=5  #set window function.


        
        if(self.tp=='bruk'):
            dic,doba=ng.bruker.read('./')
            self.PhaseData = ng.bruker.remove_digital_filter(dic,doba)
        else:
            dic,self.PhaseData = ng.pipe.read('test.fid') #read fids

        
        self.PhaseData[0]*=0.5 #pre-multiply first point by 0.5

        #window function.
        apod = numpy.exp(-numpy.pi * numpy.arange(len(self.PhaseData))*1. * self.LB /self.sw)#.astype(data.dtype)
        self.sog=apod*self.PhaseData #apodised FID.
        self.phaseArr=(numpy.arange(len(self.sog))*1. / (len(self.sog)*1.)) #will use this to set phase.



        import matplotlib
        Frq=numpy.fft.fftfreq(len(self.PhaseData),d=1/self.sw*self.sfrq) #*-1 #Get frequency range of fourier
        self.Frq=numpy.fft.fftshift(Frq)+self.waterppm                         #set indirect for plotting/

        """
        self.LB=200
        size = len(self.PhaseData)
        apod = numpy.exp(-numpy.pi * numpy.arange(size)*1. * self.LB /self.sw)#.astype(data.dtype)
        sog=apod*self.PhaseData
        Fft=numpy.fft.fft(sog)                           #Fourier Transform
        Fft=numpy.fft.fftshift(Fft)
        ndata,(p0,p1)=ng.process.proc_autophase.autops(Fft,'acme',p0=0.0,p1=1024,return_phases=True,peak_width=400)


        #ng.process.proc_autophase.manual_ps(Fft)

        print ndata,p0,p1
        
        outy=open('test.junk','w')
        for i in range(len(self.PhaseData)):
            outy.write('%e\t%e\n' % (self.Frq[i],ndata[i]))
        outy.close()
        sys.exit(100)
        """
        
        #phaseMask needs to include all intensity, including dispersive intensity.
        #basemask needs to be as out of the way as possible, exclude smiles if there are any.

        
        #self.PhaseMask=((self.Frq>self.waterppm+0.1)+(self.Frq<self.waterppm-0.1))*(self.Frq>-1)*(self.Frq<10)
        #self.baseMask=((self.Frq>self.waterppm+0.1)+(self.Frq<self.waterppm-0.1))*(self.Frq>-2)*(self.Frq<12)

        #screen the centre, and take a region that only has resonances for phasing
        

        if(self.centreFilt!=False):
            self.PhaseMask=((self.Frq>self.waterppm+self.centreFilt)+(self.Frq<self.waterppm-self.centreFilt))*(self.Frq>numpy.max(self.negNoise))*(self.Frq<numpy.min(self.posNoise))
        else:
            self.PhaseMask=(self.Frq>numpy.max(self.negNoise))*(self.Frq<numpy.min(self.posNoise))
        #self.baseMask=((self.Frq>self.waterppm+self.centreFilt)+(self.Frq<self.waterppm-self.centreFilt))*(self.Frq>-2)*(self.Frq<12) #no longer used. 


        self.DoAutoPhase() #run autophasing.
        

        
        
    def nmrPipe(self):

        outy=open(self.outdir+'/nmrproc.test.com','w')
        outy.write('#!/bin/csh\n')
        outy.write('nmrPipe -in test.fid \\\n')
        #outy.write('#| nmrPipe -fn SOL                            \\\n')
        outy.write('| nmrPipe  -fn EM  -lb %f -c 0.5              \\\n' % (self.wf))
        outy.write('| nmrPipe  -fn ZF -auto                        \\\n')
        outy.write('| nmrPipe  -fn FT -auto                        \\\n')
        outy.write('| nmrPipe  -fn PS -p0 %f -p1 %f -di -verb  \\\n' % (self.phase, self.p1))
        if(self.base=='poly'):
            outy.write('| nmrPipe -fn POLY -auto                     \\\n')
        elif(self.base=='lin'):
            outy.write('| nmrPipe -fn EXT -x1 %fppm -xn %fppm -sw     \\\n' % (numpy.min(self.negNoise),numpy.max(self.posNoise)))
            if(self.midNoise==False):
                outy.write('| nmrPipe -fn BASE -nw 5 -nl %fppm %fppm %fppm %fppm             \\\n' % (numpy.min(self.negNoise),numpy.max(self.negNoise),numpy.min(self.posNoise),numpy.max(self.posNoise)   ))
            else:
                baseStr='%fppm %fppm ' % (numpy.max(self.posNoise),numpy.min(self.posNoise))
                for i in range(len(self.midNoise)):
                    baseStr+='%fppm %fppm ' % (numpy.max(self.midNoise[i]),numpy.min(self.midNoise[i]))
                baseStr+='%fppm %fppm ' % (numpy.max(self.negNoise),numpy.min(self.negNoise))
                    
                outy.write('| nmrPipe -fn BASE -nw 5 -nl %s  \\\n' % (baseStr))

            #outy.write('| nmrPipe -fn BASE -nw 2 -nl 0%s 5%s 95%s 100%s                \\\n' % ('%','%','%','%'))
            #outy.write('| nmrPipe -fn BASE -nw 5 -nl -2ppm -1ppm 9.5ppm 10ppm              \\\n')

            #outy.write('#| nmrPipe  -fn TP                            \\\n')
        #outy.write('#| nmrPipe  -fn LP -fb                        \\\n')
        #outy.write('#| nmrPipe  -fn EM  -lb 0.0 -c 1.0            \\\n')
        #outy.write('#| nmrPipe  -fn ZF -auto                      \\\n')
        #outy.write('#| nmrPipe  -fn FT -auto                      \\\n')
        #outy.write('#| nmrPipe  -fn PS -p0 0 -p1 0.00 -di -verb   \\\n')
        #outy.write('#| nmrPipe -fn TP                             \\\n')
        #outy.write('#| nmrPipe -fn POLY -auto                     \\\n')
        outy.write('   -ov -out test.ft2\n')
        outy.close()


    def GetTempBruk(self):
        inny=open(self.outdir+'/vtc_pid_settings')
        for line in inny.readlines():
            if('temperature' in line):
                self.temp=float(line.split('"')[1])-273.19

    def Convert(self):
        if(self.tp=='var'):
            self.ConvertVarian()
        elif(self.tp=='bruk'):
            self.ConvertBruker()
        elif(self.tp=='omeg'):
            self.ConvertOmega()

    # Based on Patrik Lundstrom 011126
    #take water and sfrq, calc ppms of C and N
    def shift(self,dfrq,nuc='C13'):
        if(nuc=='H1'):
            CONV=1.
        if(nuc=='C13'):
            CONV=0.251449530
        if(nuc=='N15'):
            CONV=0.101329118
        if(nuc=='P31'):
            CONV=0.4048064954
        if(nuc=='F19'):
            CONV=0.9412866605363297

        sfrq0  = self.sfrq / (1.0 + self.waterppmTOF*1e-6);
        dfrq0 = sfrq0*CONV;
        ppm = (dfrq-dfrq0)/dfrq0*1e6;
        return ppm



    def GetNuc(self,lab):
        """
        if(nuc=='H1'):
            CONV=1.
        if(nuc=='C13'):
            CONV=0.251449530
        if(nuc=='N15'):
            CONV=0.101329118
        if(nuc=='P31'):
            CONV=0.4048064954
        if(nuc=='F19'):
            CONV=0.9412866605363297
        """
        print lab,lab[0],self.dfrq/self.sfrq,self.dfrq/self.sfrq-0.25

        if(lab[0]=='C'):
            if(numpy.fabs(self.dfrq/self.sfrq-0.25)<1E-2):
                return self.dfrq,self.dn
            if(numpy.fabs(self.dfrq2/self.sfrq-0.25)<1E-2):
                return self.dfrq2,self.dn2
        if(lab[0]=='N'):
            if(numpy.fabs(self.dfrq/self.sfrq-0.1)<1E-2):
                return self.dfrq,self.dn
            if(numpy.fabs(self.dfrq2/self.sfrq-0.1)<1E-2):
                return self.dfrq2,self.dn2
        if(lab[0]=='H'):
            if(numpy.fabs(self.sfrq/self.sfrq-1)<1E-2):
                return self.sfrq,self.tn
        print 'neither sfrq,dfrq nor dfrq2 seem to correspond to label',lab
        self.abort=1
        return -1

    def ConvertBruker(self):

        self.ns=self.GetParBruk('acqus',('','NS'))[0]
        #N=self.np,self.ni*2,self.ni2*2,self.ni3*2
        #T=self.np/2,self.ni,self.ni2,self.ni3
        self.np=int(self.GetParBruk('acqus',('','TD'))[0])

        
        if(self.dim>=3):
            if(os.path.exists(os.path.join(self.outdir,'acqu3s'))==0):
                print 'No acqu3s'
                aq3=False
            else:
                aq3=True
        else:
            aq3=True
               
        if(self.nuslist==''):
            if(aq3): #get from acquX or L345
                if(self.dim>1):
                    self.ni=int(self.GetParBruk('acqu2s',('','TD'))[0])/2
                    self.ni_TD=int(self.GetParBruk('acqu2s',('','TD'))[0])
                    if(self.dim>=3):
                        self.ni2=int(self.GetParBruk('acqu3s',('','TD'))[0])/2
                        if(self.dim==4):
                            self.ni3=int(self.GetParBruk('acqu4s',('','TD'))[0])/2
            else: #to handle Marius Clores' 4D
                self.ni=int(self.GetParBruk('acqus',('','L3'))[0])
                self.ni2=int(self.GetParBruk('acqus',('','L4'))[0])
                if(self.dim==4):
                    self.ni3=int(self.GetParBruk('acqus',('','L5'))[0])
                   
        else:
            self.ni=int(self.GetParBruk('acqu2s',('','NusTD'))[0])/2
            if(self.dim>=3):
                self.ni2=int(self.GetParBruk('acqu3s',('','NusTD'))[0])/2
                if(self.dim==4):
                    self.ni3=int(self.GetParBruk('acqu4s',('','NusTD'))[0])/2


        if(self.dim>1):
            self.yMode=self.GetParBruk('acqu2s',('','FnMODE'))[0]
            if(self.dim>=3):
                if(aq3):
                    self.zMode=self.GetParBruk('acqu3s',('','FnMODE'))[0]
                    if(self.dim==4):
                        self.aMode=self.GetParBruk('acqu4s',('','FnMODE'))[0]
                else: #deal with Marius Clore's 4D
                    self.zMode=str(0)
                    if(self.dim==4):
                        self.aMode=str(0)
                
        if(self.dim>1):
            self.aqseq=self.GetAcqseq()
            print 'aqseq:',self.aqseq
            if(self.dim==3):
                if(self.aqseq=='312'):
                    print 'Swapping ni2 and ni'
                    self.ni2,self.ni=self.ni,self.ni2
                    self.yMode,self.zMode=self.zMode,self.yMode
                else: #self.aqseq='321'
                    pass

            self.mode=[]
            self.mode.append(self.yMode)
            if(self.dim>=3):
                self.mode.append(self.zMode)
                if(self.dim==4):
                    self.mode.append(self.aMode)


        #sw=self.sw,self.sw1,self.sw2,self.sw3
        #try:
        #    self.sw=float(self.GetParBruk('acqus',('','SW_h'))[0])

        #except:
        self.sw=float(self.GetParBruk('acqus',('','SW_h'))[0])
        self.SFO1=float(self.GetParBruk('acqus',('','SFO1'))[0])  #typically H
        self.BF1=float(self.GetParBruk('acqus',('','BF1'))[0])  #typically H
        self.O1=float(self.GetParBruk('acqus',('','O1'))[0])  #typically H
        NUC1=self.GetParBruk('acqus',('','NUC1'))[0].split('<')[1].split('>')[0]


        self.NUC1=NUC1[-1:]+NUC1[:-1]

        
        if(self.NUC1=='F19'):
            self.labb='F19'
        #print self.NUC1,NUC1
        #sys.exit(100)

        self.frq=self.SFO1

        if(self.dim>1):
            self.sw1=(float(self.GetParBruk('acqu2s',('','SW_h'))[0]))
            self.SFO2=float(self.GetParBruk('acqu2s',('','SFO1'))[0])  #typically C
            self.BF2=float(self.GetParBruk('acqu2s',('','BF1'))[0])  #typically C
            self.O2=float(self.GetParBruk('acqu2s',('','O1'))[0])  #typically C
            NUC2=self.GetParBruk('acqu2s',('','NUC1'))[0].split('<')[1].split('>')[0]
            self.NUC2=NUC2[-1:]+NUC2[:-1]
            self.frq1=self.SFO2
        
        if(self.dim>=3):
            if(aq3):
                self.sw2=(float(self.GetParBruk('acqu3s',('','SW_h'))[0]))
                self.SFO3=float(self.GetParBruk('acqu3s',('','SFO1'))[0]) #typically N
                self.BF3=float(self.GetParBruk('acqu3s',('','BF1'))[0]) #typically N
                self.O3=float(self.GetParBruk('acqu3s',('','O1'))[0]) #typically N
                NUC3=self.GetParBruk('acqu3s',('','NUC1'))[0].split('<')[1].split('>')[0]            
                self.NUC3=NUC3[-1:]+NUC3[:-1]
                self.frq2=self.SFO3
                if(self.dim==4):
                    self.sw3=(float(self.GetParBruk('acqu4s',('','SW_h'))[0]))
                    self.SFO4=float(self.GetParBruk('acqu4s',('','SFO1'))[0]) #typically N
                    self.BF4=float(self.GetParBruk('acqu4s',('','BF1'))[0]) #typically N
                    self.O4=float(self.GetParBruk('acqu4s',('','O1'))[0]) #typically N
                    NUC4=self.GetParBruk('acqu4s',('','NUC1'))[0].split('<')[1].split('>')[0]
                    self.NUC4=NUC4[-1:]+NUC4[:-1]
                    self.frq3=self.SFO4        
            else: #deal with Marius Clore's 4D
                self.sw1=1/(2*float(self.GetParBruk('acqus',('','IN0'))[0]))
                self.frq1,self.f1ppm=self.GetShiftBruk(self.labb[1])
                if(self.dim>=3):
                    self.sw2=1/(2*float(self.GetParBruk('acqus',('','IN8'))[0]))
                    self.frq2,self.f2ppm=self.GetShiftBruk(self.labb[2])
                    if(self.labb[2][0]==self.labb[0][0]): #if the labelled nucleus for Z and X are the same...
                        self.SFO3=float(self.GetParBruk('acqus',('','SFO1'))[0]) #typically N
                        self.BF3=float(self.GetParBruk('acqus',('','BF1'))[0]) #typically N
                        self.O3=float(self.GetParBruk('acqus',('','O1'))[0]) #typically N
                    elif(self.labb[2][0]==self.labb[1][0]):  #if the labelled nucleus for Z and Y are the same...
                        self.SFO3=float(self.GetParBruk('acqu2s',('','SFO1'))[0]) #typically N
                        self.BF3=float(self.GetParBruk('acqu2s',('','BF1'))[0]) #typically N
                        self.O3=float(self.GetParBruk('acqu2s',('','O1'))[0]) #typically N
                        
                    if(self.labb[2][0]=='C'):
                        self.NUC3='C13'
                    elif(self.labb[2][0]=='H'):
                        self.NUC3='H1'
                    elif(self.labb[2][0]=='N'):
                        self.NUC3='N15'
                    if(self.dim==4):
                        self.sw3=1/(2*float(self.GetParBruk('acqus',('','IN19'))[0]))
                        self.frq3,self.f3ppm=self.GetShiftBruk(self.labb[3])

                        if(self.labb[3][0]==self.labb[0][0]):
                            self.SFO4=float(self.GetParBruk('acqus',('','SFO1'))[0]) #typically N
                            self.BF4=float(self.GetParBruk('acqus',('','BF1'))[0]) #typically N
                            self.O4=float(self.GetParBruk('acqus',('','O1'))[0]) #typically N
                        elif(self.labb[3][0]==self.labb[1][0]):
                            self.SFO4=float(self.GetParBruk('acqu2s',('','SFO1'))[0]) #typically N
                            self.BF4=float(self.GetParBruk('acqu2s',('','BF1'))[0]) #typically N
                            self.O4=float(self.GetParBruk('acqu2s',('','O1'))[0]) #typically N

                        if(self.labb[3][0]=='C'):
                            self.NUC4='C13'
                        elif(self.labb[3][0]=='H'):
                            self.NUC4='H1'
                        elif(self.labb[3][0]=='N'):
                            self.NUC4='N15'
                    
        if(self.dim==3):
            print self.NUC1,self.NUC2,self.NUC3
        
        #print self.SFO1,self.SFO2,self.SFO3
        #print self.O1,self.O2,self.O3
        
        #O=self.sfrq,self.frq1,self.frq2,self.frq3
        #C=self.waterppm,self.f1ppm,self.f2ppm,self.f3ppm

        self.sfrq,self.waterppm=self.GetShiftBruk(self.labb[0])

        if(self.o1p=='temp'):
            try:
                self.GetTempBruk()
                self.waterppm=self.WaterPPM()
            except:
                pass

            if(self.dim>1):
                self.f1ppm=self.O2/self.BF2        

        elif(self.o1p=='guess'):
            print self.labb
            #print self.NUC1,self.NUC2,self.NUC3
            #self.waterppm=self.o1p
            if(self.dim>1):
                self.f1ppm=self.O2/self.BF2
                #self.f1ppm=self.shift(self.frq1,nuc=self.NUC2)
                #if(self.dim>=3):
                #    self.f2ppm=self.shift(self.frq2,nuc=self.NUC3)
                #    print '2',self.frq1,self.NUC2,self.f1ppm
                #    print '3',self.frq2,self.NUC3,self.f2ppm
                #    if(self.dim==4):
                #        self.f3ppm=self.shift(self.frq3,nuc=self.NUC4)
        else:

            print self.labb
            #print self.NUC1,self.NUC2,self.NUC3
            self.waterppm=self.o1p
            if(self.dim>1):
                self.f1ppm=self.O2/self.BF2
        """
        else:
            
            try:
                self.GetTempBruk()
                self.waterppm=self.WaterPPM()
            except:
                pass

            if(self.dim>1):
                self.f1ppm=self.O2/self.BF2
            if(self.dim>=3):
                self.f2ppm=self.O3/self.BF3
                if(self.dim==4): #flemming's data implies we can't trust o1 in acqu4s
                    if(self.NUC4==self.NUC1):
                        self.f3ppm=self.waterppm
                        self.frq3=self.sfrq
                    elif(self.NUC4==self.NUC2):
                        self.f3ppm=self.f1ppm
                        self.frq3=self.frq1
                    elif(self.NUC4==self.NUC3):
                        self.f3ppm=self.f2ppm
                        self.frq3=self.frq2
                    else:
                        self.f3ppm=self.O4/self.BF4


            #print self.O1,self.BF1,self.SFO1
            #print self.O4,self.BF4,self.SFO4
                                    
            #print self.waterppm
            #sys.exit(100)
        """


            
        #self.temp=float(self.GetParBruk('acqus',('','TE'))[0])
        #print self.temp
        #try:
        #    self.waterppm=self.WaterPPM()
        #    print 'temp:',self.temp,'water ppm:',self.waterppm
        #except:
        #    print 'cannot find temperature: guessing carrier as 4.7ppm'
        #    self.waterppm=4.7
        #self.waterppm=float(self.GetParBruk('acqus',('','CNST20'))[0]) #typically N

        """
        #if(self.tn!='H1'):
        #    print 'direct nucleus isnt proton: you will need to write your own conversion script.'
        #    self.abort=1
        #    return 
        print self.labb,self.labb[0]
        self.frq1,self.n1=self.GetNuc(self.labb[1])
        self.frq2,self.n2=self.GetNuc(self.labb[2])

        
        if(self.frq1==-1 or self.frq2==-1):
            return

        if(self.dim==4):
            self.frq3,self.n3=self.GetNuc(self.labb[3])
            if(self.frq3==-1):
                return

        self.f1ppm=self.shift(self.frq1,self.n1)
        self.f2ppm=self.shift(self.frq2,self.n2)
        if(self.dim==4):
            self.f3ppm=self.shift(self.frq3,self.n3)

        #npAdj=BrukFidAdjust(np)
        """


    def GetShiftBruk(self,lab):
        if(lab[0]==self.NUC1[0]):
            return self.SFO1,self.O1/self.BF1
        elif(lab[0]==self.NUC2[0]):
            return self.SFO2,self.O2/self.BF2
        elif(lab[0]==self.NUC3[0]):
            return self.SFO3,self.O3/self.BF3
        else:
            print 'Cannot figure out shifts'
            return -1,-1
                






    def ConvertVarian(self):
        self.ni=int(self.GetParVarian(('','ni'))[0])
        self.np=int(self.GetParVarian(('','np'))[0])


        #nz=len(vpar.GetParVarian('./procpar','n',('','ncyc_cp')))

        if(self.dim==3):
            self.ni2=int(self.GetParVarian(('','ni2'))[0])
            print 'ni:',self.ni,'np:',self.np,'ni2:',self.ni2
        if(self.dim==4):
            self.ni2=int(self.GetParVarian(('','ni2'))[0])
            self.ni3=int(self.GetParVarian(('','ni3'))[0])
            print 'ni:',self.ni,'np:',self.np,'ni2:',self.ni2,'ni3:',self.ni3

        print self.GetParVarian(('','array'))[0]
        array=self.GetParVarian(('','array'))[0].split('"')[1].split(',')
        if(self.dim==2):
            self.acqORD=''

        if(self.dim==3):
            if(array[0]=='phase'): #phase,phase2
                self.acqORD=1
            else: #phase2,phase
                self.acqORD=0
        elif(self.dim==4):
            if(array[0]=='phase' and array[1]=='phase2'): #phase,phase2,phase3
                self.acqORD=1
            elif(array[0]=='phase3' and array[1]=='phase2'): #phase3,phase2,phase
                self.acqORD=0  

        #-aqORD  aqOrd  [0] Acquisition Order Code:
        #0 = 2D d2,phase
        #0 = 3D d3,d2,phase2,phase
        #0 = 4D d4,d3,d2,phase3,phase2,phase
        #1 = 3D d3,d2,phase,phase2
        #1 = 4D d4,d3,d2,phase,phase2,phase3
        #2 = 4D d3,d2,d4,phase3,phase2,phase  #NOT COVERED




        self.sw=float(self.GetParVarian(('','sw'))[0])
        self.sw1=float(self.GetParVarian(('','sw1'))[0])
        self.sw2=float(self.GetParVarian(('','sw2'))[0])
        if(self.dim==4):
            self.sw3=float(self.GetParVarian(('','sw3'))[0])
        
        self.sfrq=float(self.GetParVarian(('','sfrq'))[0])
        self.dfrq=float(self.GetParVarian(('','dfrq'))[0])
        self.dfrq2=float(self.GetParVarian(('','dfrq2'))[0])
        self.f1180_flg=self.GetParVarian(('','f1180'))[0]
        self.f2180_flg=self.GetParVarian(('','f2180'))[0]

        print self.labb

        if(self.o1p=='auto'):
            self.temp=float(self.GetParVarian(('','temp'))[0])
            print self.temp
            try:
                self.waterppm=self.WaterPPM()
                self.waterppmTOF=self.waterppm
                print 'temp:',self.temp,'water ppm:',self.waterppm
            except:
                print 'cannot find temperature: guessing carrier as 4.7ppm'
                self.waterppm=4.7
                self.waterppmTOF=self.waterppm
            try:
                print 'Found tof_me: moving carrier'
                self.tof=float(self.GetParVarian(('','tof'))[0])
                self.tof_me=float(self.GetParVarian(('','tof_me'))[0])
                self.waterppm+=-(self.tof-self.tof_me)/self.sfrq#adjust carrier, but not used for referencing other dimensions
            except:
                pass
        else:
            self.waterppm=self.o1p
            self.waterppmTOF=self.waterppm
            
        print self.waterppm
        print 'first nuclei (assuming proton):',self.labb[0],self.sfrq
        print 'second nuclei :',self.labb[1]

        self.tn=self.GetParVarian(('','tn'))[0].split('"')[1]
        self.dn=self.GetParVarian(('','dn'))[0].split('"')[1]
        self.dn2=self.GetParVarian(('','dn2'))[0].split('"')[1]
        print self.tn,self.dn,self.dn2

        if(self.tn!='H1'):
            print 'direct nucleus isnt proton: you will need to write your own conversion script.'
            self.abort=1
            return 
        print self.labb,self.labb[0]
        self.frq1,self.n1=self.GetNuc(self.labb[1])
        if(self.dim==3):
            self.frq2,self.n2=self.GetNuc(self.labb[2])

        
        if(self.frq1==-1):
            return
        if(self.dim>=3):
            if(self.frq2==-1):
                return
        
        if(self.dim==4):
            self.frq3,self.n3=self.GetNuc(self.labb[3])
            if(self.frq3==-1):
                return

        self.f1ppm=self.shift(self.frq1,self.n1)
        if(self.dim>=3):
            self.f2ppm=self.shift(self.frq2,self.n2)
        if(self.dim==4):
            self.f3ppm=self.shift(self.frq3,self.n3)


    def GetNUSsamp(self):
        inny=open(self.outdir+'/'+self.nuslist)
        samp=[]
        for line in inny.readlines():
            test=line.split()
            if(len(test)>0):
                row=[]
                for i in range(len(test)):
                    row.append(int(test[i]))
                samp.append(row)
        samp=numpy.array(samp)

        maxy=numpy.max(samp,axis=0)
        print 'maxy',maxy
        print 'Points in sampling schedule:',samp
        self.ni=maxy[0]+1
        self.ni2=maxy[1]+1
        if(self.dim==4):
            self.ni3=maxy[2]+1

        return len(samp) #return number of entries in sampling schedule        

    def PipeParse(self):    

        spa=10               #spacing in output file
        axis='x','y','z','a' #names of axes

        outy=open(self.outdir+'/fid.test.com','w')
        outy.write('#!/bin/csh\n')

        if(self.GZ=='ser'):
            outy.write('gunzip ser.gz\n')
        elif(self.GZ=='fid'):
            outy.write('gunzip fid.gz\n')
        else:
            pass

        if(self.nuslist!=''):
            samp=self.GetNUSsamp()
            if(self.tp=='var'):

                outy.write('nusExpand.tcl -mode varian -sampleCount %i -off 0 \\\n' % (samp))
                outy.write('-in %s/fid -out %s/fid_full -sample %s/%s -procpar %s/procpar\n' % (self.outdir,self.outdir,self.outdir,self.nuslist,self.outdir))
            elif(self.tp=='bruk'):
                #
                outy.write('cp %s/acqus ./\n' % self.outdir) #dodgy script. needs acqus
                outy.write('nusExpand.tcl -mode bruker -sampleCount %i -off 0 \\\n' % (samp))
                outy.write('-in %s/ser -out %s/ser_full -sample %s/%s\n' % (self.outdir,self.outdir,self.outdir,self.nuslist))
                outy.write('rm ./acqus\n')

        if(self.dim==2):
            outy.write('%s\n' % ('set ft4trec='+self.outdir+'/test.fid'))
            outy.write('if( -e %s) rm -rf %s\n' %(self.outdir+'/test.fid',self.outdir+'/test.fid'))
        if(self.dim==3):
            outy.write('%s\n' % ('set ft4trec='+self.outdir+'/fids/test%03d.fid'))
            outy.write('if( -e %s) rm -rf %s\n' %(self.outdir+'/fids/test001.fid',self.outdir+'/fids'))
        if(self.dim==4):
            outy.write('%s\n' % ('set ft4trec='+self.outdir+'/fids/test%03d%03d.fid'))
            outy.write('if( -e %s) rm -rf %s\n' %(self.outdir+'/fids/test001001.fid',self.outdir+'/fids'))

        if(self.tp=='omega'):
            pass
            #outy.write('bin2pipe -in %s -ge -neg \\\n' % (infile))
        elif(self.tp=='bruk'):
            if('ser' in os.listdir(self.outdir) or 'ser.gz' in os.listdir(self.outdir)):
                infile='ser'
            else:
                infile='fid'
            if(self.nuslist!=''):
                infile='ser_full'
            outy.write('bruk2pipe -in %s/%s  \\\n' % (self.outdir,infile))
            GRPDLY=float(self.GetParBruk('acqus',('','GRPDLY'))[0])
            DSPFVS=int(self.GetParBruk('acqus',('','DSPFVS'))[0])
            DECIM=float(self.GetParBruk('acqus',('','DECIM'))[0])

            BYTORDA=int(self.GetParBruk('acqus',('','BYTORDA'))[0])
            if(BYTORDA==1):
                flag='noaswap'
            else:
                flag='aswap'

            if(self.amx==True):
                MX='AMX'
            else:
                MX='DMX'

            if(self.ws==4):
                outy.write('-bad 0.0 -%s -%s -decim %f -dspfvs %i -grpdly %f \\\n' % (flag,MX,DECIM,DSPFVS,GRPDLY))
            else:
                outy.write('-bad 0.0 -ext -%s -%s -decim %i -dspfvs %i -grpdly %i -ws 8 -noi2f \\\n' % (flag,MX,DECIM,DSPFVS,GRPDLY))
        elif(self.tp=='var'):
            infile='fid'
            if(self.nuslist!=''):
                infile='fid_full'
            outy.write('var2pipe -in %s/%s \\\n' % (self.outdir,infile))

            outy.write(' -noaswap ')
            if(self.acqORD!=''):
                outy.write(' -aqORD %i \\\n' % (self.acqORD))
            outy.write('\\\n')

        M=[]
        if(self.tp=='bruk'):
            M.append('DQD')
        else:
            M.append('Complex')

        self.modDict={}
        self.modDict['0']='Complex'
        self.modDict['1']='QF'
        self.modDict['2']='QSEQ'
        self.modDict['3']='TPPI'
        self.modDict['4']='States'
        self.modDict['5']='States-TPPI'
        self.modDict['6']='Echo-Antiecho'
        
        for i in range(len(self.rk)):
            if(self.quant==False):
                if(self.rk[i]==0):
                    if(self.tp=='bruk'):
                        print self.mode[i]
                        print self.modDict[self.mode[i]]
                        M.append(self.modDict[self.mode[i]])
                    else:
                        M.append('Complex')
                else:
                    M.append('Rance-Kay')
            else:
                M.append('Real')

        if(self.dim==1):
            N=self.np,
            T=self.np/2,
            #M='Complex','Complex','Complex'
            sw=self.sw,
            O=self.sfrq,
            C=self.waterppm,

        if(self.dim==2):
            if(self.quant):
                N=self.np,self.ni_TD
                T=self.np/2,self.ni_TD
            else:
                N=self.np,self.ni*2
                T=self.np/2,self.ni
            #M='Complex','Complex','Complex'
            sw=self.sw,self.sw1
            O=self.sfrq,self.frq1
            C=self.waterppm,self.f1ppm

        if(self.dim==3):
            N=self.np,self.ni*2,self.ni2*2
            T=self.np/2,self.ni,self.ni2
            #M='Complex','Complex','Complex'
            sw=self.sw,self.sw1,self.sw2
            O=self.sfrq,self.frq1,self.frq2
            C=self.waterppm,self.f1ppm,self.f2ppm

        elif(self.dim==4):
            N=self.np,self.ni*2,self.ni2*2,self.ni3*2
            T=self.np/2,self.ni,self.ni2,self.ni3
            #M='Complex','Complex','Complex','Complex'
            sw=self.sw,self.sw1,self.sw2,self.sw3
            O=self.sfrq,self.frq1,self.frq2,self.frq3
            C=self.waterppm,self.f1ppm,self.f2ppm,self.f3ppm

        self.AddPipe(outy,axis,'N',N,spa)
        self.AddPipe(outy,axis,'T',T,spa)
        self.AddPipe(outy,axis,'MODE',M,spa)
        self.AddPipe(outy,axis,'SW',sw,spa)
        self.AddPipe(outy,axis,'OBS',O,spa)
        self.AddPipe(outy,axis,'CAR',C,spa)
        self.AddPipe(outy,axis,'LAB',self.labb,spa)

        outy.write(' -ndim  %s -aq2D  %s \\\n' % (str(self.dim).ljust(spa),'States'.ljust(spa)))
        #outy.write('| nmrPipe -fn MULT -c 6.10352e-02 \\\n')
        #outy.write(' -ndim  %s  \\\n' % (str(self.dim).ljust(spa)))
        #outy.write(' -ndim  %s \\\n' % (str(self.dim).ljust(spa),))

        if(self.dim==1):
            outy.write('  -out test.fid -verb -ov\n')
        if(self.dim==2):
            outy.write('  -out $ft4trec -verb -ov\n')
        if(self.dim==3):
            outy.write('  -out $ft4trec -verb -ov\n')
        if(self.dim==4):
            outy.write('| pipe2xyz -x -out $ft4trec -verb -ov -to 0\n')

        if(self.nuslist!=''):
            if(self.tp=='var'):
                infile='fid_full'
            elif(self.tp=='bruk'):
                infile='ser_full'
            else:
                infile='a'
            outy.write('rm %s/%s' % (self.outdir,infile))

        if(self.GZ=='ser'):
            outy.write('gzip ser\n')
        elif(self.GZ=='fid'):
            outy.write('gzip fid\n')
        else:
            pass

        #outy.write('| nmrPipe -ov -verb -out test.fid\n') #spit into a giant fid

        return


    def AddPipeLine(self,outy,lab,par,val,spa):
        outy.write(' %s%s ' % ('-'+lab+par.ljust(5),str(val).ljust(spa))) 

    def EndPipeLine(self,outy):
        outy.write('\\\n')

    def AddPipe(self,outy,axis,par,vals,spa):
        for i in range(self.dim):
            self.AddPipeLine(outy,axis[i],par,vals[i],spa)
        self.EndPipeLine(outy)


    #return water chemical shift in range 0-100oC
    def WaterPPM(self):
        return 5.060 - 0.0122*self.temp + (2.11E-5)*self.temp**2.

    def GetSpectrometerType(self,path='./'):
        print path
        if(os.path.exists(path+'/acqu')==1 or os.path.exists(path+'/acqu2')==1):
            sys.stdout.write('Found acqu: we are Bruker!\n')
            self.tp='bruk'
        elif(os.path.exists(path+'/procpar')==1):
            sys.stdout.write('Found procpar: we are varian!\n')
            self.tp='var'
            self.parfile=path+'/procpar'
        else:
            sys.stdout.write('Neither bruker nor varian - guessing GE!\n')
            self.tp='omega'

    def GetAcqseq(self):
        inny=open(self.outdir+'/'+'pulseprogram')
        for line in inny.readlines():
            test=line.split()
            if(len(test)>1):
                if(test[0]=='aqseq'):
                    return test[1]
        print 'Cannot find acqseq in pulseprogram'
        return -1
            
    def GetSequence(self):
        if(self.tp=='var'):
            seqfil=self.GetParVarian(('','seqfil'))[0].split('"')[1]
        elif(self.tp=='bruk'):
            seqfil=self.GetParBruk('acqu',('','PULPROG',))[0].split('<')[1].split('>')[0]
        elif(self.tp=='omeg'):
            parfile=self.GetOmegaParFile()
            test=self.GetParOmega(parfile,'n',('','seq_source',))[0].split('/')
            seqfil=test[len(test)-1]
        print 'pulse sequence name:',seqfil        

    def readfile(self,infile):
        peak=[]
        peakfile=open(infile,'r')
        for line in peakfile.readlines():
            linetosave=line.split()
            peak.append(linetosave)
        peakfile.close()
        return peak


    #analyse either acqu and acqu2
    def GetParBruk(self,infile,argv,verb='n',):
        verb='y'
        args=[]
        procpar=self.readfile(self.outdir+'/'+infile)
        for i in range(len(argv)-1):
            param=argv[i+1]
            tick=0
            for j in range(len(procpar)):
                #print procpar[j]
                test=procpar[j][0].split('##$')
                
                if(len(test)>1):
                    test2=test[1].split('=')[0]
                    if(test2==param):
                        if(verb=='y'):
                            sys.stdout.write('%s: %s\n' % (param,procpar[j][1]))
                        args.append(procpar[j][1])
                        tick=1
                else:
                    #we have a line of zeros
                    #is the previous line what we're after?
                    test=procpar[j-1][0].split('##$')
                    if(len(test)>1):
                        test2=test[1].split('=')[0]
                        for i in range(100):
                            parT=test2+str(i)
                            if(parT==param):
                                if(len(param.split(test2))>1):
                                    if(verb=='y'):
                                        sys.stdout.write('Param %s found. Range: %s\n' % (param,procpar[j-1][1]))
                                    #parameters are in rows in j,j+1,j+2...
                                    go=0
                                    cnt=0
                                    while(go==0):
                                        if(i<len(procpar[j+cnt])):
                                            val=procpar[j+cnt][i]
                                            go=1
                                        else:
                                            i-=len(procpar[j+cnt])
                                            cnt+=1

                                    args.append(val)
                                    if(verb=='y'):
                                        sys.stdout.write('%s: %s\n' % (param,val))
                                    tick=1
                #sys.exit(100)
            if(tick==0):
                if(verb=='y'):
                    sys.stdout.write('Could not find param %s in %s\n' % (param,infile))
                return 'fail'
            else:
                return args

    def GetParVarian(self,argv,verb='n'):
        verb='y'
        args=[]
        procpar=self.readfile(self.parfile)
        for i in range(len(argv)-1):
            param=argv[i+1]
            tick=0
            for j in range(len(procpar)):
                if(procpar[j][0]==param):
                    if(verb=='y'):
                        sys.stdout.write('%s: %i argument' % (param,int(procpar[j+1][0])))
                    if(int(procpar[j+1][0])>1):
                        sys.stdout.write('s')
                    sys.stdout.write('\n')
                    for k in range(int(procpar[j+1][0])):
                        try:
                            if(verb=='y'):
                                sys.stdout.write('%s ' % (procpar[j+1][k+1]))
                            args.append(procpar[j+1][k+1])
                        except:
                            if(verb=='y'):
                                sys.stdout.write('%s ' % (procpar[j+1+k][0]))
                            args.append(procpar[j+1+k][0])
                    if(verb=='y'):
                        sys.stdout.write('\n')
                    tick=1
        if(tick==0):
            if(verb=='y'):
                sys.stdout.write('Could not find param %s in procpar\n' % (param))
            return 'fail'
        else:
            return args


    def GetParOmega(self,infile,verb,argv):
        args=[]
        procpar=self.readfile(infile)
        for i in range(len(argv)-1):
            param=argv[i+1]
            tick=0
            for j in range(len(procpar)):
                if(procpar[j][0]==param):
                    if(verb=='y'):
                        sys.stdout.write('%s: ' % (param))
                    if(len(procpar[j])>2):
                        for k in range(len(procpar[j])-2):
                            if(verb=='y'):
                                sys.stdout.write('%s ' % (procpar[j][k+1]))
                            args.append(procpar[j][k+1])
                    else:
                        for k in range(len(procpar[j])-1):
                            if(verb=='y'):
                                sys.stdout.write('%s ' % (procpar[j][k+1]))
                            args.append(procpar[j][k+1])


                    if(verb=='y'):
                        sys.stdout.write('\n')
                    tick=1
        if(tick==0):
            if(verb=='y'):
                sys.stdout.write('Could not find param %s in %s\n' % (param,infile))
            return 'fail'
        else:
            return args


    def GetOmegaVal(self,infile,param):
        if(self.GetParOmega(infile,'n',('',param))!='Fail'):
            return self.GetParOmega(infile,'n',('',param))


    def GetBrukVal(self,infile,param):
        test=self.GetParBruk(infile,'n',('',param))
        if(test!='Fail'):
            return test

    def GetParBrukFile(self,filey):
        inny=open(filey)
        pars=[]
        for line in inny.readlines():
            test=line.split()
            if(len(test)>0):
                pars.append(test)
        return pars
    
