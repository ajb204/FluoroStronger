#!/usr/bin/python
###########################
#A.Baldwin July 2011
#(c) University of Oxford

from __future__ import absolute_import
from __future__ import print_function
import sys,os,numpy,copy
from scipy.optimize import leastsq
from six.moves import range

class fitter():
    def __init__(self,infile,tag,rng,A,Aqs,Ass,offqs,offss,Rqs,Rss,Jqs,frq,fit=True,jiggles=True,g=0.5E11):
        self.infile=infile
        self.tag=tag
        
        self.frq=frq #564.680 

        self.readfile() #read raw data text file

        #copy variables into the class
        self.g=g
        self.A1=1.

        #initial fitting parameters
        self.A=A  #intensity of first quartet
        self.Aqs=Aqs #intensity of remaining quartets
        self.Ass=Ass #intensity of singles

        self.offPPMqs=offqs
        self.offPPMss=offss

        self.Rqs=Rqs
        self.Rss=Rss

        self.Jqs=Jqs
        
        self.rng=rng

        #jiggler settings.
        self.jiggles=20
        self.sigma=0.1

        
        ######## turn everything into numpy for quick calcs...
        self.offPPMqs=numpy.array(self.offPPMqs)*1.
        self.offPPMss=numpy.array(self.offPPMss)*1.
        self.Jqs=numpy.array(self.Jqs)*1.
        self.Ass=numpy.array(self.Ass)*1.
        self.Aqs=numpy.array(self.Aqs)*1.
        self.Rss=numpy.array(self.Rss)*1.
        self.Rqs=numpy.array(self.Rqs)*1.


        for i in range(len(self.offPPMqs)):
            self.offPPMqs[i]=numpy.sort(self.offPPMqs[i])

        ##########
        
        
        self.StrongCalc() #simulate specrtum.

        if(fit):  #if fitting, fit pars
            self.Fit()  
        if(jiggles):  #if wanting to do jiggles to get out of local minima...
            self.Jiggler()
        if(fit or jiggles):  #when done, repeat calculation
            self.StrongCalc() 
        self.printfile()  #dump output file
        print('Finished! Have a nice day.')

    ############################################
    def pack(self): #pack variables into fitting vector
        x=[]
        for i in range(len(self.offPPMqs)):
            x.append(self.offPPMqs[i,0])
            x.append(self.offPPMqs[i,1])
        for i in range(len(self.offPPMss)):
            x.append(self.offPPMss[i])


        x.append(numpy.fabs(self.A))
        for i in range(len(self.Aqs)):
            x.append(self.Aqs[i])
        for i in range(len(self.Ass)):
            x.append(self.Ass[i])
        for i in range(len(self.Rqs)):
            x.append(self.Rqs[i])
        for i in range(len(self.Rss)):
            x.append(self.Rss[i])
        for i in range(len(self.Jqs)):
            x.append(self.Jqs[i])
        return x
    def unpack(self,x): #unpack fitting vector to class 
        cnt=0
        for i in range(len(self.offPPMqs)):
            self.offPPMqs[i,0]=x[cnt];cnt+=1
            self.offPPMqs[i,1]=x[cnt];cnt+=1
        for i in range(len(self.offPPMss)):
            self.offPPMss[i]=x[cnt];cnt+=1

        
        self.A=numpy.fabs(x[cnt]);cnt+=1
        for i in range(len(self.Aqs)):
            self.Aqs[i]=numpy.fabs(x[cnt]);cnt+=1
        for i in range(len(self.Ass)):
            self.Ass[i]=numpy.fabs(x[cnt]);cnt+=1
        for i in range(len(self.Rqs)):
            self.Rqs[i]=numpy.fabs(x[cnt]);cnt+=1
        for i in range(len(self.Rss)):
            self.Rss[i]=numpy.fabs(x[cnt]);cnt+=1
        for i in range(len(self.Jqs)):
            self.Jqs[i]=numpy.fabs(x[cnt]);cnt+=1

    def chi(self,x):  #calculate chi
        self.unpack(x)
        self.StrongCalc()
        return self.raw-self.sim
    def Fit(self,verb='y'): #do the fit
        xinit=self.pack()
        x0=leastsq(self.chi,xinit)
        if(verb=='y'):
            print(x0[0])
    def chi2(self):  #return the chi2 value
        return numpy.sum((self.sim-self.raw)**2.)/len(self.ppm)/numpy.max(self.raw)
        
    def addnoise(self,pars):  #for jiggler, add noise to pars

        parsNew=pars*numpy.random.normal(1.,self.sigma,len(pars))

        #be gentle with ppms
        parsNew[:len(self.offPPMqs)*2+len(self.offPPMss)]=pars[:len(self.offPPMqs)*2+len(self.offPPMss)]*numpy.random.normal(1.,self.sigma/50,len(self.offPPMqs)*2+len(self.offPPMss))
        self.unpack(parsNew)

    

    def Jiggler(self): #run jiggler
        
        parsBest=self.pack()
        chi2best=self.chi2()
        print('starting chi2:',chi2best)
        go=0
        cnt=0
        while(go==0):
            #print
            #print parsBest
            #print 
            self.addnoise(parsBest)
            #print self.pack()
            self.Fit(verb='n')
            chi2curr=self.chi2()
            print('best:',chi2best,'curr:',chi2curr,'cnt:',cnt)


            if(chi2curr<chi2best):
                print('   yay! keeping!')
                parsBest=self.pack()
                chi2best=chi2curr
                cnt=0
            cnt+=1
            if(cnt==self.jiggles):
                go=1
        
        self.unpack(parsBest)
        self.Fit()
        

    def readfile(self): #read raw data
        print('Reading ',self.infile)
        dat=[]
        inny=open(self.infile)
        for line in inny.readlines():
            test=line.split()
            if(len(test)>1):
                dat.append((float(test[0]),float(test[1])))
        self.dat=numpy.array(dat)
        self.ppm=self.dat[:,0]
        self.raw=self.dat[:,1]
        self.sw_h=numpy.fabs(self.ppm[-1]-self.ppm[0])*self.frq
        self.hz=self.ppm*self.frq
        

    def printfile(self):  #make nice output
        outy=open('outy.out','w')
        for i in range(len(self.dat)):
            outy.write('%f\t%f\t%f\n' % (self.ppm[i],self.raw[i],self.sim[i]))
        outy.write('\n\n')

        self.A1tmp=self.A1
        self.Aqstmp=copy.deepcopy(self.Aqs)
        self.Asstmp=copy.deepcopy(self.Ass)

        
        for k in range(len(self.offPPMqs)):
            for i in range(len(self.Ass)):
                self.Ass[i]=0
            if(k==0):
                self.A1=1
                for i in range(len(self.Aqs)):
                    self.Aqs[i]=0
            else:
                self.A1=0
                for i in range(len(self.Aqs)):
                    if(i==k-1): #when k=1, i=0
                        self.Aqs[i]=self.Aqstmp[i]
                    else:
                        self.Aqs[i]=0
            self.StrongCalc()
            for i in range(len(self.dat)):
                outy.write('%f\t%f\t%f\n' % (self.ppm[i],self.raw[i],self.sim[i]))
            outy.write('\n\n')



        for k in range(len(self.offPPMss)):
            self.A1=0
            for i in range(len(self.Aqs)):
                self.Aqs[i]=0
            for i in range(len(self.Ass)):
                if(i==k): #when k=1, i=0
                    self.Ass[i]=self.Asstmp[i]
                else:
                    self.Ass[i]=0
            self.StrongCalc()
            for i in range(len(self.dat)):
                outy.write('%f\t%f\t%f\n' % (self.ppm[i],self.raw[i],self.sim[i]))
            outy.write('\n\n')

        outy.close()

        self.Aqs=self.Aqstmp
        self.Ass=self.Asstmp
        self.A1=1
        
        gnu=open('gnu.gp','w')  #make gnuplot script
        gnu.write('set term post eps enh color solid\n')
        gnu.write('set output \'fit.'+self.tag+'.eps\'\n')
        gnu.write('set title \''+self.tag+'\'\n')
        gnu.write('set xrange[%f:%f]\n' % (numpy.max(self.rng),numpy.min(self.rng)))

        gnu.write('unset ytics\n')
        #gnu.write('unset tics top\n')
        pos=0.35
        step=0.05
        gnu.write('set xlabel\'ppm\'\n')
        gnu.write('set format y ""\n')
        gnu.write('set border 1\n')

        for i in range(len(self.offPPMqs)):
            if(i==0):
                I=self.A1
            else:
                I=self.Aqs[i-1]
            gnu.write('set label sprintf("ppm: %.2f,%.2f R2: %.0f J: %.0f I: %.2f") at graph 0.02,%f font "Arial,10"\n' % (self.offPPMqs[i,0],self.offPPMqs[i,1],self.Rqs[i],self.Jqs[i],I,pos));pos-=step

        for i in range(len(self.offPPMss)):
            gnu.write('set label sprintf("ppm: %.2f R2: %.0f I: %.2f") at graph 0.02,%f font "Arial,10"\n' % (self.offPPMss[i],self.Rss[i],self.Ass[i],pos));pos-=step


        gnu.write('plot \'outy.out\' i 0 ti \'raw\' w li,\'\' i 0 u 1:3 ti \'fit\'w li')
        cnt=1
        for i in range(len(self.offPPMqs)):
            gnu.write(',\'\' i %i u 1:($3-%f) noti w li lc 3' % (cnt,self.g*cnt));cnt+=1
        for i in range(len(self.offPPMss)):
            gnu.write(',\'\' i %i u 1:($3-%f) noti w li lc 3' % (cnt,self.g*cnt));cnt+=1

        gnu.close()
        os.system('gnuplot gnu.gp')
        
        
    def Lorentz(self,x,R2,x0):  #calculate lorenztian lineshape
        return self.A*2.*R2/( (x-x0)**2.+R2**2.) *(self.sw_h)/ ( 5.*numpy.sqrt(2*numpy.pi))
    
    def StrongCalc(self): #run strong couping calculation

        self.offqs=self.offPPMqs*self.frq
        self.offss=self.offPPMss*self.frq

        Zcc=(self.Jqs**2.+((self.offqs[:,0]-self.offqs[:,1]))**2.)**0.5
        Ycc=(self.offqs[:,1]+self.offqs[:,0])    

        out=(1-numpy.sin( self.Jqs/(self.offqs[:,1]-self.offqs[:,0])))
        inn=(1+numpy.sin( self.Jqs/(self.offqs[:,1]-self.offqs[:,0])))


        
        self.sim=numpy.zeros_like(self.raw)
        for i in range(len(self.offqs)):
            if(i==0):
                I=self.A1
            else:
                I=self.Aqs[i-1]
            self.sim+=self.Lorentz(self.hz,self.Rqs[i]/numpy.pi/2., 0.5*(self.Jqs[i]+Ycc[i]-Zcc[i]))*I*inn[i]/4.#innerA
            self.sim+=self.Lorentz(self.hz,self.Rqs[i]/numpy.pi/2.,-0.5*(self.Jqs[i]-Ycc[i]-Zcc[i]))*I*inn[i]/4.#innerB
            self.sim+=self.Lorentz(self.hz,self.Rqs[i]/numpy.pi/2., 0.5*(self.Jqs[i]+Ycc[i]+Zcc[i]))*I*out[i]/4.#outerA
            self.sim+=self.Lorentz(self.hz,self.Rqs[i]/numpy.pi/2.,-0.5*(self.Jqs[i]-Ycc[i]+Zcc[i]))*I*out[i]/4.#outerB

        for i in range(len(self.offss)):
            self.sim+=self.Lorentz(self.hz,self.Rss[i]/numpy.pi/2.,self.offss[i])*self.Ass[i]

        #spectra were fit with user entered second order (Nq) and singlet contributing (Ns) species.
        #for each second order species, an intensity, scalar coupling constant, linewidth and two oscillatior frequencies need to be provided.
        #for each singlet, an intensity, a linewidth and a peak position need to be specified.
        #the program optimises the set of Nq*4+Ns*2 parameters to get the best fit.
        #the second order spectrum was calclated to be a sum of 4 lorentz functions of peak positions +/- 0.5(J +/- Y - Z) (inner) and +/- 0.5 (J+/- Y +Z ) (outer), scaled by an intensity
        #factor of A(1+ sin (J/(deltaO)))/4 for the inner lines and A(1-sin(J/deltaO)/4 for the outer lines. 
        
        #in/2 + out/2 = 1 for quartet to have same intensity as singlet
            
        #print numpy.max(self.spec)
        #print numpy.max(spec)

        #outy=open('simulations/out/strong.out','w')
        #for i in range(len(self.FrqDirect)):
        #    outy.write('%f\t%f\n' % (self.FrqDirect[i],spec[0,i].real))
        #outy.close()

