#!/usr/bin/python

from vpar_decon import vpar

import os,numpy,sys,re
import nmrglue as ng

#######################################################
# Core class for processing. Reads acqus to get
# key parameters and constructs fid.com and nmrproc.com
# processing scripts. Requires gnuplot and latex to be
# in the system's path
#
# A.Baldwin July 2011
# (c) University of Oxford

class process():
    def __init__(self,pathy,tag,start,finish,integ,phase,p1,wf,ws=4,O1='auto',base=False,freqFile='fq2list',amx=False,series=False,basePars=((9.5,10,-3.,-2),5,0.1)):

        self.O1=O1
        self.tag=tag
        self.integ=integ
        self.start=start
        self.finish=finish
        self.phase=phase
        self.p1=p1
        self.wf=wf
        self.ws=ws
        self.base=base
        self.basePars=basePars
        print basePars
        self.STD=False
        self.freqFile=freqFile
        self.amx=amx
        self.series=series

        self.before=os.getcwd() #save starting path
        os.chdir(self.before+'/'+pathy)

        self.Anal(self.start,self.finish,self.tag,self.integ)
        if(self.STD==True or self.STD=='2D'):
            self.Analyse()
        os.chdir(self.before)  #return to where we started

    def Anal(self,start,finish,tag,integ):
        files=os.listdir('./')

        #import os
        if(start=='find'):
            #go looking for paths with an fid in.
            rawfiles=[os.path.join(dp, f) for dp, dn, fn in os.walk(os.path.expanduser("./")) for f in fn]
            files=[] 
            for file in rawfiles:
                if('fid' in file and 'test.fid' not in file and 'fid.test.com' not in file):
                    #print file.split('/')[:-1]
                    #print os.path.join(file.split('/')[:-1])
                    #files.append(os.path.join(file.split('/')[:-1]))
                    files.append(file.replace('fid',''))
                    
            #print files
            #sys.exit(100)
        else:
            files=os.listdir('./')
            
        cnt=1
        outy2=open('integration.'+tag+'.out','w')
        death=[]
        outy=open('raw.'+tag+'.data','w');outy.close()
        self.STD=False
        for file in files:
            #print file,files
            if(start!='find'):
                try:
                    dirno=int(file)
                except:
                    continue
                if('.' in file or 'figs' in file or int(file)<start or int(file)>=finish):
                    continue
            
            if(1==1):
                bef=os.getcwd()
                print file
                #if(file=='./Round 1/Sample 2 - CF2 Met(O)/7 h 19F + f_b/7 h 19F main/'):
                #    pass
                #else:
                #    continue

                    #sys.exit(100)
                try:
                    os.chdir(bef+'/'+file)
                except:
                    continue


                if(os.path.exists('./pdata')):
                    os.system('rm -rf pdata')

                print
                print
                print 'Looking at expNo:',file
                
                folly=os.listdir('./')
                print folly

                if('ser' in folly or 'ser.gz' in folly): #we are a 2D or pseudo 2D


                    ##CB
                    #cleanup excite raw files.
                    os.chdir('../')
                    remove_string= "ls | grep -P \'raw."+tag+"\d.data$\' | xargs -d\'\\n\' rm"
                    os.system(remove_string)
                    os.chdir(file)
                    ##CB


                    inst=vpar('.',2,'H1',(0,0,),self.phase,self.p1,self.wf,quant=True,ws=self.ws,o1p=self.O1,base=self.base,amx=self.amx)
                    print 'file:',file
                    d20=float(inst.GetParBruk('acqus',('','D20'))[0])
                    pl10=float(inst.GetParBruk('acqus',('','CNST62'))[0])
                    pars=(inst.GetParBrukFile(self.freqFile))
                    print self.freqFile,':',pars
                    excite=[]
                    for par in pars:
                        try:
                            excite.append( (float(par[0])-inst.O1)/inst.frq+inst.waterppm)
                        except:
                            pass
                    print 'excitation ppms:',excite
                    outy=open('../raw.'+tag+'.data.excite','w') #print excitation ppm values (for use later)
                    for e in excite:
                        outy.write('%f\n' % e)
                    outy.close()
                    excite=numpy.array(excite) #turn excitation frequencies into numpy array
                    argyExMax=numpy.argmax(numpy.fabs(excite))                           #find maximum number in exciation (off resonance)
                    
                    dic,data = ng.pipe.read('test.ft2') #read fids
                    Size=data.shape
                    #now lets read in the data and make some plots.
                    uc0 = ng.pipe.make_uc(dic,data,dim=0)
                    uc1 = ng.pipe.make_uc(dic,data,dim=1)
                    index=[]
                    for i in range((Size[1])):
                        index.append((uc1.ppm(0)-i*(-uc1.ppm(Size[1]-1)+uc1.ppm(0))/(Size[1]-1)))
                    index=numpy.array(index) #ppms


                    if('vdlist' in folly):
                        #read in mixing times, and consolidate data to unique mixing times
                        pars2=(inst.GetParBrukFile('vdlist'))
                        print 'vdlist',':',pars2 
                        data=data.reshape((len(pars2),len(pars),Size[1])) #reshape the data
                        frq=numpy.zeros_like(data)
                        mix=numpy.zeros_like(data)
                        mixTotal=[]
                        for i,p in enumerate(pars2): #go through vdlist, get mixing times
                            val=float(re.findall(r"\d+\.?\d*",p[0])[0]) #get numbers from file
                            if p[0][-1:] == 'm':  #if there's an 'm', divide by 1000
                                val = val/1000.
                            mix[i,:,:]=val  #setup 3D array with mixing times
                            mixTotal.append(val) #append time
                        for i,p in enumerate(pars):
                            val=excite[i]
                            frq[:,i,:]=val  #setup 3D array with frequencies
                            print val
                        mixTotal=numpy.array(mixTotal) #turn to numpy
                        mix=numpy.unique(mix)  #get unique mixing times.
                        print mix       #unique mixing times
                        print mixTotal  #tota mixing times
                        datNew=numpy.zeros((len(mix),len(pars),Size[1])) #get new array for data
                        for m,mi in enumerate(mix):  #for each unique mixing time
                            for e,ex in enumerate(excite):  #for each unique excitation time
                                mask=(numpy.absolute(mixTotal-mi)<0.000001)  #get indicies for mixing times that are aligned with unique mixing times
                                datNew[m,e,:]=numpy.sum(data[mask,e,:],axis=0)  #sum mixing times to get unique mixing time array only.
                        data=datNew      #copy datNew into data
                        Size=data.shape  #get new shape

                        #now lets read in the data and make some plots.
                        for k in range(Size[1]): #frq
                            if(k!=argyExMax): #don't print out the argmax.
                                if(Size[1]!=2): #if multiple excitations
                                    outy=open('../raw.'+tag+str(k)+'.data','w')
                                    self.STD='2D'
                                else: #if one excitation
                                    outy=open('../raw.'+tag+'.data','w')
                                    self.STD=True
                                #print 'raw.'+tag+str(k)+'.data'  #filename indexed by excite tag.
                                for j in range(Size[0]): #mix
                                    for i in range((Size[2])): #ppm              #file, mixing time, power, ppm, pulseon,pulseoff
                                        outy.write('%s\t%f\t%f\t%f\t%e\t%e\n' % (file,mix[j],pl10,index[i],data[j,k,i],data[j,argyExMax,i]))
                                    outy.write('\n\n')
                                outy.close()
                        #do some pointless integration
                        vals=[]  #do some integration (no longer really neccessary)
                        for inty in integ:
                            argy=numpy.abs(index-inty).argmin()
                            #print argy,index[argy],inty
                            #vals.append(data[1,argy]-data[0,argy])
                            vals.append(data[-1,-1,argy]-data[0,-1,argy])
                        print 'vals:',vals


                    else:
                        for k in range(Size[0]): #frq
                            if(k!=argyExMax): #don't print out the argmax.
                                if(Size[0]!=2):
                                    fil='../raw.'+tag+str(k)+'.data'
                                    self.STD="2D"
                                else: #original structure
                                    fil='../raw.'+tag+'.data'
                                    self.STD=True
                                if os.path.exists(fil)==0:
                                    red='w'
                                else:
                                    red='a'
                                outy=open(fil,red)
                                for i in range((Size[1])): #ppm              #file, mixing time, power, ppm, pulseon,pulseoff
                                    outy.write('%s\t%f\t%f\t%f\t%e\t%e\n' % (file,d20,pl10,index[i],data[k,i],data[argyExMax,i]))
                                outy.write('\n\n')
                                outy.close()

                        vals=[]
                        for inty in integ:
                            argy=numpy.abs(index-inty).argmin()
                            #print argy,index[argy],inty
                            vals.append(data[1,argy]-data[0,argy])
                        print 'vals:',vals


                    print 'index:',index.shape
                    print "Spectrum dimensions (pts): ",Size   #print the spectral dimensions
                    print "direct dimension limits (ppm): ", numpy.max(index),numpy.min(index)  #direct 


                    outy2.write('%f\t%f\t%e\t%e\n' % (d20,pl10,vals[0],vals[1]))
                    os.system('rm test.fid') #cleanup
                    os.system('rm test.ft2') #cleanup


                else:
                    print 'hello! I am a 1D NMR spectrum'

                    #inst=vpar('.',2,'H1',(0,0,),self.phase,self.wf,quant=True)                
                    #if('test.ft2' not in os.listdir('./')):
                    if(1==1):
                        inst=vpar('.',1,'H1',(0,),self.phase,self.p1,self.wf,quant=True,ws=self.ws,o1p=self.O1,base=self.base,amx=self.amx,basePars=self.basePars)
                        dic,data = ng.pipe.read('test.ft2') #read fids
                        Size=data.shape
                        print "Spectrum dimensions (pts): ",Size   #print the spectral dimensions        
                        index=[]
                        #for i in range((Size[0])):
                        #    index.append((uc0.ppm(0)-i*(-uc0.ppm(Size[0]-1)+uc0.ppm(0))/(Size[0]-1)))
                        #    outy.write('%s\t%f\t%f\t%f\t%e\n' % (file,d20,pl10,index[i],data[i]))
                        
                    #print file
                    if(self.series):
                        death.append(int(file))
                    else:
                        death.append(file)
                print os.getcwd()
                os.chdir(bef)
                print bef
                cnt+=1

        print death
        outy2.close()

        if(os.path.exists('figs')==0):
            os.system('mkdir figs')
        else:
            os.system('rm figs/*')

        if(self.STD==True):
            gnu=open('gnu.gp','w')
            gnu.write('set term post eps enh color solid\n')
            gnu.write('set output \'figs/fig3.'+tag+'.eps\'\n')
            gnu.write('set size square\n')
            gnu.write('unset key\n')
            gnu.write('set ticslevel 0\n')
            gnu.write('set view 82,255\n')
            gnu.write('set title \'%s difference excite: %.2f, %.2f ppm\'\n' % (tag,excite[0],excite[1]))
            gnu.write('set ylabel \'ppm\'\n')
            gnu.write('set xlabel \'delay(s)\'\n')
            gnu.write('set yrange[5:0]\n')
            gnu.write('splot \'raw.'+tag+'.data\' u 2:4:($6-$5):2 w li lc palette\n')

            gnu.write('set output \'figs/fig1.'+tag+'.eps\'\n')
            gnu.write('set title \'%s raw excite: %.2f ppm\'\n' % (tag,excite[1]))
            gnu.write('splot \'raw.'+tag+'.data\' u 2:4:($6):2 w li lc palette\n')
        
            gnu.write('set output \'figs/fig2.'+tag+'.eps\'\n')
            gnu.write('set title \'%s raw excite: %.2f ppm\'\n' % (tag,excite[0]))
            gnu.write('splot \'raw.'+tag+'.data\' u 2:4:($5):2 w li lc palette\n')
            
            gnu.close()
            os.system('gnuplot gnu.gp')
            print 'deleting:',death
            for de in death:
                os.system('rm -rf '+de)
        elif(self.STD=='2D'):
            gnu=open('gnu.gp','w')
            gnu.write('set term post eps enh color solid\n')

            gnu.write('set size square\n')
            gnu.write('unset key\n')
            gnu.write('set ticslevel 0\n')
            gnu.write('set view 82,255\n')
            for k in range(len(excite)):
                if(excite[k]!=excite[argyExMax]):
                    if(len(excite)==2):
                        filtag=tag
                    else:
                        filtag=tag+str(k)

                    gnu.write('set output \'figs/fig3.'+filtag+'.eps\'\n')
                    gnu.write('set title \'%s difference excite: %.2f, %.2f ppm\'\n' % (tag,excite[-1],excite[k]))
                    gnu.write('set ylabel \'ppm\'\n')
                    gnu.write('set xlabel \'delay(s)\'\n')
                    gnu.write('set yrange[5:0]\n')
                    gnu.write('splot \'raw.'+filtag+'.data\' u 2:4:($6-$5):2 w li lc palette\n')

                    gnu.write('set output \'figs/fig1.'+tag+'.eps\'\n')
                    gnu.write('set title \'%s raw excite: %.2f ppm\'\n' % (tag,excite[-1]))
                    gnu.write('splot \'raw.'+filtag+'.data\' u 2:4:($6):2 w li lc palette\n')
                    
                    gnu.write('set output \'figs/fig2.'+tag+'.eps\'\n')
                    gnu.write('set title \'%s raw excite: %.2f ppm\'\n' % (tag,excite[k]))
                    gnu.write('splot \'raw.'+filtag+'.data\' u 2:4:($5):2 w li lc palette\n')
            
            gnu.close()
            os.system('gnuplot gnu.gp')
            print 'deleting:',death
            for de in death:
                os.system('rm -rf '+de)


        elif(self.series):

            death=numpy.array(death)
            death=numpy.sort(death)

            Dic,Data = ng.pipe.read(str(death[0])+'/test.ft2') #read fids
            #adjust to make the correct size with headers
            #Dat,Dita = ng.pipe.read(str(death[0])+'/test.fid') #read fids
            
            NewData=numpy.zeros((len(death),len(Data)),dtype='float32')
            print NewData.shape

            #fidData=numpy.zeros((len(Dita)),dtype='complex64')
            
            for i,de in enumerate(death):
                dic,data = ng.pipe.read(str(de)+'/test.ft2') #read fids
                print de,data.shape
                NewData[i,:]=data

                #dac,dita = ng.pipe.read(str(de)+'/test.fid') #read fids
                #fidData+=dita


            #ng.pipe.write('test.fid',dac,fidData,overwrite=True)
            
            #write
            #doc,dota=ng.pipe.read('/Users/andrewbaldwin/Desktop/OldDesktop/mab/study/ft/076-8543-007.ft2')
            #print dota.shape

            #convert to 2D
            Dic['FDDIMCOUNT']=2
            Dic['FDF1TDSIZE']= len(death)
            Dic['FDF1FTSIZE']= len(death)
            Dic['FDSLICECOUNT']=len(death)
            Dic['FDREALSIZE']=len(death)
            Dic['FDSPECNUM']=len(death)

            Dic['FDF1CAR']=len(death)/2.
            Dic['FDF1CENTER']=len(death)/2.
            Dic['FDF1FTFLAG']=1
            Dic['FDF1OBS']=1
            Dic['FDF1SW']=60.

            #need carrier: sweep width, obs , min and max
            #sw      car    obs   min(ppm)   max(ppm)   
            #1 Y 60 0 nan 0 nan nan

            #set xrange for dim12.

            print os.getcwd()
            #for key,val in Dic.items():
            #    if(val!=doc[key]):
            #        print key,val,doc[key]


            ng.pipe.write('test.ft2',Dic,NewData,overwrite=True) #read fids
            #read back in to check

            Dic,Data = ng.pipe.read('test.ft2') #read fids
            Size=Data.shape
        
            uc0 = ng.pipe.make_uc(Dic,Data,dim=0)
            uc1 = ng.pipe.make_uc(Dic,Data,dim=1)
            index1=[]
            index0=[]
            for i in range((Size[1])):
                index1.append((uc1.ppm(0)-i*(-uc1.ppm(Size[1]-1)+uc1.ppm(0))/(Size[1]-1)))
            for i in range((Size[0])):
                index0.append((uc0.ppm(0)-i*(-uc0.ppm(Size[0]-1)+uc0.ppm(0))/(Size[0]-1)))
                

            outy=open('test3D.txt','w')
            for j in range(Size[0]):
                for i in range(Size[1]): #ppms.
                    outy.write('%f\t%f\t%e\n' % (index1[i],index0[j],Data[j,i]))
                outy.write('\n')
            outy.close()

            outy=open('testP1.txt','w')
            for i in range(Size[1]): #ppms.
                outy.write('%f\t%e\n' % (index1[i],numpy.sum(Data[:,i])))
            outy.close()

            outy=open('testP2.txt','w')
            for j in range(Size[0]): #ppms.
                outy.write('%f\t%e\n' % (index0[j],numpy.sum(Data[j,:])))
            outy.close()

            gnu=open('figs/gnu.gp','w')
            gnu.write('set term post eps enh color solid\n')


            #plot ppm projection.
            gnu.write('set output \'figs/projPpm.eps\'\n')
            gnu.write('unset key\n')
            gnu.write('set xlabel \'ppm\'\n')
            gnu.write('set ylabel \'\'\n')
            gnu.write('set xrange [%f:%f]\n' % (numpy.max(index1),numpy.min(index1)))
            gnu.write('plot \'testP1.txt\' u 1:2 w li\n')

            #plot fraction projection
            gnu.write('set output \'figs/projFrac.eps\'\n')
            gnu.write('set  xlabel \'fraction\'\n')
            gnu.write('set ylabel \'\'\n')
            gnu.write('set xrange [%f:%f]\n' % (numpy.min(index0),numpy.max(index0)))
            gnu.write('plot \'testP2.txt\' u 1:2 w li\n')

            """
            #plot 3D
            gnu.write('set output \'figs/topDown.eps\'\n')
            gnu.write('set pm3d map\n')
            gnu.write('set xlabel \'ppm\'\n')
            gnu.write('set ylabel \'fraction\'\n')
            gnu.write('set logscale cb\n')
            gnu.write('set xrange [%f:%f]\n' % (numpy.max(index1),numpy.min(index1)))
            gnu.write('set zrange[1E5:1E10]\n')
            gnu.write('set palette defined (0\'white\',0.1\'blue\',0.2\'green\',0.3\'orange\',0.4\'red\')\n')
            gnu.write('splot \'%s\' u 1:2:3\n' % ('test3D.txt'))
            """
            gnu.close()
            os.system('gnuplot figs/gnu.gp')

            #cleanup.
            os.system('rm test3D.txt')
            os.system('rm testP1.txt')
            os.system('rm testP2.txt')
            
            print 'deleting:',death
            for de in death:
                print de
                if(os.path.exists(str(de)+'/test.ft2')==True):
                    os.system('rm -rf '+str(de)+'/test.ft2')
                    os.system('rm -rf '+str(de)+'/test.fid')
                    os.system('rm -rf '+str(de)+'/test.out')

            


        else:
            os.system('rm integration.'+tag+'.out')
            os.system('rm raw.'+tag+'.data')

    def Analyse(self):
        #for i in range(len(self.tags)):
        #    Anal(filRan[i][0],filRan[i][1],tags[i],integ[i])

        gnu=open('gnu.gp','w')
        gnu.write('set term post eps enh color solid\n')
        gnu.write('set output \'figs/fig4.eps\'\n')
        gnu.write('set size square\n')
        gnu.write('set key outside\n')
        gnu.write('set title \'Intensities\'\n')
        gnu.write('set ylabel \'intensity\'\n')
        gnu.write('set xlabel \'delay(s)\'\n')
        gnu.write('set xrange[0:*]\n')
        #gnu.write('set yrange[0:*]\n')
        gnu.write('plot ')
        #for i,tag in enumerate(self.tag):
        #    if(i!=0):
        #        gnu.write(',')
        gnu.write('\'integration.%s.out\' u 1:3 ti\'%s\' w points lc %i,\'\' u 1:4 ti \'%s\' w points lc %i' % (self.tag,self.tag+' '+str(self.integ[0]),1,self.tag+' '+str(self.integ[1]),1))
        gnu.close()
        os.system('gnuplot gnu.gp')
        print self.before+'/py/arraygraph.py 2 4 0 0 0 0 `ls figs/*.eps`'
        os.system(self.before+'/py/arraygraph.py 2 4 0 0 0 0 `ls figs/*.eps`')
        os.system('mv summary.pdf '+self.tag+'.pdf') #move pdf new place


############################
