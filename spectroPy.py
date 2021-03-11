#/usr/local/bin/python3
#spectroPy v0.1
#(c) Pavol Gajdos, 11.3.2021

import sys
import time
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as mpl
from astropy.io import fits

name=sys.argv[1]
dat=[]
dat0=[]
res=1

colors='rgb'

def plot_spektrum():
    '''plot spectrum in 1-3 channels'''
    if len(dat)==0:
        print('Calibrate spectrum!')
        return

    mpl.figure()
    x=dat[:,0]
    if dat.shape[1]==2: ch=1
    else: ch=dat.shape[1]-1

    if ch==1: mpl.plot(x,dat[:,1],'k-')
    else:
        for i in range(ch):
            mpl.plot(x,dat[:,i+1],colors[i]+'-',label='Channel '+str(i+1))
            mpl.legend()
    mpl.xlabel('Wavelength (A)')
    mpl.ylabel('Intensity')
    mpl.show()

def save():
    '''save spectrum in 1-3 channels'''
    if len(dat)==0:
        print('Calibrate spectrum!')
        return
    ans=input('Save to file: csv (c) / text (t) / fits (f) ')
    name1=name[:name.rfind('.')]+'-spectrum'
    if dat.shape[1]==2: ch=1
    else: ch=dat.shape[1]-1
    if ans=='c':
        ans=input('Save to file "'+name1+'.csv"? (y/n) ')
        if ans=='y': name1+='.csv'
        else:
            name1=input('Save to file: ').strip()
            if len(name1)==0: return
        if ch==1: head=',Intensity'
        else:
            head=''
            for i in range(ch): head+=',Intensity_in_channel_'+str(i+1)
        np.savetxt(name1,dat,fmt='%f',delimiter=',',header='Wavelength(A)'+head,comments='')
        print('Saved')


    elif ans=='t':
        ans=input('Save to file "'+name1+'.dat"? (y/n) ')
        if ans=='y': name1+='.dat'
        else:
            name1=input('Save to file: ').strip()
            if len(name1)==0: return
        if ch==1: head='     Intensity'
        else:
            head=''
            for i in range(ch): head+='     Intensity in channel '+str(i+1)
        np.savetxt(name1,dat,fmt='%f',delimiter='     ',header='Wavelength (A)'+head)
        print('Saved')

    elif ans=='f':
        ans=input('Save to file "'+name1+'.fits"? (y/n) ')
        if ans=='y': name1+='.fits'
        else:
            name1=input('Save to file: ').strip()
            if len(name1)==0: return

        hdu=fits.PrimaryHDU(dat[:,1])
        hdu.header['CTYPE1']='Wavelength'
        hdu.header.comments['CTYPE1']='Axis Type'
        hdu.header['CUNIT1']='Angstrom'
        hdu.header.comments['CUNIT1']='Wavelength unit'
        hdu.header['CDELT1']=res
        hdu.header.comments['CDELT1']='Coordinate increment'
        hdu.header['CRPIX1']=1
        hdu.header.comments['CRPIX1']='Reference pixel'
        hdu.header['CRVAL1']=dat[0,0]
        hdu.header.comments['CRVAL1']='Coordinate at reference pixel'
        hdu.header['SWCREATE']='spectroPy'
        hdu.name='Intensity'
        if ch==1: hdu.writeto(name1,overwrite=True)
        else:
            hdu.name='Channel_1'
            hdu1=fits.ImageHDU(dat[:,2])
            hdu1.header['CTYPE1']='Wavelength'
            hdu1.header.comments['CTYPE1']='Axis Type'
            hdu1.header['CUNIT1']='Angstrom'
            hdu1.header.comments['CUNIT1']='Wavelength unit'
            hdu1.header['CDELT1']=res
            hdu1.header.comments['CDELT1']='Coordinate increment'
            hdu1.header['CRPIX1']=1
            hdu1.header.comments['CRPIX1']='Reference pixel'
            hdu1.header['CRVAL1']=dat[0,0]
            hdu1.header.comments['CRVAL1']='Coordinate at reference pixel'
            hdu1.header['SWCREATE']='spectroPy'
            hdu1.name='Channel_2'
            if ch==2:
                hdul=fits.HDUList([hdu,hdu1])
            else:
                hdu2=fits.ImageHDU(dat[:,3])
                hdu2.header['CTYPE1']='Wavelength'
                hdu2.header.comments['CTYPE1']='Axis Type'
                hdu2.header['CUNIT1']='Angstrom'
                hdu2.header.comments['CUNIT1']='Wavelength unit'
                hdu2.header['CDELT1']=res
                hdu2.header.comments['CDELT1']='Coordinate increment'
                hdu2.header['CRPIX1']=1
                hdu2.header.comments['CRPIX1']='Reference pixel'
                hdu2.header['CRVAL1']=dat[0,0]
                hdu2.header.comments['CRVAL1']='Coordinate at reference pixel'
                hdu2.header['SWCREATE']='spectroPy'
                hdu2.name='Channel_3'
                hdul=fits.HDUList([hdu,hdu1,hdu2])
            hdul.writeto(name1,overwrite=True)
        print('Saved')

def fit_gauss(x,y,sgn):
    '''fit by Gauss profile'''
    def gauss(x,amp,x0,sigma,a,b):
        return amp*np.exp(-(x-x0)**2/(2*sigma**2))+a*x+b

    i=list(range(10))+list(range(-10,0,1))
    p=np.polyfit(x[i],y[i],1)
    yy=y-np.polyval(p,x)

    if sgn>0:
        amp=max(yy)-min(yy)
        i=np.argmax(yy)
    else:
        amp=min(yy)-max(yy)
        i=np.argmin(yy)

    mean=x[i]
    sigma=np.sqrt(sum(y*(x-mean)**2)/sum(y))
    try:
        popt,pcov=curve_fit(gauss,x,y,p0=[amp,mean,sigma,p[0],p[1]])
    except RuntimeError: return False

    mpl.figure()
    mpl.plot(x,y,'b.')
    x1=np.arange(x[0],x[-1],0.01)
    mpl.plot(x1,gauss(x1,*popt),'r-')
    mpl.show()

    return popt[1]

def plot(data,chI=None,title=None,manipulate=True):
    '''function for ploting for crop and calibrate'''
    MAX_CLICK = 0.5 # in seconds; anything longer is a drag motion

    global xmin,xmax,time_onclick,figs

    def onclick(event):
        global time_onclick
        time_onclick=time.time()

    def onNoclick(event):
        global xmin,xmax,figs

        if time.time()-time_onclick>MAX_CLICK: return
        xx=event.xdata
        xlim=mpl.xlim()
        ylim=mpl.ylim()
        if event.button==1: xmin=xx
        elif event.button==3: xmax=xx

        for i in figs: i.pop(0).remove()

        figs=[]
        if xmin>-1e10:
            figs.append(mpl.plot([xmin,xmin],[ylim[0],ylim[1]],'r--'))
        if xmax>-1e10:
            figs.append(mpl.plot([xmax,xmax],[ylim[0],ylim[1]],'r--'))
        if xmax>xmin and xmin>-1e10:
            figs.append(mpl.plot([xmin,xmax],[(ylim[0]+ylim[1])/2,(ylim[0]+ylim[1])/2],'r--'))
        #print x[ind][0],y[ind][0]
        fig.canvas.draw()
        mpl.xlim(xlim)
        mpl.ylim(ylim)

    xmin=-1e10
    xmax=-1e10
    fig=mpl.figure()

    figs=[]

    if data.shape[1]==2: ch=1
    else: ch=data.shape[1]-1

    if chI is None:
        if ch==1: mpl.plot(data[:,0],data[:,1],'k-')
        else:
            for i in range(ch): mpl.plot(data[:,0],data[:,i+1],colors[i]+'-',label='Channel '+str(i+1))
            mpl.legend()
    else: mpl.plot(data[:,0],data[:,chI],'k-')

    if manipulate:
        fig.canvas.mpl_connect('button_press_event',onclick)
        fig.canvas.mpl_connect('button_release_event',onNoclick)

    if title is not None: mpl.title(title)
    mpl.show()

    return xmin,xmax

def crop():
    '''crop spectrum'''
    global dat

    if len(dat)==0:
        print('Calibrate spectrum!')
        return

    print('\nSelect region to crop (keep it):')
    print('Left click - left border')
    print('Right click - right border')

    i1,i2=plot(dat,title='Crop spectrum')

    ans=input('Already crop spectrum? (y/n) ')
    if ans=='y':
        x=dat[:,0]

        if i1>=i2: return
        i=np.where((x>i1)*(x<i2))[0]

        dat=dat[i,:]

def calibrate(n=1):
    '''calibrate using 1/2 points'''
    global dat,res

    x=dat0[:,0]
    if dat0.shape[1]==2: ch=1
    elif dat0.shape[1]<=4: ch=dat0.shape[1]-1
    else: ch=3

    plot(dat0,title='Input data',manipulate=False)

    if ch==1:
        chI=1
        chs=[1]
    else:
        chs=' / '.join([str(x+1) for x in range(ch)])
        if ch>1: chs+=' / a (all)'
        ans=input('Color channel to fit - '+chs+': ')
        if ans=='a': chs=list(range(1,ch+1))
        else: chs=[int(ans)]

    pos0=[]
    for chI in chs:
        ans1='y'  #repeat?
        while ans1=='y':
            print('\nSelect region with 0th order spectra:')
            print('Channel',chI)
            print('Left click - left border')
            print('Right click - right border')
            i1,i2=plot(dat0,chI,'Select region with 0th order spectra for channel '+str(chI))
            if i1>=i2: return
            i=np.where((x>i1)*(x<i2))[0]

            pos0i=fit_gauss(x[i],dat0[i,chI],1)

            if pos0i:
                ans=input('Is fitted profile good? (y/n) ')
                if ans=='y':
                    pos0.append(pos0i)
                    print('Peak possition:',round(pos0i,2))
                    ans1='n'
                else:
                    ans1=input('Repeat fitting? (y/n) ')
                    if ans1=='n':return
            else: input('Problem with fitting. Select different region.')
    pos0=np.average(pos0)
    if len(chs)>1: print('Peak possition:',round(pos0,2))

    if n==1:
        ans=input('Spectral resolution (A/px): ')
        res=float(ans)

    if n==2:
        pos1=[]
        for chI in chs:
            ans1='y'  #repeat?
            while ans1=='y':
                print('\nSelect region with known spectral line:')
                print('Channel',chI)
                print('Left click - left border')
                print('Right click - right border')
                i1,i2=plot(dat0,chI,'Select region with known spectral line for channel '+str(chI))
                if i1>=i2: return
                i=np.where((x>i1)*(x<i2))[0]

                pos1i=fit_gauss(x[i],dat0[i,chI],-1)

                if pos1i:
                    ans=input('Is fitted profile good? (y/n) ')
                    if ans=='y':
                        pos1.append(pos1i)
                        print('Peak possition:',round(pos1i,2))
                        ans1='n'
                    else:
                        ans1=input('Repeat fitting? (y/n) ')
                        if ans1=='n':return
                else: input('Problem with fitting. Select different region.')
        pos1=np.average(pos1)
        if len(chs)>1: print('Peak possition:',round(pos1,2))

        ans=input('Wavelength of line (A): ')
        res=float(ans)/(pos1-pos0)
        print('Spectral resolution (A/px):',round(res,3))

    lam=(x-pos0)*res
    dat=np.zeros(dat0.shape)
    dat[:,0]=lam
    for i in range(ch): dat[:,i+1]=dat0[:,i+1]

    print('Calibrated')
    return lam

def menu():
    print('calibrate - 1 point: 1')
    print('calibrate - 2 points: 2')
    print('crop spectrum: 3')
    print('plot spectrum: 4')
    print('save spectrum: 5')
    print('exit: 0')
    print('========================')
    opt=input('option: ')
    if int(opt)==0:
        #exit()
        return
    elif int(opt)==1: calibrate(1)
    elif int(opt)==2: calibrate(2)
    elif int(opt)==3: crop()
    elif int(opt)==4: plot_spektrum()
    elif int(opt)==5: save()
    print('========================\n')
    menu()

print('spectroPy v0.1 - (c) Pavol Gajdos 2021\n')

if '.fit' in name:
    #load fits image
    f=fits.open(name)
    #image/spectrum?
    #TODO...

else:
    #load spectrum
    f=open(name,'r')
    lines=f.readlines()
    f.close()

    #detect header
    header=0
    for l in lines:
        if l.strip()[0] not in '-+0123456789.#': header+=1
        else: break

    #detect separator type
    l=lines[20]
    sep=''
    for x in l.strip():
        if x in '-+0123456789.':
            if len(sep)>0: break
        else: sep+=x
    if len(sep)==0: sep=' '

    dat0=np.loadtxt(name,delimiter=sep,skiprows=header)
    if len(dat0.shape)==1:
        #only intensity
        dat0=np.column_stack((dat0,dat0))
        dat0[:,0]=np.arange(0,len(dat0[:,1]),1)
    if min(dat0[:,0])==0 and max(dat0[:,0])==0:
        #input from vSpec
        dat0[:,0]=np.arange(0,len(dat0[:,1]),1)

menu()
