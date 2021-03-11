import sys
import time
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as mpl
from astropy.io import fits as pyfits

name=sys.argv[1]
dat=[]
dat0=[]

colors='rgb'

def plot_spektrum():
    '''plot spectrum in 1-3 channels'''
    mpl.figure()
    x=dat[:,0]
    if dat.shape[1]==2: ch=1
    else: ch=dat0.shape[1]-1

    if ch==1: mpl.plot(x,dat[:,1],'k-')
    else:
        for i in range(ch): mpl.plot(x,dat[:,i+1],colors[i]+'-')
    mpl.xlabel('Wavelength (A)')
    mpl.show()

def save():
    '''save spectrum in 1-3 channels'''
    return

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

def calibrate(n=1):
    '''calibrate using 1/2 points'''
    global dat

    def plot(chI=None,title=None):
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

        if chI is None:
            if ch==1: mpl.plot(x,dat0[:,1],'k-')
            else:
                for i in range(ch): mpl.plot(x,dat0[:,i+1],colors[i]+'-')
        else:
            mpl.plot(x,dat0[:,chI],'k-')
            fig.canvas.mpl_connect('button_press_event',onclick)
            fig.canvas.mpl_connect('button_release_event',onNoclick)

        if title is not None: mpl.title(title)
        mpl.show()

        return xmin,xmax

    x=dat0[:,0]
    if dat0.shape[1]==2: ch=1
    elif dat0.shape[1]<=4: ch=dat0.shape[1]-1
    else: ch=3

    plot(title='Input data')

    if ch==1: chI=1
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
            i1,i2=plot(chI,'Select region with 0th order spectra for channel '+str(chI))
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
    print('Peak possition:',round(pos0,2))

    if n==1:
        ans=input('Spectral resolution (A/px): ')
        res=float(ans)

    if n==2:
        pos1=[]
        for chI in chs:
            ans1='y'  #repeat?
            while ans1=='y':
                print('\nSelect region with known spectral line:')
                print('Left click - left border')
                print('Right click - right border')
                i1,i2=plot(chI,'Select region with known spectral line for channel '+str(chI))
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
        print('Peak possition:',round(pos1,2))

        ans=input('Wavelength of line (A): ')
        res=float(ans)/(pos1-pos0)
        print('Spectral resolution (A/px):',round(res,3))

    lam=(x-pos0)*res
    dat=np.zeros(dat0.shape)
    dat[:,0]=lam
    for i in range(ch): dat[:,i+1]=dat0[:,i+1]

    return lam

def menu():
    print('calibrate - 1 point: 1')
    print('calibrate - 2 points: 2')
    print('save spectrum: 3')
    print('plot spectrum: 4')
    print('exit: 0')
    print('========================')
    opt=input('option: ')
    if int(opt)==0:
        #exit()
        return
    elif int(opt)==1: calibrate(1)
    elif int(opt)==2: calibrate(2)
    elif int(opt)==3: save()
    elif int(opt)==4: plot_spektrum()
    print('========================\n')
    menu()


if '.fit' in name:
    #load fits image
    fits=pyfits.open(name)
    #image/spectrum?

else:
    #load spectrum
    f=open(name,'r')
    lines=f.readlines()
    f.close()

    #detect separator type
    l=lines[20]
    sep=''
    for x in l.strip():
        if x in '0123456789.':
            if len(sep)>0: break
        else: sep+=x

    #detect header
    header=0
    for l in lines:
        if l.strip()[0] not in '0123456789.#': header+=1
        else: break

    dat0=np.loadtxt(name,delimiter=sep,skiprows=header)

menu()
