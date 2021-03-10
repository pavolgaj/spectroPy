import sys
import numpy as np
import scipy
import matplotlib.pyplot as mpl
from astropy.io import fits as pyfits

name=sys.argv[1]
dat=[]

def plot_spektrum(dat):
    '''plot spektrum v 1-3 kanaloch'''
    return
    
def save(dat):
    '''ulozit spektrum v 1-3 kanaloch'''
    return

def fit_gauss(x,y,A):
    '''fitovanie gaussovym profilom'''
    return
    
def calibrate(dat,n=1):
    '''kalibracia pomocou 1/2 bodov'''
    
    fit_gauss(x,y,A)
    return
    
def menu():
    print('calibrate - 1 point: 1')
    print('calibrate - 2 points: 2')
    print('save spectrum: 3')
    print('plot spectrum: 4')
    print('exit: 0')
    print('========================')
    volba=input('option: ')
    if int(volba)==0: exit()
    elif int(volba)==1: calibrate(dat,1)
    elif int(volba)==2: calibrate(dat,2)
    elif int(volba)==3: save(dat)
    elif int(volba)==3: plot_spektrum(dat,3)
    print('========================\n')
    menu()    
    

if '.fit' in name:
    #load fits image
    fits=pyfits.open(name)
else:
    #load spektrum
    f=open(name,'r')
    lines=f.readlines()
    f.close()
    
menu()    
