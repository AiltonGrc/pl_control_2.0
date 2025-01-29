# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 12:00:09 2019

@author: Schimpf
"""

from matplotlib import pyplot as plt
import numpy as np
from scipy import signal


def PeakAreas(xvals,yvals,peak_bin=2000,resolution=128,num_sidepeaks=5,period=12500):
    """
    Calculated peak areas of given 1D histogram.
    ATTENTION: Middle peak has to be at x=0
    
    Returns:
    --------
    (middle peak area,
    average side peak areas,
    index of middle peak;
    indices of side peaks)

    Required Parameters:
    --------------------
    xvals ... x-values of histogram
    yvals ... y-values of histogram
    
    Optional Parameters:
    --------------------
    peak_bin = 2000 # area binning in ps
    resolution = 128 # resolution in ps
    num_sidepeaks = 5 # number of side peaks to consider
    period = 12500 # repetition period in ps
    """

    index_step=period/resolution
    
    indices_bin = peak_bin//(2*resolution)+1
    middlePeakIndex = np.argmin(np.abs(np.array(xvals)))
    middlePeak_Area = yvals[middlePeakIndex-indices_bin:middlePeakIndex+indices_bin].sum()
    
    """ OLD: find sidepeak indices with peakfinder
    peaks = signal.find_peaks(yvals,distance=(2/3 * index_step))
    mask= np.logical_and((np.abs(middlePeakIndex-peaks[0]) > 2*indices_bin),(np.abs(middlePeakIndex-peaks[0]) < num_sidepeaks*period/resolution))
    sidepeak_indices = peaks[0][mask]
   
    """
    
    """ NEW: Find sidepeaks by repetition rate """
    sidepeak_indices=np.array([middlePeakIndex + index_step*i for i in range(-num_sidepeaks,num_sidepeaks+1) if i not in [0]], dtype=int)
            
    sidepeak_areas = np.array(list(map(lambda x: yvals[x-indices_bin:x+indices_bin].sum(),sidepeak_indices)))

    return middlePeak_Area, np.average(sidepeak_areas), middlePeakIndex, sidepeak_indices

def g2(A_middle,A_side):
    """
    Calculates g2 and error based on poisson distribution
    
    Returns:
    --------
    (g2, stDev)
    
    Required Parameters:
    --------------------
    A_middle ... Middle peak area of histogram
    A_side ... Average side peak area of histogram       
    """
    
    g2 = A_middle/A_side
    stDev = np.sqrt( (np.sqrt(A_middle)/A_side)**2 + (np.sqrt(A_side)*A_middle/A_side**2)**2)
    
    return g2,stDev

def PairGenEfficiency(A0,As):
    """
    Calculated the pair generation efficiency from a set of X/XX cross-correlation measurement
    
    Returns:
    --------
    (Efficiency, Error)
    
    Parameters:
    -----------
    A0 ... middle peak area
    As ... average side peak area
    """
    eff = 1/(A0/As).mean()
    vsum = np.array([A0[i]/As[i]**2 + A0[i]**2/As[i]**3 for i in range(len(A0))]).sum()
    effErr = eff**2/len(A0) * np.sqrt(vsum)    
    
    return eff,effErr

def PlotCorr(xvals,yvals,resolution=512,timebin=2000,plot=True):
    
    middlePeak_Area, av_sidepeak_areas, middlepeak_index, sidepeak_indices  = PeakAreas(xvals,yvals,resolution=resolution,peak_bin=timebin)
    g2val, g2err = g2(middlePeak_Area,av_sidepeak_areas)
    
    print(f'Middlepeak area: {middlePeak_Area}')
    print(f"Average sidepeak area: {av_sidepeak_areas}")
    
    print(f"g2val: {g2val:.3f}({g2err:.3f}) ")
    
    
    plt.rcParams['figure.dpi'] = 600
    plt.rcParams['savefig.dpi'] = 600

    if plot == True:
        f = plt.figure(figsize=(6,4))
        ax = f.add_subplot(111)
        # titel
        plt.title(label=title,fontsize=10)
        ax.plot(xvals/1E3,yvals,color='blue')
        ax.vlines(xvals[sidepeak_indices]/1E3,0,max(yvals), linestyle='--')
        ax.vlines(np.array([xvals[middlepeak_index]-timebin/2,xvals[middlepeak_index]+timebin/2])/1E3,0,max(yvals), linestyle='--', color="red")
        ax.set_xlabel('Time delay (ns)')
        ax.set_ylabel('Coincidences')
        f.show()


    return middlePeak_Area, av_sidepeak_areas, g2val, g2err, middlepeak_index, sidepeak_indices


hydraharp = True
resolution = 512  # Resolution in ps


if hydraharp:
    path = r'D:\Users\AK124573\Desktop\InGaAs_Paper\Correlation\AutoCorrelation\SA0919_qd2_v2_X_autocorr_green_0.94mW_4,7K_512ps_bin_785,12nm_1200g_80µm.dat'
    xrange = 200  # range in ns
    yvals = np.loadtxt(path, skiprows=10, usecols=[0])
    # xvals=np.loadtxt(path,delimiter=",", skiprows=0, usecols=[0])

    offs_optimum = 90000
    timebin_optimum = 1024
    g2_init = 0.9

    for _ in np.arange(7.5e4,1e5,1e2):

        print("offset: ",_)
        offs = _
        timebin = 1024  # Timebin in ps
        yvals = yvals[:int(xrange / resolution * 1000)]
        xvals = np.array([resolution * i - offs for i in range(len(yvals))])
        split = path.split("/")
        title = split[len(split) - 1]

        if _ >= 1e5: Plot = True
        else: Plot = False

        middlePeak_Area, av_sidepeak_areas, g2val, g2err, middlepeak_index, sidepeak_indices = \
            PlotCorr(xvals, yvals, resolution, timebin, Plot)
        if g2val<g2_init:
            g2_init = g2val
            offs_optimum = offs

    for _ in np.arange(1024,2**10,128):

        print("timebin: ", _)
        offs = offs_optimum
        timebin = _  # Timebin in ps
        yvals = yvals[:int(xrange / resolution * 1000)]
        xvals = np.array([resolution * i - offs for i in range(len(yvals))])
        split = path.split("/")
        title = split[len(split) - 1]

        if _ >= 2048: Plot = True
        else: Plot = False

        middlePeak_Area, av_sidepeak_areas, g2val, g2err, middlepeak_index, sidepeak_indices = \
            PlotCorr(xvals, yvals,resolution,timebin,Plot)
        if g2val < g2_init:
            g2_init = g2val
            timebin_optimum = timebin

    offs = offs_optimum
    timebin = timebin_optimum  # Timebin in ps
    yvals = yvals[:int(xrange / resolution * 1000)]
    xvals = np.array([resolution * i - offs for i in range(len(yvals))])
    split = path.split("/")
    title = split[len(split) - 1]
    middlePeak_Area, av_sidepeak_areas, g2val, g2err, middlepeak_index, sidepeak_indices = \
        PlotCorr(xvals, yvals,resolution,timebin,True)

    print('Optimum values: {} offset, {} binning'.format(offs,timebin))
    print("g2 max: ",g2val)