# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 15:33:18 2022

@author: Gabriel Undeutsch

read in spe2.x(BigLab) files and spe3.x(MagnetoLab) files
edited version of https://github.com/HydraSpex/SPE_Handler : Read in of the files

class Spectra
class variables include intensity, wavelength, nr. of pixel,...

functions:
    
    

"""

#can read spe2.x files and spe3.0 files
#edited version of https://github.com/HydraSpex/SPE_Handler


import struct
import numpy as np
import matplotlib.pylab as plt
from matplotlib import ticker
from matplotlib.widgets import Slider
import time
from matplotlib import animation

#The following 3 functions are from https://github.com/HydraSpex/SPE_Handler and are slightly eddited
#They read in the information from the spe files, which incluedes the header
#Reading Bytes
def from_bytes(b, format, offset):
    calcsize = struct.calcsize(format)
    return struct.unpack(format, b[offset:offset+calcsize])[0]

#Reading the XML Line of SPE3.x Files
def openXMLline(filename):
    num_lines = sum(1 for line in open(filename, "rb"))
    f = open(filename, "rb")
    lines = f.readlines()
    
    Txt_Dai = open("Dai_new.txt", "w")  

    f = open(filename, "rb")
    lines = f.readlines()

    i = 0
    while i < num_lines:
        Txt_Dai.write(str(lines[i]) + "\n")
        i += 1
    Txt_Dai.close()

    XMLline = str(lines[num_lines-1])
    startXML = XMLline.find("<SpeFormat")
    XMLline = XMLline[startXML:]
    
    PosVersion = XMLline.find("version") + 9
    Version = XMLline[PosVersion:PosVersion+3]
   
    PosFrame = XMLline.find("Frame") + 14
    Frame = XMLline[PosFrame:]
    PosFrameEnd = Frame.find("\"")
    Frame = Frame[:PosFrameEnd]
    
    PosWidth = XMLline.find("width") + 7
    Width = XMLline[PosWidth:]
    PosWidthEnd = Width.find("\"")
    Width = Width[:PosWidthEnd]
    
    PosHeight = XMLline.find("height") + 8
    Height = XMLline[PosHeight:]
    PosHeightEnd = Height.find("\"")
    Height = Height[:PosHeightEnd]
    
    PosLaser = XMLline.find("WavelengthLaserLine")
    Laser = XMLline[PosLaser:]
    PosLaser = Laser.find("\">") + 2
    Laser = Laser[PosLaser:]
    PosLaserEnd = Laser.find("<")
    Laser = Laser[:PosLaserEnd]
    
    PosCreate = XMLline.find("created") + 9
    Created = XMLline[PosCreate:]
    PosTime = Created.find("T")
    Date = Created[:PosTime]
    PosTimeEnd = Created.find(".")
    Time = Created[PosTime+1:PosTimeEnd]
    
    PosWave = XMLline.find("<Wavelength ")
    PosWaveEnd = XMLline.find("</Wavelength>")
    WaveData = XMLline[PosWave:PosWaveEnd]
    WaveData = WaveData[(WaveData.find("\">")+2):]
    
    Wavedata = []
    WavedataRound = []
    i = 0
    
    while i < int(Width):
        PosNext = WaveData.find(",")
        Wavedata.append(WaveData[:PosNext])
        WavedataRound.append(round(float(WaveData[:PosNext]), 3))
        WaveData = WaveData[PosNext+1:]
        i += 1

    PosExpTime = XMLline.find("<ExposureTime")
    PosExpTimeEnd = XMLline.find("</ExposureTime>")
    ExpTimeData = XMLline[PosExpTime:PosExpTimeEnd]
    ExpTime = ExpTimeData[(ExpTimeData.find("\">")+2):]
    ExpTime = int(float(ExpTime)/1000)
    
    PosCWL = XMLline.find("<CenterWavelength")
    PosCWLEnd = XMLline.find("</CenterWavelength>")
    CWLData = XMLline[PosCWL:PosCWLEnd]
    CWL = CWLData[(CWLData.find("\">")+2):]

    PosBG = (XMLline.find("<BackgroundCorrection><Enabled") + 25)
    BGData = XMLline[PosBG:]
    PosBGBegin = (BGData.find(">")+1)
    PosBGEnd = BGData.find("<")
    BG = BGData[PosBGBegin:PosBGEnd]

    PosGrating = XMLline.find("<Grating><Selected")
    GratingData = XMLline[PosGrating:]
    PosGratingBegin = GratingData.find("[")
    PosGratingEnd = (GratingData.find("]")+1)
    Grating = GratingData[PosGratingBegin:PosGratingEnd]

    return Version, int(Frame), int(Width), int(Height), Laser, Date, Time, ExpTime, CWL, Grating, BG, Wavedata, WavedataRound


def getData(bytes,filename):
    # file = open(filename, "rb")
    # bytes = file.read()

    SPEVersion = round(from_bytes(bytes, "f", 1992),1)

    datatype = from_bytes(bytes, "h", 108)
    frame_width = from_bytes(bytes, "H", 42)
    frame_height = from_bytes(bytes, "H", 656)
    num_frames = from_bytes(bytes, "i", 1446)
    to_np_type = [np.float32, np.int32, np.int16, np.uint16, None, np.float64, np.uint8, None, np.uint32]
    np_type = to_np_type[datatype]
    itemsize = np.dtype(np_type).itemsize
    XMLOffset = from_bytes(bytes, "64Q", 678)

    Count = frame_width * frame_height

    if SPEVersion >= 3:
        Version, Frame, Width, Height, Laser, Date, Time, ExpTime, CWL, Grating, BG, Wavedata, WavedataRound = openXMLline(filename)
    else:
        Version = SPEVersion
        Frame = num_frames
        Width = frame_width
        Height = frame_height
        Laser = from_bytes(bytes, "d", 3311)
        LocalDate = from_bytes(bytes, "16s", 20)
        LocalDate = LocalDate.decode("utf-8", "ignore")
        Day = LocalDate[:2]
        Year = LocalDate[5:9]
        Month = LocalDate[2:5]
        if Month == "Jan":
            Month = "01"
        elif Month == "Feb":
            Month = "02"
        elif Month == "Mar":
            Month = "03"
        elif Month == "Apr":
            Month = "04"
        elif Month == "May":
            Month = "05"
        elif Month == "Jun":
            Month = "06"
        elif Month == "Jul":
            Month = "07"
        elif Month == "Aug":
            Month = "08"
        elif Month == "Sep":
            Month = "09"
        elif Month == "Oct":
            Month = "10"
        elif Month == "Nov":
            Month = "11"
        elif Month == "Dec":
            Month = "12"
        else:
            Month = "00"
        Date = Year + "-" + Month + "-" + Day
        LocalTime = from_bytes(bytes, "6s", 172)
        LocalTime = LocalTime.decode("utf-8", "ignore")
        Time = LocalTime[:2] + ":" + LocalTime[2:4] + ":" + LocalTime[4:]
        UTCTime = from_bytes(bytes, "6s", 179)
        UTCTime = UTCTime.decode("utf-8", "ignore")
        ExpTime = from_bytes(bytes, "f", 10)
        CWL = from_bytes(bytes, "f", 72)
        Grating = from_bytes(bytes, "32f", 650)
        BG = from_bytes(bytes, "i", 150)
        XStartNM = from_bytes(bytes, "d", 3183)
        XCenNM = from_bytes(bytes, "d", 3199)
        
        Wavedata = np.linspace(XStartNM,(XCenNM-XStartNM)*2+XStartNM,Width)
        WavedataRound = np.round(Wavedata,3)
        

    return np_type, itemsize, Count, Version, Frame, Width, Height, Laser, Date, Time, ExpTime, CWL, Grating, BG, Wavedata, WavedataRound



class Spectra:
    '''
    Functions:
    -----------------------------
    plot_spec(i=0): 
    Plot the i-th spectrum
    
    plot_specs(cmap = 'bone',log=False):
    Plot a map of the spectra (IV_trace, Polmap,Powermap,...), optionally in log scale. Name must be as is for XRSP to show the axis properly.
    
    plot_image(cmap = 'bone', log = False):
    Plot image of the CCD recorded image, optinally in log scale. Can also be used for 0-th order image.
    
    animate_images(i=0, fname = None, save = True, log=False,cmap = 'GnBu',y_offset = 0,y_range=[0.45,0.55],x_range=[0,1],dpi = 600)
    Create and save an animation of several images (Polmap, Powermap,...)
    '''
    
    def __init__(self,filename,pixelsize = 20,adoffset = 600):
        '''
        reads in the spe file from the given filename
       
        Parameters:
           filename(str): define the filename of the data
           pixelsize(float)(default = 20): pixelsize of the CCD camera in µm (Pylon: 20µm, iDus: 15µm)
           adoffset(int)(default = 600): Offset of the AD Converter of the CCD
        Returns:
           None
        ''' 
        #open and read in the binary spe file
        file = open(filename, "rb")
        bytes = file.read()
        
        #read in the data
        np_type, Itemsize, Count, Version, Frame, Width, Height, Laser, Date, Time, ExpTime, CWL, Grating, BG, Wavedata, WavedataRound = getData(bytes,filename)   
        Itemsize = int(Itemsize)
        
        #read in the intensities
        intensity = []
    
        for i in range(0,Frame):
            offset = 4100 + i*Count*Itemsize
            inte = np.frombuffer(bytes, dtype=np_type, count=Count, offset = offset)
            inte = np.uint16(inte)
            if Height > 1:
                inte = inte.reshape(Height,Width)
            intensity.append(inte)
        
        #if the file consists of several spectra:
        if Frame>1:
            intensity = np.array(intensity)
        else:
            intensity = np.array(intensity)[0]
        
        self.filename = filename
        
        #make sure there is no underflow since the intenstity is in uint16
        intensity[intensity<adoffset]=adoffset
        self.intensity = intensity-adoffset
        self.wavelength = WavedataRound
        
        self.adoffset = adoffset
        self.pixelsize = pixelsize 
        
        self.np_type = np_type
        self.Itemsize = Itemsize
        self.Count = Count
        self.Version = Version
        self.Frame = Frame
        self.Width = Width
        self.Height = Height
        self.Laser = Laser
        self.Data = Date
        self.Time = Time
        self.ExpTime = ExpTime
        self.CWL = CWL
        self.Grating = Grating
        self.BG = BG
        
    
    def plot_spec(self,i=0):
        '''
        Plot a spectrum
       
        Parameters:
           i(int)(default = 0): plot the i-th spectrum
        Returns:
           None
        '''
        
        #check if i is set to a too high value
        if i>self.Frame-1:
            print('Error: i out of range')
            return
        #set x and y axis
        x = self.wavelength
        if self.Frame == 1:
            y = self.intensity
        else:
            y = self.intensity[i]
            
        #plot spectrum
        fig,ax = plt.subplots(figsize=(15,9))
    
        ax.plot(x,y, c='k')
        
        ax.set_xlabel('Wavelength[nm]', fontsize = 28)
        ax.set_ylabel('Counts', fontsize = 28)
        ax.tick_params(axis='x', labelsize=24)
        ax.tick_params(axis='y', labelsize=24)
        ax.set_xlim(min(self.wavelength),max(self.wavelength))
        ax.set_ylim(max(y)*-0.03,max(y)*1.1)
        ax.grid()
        
        return
    
    def plot_specs(self,cmap = 'bone',log=False,switchxy = False):
        '''
        Plot a map of spectra
       
        Parameters:
           cmap(str)(default = 'bone'): Choose a colormap for the map
           log(bool)(default = False): Choose whether to plot the intensities (Z-values) in logscale or linscale
           switchxy(bool)(default = False): Choose whether to switch the x and y axis around
        Returns:
           None
        ''' 
        #set x values
        x = self.wavelength
        #set y values from file name
        try:
            r = self.filename.split('#')
            r = r[-1].split('_')
            r1 = float(r[1])
            r2 = float(r[2].split('.')[0])
            y=np.linspace(r1,r2,self.Frame)
        except IndexError:
            print('y-axis could not be built from filename, set to default values.')
            y=np.linspace(0,self.Frame-1,self.Frame)
        
        #plot map
        if switchxy:
            X, Y = np.meshgrid(y,x)
        else:
            X, Y = np.meshgrid(x,y)
        
        Z = self.intensity
        
        if switchxy:
            Z = Z.T
    
        fig = plt.figure(figsize=(28,15))
        ax = fig.add_subplot(111)
        
        if log:
            plt.contourf(X, Y, Z, 100,locator=ticker.LogLocator(), cmap = cmap)
        else: 
            plt.contourf(X, Y, Z, 100, cmap = cmap)
        ax.tick_params(axis='x', labelsize=40)
        ax.tick_params(axis='y', labelsize=40)
        
        if switchxy:
            ax.set_xlabel(r[0], fontsize = 60)
            ax.set_ylabel('Wavelength[nm]', fontsize = 60)
        else:
            ax.set_ylabel(r[0], fontsize = 60)
            ax.set_xlabel('Wavelength[nm]', fontsize = 60)
        
        cbar = plt.colorbar(format='%.2f')
        cbar.ax.tick_params(labelsize = 30)
        cbar.set_label('Counts', fontsize = 60)
        plt.tight_layout()
        
        return
        
    def plot_image(self, cmap = 'bone', log = False, i=0):
        '''
        Plot an image taken from a CCD Camera. Automatically checks if the image is taken in 0-th order mode or not. Depending on that, it displays the y-axis in µm or in wavelength 
       
        Parameters:
           cmap(str)(default = 'bone'): Choose a colormap for the map
           log(bool)(default = False): Choose whether to plot the intensities (Z-values) in logscale or linscale
           i(int)(default = 0): plot the i-th image
        Returns:
           None
        '''  
        #check if i is set to a too high value
        if i>self.Frame-1:
            print('Error: i out of range')
            return
        #check if the image is 0-th order or not
        if sum(self.wavelength) > 10_000:
            zero = False
        else:
            zero = True
            
        #set x value
        #if 0-th order, set to size in µm
        if zero:
            x = np.linspace(-(len(self.wavelength)/2-1)*self.pixelsize,len(self.wavelength)/2*self.pixelsize,len(self.wavelength)) *1e-3
        #else: set to wavelength
        else:
            x=self.wavelength
        #set y axis to size in µm
        y=np.linspace(0,(self.Height-1)*self.pixelsize,self.Height)*1e-3
        X, Y = np.meshgrid(x,y)
        
        #set z values
        if self.Frame == 1:
            Z = self.intensity
        else:
            Z = self.intensity[i]
        
        #plot image
        fig = plt.figure(figsize=(28,15))
        ax = fig.add_subplot(111)
        
        if log:
            plt.contourf(X, Y, Z, 100,locator=ticker.LogLocator(), cmap = cmap)
        else: 
            plt.contourf(X, Y, Z, 100, cmap = cmap)
        ax.tick_params(axis='x', labelsize=40)
        ax.tick_params(axis='y', labelsize=40)
        ax.set_ylabel('height[mm]', fontsize = 60)
        if zero:
            ax.set_xlabel('x[mm]', fontsize = 60)
        else:
            ax.set_xlabel('Wavelength[nm]', fontsize = 60)
        
        cbar = plt.colorbar(format='%.2f')
        cbar.ax.tick_params(labelsize = 30)
        cbar.set_label('Counts', fontsize = 60)
        plt.tight_layout()
        
        return
        
    
    
    def animate_images(self,i=0, save = True, fname = None, log=False,cmap = 'GnBu',y_range=[0.45,0.55],x_range=[0,1],dpi = 600):        
        '''
        Show and save the animation of several images. If possible choose a smaller size of the plotted range to save computation time.
        
        Parameters:
           save(bool)(default = True): Choose if the animation should be saved
           fname(str)(default = None): If left to None, Animation is saved to same folder, where the spe file is, otherwise to the specified location. Ends with '.gif'
           log(bool)(default = False): Choose whether to plot the intensities (Z-values) in logscale or linscale
           cmap(str)(default = 'GnBu'): Choose a colormap for the map
           y_range(float,float)(default = (0.45,0.55)): Choose the range in y in which the plots are shown.
           x_range(float,float)(default = (0,1)): Choose the range in x in which the plots are shown.
           dpi(float)(default = 600): Choose a dpi for the resolution. If set too high without limiting the y and x range as well as high number of images leads to high computational effort
        Returns:
           None
        ''' 
        #check if the image is 0-th order or not          
        if sum(self.wavelength) > 10_000:
            zero = False
        else:
            zero = True
        #if 0-th order, set to size in µm, with x_range   
        if zero:
            x = np.linspace(-(len(self.wavelength)/2-1)*self.pixelsize,len(self.wavelength)/2*self.pixelsize,len(self.wavelength)) *1e-3
            x = x[int(len(x)*x_range[0]):int(len(x)*x_range[0])]
        #else: set to wavelength, with x_range
        else:
            x=self.wavelength[int(self.Width*x_range[0]):int(self.Width*x_range[1])]
        #set y axis to size in µm, with y_range  
        y=np.linspace(0,(self.Height-1)*self.pixelsize,self.Height)*1e-3
        y=y[int(self.Height*y_range[0]):int(self.Height*y_range[1])]
        X, Y = np.meshgrid(x,y)
        
        #plot figure for animation
        fig = plt.figure(figsize=(20,5))

        #function that creates the animation
        def animate(i):
            
            #clear the last figure and plot the next
            fig.clear()
            ax = fig.add_subplot(111)
            #set the Z values of the i-th image
            Z = self.intensity[i]
            Z = Z[int(self.Height*y_range[0]):int(self.Height*y_range[1])]
            Z = Z[:,int(self.Width*x_range[0]):int(self.Width*x_range[1])]
            if log:
                self.cont = ax.contourf(X, Y, Z, 100,locator=ticker.LogLocator(), cmap = cmap)
            else: 
                self.cont = ax.contourf(X, Y, Z, 100, cmap = cmap)
            #change the counter in the bottom right
            ax.text(min(x)+0.1,min(y)+0.1,i,fontsize = 36)
            ax.tick_params(axis='x', labelsize=28)
            ax.tick_params(axis='y', labelsize=28)
            ax.set_ylabel('height[mm]', fontsize = 36)
            if zero:
                ax.set_xlabel('x[mm]', fontsize = 36)
            else:
                ax.set_xlabel('Wavelength[nm]', fontsize = 36)
            # change the cbar 
            cbar = fig.colorbar(self.cont,format='%.d')
            cbar.ax.tick_params(labelsize = 22)
            cbar.set_label('Counts', fontsize = 36)
            return ax 
        
        #call the animation function
        ani = animation.FuncAnimation(fig, animate, self.Frame)
        # if save, save the animation to the folder
        if save:
            if fname == None:
                fname = self.filename[:-3]+'gif'
            ani.save(fname, writer='Pillow',dpi=dpi)   
        return
        
        