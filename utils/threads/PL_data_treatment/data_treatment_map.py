import os.path
import time

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.widgets import Slider, Button
import matplotlib as mpl
from scipy.interpolate import griddata

from spe_reader import openXMLline,getData,Spectra

from lmfit.models import SineModel, GaussianModel, ConstantModel,LorentzianModel, QuadraticModel, LinearModel
from lmfit import report_fit

from matplotlib import use
from matplotlib.widgets import Slider

use('TkAgg')

import json

def nm_to_ev(wl):
    if wl == 0:return 0
    return 1239.8/wl

def fwhm_ev(wl,fwhm):
    return nm_to_ev(wl-fwhm)-nm_to_ev(wl+fwhm)

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


def Gauss_fit_background(root = None, filename = None, flag_save = True,show_plot=False,show_every=100):

    points_to_remove = 40

    result_dict = {"center": [],
              "fwhm":[],
              "fwhm_ev":[],
              "fwhm_reduced_csi_square":[],
              "intensity": [],
              "wl":[],
              "fit":[]
    }

    if root == "test":
        # root = r"I:\public\NANOSCALE SEMICONDUCTOR GROUP\1. DATA\SMALL-LAB\2022\12\10\SA0952"
        root = Path.cwd()#r"C:\Users\Ailton Garcia\My Drive\Desktop Notebook\PL_data_treatment"
        file = "SA0952_10xobj_200µmslit_1200grating_GreenLaser_FEL600_1s_10K_765nm_#spacemap#(µm)_80_500.spe"
        filename = os.path.join(root,file)
        flag_plot = True
        if flag_save == None: flag_save = False

    if flag_save == None: flag_save = False

    df = Spectra(filename)


    settings, measurement,meas_settings = file.replace(".spe","").split("#")
    measurement_start, measurement_end = float(meas_settings.split("_")[1]),float(meas_settings.split("_")[2])

    angles = np.arange(measurement_start, measurement_end, (measurement_end - measurement_start) / len(df.intensity))
    angles_rad = [np.deg2rad(_) for _ in angles]

    if len(settings.split("_"))==6:
        sample_name,slit,grating,source,filter,acquisition,temperature,wl,(_) = settings.split("_")
    elif len(settings.split("_"))>6:
        sample_name,comment,slit,grating,source,filter,acquisition,temperature,wl,(_) = settings.split("_")

    if os.path.isfile(filename.replace(".spe",".json")):
        try:
            with open(filename.replace(".spe",".json"), 'r') as f:
                result_dict = json.load(f)

            last_spec = len(result_dict['intensity'])
        except Exception as e:
            print("ERROR in file load: {}".format(e))
            last_spec=0
    else:last_spec=0

    for spec in np.arange(last_spec,len(df.intensity),1):

        print(spec)

        polyn_fit = QuadraticModel()

        intensity = df.intensity[spec][points_to_remove:]

        for _ in np.arange(len(intensity)):
            intensity[_] -= 600

        wl = df.wavelength[points_to_remove:]

        result_dict['wl'] = wl

        back_fit = polyn_fit.fit(intensity,x=result_dict['wl'])

        intensity = intensity - back_fit.best_fit
        intensity = intensity.tolist()

        result_dict['intensity'].append(intensity)

        # plt.figure(1)
        # plt.plot(result_dict['wl'],result_dict['intensity'][spec])
        # plt.show()
        # plt.pause(.1)

        center_point = np.where(np.array(intensity)==max(intensity))[0][0]

        center = result_dict['wl'][center_point]

        if max(intensity)-abs(min(intensity)) < 42:
            lin = LinearModel()

            result = lin.fit(result_dict['intensity'][spec], \
                             x=result_dict['wl'])

            intensity = result.best_fit

            intensity = intensity.tolist()

            result_dict['fit'].append(intensity)
            result_dict['center'].append(int(0))
            result_dict['fwhm'].append(int(500))
            result_dict['fwhm_reduced_csi_square'].append(int(5000))
            result.values['amplitude'] = 0

            # plt.figure(1)
            # plt.plot(result_dict['wl'], result_dict['intensity'][spec])
            # plt.plot(result_dict['wl'], result.best_fit)
            # plt.show()
            # plt.pause(.1)
            #
            # time.sleep(1)
            #
            # plt.close()

            continue

        intensity_center_fit = 1

        try:

            peak_lorentz = LorentzianModel()+ ConstantModel()

            result_lorentz = peak_lorentz.fit(result_dict['intensity'][spec],x=result_dict['wl']\
                              ,center = result_dict['wl'][center_point],amplitude=result_dict['intensity'][spec][center_point]*1.5,max_nfev=10)

            peak_gauss = LorentzianModel() + ConstantModel()

            result_gauss = peak_gauss.fit(result_dict['intensity'][spec],
                                      x=result_dict['wl'],\
                                      center=result_dict['wl'][center_point], amplitude=result_dict['intensity'][spec][center_point] * 1.5,max_nfev=10)

            if result_gauss.redchi<result_lorentz.redchi:
                result = result_gauss
                # print("gauss")
            else:
                result = result_lorentz
                # print("lorentz")

        except:
            lin = LinearModel()

            result = lin.fit(result_dict['intensity'][spec],\
                                 x=result_dict['wl'])

            intensity = result.best_fit

            intensity = intensity.tolist()

            result_dict['fit'].append(intensity)
            result_dict['center'].append(int(0))
            result_dict['fwhm'].append(int(500))
            result_dict['fwhm_reduced_csi_square'].append(int(5000))
            result.values['amplitude'] = 0

        if result.values['fwhm'] / 2 < 0.002 or intensity_center_fit<0:


            lin = LinearModel()

            result = lin.fit(result_dict['intensity'][spec], \
                             x=result_dict['wl'])
            intensity = result.best_fit

            intensity = intensity.tolist()

            result_dict['fit'].append(intensity)
            result_dict['center'].append(int(0))
            result_dict['fwhm'].append(int(0))
            result_dict['fwhm_reduced_csi_square'].append(int(5000))
            result.values['amplitude'] = 0

        else:
            intensity = result.best_fit

            intensity = intensity.tolist()

            result_dict['fit'].append(intensity)
            result_dict['center'].append(result.values['center'])
            result_dict['fwhm'].append(result.values['fwhm'] / 2)
            result_dict['fwhm_reduced_csi_square'].append(result.redchi)

        if show_plot:
            if spec%show_every == 0:
                plt.figure(1)
                plt.plot(result_dict['wl'],result_dict['intensity'][spec])
                plt.plot(result_dict['wl'], result.best_fit)
                plt.draw()
                plt.pause(.1)
                time.sleep(1)
                plt.close('all')

        try:
            with open(filename.replace(".spe", ".json"), 'w+') as f:
                json.dump(result_dict, f,indent=4)
        except Exception as e: print("Error with json dump: {}".format(e))

    result_dict['center_wl_average'] = np.average(result_dict['center'], weights=result_dict['fwhm_reduced_csi_square'])
    result_dict['FWHM'] = np.average(result_dict['fwhm'], weights=result_dict['fwhm_reduced_csi_square'])
    result_dict['FWHM_e'] = fwhm_ev(result_dict['center_wl_average'],result_dict['FWHM'])

    result_dict['fwhm_ev'] = [fwhm_ev(result_dict['center_wl_average'],_) for _ in result_dict['fwhm']]
    result_dict['center_ev'] = [nm_to_ev(_) for _ in result_dict['center']]

    result_dict['dev_center'] = np.std(result_dict['center'])
    result_dict['dev_fwhm_ev'] = np.std(result_dict['fwhm_ev'])

    # best_fit = np.where(result_dict['fwhm_reduced_csi_square'] == min(result_dict['fwhm_reduced_csi_square']))[0][0]

    axs0_text = "FWHM: {0:2.2f} \u00B1 {1:2.2f} \u03BCeV\nCenter wavelentgh: {2:2.2f} \u00B1 {3:2.1E} nm" \
                "".format(result_dict['FWHM_e']*1e6,result_dict['dev_fwhm_ev']*1e6,result_dict['center_wl_average'],result_dict['dev_center'])


    fig, axs = plt.subplots(2,figsize=(8,8))
    fig.suptitle('Space Map Measurement')

    # Make a horizontal slider to control the frequency.
    axfreq = fig.add_axes([0.25, 0.1, 0.65, 0.03])
    freq_slider = Slider(
        ax=axfreq,
        label='Spec #',
        valmin=0,
        valmax=len(df.intensity),
        valinit=0,
    )

    line_fit, = axs[0].plot(result_dict['wl'], result_dict['fit'][0])
    line_data, = axs[0].plot(result_dict['wl'], df.intensity[0][points_to_remove:])

    axs[0].annotate(axs0_text,xy=(.50,.85),xycoords="axes fraction")
    axs[0].set_xlabel("Wavelength (nm)")
    axs[0].set_ylabel("Intensity (arb. un.)")
    filename = filename.replace(".spe",".png")

    print(filename)

    def update(val):
        line_fit, = axs[0].plot(result_dict['wl'], result_dict['fit'][spec])
        line_data, = axs[0].plot(result_dict['wl'], df.intensity[val])
        fig.canvas.draw_idle()

    if flag_plot:fig.show()
    if flag_save:
        print("savign figure at", filename)
        fig.savefig(filename.replace(".spe",".png"),dpi=300)

    return df, result_dict,axs0_text

def plot_data(filename=None):

    if filename == None:
        filename = "SA0952_10xobj_200µmslit_1200grating_GreenLaser_FEL600_1s_10K_765nm_#spacemap#(µm)_80_500_latest.spe"

    with open(filename.replace(".spe", ".json"), 'r') as f:
        result_dict = json.load(f)

    settings, measurement, meas_settings = filename.replace(".spe", "").split("#")
    column_number, step = float(meas_settings.split("_")[1]), float(meas_settings.split("_")[2])

    sample_name, comment, slit, grating, source, filter, acquisition, temperature, wl, (_) = settings.split("_")

    lines_numb = len(result_dict['center'])/column_number

    map = [[]]

    coord = []

    print(lines_numb, column_number, step)


    xg = np.arange(0,int((column_number*step)),step)
    yg = np.arange(0,int((lines_numb*step)),step)

    X,Y = np.meshgrid(xg,yg)

    Z = []

    # zg = griddata(result_dict['center'], tempratio, (Xg, Yg), method='linear')

    j=0

    for line in np.arange(0,int(lines_numb),1):

        for column in np.arange(0,int(column_number),1):

            coord.append(((column)*step,(line)*step,result_dict['center'][j]))
            Z.append(result_dict['center'][j])

            print("line,column,spectra,center\n",line,column,j,result_dict['center'][j])

            map[line].append(result_dict['center'][j])

            j+=1

        map.append([])

    Z = np.array(Z).reshape(52,80)

    fig, ax = plt.subplots()

    levels = np.arange(755,775,1)

    cmap = mpl.cm.get_cmap("plasma").copy()

    # cs = ax.contourf(X, Y, Z, levels,cmap=cmap)

    cs = ax.pcolormesh(X, Y, Z, cmap=cmap,vmin=755,vmax=775)

    cs.cmap.set_under('black')
    cs.cmap.set_over('white')

    fig.colorbar(cs)

    plt.show()

def plot_single(filename=None):

    if filename == None:
        filename = "SA0952_10xobj_200µmslit_1200grating_GreenLaser_FEL600_1s_10K_765nm_#spacemap#(µm)_80_500_latest.spe"

    with open(filename.replace(".spe", ".json"), 'r') as f:
        result_dict = json.load(f)

    settings, measurement, meas_settings = filename.replace(".spe", "").split("#")
    column_number, step = float(meas_settings.split("_")[1]), float(meas_settings.split("_")[2])

    sample_name, comment, slit, grating, source, filter, acquisition, temperature, wl, (_) = settings.split("_")

    lines_numb = len(result_dict['center'])/column_number

    map = [[]]

    coord  = []

    print(lines_numb, column_number, step)

    fig, ax = plt.subplots()
    line_fit, = ax.plot(result_dict['wl'], result_dict['fit'][0])
    line_data, = ax.plot(result_dict['wl'], result_dict['intensity'][0])
    ax.set_xlim([result_dict['wl'][0], result_dict['wl'][-1]])
    ax.set_ylim([min(result_dict['intensity'][0]), max(result_dict['intensity'][0])])
    ax_text = fig.text(.01,.90,"Center wavelength (nm):\n {}".format(result_dict['center'][0]))
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Intensity (a.u.)")
    plt.tight_layout()

    # adjust the main plot to make room for the sliders
    fig.subplots_adjust(left=0.25, bottom=0.25)

    # Make a horizontal slider to control the frequency.
    axfreq = fig.add_axes([0.25, 0.1, 0.65, 0.03])
    freq_slider = Slider(
        ax=axfreq,
        label='Spectra',
        valmin=1,
        valmax=int(lines_numb*column_number),
        valinit=1,
    )

    def update(val):
        val = int(val)-1
        line_fit.set_ydata(result_dict['fit'][val])
        line_data.set_ydata(result_dict['intensity'][val])
        ax.set_ylim([min(result_dict['intensity'][val]), max(result_dict['intensity'][val])])
        ax_text.set_text("Center wavelength (nm):\n {}".format(result_dict['center'][val]))
        fig.canvas.draw_idle()

    freq_slider.on_changed(update)

    plt.show()

# Gauss_fit_background("test")
plot_single()
# plot_data()