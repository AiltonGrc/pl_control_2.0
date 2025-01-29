import os.path
import time

import lmfit
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

from .spe_reader import openXMLline,getData,Spectra

from lmfit.models import SineModel, GaussianModel, ConstantModel,LorentzianModel,LinearModel
from lmfit import report_fit


def nm_to_ev(wl):
    return 1239.8/wl

def fwhm_ev(wl,fwhm):
    return nm_to_ev(wl-fwhm)-nm_to_ev(wl+fwhm)

def resid(params, x,model, data=None):

    if model == 'cossinus':

        parvals = params.valuesdict()

        frequency = parvals['frequency']
        slope = parvals['slope']
        intercept = parvals['intercept']
        amplitude = parvals['amplitude']
        shift = parvals['shift']

        model = amplitude * np.cos(x * 2 * np.pi * frequency + shift) + intercept

    if data != None: weight = [1 / _ for _ in data]

    if data == None:
        return model
    else: return model - data

def fit_higher_peak(root = None, filename = None, flag_save = None,flag_plot=False):

    result_dict = {"sample_name"            :'empty',
                   "center_point"           : 5,
                   "center"                 : [],
                   "fwhm"                   : [],
                   "fwhm_ev"                : [],
                   "fwhm_reduced_csi_square": [],
                   "FSS"                    : [],
                   "Polarization"           : [],
                   "Phase"                  : [],
                   "second_x_axis"          : [],
                   "df"                     : []
                   }

    if root == "test" or (root == None and filename in (None,'test')):
        root = "I:\\public\\NANOSCALE SEMICONDUCTOR GROUP\\1. DATA\\SMALL-LAB\\2022\\11\\24\\SA0484"
        file = "Sa0484_QD1_30µmslit_150grating_White Light_FEL600_1s_10K_#Pol(deg)_0.0_360.0.spe"
        filename = os.path.join(root, file)
        flag_plot = True

    else:
        if root == None:
            file = filename.split('\\')[-1]
        else:
            file = filename
            filename = os.path.join(root, filename)

    if flag_save == None: flag_save = False

    df = Spectra(os.path.join(filename))
    result_dict['df'] = df

    settings, measurement = file.replace(".spe", "").split("#")
    try:
        measurement_start, measurement_end = float(measurement.split("_")[1]), float(measurement.split("_")[2])
    except Exception as e:
        print("Exception: ", e)
        print("Could not process file")
        return None, None

    result_dict['second_x_axis'] = np.arange(measurement_start, measurement_end, (measurement_end - measurement_start) / len(df.intensity))

    result_dict['sample_name'] = settings.split("_")[:6]

    for _ in np.arange(len(df.intensity)):

        maximum = max(df.intensity[0])
        center_point = np.where(df.intensity[0] == maximum)[0][0]
        result_dict['center_point'] = center_point
        result_dict['peak_x_coord'] = df.wavelength[center_point - 15:center_point + 15]
        center = df.wavelength[center_point]

        peak_lorentz = LorentzianModel() + ConstantModel()

        result_lorentz = peak_lorentz.fit(df.intensity[_][center_point - 15:center_point + 15],
                                          x=df.wavelength[center_point - 15:center_point + 15] \
                                          , center=df.wavelength[center_point],
                                          amplitude=df.intensity[_][center_point] * 1.5)

        peak_gauss = LorentzianModel() + ConstantModel()

        result_gauss = peak_gauss.fit(df.intensity[_][center_point - 15:center_point + 15],
                                      x=df.wavelength[center_point - 15:center_point + 15] \
                                      , center=df.wavelength[center_point],
                                      amplitude=df.intensity[_][center_point] * 1.5)

        if result_gauss.redchi < result_lorentz.redchi:
            result = result_gauss
            # print("gauss")
        else:
            result = result_lorentz
            # print("lorentz")

        result_dict['center'].append(result.values['center'])
        result_dict['fwhm'].append(result.values['fwhm'] / 2)
        result_dict['fwhm_reduced_csi_square'].append(result.redchi)

    result_dict['center_wl_average'] = np.average(result_dict['center'], weights=result_dict['fwhm_reduced_csi_square'])
    result_dict['FWHM'] = np.average(result_dict['fwhm'], weights=result_dict['fwhm_reduced_csi_square'])
    result_dict['FWHM_e'] = fwhm_ev(result_dict['center_wl_average'], result_dict['fwhm'])

    result_dict['fwhm_ev'] = [fwhm_ev(result_dict['center_wl_average'], _) for _ in result_dict['fwhm']]
    result_dict['center_ev'] = [nm_to_ev(_) for _ in result_dict['center']]
    result_dict['exciton_fit'] = result.best_fit


    return result_dict,result

def plot_graph(result_dict,fit_result):

        text = text = lmfit.fit_report(fit_result)

        spectra = int(np.round(len(result_dict['df'].intensity)/2,0))

        fig, (ax0,ax1) = plt.subplots(2,1,sharex=False,sharey=False)
        ax0.plot(result_dict['df'].wavelength, result_dict['df'].intensity[spectra])
        ax0.plot(result_dict['df'].wavelength[result_dict['center_point'] - 15:result_dict['center_point'] + 15]
                 ,fit_result.best_fit,'--r')
        ax1.plot(result_dict['second_x_axis'],result_dict['center_ev'],label='data')
        # ax1.plot(angles_rad,resid(result_fss.init_vals,angles),'-g',label='start')
        ax1.plot(result_dict['second_x_axis'], result_dict['fit_data'],'--r',label='best fit')
        ax1.text(.1,.4,text ,transform = ax1.transAxes)
        plt.legend()
        fig.canvas.draw_idle()
        plt.show()
        plt.pause(.1)

def Fss_fit(root = None, filename = None, flag_plot=False):

    result_dict, result_exciton = fit_higher_peak(root, filename, flag_plot)

    amplitude = abs(max(result_dict['center_ev']) - min(result_dict['center_ev']))/2

    period = abs(result_dict['center'].index(max(result_dict['center'])) - \
                    result_dict['center'].index(min(result_dict['center'])))
    period *= (result_dict['second_x_axis'][-1] - result_dict['second_x_axis'][0]) / len(result_dict['df'].intensity)

    phase = max(result_dict['center_ev'], key=lambda x:\
        abs(x - np.average(result_dict['center_ev'])))
    # phase *= (angles[-1] - angles[0]) / len(df.intensity)

    # print(frequency,Amplitude,phase)

    intercept = np.average(result_dict['center_ev'])
    slope = abs(result_dict['center_ev'][-1]-result_dict['center_ev'][0])\
            /(abs(result_dict['second_x_axis'][-1] - result_dict['second_x_axis'][0]))

    params = lmfit.Parameters()

    # fss_model.guess(result_dict['center_ev'], x=angles_rad)
    params.add('amplitude',  value=amplitude,min=amplitude*.75, max=amplitude)
    params.add('frequency', value=1/180,min=0,max=270)
    params.add('shift', value=phase,min=0,max=np.deg2rad(180))
    params.add('intercept', value=intercept,min=1e-8,max=2)
    params.add('slope',value=slope,min=-1e-8,max=1e-8)

    result_dict['model'] = 'cossinus'

    result_fss_minimizer = lmfit.Minimizer(resid,params,fcn_args=(result_dict['second_x_axis'],result_dict['model']
                                                   ,result_dict['center_ev']))

    result_fss = result_fss_minimizer.minimize()

    result_dict['fit_data'] = result_dict['center_ev'] - result_fss.residual

    result_dict['FSS'].append(2*result_fss.params['amplitude'].value*1e6)
    result_dict['Polarization'].append(np.rad2deg(result_fss.params['frequency'].value))
    result_dict['Phase'].append(np.rad2deg(result_fss.params['shift'].value))
    result_dict['plot_text'] = lmfit.fit_report(result_fss)

    # print(result_dict['fwhm_ev'])
    try: result_dict['dev_center'] = np.std(result_dict['center'])
    except: result_dict['dev_center'] = 0
    try: result_dict['dev_fwhm_ev'] = np.std(result_dict['fwhm_ev'])
    except: result_dict['dev_center'] = 0
    try: result_dict['dev_FSS'] = result_fss.params['amplitude'].stderr
    except: result_dict['dev_FSS'] = 0
    try: result_dict['dev_pol'] = np.rad2deg(result_fss.params['frequency'].stderr)
    except: result_dict['dev_pol'] = 0
    try: result_dict['dev_phase'] = np.rad2deg(result_fss.params['shift'].stderr)
    except: result_dict['dev_phase'] = 0

    if flag_plot:
        plot_graph(result_dict,result_exciton)

    return result_dict

def Pow_fit(root = None, filename = None, flag_plot=False):

    result_dict, result_exciton = fit_higher_peak(root, filename,flag_plot)

    return result_dict

def IV_fit(root = None, filename = None, flag_plot=False):

    result_dict, result_exciton = fit_higher_peak(root, filename, flag_plot)

    return result_dict