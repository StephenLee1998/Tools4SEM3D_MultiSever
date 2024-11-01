from obspy import Trace, UTCDateTime
from scipy.interpolate import interp1d
import numpy as np

### Bandpass the waveforms with obspy
def bpfilter(data, fs, trange):
    tr = Trace(data=data)
    tr.stats.starttime = UTCDateTime()  # Assuming the current time for simplicity
    tr.stats.sampling_rate = fs
    tr.filter('bandpass', freqmin=1/trange[1], freqmax=1/trange[0], corners=4, zerophase=True)
    return tr.data

def interpolate_data(t_old, data_old, t_new):
    #interpolate the data in new time series
    #t_old: raw time series 
    #data_old: raw data series
    #t_new: new data series
    
    interpolation_func = interp1d(t_old, data_old, bounds_error=False, fill_value="extrapolate")
    data_new = interpolation_func(t_new)
    return data_new

### taper the waveforms part from Pyadjoint
def window_taper(signal, taper_percentage, taper_type):
    """
    window taper function.

    :param signal: time series
    :type signal: ndarray(float)

    :param taper_percentage: total percentage of taper in decimal
    :type taper_percentage: float

    return : tapered input ndarray

    taper_type:
    1, cos
    2, cos_p10
    3, hann
    4, hamming

    To do:
    with options of more tapers
    """
    taper_collection = ('cos', 'cos_p10', 'hann', "hamming")

    if taper_type not in taper_collection:
        raise ValueError("Window taper not supported")

    if taper_percentage < 0 or taper_percentage > 1:
        raise ValueError("Wrong taper percentage")

    npts = len(signal)

    if taper_percentage == 0.0 or taper_percentage == 1.0:
        frac = int(npts*taper_percentage / 2.0)
    else:
        frac = int(npts*taper_percentage / 2.0 + 0.5)

    idx1 = frac
    idx2 = npts - frac

    if taper_type == 'hann':
        signal[:idx1] *=\
            (0.5 - 0.5 * np.cos(2.0 * np.pi * np.arange(0, frac) /
                                (2 * frac - 1)))
        signal[idx2:] *=\
            (0.5 - 0.5 * np.cos(2.0 * np.pi * np.arange(frac, 2 * frac) /
                                (2 * frac - 1)))

    if taper_type == 'hamming':
        signal[:idx1] *=\
            (0.54 - 0.46 * np.cos(2.0 * np.pi * np.arange(0, frac) /
                                  (2 * frac - 1)))
        signal[idx2:] *=\
            (0.54 - 0.46 * np.cos(2.0 * np.pi * np.arange(frac, 2 * frac) /
                                  (2 * frac - 1)))

    if taper_type == 'cos':
        power = 1.
        signal[:idx1] *= np.cos(np.pi * np.arange(0, frac) /
                                (2 * frac - 1) - np.pi / 2.0) ** power
        signal[idx2:] *= np.cos(np.pi * np.arange(frac, 2 * frac) /
                                (2 * frac - 1) - np.pi / 2.0) ** power

    if taper_type == 'cos_p10':
        power = 10.
        signal[:idx1] *= 1. - np.cos(np.pi * np.arange(0, frac) /
                                     (2 * frac - 1)) ** power
        signal[idx2:] *= 1. - np.cos(np.pi * np.arange(frac, 2 * frac) /
                                     (2 * frac - 1)) ** power
    return signal