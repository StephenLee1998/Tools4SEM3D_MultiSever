import numpy as np
from scipy.integrate import simps
from Utils.Signals import window_taper

def _xcorr_shift(d, s):
    cc = np.correlate(d, s, mode="full")
    time_shift = cc.argmax() - len(d) + 1
    return time_shift


def cc_error(d1, d2, deltat, cc_shift, cc_dlna, sigma_dt_min, sigma_dlna_min):
    """
    Estimate error for dt and dlna with uncorrelation assumption
    """
    nlen_t = len(d1)

    d2_cc_dt = np.zeros(nlen_t)
    d2_cc_dtdlna = np.zeros(nlen_t)

    for index in range(0, nlen_t):
        index_shift = index - cc_shift

        if 0 <= index_shift < nlen_t:
            # corrected by c.c. shift
            d2_cc_dt[index] = d2[index_shift]

            # corrected by c.c. shift and amplitude
            d2_cc_dtdlna[index] = np.exp(cc_dlna) * d2[index_shift]

    # time derivative of d2_cc (velocity)
    d2_cc_vel = np.gradient(d2_cc_dtdlna, deltat)

    # the estimated error for dt and dlna with uncorrelation assumption
    sigma_dt_top = np.sum((d1 - d2_cc_dtdlna)**2)
    sigma_dt_bot = np.sum(d2_cc_vel**2)

    sigma_dlna_top = sigma_dt_top
    sigma_dlna_bot = np.sum(d2_cc_dt**2)

    sigma_dt = np.sqrt(sigma_dt_top / sigma_dt_bot)
    sigma_dlna = np.sqrt(sigma_dlna_top / sigma_dlna_bot)

    if sigma_dt < sigma_dt_min:
        sigma_dt = sigma_dt_min

    if sigma_dlna < sigma_dlna_min:
        sigma_dlna = sigma_dlna_min

    return sigma_dt, sigma_dlna

def cc_traveltime_misfit(observed, synthetic, deltat, window, config):  # NOQA
    ret_val_p = {}
    ret_val_q = {}

    nlen_data = len(synthetic)

    fp = np.zeros(nlen_data)
    fq = np.zeros(nlen_data)

    misfit_sum_p = 0.0
    misfit_sum_q = 0.0
    deltat = float(deltat)

    for wins in window:
        left_window_border = wins[0]
        right_window_border = wins[1]

        left_sample = int(np.floor(left_window_border / deltat)) + 1
        if left_sample < 0:
            left_sample = 0
        nlen = int(np.floor((right_window_border -
                             left_window_border) / deltat)) + 1
        right_sample = left_sample + nlen

        d = np.zeros(nlen)
        s = np.zeros(nlen)

        d[0:nlen] = observed[left_sample:right_sample]
        s[0:nlen] = synthetic[left_sample:right_sample]

        # All adjoint sources will need some kind of windowing taper
        d = window_taper(d, taper_percentage=np.float(config.taper_percentage),
                     taper_type=config.taper_type)
        s = window_taper(s, taper_percentage=np.float(config.taper_percentage),
                     taper_type=config.taper_type)

        i_shift = _xcorr_shift(d, s)
        t_shift = i_shift * deltat

        cc_dlna = 0.5 * np.log(sum(d[0:nlen]*d[0:nlen]) /
                               sum(s[0:nlen]*s[0:nlen]))

        sigma_dt, sigma_dlna = cc_error(d, s, deltat, i_shift, cc_dlna,
                                        np.float(config.dt_sigma_min),
                                        np.float(config.dlna_sigma_min))

        misfit_sum_p += 0.5 * (t_shift/sigma_dt) ** 2
        misfit_sum_q += 0.5 * (cc_dlna/sigma_dlna) ** 2

        dsdt = np.gradient(s, deltat)
        nnorm = simps(y=dsdt*dsdt, dx=deltat)
        fp[left_sample:right_sample] = dsdt[:] * t_shift / nnorm / sigma_dt**2

        mnorm = simps(y=s*s, dx=deltat)
        fq[left_sample:right_sample] =\
            -1.0 * s[:] * cc_dlna / mnorm / sigma_dlna ** 2

    ret_val_p["misfit"] = misfit_sum_p
    ret_val_q["misfit"] = misfit_sum_q

    
    ret_val_p["adjoint_source"] = fp[::-1]
    ret_val_q["adjoint_source"] = fq[::-1]

    if config.measure_type == "dt":
        return ret_val_p, t_shift

    if config.measure_type == "am":
        return ret_val_q, t_shift
    

    
