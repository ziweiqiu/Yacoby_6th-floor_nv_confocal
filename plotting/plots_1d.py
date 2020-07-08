"""
    This file is part of b26_toolkit, a pylabcontrol add-on for experiments in Harvard LISE B26.
    Copyright (C) <2016>  Arthur Safira, Jan Gieseler, Aaron Kabcenell

    b26_toolkit is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    b26_toolkit is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with b26_toolkit.  If not, see <http://www.gnu.org/licenses/>.
"""

import matplotlib.patches as patches
import numpy as np
from matplotlib.collections import PatchCollection
from b26_toolkit.data_processing.fit_functions import lorentzian, double_lorentzian
import math

def plot_psd(freq, psd, axes, y_scaling = 'log', x_scaling = 'lin'):
    '''
    plots the power spectral density on to the canvas axes
    :param freq: x-values array of length N
    :param psd: y-values array of length N
    :param axes: target axes object
    :return: None
    '''
    unit = 'Hz'
    c_unit = 1.0
    if x_scaling == 'lin':
        if np.mean(freq) > 1e6:
            c_unit = 1e-6
            unit = 'MHz'
        elif np.mean(freq) > 1e3:
            c_unit = 1e-3
            unit = 'kHz'


    if y_scaling == 'log' and x_scaling == 'log':
        axes.loglog(c_unit*freq, psd, 'b')
    elif y_scaling == 'log' and x_scaling == 'lin':
        axes.semilogy(c_unit*freq, psd, 'b')
    elif y_scaling == 'lin' and x_scaling == 'log':
        axes.semilogx(c_unit*freq, psd, 'b')
    elif y_scaling == 'lin' and x_scaling == 'lin':
        axes.plot(c_unit*freq, psd, 'b')

    axes.set_xlabel('frequency ({:s})'.format(unit))

    axes.set_xlim([min(c_unit*freq), max(c_unit*freq)])

def plot_esr(axes, frequency, counts, fit_params=None, plot_marker_data = 'b', plot_marker_fit = 'r', linestyle = '-', marker = '.', avg_counts = 0, mw_power = 0, D = 2.87E9, gama = 2.8028E6, err = None):
    """
    plots the esr
    Args:
        axes: axes object
        fit_params: array with fitparameters either length 4 for single peak or length 6 for double peak
        frequency: mw frequency (array)
        counts: counts (array)
        plot_marker:  (str) plot marker
        D: zero-field splitting ~ 2.87 GHz
        gama: gyromagnetic ratio ~ 2.8 MHz/Gauss

    Returns:

    """
    #  ======== plot data =========
    def B_field(f1,f2,D,gama):
        # this is to calculate the magnitude of the magnetic field and its angle to the NV axis
        B_parallel_sqrd = -(D + f2 - 2 * f1) * (D + f1 - 2 * f2) * (D + f1 + f2)
        B_perpendicular_sqrd = -(2 * D - f2 - f1) * (2 * D + 2* f1 - f2) * (2 * D - f1 + 2 * f2)

        if B_parallel_sqrd < 0 or B_perpendicular_sqrd < 0:
            print('check zero-field splitting')
            return -1, -1
        else:
            B_parallel = np.sqrt(B_parallel_sqrd) / (3 * gama * np.sqrt(3 * D))
            B_perpendicular = np.sqrt(B_perpendicular_sqrd) / (3 * gama * np.sqrt(3 * D))
            B_mag = np.sqrt(B_parallel * B_parallel + B_perpendicular * B_perpendicular)
            angle = np.arctan(B_perpendicular / B_parallel) / math.pi * 180
            return B_mag, angle

    axes.clear() # ER 20181012 - matplotlib axes.hold() removed in update to 3.0.0

    if len(counts)<=len(frequency):
        len1 = len(counts)
    else:
        len1 = len(frequency)

    if err is not None:
        counts_err = err
    else:
        counts_err = np.zeros(int(len1))

    axes.errorbar(frequency[0:len1], counts[0:len1], counts_err[0:len1],linestyle=linestyle, marker=marker)

    title = 'ESR'
    fit_data = None

    #  ======== plot fit =========
    if fit_params is not None and len(fit_params) and fit_params[0] != -1:  # check if fit valid
        if len(fit_params) == 4: # single peak
            fit_data = lorentzian(frequency, *fit_params)
            title = 'ESR, averaged counts: {:0.2f}kcps \n fo = {:0.4f} GHz, wo = {:0.2f} MHz, a0 = {:.2f} % \n microwave power: {:.2f}dBm'.format(avg_counts, fit_params[2]/1E9,
                                                                                             fit_params[3]/1E6, fit_params[1]*100, mw_power)
        elif len(fit_params) == 6: # double peak
            B_mag, angle = B_field(fit_params[4],fit_params[5],D,gama)
            fit_data = double_lorentzian(frequency, *fit_params)
            title = 'ESR, averaged counts: {:.2f}kcps  \n f1 = {:0.4f} GHz, f2 = {:0.4f} GHz, wo = {:0.2f} MHz \n a1 = {:.2f} %, a2 = {:.2f} %  \n microwave power: {:.2f}dBm, |B| ~ {:.2f} G, θ ~ {:.2f} deg'.format(
                avg_counts, fit_params[4]/1E9, fit_params[5]/1E9, fit_params[1]/1E6, fit_params[2]*100, fit_params[3]*100, mw_power, float(B_mag), float(angle))
    else:
        title = 'ESR, averaged counts: {:.2f}kcps \n microwave power: {:.2f}dBm'.format(avg_counts, mw_power)

    if fit_data is not None:
        axes.plot(frequency, fit_data)


    axes.set_title(title)
    axes.set_xlabel('Frequency (Hz)')
    axes.set_ylabel('Fluorescence')

    # return lines

def plot_magnet_sweep1D_ESR(axis, x, y1, y2, lbl, x_r = None, y1_r  = None , y2_r  = None):
    """
        plots ESR f0 and wo during magnet sweep

        Args:
            axis: the axis to draw the plot
            x (1d array): magnet position (forward)
            y1 (1d array): fo data (forward)
            y2 (1d array): wo data (forward)
            x_r (1d array): magnet position (backward)
            y1_r (1d array): fo data (backward)
            y2_r (1d array): wo data (backward)
            lbl (1d array): axis labels and title

        Returns: (none)
    """

    # Clear any previous plots on the axis
    axLeft = axis[0]
    axLeft.clear()
    axRight = axis[1]
    axRight.clear()

    # Plot the forward sweep
    if len(y1) < len(x):
        len1 = len(y1)
    else:
        len1 = len(x)

    if len(y2) < len(x):
        len2 = len(y2)
    else:
        len2 = len(x)

    axLeft.plot(x[0:len1], y1[0:len1], color = '#2379af', marker = '.')
    axRight.plot(x[0:len2], y2[0:len2], color = '#2379af', marker = '.')

    # Plot the backward sweep if it exists
    if x_r is not None and y1_r is not None:
        if len(y1_r) < len(x_r):
            len1_r = len(y1_r)
        else:
            len1_r = len(x_r)
        axLeft.plot(x_r[0:len1_r], y1_r[0:len1_r], color = '#d67b2c', marker='.')

    if x_r is not None and y2_r is not None:
        if len(y2_r) < len(x_r):
            len2_r = len(y2_r)
        else:
            len2_r = len(x_r)
        axRight.plot(x_r[0:len2_r], y2_r[0:len2_r], color = '#d67b2c', marker='.')

    if not (lbl == None):

        axLeft.set_xlabel(lbl[0])
        axLeft.set_ylabel(lbl[1])
        axLeft.set_title(lbl[3])
        axRight.set_xlabel(lbl[0])
        axRight.set_ylabel(lbl[2])
        axRight.set_title(lbl[3])

def plot_magnet_sweep1D_Fluor(axis, x, y1, lbl, x_r = None, y1_r  = None):
    """
       plots ESR f0 and wo during magnet sweep

        Args:
        axis: the axis to draw the plot
        x (1d array): magnet position
        y1 (1d array): fluorescence data

        lbl (1d array): axis labels

            Returns: (none)
    """

    # Clear any previous plots on the axis
    axLeft = axis[0]
    axLeft.clear()

    # Plot the forward sweep
    if len(y1)<len(x):
        len1 = len(y1)
    else:
        len1 = len(x)

    axLeft.plot(x[0:len1], y1[0:len1], color = '#2379af', marker = '.')

    # Plot the backward sweep if it exists
    if x_r is not None and y1_r is not None:
        if len(y1_r) < len(x_r):
            len1_r = len(y1_r)
        else:
            len1_r = len(x_r)

        axLeft.plot(x_r[0:len1_r], y1_r[0:len1_r], color='#d67b2c', marker='.')

    if not (lbl == None):
        axLeft.set_xlabel(lbl[0])
        axLeft.set_ylabel(lbl[1])
        axLeft.set_title(lbl[2])

def plot_pulses(axis, pulse_collection, pulse_colors=None, pulse_tag = None):
    """
    creates a visualization of pulses (in pulse_collection) on a matplotlib axis (axis)

    Args:
        axis: The axis for the matplotlib plot
        pulse_collection: a collection of pulses, named tuples (channel_id, start_time, duration)
        pulse_colors: a dictionary of {channel_id:matplotlib_color} that maps channels to colors

    Returns:

    """

    # create a list of unique instruments from the pulses
    instrument_names = sorted(list(set([pulse.channel_id for pulse in pulse_collection])))

    # assign colors for certain specific channels
    if pulse_colors is None:
        pulse_colors = {'laser': '#50FF00', 'microwave_i': 'r', 'apd_readout': '#fb0b03',
                        'microwave_switch': '#000000', 'microwave_switch_I': '#0000FF','microwave_switch_II': '#551A8B', 'apd_switch': '#feccca'}

    # find the maximum time from the list of pulses
    max_time = max([pulse.start_time + pulse.duration for pulse in pulse_collection])

    axis.clear()

    # set axis boundaries
    axis.set_ylim(-0.75, len(instrument_names) - .25)
    axis.set_xlim(0, max_time)

    # label y axis with pulse names
    axis.set_yticks(list(range(len(instrument_names))))
    axis.set_yticklabels(instrument_names)

    # create horizontal lines for each pulse
    for pulse_plot_y_position in range(0, len(instrument_names)):
        axis.axhline(pulse_plot_y_position - .25, 0.0, max_time, color='k')
    axis.tick_params(axis='y', which='both', length=0)  # remove tick marks on y axis

    # create a vertical line denoting the end of the pulse sequence loop
    # axis.axvline(max_time, -0.5, len(instrument_names), color='r')

    # create rectangles for the pulses
    patch_list = []
    for pulse in pulse_collection:
        patch_list.append(
            patches.Rectangle((pulse.start_time, instrument_names.index(pulse.channel_id) - .25), pulse.duration, 0.5,
                              fc=pulse_colors.get(pulse.channel_id, 'b')))

    patch_collection = PatchCollection(patch_list, match_original=True)
    axis.add_collection(patch_collection)

    # label the axis
    if pulse_tag is not None:
        axis.set_title('Pulse Visualization' + pulse_tag)
    else:
        axis.set_title('Pulse Visualization')
    axis.set_ylabel('pulse destination')
    axis.set_xlabel('time [ns]')

    xticks = np.array(axis.get_xticks())

    if np.sum(xticks > 1E9) > 0:
        axis.set_xticklabels(xticks/1e9)
        axis.set_xlabel('time [s]')

    elif np.sum(xticks > 1E6) > 0:
        axis.set_xticklabels(xticks/1e6)
        axis.set_xlabel('time [ms]')

    elif np.sum(xticks > 1E3) > 0:
        axis.set_xticklabels(xticks/1e3)
        axis.set_xlabel('time [us]')

def update_pulse_plot(axis, pulse_collection, pulse_colors=None, pulse_tag = None):
    """
    updates a previously created plot of pulses, removing the previous ones and adding ones corresponding to
    pulse_collection. The new pulse collection must only contain channel_ids already present on the passed axis

    Args:
        axis: The axis for the matplotlib plot
        pulse_collection: a collection of pulses, named tuples (channel_id, start_time, duration)
        pulse_colors: a dictionary of {channel_id:matplotlib_color} that maps channels to colors

    Returns:

    """

    # assign colors for certain specific channels
    if pulse_colors is None:
        # pulse_colors = {'laser': '#50FF00', 'microwave_i': 'r', 'apd_readout': 'k'}
        pulse_colors = {'laser': '#50FF00', 'microwave_i': 'r', 'apd_readout': '#fb0b03',
                        'microwave_switch': '#000000', 'microwave_switch_I': '#0000FF',
                        'microwave_switch_II': '#551A8B', 'apd_switch': '#feccca'}

    # get a list of unique instruments from the pulses
    instrument_names_old = [str(label.get_text()) for label in axis.get_yticklabels()]
    instrument_names = sorted(list(set([pulse.channel_id for pulse in pulse_collection])))
    if pulse_tag is not None:
        axis.set_title('Pulse Visualization' + pulse_tag)
    else:
        axis.set_title('Pulse Visualization')

    # if the number of pulses has changed call plot instead of update
    if set(instrument_names_old) != set(instrument_names):
        plot_pulses(axis, pulse_collection, pulse_colors = pulse_colors)
    else:
        # find the maximum time from the list of pulses
        max_time = max([pulse.start_time + pulse.duration for pulse in pulse_collection])

        axis.set_xlim(0, max_time)

        # remove the previous pulses
        [child.remove() for child in axis.get_children() if isinstance(child, PatchCollection)]

        # create rectangles for the pulses
        patch_list = []
        for pulse in pulse_collection:
            patch_list.append(
                patches.Rectangle((pulse.start_time, instrument_names.index(pulse.channel_id) - .25), pulse.duration, 0.5,
                                  fc=pulse_colors.get(pulse.channel_id, 'b')))

        patch_collection = PatchCollection(patch_list, match_original=True)
        axis.add_collection(patch_collection)


        xticks = np.array(axis.get_xticks())

        if np.sum(xticks > 1E9) > 0:
            axis.set_xticklabels(xticks/1e9)
            axis.set_xlabel('time [s]')

        elif np.sum(xticks > 1E6) > 0:
            axis.set_xticklabels(xticks/1e6)
            axis.set_xlabel('time [ms]')

        elif np.sum(xticks > 1E3) > 0:
            axis.set_xticklabels(xticks/1e3)
            axis.set_xlabel('time [us]')

def plot_counts(axis, data, axes_labels = None):
    """
    plots APD timeseries data

    Args:
        axis: the axis to draw the plot
        data (2d array): APD count timeseries data

    Returns: (none)
    """
    # print('start plot_counts')
    # print(data)

    # axis.plot(data[0],data[1], linewidth=2.0)
    axis.plot(data, linewidth=2.0)
    # axis.hold(False)
    # print('ok plot_counts')

    if axes_labels is None:
        x_label = 'time'

    else:
        x_label = axes_labels

    axis.set_xlabel(x_label)
    axis.set_ylabel('kCounts/sec')

def update_counts(axis, data):
    if data == None:
        return

    axis.lines[0].set_ydata(data)
    axis.lines[0].set_xdata(range(0,len(data)))
    axis.relim()
    axis.autoscale_view()

def plot_voltage(axis, data):
    """
    plots APD timeseries data

    Args:
        axis: the axis to draw the plot
        data (2d array): APD count timeseries data

    Returns: (none)
    """

    axis.plot(data, linewidth=2.0)
    # axis.hold(False)

    axis.set_xlabel('time')
    axis.set_ylabel('voltage (V)')

def plot_temperature(axis, data, sample_rate):
    """
    plots the temperature

    Args:
        axis:
        data:
        sample_rate: at which data has been acquired
    Returns:

    """

    time = np.arange(len(data))/float(sample_rate)

    label = 'time (s)'
    if max(time)>60:
        time /= 60.
        label = 'time (min)'
    if max(time)>60:
        time /= 60.
        label = 'time (h)'
    axis.plot(time, data)
    # axis.hold(False)

    axis.set_xlabel(label)
    axis.set_ylabel('temperature (K)')

def plot_1d_simple_timetrace_ns(axis, times, data_list, y_label='kCounts/sec', title=None, data_err_list = None):
    """
    plots a time trace for a list of data assuming that the times are give in ns
    Args:
        axis: axis object on which to plot
        times: times in ns (list or array of length N)
        data_list: list of data (size MxN)
        y_label: (optional) label for y axis
        title:  (optional) title

    """
    axis.clear()

    times = 1.*np.array(times) # cast onto numpy in case we got a list, which breaks some of the commands below

    x_label = 'time'
    if max(times) < 1e3:
        x_label = 'time [ns]'
    elif max(times) < 1e6:
        x_label = 'time [us]'
        times *= 1e-3
    elif max(times) < 1e9:
        x_label = 'time [ms]'
        times *= 1e-6
    elif max(times) < 1e12:
        x_label = 'time [s]'
        times *= 1e-9

    for counts in data_list:
        axis.plot(times, counts)
    # if data_err_list is None:
    #     for counts in data_list:
    #         axis.plot(times, counts)
    # else:
    #     i = 0
    #     for counts in data_list:
    #         counts_err = data_err_list[i]
    #         axis.errorbar(times, counts, yerr = counts_err)
    #         i += 1

    # axis.hold(False)


    axis.set_xlabel(x_label)
    axis.set_ylabel(y_label)
    if title:
        axis.set_title(title)
    axis.set_xlim([min(times), max(times)])

def update_1d_simple(axis, times, counts_list, counts_err_list = None):
    """
    Args:
        axis: axes object
        times: JG: THIS IS NOT USED! WHAT IS IT? => add comment, e.g. for future purpose or delete!
        counts_list: list of

    Returns:
    """

    if len(np.shape(counts_list)) == 1:
        counts_list = [counts_list]
    if len(axis.lines) != len(counts_list):
        counts_list = np.transpose(counts_list)


    #print('===>> JG 20181128: number of lines', len(axis.lines),' number of list ', len(counts_list), ' (they should be the same!!)')
    # if len(axis.lines) != len(counts_list):
    #     print('UUUUUU axes.lines:', len(axis.lines), 'len counts:', len(counts_list))

    #assert len(axis.lines) == len(counts_list)
    if len(axis.lines) != len(counts_list): # don't update the plot if the number of lines to plot isn't equal to the number of counts lists
        print('number of lines to plot is not equal to number of counts lists!!!')
        print('===>> ER 20181201: number of lines', len(axis.lines), ' number of list ', len(counts_list),
              ' (they should be the same!!)')
        print('counts_list: ', counts_list)
        print('axis.lines: ', axis.lines)
        return False

    else:
        for index, counts in enumerate(counts_list):
            axis.lines[index].set_ydata(counts)
        axis.relim()
        axis.autoscale_view()
        return True
        # if counts_err_list is None:
        #     for index, counts in enumerate(counts_list):
        #         axis.lines[index].set_ydata(counts)
        #     axis.relim()
        #     axis.autoscale_view()
        #     return True
        # else:
        #     pass
    # if len(axis.lines) != len(counts_list): # don't update the plot if the number of lines to plot isn't equal to the number of counts lists
    #     # print('number of lines to plot is not equal to number of counts lists!!!')
    #     # print('===>> ER 20181201: number of lines', len(axis.lines), ' number of list ', len(counts_list),
    #     #       ' (they should be the same!!)')
    #     # print('counts_list: ', counts_list)
    #     # print('axis.lines: ', axis.lines)
    #     print('WARNING: len(axis.lines) != len(counts_list) --> Updated plotting might be messed up')
    # for index, counts in enumerate(counts_list):
    #
    #     axis.lines[index].set_ydata(counts)
    # axis.relim()
    # axis.autoscale_view()

def plot_counts_vs_pos(axis, data, pos, x_label = None, y_label = None, title = None, marker = None):
    """
    plots magnet position vs. counts for aligning the field with fluorescence

    """
    axis.clear()
    if x_label is None:
        axis.set_xlabel('position (mm)')
    else:
        axis.set_xlabel(x_label)

    if y_label is None:
        axis.set_ylabel('kCounts/sec')
    else:
        axis.set_ylabel(y_label)

    if title:
        axis.set_title(title)

    if marker:
        axis.plot(pos, data, marker, linewidth=2.0, color='#2379af')
    else:
        axis.plot(pos, data, linewidth=2.0, color = '#2379af')

def update_counts_vs_pos(axis, data, pos):
    """

    updates the plot of the position vs counts when scanning the magnet

    """

    if data is None:
        return

    axis.lines[0].set_ydata(data)
    axis.lines[0].set_xdata(pos)
    axis.relim()
    #axis.set
    axis.autoscale_view()
