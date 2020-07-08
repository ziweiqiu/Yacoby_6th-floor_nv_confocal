# """
#     This file is part of b26_toolkit, a pylabcontrol add-on for experiments in Harvard LISE B26.
#     Copyright (C) <2016>  Arthur Safira, Jan Gieseler, Aaron Kabcenell
#
#     b26_toolkit is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     b26_toolkit is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with b26_toolkit.  If not, see <http://www.gnu.org/licenses/>.
# """

from pylabcontrol.core import Script, Parameter

# import standard libraries
import numpy as np
from b26_toolkit.instruments import AgilentMicrowaveGenerator, NI6353, LISE607RTPulseBlaster
# from b26_toolkit.scripts.spec_analyzer_get_spectrum import SpecAnalyzerGetSpectrum
from b26_toolkit.plotting.plots_1d import plot_esr

# from b26_toolkit.plotting.plots_1d import plot_diff_freq_vs_freq
from b26_toolkit.data_processing.esr_signal_processing import fit_esr
import time

class ESR_N9310A(Script):
    """

    This class runs ESR on an NV center, outputing microwaves using Agilent N9310A MicrowaveGenerator in the frequency sweep mode and reading in NV counts using a DAQ. This is relatively fast.
    This script will repeat ESR until a good fit is found, however the minimum and maximum repetition numbers are specified.
    Manual: https://scdn.rohde-schwarz.com/ur/pws/dl_downloads/dl_common_library/dl_manuals/gb_1/s/smb/SMB100A_OperatingManual_en_21.pdf

    --> written by Ziwei Qiu 3/18/2019

    """

    _DEFAULT_SETTINGS = [
        Parameter('power_out', -8.0, float, 'output power (dBm)'),
        Parameter('esr_avg_min', 30, int, 'minimum number of esr averages'),
        Parameter('esr_avg_max', 50, int, 'maximum number of esr averages'),
        Parameter('freq_start', 2.88e9, float, 'start frequency of scan [Hz]'),
        Parameter('freq_stop', 1e8, float, 'end frequency of scan [Hz]'),
        Parameter('range_type', 'center_range', ['start_stop', 'center_range'],
                  'start_stop: freq. range from freq_start to freq_stop. center_range: centered at freq_start and width freq_stop'),
        Parameter('freq_points', 100, int, 'number of frequencies in scan'),
        Parameter('reverse_sweep', False, bool, 'select to reverse the sweep direction'),
        Parameter('time_per_pt', 0.02, [0.01, 0.02, 0.04, 0.06, 0.08, 0.1,0.15,0.2,0.4,0.6,0.8,1], 'time per frequency point [s] (from 10ms to 10s)'),
        Parameter('turn_off_after', True, bool, 'if true MW output is turned off after the measurement'),
        Parameter('norm_to_ref', True, bool, 'If true normalize each frequency sweep by the average counts.'),
        Parameter('save_full_esr', True, bool, 'If true save all the esr traces individually'),
        Parameter('daq_type', 'PCI', ['PCI'], 'Type of daq to use for scan'),
        Parameter('fit_constants',
                  [
                      Parameter('num_of_peaks', -1, [-1, 1, 2],
                                'specify number of peaks for fitting. if not specifying the number of peaks, choose -1'),
                      Parameter('minimum_counts', 0.9, float, 'minumum counts for an ESR to not be considered noise'),
                      Parameter('contrast_factor', 3.0, float, 'minimum contrast for an ESR to not be considered noise'),
                      Parameter('zfs', 2.87E9, float, 'zero-field splitting [Hz]'),
                      Parameter('gama', 2.8028E6, float, 'NV spin gyromagnetic ratio [Hz/Gauss]'),
                  ]),
        Parameter('track_laser_power',
                  [
                      Parameter('on/off', False, bool,
                                'If true, measure and normalize out laser power drifts during esr'),
                      Parameter('ai_channel', 'ai4', ['ai0', 'ai1', 'ai2', 'ai3', 'ai4'],
                                'channel to use for analog input, to which the photodiode is connected')
                  ]),
    ]

    _INSTRUMENTS = {'microwave_generator_a': AgilentMicrowaveGenerator, 'NI6353': NI6353, 'PB': LISE607RTPulseBlaster
    }

    _SCRIPTS = {}

    def __init__(self, instruments, scripts=None, name=None, settings=None, log_function=None, data_path=None):
        Script.__init__(self, name, settings=settings, scripts=scripts, instruments=instruments,
                        log_function=log_function, data_path=data_path)

        self.daq_in = self.instruments['NI6353']['instance']

    def setup_microwave_gen(self):
        '''
        set relevant parameters on the MW generator.
        '''
        self.instruments['microwave_generator_a']['instance'].update({'amplitude': self.settings['power_out']})

        # There will be no counts recorded for this starting frequency. The MW generator will just wait for the pulse at this freq...
        freq_start = self.x_array[0] - (self.x_array[1] - self.x_array[0])
        print('ESR start frequency: {:0.2f} [Hz]'.format(float(self.x_array[0])))
        print('ESR stop frequency: {:0.2f} [Hz]'.format(float(self.x_array[-1])))

        self.instruments['microwave_generator_a']['instance'].update({'enable_modulation': False})
        self.instruments['microwave_generator_a']['instance'].update({'enable_IQ': False})

        self.instruments['microwave_generator_a']['instance'].update({'power_mode': 'CW'})
        self.instruments['microwave_generator_a']['instance'].update({'enable_output': True})

        # Avoid setting conflict issue!
        self.instruments['microwave_generator_a']['instance'].update({'freq_start': float(freq_start)})
        time.sleep(0.1)
        self.instruments['microwave_generator_a']['instance'].update({'freq_stop': float(self.x_array[-1])})
        time.sleep(0.1)
        self.instruments['microwave_generator_a']['instance'].update({'freq_stop': float(self.x_array[-1])})
        time.sleep(0.1)
        self.instruments['microwave_generator_a']['instance'].update({'freq_start': float(freq_start)})
        time.sleep(0.1)

        self.instruments['microwave_generator_a']['instance'].update({'freq_pts': self.settings['freq_points'] + 1})
        if self.settings['reverse_sweep']:
            self.instruments['microwave_generator_a']['instance'].update({'swp_direction': 'DOWN'})
        else:
            self.instruments['microwave_generator_a']['instance'].update({'swp_direction': 'UP'})

        self.instruments['microwave_generator_a']['instance'].update({'freq_mode': 'Sweep'})


        # print('just setup_microwave_gen')

    def setup_daq(self):
        '''
        Initialize the relevant DAQ and set the sample rate.
        '''

        self.daq_in = self.instruments['NI6353']['instance']
        # DAQ minimum buffer size is 2, so we break the integration time in half
        # sample_rate = float(1) / (self.settings['integration_time'] / self.settings['num_samps_per_pt'])
        sample_rate = float(1) / self.time_per_pt # self.time_per_pt is at least 10ms here

        self.daq_in.settings['digital_input']['ctr0']['sample_rate'] = sample_rate

    def setup_pb(self):  # ER 20181017
        '''
        Setup the channels on the PB card.
        '''

        self.instruments['PB']['instance'].update({'microwave_switch': {'status': True}})
        # if self.instruments['microwave_generator']['instance'].amplitude < -5.0:
        #     self.instruments['PB']['instance'].update({'microwave_switch': {'status': True}})
        # else:
        #     self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})

    def get_freq_array(self):
        '''
        Construct the frequency array based on the setting parameters.

        Returns:
            freq_values: array of the frequencies to be tested
            freq_range: the range of the frequency to be tested, i.e. maximum frequency - minumum frequency
        '''

        # construct the frequency array
        if self.settings['range_type'] == 'start_stop':
            if self.settings['freq_start'] > self.settings['freq_stop']:
                self.log('end freq. must be larger than start freq when range_type is start_stop. Abort script')
                self._abort = True

            if self.settings['freq_start'] < 9.0E3 or self.settings['freq_stop'] > 3.20E9:  # freq range of the RnS
                self.log('start or stop frequency out of bounds')
                self._abort = True

            freq_values = np.linspace(self.settings['freq_start'], self.settings['freq_stop'],
                                      self.settings['freq_points'])
            freq_range = max(freq_values) - min(freq_values)

        elif self.settings['range_type'] == 'center_range':
            if self.settings['freq_start'] < 2 * self.settings['freq_stop']:
                self.log(
                    'end freq. (range) must be smaller than 2x start freq (center) when range_type is center_range. Abort script')
                self._abort = True
            freq_values = np.linspace(self.settings['freq_start'] - self.settings['freq_stop'] / 2,
                                      self.settings['freq_start'] + self.settings['freq_stop'] / 2,
                                      self.settings['freq_points'])
            freq_range = max(freq_values) - min(freq_values)

            if self.settings['freq_stop'] > 1e9:
                self.log('freq_stop (range) is quite large --- did you mean to set \'range_type\' to \'start_stop\'? ')
        else:
            self.log('unknown range parameter. Abort script')
            self._abort = True

        return freq_values, freq_range

    def run_sweep(self):

        # initialize APD thread
        ctrtask = self.daq_in.setup_counter("ctr0",len(self.x_array) + 1)

        # Resets the sweep. With the next trigger event, the sweep starts with at the initial value.
        #self.instruments['microwave_generator_a']['instance'].reset_sweep()

        # Start counter and scanning sequence
        self.daq_in.run(ctrtask)
        xLineData, _ = self.daq_in.read(ctrtask)
        self.daq_in.stop(ctrtask)
        single_sweep_data = np.diff(xLineData)

        return single_sweep_data # return raw data

    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """

        # N9310A dwell time for freq sweep is 10ms - 10sec

        self.time_per_pt = self.settings['time_per_pt']
        if self.time_per_pt < 0.01:
            self.log('Minimum dwell time is 10ms.')
            self.time_per_pt = 0.01

        # turn on laser and apd_switch
        self.instruments['PB']['instance'].update({'laser': {'status': True}})
        self.instruments['PB']['instance'].update({'apd_switch': {'status': True}})

        # if tracking laser power drifts, check for PCI daq
        self.settings['track_laser_power']['on/off'] = False # Laser power tracking is not implenmented. Set it to be false for now.

        if self.settings['track_laser_power']['on/off'] and not self.settings['daq_type'] == 'PCI':
            print("tracking laser power drifts only enabled for PCI daq")
            self._abort = True

        self.lines = []
        take_ref = self.settings['norm_to_ref']

        # get the frequencices of the sweep
        freq_values, freq_range = self.get_freq_array()
        self.x_array = freq_values

        # initialize some of the fields in self.data
        self.data = {'frequency': freq_values, 'data': [], 'data_err': None, 'fit_params': None, 'avrg_counts': None}
        self.data['current_avg_cnts'] = 0.
        self.data['total_avg_cnts'] = 0.
        # initialize data arrays
        esr_data = np.zeros((self.settings['esr_avg_max'], len(freq_values)))  # for the raw esr data
        # laser_data = np.zeros((self.settings['esr_avg'], len(freq_values)))  # for the raw photodiode data
        avrg_counts = np.zeros(self.settings['esr_avg_max'])  # average counts for EACH average to normalize the plot if take_ref is true

        # setup the microwave generator
        self.setup_microwave_gen()

        # setup the daq
        self.setup_daq()

        # setup the pulseblaster card (i.e., turn mw_switch on)
        self.setup_pb()

        # run sweeps
        for scan_num in range(0, self.settings['esr_avg_max']):

            esr_data_pos = 0
            if self._abort:
                break

            # get the data for a single sweep. These are raw data.
            # single_sweep_data, single_sweep_laser_data = self.run_sweep(freq_values)
            single_sweep_data = self.run_sweep()

            # save the single sweep data and normalize to kcounts/sec
            esr_data[scan_num, esr_data_pos:(esr_data_pos + len(single_sweep_data))] = single_sweep_data * (
                        .001 / self.settings['time_per_pt'])


            # laser_data[scan_num, esr_data_pos:(esr_data_pos + len(single_sweep_laser_data))] = single_sweep_laser_data
            esr_data_pos += len(single_sweep_data)

            # average counts of the expt
            avrg_counts[scan_num] = np.mean(esr_data[scan_num])
            # self.current_avg_cnts = np.mean(avrg_counts[scan_num])
            self.current_avg_cnts = avrg_counts[scan_num]
            # print('scan_num', scan_num)
            # print('self.current_avg_cnts is', self.current_avg_cnts)
            self.data.update({'current_avg_cnts': self.current_avg_cnts}) # save as data

            if take_ref is True:
                esr_data[scan_num] /= avrg_counts[scan_num]

            # # normalize the data if track_laser
            # if self.settings['track_laser_power']['on/off']:
            #     laser_norm_data = np.divide(esr_data, laser_data)
            #     # average of the normalized data for the number of averages completed so far, to plot and fit to if laser power tracking is on
            #     data_laser_norm = (np.mean(laser_norm_data[0:(scan_num + 1)], axis=0))  # *tmp_laser

            # current non-normalized averaged data to plot and fit to if laser power tracking is off
            esr_avg = np.mean(esr_data[0:(scan_num + 1)], axis=0)
            self.total_avg_cnts = np.mean(avrg_counts[0:(scan_num + 1)], axis=0)
            self.data.update({'total_avg_cnts': self.total_avg_cnts})

            # if not self.settings['track_laser_power']['on/off']:
            #     # fit to the data
            #     fit_params = fit_esr(freq_values, esr_avg, min_counts=self.settings['fit_constants']['minimum_counts'],
            #                          contrast_factor=self.settings['fit_constants']['contrast_factor'])
            # elif self.settings['track_laser_power']['on/off']:
            #     # fit to the data
            #     fit_params = fit_esr(freq_values, data_laser_norm,
            #                          min_counts=self.settings['fit_constants']['minimum_counts'],
            #                          contrast_factor=self.settings['fit_constants']['contrast_factor'])
            #
            #     # save the data
            #     self.data.update({'laser_data': laser_data})
            #     self.data.update({'data_laser_norm': data_laser_norm})

            # do fitting
            fit_params = fit_esr(freq_values, esr_avg, min_counts=self.settings['fit_constants']['minimum_counts'],
                                 contrast_factor=self.settings['fit_constants']['contrast_factor'],
                                 num_of_peaks=self.settings['fit_constants']['num_of_peaks'])
            # if fit_params is not None:
            #     print('fit para number', len(fit_params))
            #     print('fit para ', fit_params)

            self.data.update({'data': esr_avg, 'fit_params': fit_params, 'avrg_counts': avrg_counts})

            if scan_num > 0:
                esr_avg_err = np.std(esr_data[0:(scan_num + 1)], axis=0)/np.sqrt(scan_num)
                self.data.update({'data_err': esr_avg_err})
            #     self.data.update({'data': avrg_counts})  # ER 20181022

            if self.settings['save_full_esr']:
                self.data.update({'esr_data': esr_data[0:(scan_num + 1), :]})

            # Break out of the loop if # of averages is enough and a good fit has been found.
            if scan_num >= self.settings['esr_avg_min'] and fit_params is not None:
                break

            progress = self._calc_progress(scan_num)
            self.updateProgress.emit(progress)

        # Turn off the microwave switch channel
        self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})

        # Turn off laser and apd_switch
        self.instruments['PB']['instance'].update({'laser': {'status': False}})
        self.instruments['PB']['instance'].update({'apd_switch': {'status': False}})

        # This is used to flag when experiment is done and it allows proper saving to datasets
        self.data['exp_finished'] = 1

        # Bring the MW generator to the normal state
        self.instruments['microwave_generator_a']['instance'].update({'power_mode': 'CW'})
        self.instruments['microwave_generator_a']['instance'].update({'freq_mode': 'CW'})
        # time.sleep(0.1)
        if self.settings['turn_off_after']:
            self.instruments['microwave_generator_a']['instance'].update({'enable_output': False})

    def _calc_progress(self, scan_num):
        # COMMENT_ME
        progress = float(scan_num) / self.settings['esr_avg_min'] * 100.
        self.progress = progress
        return int(progress)

    def _plot(self, axes_list, data=None):
        """
        plotting function for esr
        Args:
            axes_list: list of axes objects on which to plot plots the esr on the first axes object
            data: data (dictionary that contains keys frequency, data and fit_params) if not provided use self.data
        Returns:

        """

        if data is None:
            data = self.data

        if 'exp_finished' in data.keys() is not None and data['exp_finished'] == 1: # the final plot after the exp is done
            mean_counts = data['total_avg_cnts']  # averaged counts of all the scans
        else:
            mean_counts = data['current_avg_cnts']  # averaged counts of the current scan

        if not self.settings['track_laser_power']['on/off']:

            plot_esr(axes_list[0], data['frequency'], data['data'], data['fit_params'], avg_counts=mean_counts,
                     mw_power=self.settings['power_out'], D=self.settings['fit_constants']['zfs'],
                     gama=self.settings['fit_constants']['gama'], err=data['data_err'])
        elif self.settings['track_laser_power']['on/off'] and self.settings['daq_type'] == 'PCI':

            plot_esr(axes_list[0], data['frequency'], data['data_laser_norm'], data['fit_params'],
                     avg_counts=mean_counts, mw_power=self.settings['power_out'],
                     D=self.settings['fit_constants']['zfs'],
                     gama=self.settings['fit_constants']['gama'], err=data['data_err'])

    def get_axes_layout(self, figure_list):
        """
        returns the axes objects the script needs to plot its data
        the default creates a single axes object on each figure
        This can/should be overwritten in a child script if more axes objects are needed
        Args:
            figure_list: a list of figure objects
        Returns:
            axes_list: a list of axes objects

        """
        new_figure_list = [figure_list[1]]
        # if 'exp_finished' in self.data.keys() is not None and self.data['exp_finished'] == 1:
        #     pass
        # else:
        #     return super(ESR_RnS, self).get_axes_layout(new_figure_list)

        return super(ESR_N9310A, self).get_axes_layout(new_figure_list)

if __name__ == '__main__':
    script = {}
    instr = {}
    script, failed, instr = Script.load_and_append({'ESR_N9310A': 'ESR_N9310A'}, script, instr)

    print(script)
    print(failed)
    print(instr)