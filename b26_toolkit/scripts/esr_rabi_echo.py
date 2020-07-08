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

# from b26_toolkit.scripts import ESR, Rabi
from b26_toolkit.scripts.esr_RnS import ESR_FastSwp_RnS_FitGuaranteed
from b26_toolkit.scripts.pulse_sequences.rabi import Rabi_RnS
from b26_toolkit.scripts.pulse_sequences.hahn_echo import PDD_RnS
from pylabcontrol.core import Script, Parameter
from collections import deque
import numpy as np
# from b26_toolkit.instruments import R8SMicrowaveGenerator, NI6353, LISE607RTPulseBlaster
# import time
# from b26_toolkit.scripts.daq_read_counter import Daq_Read_Counter





class ESR_Rabi_RnS(Script):
    """
    Does both an ESR experiment and a Rabi experiment on an NV, using the reference frequency from the esr data.
    Experiments are done with Rhode&Schwartz microwave generator SMB100A and no IQ modulation is available on it.

    --> Last edited by Ziwei Qiu 3/29/2019

    """

    _DEFAULT_SETTINGS = [

        Parameter('esr_settings', [
            Parameter('esr_mw_pwr', -10, float, 'microwave power for ESR scan'),
            Parameter('esr_avg_min', 12, int, 'minimum number of esr averages'),
            Parameter('esr_avg_max', 50, int, 'maximum number of esr averages'),
            Parameter('esr_cntr_freq', 2.82e9, float, 'center frequency for ESR scan'),
            Parameter('esr_freq_range', 8.5e7, float, 'frequency range for ESR scan (suggest 6e7 - 9e7)'),
            Parameter('esr_num_of_pts', 65, int, 'number of frequency points for ESR scan'),
            Parameter('esr_time_per_pt', 0.02, [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.1],
                      'integration time for each point in the fast ESR scan (suggest 0.02-0.04)'),
            Parameter('esrfit_num_of_peaks', 1, [-1, 1, 2],
                      'specify number of peaks for fitting. if not specifying the number of peaks, choose -1'),
            Parameter('esrfit_minimum_counts', 0.9, float,
                      'minumum counts for an ESR to not be considered noise (suggest 0.8 - 1.01 if esr is normalized)'),
            Parameter('esrfit_contrast_factor', 3.0, float,
                      'minimum contrast for an ESR to not be considered noise (suggest 3.0-4.0)')
        ]),
        Parameter('rabi_settings', [
            Parameter('tau_min', 15, float, 'minimum time for rabi oscillations (in ns)'),
            Parameter('tau_max', 1000, float, 'total time of rabi oscillations (in ns)'),
            Parameter('tau_step', 20., [2.5, 5., 10., 20., 50., 100., 200., 500., 1000., 10000., 100000., 500000.],
                      'time step increment of rabi pulse duration (in ns)'),
            Parameter('num_averages', 500000, int, 'number of averages (>100000)'),
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 100, float, 'apd readout duration (in ns)'),
            Parameter('nv_reset_time', 2000, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 600, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 530, int,
                      '[ns] delay between laser on and APD readout (given by spontaneous decay rate)')
        ])
    ]

    _INSTRUMENTS = { }

    _SCRIPTS = {'esr': ESR_FastSwp_RnS_FitGuaranteed, 'rabi': Rabi_RnS}

    def __init__(self, instruments = None, scripts = None, name = None, settings = None, log_function = None, data_path = None):
        """
        Example of a script that makes use of an instrument
        Args:
            instruments: instruments the script will make use of
            name (optional): name of script, if empty same as class name
            settings (optional): settings for this script, if empty same as default settings
        """

        # call init of superclass
        Script.__init__(self, name, settings=settings, instruments=instruments, scripts=scripts,
                        log_function=log_function, data_path=data_path)

    def _function(self):

        # need to create the data, even though it is empty, in order for the plot to refresh!!
        self.data = {'esr_fo': deque()}

        # Set the right parameters for the ESR scan
        self.scripts['esr'].settings['power_out'] = self.settings['esr_settings']['esr_mw_pwr']
        self.scripts['esr'].settings['esr_avg_min'] = self.settings['esr_settings']['esr_avg_min']
        self.scripts['esr'].settings['esr_avg_max'] = self.settings['esr_settings']['esr_avg_max']
        self.scripts['esr'].settings['freq_start'] = self.settings['esr_settings']['esr_cntr_freq']
        self.scripts['esr'].settings['freq_stop'] = self.settings['esr_settings']['esr_freq_range']
        self.scripts['esr'].settings['range_type'] = 'center_range'
        self.scripts['esr'].settings['freq_points'] = self.settings['esr_settings']['esr_num_of_pts']
        self.scripts['esr'].settings['time_per_pt'] = self.settings['esr_settings']['esr_time_per_pt']
        self.scripts['esr'].settings['fit_constants']['num_of_peaks'] = self.settings['esr_settings'][
            'esrfit_num_of_peaks']
        self.scripts['esr'].settings['fit_constants']['minimum_counts'] = self.settings['esr_settings'][
            'esrfit_minimum_counts']
        self.scripts['esr'].settings['fit_constants']['contrast_factor'] = self.settings['esr_settings'][
            'esrfit_contrast_factor']

        # Set the right parameters for Rabi
        self.scripts['rabi'].settings['tau_times']['min_time'] = self.settings['rabi_settings']['tau_min']
        self.scripts['rabi'].settings['tau_times']['max_time'] = self.settings['rabi_settings']['tau_max']
        self.scripts['rabi'].settings['tau_times']['time_step'] = self.settings['rabi_settings']['tau_step']
        self.scripts['rabi'].settings['num_averages'] = self.settings['rabi_settings']['num_averages']
        self.scripts['rabi'].settings['read_out']['meas_time'] = self.settings['read_out']['meas_time']
        self.scripts['rabi'].settings['read_out']['nv_reset_time'] = self.settings['read_out']['nv_reset_time']
        self.scripts['rabi'].settings['read_out']['laser_off_time'] = self.settings['read_out']['laser_off_time']
        self.scripts['rabi'].settings['read_out']['delay_mw_readout'] = self.settings['read_out']['delay_mw_readout']
        self.scripts['rabi'].settings['read_out']['delay_readout'] = self.settings['read_out']['delay_readout']


        # Start ESR
        print('==> Start measuring ESR...')
        # print('self._plot_refresh', self._plot_refresh)
        self.scripts['esr'].run()
        esr_fit_data = self.scripts['esr'].data['fit_params']
        print('     len(esr_fit_data) =  ', esr_fit_data)

        if self.scripts['esr'].data['fit_params'] is not None:
            if len(self.scripts['esr'].data['fit_params']) == 4:
                self.rabi_frequency = self.scripts['esr'].data['fit_params'][2]
            elif len(self.scripts['esr'].data['fit_params']) == 6:
                self.rabi_frequency = self.scripts['esr'].data['fit_params'][4]
            else:
                self.log('Could not get fit parameters from esr script. Experiment aborted.')
                print('Could not get fit parameters from esr script. Experiment aborted.')
                self._abort = True
                # raise RuntimeError('Could not get fit parameters from esr script')

            # if self.rabi_frequency < self.scripts['esr'].settings['freq_start']:
            if self.rabi_frequency < (self.scripts['esr'].settings['freq_start'] - 0.5 * self.scripts['esr'].settings['freq_stop']):
                self.log('Resonance frequency found ({:0.2e}) Hz was below esr sweep range, aborting rabi attempt'.format(self.rabi_frequency))
            # elif self.rabi_frequency > self.scripts['esr'].settings['freq_stop']:
            elif self.rabi_frequency > (self.scripts['esr'].settings['freq_start'] + 0.5 * self.scripts['esr'].settings['freq_stop']):
                self.log('Resonance frequency found ({:0.2e}) Hz was above esr sweep range, aborting rabi attempt'.format(self.rabi_frequency))
            else:
                self.log(
                    'Starting rabi with frequency {:.4f} GHz and power {:.2f} dBm'.format(self.rabi_frequency / 1E9,
                                                                                          self.settings['esr_settings'][
                                                                                              'esr_mw_pwr']))
                print('Starting rabi with frequency {:.4f} GHz and power {:.2f} dBm'.format(self.rabi_frequency / 1E9,
                                                                                          self.settings['esr_settings'][
                                                                                              'esr_mw_pwr']))
                self.scripts['rabi'].settings['mw_pulses']['mw_power'] = self.settings['esr_settings']['esr_mw_pwr']
                self.scripts['rabi'].settings['mw_pulses']['mw_frequency'] = float(self.rabi_frequency)

                self.scripts['rabi'].run()
        else:
            self.log('No resonance frequency found skipping rabi attempt')

    def _plot(self, axes_list, data=None):

        # print('(plot)')
        if self._current_subscript_stage['current_subscript'] is self.scripts['esr'] and self.scripts['esr'].is_running:
            self.scripts['esr']._plot([axes_list[1]])
        elif self._current_subscript_stage['current_subscript'] is self.scripts['rabi'] and self.scripts['rabi'].is_running:
            self.scripts['rabi']._plot(axes_list)

    def _update_plot(self, axes_list):
        # print('(update plot)')

        if self._current_subscript_stage['current_subscript'] is self.scripts['esr'] and self.scripts['esr'].is_running:
            self.scripts['esr']._update_plot([axes_list[1]])
        elif self._current_subscript_stage['current_subscript'] is self.scripts['rabi'] and self.scripts['rabi'].is_running:
            self.scripts['rabi']._update_plot(axes_list)


class ESR_Rabi_Echo_RnS(Script):
    """
        Does an ESR experiment, a Rabi experiment and a Hahn-echo experiment on an NV.
        Experiments are done with Rhode&Schwartz microwave generator SMB100A and no IQ modulation is available on it.

        --> Last edited by Ziwei Qiu 3/30/2019
    """

    _DEFAULT_SETTINGS = [

        Parameter('to-do', 'esr_rabi_echo', ['esr', 'esr_rabi','esr_rabi_echo'], 'Choose what experiments to do'),
        Parameter('esr_settings', [
            Parameter('esr_mw_pwr', -12, float, 'microwave power for ESR scan'),
            Parameter('esr_avg_min', 25, int, 'minimum number of esr averages'),
            Parameter('esr_avg_max', 50, int, 'maximum number of esr averages'),
            Parameter('esr_cntr_freq', 2.92e9, float, 'center frequency for ESR scan'),
            Parameter('esr_freq_range', 8e7, float, 'frequency range for ESR scan (suggest 6e7 - 9e7)'),
            Parameter('esr_num_of_pts', 121, int, 'number of frequency points for ESR scan'),
            Parameter('esr_time_per_pt', 0.02, [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.1],
                      'integration time for each point in the fast ESR scan (suggest 0.02-0.04)'),
            Parameter('esrfit_num_of_peaks', 1, [-1, 1, 2],
                      'specify number of peaks for fitting. if not specifying the number of peaks, choose -1'),
            Parameter('esrfit_minimum_counts', 0.9, float,
                      'minumum counts for an ESR to not be considered noise (suggest 0.8 - 1.01 if esr is normalized)'),
            Parameter('esrfit_contrast_factor', 3.0, float,
                      'minimum contrast for an ESR to not be considered noise (suggest 3.0-4.0)')
        ]),
        Parameter('rabi_settings', [
            Parameter('rabi_mw_pwr', -8, float, 'microwave power for ESR scan'),
            Parameter('tau_min', 15, float, 'minimum time for rabi oscillations (in ns)'),
            Parameter('tau_max', 1000, float, 'total time of rabi oscillations (in ns)'),
            Parameter('tau_step', 20., [2.5, 5., 10., 20., 50., 100., 200., 500., 1000., 10000., 100000., 500000.],
                      'time step increment of rabi pulse duration (in ns)'),
            Parameter('num_averages', 500000, int, 'number of averages (>100000), 1 block = 100000'),
        ]),
        Parameter('echo_settings', [
            Parameter('tau_min', 200, float, 'minimum time for rabi oscillations (in ns)'),
            Parameter('tau_max', 10000, float, 'total time of rabi oscillations (in ns)'),
            Parameter('tau_step', 100., [2.5, 5., 10., 20., 50., 100., 200., 500., 1000., 10000., 100000., 500000.],
                      'time step increment of rabi pulse duration (in ns)'),
            Parameter('num_averages', 1000000, int, 'number of averages (>100000), 1 block = 100000'),
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 100, float, 'apd readout duration (in ns)'),
            Parameter('nv_reset_time', 2000, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 600, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 530, int,
                      '[ns] delay between laser on and APD readout (given by spontaneous decay rate)')
        ]),
        Parameter('Tracking', [
            Parameter('on/off', True, bool, 'used to turn on tracking for both rabi and echo'),
            Parameter('threshold', 0.7, float, 'threshold for tracking'),
            Parameter('init_fluor', 80, float, 'initial fluor of NV to compare to, in kcps')
        ])
    ]

    _INSTRUMENTS = {}

    _SCRIPTS = {'esr': ESR_FastSwp_RnS_FitGuaranteed, 'rabi': Rabi_RnS, 'echo': PDD_RnS}

    def __init__(self, instruments = None, scripts = None, name = None, settings = None, log_function = None, data_path = None):
        """
        Example of a script that makes use of an instrument
        Args:
            instruments: instruments the script will make use of
            name (optional): name of script, if empty same as class name
            settings (optional): settings for this script, if empty same as default settings
        """

        # call init of superclass
        Script.__init__(self, name, settings=settings, instruments=instruments, scripts=scripts,
                        log_function=log_function, data_path=data_path)

    def _function(self):

        # need to create the data, even though it is empty, in order for the plot to refresh!!
        self.data = {'mw_pulses': deque()}

        # Set the right parameters for the ESR scan
        self.scripts['esr'].settings['power_out'] = self.settings['esr_settings']['esr_mw_pwr']
        self.scripts['esr'].settings['esr_avg_min'] = self.settings['esr_settings']['esr_avg_min']
        self.scripts['esr'].settings['esr_avg_max'] = self.settings['esr_settings']['esr_avg_max']
        self.scripts['esr'].settings['freq_start'] = self.settings['esr_settings']['esr_cntr_freq']
        self.scripts['esr'].settings['freq_stop'] = self.settings['esr_settings']['esr_freq_range']
        self.scripts['esr'].settings['range_type'] = 'center_range'
        self.scripts['esr'].settings['freq_points'] = self.settings['esr_settings']['esr_num_of_pts']
        self.scripts['esr'].settings['time_per_pt'] = self.settings['esr_settings']['esr_time_per_pt']
        self.scripts['esr'].settings['fit_constants']['num_of_peaks'] = self.settings['esr_settings'][
            'esrfit_num_of_peaks']
        self.scripts['esr'].settings['fit_constants']['minimum_counts'] = self.settings['esr_settings'][
            'esrfit_minimum_counts']
        self.scripts['esr'].settings['fit_constants']['contrast_factor'] = self.settings['esr_settings'][
            'esrfit_contrast_factor']

        # Set the right parameters for Rabi
        self.scripts['rabi'].settings['mw_pulses']['mw_power'] = self.settings['rabi_settings']['rabi_mw_pwr']
        self.scripts['rabi'].settings['tau_times']['min_time'] = self.settings['rabi_settings']['tau_min']
        self.scripts['rabi'].settings['tau_times']['max_time'] = self.settings['rabi_settings']['tau_max']
        self.scripts['rabi'].settings['tau_times']['time_step'] = self.settings['rabi_settings']['tau_step']
        self.scripts['rabi'].settings['num_averages'] = self.settings['rabi_settings']['num_averages']
        self.scripts['rabi'].settings['read_out']['meas_time'] = self.settings['read_out']['meas_time']
        self.scripts['rabi'].settings['read_out']['nv_reset_time'] = self.settings['read_out']['nv_reset_time']
        self.scripts['rabi'].settings['read_out']['laser_off_time'] = self.settings['read_out']['laser_off_time']
        self.scripts['rabi'].settings['read_out']['delay_mw_readout'] = self.settings['read_out']['delay_mw_readout']
        self.scripts['rabi'].settings['read_out']['delay_readout'] = self.settings['read_out']['delay_readout']
        self.scripts['rabi'].settings['Tracking']['on/off'] = self.settings['Tracking']['on/off']
        self.scripts['rabi'].settings['Tracking']['threshold'] = self.settings['Tracking']['threshold']
        self.scripts['rabi'].settings['Tracking']['init_fluor'] = self.settings['Tracking']['init_fluor']

        # Set the right parameters for Echo
        self.scripts['echo'].settings['mw_pulses']['mw_power'] = self.settings['rabi_settings']['rabi_mw_pwr']
        self.scripts['echo'].settings['tau_times']['min_time'] = self.settings['echo_settings']['tau_min']
        self.scripts['echo'].settings['tau_times']['max_time'] = self.settings['echo_settings']['tau_max']
        self.scripts['echo'].settings['tau_times']['time_step'] = self.settings['echo_settings']['tau_step']
        self.scripts['echo'].settings['num_averages'] = self.settings['echo_settings']['num_averages']
        self.scripts['echo'].settings['read_out']['meas_time'] = self.settings['read_out']['meas_time']
        self.scripts['echo'].settings['read_out']['nv_reset_time'] = self.settings['read_out']['nv_reset_time']
        self.scripts['echo'].settings['read_out']['laser_off_time'] = self.settings['read_out']['laser_off_time']
        self.scripts['echo'].settings['read_out']['delay_mw_readout'] = self.settings['read_out']['delay_mw_readout']
        self.scripts['echo'].settings['read_out']['delay_readout'] = self.settings['read_out']['delay_readout']
        self.scripts['echo'].settings['Tracking']['on/off'] = self.settings['Tracking']['on/off']
        self.scripts['echo'].settings['Tracking']['threshold'] = self.settings['Tracking']['threshold']
        self.scripts['echo'].settings['Tracking']['init_fluor'] = self.settings['Tracking']['init_fluor']
        self.scripts['echo'].settings['decoupling_seq']['type'] = 'spin_echo'
        self.scripts['echo'].settings['decoupling_seq']['num_of_pulse_blocks'] = 1


        # Start ESR
        print('==> Start measuring ESR...')
        # print('self._plot_refresh', self._plot_refresh)
        self.scripts['esr'].run()
        esr_fit_data = self.scripts['esr'].data['fit_params']
        print('     len(esr_fit_data) =  ', esr_fit_data)

        # Start Rabi
        if self.settings['to-do'] ==  'esr_rabi' or  'esr_rabi_echo':
            print('==> Start measuring Rabi...')
            if self.scripts['esr'].data['fit_params'] is not None:
                if len(self.scripts['esr'].data['fit_params']) == 4:
                    self.rabi_frequency = self.scripts['esr'].data['fit_params'][2]
                elif len(self.scripts['esr'].data['fit_params']) == 6:
                    self.rabi_frequency = self.scripts['esr'].data['fit_params'][4]
                else:
                    self.log('Could not get fit parameters from esr script. Experiment aborted.')
                    print('Could not get fit parameters from esr script. Experiment aborted.')
                    self._abort = True
                    # raise RuntimeError('Could not get fit parameters from esr script')

                # if self.rabi_frequency < self.scripts['esr'].settings['freq_start']:
                if self.rabi_frequency < (self.scripts['esr'].settings['freq_start'] - 0.5 * self.scripts['esr'].settings['freq_stop']):
                    self.log(
                        'Resonance frequency found ({:0.2e}) Hz was below esr sweep range, aborting rabi attempt (if any)'.format(
                            self.rabi_frequency))
                    print(
                        'Resonance frequency found ({:0.2e}) Hz was below esr sweep range, aborting rabi attempt (if any)'.format(
                            self.rabi_frequency))
                    self._abort = True
                # elif self.rabi_frequency > self.scripts['esr'].settings['freq_stop']:
                elif self.rabi_frequency > (self.scripts['esr'].settings['freq_start'] + 0.5 * self.scripts['esr'].settings['freq_stop']):
                    self.log(
                        'Resonance frequency found ({:0.2e}) Hz was above esr sweep range, aborting rabi attempt (if any)'.format(
                            self.rabi_frequency))
                    print(
                        'Resonance frequency found ({:0.2e}) Hz was above esr sweep range, aborting rabi attempt (if any)'.format(
                            self.rabi_frequency))
                    self._abort = True
                else:
                    self.log(
                        'Starting rabi with frequency {:.4f} GHz and power {:.2f} dBm'.format(self.rabi_frequency / 1E9,
                                                                                              self.settings['rabi_settings'][
                                                                                                  'rabi_mw_pwr']))
                    print('Starting rabi with frequency {:.4f} GHz and power {:.2f} dBm'.format(self.rabi_frequency / 1E9,
                                                                                              self.settings['rabi_settings'][
                                                                                                  'rabi_mw_pwr']))
                    self.scripts['rabi'].settings['mw_pulses']['mw_frequency'] = float(self.rabi_frequency)
                    self.data['mw_pulses'].append(float(self.settings['rabi_settings']['rabi_mw_pwr']))
                    self.data['mw_pulses'].append(float(self.rabi_frequency))
                    self.flag_rabi_plot = True
                    self.scripts['rabi'].run()
            else:
                self.log('No resonance frequency found skipping rabi attempt (if any)')
                print('No resonance frequency found skipping rabi attempt (if any)')
                self._abort = True

        # Start Echo
        if self.settings['to-do'] == 'esr_rabi_echo':

            if 'fits' in self.scripts['rabi'].data.keys() and self.scripts['rabi'].data['fits'] is not None:
                self.pi_time = self.scripts['rabi'].data['pi_time']
                self.pi_half_time = self.scripts['rabi'].data['pi_half_time']
                self.three_pi_half_time = self.scripts['rabi'].data['three_pi_half_time']

                if self.three_pi_half_time >  self.pi_time > self.pi_half_time > 20:

                    self.scripts['echo'].settings['mw_pulses']['mw_frequency'] = float(self.rabi_frequency)
                    self.scripts['echo'].settings['mw_pulses']['pi_pulse_time'] = float(self.pi_time)
                    self.scripts['echo'].settings['mw_pulses']['pi_half_pulse_time'] = float(self.pi_half_time)
                    self.scripts['echo'].settings['mw_pulses']['3pi_half_pulse_time'] = float(self.three_pi_half_time)
                    self.data['mw_pulses'].append(float(self.pi_time))
                    self.data['mw_pulses'].append(float(self.pi_half_time))
                    self.data['mw_pulses'].append(float(self.three_pi_half_time))
                    self.data['rabi_T2_star'] = float(self.scripts['rabi'].data['T2_star'])
                    self.log('Starting Echo experiment.')
                    print('Starting Echo experiment.')
                    self.flag_echo_plot = True
                    self.scripts['echo'].run()
                    if 'fits' in self.scripts['echo'].data.keys() and self.scripts['echo'].data['fits'] is not None:
                        self.data['echo_T2'] = float(self.scripts['echo'].data['fits'][1])

                else:
                    self.log('Rabi fits are not real, aborting echo attempt (if any)')
                    print('Rabi fits are not real, aborting echo attempt (if any)')
                    self._abort = True
            else:
                self.log('No fitting from the rabi, skipping echo attempt (if any)')
                print('No fitting from the rabi, skipping echo attempt (if any)')
                self._abort = True

        # for proper saving, deque objects need to be converted to np.arrays or list
        if 'mw_pulses' in self.data.keys() is not None:
            # self.data['counts'] = list(self.data['counts'])
            self.data['mw_pulses'] = np.asarray(self.data['mw_pulses'])

    def _plot(self, axes_list, data=None):

        if self._current_subscript_stage['current_subscript'] is self.scripts['esr'] and self.scripts['esr'].is_running:
            self.scripts['esr']._plot([axes_list[1]])
        elif self._current_subscript_stage['current_subscript'] is self.scripts['rabi'] and self.scripts['rabi'].is_running:
            self.scripts['rabi']._plot(axes_list)
        elif self._current_subscript_stage['current_subscript'] is self.scripts['echo'] and self.scripts['echo'].is_running:
            self.scripts['echo']._plot(axes_list)

    def _update_plot(self, axes_list):

        if self._current_subscript_stage['current_subscript'] is self.scripts['esr'] and self.scripts['esr'].is_running:
            self.scripts['esr']._update_plot([axes_list[1]])
        elif self._current_subscript_stage['current_subscript'] is self.scripts['rabi'] and self.scripts['rabi'].is_running:
            # self.scripts['rabi']._update_plot(axes_list)
            if self.flag_rabi_plot:

                self.scripts['rabi']._plot([axes_list[0],axes_list[1]])
                self.flag_rabi_plot = False
            else:

                self.scripts['rabi']._update_plot([axes_list[0], axes_list[1]])
        elif self._current_subscript_stage['current_subscript'] is self.scripts['echo'] and self.scripts['echo'].is_running:
            if self.flag_echo_plot:
                self.scripts['echo']._plot([axes_list[0],axes_list[1]])
                self.flag_echo_plot = False
            else:
                self.scripts['echo']._update_plot([axes_list[0],axes_list[1]])

    def get_axes_layout(self, figure_list):
        """
        returns the axes objects the script needs to plot its data
        this overwrites the default get_axis_layout in PyLabControl.src.core.scripts
        Args:
            figure_list: a list of figure objects
        Returns:
            axes_list: a list of axes objects

        """
        axes_list = []
        if self._plot_refresh is True:
            for fig in figure_list:
                fig.clf()
            axes_list.append(figure_list[0].add_subplot(111))  # axes_list[0]
            axes_list.append(figure_list[1].add_subplot(111))  # axes_list[1]
        else:
            axes_list.append(figure_list[0].axes[0])
            axes_list.append(figure_list[1].axes[0])
        return axes_list



if __name__ == '__main__':
    script = {}
    instr = {}
    script, failed, instr = Script.load_and_append({'ESR_Rabi_RnS': 'ESR_Rabi_RnS'}, script, instr)

    print(script)
    print(failed)
    print(instr)


