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
import numpy as np
from b26_toolkit.scripts.pulse_sequences.pulsed_experiment_base_script import PulsedExperimentBaseScript
# from b26_toolkit.instruments import NI6259, NI9402, B26PulseBlaster, MicrowaveGenerator, Pulse
from b26_toolkit.instruments import NI6353, LISE607RTPulseBlaster, AgilentMicrowaveGenerator, R8SMicrowaveGenerator, Pulse, AgilentMicrowaveGeneratorII
from pylabcontrol.core import Parameter
from b26_toolkit.data_processing.fit_functions import cose_with_decay, fit_exp_decay
from b26_toolkit.plotting.plots_1d import plot_1d_simple_timetrace_ns, plot_pulses

MIN_DURATION = 15 # minimum allowed duration in [ns] for pulse blaster pulses

# class GrDelayMeas(PulsedExperimentBaseScript):
#     """
#         This script measures the green laser delay during AOM turn ON and turn OFF
#         ==> Last edited by ZQ 2/26/2019 7:05pm
#     """
#     _DEFAULT_SETTINGS = [
#         Parameter('tau_times', [
#             Parameter('min_time', 20, float, 'minimum green delay (in ns)'),
#             Parameter('max_time', 800, float, 'max green delay (in ns)'),
#             Parameter('time_step', 10, [10,20,30,40,50,80,100],'time step increment of green delay (in ns)')
#         ]),
#         Parameter('read_out', [
#             Parameter('meas_time', 30, int, '[ns] APD window to count photons'),
#             Parameter('green_time', 500, float, '[ns] duration of green pulse '),
#             Parameter('laser_off_time', 1000, float, '[ns] dark time before AOM turn on')
#         ]),
#         Parameter('num_averages', 100000, int, 'number of averages'),
#         # Parameter('APD_type', 'SPCM', ['SPCM', 'MPPC'], 'SPCM = single-photon counting module (digital input on DAQ); MPPC = multi-pixel photon counter (analog input on DAQ)'),
#         # Parameter('MPPC_ai_channel', 'ai4', ['ai0', 'ai1', 'ai2', 'ai3', 'ai4', 'ai5', 'ai6', 'ai7', 'ai8', 'ai9'],
#         #           'NI DAQ analog input channel for MPPC voltage'),
#         # Parameter('skip_invalid_sequences', True, bool, 'Skips any sequences with <MIN_DURATIONns commands')
#     ]
#
#     _INSTRUMENTS = {'daq': NI6353, 'PB': LISE607RTPulseBlaster}
#
#     def _function(self):
#         super(GrDelayMeas, self)._function(self.data)
#         # counts = self.data['counts']
#         # tau = self.data['tau']
#
#         ## Turn off green light (the pulse blaster will pulse it on when needed)
#         self.instruments['PB']['instance'].update({'laser': {'status': False}})
#
#         self.data[
#             'exp_finished'] = 1  # this is to flag when the experiment is done or stopped and this allows proper saving even fit fails
#
#         # turn off laser, apd_switch and MW (just in case)
#         self.instruments['PB']['instance'].update({'laser': {'status': False}})
#         self.instruments['PB']['instance'].update({'apd_switch': {'status': False}})
#
#
#     def _create_pulse_sequences(self):
#         '''
#
#         Returns: pulse_sequences, num_averages, tau_list, meas_time
#             pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
#             scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
#             sequence must have the same number of daq read pulses
#             num_averages: the number of times to repeat each pulse sequence
#             tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
#             meas_time: the width (in ns) of the daq measurement
#
#         '''
#         pulse_sequences = []
#
#         # JG 16-08-25 changed (MIN_DURATION min spacing is taken care of later):
#         # tau_list = range(0, int(self.settings['tau_times']['max_time']), self.settings['tau_times']['time_step'])
#         # tau_list = range(int(self.settings['tau_times']['min_time']), int(self.settings['tau_times']['max_time']),self.settings['tau_times']['time_step'])
#
#         max_range = int(np.floor((self.settings['tau_times']['max_time'] - self.settings['tau_times']['min_time']) /
#                                  self.settings['tau_times']['time_step']))
#         tau_list = np.array(
#             [self.settings['tau_times']['min_time'] + i * self.settings['tau_times']['time_step'] for i in
#              range(max_range)])
#
#         # ignore the sequence if the mw-pulse is shorter than MIN_DURATION (0 is ok because there is no mw pulse!)
#         tau_list = [x for x in tau_list if x == 0 or x >= MIN_DURATION]
#         print('tau_list', tau_list)
#
#         meas_time = self.settings['read_out']['meas_time']
#         green_time = self.settings['read_out']['green_time']
#         laser_off_time = self.settings['read_out']['laser_off_time']
#
#         for tau in tau_list:
#             pulse_sequence = \
#                 [Pulse('laser', laser_off_time, green_time),
#                  Pulse('apd_readout', laser_off_time - green_time, meas_time) # dark reference
#                  ]
#             pulse_sequence.append(Pulse('apd_readout', laser_off_time + tau, meas_time))
#             pulse_sequence.append(Pulse('apd_switch', 30, laser_off_time + green_time + green_time))
#
#             pulse_sequences.append(pulse_sequence)
#
#         print('number of sequences before validation ', len(pulse_sequences))
#         # return pulse_sequences, self.settings['num_averages'], tau_list, meas_time
#
#         return pulse_sequences, tau_list, meas_time
#
#
#     def _plot(self, axislist, data = None):
#         '''
#         Plot 1: self.data['tau'], the list of times specified for a given experiment, verses self.data['counts'], the data
#         received for each time
#         Plot 2: the pulse sequence performed at the current time (or if plotted statically, the last pulse sequence
#         performed
#
#         Args:
#             axes_list: list of axes to write plots to (uses first 2)
#             data (optional) dataset to plot (dictionary that contains keys counts, tau, fits), if not provided use self.data
#         '''
#
#
#
#         super(GrDelayMeas, self)._plot(axislist)
#         axislist[0].set_title('AOM delay calibration')
#
#         if data is None:
#             data = self.data
#
#         if 'exp_finished' in data.keys() is not None and data['exp_finished'] == 1:
#
#             axislist[0].plot(data['tau'], data['counts'])
#             axislist[0].set_title('(final plot) AOM delay calibration')
#
#         else:
#
#             super(GrDelayMeas, self)._plot(axislist)
#             axislist[0].set_title('(initial plot) AOM delay calibration')
#
#
#     def _update_plot(self, axislist):
#             if len(axislist[0].lines) == 0:
#                 self._plot(axislist)
#                 return
#             # axislist[0].lines = axislist[0].lines[0]
#             super(GrDelayMeas, self)._update_plot(axislist)
#
#             axislist[0].set_title('(updating plot) AOM delay calibration')

class GrDelayMeas(PulsedExperimentBaseScript):
    """
#         This script measures the green laser delay during AOM turn ON and turn OFF
#         ==> Last edited by ZQ 2/26/2019 7:05pm
#     """
    _DEFAULT_SETTINGS = [
        Parameter('tau_times', [
            Parameter('min_time', 20, float, 'minimum green delay (in ns)'),
            Parameter('max_time', 800, float, 'max green delay (in ns)'),
            Parameter('time_step', 10, [10, 20, 30, 40, 50, 80, 100], 'time step increment of green delay (in ns)')
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 30, int, '[ns] APD window to count photons'),
            Parameter('green_time', 500, float, '[ns] duration of green pulse '),
            Parameter('laser_off_time', 1000, float, '[ns] dark time before AOM turn on')
        ]),
        Parameter('num_averages', 100000, int, 'number of averages'),
        # Parameter('APD_type', 'SPCM', ['SPCM', 'MPPC'], 'SPCM = single-photon counting module (digital input on DAQ); MPPC = multi-pixel photon counter (analog input on DAQ)'),
        # Parameter('MPPC_ai_channel', 'ai4', ['ai0', 'ai1', 'ai2', 'ai3', 'ai4', 'ai5', 'ai6', 'ai7', 'ai8', 'ai9'],
        #           'NI DAQ analog input channel for MPPC voltage'),
        # Parameter('skip_invalid_sequences', True, bool, 'Skips any sequences with <MIN_DURATIONns commands')
    ]

    _INSTRUMENTS = {'daq': NI6353, 'PB': LISE607RTPulseBlaster, 'mw_gen': R8SMicrowaveGenerator}

    def _function(self):


        ## Turn off green light (the pulse blaster will pulse it on when needed)
        self.instruments['PB']['instance'].update({'laser': {'status': False}})


        super(GrDelayMeas, self)._function(self.data)

        ## Turn off green light (the pulse blaster will pulse it on when needed)
        self.instruments['PB']['instance'].update({'laser': {'status': False}})

        self.data[
            'exp_finished'] = 1  # this is to flag when the experiment is done or stopped and this allows proper saving even fit fails

        # turn off laser, apd_switch and MW (just in case)
        self.instruments['PB']['instance'].update({'laser': {'status': False}})
        self.instruments['PB']['instance'].update({'apd_switch': {'status': False}})

    def _create_pulse_sequences(self):
        '''

        Returns: pulse_sequences, num_averages, tau_list, meas_time
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
            meas_time: the width (in ns) of the daq measurement

        '''
        pulse_sequences = []

        # JG 16-08-25 changed (MIN_DURATION min spacing is taken care of later):
        # tau_list = range(0, int(self.settings['tau_times']['max_time']), self.settings['tau_times']['time_step'])
        # tau_list = range(int(self.settings['tau_times']['min_time']), int(self.settings['tau_times']['max_time']),self.settings['tau_times']['time_step'])

        max_range = int(np.floor((self.settings['tau_times']['max_time'] - self.settings['tau_times']['min_time']) /
                                 self.settings['tau_times']['time_step']))
        tau_list = np.array(
            [self.settings['tau_times']['min_time'] + i * self.settings['tau_times']['time_step'] for i in
             range(max_range)])

        # ignore the sequence if the mw-pulse is shorter than MIN_DURATION (0 is ok because there is no mw pulse!)
        tau_list = [x for x in tau_list if x == 0 or x >= MIN_DURATION]
        print('tau_list', tau_list)

        meas_time = self.settings['read_out']['meas_time']
        green_time = self.settings['read_out']['green_time']
        laser_off_time = self.settings['read_out']['laser_off_time']

        for tau in tau_list:
            pulse_sequence = \
                [Pulse('laser', laser_off_time, green_time),
                 Pulse('apd_readout', laser_off_time - green_time, meas_time) # dark reference
                 ]
            pulse_sequence.append(Pulse('apd_readout', laser_off_time + tau, meas_time))
            pulse_sequence.append(Pulse('apd_switch', 30, laser_off_time + green_time + green_time))

            pulse_sequences.append(pulse_sequence)

        print('number of sequences before validation ', len(pulse_sequences))
        # return pulse_sequences, self.settings['num_averages'], tau_list, meas_time

        return pulse_sequences, tau_list, meas_time

    def _plot(self, axislist, data = None):
        '''
        Plot 1: self.data['tau'], the list of times specified for a given experiment, verses self.data['counts'], the data
        received for each time
        Plot 2: the pulse sequence performed at the current time (or if plotted statically, the last pulse sequence
        performed

        Args:
            axes_list: list of axes to write plots to (uses first 2)
            data (optional) dataset to plot (dictionary that contains keys counts, tau, fits), if not provided use self.data
        '''

        if data is None:
            data = self.data

        if 'exp_finished' in data.keys() is not None and data['exp_finished'] == 1:

            axislist[0].plot(data['tau'], data['counts'][:, 0])
            axislist[0].plot(data['tau'], data['counts'][:, 1])

            axislist[0].set_xlabel('AOM delay [ns]')
            axislist[0].set_ylabel('kCounts/sec')

            axislist[0].set_title('(final plot) AOM delay calibration')
            axislist[0].legend(labels=('Dark Fluorescence', 'Data'), fontsize=8)

        else:

            super(GrDelayMeas, self)._plot(axislist)
            axislist[0].set_title('(initial plot) AOM delay calibration')
            axislist[0].legend(labels=('Dark Fluorescence', 'Data'), fontsize=8)

    def _update_plot(self, axislist):
            # self._plot(axislist)

            if len(axislist[0].lines) == 0:
                self._plot(axislist)
                return
            super(GrDelayMeas, self)._update_plot(axislist)

            axislist[0].set_title('(updating plot) AOM delay calibration')
            axislist[0].legend(labels=('Dark Fluorescence', 'Data'), fontsize=8)


class IQCalibration_N9310A(PulsedExperimentBaseScript):
    """
    This script calibrates the MW amplitude and phase between I and Q using a modified Hahn-echo sequence.
    Keysight N9310A generator is used and it has IQ modulation.
    ==> last edited by Ziwei Qiu on 4/18/2019
    """

    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_power', -10.0, float, 'microwave power in dB'),
            Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
            Parameter('microwave_channel', 'i', ['i', 'q'], 'Channel to use for mw pulses'),
            Parameter('pi_pulse_time', 50.0, float, 'time duration of a pi pulse (in ns)'),
            Parameter('pi_half_pulse_time', 25.0, float, 'time duration of a pi/2 pulse (in ns)'),
            Parameter('3pi_half_pulse_time', 75.0, float, 'time duration of a 3pi/2 pulse (in ns)')
        ]),
        Parameter('tau_times', [
            Parameter('min_time', 200, float, 'minimum time between pi pulses'),
            Parameter('max_time', 10000, float, 'maximum time between pi pulses'),
            Parameter('time_step', 100.,
                      [2.5, 5., 10., 20., 50., 100., 200., 300., 400., 500., 600., 800., 1000., 2000., 10000., 100000.,
                       500000.],
                      'time step increment of time between pi pulses (in ns)')
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 100, float, 'measurement time after rabi sequence (in ns)'),
            Parameter('nv_reset_time', 2000, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 600, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 540, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('min_pulse_dur', 15, [10, 15, 50], 'Minimum allowed pulse duration (ns)'),
        Parameter('num_averages', 500000, int, 'number of averages (>100000)'),
    ]

    _INSTRUMENTS = {'daq': NI6353, 'PB': LISE607RTPulseBlaster, 'mw_gen_a': AgilentMicrowaveGenerator}

    def _function(self):

        ## Turn off green light (the pulse blaster will pulse it on when needed)
        self.instruments['PB']['instance'].update({'laser': {'status': False}})
        self.instruments['PB']['instance'].update({'apd_switch': {'status': False}})

        # self.instruments['mw_gen_a']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen_a']['instance'].update({'enable_IQ': True})
        self.instruments['mw_gen_a']['instance'].update({'enable_modulation': True})
        self.instruments['mw_gen_a']['instance'].update({'freq_mode': 'CW'})
        self.instruments['mw_gen_a']['instance'].update({'power_mode': 'CW'})
        self.instruments['mw_gen_a']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power']})
        self.instruments['mw_gen_a']['instance'].update({'frequency': self.settings['mw_pulses']['mw_frequency']})
        self.instruments['mw_gen_a']['instance'].update({'enable_output': True})

        self.avg_block_number = super(IQCalibration_N9310A, self)._function()

        # this is to flag when the experiment is done or stopped and this allows proper saving even fit fails
        self.data['exp_finished'] = 1

        # turn off laser, apd_switch and MW (just in case)
        self.instruments['PB']['instance'].update({'laser': {'status': False}})
        self.instruments['PB']['instance'].update({'apd_switch': {'status': False}})
        self.instruments['mw_gen_a']['instance'].update({'enable_output': False})
        self.instruments['mw_gen_a']['instance'].update({'enable_IQ': False})
        self.instruments['mw_gen_a']['instance'].update({'enable_modulation': False})

        if 'counts' in self.data.keys() and 'tau' in self.data.keys() and 0 not in (
                self.data['counts'][:, 1] + self.data['counts'][:, 0]) and 0 not in (
                self.data['counts'][:, 2] + self.data['counts'][:, 3]) and 0 not in (
                self.data['counts'][:, 5] + self.data['counts'][:, 4]) and 0 not in (
                self.data['counts'][:, 6] + self.data['counts'][:, 7]):

            self.data['XYX'] = 2. * (self.data['counts'][:, 1] - self.data['counts'][:, 0]) / (
                        self.data['counts'][:, 1] + self.data['counts'][:, 0])
            self.data['XXX'] = 2. * (self.data['counts'][:, 2] - self.data['counts'][:, 3]) / (
                        self.data['counts'][:, 2] + self.data['counts'][:, 3])
            self.data['YXY'] = 2. * (self.data['counts'][:, 5] - self.data['counts'][:, 4]) / (
                        self.data['counts'][:, 5] + self.data['counts'][:, 4])
            self.data['YYY'] = 2. * (self.data['counts'][:, 6] - self.data['counts'][:, 7]) / (
                        self.data['counts'][:, 6] + self.data['counts'][:, 7])

    def _create_pulse_sequences(self):
        '''
        Returns: pulse_sequences, num_averages, tau_list, meas_time
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
            meas_time: the width (in ns) of the daq measurement
        '''
        MIN_DURATION = self.settings['min_pulse_dur']
        pulse_sequences = []

        max_range = int(np.floor((self.settings['tau_times']['max_time'] - self.settings['tau_times']['min_time']) /
                                 self.settings['tau_times']['time_step']))
        tau_list = np.array(
            [self.settings['tau_times']['min_time'] + i * self.settings['tau_times']['time_step'] for i in
             range(max_range)])

        # ignore the sequence if the mw-pulse is shorter than 15ns (0 is ok because there is no mw pulse!)
        tau_list = [x for x in tau_list if x == 0 or x >= MIN_DURATION]

        nv_reset_time = self.settings['read_out']['nv_reset_time']
        delay_readout = self.settings['read_out']['delay_readout']
        if self.settings['mw_pulses']['microwave_channel'] == 'i':
            microwave_channel = 'microwave_i'
            microwave_channel_2 = 'microwave_q'
        else:
            microwave_channel = 'microwave_q'
            microwave_channel_2 = 'microwave_i'
        pi_time = self.settings['mw_pulses']['pi_pulse_time']
        pi_half_time = self.settings['mw_pulses']['pi_half_pulse_time']
        three_pi_half_time = self.settings['mw_pulses']['3pi_half_pulse_time']

        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']


        for tau_total in tau_list:
            tau = tau_total
            pulse_sequence = []

            # Two XYX:
            pulse_sequence = [Pulse(microwave_channel, laser_off_time, pi_half_time)]
            pulse_sequence.append(
                Pulse(microwave_channel_2, laser_off_time + pi_half_time / 2. + tau - pi_time / 2., pi_time))
            pulse_sequence.append(
                Pulse(microwave_channel, laser_off_time + pi_half_time / 2. + tau + tau - pi_half_time / 2.,
                      pi_half_time))

            end_of_first_XYX = laser_off_time + pi_half_time/2. + tau + tau - pi_half_time/2. + pi_half_time

            pulse_sequence.append(Pulse('laser', end_of_first_XYX + delay_mw_readout, nv_reset_time))
            pulse_sequence.append(Pulse('apd_readout', end_of_first_XYX + delay_mw_readout + delay_readout, meas_time))

            start_of_second_XYX = end_of_first_XYX + delay_mw_readout + nv_reset_time + laser_off_time

            pulse_sequence.append(Pulse(microwave_channel, start_of_second_XYX, pi_half_time))
            pulse_sequence.append(Pulse(microwave_channel_2, start_of_second_XYX + pi_half_time/2. + tau - pi_time/2., pi_time))
            pulse_sequence.append(Pulse(microwave_channel, start_of_second_XYX + pi_half_time/2. + tau + tau - pi_half_time/2., three_pi_half_time))

            end_of_second_XYX = start_of_second_XYX + pi_half_time/2. + tau + tau - pi_half_time/2. + three_pi_half_time

            pulse_sequence.append(Pulse('laser', end_of_second_XYX + delay_mw_readout, nv_reset_time))
            pulse_sequence.append(Pulse('apd_readout', end_of_second_XYX + delay_mw_readout + delay_readout, meas_time))

            #Two XXX:

            start_of_first_XXX = end_of_second_XYX + delay_mw_readout + nv_reset_time + laser_off_time

            pulse_sequence.append(Pulse(microwave_channel, start_of_first_XXX, pi_half_time))
            pulse_sequence.append(
                Pulse(microwave_channel, start_of_first_XXX + pi_half_time / 2. + tau - pi_time / 2., pi_time))
            pulse_sequence.append(
                Pulse(microwave_channel, start_of_first_XXX + pi_half_time / 2. + tau + tau - pi_half_time / 2.,
                      pi_half_time))

            end_of_first_XXX =  start_of_first_XXX + pi_half_time/2. + tau + tau - pi_half_time/2. + pi_half_time

            pulse_sequence.append(Pulse('laser', end_of_first_XXX + delay_mw_readout, nv_reset_time))
            pulse_sequence.append(Pulse('apd_readout', end_of_first_XXX + delay_mw_readout + delay_readout, meas_time))

            start_of_second_XXX = end_of_first_XXX + delay_mw_readout + nv_reset_time + laser_off_time

            pulse_sequence.append(Pulse(microwave_channel, start_of_second_XXX, pi_half_time))
            pulse_sequence.append(
                Pulse(microwave_channel, start_of_second_XXX + pi_half_time / 2. + tau - pi_time / 2., pi_time))
            pulse_sequence.append(
                Pulse(microwave_channel, start_of_second_XXX + pi_half_time / 2. + tau + tau - pi_half_time / 2.,
                      three_pi_half_time))

            end_of_second_XXX = start_of_second_XXX + pi_half_time / 2. + tau + tau - pi_half_time / 2. + three_pi_half_time

            pulse_sequence.append(Pulse('laser', end_of_second_XXX + delay_mw_readout, nv_reset_time))
            pulse_sequence.append(Pulse('apd_readout', end_of_second_XXX + delay_mw_readout + delay_readout, meas_time))

            # Two YXY:
            start_of_first_YXY = end_of_second_XXX + delay_mw_readout + nv_reset_time + laser_off_time

            pulse_sequence.append(Pulse(microwave_channel_2, start_of_first_YXY, pi_half_time))
            pulse_sequence.append(
                Pulse(microwave_channel, start_of_first_YXY + pi_half_time / 2. + tau - pi_time / 2., pi_time))
            pulse_sequence.append(
                Pulse(microwave_channel_2, start_of_first_YXY + pi_half_time / 2. + tau + tau - pi_half_time / 2.,
                      pi_half_time))

            end_of_first_YXY = start_of_first_YXY + pi_half_time / 2. + tau + tau - pi_half_time / 2. + pi_half_time

            pulse_sequence.append(Pulse('laser', end_of_first_YXY + delay_mw_readout, nv_reset_time))
            pulse_sequence.append(Pulse('apd_readout', end_of_first_YXY + delay_mw_readout + delay_readout, meas_time))

            start_of_second_YXY = end_of_first_YXY + delay_mw_readout + nv_reset_time + laser_off_time

            pulse_sequence.append(Pulse(microwave_channel_2, start_of_second_YXY, pi_half_time))
            pulse_sequence.append(
                Pulse(microwave_channel, start_of_second_YXY + pi_half_time / 2. + tau - pi_time / 2., pi_time))
            pulse_sequence.append(
                Pulse(microwave_channel_2, start_of_second_YXY + pi_half_time / 2. + tau + tau - pi_half_time / 2.,
                      three_pi_half_time))

            end_of_second_YXY = start_of_second_YXY + pi_half_time / 2. + tau + tau - pi_half_time / 2. + three_pi_half_time

            pulse_sequence.append(Pulse('laser', end_of_second_YXY + delay_mw_readout, nv_reset_time))
            pulse_sequence.append(Pulse('apd_readout', end_of_second_YXY + delay_mw_readout + delay_readout, meas_time))

            # Two YYY:

            start_of_first_YYY = end_of_second_YXY + delay_mw_readout + nv_reset_time + laser_off_time

            pulse_sequence.append(Pulse(microwave_channel_2, start_of_first_YYY, pi_half_time))
            pulse_sequence.append(
                Pulse(microwave_channel_2, start_of_first_YYY + pi_half_time / 2. + tau - pi_time / 2., pi_time))
            pulse_sequence.append(
                Pulse(microwave_channel_2, start_of_first_YYY + pi_half_time / 2. + tau + tau - pi_half_time / 2.,
                      pi_half_time))

            end_of_first_YYY = start_of_first_YYY + pi_half_time / 2. + tau + tau - pi_half_time / 2. + pi_half_time

            pulse_sequence.append(Pulse('laser', end_of_first_YYY + delay_mw_readout, nv_reset_time))
            pulse_sequence.append(Pulse('apd_readout', end_of_first_YYY + delay_mw_readout + delay_readout, meas_time))

            start_of_second_YYY = end_of_first_YYY + delay_mw_readout + nv_reset_time + laser_off_time

            pulse_sequence.append(Pulse(microwave_channel_2, start_of_second_YYY, pi_half_time))
            pulse_sequence.append(
                Pulse(microwave_channel_2, start_of_second_YYY + pi_half_time / 2. + tau - pi_time / 2., pi_time))
            pulse_sequence.append(
                Pulse(microwave_channel_2, start_of_second_YYY + pi_half_time / 2. + tau + tau - pi_half_time / 2.,
                      three_pi_half_time))

            end_of_second_YYY = start_of_second_YYY + pi_half_time / 2. + tau + tau - pi_half_time / 2. + three_pi_half_time

            pulse_sequence.append(Pulse('laser', end_of_second_YYY + delay_mw_readout, nv_reset_time))
            pulse_sequence.append(Pulse('apd_readout', end_of_second_YYY + delay_mw_readout + delay_readout, meas_time))

            # Turn on APD switch all the time
            pulse_sequence.append(Pulse('apd_switch', 0,
                                        end_of_second_YYY + delay_mw_readout + nv_reset_time))

            # Append this singlge sequence to the total sequnce
            pulse_sequences.append(pulse_sequence)

        print('number of sequences before validation ', len(pulse_sequences))
        return pulse_sequences, tau_list, meas_time

    def _plot(self, axislist, data=None):

        if data is None:
            data = self.data
        if 'exp_finished' in data.keys() is not None and data['exp_finished'] == 1:
            if 'XYX' in data.keys() is not None and 'YXY' in data.keys() is not None and 'XXX' in data.keys() is not None and 'YYY' in data.keys() is not None:

                axislist[0].plot(data['tau'], data['XYX'])
                axislist[0].plot(data['tau'], data['XXX'])
                axislist[0].plot(data['tau'], data['YXY'])
                axislist[0].plot(data['tau'], data['YYY'])
                axislist[0].set_ylabel('normalized fluorescence')
                axislist[0].set_xlabel('tau [ns]')

                axislist[0].set_title('(final plot) Keysight N9310A - IQ_calibration \n mw-power:{:.2f}dBm, mw_freq:{:.3f} GHz'.format(
                    self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency'] * 1e-9))
                axislist[0].legend(labels=('XYX', 'XXX', 'YXY', 'YYY'), fontsize=10)
            else:
                axislist[0].plot(data['tau'], data['counts'][:, 0])
                axislist[0].plot(data['tau'], data['counts'][:, 1])
                axislist[0].plot(data['tau'], data['counts'][:, 2])
                axislist[0].plot(data['tau'], data['counts'][:, 3])
                axislist[0].plot(data['tau'], data['counts'][:, 4])
                axislist[0].plot(data['tau'], data['counts'][:, 5])
                axislist[0].plot(data['tau'], data['counts'][:, 6])
                axislist[0].plot(data['tau'], data['counts'][:, 7])
                axislist[0].set_ylabel('fluorescence [kcps]')
                axislist[0].set_xlabel('tau [ns]')
                axislist[0].legend(labels=('XYX_pi/2 {:.0f}kcps'.format(np.mean(data['counts'][:, 0])),
                                           'XYX_3pi/2 {:.0f}kcps'.format(np.mean(data['counts'][:, 1])),
                                           'XXX_pi/2 {:.0f}kcps'.format(np.mean(data['counts'][:, 2])),
                                           'XXX_3pi/2 {:.0f}kcps'.format(np.mean(data['counts'][:, 3])),
                                           'YXY_pi/2 {:.0f}kcps'.format(np.mean(data['counts'][:, 4])),
                                           'YXY_3pi/2 {:.0f}kcps'.format(np.mean(data['counts'][:, 5])),
                                           'YYY_pi/2 {:.0f}kcps'.format(np.mean(data['counts'][:, 6])),
                                           'YYY_3pi/2 {:.0f}kcps'.format(np.mean(data['counts'][:, 7]))), fontsize=8)

                axislist[0].set_title('(final plot) Keysight N9310A - IQ_calibration \n mw-power:{:.2f}dBm, mw_freq:{:.3f} GHz'.format(
                    self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency'] * 1e-9))

        else:
            super(IQCalibration_N9310A, self)._plot(axislist)
            axislist[0].set_ylabel('fluorescence [kcps]')
            axislist[0].set_xlabel('tau [ns]')
            axislist[0].legend(labels=('XYX_pi/2', 'XYX_3pi/2', 'XXX_pi/2 ', 'XXX_3pi/2 ', 'YXY_pi/2 ', 'YXY_3pi/2 ', 'YYY_pi/2 ', 'YYY_3pi/2 '), fontsize = 8)
            axislist[0].set_title(
                '(initial plot) Keysight N9310A - IQ_calibration \n mw-power:{:.2f}dBm, mw_freq:{:.3f} GHz'.format(
                    self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency'] * 1e-9))

    def _update_plot(self, axislist):

        if len(axislist[0].lines) == 0:
            self._plot(axislist)
            return
        super(IQCalibration_N9310A, self)._update_plot(axislist)

        data = self.data
        axislist[0].set_title(
            '(updating plot) Keysight N9310A - IQ_calibration \n mw-power:{:.2f}dBm, mw_freq:{:.3f} GHz'.format(
                self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency'] * 1e-9))
        axislist[0].set_ylabel('fluorescence [kcps]')
        axislist[0].set_xlabel('tau [ns]')
        axislist[0].legend(labels=('XYX_pi/2 {:.0f}kcps'.format(np.mean(data['counts'][:, 0])),
                                   'XYX_3pi/2 {:.0f}kcps'.format(np.mean(data['counts'][:, 1])),
                                   'XXX_pi/2 {:.0f}kcps'.format(np.mean(data['counts'][:, 2])),
                                   'XXX_3pi/2 {:.0f}kcps'.format(np.mean(data['counts'][:, 3])),
                                   'YXY_pi/2 {:.0f}kcps'.format(np.mean(data['counts'][:, 4])),
                                   'YXY_3pi/2 {:.0f}kcps'.format(np.mean(data['counts'][:, 5])),
                                   'YYY_pi/2 {:.0f}kcps'.format(np.mean(data['counts'][:, 6])),
                                   'YYY_3pi/2 {:.0f}kcps'.format(np.mean(data['counts'][:, 7]))), fontsize=8)
# class ReadoutStartTime(PulsedExperimentBaseScript):  # ER 10.21.2017
#     """
#     This script sweeps the start time of the APD readout window, using a fixed-duration readout pulse.
#     The goal is to figure out when to turn on the readout window to maximize SNR.
#     To maximize sensitivity and protect against laser power drift noise, we use the 'double-init' scheme.
#     """
#     _DEFAULT_SETTINGS = [
#         Parameter('mw_pulse', [
#             Parameter('mw_power', -45.0, float, 'microwave power in dB'),
#             Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
#             Parameter('microwave_channel', 'i', ['i', 'q'], 'Channel to use for mw pulses'),
#             Parameter('pi_time', 30.0, float, 'pi time in ns')
#         ]),
#         Parameter('tau_times', [
#             Parameter('min_time', 15, float, 'minimum time for rabi oscillations (in ns)'),
#             Parameter('max_time', 200, float, 'total time of rabi oscillations (in ns)'),
#             Parameter('time_step', 5, [5, 10, 20, 50, 100, 200, 500, 1000, 10000, 100000, 500000],
#                       'time step increment of readout pulse duration (in ns)')
#         ]),
#         Parameter('read_out', [
#             Parameter('nv_reset_time', 1750, int, 'time with laser on to reset state'),
#             Parameter('laser_off_time', 1000, int,
#                       'minimum laser off time before taking measurements (ns)'),
#             Parameter('delay_mw_readout', 100, int, 'delay between mw and readout (in ns)'),
#             Parameter('delay_readout', 30, int, 'delay between laser on and readout (given by spontaneous decay rate)'),
#             Parameter('readout_window', 300, int, 'length of readout window')
#         ]),
#         Parameter('num_averages', 100000, int, 'number of averages'),
#     ]
#
#     _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator}
#
#     def _function(self):
#         # COMMENT_ME
#
#         self.data['fits'] = None
#         self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
#         self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulse']['mw_power']})
#         self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulse']['mw_frequency']})
#         super(ReadoutStartTime, self)._function(self.data)
#
#
#         counts = self.data['counts'][:, 1] # / self.data['counts'][:, 0]
#         tau = self.data['tau']
#
#         try:
#             fits = fit_exp_decay(tau, counts, varibale_phase=True)
#             self.data['fits'] = fits
#         except:
#             self.data['fits'] = None
#             self.log('fit failed')
#
#     def _create_pulse_sequences(self):
#         '''
#
#         Returns: pulse_sequences, num_averages, tau_list, meas_time
#             pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
#             scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
#             sequence must have the same number of daq read pulses
#             num_averages: the number of times to repeat each pulse sequence
#             tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
#             meas_time: the width (in ns) of the daq measurement
#
#         '''
#         pulse_sequences = []
#         # tau_list = range(int(max(15, self.settings['tau_times']['time_step'])), int(self.settings['tau_times']['max_time'] + 15),
#         #                  self.settings['tau_times']['time_step'])
#         # JG 16-08-25 changed (15ns min spacing is taken care of later):
#         tau_list = list(range(int(self.settings['tau_times']['min_time']), int(self.settings['tau_times']['max_time']),
#                          self.settings['tau_times']['time_step']))
#
#         # ignore the sequence if the mw-pulse is shorter than 15ns (0 is ok because there is no mw pulse!)
#         tau_list = [x for x in tau_list if x == 0 or x >= 15]
#         nv_reset_time = self.settings['read_out']['nv_reset_time']
#         delay_readout = self.settings['read_out']['delay_readout']
#         microwave_channel = 'microwave_' + self.settings['mw_pulse']['microwave_channel']
#
#         laser_off_time = self.settings['read_out']['laser_off_time']
#         delay_mw_readout = self.settings['read_out']['delay_mw_readout']
#         pi_time = self.settings['mw_pulse']['pi_time']
#         meas_time = self.settings['read_out']['readout_window']
#
#
#         for tau in tau_list:
#             pulse_sequence = \
#                 [Pulse('laser', laser_off_time, nv_reset_time),
#                  Pulse('apd_readout', laser_off_time+ delay_readout + tau, meas_time),
#                  ]
#             # if tau is 0 there is actually no mw pulse
#             if tau > 0:
#                 pulse_sequence += [Pulse(microwave_channel, laser_off_time + nv_reset_time + laser_off_time, pi_time)]
#
#             pulse_sequence += [
#                 Pulse('laser', laser_off_time + nv_reset_time + laser_off_time + delay_mw_readout, nv_reset_time),
#                 Pulse('apd_readout', laser_off_time + nv_reset_time + laser_off_time + delay_mw_readout + delay_readout + tau, meas_time)
#             ]
#             # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
#             # if tau == 0 or tau>=15:
#             pulse_sequences.append(pulse_sequence)
#
#         return pulse_sequences, tau_list, meas_time
#
#     def _plot(self, axislist, data=None):
#         '''
#         Plot 1: self.data['tau'], the list of times specified for a given experiment, verses self.data['counts'], the data
#         received for each time
#         Plot 2: the pulse sequence performed at the current time (or if plotted statically, the last pulse sequence
#         performed
#
#         Args:
#             axes_list: list of axes to write plots to (uses first 2)
#             data (optional) dataset to plot (dictionary that contains keys counts, tau, fits), if not provided use self.data
#         '''
#
#         if data is None:
#             data = self.data
#
#         if data['fits'] is not None:
#             counts = data['counts'][:, 1] / data['counts'][:, 0]
#             tau = data['tau']
#             fits = data['fits']  # amplitude, frequency, phase, offset
#
#             axislist[0].plot(tau, counts, 'b')
#             axislist[0].hold(True)
#
#             axislist[0].plot(tau, cose_with_decay(tau, *fits), 'k', lw=3)
#             # pi_time = 2*np.pi / fits[1] / 2
#             pi_time = (np.pi - fits[2]) / fits[1]
#             pi_half_time = (np.pi / 2 - fits[2]) / fits[1]
#             three_pi_half_time = (3 * np.pi / 2 - fits[2]) / fits[1]
#             rabi_freq = 1000 * fits[1] / (2 * np.pi)
#             #   axislist[0].set_title('Rabi mw-power:{:0.1f}dBm, mw_freq:{:0.3f} GHz, pi-time: {:2.1f}ns, pi-half-time: {:2.1f}ns, 3pi_half_time: {:2.1f}ns, Rabi freq: {2.1f}MHz'.format(self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency']*1e-9, pi_time, pi_half_time, three_pi_half_time, rabi_freq))
#             axislist[0].set_title('Readout pulse width counts')
#         else:
#             super(ReadoutStartTime, self)._plot(axislist)
#             axislist[0].set_title('Readout pulse width counts')
#             axislist[0].legend(labels=('Ref Fluorescence', 'Pi pulse Data'), fontsize=8)
#
# class ReadoutDuration(PulsedExperimentBaseScript):
#     """
#   This script sweeps the readout pulse duration. To symmetrize the sequence between the 0 and +/-1 state we reinitialize every time
#       """
#     _DEFAULT_SETTINGS = [
#         Parameter('mw_pulse', [
#             Parameter('mw_power', -45.0, float, 'microwave power in dB'),
#             Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
#             Parameter('microwave_channel', 'i', ['i', 'q'], 'Channel to use for mw pulses'),
#             Parameter('pi_time', 30.0, float, 'pi time in ns')
#         ]),
#         Parameter('tau_times', [
#             Parameter('min_time', 15, float, 'minimum time for rabi oscillations (in ns)'),
#             Parameter('max_time', 200, float, 'total time of rabi oscillations (in ns)'),
#             Parameter('time_step', 5, [5, 10, 20, 50, 100, 200, 500, 1000, 10000, 100000, 500000],
#                       'time step increment of readout pulse duration (in ns)')
#         ]),
#         Parameter('read_out', [
#             Parameter('nv_reset_time', 1750, int, 'time with laser on to reset state'),
#             Parameter('laser_off_time', 1000, int,
#                       'minimum laser off time before taking measurements (ns)'),
#             Parameter('delay_mw_readout', 100, int, 'delay between mw and readout (in ns)'),
#             Parameter('delay_readout', 30, int, 'delay between laser on and readout (given by spontaneous decay rate)')
#         ]),
#         Parameter('num_averages', 100000, int, 'number of averages'),
#     ]
#
#     _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator}
#
#     def _function(self):
#         # COMMENT_ME
#
#         self.data['fits'] = None
#         self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
#         self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulse']['mw_power']})
#         self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulse']['mw_frequency']})
#         super(ReadoutDuration, self)._function(self.data)
#
#         counts = self.data['counts'][:, 1]  # / self.data['counts'][:, 0]
#         tau = self.data['tau']
#
#         try:
#             fits = fit_exp_decay(tau, counts, varibale_phase=True)
#             self.data['fits'] = fits
#         except:
#             self.data['fits'] = None
#             self.log('fit failed')
#
#     def _create_pulse_sequences(self):
#         '''
#         Returns: pulse_sequences, num_averages, tau_list, meas_time
#             pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
#             scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
#             sequence must have the same number of daq read pulses
#             num_averages: the number of times to repeat each pulse sequence
#             tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
#             meas_time: the width (in ns) of the daq measurement
#         '''
#         pulse_sequences = []
#         # tau_list = range(int(max(15, self.settings['tau_times']['time_step'])), int(self.settings['tau_times']['max_time'] + 15),
#         #                  self.settings['tau_times']['time_step'])
#         # JG 16-08-25 changed (15ns min spacing is taken care of later):
#         tau_list = range(int(self.settings['tau_times']['min_time']), int(self.settings['tau_times']['max_time']),
#                          self.settings['tau_times']['time_step'])
#
#         # ignore the sequence if the mw-pulse is shorter than 15ns (0 is ok because there is no mw pulse!)
#         tau_list = [x for x in tau_list if x == 0 or x >= 15]
#
#         nv_reset_time = self.settings['read_out']['nv_reset_time']
#         delay_readout = self.settings['read_out']['delay_readout']
#         microwave_channel = 'microwave_' + self.settings['mw_pulse']['microwave_channel']
#
#         laser_off_time = self.settings['read_out']['laser_off_time']
#         delay_mw_readout = self.settings['read_out']['delay_mw_readout']
#         pi_time = self.settings['mw_pulse']['pi_time']
#         meas_time = 100
#
#         for tau in tau_list:
#             pulse_sequence = \
#                 [Pulse('laser', laser_off_time, nv_reset_time),
#                  Pulse('apd_readout', laser_off_time + delay_readout, tau),
#                  ]
#             # if tau is 0 there is actually no mw pulse
#             if tau > 0:
#                 pulse_sequence += [Pulse(microwave_channel, laser_off_time + nv_reset_time + laser_off_time, pi_time)]
#
#             pulse_sequence += [
#                 Pulse('laser', laser_off_time + nv_reset_time + laser_off_time + delay_mw_readout, nv_reset_time),
#                 Pulse('apd_readout', laser_off_time + nv_reset_time + laser_off_time + delay_mw_readout + delay_readout,
#                       tau)
#             ]
#             # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
#             # if tau == 0 or tau>=15:
#             pulse_sequences.append(pulse_sequence)
#
#         return pulse_sequences, tau_list, meas_time
#
#     def _plot(self, axislist, data=None):
#         """
#         Plot 1: self.data['tau'], the list of times specified for a given experiment, verses self.data['counts'], the data
#         received for each time
#         Plot 2: the pulse sequence performed at the current time (or if plotted statically, the last pulse sequence
#         performed
#         Args:
#             axes_list: list of axes to write plots to (uses first 2)
#             data (optional) dataset to plot (dictionary that contains keys counts, tau, fits), if not provided use self.data
#         """
#
#         if data is None:
#             data = self.data
#
#         if data['fits'] is not None:
#             counts = data['counts'][:, 1] / data['counts'][:, 0]
#             tau = data['tau']
#             fits = data['fits']  # amplitude, frequency, phase, offset
#
#             axislist[0].plot(tau, counts, 'b')
#             axislist[0].hold(True)
#
#             axislist[0].plot(tau, cose_with_decay(tau, *fits), 'k', lw=3)
#             # pi_time = 2*np.pi / fits[1] / 2
#             pi_time = (np.pi - fits[2]) / fits[1]
#             pi_half_time = (np.pi / 2 - fits[2]) / fits[1]
#             three_pi_half_time = (3 * np.pi / 2 - fits[2]) / fits[1]
#             rabi_freq = 1000 * fits[1] / (2 * np.pi)
#             #   axislist[0].set_title('Rabi mw-power:{:0.1f}dBm, mw_freq:{:0.3f} GHz, pi-time: {:2.1f}ns, pi-half-time: {:2.1f}ns, 3pi_half_time: {:2.1f}ns, Rabi freq: {2.1f}MHz'.format(self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency']*1e-9, pi_time, pi_half_time, three_pi_half_time, rabi_freq))
#             axislist[0].set_title('Readout pulse width counts')
#         else:
#             super(ReadoutDuration, self)._plot(axislist)
#             axislist[0].set_title('Readout pulse width counts')
#             axislist[0].legend(labels=('Ref Fluorescence', 'Pi pulse Data'), fontsize=8)
#
# class ReadoutStartTimeWithoutMW(PulsedExperimentBaseScript):
#     """
#
#     This script sweeps the start time of the APD readout window, using a fixed-duration readout pulse.
#     The goal is to figure out when to turn on the readout window to maximize SNR.
#     Does not include a MW pulse (so no double-init scheme)
#
#     """
#     _DEFAULT_SETTINGS = [
#         Parameter('count_source_pulse_width', 10000, int, 'How long to pulse the count source (in ns)'),
#         Parameter('measurement_gate_pulse_width', 15, int, 'How long to have the DAQ acquire data (in ns)'),
#         Parameter('min_delay', 0, int, 'minimum delay over which to scan'),
#         Parameter('max_delay', 1000, int, 'maximum delay over which to scan'),
#         Parameter('delay_interval_step_size', 15, int, 'Amount delay is increased for each new run'),
#         Parameter('num_averages', 1000, int, 'number of times to average for each delay'),
#         Parameter('reset_time', 10000, int, 'How long to wait for laser to turn off and reach steady state'),
#     ]
#
#     def _create_pulse_sequences(self):
#         '''
#         Creates a pulse sequence with no pulses for a reset time, then the laser on for a count_source_pulse_width time.
#         The daq measurement window is then swept across the laser pulse window to measure counts from min_delay
#         (can be negative) to max_delay, which are both defined with the laser on at 0.
#
#         Returns: pulse_sequences, num_averages, tau_list, measurement_gate_width
#             pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
#             scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
#             sequence must have the same number of daq read pulses
#             num_averages: the number of times to repeat each pulse sequence
#             tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
#             measurement_gate_width: the width (in ns) of the daq measurement
#
#
#         '''
#         pulse_sequences = []
#         gate_delays = list(range(self.settings['min_delay'], self.settings['max_delay'], self.settings['delay_interval_step_size']))
#         reset_time = self.settings['reset_time']
#         for delay in gate_delays:
#             pulse_sequences.append([Pulse('laser', reset_time, self.settings['count_source_pulse_width']),
#                                     Pulse('apd_readout', delay + reset_time,
#                                            self.settings['measurement_gate_pulse_width'])
#                                     ])
#         return pulse_sequences, gate_delays, self.settings['measurement_gate_pulse_width']
#
#     def _plot(self, axes_list, data=None):
#         """
#         Very similar to plot for PulseBlasterBaseScript but here deals with case where first plot has counts=[], using
#         PulseBlasterBaseScript plot leads to errors
#         Plot 1: self.data['tau'], the list of times specified for a given experiment, verses self.data['counts'], the data
#         received for each time
#         Plot 2: the pulse sequence performed at the current time (or if plotted statically, the last pulse sequence
#         performed
#
#         Args:
#             axes_list: list of axes to write plots to (uses first 2)
#             data (optional) dataset to plot (dictionary that contains keys counts, tau), if not provided use self.data
#         """
#         if data is None:
#             data = self.data
#
#         counts = data['counts']
#         x_data = data['tau']
#         axis1 = axes_list[0]
#         plot_1d_simple_timetrace_ns(axis1, x_data, [counts])
#         axis2 = axes_list[1]
#         plot_pulses(axis2, self.pulse_sequences[self.sequence_index])



if __name__ == '__main__':

    from pylabcontrol.core.scripts import Script # import script, AS and ER 20180426

    script = {}
    instr = {}
    script, failed, instr = Script.load_and_append({'IQCalibration_N9310A': 'IQCalibration_N9310A'}, script, instr)

    print(script)
    print(('failed', failed))
    print(instr)

