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
from b26_toolkit.instruments import NI6353, LISE607RTPulseBlaster, AgilentMicrowaveGenerator, R8SMicrowaveGenerator, Pulse, AgilentMicrowaveGeneratorII, Agilent33120A
from b26_toolkit.plotting.plots_1d import plot_pulses, update_pulse_plot, plot_1d_simple_timetrace_ns, update_1d_simple
from pylabcontrol.core import Parameter
from b26_toolkit.data_processing.fit_functions import fit_rabi_decay, cose_with_decay
import time

# MIN_DURATION = 15 # minimum allowed duration in [ns] for pulse blaster pulses

class Rabi_RnS(PulsedExperimentBaseScript): # ER 5.25.2017
    """
        This script applies a microwave pulse at fixed power for varying durations to measure Rabi oscillations.
        Uses a double_init scheme.

        The script uses Rhode&Schwarz signal generator SMB100A.
        Note that the MW generator used in this script doesn't have I/Q modulation option.

        ==> Last edited by Ziwei 4/17/2019

    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_power', -10.0, float, 'microwave power in dBm'),
            Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
            # Parameter('microwave_channel', 'i', ['i', 'q'], 'Channel to use for mw pulses')
            Parameter('add_initial_pulse', False, bool, 'check to add an initial pulse in the beginning'),
            Parameter('pulse_time', 50.0, float, 'time duration of an initial pulse (in ns)')
        ]),
        Parameter('tau_times', [
            Parameter('min_time', 15, float, 'minimum time for rabi oscillations (in ns)'),
            Parameter('max_time', 1000, float, 'total time of rabi oscillations (in ns)'),
            Parameter('time_step', 20., [2.5, 5., 10., 15., 20., 50., 100., 200., 500., 1000., 10000., 100000., 500000.],
                  'time step increment of rabi pulse duration (in ns)')
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 235, float, 'measurement time after rabi sequence (in ns)'),
            Parameter('nv_reset_time', 2000, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 600, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 550, int, '[ns] delay between laser on and APD readout (given by spontaneous decay rate)'),
            # Parameter('APD_switch_delay', 100, int, '[ns] delay between laser on and APD switch'),
            # Parameter('mw_switch_extra_time', 0, [0, 10, 15, 20, 30, 40],
            #           '[ns] buffer time of the MW switch window on both sides of MW_i or MW_q pulses. choose 0 if you are only using MW switches to control pulses')
        ]),
        Parameter('min_pulse_dur', 15, [10, 15, 50], 'Minimum allowed pulse duration (ns)'),
        Parameter('num_averages', 500000, int, 'number of averages (>100000)'),
    ]

    _INSTRUMENTS = {'daq': NI6353, 'PB': LISE607RTPulseBlaster, 'mw_gen': R8SMicrowaveGenerator}

    def _function(self):


        ## Turn off green light (the pulse blaster will pulse it on when needed)
        self.instruments['PB']['instance'].update({'laser': {'status': False}})

        # Set up the microwave generator
        self.instruments['mw_gen']['instance'].update({'enable_output': True})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulses']['mw_frequency']})
        self.instruments['mw_gen']['instance'].update({'freq_mode': 'CW'})
        self.instruments['mw_gen']['instance'].update({'power_mode': 'CW'})

        super(Rabi_RnS, self)._function()

        self.data['exp_finished'] = 1 # this is to flag when the experiment is done or stopped and this allows proper saving even fit fails

        # turn off laser, apd_switch and MW (just in case)
        self.instruments['PB']['instance'].update({'laser': {'status': False}})
        self.instruments['PB']['instance'].update({'apd_switch': {'status': False}})
        self.instruments['mw_gen']['instance'].update({'enable_output': False})

        # start fitting
        if 'counts' in self.data.keys() and 'tau' in self.data.keys() and not 0 in self.data['counts'][:, 0]:
            counts = self.data['counts'][:, 1] / self.data['counts'][:, 0]
            tau = self.data['tau']
            fit_success = False

            try:
                rabi_fits = fit_rabi_decay(tau, counts, variable_phase=True)
                self.data['fits'] = rabi_fits[0]
                fit_success = True
                # if rabi_fits[1] == True:
                #     self.data['fits'] = rabi_fits[0]
                #     fit_success = True
                # else:
                #     # self.data['fits'] = None
                #     fit_success = False
                #     print('rabi fit failed')
                #     self.log('rabi fit failed')
            except:
                # self.data['fits'] = None
                fit_success = False
                print('no rabi fit found')
                self.log('no rabi fit found')


            if fit_success:
                counts = self.data['counts'][:,1]/ self.data['counts'][:,0]
                tau =self.data['tau']
                fits = self.data['fits']

                RabiT = 2*np.pi / fits[1]
                # RabiT = (2 * np.pi / fits[1]) # avoid negative rabi period ZQ 8/25/2017
                T2_star = fits[4]
                phaseoffs = fits[2]
                pi_time = RabiT / 2 - phaseoffs * RabiT / (2 * np.pi)
                pi_half_time = RabiT / 4 - phaseoffs * RabiT / (2 * np.pi)
                three_pi_half_time = 3 * RabiT / 4 - phaseoffs * RabiT / (2 * np.pi)
                # pi_time = RabiT / 2
                # pi_half_time = RabiT / 4
                # three_pi_half_time = 3 * RabiT / 4
                self.data['phaseoffs'] = phaseoffs * RabiT / (2 * np.pi)
                self.data['pi_time'] = pi_time
                self.data['pi_half_time'] = pi_half_time
                self.data['three_pi_half_time'] = three_pi_half_time
                self.data['T2_star'] = T2_star
                self.data['rabi_freq'] = 1000 * fits[1] / (2 * np.pi) # Rabi frequency in [MHz]

    def _create_pulse_sequences(self):
        """

        Returns: pulse_sequences, num_averages, tau_list, meas_time
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
            meas_time: the width (in ns) of the daq measurement

        """
        MIN_DURATION = self.settings['min_pulse_dur']
        pulse_sequences = []
        # tau_list = range(int(max(15, self.settings['tau_times']['time_step'])), int(self.settings['tau_times']['max_time'] + 15),
        #                  self.settings['tau_times']['time_step'])
        # JG 16-08-25 changed (15ns min spacing is taken care of later):
       # tau_list = list(range(int(self.settings['tau_times']['min_time']),
                            #  int(self.settings['tau_times']['max_time']),
                            #  self.settings['tau_times']['time_step']))

        max_range = int(np.floor((self.settings['tau_times']['max_time']-self.settings['tau_times']['min_time'])/self.settings['tau_times']['time_step']))
        tau_list = np.array([self.settings['tau_times']['min_time'] + i*self.settings['tau_times']['time_step'] for i in range(max_range)])

        # ignore the sequence if the mw-pulse is shorter than MIN_DURATION (0 is ok because there is no mw pulse!)
        tau_list = [x for x in tau_list if x == 0 or x >= MIN_DURATION]
        print('tau_list in [ns]:', tau_list)

        meas_time = self.settings['read_out']['meas_time']
        nv_reset_time = self.settings['read_out']['nv_reset_time']
        laser_off_time = self.settings['read_out']['laser_off_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        delay_readout = self.settings['read_out']['delay_readout']
        initial_pulse_time = self.settings['mw_pulses']['pulse_time']

        # mw_sw_buffer = self.settings['read_out']['mw_switch_extra_time']
        # APD_sw_delay = self.settings['read_out']['APD_switch_delay']
        # microwave_channel = 'microwave_' + self.settings['mw_pulses']['microwave_channel']

        # print('start the for loop')

        for tau in tau_list:

            if self.settings['mw_pulses']['add_initial_pulse']:
                pulse_sequence = [Pulse('laser', laser_off_time + initial_pulse_time + delay_mw_readout + tau, nv_reset_time),
                                  Pulse('apd_readout', laser_off_time + initial_pulse_time + delay_mw_readout + tau + delay_readout, meas_time)]
                beginning_of_rabi = laser_off_time + initial_pulse_time + delay_mw_readout + tau + nv_reset_time
            else:
                pulse_sequence = [Pulse('laser', laser_off_time + tau, nv_reset_time),
                                  Pulse('apd_readout', laser_off_time + tau + delay_readout, meas_time)]
                beginning_of_rabi = laser_off_time + tau + nv_reset_time

            # if tau is 0 there is actually no mw pulse
            if tau > 0:
                if self.settings['mw_pulses']['add_initial_pulse']:
                    pulse_sequence.append(Pulse('microwave_switch_II', beginning_of_rabi + laser_off_time, initial_pulse_time))
                    pulse_sequence.append(
                        Pulse('microwave_switch_II', beginning_of_rabi + laser_off_time + initial_pulse_time + delay_mw_readout, tau))
                    end_of_rabi = beginning_of_rabi + laser_off_time + initial_pulse_time + delay_mw_readout + tau

                else:
                    pulse_sequence.append(Pulse('microwave_switch_II', beginning_of_rabi + laser_off_time, tau))
                    end_of_rabi = beginning_of_rabi + laser_off_time + initial_pulse_time + tau

                pulse_sequence.append(Pulse('laser', end_of_rabi + delay_mw_readout, nv_reset_time))
                pulse_sequence.append(Pulse('apd_readout',end_of_rabi  + delay_mw_readout+ delay_readout, meas_time))
                pulse_sequence.append(Pulse('apd_switch', 0, end_of_rabi  + delay_mw_readout+ nv_reset_time))

            # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
            # if tau == 0 or tau>=15:
            pulse_sequences.append(pulse_sequence)

        print('number of sequences before validation: ', len(pulse_sequences))
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

        if 'fits' in data.keys() and data['fits'] is not None:
            counts = data['counts'][:,1]/ data['counts'][:,0]
            mean_fluor = np.mean(data['counts'][:,0])
            tau = data['tau']
            fits = data['fits'] # amplitude, frequency, phase, offset

            # axislist[0].plot(tau, counts, 'b')
            axislist[0].plot(tau, counts, '.-')

            #axislist[0].hold(True) ER 20181012

            # axislist[0].plot(tau, cose_with_decay(tau, *fits), 'k', lw=3)
            axislist[0].plot(tau, cose_with_decay(tau, *fits), lw=3)


            # RabiT = 2 * np.pi / fits[1]
            T2_star = self.data['T2_star']
            phaseoffs = self.data['phaseoffs']
            pi_time = self.data['pi_time']
            pi_half_time = self.data['pi_half_time']
            three_pi_half_time = self.data['three_pi_half_time']
            # rabi_freq = 1000 * fits[1] / (2 * np.pi)
            rabi_freq = self.data['rabi_freq'] # Rabi frequency in [MHz]


            axislist[0].plot(pi_time, cose_with_decay(pi_time, *fits), 'ro', lw=3)
            axislist[0].annotate('$\pi$={:0.1f} ns'.format(pi_time), xy=(pi_time, cose_with_decay(pi_time, *fits)),
                                 xytext=(pi_time - 10., cose_with_decay(pi_time, *fits) - .01), xycoords='data')
            axislist[0].plot(pi_half_time, cose_with_decay(pi_half_time, *fits), 'ro', lw=3)
            axislist[0].annotate('$\pi/2$={:0.1f} ns'.format(pi_half_time),
                                 xy=(pi_half_time, cose_with_decay(pi_half_time, *fits)),
                                 xytext=(pi_half_time + 10., cose_with_decay(pi_half_time, *fits)), xycoords='data')
            axislist[0].plot(three_pi_half_time, cose_with_decay(three_pi_half_time, *fits), 'ro', lw=3)
            axislist[0].annotate('$3\pi/2$={:0.1f} ns'.format(three_pi_half_time),
                                 xy=(three_pi_half_time, cose_with_decay(three_pi_half_time, *fits)),
                                 xytext=(three_pi_half_time + 10., cose_with_decay(three_pi_half_time, *fits)),
                                 xycoords='data')
            axislist[0].plot(-phaseoffs, cose_with_decay(-phaseoffs, *fits), 'gd', lw=3)
            axislist[0].annotate('$start$={:0.1f} ns'.format(-phaseoffs),
                                 xy=(-phaseoffs, cose_with_decay(-phaseoffs, *fits)),
                                 xytext=(-phaseoffs + 10., cose_with_decay(-phaseoffs, *fits)), xycoords='data')

            #pi_time = 2*np.pi / fits[1] / 2
            # pi_time = (np.pi - fits[2])/fits[1]
            # pi_half_time = (np.pi/2 - fits[2])/fits[1]
            # three_pi_half_time = (3*np.pi/2 - fits[2])/fits[1]

         #   axislist[0].set_title('Rabi mw-power:{:0.1f}dBm, mw_freq:{:0.3f} GHz, pi-time: {:2.1f}ns, pi-half-time: {:2.1f}ns, 3pi_half_time: {:2.1f}ns, Rabi freq: {2.1f}MHz'.format(self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency']*1e-9, pi_time, pi_half_time, three_pi_half_time, rabi_freq))
            axislist[0].set_title(
                'R&S SMB100A \n Rabi: {:0.2f}MHz, pi-half time: {:2.1f}ns, pi-time: {:2.1f}ns, 3pi-half time: {:2.1f}ns \n T2*: {:2.1f}ns, Ref fluor: {:0.1f}kcps, mw-power: {:0.1f}dBm, mw_freq: {:0.4f} GHz'.format(
                    rabi_freq, pi_half_time, pi_time, three_pi_half_time, T2_star, mean_fluor,
                    self.settings['mw_pulses']['mw_power'],
                    self.settings['mw_pulses'][
                        'mw_frequency'] * 1e-9))
            axislist[0].set_xlabel('Rabi tau [ns]')
            axislist[0].set_ylabel('Normalized Fluorescence')

        elif 'exp_finished' in data.keys() is not None and data['exp_finished'] == 1:

            if 0 not in self.data['counts'][:, 0]:
                axislist[0].plot(data['tau'], data['counts'][:, 1] / data['counts'][:, 0])
                axislist[0].set_xlabel('Rabi tau [ns]')
                axislist[0].set_ylabel('Normalized Fluorescence')
            else:
                axislist[0].plot(data['tau'], data['counts'][:, 0])
                axislist[0].plot(data['tau'], data['counts'][:, 1])
                axislist[0].set_xlabel('Rabi tau [ns]')
                axislist[0].set_ylabel('Fluorescence [kcps]')
                axislist[0].legend(labels=('Ref Fluorescence', 'Rabi Data'), fontsize=8)

            axislist[0].set_title(
                '(final plot) R&S SMB100A \n Rabi mw-power:{:0.1f}dBm, mw_freq:{:0.4f} GHz'.format(self.settings['mw_pulses']['mw_power'],
                                                                       self.settings['mw_pulses'][
                                                                           'mw_frequency'] * 1e-9))
        else:

            super(Rabi_RnS, self)._plot(axislist)
            axislist[0].set_title('(initial plot) R&S SMB100A \n Rabi mw-power:{:0.1f}dBm, mw_freq:{:0.4f} GHz'.format(
                self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency'] * 1e-9))

            axislist[0].legend(labels=('Ref Fluorescence', 'Rabi Data'), fontsize=8)

    def _update_plot(self, axislist):
            # self._plot(axislist)

        if len(axislist[0].lines) == 0:
            self._plot(axislist)
            return
        super(Rabi_RnS, self)._update_plot(axislist)

        axislist[0].set_title('(updating plot) R&S SMB100A \n Rabi mw-power:{:0.1f}dBm, mw_freq:{:0.4f} GHz'.format(self.settings['mw_pulses']['mw_power'],
                                                                                         self.settings['mw_pulses'][
                                                                                             'mw_frequency'] * 1e-9))
        axislist[0].legend(labels=('Ref Fluorescence', 'Rabi Data'), fontsize=8)

class Rabi_N9310A(PulsedExperimentBaseScript): # ER 5.25.2017
    """
        This script applies a microwave pulse at fixed power for varying durations to measure Rabi oscillations.
        Uses a double_init scheme.

        The script uses Keysight (Agilent) signal generator N9310A.
        Note that the MW generator used in this script has I/Q modulation option.

        ==> Last edited by Ziwei 4/17/2019

    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_power', -10.0, float, 'microwave power in dBm'),
            Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
            Parameter('microwave_channel', 'i', ['i', 'q'], 'Channel to use for mw pulses'),
            Parameter('add_initial_pulse', False, bool, 'check to add an initial pulse in the beginning'),
            Parameter('pulse_time', 50.0, float, 'time duration of an initial pulse (in ns)')
        ]),
        Parameter('tau_times', [
            Parameter('min_time', 15, float, 'minimum time for rabi oscillations (in ns)'),
            Parameter('max_time', 1000, float, 'total time of rabi oscillations (in ns)'),
            Parameter('time_step', 20., [2.5, 5., 10., 15., 20., 50., 100., 200., 500., 1000., 10000., 100000., 500000.],
                  'time step increment of rabi pulse duration (in ns)')
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 200, float, 'measurement time after rabi sequence (in ns)'),
            Parameter('nv_reset_time', 2000, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 600, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 550, int, '[ns] delay between laser on and APD readout (given by spontaneous decay rate)'),
        ]),
        Parameter('min_pulse_dur', 15, [10, 15, 50], 'Minimum allowed pulse duration (ns)'),
        Parameter('num_averages', 500000, int, 'number of averages (>100000)'),
    ]

    _INSTRUMENTS = {'daq': NI6353, 'PB': LISE607RTPulseBlaster, 'mw_gen_a': AgilentMicrowaveGenerator}

    def _function(self):


        ## Turn off green light (the pulse blaster will pulse it on when needed)
        self.instruments['PB']['instance'].update({'laser': {'status': False}})
        self.instruments['PB']['instance'].update({'apd_switch': {'status': False}})

        self.instruments['mw_gen_a']['instance'].update({'enable_IQ': True})
        self.instruments['mw_gen_a']['instance'].update({'enable_modulation': True})
        self.instruments['mw_gen_a']['instance'].update({'freq_mode': 'CW'})
        self.instruments['mw_gen_a']['instance'].update({'power_mode': 'CW'})
        self.instruments['mw_gen_a']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power']})
        self.instruments['mw_gen_a']['instance'].update({'frequency': self.settings['mw_pulses']['mw_frequency']})
        self.instruments['mw_gen_a']['instance'].update({'enable_output': True})

        self.avg_block_number = super(Rabi_N9310A, self)._function()

        self.data['exp_finished'] = 1 # this is to flag when the experiment is done or stopped and this allows proper saving even fit fails
        self.data['avg_block_number'] = self.avg_block_number

        # turn off laser, apd_switch and MW (just in case)
        self.instruments['PB']['instance'].update({'laser': {'status': False}})
        self.instruments['PB']['instance'].update({'apd_switch': {'status': False}})
        self.instruments['mw_gen_a']['instance'].update({'enable_output': False})
        self.instruments['mw_gen_a']['instance'].update({'enable_IQ': False})
        self.instruments['mw_gen_a']['instance'].update({'enable_modulation': False})

        # start fitting
        if 'counts' in self.data.keys() and 'tau' in self.data.keys() and not 0 in self.data['counts'][:, 0]:
            counts = self.data['counts'][:, 1] / self.data['counts'][:, 0]
            tau = self.data['tau']
            fit_success = False

            try:
                rabi_fits = fit_rabi_decay(tau, counts, variable_phase=True)
                self.data['fits'] = rabi_fits[0]
                fit_success = True
                # if rabi_fits[1] == True:
                #     self.data['fits'] = rabi_fits[0]
                #     fit_success = True
                # else:
                #     # self.data['fits'] = None
                #     fit_success = False
                #     print('rabi fit failed')
                #     self.log('rabi fit failed')
            except:
                # self.data['fits'] = None
                fit_success = False
                print('no rabi fit found')
                self.log('no rabi fit found')


            if fit_success:
                counts = self.data['counts'][:,1]/ self.data['counts'][:,0]
                tau =self.data['tau']
                fits = self.data['fits']

                RabiT = 2*np.pi / fits[1]
                T2_star = fits[4]
                phaseoffs = fits[2]
                pi_time = RabiT / 2 - phaseoffs * RabiT / (2 * np.pi)
                pi_half_time = RabiT / 4 - phaseoffs * RabiT / (2 * np.pi)
                three_pi_half_time = 3 * RabiT / 4 - phaseoffs * RabiT / (2 * np.pi)
                self.data['phaseoffs'] = phaseoffs * RabiT / (2 * np.pi)
                self.data['pi_time'] = pi_time
                self.data['pi_half_time'] = pi_half_time
                self.data['three_pi_half_time'] = three_pi_half_time
                self.data['T2_star'] = T2_star
                self.data['rabi_freq'] = 1000 * fits[1] / (2 * np.pi) # Rabi frequency in [MHz]

    def _create_pulse_sequences(self):
        """

        Returns: pulse_sequences, num_averages, tau_list, meas_time
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
            meas_time: the width (in ns) of the daq measurement

        """

        MIN_DURATION = self.settings['min_pulse_dur']
        pulse_sequences = []
        max_range = int(np.floor((self.settings['tau_times']['max_time']-self.settings['tau_times']['min_time'])/self.settings['tau_times']['time_step']))
        tau_list = np.array([self.settings['tau_times']['min_time'] + i*self.settings['tau_times']['time_step'] for i in range(max_range)])

        # ignore the sequence if the mw-pulse is shorter than MIN_DURATION (0 is ok because there is no mw pulse!)
        tau_list = [x for x in tau_list if x == 0 or x >= MIN_DURATION]
        print('tau_list in [ns]:', tau_list)

        meas_time = self.settings['read_out']['meas_time']
        nv_reset_time = self.settings['read_out']['nv_reset_time']
        laser_off_time = self.settings['read_out']['laser_off_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        delay_readout = self.settings['read_out']['delay_readout']
        initial_pulse_time = self.settings['mw_pulses']['pulse_time']
        microwave_channel = 'microwave_' + self.settings['mw_pulses']['microwave_channel']

        for tau in tau_list:

            if self.settings['mw_pulses']['add_initial_pulse']:
                pulse_sequence = [Pulse('laser', laser_off_time + initial_pulse_time + delay_mw_readout + tau, nv_reset_time),
                                  Pulse('apd_readout', laser_off_time + initial_pulse_time + delay_mw_readout + tau + delay_readout, meas_time)]
                beginning_of_rabi = laser_off_time + initial_pulse_time + delay_mw_readout + tau + nv_reset_time
            else:
                pulse_sequence = [Pulse('laser', laser_off_time + tau, nv_reset_time),
                                  Pulse('apd_readout', laser_off_time + tau + delay_readout, meas_time)]
                beginning_of_rabi = laser_off_time + tau + nv_reset_time

            # if tau is 0 there is actually no mw pulse
            if tau > 0:
                if self.settings['mw_pulses']['add_initial_pulse']:
                    pulse_sequence.append(Pulse(microwave_channel, beginning_of_rabi + laser_off_time, initial_pulse_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel, beginning_of_rabi + laser_off_time + initial_pulse_time + delay_mw_readout, tau))
                    end_of_rabi = beginning_of_rabi + laser_off_time + initial_pulse_time + delay_mw_readout + tau

                else:
                    pulse_sequence.append(Pulse(microwave_channel, beginning_of_rabi + laser_off_time, tau))
                    end_of_rabi = beginning_of_rabi + laser_off_time + initial_pulse_time + tau

                pulse_sequence.append(Pulse('laser', end_of_rabi + delay_mw_readout, nv_reset_time))
                pulse_sequence.append(Pulse('apd_readout',end_of_rabi  + delay_mw_readout+ delay_readout, meas_time))
                pulse_sequence.append(Pulse('apd_switch', 0, end_of_rabi  + delay_mw_readout+ nv_reset_time))

            # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
            # if tau == 0 or tau>=15:
            pulse_sequences.append(pulse_sequence)

        print('number of sequences before validation: ', len(pulse_sequences))
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

        if 'avg_block_number' in data.keys() is not None and data['avg_block_number'] is not None:
            avg_block_number = data['avg_block_number']
        else:
            avg_block_number = 0

        if 'fits' in data.keys() and data['fits'] is not None:
            counts = data['counts'][:,1]/ data['counts'][:,0]
            mean_fluor = np.mean(data['counts'][:,0])
            tau = data['tau']
            fits = data['fits'] # amplitude, frequency, phase, offset

            axislist[0].plot(tau, counts, '.-')
            tauinterp = np.linspace(np.min(tau), np.max(tau), 100)
            # axislist[0].plot(tau, cose_with_decay(tau, *fits), lw=3)
            axislist[0].plot(tauinterp, cose_with_decay(tauinterp, *fits), lw=3)
            T2_star = self.data['T2_star']
            phaseoffs = self.data['phaseoffs']
            pi_time = self.data['pi_time']
            pi_half_time = self.data['pi_half_time']
            three_pi_half_time = self.data['three_pi_half_time']
            rabi_freq = self.data['rabi_freq'] # Rabi frequency in [MHz]


            axislist[0].plot(pi_time, cose_with_decay(pi_time, *fits), 'ro', lw=3)
            axislist[0].annotate('$\pi$={:0.1f} ns'.format(pi_time), xy=(pi_time, cose_with_decay(pi_time, *fits)),
                                 xytext=(pi_time - 10., cose_with_decay(pi_time, *fits) - .01), xycoords='data')
            axislist[0].plot(pi_half_time, cose_with_decay(pi_half_time, *fits), 'ro', lw=3)
            axislist[0].annotate('$\pi/2$={:0.1f} ns'.format(pi_half_time),
                                 xy=(pi_half_time, cose_with_decay(pi_half_time, *fits)),
                                 xytext=(pi_half_time + 10., cose_with_decay(pi_half_time, *fits)), xycoords='data')
            axislist[0].plot(three_pi_half_time, cose_with_decay(three_pi_half_time, *fits), 'ro', lw=3)
            axislist[0].annotate('$3\pi/2$={:0.1f} ns'.format(three_pi_half_time),
                                 xy=(three_pi_half_time, cose_with_decay(three_pi_half_time, *fits)),
                                 xytext=(three_pi_half_time + 10., cose_with_decay(three_pi_half_time, *fits)),
                                 xycoords='data')
            axislist[0].plot(-phaseoffs, cose_with_decay(-phaseoffs, *fits), 'gd', lw=3)
            axislist[0].annotate('$start$={:0.1f} ns'.format(-phaseoffs),
                                 xy=(-phaseoffs, cose_with_decay(-phaseoffs, *fits)),
                                 xytext=(-phaseoffs + 10., cose_with_decay(-phaseoffs, *fits)), xycoords='data')


            axislist[0].set_title(
                'Keysight N9310A, pulse carved by ' + self.settings['mw_pulses'][
                    'microwave_channel'] + '\nRabi: {:0.2f}MHz, pi-half time: {:2.1f}ns, pi-time: {:2.1f}ns, 3pi-half time: {:2.1f}ns \n after {:d} averages \n T2*: {:2.1f}ns, Ref fluor: {:0.1f}kcps \n mw-power: {:0.2f}dBm, mw_freq: {:0.4f} GHz'.format(
                    rabi_freq, pi_half_time, pi_time, three_pi_half_time, avg_block_number, T2_star, mean_fluor,
                    self.settings['mw_pulses']['mw_power'],
                    self.settings['mw_pulses'][
                        'mw_frequency'] * 1e-9))
            axislist[0].set_xlabel('Rabi tau [ns]')
            axislist[0].set_ylabel('Normalized Fluorescence')

        elif 'exp_finished' in data.keys() is not None and data['exp_finished'] == 1:

            if 0 not in self.data['counts'][:, 0]:
                axislist[0].plot(data['tau'], data['counts'][:, 1] / data['counts'][:, 0])
                axislist[0].set_xlabel('Rabi tau [ns]')
                axislist[0].set_ylabel('Normalized Fluorescence')
            else:
                axislist[0].plot(data['tau'], data['counts'][:, 0])
                axislist[0].plot(data['tau'], data['counts'][:, 1])
                axislist[0].set_xlabel('Rabi tau [ns]')
                axislist[0].set_ylabel('Fluorescence [kcps]')
                axislist[0].legend(labels=('Ref Fluorescence', 'Rabi Data'), fontsize=8)

            axislist[0].set_title(
                '(final plot) Keysight N9310A, pulse carved by ' + self.settings['mw_pulses'][
                    'microwave_channel'] + '\n after {:d} averages \nRabi: mw-power:{:0.2f}dBm, mw_freq:{:0.4f} GHz'.format(
                    avg_block_number, self.settings['mw_pulses']['mw_power'],
                    self.settings['mw_pulses'][
                        'mw_frequency'] * 1e-9))
        else:

            super(Rabi_N9310A, self)._plot(axislist)
            axislist[0].set_title('(initial plot) Keysight N9310A, pulse carved by ' + self.settings['mw_pulses'][
                'microwave_channel'] + '\n after {:d} averages \nRabi: mw-power:{:0.2f}dBm, mw_freq:{:0.4f} GHz'.format(
                avg_block_number, self.settings['mw_pulses']['mw_power'],
                self.settings['mw_pulses']['mw_frequency'] * 1e-9))

            axislist[0].legend(labels=('Ref Fluorescence', 'Rabi Data'), fontsize=8)

    def _update_plot(self, axislist):
            # self._plot(axislist)

        if len(axislist[0].lines) == 0:
            self._plot(axislist)
            return
        super(Rabi_N9310A, self)._update_plot(axislist)

        axislist[0].set_title('(updating plot) Keysight N9310A, pulse carved by ' + self.settings['mw_pulses']['microwave_channel'] + '\nRabi: mw-power:{:0.2f}dBm, mw_freq:{:0.4f} GHz'.format(self.settings['mw_pulses']['mw_power'],
                                                                                         self.settings['mw_pulses'][
                                                                                             'mw_frequency'] * 1e-9))
        axislist[0].legend(labels=('Ref Fluorescence', 'Rabi Data'), fontsize=8)

class Rabi_RnSIQ(PulsedExperimentBaseScript): # ER 5.25.2017
    """
        This script applies a microwave pulse at fixed power for varying durations to measure Rabi oscillations.
        The script uses Rohde&Schwarz signal generator SMB100A.
        The phase is controlled by an IQ mixer.

        ==> Last edited by Ziwei 5/11/2019

    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_power', -10.0, float, 'microwave power in dBm'),
            Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
            Parameter('microwave_channel', 'i', ['i', 'q'], 'Channel to use for mw pulses'),
            Parameter('add_initial_pulse', False, bool, 'check to add an initial pulse in the beginning'),
            Parameter('pulse_time', 50.0, float, 'time duration of an initial pulse (in ns)')
        ]),
        Parameter('tau_times', [
            Parameter('min_time', 15, float, 'minimum time for rabi oscillations (in ns)'),
            Parameter('max_time', 1000, float, 'total time of rabi oscillations (in ns)'),
            Parameter('time_step', 20., [2.5, 5., 10., 20., 50., 100., 200., 500., 1000., 10000., 100000., 500000.],
                  'time step increment of rabi pulse duration (in ns)')
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 235, float, 'measurement time after rabi sequence (in ns)'),
            Parameter('nv_reset_time', 2000, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 600, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 550, int, '[ns] delay between laser on and APD readout (given by spontaneous decay rate)'),
        ]),
        Parameter('min_pulse_dur', 15, [10, 15, 50], 'Minimum allowed pulse duration (ns)'),
        Parameter('num_averages', 500000, int, 'number of averages (>100000)'),
    ]

    _INSTRUMENTS = {'daq': NI6353, 'PB': LISE607RTPulseBlaster, 'mw_gen': R8SMicrowaveGenerator}

    def _function(self):


        ## Turn off green light (the pulse blaster will pulse it on when needed)
        self.instruments['PB']['instance'].update({'laser': {'status': False}})
        self.instruments['PB']['instance'].update({'apd_switch': {'status': False}})

        # Set up the microwave generator SMB100A
        self.instruments['mw_gen']['instance'].update({'freq_mode': 'CW'})
        self.instruments['mw_gen']['instance'].update({'power_mode': 'CW'})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulses']['mw_frequency']})
        self.instruments['mw_gen']['instance'].update({'enable_output': True})
        # This lets the RF go through the IQ mixer
        self.instruments['PB']['instance'].update({'microwave_switch_II': {'status': False}})

        super(Rabi_RnSIQ, self)._function()

        self.data['exp_finished'] = 1 # this is to flag when the experiment is done or stopped and this allows proper saving even fit fails

        # turn off laser, apd_switch and MW (just in case)
        self.instruments['mw_gen']['instance'].update({'enable_output': False})
        self.instruments['PB']['instance'].update({'laser': {'status': False}})
        self.instruments['PB']['instance'].update({'apd_switch': {'status': False}})

        # start fitting
        if 'counts' in self.data.keys() and 'tau' in self.data.keys() and not 0 in self.data['counts'][:, 0]:
            counts = self.data['counts'][:, 1] / self.data['counts'][:, 0]
            tau = self.data['tau']
            fit_success = False

            try:
                rabi_fits = fit_rabi_decay(tau, counts, variable_phase=True)
                self.data['fits'] = rabi_fits[0]
                fit_success = True
                # if rabi_fits[1] == True:
                #     self.data['fits'] = rabi_fits[0]
                #     fit_success = True
                # else:
                #     # self.data['fits'] = None
                #     fit_success = False
                #     print('rabi fit failed')
                #     self.log('rabi fit failed')
            except:
                # self.data['fits'] = None
                fit_success = False
                print('no rabi fit found')
                self.log('no rabi fit found')


            if fit_success:
                counts = self.data['counts'][:,1]/ self.data['counts'][:,0]
                tau =self.data['tau']
                fits = self.data['fits']

                RabiT = 2*np.pi / fits[1]
                T2_star = fits[4]
                phaseoffs = fits[2]
                pi_time = RabiT / 2 - phaseoffs * RabiT / (2 * np.pi)
                pi_half_time = RabiT / 4 - phaseoffs * RabiT / (2 * np.pi)
                three_pi_half_time = 3 * RabiT / 4 - phaseoffs * RabiT / (2 * np.pi)
                self.data['phaseoffs'] = phaseoffs * RabiT / (2 * np.pi)
                self.data['pi_time'] = pi_time
                self.data['pi_half_time'] = pi_half_time
                self.data['three_pi_half_time'] = three_pi_half_time
                self.data['T2_star'] = T2_star
                self.data['rabi_freq'] = 1000 * fits[1] / (2 * np.pi) # Rabi frequency in [MHz]

    def _create_pulse_sequences(self):
        """

        Returns: pulse_sequences, num_averages, tau_list, meas_time
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
            meas_time: the width (in ns) of the daq measurement

        """

        MIN_DURATION = self.settings['min_pulse_dur']
        pulse_sequences = []
        max_range = int(np.floor((self.settings['tau_times']['max_time']-self.settings['tau_times']['min_time'])/self.settings['tau_times']['time_step']))
        tau_list = np.array([self.settings['tau_times']['min_time'] + i*self.settings['tau_times']['time_step'] for i in range(max_range)])

        # ignore the sequence if the mw-pulse is shorter than MIN_DURATION (0 is ok because there is no mw pulse!)
        tau_list = [x for x in tau_list if x == 0 or x >= MIN_DURATION]
        print('tau_list in [ns]:', tau_list)

        meas_time = self.settings['read_out']['meas_time']
        nv_reset_time = self.settings['read_out']['nv_reset_time']
        laser_off_time = self.settings['read_out']['laser_off_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        delay_readout = self.settings['read_out']['delay_readout']
        initial_pulse_time = self.settings['mw_pulses']['pulse_time']
        microwave_channel = 'microwave_' + self.settings['mw_pulses']['microwave_channel']+'_II'

        for tau in tau_list:

            if self.settings['mw_pulses']['add_initial_pulse']:
                pulse_sequence = [Pulse('laser', laser_off_time + initial_pulse_time + delay_mw_readout + tau, nv_reset_time),
                                  Pulse('apd_readout', laser_off_time + initial_pulse_time + delay_mw_readout + tau + delay_readout, meas_time)]
                beginning_of_rabi = laser_off_time + initial_pulse_time + delay_mw_readout + tau + nv_reset_time
            else:
                pulse_sequence = [Pulse('laser', laser_off_time + tau, nv_reset_time),
                                  Pulse('apd_readout', laser_off_time + tau + delay_readout, meas_time)]
                beginning_of_rabi = laser_off_time + tau + nv_reset_time

            # if tau is 0 there is actually no mw pulse
            if tau > 0:
                if self.settings['mw_pulses']['add_initial_pulse']:
                    pulse_sequence.append(Pulse(microwave_channel, beginning_of_rabi + laser_off_time, initial_pulse_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel, beginning_of_rabi + laser_off_time + initial_pulse_time + delay_mw_readout, tau))
                    end_of_rabi = beginning_of_rabi + laser_off_time + initial_pulse_time + delay_mw_readout + tau

                else:
                    pulse_sequence.append(Pulse(microwave_channel, beginning_of_rabi + laser_off_time, tau))
                    end_of_rabi = beginning_of_rabi + laser_off_time + initial_pulse_time + tau

                pulse_sequence.append(Pulse('laser', end_of_rabi + delay_mw_readout, nv_reset_time))
                pulse_sequence.append(Pulse('apd_readout',end_of_rabi  + delay_mw_readout+ delay_readout, meas_time))
                pulse_sequence.append(Pulse('apd_switch', 0, end_of_rabi  + delay_mw_readout+ nv_reset_time))

            # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
            # if tau == 0 or tau>=15:
            pulse_sequences.append(pulse_sequence)

        print('number of sequences before validation: ', len(pulse_sequences))
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

        if 'fits' in data.keys() and data['fits'] is not None:
            counts = data['counts'][:,1]/ data['counts'][:,0]
            mean_fluor = np.mean(data['counts'][:,0])
            tau = data['tau']
            fits = data['fits'] # amplitude, frequency, phase, offset

            axislist[0].plot(tau, counts, '.-')
            axislist[0].plot(tau, cose_with_decay(tau, *fits), lw=3)
            T2_star = self.data['T2_star']
            phaseoffs = self.data['phaseoffs']
            pi_time = self.data['pi_time']
            pi_half_time = self.data['pi_half_time']
            three_pi_half_time = self.data['three_pi_half_time']
            rabi_freq = self.data['rabi_freq'] # Rabi frequency in [MHz]


            axislist[0].plot(pi_time, cose_with_decay(pi_time, *fits), 'ro', lw=3)
            axislist[0].annotate('$\pi$={:0.1f} ns'.format(pi_time), xy=(pi_time, cose_with_decay(pi_time, *fits)),
                                 xytext=(pi_time - 10., cose_with_decay(pi_time, *fits) - .01), xycoords='data')
            axislist[0].plot(pi_half_time, cose_with_decay(pi_half_time, *fits), 'ro', lw=3)
            axislist[0].annotate('$\pi/2$={:0.1f} ns'.format(pi_half_time),
                                 xy=(pi_half_time, cose_with_decay(pi_half_time, *fits)),
                                 xytext=(pi_half_time + 10., cose_with_decay(pi_half_time, *fits)), xycoords='data')
            axislist[0].plot(three_pi_half_time, cose_with_decay(three_pi_half_time, *fits), 'ro', lw=3)
            axislist[0].annotate('$3\pi/2$={:0.1f} ns'.format(three_pi_half_time),
                                 xy=(three_pi_half_time, cose_with_decay(three_pi_half_time, *fits)),
                                 xytext=(three_pi_half_time + 10., cose_with_decay(three_pi_half_time, *fits)),
                                 xycoords='data')
            axislist[0].plot(-phaseoffs, cose_with_decay(-phaseoffs, *fits), 'gd', lw=3)
            axislist[0].annotate('$start$={:0.1f} ns'.format(-phaseoffs),
                                 xy=(-phaseoffs, cose_with_decay(-phaseoffs, *fits)),
                                 xytext=(-phaseoffs + 10., cose_with_decay(-phaseoffs, *fits)), xycoords='data')
            axislist[0].set_title(
                'R&S SMB100A, pulse carved by ' + self.settings['mw_pulses'][
                    'microwave_channel'] + '\nRabi: {:0.2f}MHz, pi-half time: {:2.1f}ns, pi-time: {:2.1f}ns, 3pi-half time: {:2.1f}ns \n T2*: {:2.1f}ns, Ref fluor: {:0.1f}kcps, mw-power: {:0.2f}dBm, mw_freq: {:0.4f} GHz'.format(
                    rabi_freq, pi_half_time, pi_time, three_pi_half_time, T2_star, mean_fluor,
                    self.settings['mw_pulses']['mw_power'],
                    self.settings['mw_pulses'][
                        'mw_frequency'] * 1e-9))
            axislist[0].set_xlabel('Rabi tau [ns]')
            axislist[0].set_ylabel('Normalized Fluorescence')

        elif 'exp_finished' in data.keys() is not None and data['exp_finished'] == 1:

            if 0 not in self.data['counts'][:, 0]:
                axislist[0].plot(data['tau'], data['counts'][:, 1] / data['counts'][:, 0])
                axislist[0].set_xlabel('Rabi tau [ns]')
                axislist[0].set_ylabel('Normalized Fluorescence')
            else:
                axislist[0].plot(data['tau'], data['counts'][:, 0])
                axislist[0].plot(data['tau'], data['counts'][:, 1])
                axislist[0].set_xlabel('Rabi tau [ns]')
                axislist[0].set_ylabel('Fluorescence [kcps]')
                axislist[0].legend(labels=('Ref Fluorescence', 'Rabi Data'), fontsize=8)

            axislist[0].set_title(
                '(final plot) R&S SMB100A, pulse carved by ' + self.settings['mw_pulses']['microwave_channel'] + '\nRabi: mw-power:{:0.2f}dBm, mw_freq:{:0.4f} GHz'.format(self.settings['mw_pulses']['mw_power'],
                                                                       self.settings['mw_pulses'][
                                                                           'mw_frequency'] * 1e-9))
        else:

            super(Rabi_RnSIQ, self)._plot(axislist)
            axislist[0].set_title('(initial plot) R&S SMB100A, pulse carved by ' + self.settings['mw_pulses']['microwave_channel'] + '\nRabi: mw-power:{:0.2f}dBm, mw_freq:{:0.4f} GHz'.format(
                self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency'] * 1e-9))

            axislist[0].legend(labels=('Ref Fluorescence', 'Rabi Data'), fontsize=8)

    def _update_plot(self, axislist):
            # self._plot(axislist)

        if len(axislist[0].lines) == 0:
            self._plot(axislist)
            return
        super(Rabi_RnSIQ, self)._update_plot(axislist)

        axislist[0].set_title('(updating plot) R&S SMB100A, pulse carved by ' + self.settings['mw_pulses']['microwave_channel'] + '\nRabi: mw-power:{:0.2f}dBm, mw_freq:{:0.4f} GHz'.format(self.settings['mw_pulses']['mw_power'],
                                                                                         self.settings['mw_pulses'][
                                                                                             'mw_frequency'] * 1e-9))
        axislist[0].legend(labels=('Ref Fluorescence', 'Rabi Data'), fontsize=8)

class Ramsey_RnS(PulsedExperimentBaseScript):
    """
        This script measures the FID signal of NV with Ramsey sequence to measure its dephasing time.
        Rhode&Schwarz signal generator SMB100A is used.
        Note that the MW generator used in this script doesn't have I/Q modulation option.

        ==> Last edited by Ziwei 4/17/2019

    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_power', -45.0, float, 'microwave power in dBm'),
            Parameter('resonant_freq', 2.87e9, float, 'resonant frequency in Hz'),
            Parameter('detuning', -50e6, float, 'detuning for Ramsey experiment in Hz]'),
            Parameter('pi_half_pulse_time', 50.0, float, 'time duration of a pi/2 pulse (in ns)'),
            Parameter('three_pi_half_pulse_time', 135.0, float, 'time duration of a 3pi/2 pulse (in ns)'),
        ]),
        Parameter('tau_times', [
            Parameter('min_time', 20, float, 'minimum time between the two pi-half pulses (in ns)'),
            Parameter('max_time', 2000, float, 'maximum time between the two pi-half pulses (in ns)'),
            Parameter('time_step', 20., [2.5, 5., 10., 20., 50., 100., 200., 500., 1000., 10000., 100000., 500000.],
                  'time step increment of rabi pulse duration (in ns)')
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 235, float, 'measurement time after rabi sequence (in ns)'),
            Parameter('nv_reset_time', 2000, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 600, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 550, int, '[ns] delay between laser on and APD readout (given by spontaneous decay rate)')
        ]),
        Parameter('min_pulse_dur', 15, [10, 15, 50], 'Minimum allowed pulse duration (ns) > 10ns'),
        Parameter('num_averages', 500000, int, 'number of averages (>100000)'),
    ]

    _INSTRUMENTS = {'daq': NI6353, 'PB': LISE607RTPulseBlaster, 'mw_gen': R8SMicrowaveGenerator}

    def _function(self):

        freq_to_do = self.settings['mw_pulses']['resonant_freq'] + self.settings['mw_pulses']['detuning']
        ## Turn off green light (the pulse blaster will pulse it on when needed)
        self.instruments['PB']['instance'].update({'laser': {'status': False}})
        self.instruments['PB']['instance'].update({'apd_switch': {'status': False}})
        self.instruments['mw_gen']['instance'].update({'enable_output': True})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power']})
        self.instruments['mw_gen']['instance'].update({'frequency': freq_to_do})


        super(Ramsey_RnS, self)._function()

        self.data['exp_finished'] = 1 # this is to flag when the experiment is done or stopped. this allows proper saving even fit fails

        # turn off laser, apd_switch and MW (just in case)
        self.instruments['PB']['instance'].update({'laser': {'status': False}})
        self.instruments['PB']['instance'].update({'apd_switch': {'status': False}})
        self.instruments['mw_gen']['instance'].update({'enable_output': False})

        # start fitting
        # if 'counts' in self.data.keys() and 'tau' in self.data.keys() and not 0 in self.data['counts'][:, 0]:
        #     counts = self.data['counts'][:, 1] / self.data['counts'][:, 0]
        #     tau = self.data['tau']
        #     fit_success = False
            # try:
            #     rabi_fits = fit_rabi_decay(tau, counts, variable_phase=True)
            #     if rabi_fits[1] == True:
            #         self.data['fits'] = rabi_fits[0]
            #         fit_success = True
            #     else:
            #         # self.data['fits'] = None
            #         fit_success = False
            #         print('rabi fit failed')
            #         self.log('rabi fit failed')
            # except:
            #     # self.data['fits'] = None
            #     fit_success = False
            #     print('rabi fit failed')
            #     self.log('rabi fit failed')
            #
            # if fit_success:
            #     counts = self.data['counts'][:,1]/ self.data['counts'][:,0]
            #     tau =self.data['tau']
            #     fits = self.data['fits']
            #
            #     RabiT = 2*np.pi / fits[1]
            #     # RabiT = (2 * np.pi / fits[1]) # avoid negative rabi period ZQ 8/25/2017
            #     T2_star = fits[4]
            #     phaseoffs = fits[2]
            #     pi_time = RabiT / 2 - phaseoffs * RabiT / (2 * np.pi)
            #     pi_half_time = RabiT / 4 - phaseoffs * RabiT / (2 * np.pi)
            #     three_pi_half_time = 3 * RabiT / 4 - phaseoffs * RabiT / (2 * np.pi)
            #     # pi_time = RabiT / 2
            #     # pi_half_time = RabiT / 4
            #     # three_pi_half_time = 3 * RabiT / 4
            #     self.data['phaseoffs'] = phaseoffs * RabiT / (2 * np.pi)
            #     self.data['pi_time'] = pi_time
            #     self.data['pi_half_time'] = pi_half_time
            #     self.data['three_pi_half_time'] = three_pi_half_time
            #     self.data['T2_star'] = T2_star
            #     self.data['rabi_freq'] = 1000 * fits[1] / (2 * np.pi) # Rabi frequency in [MHz]

    def _create_pulse_sequences(self):
        """

        Returns: pulse_sequences, num_averages, tau_list, meas_time
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
            meas_time: the width (in ns) of the daq measurement

        """
        MIN_DURATION = self.settings['tau_times']['min_pulse_dur']
        pulse_sequences = []
        max_range = int(np.floor((self.settings['tau_times']['max_time']-self.settings['tau_times']['min_time'])/self.settings['tau_times']['time_step']))
        tau_list = np.array([self.settings['tau_times']['min_time'] + i*self.settings['tau_times']['time_step'] for i in range(max_range)])

        # ignore the sequence if the mw-pulse is shorter than MIN_DURATION (0 is ok because there is no mw pulse!)
        tau_list = [x for x in tau_list if x == 0 or x >= MIN_DURATION]
        print('tau_list in [ns]:', tau_list)

        pi_half_pulse_time = self.settings['mw_pulses']['pi_half_pulse_time']
        three_pi_half_pulse_time = self.settings['mw_pulses']['three_pi_half_pulse_time']
        meas_time = self.settings['read_out']['meas_time']
        nv_reset_time = self.settings['read_out']['nv_reset_time']
        laser_off_time = self.settings['read_out']['laser_off_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        delay_readout = self.settings['read_out']['delay_readout']


        for tau in tau_list:
            # Normal Ramsey sequence
            pulse_sequence = [Pulse('laser', laser_off_time + pi_half_pulse_time + tau + pi_half_pulse_time, nv_reset_time),
                              Pulse('apd_readout', laser_off_time + pi_half_pulse_time + tau + pi_half_pulse_time + delay_readout, meas_time),
                              ]

            end_of_first_laser = laser_off_time + pi_half_pulse_time + tau + pi_half_pulse_time + nv_reset_time

            pulse_sequence.append(Pulse('microwave_switch_II', end_of_first_laser + laser_off_time, pi_half_pulse_time))
            pulse_sequence.append(Pulse('microwave_switch_II',
                                        end_of_first_laser + laser_off_time + pi_half_pulse_time + tau, pi_half_pulse_time))
            pulse_sequence.append(Pulse('laser',
                                        end_of_first_laser + laser_off_time + pi_half_pulse_time + tau + pi_half_pulse_time + delay_mw_readout,
                                        nv_reset_time))
            pulse_sequence.append(Pulse('apd_readout',
                                        end_of_first_laser + laser_off_time + pi_half_pulse_time + tau + pi_half_pulse_time + delay_mw_readout + delay_readout,
                                        meas_time))

            pulse_sequence.append(
                Pulse('apd_switch', 0, end_of_first_laser + laser_off_time + pi_half_pulse_time + tau + pi_half_pulse_time + delay_mw_readout + nv_reset_time))

            # The sequence implementing 3pi/2 readout...
            # start_of_first_Ramsey = laser_off_time
            # # the first Ramsey with pi/2 readout
            # pulse_sequence = [Pulse('microwave_switch', start_of_first_Ramsey, pi_half_pulse_time)]
            # pulse_sequence.append(Pulse('microwave_switch', start_of_first_Ramsey + pi_half_pulse_time + tau, pi_half_pulse_time))
            # pulse_sequence.append(Pulse('laser',
            #                             start_of_first_Ramsey + pi_half_pulse_time + tau + pi_half_pulse_time + delay_mw_readout,
            #                             nv_reset_time))
            # pulse_sequence.append(Pulse('apd_readout',
            #                             start_of_first_Ramsey + pi_half_pulse_time + tau + pi_half_pulse_time + delay_mw_readout + delay_readout,
            #                             meas_time))
            #
            # start_of_second_Ramsey = start_of_first_Ramsey + pi_half_pulse_time + tau + pi_half_pulse_time + delay_mw_readout + nv_reset_time + laser_off_time
            # # the second Ramsey with 3pi/2 readout
            # pulse_sequence.append(Pulse('microwave_switch', start_of_second_Ramsey, pi_half_pulse_time))
            # pulse_sequence.append(
            #     Pulse('microwave_switch', start_of_second_Ramsey + pi_half_pulse_time + tau, three_pi_half_pulse_time))
            # pulse_sequence.append(Pulse('laser',
            #                             start_of_second_Ramsey + pi_half_pulse_time + tau + three_pi_half_pulse_time + delay_mw_readout,
            #                             nv_reset_time))
            # pulse_sequence.append(Pulse('apd_readout',
            #                             start_of_second_Ramsey + pi_half_pulse_time + tau + three_pi_half_pulse_time + delay_mw_readout + delay_readout,
            #                             meas_time))
            # pulse_sequence.append(
            #     Pulse('apd_switch', 0, start_of_second_Ramsey + pi_half_pulse_time + tau + three_pi_half_pulse_time + delay_mw_readout + nv_reset_time))


            # pulse_sequence = [
            #     Pulse('laser', laser_off_time + pi_half_pulse_time + tau + pi_half_pulse_time, nv_reset_time),
            #     Pulse('apd_readout', laser_off_time + pi_half_pulse_time + tau + pi_half_pulse_time + delay_readout,
            #           meas_time)]
            #
            # pulse_sequence.append(Pulse('microwave_switch',
            #                             laser_off_time + pi_half_pulse_time + tau + pi_half_pulse_time + nv_reset_time + laser_off_time,
            #                             pi_half_pulse_time))
            #
            # end_of_first_MW = laser_off_time + pi_half_pulse_time + tau + pi_half_pulse_time + nv_reset_time + laser_off_time + pi_half_pulse_time
            #
            # pulse_sequence.append(Pulse('microwave_switch', end_of_first_MW + tau, pi_half_pulse_time))
            #
            # end_of_second_MW = end_of_first_MW + tau + pi_half_pulse_time
            #
            # pulse_sequence.append(Pulse('laser',
            #                             end_of_second_MW + tau + pi_half_pulse_time + delay_mw_readout, nv_reset_time))
            #
            # pulse_sequence.append(
            #     Pulse('apd_readout', end_of_second_MW + tau + pi_half_pulse_time + delay_mw_readout + delay_readout,
            #           meas_time))
            # pulse_sequence.append(
            #     Pulse('apd_switch', 0, end_of_second_MW + tau + pi_half_pulse_time + delay_mw_readout + nv_reset_time))

            # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
            # if tau == 0 or tau>=15:
            pulse_sequences.append(pulse_sequence)

        print('number of sequences before validation: ', len(pulse_sequences))
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

        # if 'fits' in data.keys() and data['fits'] is not None:
        #     counts = data['counts'][:,1]/ data['counts'][:,0]
        #     mean_fluor = np.mean(data['counts'][:,0])
        #     tau = data['tau']
        #     fits = data['fits'] # amplitude, frequency, phase, offset
        #
        #     # axislist[0].plot(tau, counts, 'b')
        #     axislist[0].plot(tau, counts, '.-')
        #
        #     #axislist[0].hold(True) ER 20181012
        #
        #     # axislist[0].plot(tau, cose_with_decay(tau, *fits), 'k', lw=3)
        #     axislist[0].plot(tau, cose_with_decay(tau, *fits), lw=3)
        #
        #
        #     # RabiT = 2 * np.pi / fits[1]
        #     T2_star = self.data['T2_star']
        #     phaseoffs = self.data['phaseoffs']
        #     pi_time = self.data['pi_time']
        #     pi_half_time = self.data['pi_half_time']
        #     three_pi_half_time = self.data['three_pi_half_time']
        #     # rabi_freq = 1000 * fits[1] / (2 * np.pi)
        #     rabi_freq = self.data['rabi_freq'] # Rabi frequency in [MHz]
        #
        #
        #     axislist[0].plot(pi_time, cose_with_decay(pi_time, *fits), 'ro', lw=3)
        #     axislist[0].annotate('$\pi$={:0.1f} ns'.format(pi_time), xy=(pi_time, cose_with_decay(pi_time, *fits)),
        #                          xytext=(pi_time - 10., cose_with_decay(pi_time, *fits) - .01), xycoords='data')
        #     axislist[0].plot(pi_half_time, cose_with_decay(pi_half_time, *fits), 'ro', lw=3)
        #     axislist[0].annotate('$\pi/2$={:0.1f} ns'.format(pi_half_time),
        #                          xy=(pi_half_time, cose_with_decay(pi_half_time, *fits)),
        #                          xytext=(pi_half_time + 10., cose_with_decay(pi_half_time, *fits)), xycoords='data')
        #     axislist[0].plot(three_pi_half_time, cose_with_decay(three_pi_half_time, *fits), 'ro', lw=3)
        #     axislist[0].annotate('$3\pi/2$={:0.1f} ns'.format(three_pi_half_time),
        #                          xy=(three_pi_half_time, cose_with_decay(three_pi_half_time, *fits)),
        #                          xytext=(three_pi_half_time + 10., cose_with_decay(three_pi_half_time, *fits)),
        #                          xycoords='data')
        #     axislist[0].plot(-phaseoffs, cose_with_decay(-phaseoffs, *fits), 'gd', lw=3)
        #     axislist[0].annotate('$start$={:0.1f} ns'.format(-phaseoffs),
        #                          xy=(-phaseoffs, cose_with_decay(-phaseoffs, *fits)),
        #                          xytext=(-phaseoffs + 10., cose_with_decay(-phaseoffs, *fits)), xycoords='data')
        #
        #     #pi_time = 2*np.pi / fits[1] / 2
        #     # pi_time = (np.pi - fits[2])/fits[1]
        #     # pi_half_time = (np.pi/2 - fits[2])/fits[1]
        #     # three_pi_half_time = (3*np.pi/2 - fits[2])/fits[1]
        #
        #  #   axislist[0].set_title('Rabi mw-power:{:0.1f}dBm, mw_freq:{:0.3f} GHz, pi-time: {:2.1f}ns, pi-half-time: {:2.1f}ns, 3pi_half_time: {:2.1f}ns, Rabi freq: {2.1f}MHz'.format(self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency']*1e-9, pi_time, pi_half_time, three_pi_half_time, rabi_freq))
        #     axislist[0].set_title(
        #         'Rabi: {:0.2f}MHz, pi-half time: {:2.1f}ns, pi-time: {:2.1f}ns, 3pi-half time: {:2.1f}ns \n T2*: {:2.1f}ns, Ref fluor: {:0.1f}kcps, mw-power: {:0.1f}dBm, mw_freq: {:0.3f} GHz'.format(
        #             rabi_freq, pi_half_time, pi_time, three_pi_half_time, T2_star, mean_fluor,
        #             self.settings['mw_pulses']['mw_power'],
        #             self.settings['mw_pulses'][
        #                 'mw_frequency'] * 1e-9))
        #     axislist[0].set_xlabel('Rabi tau [ns]')
        #     axislist[0].set_ylabel('Normalized Fluorescence')

        # elif 'exp_finished' in data.keys() is not None and data['exp_finished'] == 1:
        if 'exp_finished' in data.keys() is not None and data['exp_finished'] == 1:

            if 0 not in self.data['counts'][:, 0]:
                axislist[0].plot(data['tau'], data['counts'][:, 1] / data['counts'][:, 0])
                axislist[0].set_xlabel('Ramsey tau [ns]')
                axislist[0].set_ylabel('Normalized Fluorescence')
            else:
                axislist[0].plot(data['tau'], data['counts'][:, 0])
                axislist[0].plot(data['tau'], data['counts'][:, 1])
                axislist[0].set_xlabel('Ramsey tau [ns]')
                axislist[0].set_ylabel('Fluorescence [kcps]')
                axislist[0].legend(labels=('Ref', 'Ramsey data'), fontsize=8)

            # axislist[0].plot(data['tau'], data['counts'][:, 0])
            # axislist[0].plot(data['tau'], data['counts'][:, 1])
            # axislist[0].set_xlabel('Ramsey tau [ns]')
            # axislist[0].set_ylabel('Fluorescence [kcps]')

            axislist[0].set_title(
                '(final plot) R&S SMB100A \nRamsey mw-power:{:0.2f}dBm, mw_freq:{:0.4f} GHz \n detuning: {:0.4f} MHz, π/2 time: {:0.2f} ns'.format(
                    self.settings['mw_pulses']['mw_power'],
                    self.settings['mw_pulses']['resonant_freq'] * 1e-9, self.settings['mw_pulses']['detuning'] * 1e-6,
                    self.settings['mw_pulses']['pi_half_pulse_time']))
            axislist[0].legend(labels=('Ref', 'Ramsey data'), fontsize=8)

        else:

            super(Ramsey_RnS, self)._plot(axislist)
            axislist[0].set_title(
                '(initial plot) R&S SMB100A \nRamsey mw-power:{:0.2f}dBm, mw_freq:{:0.4f} GHz \n detuning: {:0.4f} MHz, π/2 time: {:0.2f} ns'.format(
                    self.settings['mw_pulses']['mw_power'],
                    self.settings['mw_pulses']['resonant_freq'] * 1e-9, self.settings['mw_pulses']['detuning'] * 1e-6,
                    self.settings['mw_pulses']['pi_half_pulse_time']))
            axislist[0].legend(labels=('Ref', 'Ramsey data'), fontsize=8)

    def _update_plot(self, axislist):
            # self._plot(axislist)
            if len(axislist[0].lines) == 0:
                self._plot(axislist)
                return
            super(Ramsey_RnS, self)._update_plot(axislist)

            axislist[0].set_title(
                '(updating plot) R&S SMB100A \nRamsey mw-power:{:0.2f}dBm, mw_freq:{:0.4f} GHz \n detuning: {:0.4f} MHz, π/2 time: {:0.2f} ns'.format(
                    self.settings['mw_pulses']['mw_power'],
                    self.settings['mw_pulses']['resonant_freq'] * 1e-9, self.settings['mw_pulses']['detuning'] * 1e-6,
                    self.settings['mw_pulses']['pi_half_pulse_time']))
            axislist[0].legend(labels=('Ref', 'Ramsey data'), fontsize=8)

class Ramsey_N9310A(PulsedExperimentBaseScript):
    """
        This script measures the FID signal of NV with Ramsey sequence to measure its dephasing time.
        Keysight signal generator N9310A is used.
        Note that the MW generator used in this script has I/Q modulation option.

        ==> Last edited by Ziwei 4/17/2019

    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_power', -45.0, float, 'microwave power in dBm'),
            Parameter('resonant_freq', 2.87e9, float, 'resonant frequency in Hz'),
            Parameter('detuning', -50e6, float, 'detuning for Ramsey experiment in Hz]'),
            Parameter('pi_half_pulse_time', 50.0, float, 'time duration of a pi/2 pulse (in ns)'),
            Parameter('three_pi_half_pulse_time', 135.0, float, 'time duration of a 3pi/2 pulse (in ns)'),
            Parameter('microwave_channel', 'i', ['i', 'q'], 'channel to use for mw pulses')
        ]),
        Parameter('tau_times', [
            Parameter('min_time', 20, float, 'minimum time between the two pi-half pulses (in ns)'),
            Parameter('max_time', 2000, float, 'maximum time between the two pi-half pulses (in ns)'),
            Parameter('time_step', 20., [2.5, 5., 10., 15., 20., 50., 100., 200., 500., 1000., 10000., 100000., 500000.],
                  'time step increment of rabi pulse duration (in ns)')
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 200, float, 'measurement time after rabi sequence (in ns)'),
            Parameter('nv_reset_time', 2000, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 600, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 550, int, '[ns] delay between laser on and APD readout (given by spontaneous decay rate)')
        ]),
        Parameter('min_pulse_dur', 15, [10, 15, 50], 'Minimum allowed pulse duration (ns) >10 ns'),
        Parameter('num_averages', 500000, int, 'number of averages (>100000)'),
    ]

    _INSTRUMENTS = {'daq': NI6353, 'PB': LISE607RTPulseBlaster, 'mw_gen_a': AgilentMicrowaveGenerator}

    def _function(self):

        freq_to_do = self.settings['mw_pulses']['resonant_freq'] + self.settings['mw_pulses']['detuning']
        ## Turn off green light (the pulse blaster will pulse it on when needed)
        self.instruments['PB']['instance'].update({'laser': {'status': False}})
        self.instruments['PB']['instance'].update({'apd_switch': {'status': False}})

        self.instruments['mw_gen_a']['instance'].update({'enable_IQ': True})
        self.instruments['mw_gen_a']['instance'].update({'enable_modulation': True})
        self.instruments['mw_gen_a']['instance'].update({'freq_mode': 'CW'})
        self.instruments['mw_gen_a']['instance'].update({'power_mode': 'CW'})
        self.instruments['mw_gen_a']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power']})
        self.instruments['mw_gen_a']['instance'].update({'frequency': freq_to_do})
        self.instruments['mw_gen_a']['instance'].update({'enable_output': True})

        self.avg_block_number = super(Ramsey_N9310A, self)._function()

        self.data['exp_finished'] = 1 # this is to flag when the experiment is done or stopped. this allows proper saving even fit fails
        self.data['avg_block_number'] = self.avg_block_number

        # turn off laser, apd_switch and MW (just in case)
        self.instruments['PB']['instance'].update({'laser': {'status': False}})
        self.instruments['PB']['instance'].update({'apd_switch': {'status': False}})
        self.instruments['mw_gen_a']['instance'].update({'enable_output': False})
        self.instruments['mw_gen_a']['instance'].update({'enable_IQ': False})
        self.instruments['mw_gen_a']['instance'].update({'enable_modulation': False})

        # start fitting
        # if 'counts' in self.data.keys() and 'tau' in self.data.keys() and not 0 in self.data['counts'][:, 0]:
        #     counts = self.data['counts'][:, 1] / self.data['counts'][:, 0]
        #     tau = self.data['tau']
        #     fit_success = False
            # try:
            #     rabi_fits = fit_rabi_decay(tau, counts, variable_phase=True)
            #     if rabi_fits[1] == True:
            #         self.data['fits'] = rabi_fits[0]
            #         fit_success = True
            #     else:
            #         # self.data['fits'] = None
            #         fit_success = False
            #         print('rabi fit failed')
            #         self.log('rabi fit failed')
            # except:
            #     # self.data['fits'] = None
            #     fit_success = False
            #     print('rabi fit failed')
            #     self.log('rabi fit failed')
            #
            # if fit_success:
            #     counts = self.data['counts'][:,1]/ self.data['counts'][:,0]
            #     tau =self.data['tau']
            #     fits = self.data['fits']
            #
            #     RabiT = 2*np.pi / fits[1]
            #     # RabiT = (2 * np.pi / fits[1]) # avoid negative rabi period ZQ 8/25/2017
            #     T2_star = fits[4]
            #     phaseoffs = fits[2]
            #     pi_time = RabiT / 2 - phaseoffs * RabiT / (2 * np.pi)
            #     pi_half_time = RabiT / 4 - phaseoffs * RabiT / (2 * np.pi)
            #     three_pi_half_time = 3 * RabiT / 4 - phaseoffs * RabiT / (2 * np.pi)
            #     # pi_time = RabiT / 2
            #     # pi_half_time = RabiT / 4
            #     # three_pi_half_time = 3 * RabiT / 4
            #     self.data['phaseoffs'] = phaseoffs * RabiT / (2 * np.pi)
            #     self.data['pi_time'] = pi_time
            #     self.data['pi_half_time'] = pi_half_time
            #     self.data['three_pi_half_time'] = three_pi_half_time
            #     self.data['T2_star'] = T2_star
            #     self.data['rabi_freq'] = 1000 * fits[1] / (2 * np.pi) # Rabi frequency in [MHz]

    def _create_pulse_sequences(self):
        """

        Returns: pulse_sequences, num_averages, tau_list, meas_time
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
            meas_time: the width (in ns) of the daq measurement

        """
        MIN_DURATION = self.settings['min_pulse_dur']
        pulse_sequences = []
        max_range = int(np.floor((self.settings['tau_times']['max_time']-self.settings['tau_times']['min_time'])/self.settings['tau_times']['time_step']))
        tau_list = np.array([self.settings['tau_times']['min_time'] + i*self.settings['tau_times']['time_step'] for i in range(max_range)])

        # ignore the sequence if the mw-pulse is shorter than MIN_DURATION (0 is ok because there is no mw pulse!)
        tau_list = [x for x in tau_list if x == 0 or x >= MIN_DURATION]
        print('tau_list in [ns]:', tau_list)

        pi_half_pulse_time = self.settings['mw_pulses']['pi_half_pulse_time']
        three_pi_half_pulse_time = self.settings['mw_pulses']['three_pi_half_pulse_time']
        meas_time = self.settings['read_out']['meas_time']
        nv_reset_time = self.settings['read_out']['nv_reset_time']
        laser_off_time = self.settings['read_out']['laser_off_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        delay_readout = self.settings['read_out']['delay_readout']
        # mw_sw_buffer = self.settings['read_out']['mw_switch_extra_time']
        # APD_sw_delay = self.settings['read_out']['APD_switch_delay']
        microwave_channel = 'microwave_' + self.settings['mw_pulses']['microwave_channel']

        for tau in tau_list:
            # Normal Ramsey sequence
            pulse_sequence = [Pulse('laser', laser_off_time + pi_half_pulse_time + tau + pi_half_pulse_time, nv_reset_time),
                              Pulse('apd_readout', laser_off_time + pi_half_pulse_time + tau + pi_half_pulse_time + delay_readout, meas_time),
                              ]

            end_of_first_laser = laser_off_time + pi_half_pulse_time + tau + pi_half_pulse_time + nv_reset_time

            pulse_sequence.append(Pulse(microwave_channel, end_of_first_laser + laser_off_time, pi_half_pulse_time))
            pulse_sequence.append(Pulse(microwave_channel,
                                        end_of_first_laser + laser_off_time + pi_half_pulse_time + tau, pi_half_pulse_time))
            pulse_sequence.append(Pulse('laser',
                                        end_of_first_laser + laser_off_time + pi_half_pulse_time + tau + pi_half_pulse_time + delay_mw_readout,
                                        nv_reset_time))
            pulse_sequence.append(Pulse('apd_readout',
                                        end_of_first_laser + laser_off_time + pi_half_pulse_time + tau + pi_half_pulse_time + delay_mw_readout + delay_readout,
                                        meas_time))

            pulse_sequence.append(
                Pulse('apd_switch', 0, end_of_first_laser + laser_off_time + pi_half_pulse_time + tau + pi_half_pulse_time + delay_mw_readout + nv_reset_time))

            # The sequence implementing 3pi/2 readout...
            # start_of_first_Ramsey = laser_off_time
            # # the first Ramsey with pi/2 readout
            # pulse_sequence = [Pulse('microwave_switch', start_of_first_Ramsey, pi_half_pulse_time)]
            # pulse_sequence.append(Pulse('microwave_switch', start_of_first_Ramsey + pi_half_pulse_time + tau, pi_half_pulse_time))
            # pulse_sequence.append(Pulse('laser',
            #                             start_of_first_Ramsey + pi_half_pulse_time + tau + pi_half_pulse_time + delay_mw_readout,
            #                             nv_reset_time))
            # pulse_sequence.append(Pulse('apd_readout',
            #                             start_of_first_Ramsey + pi_half_pulse_time + tau + pi_half_pulse_time + delay_mw_readout + delay_readout,
            #                             meas_time))
            #
            # start_of_second_Ramsey = start_of_first_Ramsey + pi_half_pulse_time + tau + pi_half_pulse_time + delay_mw_readout + nv_reset_time + laser_off_time
            # # the second Ramsey with 3pi/2 readout
            # pulse_sequence.append(Pulse('microwave_switch', start_of_second_Ramsey, pi_half_pulse_time))
            # pulse_sequence.append(
            #     Pulse('microwave_switch', start_of_second_Ramsey + pi_half_pulse_time + tau, three_pi_half_pulse_time))
            # pulse_sequence.append(Pulse('laser',
            #                             start_of_second_Ramsey + pi_half_pulse_time + tau + three_pi_half_pulse_time + delay_mw_readout,
            #                             nv_reset_time))
            # pulse_sequence.append(Pulse('apd_readout',
            #                             start_of_second_Ramsey + pi_half_pulse_time + tau + three_pi_half_pulse_time + delay_mw_readout + delay_readout,
            #                             meas_time))
            # pulse_sequence.append(
            #     Pulse('apd_switch', 0, start_of_second_Ramsey + pi_half_pulse_time + tau + three_pi_half_pulse_time + delay_mw_readout + nv_reset_time))


            # pulse_sequence = [
            #     Pulse('laser', laser_off_time + pi_half_pulse_time + tau + pi_half_pulse_time, nv_reset_time),
            #     Pulse('apd_readout', laser_off_time + pi_half_pulse_time + tau + pi_half_pulse_time + delay_readout,
            #           meas_time)]
            #
            # pulse_sequence.append(Pulse('microwave_switch',
            #                             laser_off_time + pi_half_pulse_time + tau + pi_half_pulse_time + nv_reset_time + laser_off_time,
            #                             pi_half_pulse_time))
            #
            # end_of_first_MW = laser_off_time + pi_half_pulse_time + tau + pi_half_pulse_time + nv_reset_time + laser_off_time + pi_half_pulse_time
            #
            # pulse_sequence.append(Pulse('microwave_switch', end_of_first_MW + tau, pi_half_pulse_time))
            #
            # end_of_second_MW = end_of_first_MW + tau + pi_half_pulse_time
            #
            # pulse_sequence.append(Pulse('laser',
            #                             end_of_second_MW + tau + pi_half_pulse_time + delay_mw_readout, nv_reset_time))
            #
            # pulse_sequence.append(
            #     Pulse('apd_readout', end_of_second_MW + tau + pi_half_pulse_time + delay_mw_readout + delay_readout,
            #           meas_time))
            # pulse_sequence.append(
            #     Pulse('apd_switch', 0, end_of_second_MW + tau + pi_half_pulse_time + delay_mw_readout + nv_reset_time))

            # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
            # if tau == 0 or tau>=15:
            pulse_sequences.append(pulse_sequence)

        print('number of sequences before validation: ', len(pulse_sequences))
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

        if 'avg_block_number' in data.keys() is not None and data['avg_block_number'] is not None:
            avg_block_number = data['avg_block_number']
        else:
            avg_block_number = 0

        # if 'fits' in data.keys() and data['fits'] is not None:
        #     counts = data['counts'][:,1]/ data['counts'][:,0]
        #     mean_fluor = np.mean(data['counts'][:,0])
        #     tau = data['tau']
        #     fits = data['fits'] # amplitude, frequency, phase, offset
        #
        #     # axislist[0].plot(tau, counts, 'b')
        #     axislist[0].plot(tau, counts, '.-')
        #
        #     #axislist[0].hold(True) ER 20181012
        #
        #     # axislist[0].plot(tau, cose_with_decay(tau, *fits), 'k', lw=3)
        #     axislist[0].plot(tau, cose_with_decay(tau, *fits), lw=3)
        #
        #
        #     # RabiT = 2 * np.pi / fits[1]
        #     T2_star = self.data['T2_star']
        #     phaseoffs = self.data['phaseoffs']
        #     pi_time = self.data['pi_time']
        #     pi_half_time = self.data['pi_half_time']
        #     three_pi_half_time = self.data['three_pi_half_time']
        #     # rabi_freq = 1000 * fits[1] / (2 * np.pi)
        #     rabi_freq = self.data['rabi_freq'] # Rabi frequency in [MHz]
        #
        #
        #     axislist[0].plot(pi_time, cose_with_decay(pi_time, *fits), 'ro', lw=3)
        #     axislist[0].annotate('$\pi$={:0.1f} ns'.format(pi_time), xy=(pi_time, cose_with_decay(pi_time, *fits)),
        #                          xytext=(pi_time - 10., cose_with_decay(pi_time, *fits) - .01), xycoords='data')
        #     axislist[0].plot(pi_half_time, cose_with_decay(pi_half_time, *fits), 'ro', lw=3)
        #     axislist[0].annotate('$\pi/2$={:0.1f} ns'.format(pi_half_time),
        #                          xy=(pi_half_time, cose_with_decay(pi_half_time, *fits)),
        #                          xytext=(pi_half_time + 10., cose_with_decay(pi_half_time, *fits)), xycoords='data')
        #     axislist[0].plot(three_pi_half_time, cose_with_decay(three_pi_half_time, *fits), 'ro', lw=3)
        #     axislist[0].annotate('$3\pi/2$={:0.1f} ns'.format(three_pi_half_time),
        #                          xy=(three_pi_half_time, cose_with_decay(three_pi_half_time, *fits)),
        #                          xytext=(three_pi_half_time + 10., cose_with_decay(three_pi_half_time, *fits)),
        #                          xycoords='data')
        #     axislist[0].plot(-phaseoffs, cose_with_decay(-phaseoffs, *fits), 'gd', lw=3)
        #     axislist[0].annotate('$start$={:0.1f} ns'.format(-phaseoffs),
        #                          xy=(-phaseoffs, cose_with_decay(-phaseoffs, *fits)),
        #                          xytext=(-phaseoffs + 10., cose_with_decay(-phaseoffs, *fits)), xycoords='data')
        #
        #     #pi_time = 2*np.pi / fits[1] / 2
        #     # pi_time = (np.pi - fits[2])/fits[1]
        #     # pi_half_time = (np.pi/2 - fits[2])/fits[1]
        #     # three_pi_half_time = (3*np.pi/2 - fits[2])/fits[1]
        #
        #  #   axislist[0].set_title('Rabi mw-power:{:0.1f}dBm, mw_freq:{:0.3f} GHz, pi-time: {:2.1f}ns, pi-half-time: {:2.1f}ns, 3pi_half_time: {:2.1f}ns, Rabi freq: {2.1f}MHz'.format(self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency']*1e-9, pi_time, pi_half_time, three_pi_half_time, rabi_freq))
        #     axislist[0].set_title(
        #         'Rabi: {:0.2f}MHz, pi-half time: {:2.1f}ns, pi-time: {:2.1f}ns, 3pi-half time: {:2.1f}ns \n T2*: {:2.1f}ns, Ref fluor: {:0.1f}kcps, mw-power: {:0.1f}dBm, mw_freq: {:0.3f} GHz'.format(
        #             rabi_freq, pi_half_time, pi_time, three_pi_half_time, T2_star, mean_fluor,
        #             self.settings['mw_pulses']['mw_power'],
        #             self.settings['mw_pulses'][
        #                 'mw_frequency'] * 1e-9))
        #     axislist[0].set_xlabel('Rabi tau [ns]')
        #     axislist[0].set_ylabel('Normalized Fluorescence')

        # elif 'exp_finished' in data.keys() is not None and data['exp_finished'] == 1:
        if 'exp_finished' in data.keys() is not None and data['exp_finished'] == 1:

            if 0 not in self.data['counts'][:, 0]:
                axislist[0].plot(data['tau'], data['counts'][:, 1] / data['counts'][:, 0])
                axislist[0].set_xlabel('Ramsey tau [ns]')
                axislist[0].set_ylabel('Normalized Fluorescence')
                axislist[0].legend(labels=('Ramsey data'), fontsize=8)
            else:
                axislist[0].plot(data['tau'], data['counts'][:, 0])
                axislist[0].plot(data['tau'], data['counts'][:, 1])
                axislist[0].set_xlabel('Ramsey tau [ns]')
                axislist[0].set_ylabel('Fluorescence [kcps]')
                axislist[0].legend(labels=('Ref', 'Ramsey data'), fontsize=8)

            # axislist[0].plot(data['tau'], data['counts'][:, 0])
            # axislist[0].plot(data['tau'], data['counts'][:, 1])
            # axislist[0].set_xlabel('Ramsey tau [ns]')
            # axislist[0].set_ylabel('Fluorescence [kcps]')

            axislist[0].set_title(
                '(final plot) Keysight N9310A, pulse carved by ' + self.settings['mw_pulses'][
                    'microwave_channel'] + '\nRamsey mw-power:{:0.2f}dBm, mw_freq:{:0.4f} GHz \nafter {:d} averages\ndetuning: {:0.4f} MHz, π/2 time: {:0.2f} ns'.format(
                    self.settings['mw_pulses']['mw_power'],
                    self.settings['mw_pulses']['resonant_freq'] * 1e-9, avg_block_number, self.settings['mw_pulses']['detuning'] * 1e-6,
                    self.settings['mw_pulses']['pi_half_pulse_time']))
            axislist[0].legend(labels=('Ref', 'Ramsey data'), fontsize=8)

        else:

            super(Ramsey_N9310A, self)._plot(axislist)
            axislist[0].set_title(
                '(initial plot) Keysight N9310A, pulse carved by ' + self.settings['mw_pulses'][
                    'microwave_channel'] + '\nRamsey mw-power:{:0.2f}dBm, mw_freq:{:0.4f} GHz \nafter {:d} averages\n detuning: {:0.4f} MHz, π/2 time: {:0.2f} ns'.format(
                    self.settings['mw_pulses']['mw_power'],
                    self.settings['mw_pulses']['resonant_freq'] * 1e-9, avg_block_number, self.settings['mw_pulses']['detuning'] * 1e-6,
                    self.settings['mw_pulses']['pi_half_pulse_time']))
            axislist[0].legend(labels=('Ref', 'Ramsey data'), fontsize=8)

    def _update_plot(self, axislist):
            # self._plot(axislist)
            if len(axislist[0].lines) == 0:
                self._plot(axislist)
                return
            super(Ramsey_N9310A, self)._update_plot(axislist)

            axislist[0].set_title(
                '(updating plot) Keysight N9310A, pulse carved by ' + self.settings['mw_pulses'][
                    'microwave_channel'] + '\nRamsey mw-power:{:0.2f}dBm, mw_freq:{:0.4f} GHz \n detuning: {:0.4f} MHz, π/2 time: {:0.2f} ns'.format(
                    self.settings['mw_pulses']['mw_power'],
                    self.settings['mw_pulses']['resonant_freq'] * 1e-9, self.settings['mw_pulses']['detuning'] * 1e-6,
                    self.settings['mw_pulses']['pi_half_pulse_time']))
            axislist[0].legend(labels=('Ref', 'Ramsey data'), fontsize=8)

class Ramsey_CoilVol(PulsedExperimentBaseScript):
    """
        This script measures the FID signal of NV with Ramsey sequence to measure its dephasing time.
        Keysight signal generator N9310A is used.
        The script is essentially the same as Ramsey_N9310A, but it allows controlling the voltage on the coil
        Note that the MW generator used in this script has I/Q modulation option.

        ==> Last edited by Ziwei 11/23/2019

    """
    _DEFAULT_SETTINGS = [
        Parameter('coil_voltage', 0.0, float, 'voltage to be applied on the coil'),
        Parameter('mw_pulses', [
            Parameter('mw_power', -45.0, float, 'microwave power in dBm'),
            Parameter('resonant_freq', 2.87e9, float, 'resonant frequency in Hz'),
            Parameter('detuning', -50e6, float, 'detuning for Ramsey experiment in Hz]'),
            Parameter('pi_half_pulse_time', 50.0, float, 'time duration of a pi/2 pulse (in ns)'),
            Parameter('three_pi_half_pulse_time', 135.0, float, 'time duration of a 3pi/2 pulse (in ns)'),
            Parameter('microwave_channel', 'i', ['i', 'q'], 'channel to use for mw pulses')
        ]),
        Parameter('tau_times', [
            Parameter('min_time', 20, float, 'minimum time between the two pi-half pulses (in ns)'),
            Parameter('max_time', 2000, float, 'maximum time between the two pi-half pulses (in ns)'),
            Parameter('time_step', 20., [2.5, 5., 10., 20., 50., 100., 200., 500., 1000., 10000., 100000., 500000.],
                  'time step increment of rabi pulse duration (in ns)')
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 200, float, 'measurement time after rabi sequence (in ns)'),
            Parameter('nv_reset_time', 2000, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 600, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 550, int, '[ns] delay between laser on and APD readout (given by spontaneous decay rate)')
        ]),
        Parameter('min_pulse_dur', 15, [10, 15, 50], 'Minimum allowed pulse duration (ns) >10 ns'),
        Parameter('num_averages', 500000, int, 'number of averages (>100000)'),
    ]

    _INSTRUMENTS = {'daq': NI6353, 'PB': LISE607RTPulseBlaster, 'mw_gen_a': AgilentMicrowaveGenerator,
                    'fg': Agilent33120A}

    def _function(self):

        # Set up the function generator
        self.instruments['fg']['instance'].update({'display_on': True})
        self.instruments['fg']['instance'].update({'burst_mod': False})
        self.instruments['fg']['instance'].update({'wave_shape': 'DC'})
        print('Setting coil voltage to {:0.4f} V'.format(self.settings['coil_voltage']))
        self.instruments['fg']['instance'].update({'offset': self.settings['coil_voltage']})
        time.sleep(0.1)

        freq_to_do = self.settings['mw_pulses']['resonant_freq'] + self.settings['mw_pulses']['detuning']
        ## Turn off green light (the pulse blaster will pulse it on when needed)
        self.instruments['PB']['instance'].update({'laser': {'status': False}})
        self.instruments['PB']['instance'].update({'apd_switch': {'status': False}})

        self.instruments['mw_gen_a']['instance'].update({'enable_IQ': True})
        self.instruments['mw_gen_a']['instance'].update({'enable_modulation': True})
        self.instruments['mw_gen_a']['instance'].update({'freq_mode': 'CW'})
        self.instruments['mw_gen_a']['instance'].update({'power_mode': 'CW'})
        self.instruments['mw_gen_a']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power']})
        self.instruments['mw_gen_a']['instance'].update({'frequency': freq_to_do})
        self.instruments['mw_gen_a']['instance'].update({'enable_output': True})

        self.avg_block_number = super(Ramsey_CoilVol, self)._function()

        self.data['exp_finished'] = 1 # this is to flag when the experiment is done or stopped. this allows proper saving even fit fails
        self.data['avg_block_number'] = self.avg_block_number

        # turn off laser, apd_switch and MW (just in case)
        self.instruments['PB']['instance'].update({'laser': {'status': False}})
        self.instruments['PB']['instance'].update({'apd_switch': {'status': False}})
        self.instruments['mw_gen_a']['instance'].update({'enable_output': False})
        self.instruments['mw_gen_a']['instance'].update({'enable_IQ': False})
        self.instruments['mw_gen_a']['instance'].update({'enable_modulation': False})

        self.instruments['fg']['instance'].update({'offset': 0.0})
        print('Setting coil voltage to 0.0 V')
        time.sleep(0.1)

        # start fitting
        # if 'counts' in self.data.keys() and 'tau' in self.data.keys() and not 0 in self.data['counts'][:, 0]:
        #     counts = self.data['counts'][:, 1] / self.data['counts'][:, 0]
        #     tau = self.data['tau']
        #     fit_success = False
            # try:
            #     rabi_fits = fit_rabi_decay(tau, counts, variable_phase=True)
            #     if rabi_fits[1] == True:
            #         self.data['fits'] = rabi_fits[0]
            #         fit_success = True
            #     else:
            #         # self.data['fits'] = None
            #         fit_success = False
            #         print('rabi fit failed')
            #         self.log('rabi fit failed')
            # except:
            #     # self.data['fits'] = None
            #     fit_success = False
            #     print('rabi fit failed')
            #     self.log('rabi fit failed')
            #
            # if fit_success:
            #     counts = self.data['counts'][:,1]/ self.data['counts'][:,0]
            #     tau =self.data['tau']
            #     fits = self.data['fits']
            #
            #     RabiT = 2*np.pi / fits[1]
            #     # RabiT = (2 * np.pi / fits[1]) # avoid negative rabi period ZQ 8/25/2017
            #     T2_star = fits[4]
            #     phaseoffs = fits[2]
            #     pi_time = RabiT / 2 - phaseoffs * RabiT / (2 * np.pi)
            #     pi_half_time = RabiT / 4 - phaseoffs * RabiT / (2 * np.pi)
            #     three_pi_half_time = 3 * RabiT / 4 - phaseoffs * RabiT / (2 * np.pi)
            #     # pi_time = RabiT / 2
            #     # pi_half_time = RabiT / 4
            #     # three_pi_half_time = 3 * RabiT / 4
            #     self.data['phaseoffs'] = phaseoffs * RabiT / (2 * np.pi)
            #     self.data['pi_time'] = pi_time
            #     self.data['pi_half_time'] = pi_half_time
            #     self.data['three_pi_half_time'] = three_pi_half_time
            #     self.data['T2_star'] = T2_star
            #     self.data['rabi_freq'] = 1000 * fits[1] / (2 * np.pi) # Rabi frequency in [MHz]

    def _create_pulse_sequences(self):
        """

        Returns: pulse_sequences, num_averages, tau_list, meas_time
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
            meas_time: the width (in ns) of the daq measurement

        """
        MIN_DURATION = self.settings['min_pulse_dur']
        pulse_sequences = []
        max_range = int(np.floor((self.settings['tau_times']['max_time']-self.settings['tau_times']['min_time'])/self.settings['tau_times']['time_step']))
        tau_list = np.array([self.settings['tau_times']['min_time'] + i*self.settings['tau_times']['time_step'] for i in range(max_range)])

        # ignore the sequence if the mw-pulse is shorter than MIN_DURATION (0 is ok because there is no mw pulse!)
        tau_list = [x for x in tau_list if x == 0 or x >= MIN_DURATION]
        print('tau_list in [ns]:', tau_list)

        pi_half_pulse_time = self.settings['mw_pulses']['pi_half_pulse_time']
        three_pi_half_pulse_time = self.settings['mw_pulses']['three_pi_half_pulse_time']
        meas_time = self.settings['read_out']['meas_time']
        nv_reset_time = self.settings['read_out']['nv_reset_time']
        laser_off_time = self.settings['read_out']['laser_off_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        delay_readout = self.settings['read_out']['delay_readout']
        # mw_sw_buffer = self.settings['read_out']['mw_switch_extra_time']
        # APD_sw_delay = self.settings['read_out']['APD_switch_delay']
        microwave_channel = 'microwave_' + self.settings['mw_pulses']['microwave_channel']

        for tau in tau_list:
            # Normal Ramsey sequence
            pulse_sequence = [Pulse('laser', laser_off_time + pi_half_pulse_time + tau + pi_half_pulse_time, nv_reset_time),
                              Pulse('apd_readout', laser_off_time + pi_half_pulse_time + tau + pi_half_pulse_time + delay_readout, meas_time),
                              ]

            end_of_first_laser = laser_off_time + pi_half_pulse_time + tau + pi_half_pulse_time + nv_reset_time

            pulse_sequence.append(Pulse(microwave_channel, end_of_first_laser + laser_off_time, pi_half_pulse_time))
            pulse_sequence.append(Pulse(microwave_channel,
                                        end_of_first_laser + laser_off_time + pi_half_pulse_time + tau, pi_half_pulse_time))
            pulse_sequence.append(Pulse('laser',
                                        end_of_first_laser + laser_off_time + pi_half_pulse_time + tau + pi_half_pulse_time + delay_mw_readout,
                                        nv_reset_time))
            pulse_sequence.append(Pulse('apd_readout',
                                        end_of_first_laser + laser_off_time + pi_half_pulse_time + tau + pi_half_pulse_time + delay_mw_readout + delay_readout,
                                        meas_time))

            pulse_sequence.append(
                Pulse('apd_switch', 0, end_of_first_laser + laser_off_time + pi_half_pulse_time + tau + pi_half_pulse_time + delay_mw_readout + nv_reset_time))

            # The sequence implementing 3pi/2 readout...
            # start_of_first_Ramsey = laser_off_time
            # # the first Ramsey with pi/2 readout
            # pulse_sequence = [Pulse('microwave_switch', start_of_first_Ramsey, pi_half_pulse_time)]
            # pulse_sequence.append(Pulse('microwave_switch', start_of_first_Ramsey + pi_half_pulse_time + tau, pi_half_pulse_time))
            # pulse_sequence.append(Pulse('laser',
            #                             start_of_first_Ramsey + pi_half_pulse_time + tau + pi_half_pulse_time + delay_mw_readout,
            #                             nv_reset_time))
            # pulse_sequence.append(Pulse('apd_readout',
            #                             start_of_first_Ramsey + pi_half_pulse_time + tau + pi_half_pulse_time + delay_mw_readout + delay_readout,
            #                             meas_time))
            #
            # start_of_second_Ramsey = start_of_first_Ramsey + pi_half_pulse_time + tau + pi_half_pulse_time + delay_mw_readout + nv_reset_time + laser_off_time
            # # the second Ramsey with 3pi/2 readout
            # pulse_sequence.append(Pulse('microwave_switch', start_of_second_Ramsey, pi_half_pulse_time))
            # pulse_sequence.append(
            #     Pulse('microwave_switch', start_of_second_Ramsey + pi_half_pulse_time + tau, three_pi_half_pulse_time))
            # pulse_sequence.append(Pulse('laser',
            #                             start_of_second_Ramsey + pi_half_pulse_time + tau + three_pi_half_pulse_time + delay_mw_readout,
            #                             nv_reset_time))
            # pulse_sequence.append(Pulse('apd_readout',
            #                             start_of_second_Ramsey + pi_half_pulse_time + tau + three_pi_half_pulse_time + delay_mw_readout + delay_readout,
            #                             meas_time))
            # pulse_sequence.append(
            #     Pulse('apd_switch', 0, start_of_second_Ramsey + pi_half_pulse_time + tau + three_pi_half_pulse_time + delay_mw_readout + nv_reset_time))


            # pulse_sequence = [
            #     Pulse('laser', laser_off_time + pi_half_pulse_time + tau + pi_half_pulse_time, nv_reset_time),
            #     Pulse('apd_readout', laser_off_time + pi_half_pulse_time + tau + pi_half_pulse_time + delay_readout,
            #           meas_time)]
            #
            # pulse_sequence.append(Pulse('microwave_switch',
            #                             laser_off_time + pi_half_pulse_time + tau + pi_half_pulse_time + nv_reset_time + laser_off_time,
            #                             pi_half_pulse_time))
            #
            # end_of_first_MW = laser_off_time + pi_half_pulse_time + tau + pi_half_pulse_time + nv_reset_time + laser_off_time + pi_half_pulse_time
            #
            # pulse_sequence.append(Pulse('microwave_switch', end_of_first_MW + tau, pi_half_pulse_time))
            #
            # end_of_second_MW = end_of_first_MW + tau + pi_half_pulse_time
            #
            # pulse_sequence.append(Pulse('laser',
            #                             end_of_second_MW + tau + pi_half_pulse_time + delay_mw_readout, nv_reset_time))
            #
            # pulse_sequence.append(
            #     Pulse('apd_readout', end_of_second_MW + tau + pi_half_pulse_time + delay_mw_readout + delay_readout,
            #           meas_time))
            # pulse_sequence.append(
            #     Pulse('apd_switch', 0, end_of_second_MW + tau + pi_half_pulse_time + delay_mw_readout + nv_reset_time))

            # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
            # if tau == 0 or tau>=15:
            pulse_sequences.append(pulse_sequence)

        print('number of sequences before validation: ', len(pulse_sequences))
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

        if 'avg_block_number' in data.keys() is not None and data['avg_block_number'] is not None:
            avg_block_number = data['avg_block_number']
        else:
            avg_block_number = 0

        # if 'fits' in data.keys() and data['fits'] is not None:
        #     counts = data['counts'][:,1]/ data['counts'][:,0]
        #     mean_fluor = np.mean(data['counts'][:,0])
        #     tau = data['tau']
        #     fits = data['fits'] # amplitude, frequency, phase, offset
        #
        #     # axislist[0].plot(tau, counts, 'b')
        #     axislist[0].plot(tau, counts, '.-')
        #
        #     #axislist[0].hold(True) ER 20181012
        #
        #     # axislist[0].plot(tau, cose_with_decay(tau, *fits), 'k', lw=3)
        #     axislist[0].plot(tau, cose_with_decay(tau, *fits), lw=3)
        #
        #
        #     # RabiT = 2 * np.pi / fits[1]
        #     T2_star = self.data['T2_star']
        #     phaseoffs = self.data['phaseoffs']
        #     pi_time = self.data['pi_time']
        #     pi_half_time = self.data['pi_half_time']
        #     three_pi_half_time = self.data['three_pi_half_time']
        #     # rabi_freq = 1000 * fits[1] / (2 * np.pi)
        #     rabi_freq = self.data['rabi_freq'] # Rabi frequency in [MHz]
        #
        #
        #     axislist[0].plot(pi_time, cose_with_decay(pi_time, *fits), 'ro', lw=3)
        #     axislist[0].annotate('$\pi$={:0.1f} ns'.format(pi_time), xy=(pi_time, cose_with_decay(pi_time, *fits)),
        #                          xytext=(pi_time - 10., cose_with_decay(pi_time, *fits) - .01), xycoords='data')
        #     axislist[0].plot(pi_half_time, cose_with_decay(pi_half_time, *fits), 'ro', lw=3)
        #     axislist[0].annotate('$\pi/2$={:0.1f} ns'.format(pi_half_time),
        #                          xy=(pi_half_time, cose_with_decay(pi_half_time, *fits)),
        #                          xytext=(pi_half_time + 10., cose_with_decay(pi_half_time, *fits)), xycoords='data')
        #     axislist[0].plot(three_pi_half_time, cose_with_decay(three_pi_half_time, *fits), 'ro', lw=3)
        #     axislist[0].annotate('$3\pi/2$={:0.1f} ns'.format(three_pi_half_time),
        #                          xy=(three_pi_half_time, cose_with_decay(three_pi_half_time, *fits)),
        #                          xytext=(three_pi_half_time + 10., cose_with_decay(three_pi_half_time, *fits)),
        #                          xycoords='data')
        #     axislist[0].plot(-phaseoffs, cose_with_decay(-phaseoffs, *fits), 'gd', lw=3)
        #     axislist[0].annotate('$start$={:0.1f} ns'.format(-phaseoffs),
        #                          xy=(-phaseoffs, cose_with_decay(-phaseoffs, *fits)),
        #                          xytext=(-phaseoffs + 10., cose_with_decay(-phaseoffs, *fits)), xycoords='data')
        #
        #     #pi_time = 2*np.pi / fits[1] / 2
        #     # pi_time = (np.pi - fits[2])/fits[1]
        #     # pi_half_time = (np.pi/2 - fits[2])/fits[1]
        #     # three_pi_half_time = (3*np.pi/2 - fits[2])/fits[1]
        #
        #  #   axislist[0].set_title('Rabi mw-power:{:0.1f}dBm, mw_freq:{:0.3f} GHz, pi-time: {:2.1f}ns, pi-half-time: {:2.1f}ns, 3pi_half_time: {:2.1f}ns, Rabi freq: {2.1f}MHz'.format(self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency']*1e-9, pi_time, pi_half_time, three_pi_half_time, rabi_freq))
        #     axislist[0].set_title(
        #         'Rabi: {:0.2f}MHz, pi-half time: {:2.1f}ns, pi-time: {:2.1f}ns, 3pi-half time: {:2.1f}ns \n T2*: {:2.1f}ns, Ref fluor: {:0.1f}kcps, mw-power: {:0.1f}dBm, mw_freq: {:0.3f} GHz'.format(
        #             rabi_freq, pi_half_time, pi_time, three_pi_half_time, T2_star, mean_fluor,
        #             self.settings['mw_pulses']['mw_power'],
        #             self.settings['mw_pulses'][
        #                 'mw_frequency'] * 1e-9))
        #     axislist[0].set_xlabel('Rabi tau [ns]')
        #     axislist[0].set_ylabel('Normalized Fluorescence')

        # elif 'exp_finished' in data.keys() is not None and data['exp_finished'] == 1:
        if 'exp_finished' in data.keys() is not None and data['exp_finished'] == 1:

            if 0 not in self.data['counts'][:, 0]:
                axislist[0].plot(data['tau'], data['counts'][:, 1] / data['counts'][:, 0])
                axislist[0].set_xlabel('Ramsey tau [ns]')
                axislist[0].set_ylabel('Normalized Fluorescence')
                axislist[0].legend(labels=('Ramsey data'), fontsize=8)
            else:
                axislist[0].plot(data['tau'], data['counts'][:, 0])
                axislist[0].plot(data['tau'], data['counts'][:, 1])
                axislist[0].set_xlabel('Ramsey tau [ns]')
                axislist[0].set_ylabel('Fluorescence [kcps]')
                axislist[0].legend(labels=('Ref', 'Ramsey data'), fontsize=8)

            axislist[0].set_title(
                '(final plot) Keysight N9310A, pulse carved by ' + self.settings['mw_pulses'][
                    'microwave_channel'] + '\nRamsey mw-power:{:0.2f}dBm, mw_freq:{:0.4f} GHz \nafter {:d} averages\ndetuning: {:0.4f} MHz, π/2 time: {:0.2f} ns\n coil voltage = {:0.4f} V'.format(
                    self.settings['mw_pulses']['mw_power'],
                    self.settings['mw_pulses']['resonant_freq'] * 1e-9, avg_block_number, self.settings['mw_pulses']['detuning'] * 1e-6,
                    self.settings['mw_pulses']['pi_half_pulse_time'], self.settings['coil_voltage']))
            axislist[0].legend(labels=('Ref', 'Ramsey data'), fontsize=8)

        else:

            super(Ramsey_CoilVol, self)._plot(axislist)
            axislist[0].set_title(
                '(initial plot) Keysight N9310A, pulse carved by ' + self.settings['mw_pulses'][
                    'microwave_channel'] + '\nRamsey mw-power:{:0.2f}dBm, mw_freq:{:0.4f} GHz \nafter {:d} averages\n detuning: {:0.4f} MHz, π/2 time: {:0.2f} ns\n coil voltage = {:0.4f} V'.format(
                    self.settings['mw_pulses']['mw_power'],
                    self.settings['mw_pulses']['resonant_freq'] * 1e-9, avg_block_number, self.settings['mw_pulses']['detuning'] * 1e-6,
                    self.settings['mw_pulses']['pi_half_pulse_time'], self.settings['coil_voltage']))
            axislist[0].legend(labels=('Ref', 'Ramsey data'), fontsize=8)

    def _update_plot(self, axislist):
            # self._plot(axislist)
            if len(axislist[0].lines) == 0:
                self._plot(axislist)
                return
            super(Ramsey_CoilVol, self)._update_plot(axislist)

            axislist[0].set_title(
                '(updating plot) Keysight N9310A, pulse carved by ' + self.settings['mw_pulses'][
                    'microwave_channel'] + '\nRamsey mw-power:{:0.2f}dBm, mw_freq:{:0.4f} GHz \n detuning: {:0.4f} MHz, π/2 time: {:0.2f} ns\n coil voltage = {:0.4f} V'.format(
                    self.settings['mw_pulses']['mw_power'],
                    self.settings['mw_pulses']['resonant_freq'] * 1e-9, self.settings['mw_pulses']['detuning'] * 1e-6,
                    self.settings['mw_pulses']['pi_half_pulse_time'], self.settings['coil_voltage']))
            axislist[0].legend(labels=('Ref', 'Ramsey data'), fontsize=8)

class Ramsey_RnSIQ(PulsedExperimentBaseScript):
    """
        This script measures the FID signal of NV with Ramsey sequence to measure its dephasing time.
        The script uses Rohde&Schwarz signal generator SMB100A.
        The phase is controlled by an IQ mixer.

        ==> Last edited by Ziwei 5/11/2019

    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_power', -45.0, float, 'microwave power in dBm'),
            Parameter('resonant_freq', 2.87e9, float, 'resonant frequency in Hz'),
            Parameter('detuning', -50e6, float, 'detuning for Ramsey experiment in Hz]'),
            Parameter('pi_half_pulse_time', 50.0, float, 'time duration of a pi/2 pulse (in ns)'),
            Parameter('three_pi_half_pulse_time', 135.0, float, 'time duration of a 3pi/2 pulse (in ns)'),
            Parameter('microwave_channel', 'i', ['i', 'q'], 'channel to use for mw pulses')
        ]),
        Parameter('tau_times', [
            Parameter('min_time', 20, float, 'minimum time between the two pi-half pulses (in ns)'),
            Parameter('max_time', 2000, float, 'maximum time between the two pi-half pulses (in ns)'),
            Parameter('time_step', 20., [2.5, 5., 10., 20., 50., 100., 200., 500., 1000., 10000., 100000., 500000.],
                  'time step increment of rabi pulse duration (in ns)')
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 235, float, 'measurement time after rabi sequence (in ns)'),
            Parameter('nv_reset_time', 2000, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 600, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 550, int, '[ns] delay between laser on and APD readout (given by spontaneous decay rate)')
        ]),
        Parameter('min_pulse_dur', 15, [10, 15, 50], 'Minimum allowed pulse duration (ns) >10 ns'),
        Parameter('num_averages', 500000, int, 'number of averages (>100000)'),
    ]

    _INSTRUMENTS = {'daq': NI6353, 'PB': LISE607RTPulseBlaster, 'mw_gen': R8SMicrowaveGenerator}

    def _function(self):

        freq_to_do = self.settings['mw_pulses']['resonant_freq'] + self.settings['mw_pulses']['detuning']
        ## Turn off green light (the pulse blaster will pulse it on when needed)
        self.instruments['PB']['instance'].update({'laser': {'status': False}})
        self.instruments['PB']['instance'].update({'apd_switch': {'status': False}})

        # Set up the microwave generator SMB100A
        self.instruments['mw_gen']['instance'].update({'freq_mode': 'CW'})
        self.instruments['mw_gen']['instance'].update({'power_mode': 'CW'})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulses']['mw_frequency']})
        self.instruments['mw_gen']['instance'].update({'enable_output': True})
        # This lets the RF go through the IQ mixer
        self.instruments['PB']['instance'].update({'microwave_switch_II': {'status': False}})

        super(Ramsey_RnSIQ, self)._function()

        self.data['exp_finished'] = 1 # this is to flag when the experiment is done or stopped. this allows proper saving even fit fails

        # turn off laser, apd_switch and MW (just in case)
        self.instruments['mw_gen']['instance'].update({'enable_output': False})
        self.instruments['PB']['instance'].update({'laser': {'status': False}})
        self.instruments['PB']['instance'].update({'apd_switch': {'status': False}})

        # start fitting
        # if 'counts' in self.data.keys() and 'tau' in self.data.keys() and not 0 in self.data['counts'][:, 0]:
        #     counts = self.data['counts'][:, 1] / self.data['counts'][:, 0]
        #     tau = self.data['tau']
        #     fit_success = False
            # try:
            #     rabi_fits = fit_rabi_decay(tau, counts, variable_phase=True)
            #     if rabi_fits[1] == True:
            #         self.data['fits'] = rabi_fits[0]
            #         fit_success = True
            #     else:
            #         # self.data['fits'] = None
            #         fit_success = False
            #         print('rabi fit failed')
            #         self.log('rabi fit failed')
            # except:
            #     # self.data['fits'] = None
            #     fit_success = False
            #     print('rabi fit failed')
            #     self.log('rabi fit failed')
            #
            # if fit_success:
            #     counts = self.data['counts'][:,1]/ self.data['counts'][:,0]
            #     tau =self.data['tau']
            #     fits = self.data['fits']
            #
            #     RabiT = 2*np.pi / fits[1]
            #     # RabiT = (2 * np.pi / fits[1]) # avoid negative rabi period ZQ 8/25/2017
            #     T2_star = fits[4]
            #     phaseoffs = fits[2]
            #     pi_time = RabiT / 2 - phaseoffs * RabiT / (2 * np.pi)
            #     pi_half_time = RabiT / 4 - phaseoffs * RabiT / (2 * np.pi)
            #     three_pi_half_time = 3 * RabiT / 4 - phaseoffs * RabiT / (2 * np.pi)
            #     # pi_time = RabiT / 2
            #     # pi_half_time = RabiT / 4
            #     # three_pi_half_time = 3 * RabiT / 4
            #     self.data['phaseoffs'] = phaseoffs * RabiT / (2 * np.pi)
            #     self.data['pi_time'] = pi_time
            #     self.data['pi_half_time'] = pi_half_time
            #     self.data['three_pi_half_time'] = three_pi_half_time
            #     self.data['T2_star'] = T2_star
            #     self.data['rabi_freq'] = 1000 * fits[1] / (2 * np.pi) # Rabi frequency in [MHz]

    def _create_pulse_sequences(self):
        """

        Returns: pulse_sequences, num_averages, tau_list, meas_time
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
            meas_time: the width (in ns) of the daq measurement

        """
        MIN_DURATION = self.settings['min_pulse_dur']
        pulse_sequences = []
        max_range = int(np.floor((self.settings['tau_times']['max_time']-self.settings['tau_times']['min_time'])/self.settings['tau_times']['time_step']))
        tau_list = np.array([self.settings['tau_times']['min_time'] + i*self.settings['tau_times']['time_step'] for i in range(max_range)])

        # ignore the sequence if the mw-pulse is shorter than MIN_DURATION (0 is ok because there is no mw pulse!)
        tau_list = [x for x in tau_list if x == 0 or x >= MIN_DURATION]
        print('tau_list in [ns]:', tau_list)

        pi_half_pulse_time = self.settings['mw_pulses']['pi_half_pulse_time']
        three_pi_half_pulse_time = self.settings['mw_pulses']['three_pi_half_pulse_time']
        meas_time = self.settings['read_out']['meas_time']
        nv_reset_time = self.settings['read_out']['nv_reset_time']
        laser_off_time = self.settings['read_out']['laser_off_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        delay_readout = self.settings['read_out']['delay_readout']

        microwave_channel = 'microwave_' + self.settings['mw_pulses']['microwave_channel'] + '_II'

        for tau in tau_list:
            # Normal Ramsey sequence
            pulse_sequence = [Pulse('laser', laser_off_time + pi_half_pulse_time + tau + pi_half_pulse_time, nv_reset_time),
                              Pulse('apd_readout', laser_off_time + pi_half_pulse_time + tau + pi_half_pulse_time + delay_readout, meas_time),
                              ]

            end_of_first_laser = laser_off_time + pi_half_pulse_time + tau + pi_half_pulse_time + nv_reset_time

            pulse_sequence.append(Pulse(microwave_channel, end_of_first_laser + laser_off_time, pi_half_pulse_time))
            pulse_sequence.append(Pulse(microwave_channel,
                                        end_of_first_laser + laser_off_time + pi_half_pulse_time + tau, pi_half_pulse_time))
            pulse_sequence.append(Pulse('laser',
                                        end_of_first_laser + laser_off_time + pi_half_pulse_time + tau + pi_half_pulse_time + delay_mw_readout,
                                        nv_reset_time))
            pulse_sequence.append(Pulse('apd_readout',
                                        end_of_first_laser + laser_off_time + pi_half_pulse_time + tau + pi_half_pulse_time + delay_mw_readout + delay_readout,
                                        meas_time))

            pulse_sequence.append(
                Pulse('apd_switch', 0, end_of_first_laser + laser_off_time + pi_half_pulse_time + tau + pi_half_pulse_time + delay_mw_readout + nv_reset_time))

            # The sequence implementing 3pi/2 readout...
            # start_of_first_Ramsey = laser_off_time
            # # the first Ramsey with pi/2 readout
            # pulse_sequence = [Pulse('microwave_switch', start_of_first_Ramsey, pi_half_pulse_time)]
            # pulse_sequence.append(Pulse('microwave_switch', start_of_first_Ramsey + pi_half_pulse_time + tau, pi_half_pulse_time))
            # pulse_sequence.append(Pulse('laser',
            #                             start_of_first_Ramsey + pi_half_pulse_time + tau + pi_half_pulse_time + delay_mw_readout,
            #                             nv_reset_time))
            # pulse_sequence.append(Pulse('apd_readout',
            #                             start_of_first_Ramsey + pi_half_pulse_time + tau + pi_half_pulse_time + delay_mw_readout + delay_readout,
            #                             meas_time))
            #
            # start_of_second_Ramsey = start_of_first_Ramsey + pi_half_pulse_time + tau + pi_half_pulse_time + delay_mw_readout + nv_reset_time + laser_off_time
            # # the second Ramsey with 3pi/2 readout
            # pulse_sequence.append(Pulse('microwave_switch', start_of_second_Ramsey, pi_half_pulse_time))
            # pulse_sequence.append(
            #     Pulse('microwave_switch', start_of_second_Ramsey + pi_half_pulse_time + tau, three_pi_half_pulse_time))
            # pulse_sequence.append(Pulse('laser',
            #                             start_of_second_Ramsey + pi_half_pulse_time + tau + three_pi_half_pulse_time + delay_mw_readout,
            #                             nv_reset_time))
            # pulse_sequence.append(Pulse('apd_readout',
            #                             start_of_second_Ramsey + pi_half_pulse_time + tau + three_pi_half_pulse_time + delay_mw_readout + delay_readout,
            #                             meas_time))
            # pulse_sequence.append(
            #     Pulse('apd_switch', 0, start_of_second_Ramsey + pi_half_pulse_time + tau + three_pi_half_pulse_time + delay_mw_readout + nv_reset_time))


            # pulse_sequence = [
            #     Pulse('laser', laser_off_time + pi_half_pulse_time + tau + pi_half_pulse_time, nv_reset_time),
            #     Pulse('apd_readout', laser_off_time + pi_half_pulse_time + tau + pi_half_pulse_time + delay_readout,
            #           meas_time)]
            #
            # pulse_sequence.append(Pulse('microwave_switch',
            #                             laser_off_time + pi_half_pulse_time + tau + pi_half_pulse_time + nv_reset_time + laser_off_time,
            #                             pi_half_pulse_time))
            #
            # end_of_first_MW = laser_off_time + pi_half_pulse_time + tau + pi_half_pulse_time + nv_reset_time + laser_off_time + pi_half_pulse_time
            #
            # pulse_sequence.append(Pulse('microwave_switch', end_of_first_MW + tau, pi_half_pulse_time))
            #
            # end_of_second_MW = end_of_first_MW + tau + pi_half_pulse_time
            #
            # pulse_sequence.append(Pulse('laser',
            #                             end_of_second_MW + tau + pi_half_pulse_time + delay_mw_readout, nv_reset_time))
            #
            # pulse_sequence.append(
            #     Pulse('apd_readout', end_of_second_MW + tau + pi_half_pulse_time + delay_mw_readout + delay_readout,
            #           meas_time))
            # pulse_sequence.append(
            #     Pulse('apd_switch', 0, end_of_second_MW + tau + pi_half_pulse_time + delay_mw_readout + nv_reset_time))

            # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
            # if tau == 0 or tau>=15:
            pulse_sequences.append(pulse_sequence)

        print('number of sequences before validation: ', len(pulse_sequences))
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

        # if 'fits' in data.keys() and data['fits'] is not None:
        #     counts = data['counts'][:,1]/ data['counts'][:,0]
        #     mean_fluor = np.mean(data['counts'][:,0])
        #     tau = data['tau']
        #     fits = data['fits'] # amplitude, frequency, phase, offset
        #
        #     # axislist[0].plot(tau, counts, 'b')
        #     axislist[0].plot(tau, counts, '.-')
        #
        #     #axislist[0].hold(True) ER 20181012
        #
        #     # axislist[0].plot(tau, cose_with_decay(tau, *fits), 'k', lw=3)
        #     axislist[0].plot(tau, cose_with_decay(tau, *fits), lw=3)
        #
        #
        #     # RabiT = 2 * np.pi / fits[1]
        #     T2_star = self.data['T2_star']
        #     phaseoffs = self.data['phaseoffs']
        #     pi_time = self.data['pi_time']
        #     pi_half_time = self.data['pi_half_time']
        #     three_pi_half_time = self.data['three_pi_half_time']
        #     # rabi_freq = 1000 * fits[1] / (2 * np.pi)
        #     rabi_freq = self.data['rabi_freq'] # Rabi frequency in [MHz]
        #
        #
        #     axislist[0].plot(pi_time, cose_with_decay(pi_time, *fits), 'ro', lw=3)
        #     axislist[0].annotate('$\pi$={:0.1f} ns'.format(pi_time), xy=(pi_time, cose_with_decay(pi_time, *fits)),
        #                          xytext=(pi_time - 10., cose_with_decay(pi_time, *fits) - .01), xycoords='data')
        #     axislist[0].plot(pi_half_time, cose_with_decay(pi_half_time, *fits), 'ro', lw=3)
        #     axislist[0].annotate('$\pi/2$={:0.1f} ns'.format(pi_half_time),
        #                          xy=(pi_half_time, cose_with_decay(pi_half_time, *fits)),
        #                          xytext=(pi_half_time + 10., cose_with_decay(pi_half_time, *fits)), xycoords='data')
        #     axislist[0].plot(three_pi_half_time, cose_with_decay(three_pi_half_time, *fits), 'ro', lw=3)
        #     axislist[0].annotate('$3\pi/2$={:0.1f} ns'.format(three_pi_half_time),
        #                          xy=(three_pi_half_time, cose_with_decay(three_pi_half_time, *fits)),
        #                          xytext=(three_pi_half_time + 10., cose_with_decay(three_pi_half_time, *fits)),
        #                          xycoords='data')
        #     axislist[0].plot(-phaseoffs, cose_with_decay(-phaseoffs, *fits), 'gd', lw=3)
        #     axislist[0].annotate('$start$={:0.1f} ns'.format(-phaseoffs),
        #                          xy=(-phaseoffs, cose_with_decay(-phaseoffs, *fits)),
        #                          xytext=(-phaseoffs + 10., cose_with_decay(-phaseoffs, *fits)), xycoords='data')
        #
        #     #pi_time = 2*np.pi / fits[1] / 2
        #     # pi_time = (np.pi - fits[2])/fits[1]
        #     # pi_half_time = (np.pi/2 - fits[2])/fits[1]
        #     # three_pi_half_time = (3*np.pi/2 - fits[2])/fits[1]
        #
        #  #   axislist[0].set_title('Rabi mw-power:{:0.1f}dBm, mw_freq:{:0.3f} GHz, pi-time: {:2.1f}ns, pi-half-time: {:2.1f}ns, 3pi_half_time: {:2.1f}ns, Rabi freq: {2.1f}MHz'.format(self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency']*1e-9, pi_time, pi_half_time, three_pi_half_time, rabi_freq))
        #     axislist[0].set_title(
        #         'Rabi: {:0.2f}MHz, pi-half time: {:2.1f}ns, pi-time: {:2.1f}ns, 3pi-half time: {:2.1f}ns \n T2*: {:2.1f}ns, Ref fluor: {:0.1f}kcps, mw-power: {:0.1f}dBm, mw_freq: {:0.3f} GHz'.format(
        #             rabi_freq, pi_half_time, pi_time, three_pi_half_time, T2_star, mean_fluor,
        #             self.settings['mw_pulses']['mw_power'],
        #             self.settings['mw_pulses'][
        #                 'mw_frequency'] * 1e-9))
        #     axislist[0].set_xlabel('Rabi tau [ns]')
        #     axislist[0].set_ylabel('Normalized Fluorescence')

        # elif 'exp_finished' in data.keys() is not None and data['exp_finished'] == 1:
        if 'exp_finished' in data.keys() is not None and data['exp_finished'] == 1:

            if 0 not in self.data['counts'][:, 0]:
                axislist[0].plot(data['tau'], data['counts'][:, 1] / data['counts'][:, 0])
                axislist[0].set_xlabel('Ramsey tau [ns]')
                axislist[0].set_ylabel('Normalized Fluorescence')
                axislist[0].legend(labels=('Ramsey data'), fontsize=8)
            else:
                axislist[0].plot(data['tau'], data['counts'][:, 0])
                axislist[0].plot(data['tau'], data['counts'][:, 1])
                axislist[0].set_xlabel('Ramsey tau [ns]')
                axislist[0].set_ylabel('Fluorescence [kcps]')
                axislist[0].legend(labels=('Ref', 'Ramsey data'), fontsize=8)


            axislist[0].set_title(
                '(final plot) R&S SMB100A, pulse carved by ' + self.settings['mw_pulses'][
                    'microwave_channel'] + '\nRamsey mw-power:{:0.2f}dBm, mw_freq:{:0.4f} GHz \n detuning: {:0.4f} MHz, π/2 time: {:0.2f} ns'.format(
                    self.settings['mw_pulses']['mw_power'],
                    self.settings['mw_pulses']['resonant_freq'] * 1e-9, self.settings['mw_pulses']['detuning'] * 1e-6,
                    self.settings['mw_pulses']['pi_half_pulse_time']))
            axislist[0].legend(labels=('Ref', 'Ramsey data'), fontsize=8)

        else:

            super(Ramsey_RnSIQ, self)._plot(axislist)
            axislist[0].set_title(
                '(initial plot) R&S SMB100A, pulse carved by ' + self.settings['mw_pulses'][
                    'microwave_channel'] + '\nRamsey mw-power:{:0.2f}dBm, mw_freq:{:0.4f} GHz \n detuning: {:0.4f} MHz, π/2 time: {:0.2f} ns'.format(
                    self.settings['mw_pulses']['mw_power'],
                    self.settings['mw_pulses']['resonant_freq'] * 1e-9, self.settings['mw_pulses']['detuning'] * 1e-6,
                    self.settings['mw_pulses']['pi_half_pulse_time']))
            axislist[0].legend(labels=('Ref', 'Ramsey data'), fontsize=8)

    def _update_plot(self, axislist):
            # self._plot(axislist)
            if len(axislist[0].lines) == 0:
                self._plot(axislist)
                return
            super(Ramsey_RnSIQ, self)._update_plot(axislist)

            axislist[0].set_title(
                '(updating plot) R&S SMB100A, pulse carved by ' + self.settings['mw_pulses'][
                    'microwave_channel'] + '\nRamsey mw-power:{:0.2f}dBm, mw_freq:{:0.4f} GHz \n detuning: {:0.4f} MHz, π/2 time: {:0.2f} ns'.format(
                    self.settings['mw_pulses']['mw_power'],
                    self.settings['mw_pulses']['resonant_freq'] * 1e-9, self.settings['mw_pulses']['detuning'] * 1e-6,
                    self.settings['mw_pulses']['pi_half_pulse_time']))
            axislist[0].legend(labels=('Ref', 'Ramsey data'), fontsize=8)

class TwoPhotonRabi(PulsedExperimentBaseScript):
    """
        This script measures the two-photon rabi oscillation between |+1> and |-1> states.
        Keysight N9310A I is used to apply the pi pulse.
        Rhode&Schwartz RF generator (SMB100A) is used to drive the upper transition.
        Keysight N9310A II is used to drive the lower transition.

        ==> Last edited by ZQ 4/19/2019 5:42pm

    """
    _DEFAULT_SETTINGS = [

        Parameter('mw_pulses', [
            Parameter('mw_power', -45.0, float, 'microwave power in dBm, used to apply pi pulses'),
            Parameter('mw_frequency', 'resonance_frequency_I', ['resonance_frequency_I', 'resonance_frequency_II'],
                      'choose which resonance frequency to use'),
            Parameter('resonance_frequency_I', 2.82e9, float,
                      'transition frequency between |0> and |-1> in Hz'),
            Parameter('resonance_frequency_II', 2.92e9, float,
                      'transition frequency between |0> and |+1> in Hz'),
            Parameter('pi_pulse_time', 50.0, float, 'time duration of a pi pulse (in ns)'),
            Parameter('microwave_channel', 'i', ['i', 'q'], 'Channel to use for mw pulses'),
            Parameter('mw_power_I', -45.0, float, 'microwave power in dBm, used to induce two-photon rabi'),
            Parameter('mw_power_II', -45.0, float, 'microwave power in dBm, used to induce two-photon rabi'),
            Parameter('detuning_I', 1.0e6, float,
                      'detuning from resonance_frequency_I in Hz'),
            Parameter('detuning_II', 1.0e6, float,
                      'detuning from resonance_frequency_II in Hz'),
        ]),

        Parameter('tau_times', [
            Parameter('min_time', 20, float, 'minimum time for rabi oscillations (in ns)'),
            Parameter('max_time', 200, float, 'total time of rabi oscillations (in ns)'),
            Parameter('time_step', 20., [2.5, 5., 10., 20., 50., 100., 200., 500., 1000., 10000., 100000., 500000.],
                  'time step increment of rabi pulse duration (in ns)')
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 235, float, 'measurement time after rabi sequence (in ns)'),
            Parameter('nv_reset_time', 2000, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 600, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 550, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('min_pulse_dur', 15, [2, 5, 10, 15, 50], 'Minimum allowed pulse duration (ns)'),
        Parameter('num_averages', 500000, int, 'number of averages (>100000)'),
    ]

    _INSTRUMENTS = {'daq': NI6353, 'PB': LISE607RTPulseBlaster, 'mw_gen_a': AgilentMicrowaveGenerator, 'mw_gen': R8SMicrowaveGenerator,
                    'mw_gen_b': AgilentMicrowaveGeneratorII}

    def _function(self):


        self.resonance_frequency_I = self.settings['mw_pulses']['resonance_frequency_I']
        self.resonance_frequency_II = self.settings['mw_pulses']['resonance_frequency_II']
        self.detuned_frequency_I = self.settings['mw_pulses']['resonance_frequency_I'] + \
                                           self.settings['mw_pulses']['detuning_I']
        self.detuned_frequency_II = self.settings['mw_pulses']['resonance_frequency_II'] + \
                                            self.settings['mw_pulses']['detuning_II']
        if self.settings['mw_pulses']['mw_frequency'] == 'resonance_frequency_I':
            self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulses']['resonance_frequency_I']})
            self.resonance_frequency_used = self.settings['mw_pulses']['resonance_frequency_I']
        else:
            self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulses']['resonance_frequency_II']})
            self.resonance_frequency_used = self.settings['mw_pulses']['resonance_frequency_II']

        ## Turn off green light (the pulse blaster will pulse it on when needed)
        self.instruments['PB']['instance'].update({'laser': {'status': False}})
        self.instruments['PB']['instance'].update({'apd_switch': {'status': False}})

        ## turn on Keysight N9310A I
        self.instruments['mw_gen_a']['instance'].update({'enable_IQ': True})
        self.instruments['mw_gen_a']['instance'].update({'enable_modulation': True})
        self.instruments['mw_gen_a']['instance'].update({'freq_mode': 'CW'})
        self.instruments['mw_gen_a']['instance'].update({'power_mode': 'CW'})
        self.instruments['mw_gen_a']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power']})
        self.instruments['mw_gen_a']['instance'].update({'frequency': self.resonance_frequency_used})

        ## turn on Rhode & Schartz microwave generator and Agilent RF generator II
        # self.instruments['mw_gen_b']['instance'].update({'enable_IQ': False})
        # self.instruments['mw_gen_b']['instance'].update({'enable_modulation': False})
        self.instruments['mw_gen_b']['instance'].update({'freq_mode': 'CW'})
        self.instruments['mw_gen_b']['instance'].update({'power_mode': 'CW'})
        self.instruments['mw_gen_b']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power_I']})
        self.instruments['mw_gen_b']['instance'].update({'frequency': self.detuned_frequency_I})

        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power_II']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.detuned_frequency_II})

        self.instruments['mw_gen_a']['instance'].update({'enable_output': True})
        self.instruments['mw_gen']['instance'].update({'enable_output': True})
        self.instruments['mw_gen_b']['instance'].update({'enable_output': True})

        super(TwoPhotonRabi, self)._function() # self.data is emptied...

        self.data['exp_finished'] = 1 # this is to flag when the experiment is done or stopped and this allows proper saving even fit fails

        self.data['resonance_frequency_I'] = self.resonance_frequency_I
        self.data['resonance_frequency_II'] = self.resonance_frequency_II
        self.data['detuned_frequency_I'] = self.detuned_frequency_I
        self.data['detuned_frequency_II'] = self.detuned_frequency_II
        self.data['resonance_frequency_used'] = self.resonance_frequency_used

        # turn off laser, apd_switch and MW (just in case)

        self.instruments['mw_gen_a']['instance'].update({'enable_IQ': False})
        self.instruments['mw_gen_a']['instance'].update({'enable_modulation': False})
        self.instruments['mw_gen']['instance'].update({'enable_output': False})
        self.instruments['mw_gen_b']['instance'].update({'enable_output': False})
        self.instruments['mw_gen_a']['instance'].update({'enable_output': False})
        self.instruments['PB']['instance'].update({'laser': {'status': False}})
        self.instruments['PB']['instance'].update({'apd_switch': {'status': False}})

        # # start fitting
        #
        # if 'counts' in self.data.keys() and 'tau' in self.data.keys():
        #     counts = self.data['counts'][:, 1] / self.data['counts'][:, 0]
        #     tau = self.data['tau']
        #     fit_success = False
        #
        #     # try:
        #     #
        #     #     rabi_fits = fit_rabi_decay(tau, counts, variable_phase=True)
        #     #     if rabi_fits[1] == True:
        #     #         self.data['fits'] = rabi_fits[0]
        #     #         fit_success = True
        #     #     else:
        #     #         # self.data['fits'] = None
        #     #         fit_success = False
        #     #         print('rabi fit failed')
        #     #         self.log('rabi fit failed')
        #     # except:
        #     #     # self.data['fits'] = None
        #     #     fit_success = False
        #     #     print('rabi fit failed')
        #     #     self.log('rabi fit failed')
        # if fit_success:
        #     counts = self.data['counts'][:,1]/ self.data['counts'][:,0]
        #     tau =self.data['tau']
        #     fits = self.data['fits']
        #
        #     RabiT = 2*np.pi / fits[1]
        #     # RabiT = (2 * np.pi / fits[1]) # avoid negative rabi period ZQ 8/25/2017
        #     T2_star = fits[4]
        #     phaseoffs = fits[2]
        #     pi_time = RabiT / 2 - phaseoffs * RabiT / (2 * np.pi)
        #     pi_half_time = RabiT / 4 - phaseoffs * RabiT / (2 * np.pi)
        #     three_pi_half_time = 3 * RabiT / 4 - phaseoffs * RabiT / (2 * np.pi)
        #     # pi_time = RabiT / 2
        #     # pi_half_time = RabiT / 4
        #     # three_pi_half_time = 3 * RabiT / 4
        #     self.data['phaseoffs'] = phaseoffs * RabiT / (2 * np.pi)
        #     self.data['pi_time'] = pi_time
        #     self.data['pi_half_time'] = pi_half_time
        #     self.data['three_pi_half_time'] = three_pi_half_time
        #     self.data['T2_star'] = T2_star
        #     self.data['rabi_freq'] = 1000 * fits[1] / (2 * np.pi) # Rabi frequency in [MHz]

    def _create_pulse_sequences(self):
        """

        Returns: pulse_sequences, num_averages, tau_list, meas_time
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
            meas_time: the width (in ns) of the daq measurement

        """
        MIN_DURATION = self.settings['min_pulse_dur']
        pulse_sequences = []
        # tau_list = range(int(max(15, self.settings['tau_times']['time_step'])), int(self.settings['tau_times']['max_time'] + 15),
        #                  self.settings['tau_times']['time_step'])
        # JG 16-08-25 changed (15ns min spacing is taken care of later):
       # tau_list = list(range(int(self.settings['tau_times']['min_time']),
                            #  int(self.settings['tau_times']['max_time']),
                            #  self.settings['tau_times']['time_step']))

        max_range = int(np.floor((self.settings['tau_times']['max_time']-self.settings['tau_times']['min_time'])/self.settings['tau_times']['time_step']))
        tau_list = np.array([self.settings['tau_times']['min_time'] + i*self.settings['tau_times']['time_step'] for i in range(max_range)])

        # ignore the sequence if the mw-pulse is shorter than MIN_DURATION (0 is ok because there is no mw pulse!)
        tau_list = [x for x in tau_list if x == 0 or x >= MIN_DURATION]
        print('tau_list in [ns]:', tau_list)

        meas_time = self.settings['read_out']['meas_time']
        nv_reset_time = self.settings['read_out']['nv_reset_time']
        laser_off_time = self.settings['read_out']['laser_off_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        delay_readout = self.settings['read_out']['delay_readout']
        pi_time = self.settings['mw_pulses']['pi_pulse_time']
        microwave_channel = 'microwave_' + self.settings['mw_pulses']['microwave_channel']

        for tau in tau_list:
            pulse_sequence = [Pulse('laser', laser_off_time, nv_reset_time),
                              Pulse('apd_readout', laser_off_time + delay_readout, meas_time),
                              ]

            end_of_first_laser = laser_off_time + nv_reset_time

            pulse_sequence.append(Pulse(microwave_channel, end_of_first_laser + laser_off_time, pi_time))
            pulse_sequence.append(Pulse(microwave_channel,
                          end_of_first_laser + laser_off_time + pi_time + tau , pi_time))
            pulse_sequence.append(Pulse('laser',
                                    end_of_first_laser + laser_off_time + pi_time + tau +  pi_time + delay_mw_readout,
                                    nv_reset_time))
            pulse_sequence.append(Pulse('apd_readout',
                                    end_of_first_laser + laser_off_time + pi_time + tau + pi_time + delay_mw_readout + delay_readout,
                                    meas_time))

            end_of_second_laser = end_of_first_laser + laser_off_time + pi_time +  tau +  pi_time + delay_mw_readout + nv_reset_time

            pulse_sequence.append(Pulse(microwave_channel, end_of_second_laser + laser_off_time, pi_time))
            pulse_sequence.append(Pulse('laser',
                                    end_of_second_laser + laser_off_time + pi_time +  tau +  pi_time + delay_mw_readout,
                                    nv_reset_time))
            pulse_sequence.append(Pulse('apd_readout',
                                    end_of_second_laser + laser_off_time + pi_time +  tau + pi_time + delay_mw_readout + delay_readout,
                                    meas_time))

            end_of_third_laser = end_of_second_laser + laser_off_time + pi_time +  tau +  pi_time + delay_mw_readout + nv_reset_time

            pulse_sequence.append(Pulse(microwave_channel, end_of_third_laser + laser_off_time, pi_time))
            pulse_sequence.append(Pulse('microwave_switch_II', end_of_third_laser + laser_off_time + pi_time, tau))
            pulse_sequence.append(Pulse('microwave_switch_III', end_of_third_laser + laser_off_time + pi_time, tau))
            pulse_sequence.append(Pulse(microwave_channel,
                          end_of_third_laser + laser_off_time + pi_time  + tau, pi_time))

            pulse_sequence.append(Pulse('laser',
                                    end_of_third_laser + laser_off_time + pi_time +  tau + pi_time + delay_mw_readout,
                                    nv_reset_time))
            pulse_sequence.append(Pulse('apd_readout',
                                    end_of_third_laser + laser_off_time + pi_time + tau + pi_time + delay_mw_readout + delay_readout,
                                    meas_time))

            end_of_fourth_laser = end_of_third_laser + laser_off_time + pi_time + tau + pi_time + delay_mw_readout + nv_reset_time

            pulse_sequence.append(Pulse('apd_switch', 0, end_of_fourth_laser))
            pulse_sequences.append(pulse_sequence)

        # for tau in tau_list:
        #     # print(tau)
        #     pulse_sequence = [Pulse('laser', laser_off_time, nv_reset_time),
        #                       Pulse('apd_readout', laser_off_time + delay_readout, meas_time),
        #                       Pulse('microwave_switch', laser_off_time + nv_reset_time + laser_off_time, meas_time),
        #                       ]
        #
        #     # if tau is 0 there is actually no mw pulse
        #     if tau > 0:
        #         pulse_sequence.append(Pulse('microwave_switch',
        #                                     laser_off_time + tau + 2*mw_sw_buffer + nv_reset_time + laser_off_time,
        #                                     tau))
        #
        #     pulse_sequence.append(Pulse('laser',
        #                                 laser_off_time + tau + 2*mw_sw_buffer + nv_reset_time + laser_off_time + tau + 2*mw_sw_buffer + delay_mw_readout,
        #                                 nv_reset_time))
        #     pulse_sequence.append(Pulse('apd_readout',
        #                                 laser_off_time + tau + 2*mw_sw_buffer + nv_reset_time + laser_off_time + tau + 2*mw_sw_buffer + delay_mw_readout + delay_readout,
        #                                 meas_time))
        #     pulse_sequence.append(Pulse('apd_switch', 0,
        #                                 laser_off_time + tau + 2 * mw_sw_buffer + nv_reset_time + laser_off_time + tau + 2 * mw_sw_buffer + delay_mw_readout + nv_reset_time))
        #
        #     # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
        #     # if tau == 0 or tau>=15:
        #     pulse_sequences.append(pulse_sequence)

        print('number of sequences before validation ', len(pulse_sequences))
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

        if 'fits' in data.keys() and data['fits'] is not None:
            print('fit exists')
            counts = data['counts'][:,1]/ data['counts'][:,0]
            tau = data['tau']
            fits = data['fits'] # amplitude, frequency, phase, offset

            axislist[0].plot(tau, counts, 'b')
            #axislist[0].hold(True) ER 20181012

            axislist[0].plot(tau, cose_with_decay(tau, *fits), 'k', lw=3)

            # RabiT = 2 * np.pi / fits[1]
            T2_star = self.data['T2_star']
            phaseoffs = self.data['phaseoffs']
            pi_time = self.data['pi_time']
            pi_half_time = self.data['pi_half_time']
            three_pi_half_time = self.data['three_pi_half_time']
            # rabi_freq = 1000 * fits[1] / (2 * np.pi)
            rabi_freq = self.data['rabi_freq'] # Rabi frequency in [MHz]


            axislist[0].plot(pi_time, cose_with_decay(pi_time, *fits), 'ro', lw=3)
            axislist[0].annotate('$\pi$={:0.1f}'.format(pi_time), xy=(pi_time, cose_with_decay(pi_time, *fits)),
                                 xytext=(pi_time - 10., cose_with_decay(pi_time, *fits) - .01), xycoords='data')
            axislist[0].plot(pi_half_time, cose_with_decay(pi_half_time, *fits), 'ro', lw=3)
            axislist[0].annotate('$\pi/2$={:0.1f}'.format(pi_half_time),
                                 xy=(pi_half_time, cose_with_decay(pi_half_time, *fits)),
                                 xytext=(pi_half_time + 10., cose_with_decay(pi_half_time, *fits)), xycoords='data')
            axislist[0].plot(three_pi_half_time, cose_with_decay(three_pi_half_time, *fits), 'ro', lw=3)
            axislist[0].annotate('$3\pi/2$={:0.1f}'.format(three_pi_half_time),
                                 xy=(three_pi_half_time, cose_with_decay(three_pi_half_time, *fits)),
                                 xytext=(three_pi_half_time + 10., cose_with_decay(three_pi_half_time, *fits)),
                                 xycoords='data')
            axislist[0].plot(-phaseoffs, cose_with_decay(-phaseoffs, *fits), 'gd', lw=3)
            axislist[0].annotate('$start$={:0.1f}'.format(-phaseoffs),
                                 xy=(-phaseoffs, cose_with_decay(-phaseoffs, *fits)),
                                 xytext=(-phaseoffs + 10., cose_with_decay(-phaseoffs, *fits)), xycoords='data')

            #pi_time = 2*np.pi / fits[1] / 2
            # pi_time = (np.pi - fits[2])/fits[1]
            # pi_half_time = (np.pi/2 - fits[2])/fits[1]
            # three_pi_half_time = (3*np.pi/2 - fits[2])/fits[1]

         #   axislist[0].set_title('Rabi mw-power:{:0.1f}dBm, mw_freq:{:0.3f} GHz, pi-time: {:2.1f}ns, pi-half-time: {:2.1f}ns, 3pi_half_time: {:2.1f}ns, Rabi freq: {2.1f}MHz'.format(self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency']*1e-9, pi_time, pi_half_time, three_pi_half_time, rabi_freq))
            axislist[0].set_title(
                'Rabi: {:0.1f}MHz, pi-half time: {:2.1f}ns, pi-time: {:2.1f}ns, 3pi-half time: {:2.1f}ns'.format(
                    rabi_freq, pi_half_time, pi_time, three_pi_half_time))

        elif 'exp_finished' in data.keys() is not None and data['exp_finished'] == 1:

            axislist[0].plot(data['tau'], data['counts'][:, 0])
            axislist[0].plot(data['tau'], data['counts'][:, 1])
            axislist[0].plot(data['tau'], data['counts'][:, 2])
            axislist[0].plot(data['tau'], data['counts'][:, 3])

            axislist[0].set_xlabel('Rabi tau [ns]')
            axislist[0].set_ylabel('kCounts/sec')

            axislist[0].set_title(
                '(final plot) Two-photon Rabi \n mw_power_I:{:0.2f}dBm, resonance_freq_I:{:0.4f} GHz, detuning_I: {:0.2f} MHz \n mw_power_II:{:0.2f}dBm, resonance_freq_II:{:0.4f} GHz, detuning_II: {:0.2f} MHz\n'.format(
                    self.settings['mw_pulses']['mw_power_I'],
                    self.settings['mw_pulses']['resonance_frequency_I'] * 1e-9,
                    self.settings['mw_pulses']['detuning_I'] * 1e-6,
                    self.settings['mw_pulses']['mw_power_II'],
                    self.settings['mw_pulses']['resonance_frequency_II'] * 1e-9,
                    self.settings['mw_pulses']['detuning_II'] * 1e-6) +
                self.settings['mw_pulses']['mw_frequency'] + ' is used for π pulse')

            axislist[0].legend(
                labels=('Ref |0>', 'Ref with two π pulses', 'Ref with one π pulses', 'Two-photon Rabi Data'),
                fontsize=8)

        else:

            super(TwoPhotonRabi, self)._plot(axislist)
            axislist[0].set_title(
                '(initial plot) Two-photon Rabi \n mw_power_I:{:0.2f}dBm, resonance_freq_I:{:0.4f} GHz, detuning_I: {:0.2f} MHz \n mw_power_II:{:0.2f}dBm, resonance_freq_II:{:0.4f} GHz, detuning_II: {:0.2f} MHz\n'.format(
                    self.settings['mw_pulses']['mw_power_I'],
                    self.settings['mw_pulses']['resonance_frequency_I'] * 1e-9,
                    self.settings['mw_pulses']['detuning_I'] * 1e-6,
                    self.settings['mw_pulses']['mw_power_II'],
                    self.settings['mw_pulses']['resonance_frequency_II'] * 1e-9,
                    self.settings['mw_pulses']['detuning_II'] * 1e-6) +
                self.settings['mw_pulses']['mw_frequency'] + ' is used for π pulse')

            axislist[0].legend(
                labels=('Ref |0>', 'Ref with two π pulses', 'Ref with one π pulses', 'Two-photon Rabi Data'),
                fontsize=8)

    def _update_plot(self, axislist):
            if len(axislist[0].lines) == 0:
                self._plot(axislist)
                return
            super(TwoPhotonRabi, self)._update_plot(axislist)

            axislist[0].set_title(
                '(updating plot) Two-photon Rabi \n mw_power_I:{:0.2f}dBm, resonance_freq_I:{:0.4f} GHz, detuning_I: {:0.2f} MHz \n mw_power_II:{:0.2f}dBm, resonance_freq_II:{:0.4f} GHz, detuning_II: {:0.2f} MHz\n'.format(
                    self.settings['mw_pulses']['mw_power_I'],
                    self.settings['mw_pulses']['resonance_frequency_I'] * 1e-9,
                    self.settings['mw_pulses']['detuning_I'] * 1e-6,
                    self.settings['mw_pulses']['mw_power_II'],
                    self.settings['mw_pulses']['resonance_frequency_II'] * 1e-9,
                    self.settings['mw_pulses']['detuning_II'] * 1e-6) +
                self.settings['mw_pulses']['mw_frequency'] + ' is used for π pulse')
            axislist[0].legend(
                labels=('Ref |0>', 'Ref with two π pulses', 'Ref with one π pulses', 'Two-photon Rabi Data'),
                fontsize=8)


# class Rabi_old_fromB26(PulsedExperimentBaseScript): # ER 5.25.2017
#     """
# This script applies a microwave pulse at fixed power for varying durations to measure Rabi oscillations.
# Uses a double_init scheme
#     """
#     _DEFAULT_SETTINGS = [
#         Parameter('mw_pulses', [
#             Parameter('mw_power', -45.0, float, 'microwave power in dB'),
#             Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
#             Parameter('microwave_channel', 'i', ['i', 'q'], 'Channel to use for mw pulses')
#         ]),
#         Parameter('tau_times', [
#             Parameter('min_time', 15, float, 'minimum time for rabi oscillations (in ns)'),
#             Parameter('max_time', 200, float, 'total time of rabi oscillations (in ns)'),
#             Parameter('time_step', 5., [2.5, 5., 10., 20., 50., 100., 200., 500., 1000., 10000., 100000., 500000.],
#                   'time step increment of rabi pulse duration (in ns)')
#         ]),
#         Parameter('read_out', [
#             Parameter('meas_time', 250, float, 'measurement time after rabi sequence (in ns)'),
#             Parameter('nv_reset_time', 1750, int, 'time with laser on to reset state'),
#             Parameter('laser_off_time', 1000, int,
#                       'minimum laser off time before taking measurements (ns)'),
#             Parameter('delay_mw_readout', 100, int, 'delay between mw and readout (in ns)'),
#             Parameter('delay_readout', 30, int, 'delay between laser on and readout (given by spontaneous decay rate)')
#         ]),
#         Parameter('num_averages', 100000, int, 'number of averages'),
#     ]
#
#     _INSTRUMENTS = {'daq': NI6259, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator}
#
#     def _function(self):
#         #COMMENT_ME
#
#         self.data['fits'] = None
#         self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
#         self.instruments['mw_gen']['instance'].update({'enable_modulation': True})
#         self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power']})
#         self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulses']['mw_frequency']})
#         super(Rabi, self)._function()
#
#         if 'counts' in self.data.keys() and 'tau' in self.data.keys():
#             counts = self.data['counts'][:, 1] / self.data['counts'][:, 0]
#             tau = self.data['tau']
#
#             try:
#                 fits = fit_rabi_decay(tau, counts, variable_phase=True)
#                 self.data['fits'] = fits
#             except:
#                 self.data['fits'] = None
#                 self.log('rabi fit failed')
#
#     def _create_pulse_sequences(self):
#         """
#
#         Returns: pulse_sequences, num_averages, tau_list, meas_time
#             pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
#             scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
#             sequence must have the same number of daq read pulses
#             num_averages: the number of times to repeat each pulse sequence
#             tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
#             meas_time: the width (in ns) of the daq measurement
#
#         """
#         pulse_sequences = []
#         # tau_list = range(int(max(15, self.settings['tau_times']['time_step'])), int(self.settings['tau_times']['max_time'] + 15),
#         #                  self.settings['tau_times']['time_step'])
#         # JG 16-08-25 changed (15ns min spacing is taken care of later):
#        # tau_list = list(range(int(self.settings['tau_times']['min_time']),
#                             #  int(self.settings['tau_times']['max_time']),
#                             #  self.settings['tau_times']['time_step']))
#
#         max_range = int(np.floor((self.settings['tau_times']['max_time']-self.settings['tau_times']['min_time'])/self.settings['tau_times']['time_step']))
#         tau_list = np.array([self.settings['tau_times']['min_time'] + i*self.settings['tau_times']['time_step'] for i in range(max_range)])
#
#         # ignore the sequence if the mw-pulse is shorter than 15ns (0 is ok because there is no mw pulse!)
#         tau_list = [x for x in tau_list if x == 0 or x >= 15]
#
#         nv_reset_time = self.settings['read_out']['nv_reset_time']
#         delay_readout = self.settings['read_out']['delay_readout']
#         microwave_channel = 'microwave_' + self.settings['mw_pulses']['microwave_channel']
#
#         laser_off_time = self.settings['read_out']['laser_off_time']
#         meas_time = self.settings['read_out']['meas_time']
#         delay_mw_readout = self.settings['read_out']['delay_mw_readout']
#
#         for tau in tau_list:
#             pulse_sequence = [Pulse('laser', laser_off_time + tau + 2*40, nv_reset_time),
#                               Pulse('apd_readout', laser_off_time + tau + 2*40 + delay_readout, meas_time)]
#
#             # if tau is 0 there is actually no mw pulse
#             if tau > 0:
#                 pulse_sequence.append(Pulse(microwave_channel,
#                                             laser_off_time + tau + 2*40 + nv_reset_time + laser_off_time,
#                                             tau))
#
#             pulse_sequence.append(Pulse('laser',
#                                         laser_off_time + tau + 2*40 + nv_reset_time + laser_off_time + tau + 2*40 + delay_mw_readout,
#                                         nv_reset_time))
#             pulse_sequence.append(Pulse('apd_readout',
#                                         laser_off_time + tau + 2*40 + nv_reset_time + laser_off_time + tau + 2*40 + delay_mw_readout + delay_readout,
#                                         meas_time))
#             # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
#             # if tau == 0 or tau>=15:
#             pulse_sequences.append(pulse_sequence)
#
#         return pulse_sequences, tau_list, meas_time
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
#         if data is None:
#             data = self.data
#
#         if 'fits' in data.keys() and data['fits'] is not None:
#             counts = data['counts'][:,1]/ data['counts'][:,0]
#             tau = data['tau']
#             fits = data['fits'] # amplitude, frequency, phase, offset
#
#             axislist[0].plot(tau, counts, 'b')
#             #axislist[0].hold(True) ER 20181012
#
#             axislist[0].plot(tau, cose_with_decay(tau, *fits), 'k', lw=3)
#             #pi_time = 2*np.pi / fits[1] / 2
#             pi_time = (np.pi - fits[2])/fits[1]
#             pi_half_time = (np.pi/2 - fits[2])/fits[1]
#             three_pi_half_time = (3*np.pi/2 - fits[2])/fits[1]
#             rabi_freq = 1000*fits[1]/(2*np.pi)
#          #   axislist[0].set_title('Rabi mw-power:{:0.1f}dBm, mw_freq:{:0.3f} GHz, pi-time: {:2.1f}ns, pi-half-time: {:2.1f}ns, 3pi_half_time: {:2.1f}ns, Rabi freq: {2.1f}MHz'.format(self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency']*1e-9, pi_time, pi_half_time, three_pi_half_time, rabi_freq))
#             axislist[0].set_title('Rabi: {:0.1f}MHz, pi-half time: {:2.1f}ns, pi-time: {:2.1f}ns, 3pi-half time: {:2.1f}ns'.format(rabi_freq, pi_half_time, pi_time, three_pi_half_time))
#         else:
#             super(Rabi, self)._plot(axislist)
#             axislist[0].set_title('Rabi mw-power:{:0.1f}dBm, mw_freq:{:0.3f} GHz'.format(self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency']*1e-9))
#             axislist[0].legend(labels=('Ref Fluorescence', 'Rabi Data'), fontsize=8)

# class RabiPowerSweepSingleTau(PulsedExperimentBaseScript):
#     """
# This script applies a microwave pulse at fixed power for varying durations to measure Rabi Oscillations
# todo(emma): (write as a double_init scheme)
#     """
#     _DEFAULT_SETTINGS = [
#         Parameter('min_mw_power', -45.0, float, 'minimum microwave power in dB'),
#         Parameter('max_mw_power', -45.0, float, 'maximum microwave power in dB'),
#         Parameter('mw_power_step', 1.0, float, 'power to step by in dB'),
#         Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
#         Parameter('mw_time', 200, float, 'total time of rabi oscillations (in ns)'),
#         Parameter('meas_time', 300, float, 'measurement time after rabi sequence (in ns)'),
#         Parameter('num_averages', 1000000, int, 'number of averages'),
#         Parameter('reset_time', 10000, int, 'time with laser on at the beginning to reset state'),
#     ]
#
#     _INSTRUMENTS = {'daq': NI6259, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator}
#
#     def _function(self):
#         # COMMENT_ME
#         self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
#         self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_frequency']})
#         mw_power_values = np.arange(self.settings['min_mw_power'],
#                                     self.settings['max_mw_power'] + self.settings['mw_power_step'],
#                                     self.settings['mw_power_step'])
#
#      #   print(mw_power_values)
#         self.data = {'mw_power_values': mw_power_values, 'counts_for_mw': np.zeros(len(mw_power_values))}
#         for index, power in enumerate(mw_power_values):
#             self.instruments['mw_gen']['instance'].update({'amplitude': float(power)})
#             super(RabiPowerSweepSingleTau, self)._function(self.data)
#             self.data['counts_for_mw'][index] = self.data['counts'][0]
#
#     def _create_pulse_sequences(self):
#         '''
#
#         Returns: pulse_sequences, num_averages, tau_list
#             pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
#             scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
#             sequence must have the same number of daq read pulses
#             num_averages: the number of times to repeat each pulse sequence
#             tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
#             meas_time: the width (in ns) of the daq measurement
#
#         '''
#         pulse_sequences = []
#         reset_time = self.settings['reset_time']
#         mw_time = self.settings['mw_time']
#         pulse_sequences.append([Pulse('laser', 0, reset_time),
#                                 Pulse('microwave_i', reset_time + 200, mw_time),
#                                 Pulse('laser', reset_time + mw_time + 300, self.settings['meas_time']),
#                                 Pulse('apd_readout', reset_time + mw_time + 300, self.settings['meas_time'])
#                                 ])
#
#         end_time_max = 0
#         for pulse_sequence in pulse_sequences:
#             for pulse in pulse_sequence:
#                 end_time_max = max(end_time_max, pulse.start_time + pulse.duration)
#         for pulse_sequence in pulse_sequences:
#             pulse_sequence.append(
#                 Pulse('laser', end_time_max + 1850, 15))  # Jan Feb 1st 2017: what is 1850??? Need to comment!
#
#         return pulse_sequences, [mw_time], self.settings['meas_time']
#
#     def _plot(self, axes_list, data=None):
#         '''
#         Plot 1: self.data['tau'], the list of times specified for a given experiment, verses self.data['counts'], the data
#         received for each time
#         Plot 2: the pulse sequence performed at the current time (or if plotted statically, the last pulse sequence
#         performed
#
#         Args:
#             axes_list: list of axes to write plots to (uses first 2)
#             data (optional) dataset to plot (dictionary that contains keys counts_for_mw, mw_power_values), if not provided use self.data
#         '''
#         if data is None:
#             data = self.data
#
#         counts = data['counts_for_mw']
#         x_data = data['mw_power_values']
#         axis1 = axes_list[0]
#         if not counts == []:
#             plot_1d_simple_timetrace_ns(axis1, x_data, [counts], x_label='microwave power (dBm)')
#         axis2 = axes_list[1]
#         plot_pulses(axis2, self.pulse_sequences[self.sequence_index])
#
#     def _update_plot(self, axes_list):
#         '''
#         Updates plots specified in _plot above
#         Args:
#             axes_list: list of axes to write plots to (uses first 2)
#
#         '''
#         counts = self.data['counts_for_mw']
#         x_data = self.data['mw_power_values']
#         axis1 = axes_list[0]
#         if not counts == []:
#             update_1d_simple(axis1, x_data, [counts])
#         axis2 = axes_list[1]
#         update_pulse_plot(axis2, self.pulse_sequences[self.sequence_index])