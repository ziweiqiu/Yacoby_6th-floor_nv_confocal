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
from b26_toolkit.instruments import NI6353, LISE607RTPulseBlaster, R8SMicrowaveGenerator, Pulse, AgilentMicrowaveGenerator, AgilentMicrowaveGeneratorII, Agilent33120A
from pylabcontrol.core import Parameter, Script
from pylabcontrol.scripts import SelectPoints
from b26_toolkit.plotting.plots_1d import plot_pulses
from b26_toolkit.data_processing.fit_functions import fit_exp_decay, exp_offset
import random
# from b26_toolkit.scripts import ESR
# from .rabi import Rabi

# MIN_DURATION = 15 # minimum allowed duration in [ns] for pulse blaster pulses

class PDD_N9310A(PulsedExperimentBaseScript):
    """
    This script runs a PDD ( Periodic Dynamical Decoupling) sequence for different number of pi pulses.
    To symmetrize the sequence between the 0 and +/-1 state we reinitialize every time.
    Keysight N9310A generator is used and it has IQ modulation.
    For a single pi-pulse this is a Hahn-echo sequence. For zero pulses this is a Ramsey sequence.
    The sequence is pi/2 - tau/4 - (tau/4 - pi  - tau/4)^n - tau/4 - pi/2
    Tau/2 is the time between the center of the pulses!

    --> Ziwei Qiu 4/17/2019

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
            Parameter('time_step', 100., [2.5, 5., 10., 20., 50., 100., 200., 300., 400., 500., 600., 800., 1000.,2000.,5000., 8000., 10000.,50000., 100000., 500000.],
                      'time step increment of time between pi pulses (in ns)')
        ]),
        Parameter('decoupling_seq', [
            Parameter('type', 'spin_echo', ['spin_echo', 'CPMG', 'XY4', 'XY8'], 'type of dynamical decoupling sequences'),
            Parameter('num_of_pulse_blocks', 1, int, 'number of pulse blocks.'),
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 200, float, 'measurement time after rabi sequence (in ns)'),
            Parameter('nv_reset_time', 2000, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 600, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 550, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('min_pulse_dur', 15, [10, 15, 50], 'Minimum allowed pulse duration (ns)'),
        Parameter('num_averages', 500000, int, 'number of averages (>100000)'),
    ]

    _INSTRUMENTS = {'daq': NI6353, 'PB': LISE607RTPulseBlaster, 'mw_gen_a': AgilentMicrowaveGenerator}

    def _function(self):
        #COMMENT_ME

        self.data['fits'] = None

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

        self.avg_block_number = super(PDD_N9310A, self)._function()

        # this is to flag when the experiment is done or stopped and this allows proper saving even fit fails
        self.data['exp_finished'] = 1
        self.data['avg_block_number'] = self.avg_block_number

        # turn off laser, apd_switch and MW (just in case)
        self.instruments['PB']['instance'].update({'laser': {'status': False}})
        self.instruments['PB']['instance'].update({'apd_switch': {'status': False}})
        self.instruments['mw_gen_a']['instance'].update({'enable_output': False})
        self.instruments['mw_gen_a']['instance'].update({'enable_IQ': False})
        self.instruments['mw_gen_a']['instance'].update({'enable_modulation': False})

        # start fitting
        if 'counts' in self.data.keys() and 'tau' in self.data.keys() and not 0 in (self.data['counts'][:,1] + self.data['counts'][:, 0]):
            self.data['norm_echo'] = 2. * (- self.data['counts'][:, 1] + self.data['counts'][:, 0]) / (
                    self.data['counts'][:, 1] + self.data['counts'][:, 0])
            counts = self.data['norm_echo']

            # counts = (- self.data['counts'][:, 1] + self.data['counts'][:,0]) / (self.data['counts'][:,1] + self.data['counts'][:, 0])
            tau = self.data['tau']
            fit_success = False

            try:
                fits = fit_exp_decay(tau, counts, offset = True, verbose = True)
                self.data['fits'] = fits
                fit_success = True
            except:
                self.data['fits'] = None
                self.log('T2 fit failed')
                fit_success = False

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
            microwave_channel_1 = 'microwave_i'
            microwave_channel_2 = 'microwave_q'
        else:
            microwave_channel_1 = 'microwave_q'
            microwave_channel_2 = 'microwave_i'
        pi_time = self.settings['mw_pulses']['pi_pulse_time']
        pi_half_time = self.settings['mw_pulses']['pi_half_pulse_time']
        three_pi_half_time = self.settings['mw_pulses']['3pi_half_pulse_time']

        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']

        number_of_pulse_blocks = self.settings['decoupling_seq']['num_of_pulse_blocks']

        # for tau in tau_list:
        #     pulse_sequence = \
        #     [
        #         Pulse(microwave_channel, laser_off_time, pi_half_time),
        #         Pulse(microwave_channel, laser_off_time + pi_half_time/2. + tau - pi_time/2., pi_time),
        #         Pulse(microwave_channel, laser_off_time + pi_half_time/2. + tau + tau - pi_half_time/2., pi_half_time)
        #     ]
        #
        #     end_of_first_HE = laser_off_time + pi_half_time/2. + tau + tau - pi_half_time/2. + pi_half_time
        #
        #     pulse_sequence += [
        #          Pulse('laser', end_of_first_HE + delay_mw_readout, nv_reset_time),
        #          Pulse('apd_readout', end_of_first_HE + delay_mw_readout + delay_readout, meas_time),
        #          ]
        #
        #     start_of_second_HE = end_of_first_HE + delay_mw_readout + nv_reset_time + laser_off_time
        #
        #     pulse_sequence += \
        #     [
        #         Pulse(microwave_channel, start_of_second_HE, pi_half_time),
        #         Pulse(microwave_channel, start_of_second_HE + pi_half_time/2. + tau - pi_time/2., pi_time),
        #         Pulse(microwave_channel, start_of_second_HE + pi_half_time/2. + tau + tau - pi_half_time/2., three_pi_half_time)
        #     ]
        #
        #     end_of_second_HE = start_of_second_HE + pi_half_time/2. + tau + tau - pi_half_time/2. + pi_half_time
        #
        #     pulse_sequence += [
        #         Pulse('laser', end_of_second_HE + delay_mw_readout, nv_reset_time),
        #         Pulse('apd_readout', end_of_second_HE + delay_mw_readout + delay_readout, meas_time)
        #     ]
        #
        #     pulse_sequence.append(Pulse('apd_switch', 0, end_of_second_HE + delay_mw_readout + nv_reset_time))
        #     # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
        #     # if tau == 0 or tau>=15:
        #     pulse_sequences.append(pulse_sequence)

        if self.settings['decoupling_seq']['type'] == 'spin_echo':
            for tau_total in tau_list:
                tau = tau_total / (1 * self.settings['decoupling_seq']['num_of_pulse_blocks'])
                pulse_sequence = []

                # ECHO SEQUENCE:
                start_of_first_HE = laser_off_time
                # the first pi/2 pulse

                pulse_sequence = [Pulse(microwave_channel_1, start_of_first_HE, pi_half_time)]
                section_begin_time = start_of_first_HE + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                for i in range(0, number_of_pulse_blocks):
                    # pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                    section_begin_time += 1 * tau

                # the second pi/2 pulse and readout
                pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + tau / 2, pi_half_time))
                pulse_sequence.append(Pulse('laser', section_begin_time + tau / 2 + pi_half_time + delay_mw_readout,
                                             nv_reset_time))
                pulse_sequence.append(Pulse('apd_readout',
                                             section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + delay_readout,
                                             meas_time))
                start_of_second_HE = section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time

                # the first pi/2 pulse
                pulse_sequence.append(Pulse(microwave_channel_1, start_of_second_HE, pi_half_time))
                section_begin_time = start_of_second_HE + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                for i in range(0, number_of_pulse_blocks):
                    # pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                    section_begin_time += 1 * tau

                # the second 3*pi/2 pulse and readout
                pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + tau / 2, three_pi_half_time))
                pulse_sequence.append(Pulse('laser',
                                             section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout,
                                             nv_reset_time))
                pulse_sequence.append(Pulse('apd_readout',
                                             section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + delay_readout,
                                             meas_time))

                # Turn on APD switch all the time
                pulse_sequence.append(Pulse('apd_switch', 0,
                                            section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + nv_reset_time))

                # # DEER SEQUENCE
                # if self.settings['decoupling_seq']['do_deer']:
                #     start_of_DEER = section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time
                #     # the first pi/2 pulse
                #     pulse_sequence.extend([Pulse(microwave_channel_1, start_of_DEER, pi_half_time)])
                #     section_begin_time = start_of_DEER + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                #     for i in range(0, number_of_pulse_blocks):
                #         pulse_sequence.extend(
                #             [Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 1 * tau - RF_pi_time / 2, RF_pi_time)
                #              ])
                #         section_begin_time += 1 * tau
                #
                #     # the second pi/2 pulse and readout
                #     pulse_sequence.extend([Pulse(microwave_channel_1, section_begin_time + tau / 2, pi_half_time),
                #                            Pulse('laser',
                #                                  section_begin_time + tau / 2 + pi_half_time + delay_mw_readout,
                #                                  nv_reset_time),
                #                            Pulse('apd_readout',
                #                                  section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + delay_readout,
                #                                  meas_time)
                #                            ])
                #
                #     start_of_second_DEER = section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time
                #     pulse_sequence.extend([Pulse(microwave_channel_1, start_of_second_DEER, pi_half_time)])
                #     section_begin_time = start_of_second_DEER + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                #
                #     for i in range(0, number_of_pulse_blocks):
                #         pulse_sequence.extend(
                #             [Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 1 * tau - RF_pi_time / 2, RF_pi_time)
                #              ])
                #         section_begin_time += 1 * tau
                #
                #     # the second 3*pi/2 pulse and readout
                #     pulse_sequence.extend([Pulse(microwave_channel_1, section_begin_time + tau / 2, three_pi_half_time),
                #                            Pulse('laser',
                #                                  section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout,
                #                                  nv_reset_time),
                #                            Pulse('apd_readout',
                #                                  section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + delay_readout,
                #                                  meas_time)
                #                            ])

                # print(pulse_sequence)
                pulse_sequences.append(pulse_sequence)

            print('number of sequences before validation ', len(pulse_sequences))
            return pulse_sequences, tau_list, meas_time
            # return pulse_sequences, self.settings['num_averages'], tau_list, meas_time

        elif self.settings['decoupling_seq']['type'] == 'XY4':

            for tau_total in tau_list:
                tau = tau_total / (4 * self.settings['decoupling_seq']['num_of_pulse_blocks'])
                pulse_sequence = []

                # ECHO SEQUENCE:
                start_of_first_HE = laser_off_time
                # the first pi/2 pulse
                pulse_sequence = [Pulse(microwave_channel_1, start_of_first_HE, pi_half_time)]
                section_begin_time = start_of_first_HE + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                for i in range(0, number_of_pulse_blocks):
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 2 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time / 2, pi_time))

                    section_begin_time += 4 * tau

                # the second pi/2 pulse and readout
                pulse_sequence.append(
                    Pulse(microwave_channel_1, section_begin_time + tau / 2, three_pi_half_time))
                pulse_sequence.append(
                    Pulse('laser', section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout,
                          nv_reset_time))
                pulse_sequence.append(
                    Pulse('apd_readout',
                          section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + delay_readout,
                          meas_time))
                start_of_second_HE = section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time

                # the first pi/2 pulse
                pulse_sequence.append(
                    Pulse(microwave_channel_1, start_of_second_HE, pi_half_time))
                section_begin_time = start_of_second_HE + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                for i in range(0, number_of_pulse_blocks):
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 2 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time / 2, pi_time))
                    section_begin_time += 4 * tau

                # the second 3*pi/2 pulse and readout
                pulse_sequence.append(
                    Pulse(microwave_channel_1, section_begin_time + tau / 2, pi_half_time))
                pulse_sequence.append(
                    Pulse('laser',
                          section_begin_time + tau / 2 + pi_half_time + delay_mw_readout,
                          nv_reset_time))
                pulse_sequence.append(
                    Pulse('apd_readout',
                          section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + delay_readout,
                          meas_time))

                # # DEER SEQUENCE
                # if self.settings['decoupling_seq']['do_deer']:
                #     start_of_DEER = section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time
                #     # the first pi/2 pulse
                #     pulse_sequence.extend([Pulse(microwave_channel_1, start_of_DEER, pi_half_time)])
                #     section_begin_time = start_of_DEER + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                #     for i in range(0, number_of_pulse_blocks):
                #         pulse_sequence.extend(
                #             [Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 1 * tau - RF_pi_time / 2, RF_pi_time),
                #              Pulse(microwave_channel_1, section_begin_time + 2 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 2 * tau - RF_pi_time / 2, RF_pi_time),
                #              Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 3 * tau - RF_pi_time / 2, RF_pi_time),
                #              Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 4 * tau - RF_pi_time / 2, RF_pi_time)
                #              ])
                #         section_begin_time += 4 * tau
                #
                #     # the second pi/2 pulse and readout
                #     pulse_sequence.extend([Pulse(microwave_channel_1, section_begin_time + tau / 2, pi_half_time),
                #                            Pulse('laser',
                #                                  section_begin_time + tau / 2 + pi_half_time + delay_mw_readout,
                #                                  nv_reset_time),
                #                            Pulse('apd_readout',
                #                                  section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + delay_readout,
                #                                  meas_time)
                #                            ])
                #
                #     start_of_second_DEER = section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time
                #     pulse_sequence.extend([Pulse(microwave_channel_1, start_of_second_DEER, pi_half_time)])
                #     section_begin_time = start_of_second_DEER + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                #
                #     for i in range(0, number_of_pulse_blocks):
                #         pulse_sequence.extend(
                #             [Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 1 * tau - RF_pi_time / 2, RF_pi_time),
                #              Pulse(microwave_channel_1, section_begin_time + 2 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 2 * tau - RF_pi_time / 2, RF_pi_time),
                #              Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 3 * tau - RF_pi_time / 2, RF_pi_time),
                #              Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 4 * tau - RF_pi_time / 2, RF_pi_time)
                #              ])
                #         section_begin_time += 4 * tau
                #
                #     # the second 3*pi/2 pulse and readout
                #     pulse_sequence.extend([Pulse(microwave_channel_1, section_begin_time + tau / 2, three_pi_half_time),
                #                            Pulse('laser',
                #                                  section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout,
                #                                  nv_reset_time),
                #                            Pulse('apd_readout',
                #                                  section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + delay_readout,
                #                                  meas_time)
                #                            ])
                # Turn on APD switch all the time
                pulse_sequence.append(Pulse('apd_switch', 0,
                                            section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + nv_reset_time))
                pulse_sequences.append(pulse_sequence)

            print('number of sequences before validation ', len(pulse_sequences))
            # return pulse_sequences, self.settings['num_averages'], tau_list, meas_time
            return pulse_sequences, tau_list, meas_time

        elif self.settings['decoupling_seq']['type'] == 'XY8':

            for tau_total in tau_list:
                tau = tau_total / (8 * self.settings['decoupling_seq']['num_of_pulse_blocks'])
                pulse_sequence = []

                # ECHO SEQUENCE:
                start_of_first_HE = laser_off_time

                # the first pi/2 pulse
                pulse_sequence.append(
                    Pulse(microwave_channel_1, start_of_first_HE, pi_half_time))

                section_begin_time = start_of_first_HE + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                for i in range(0, number_of_pulse_blocks):
                    # pulse_sequence.extend(
                    #     # [Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time),
                    #     #  Pulse(microwave_channel_1, section_begin_time + 2 * tau - pi_time / 2, pi_time),
                    #     #  Pulse(microwave_channel_2, section_begin_time + 3 * tau - pi_time / 2, pi_time),
                    #     #  Pulse(microwave_channel_1, section_begin_time + 4 * tau - pi_time / 2, pi_time),
                    #     #  Pulse(microwave_channel_1, section_begin_time + 5 * tau - pi_time / 2, pi_time),
                    #     #  Pulse(microwave_channel_2, section_begin_time + 6 * tau - pi_time / 2, pi_time),
                    #     #  Pulse(microwave_channel_1, section_begin_time + 7 * tau - pi_time / 2, pi_time),
                    #     #  Pulse(microwave_channel_2, section_begin_time + 8 * tau - pi_time / 2, pi_time)
                    #     [Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time),
                    #      Pulse(microwave_channel_2, section_begin_time + 2 * tau - pi_time / 2, pi_time),
                    #      Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time / 2, pi_time),
                    #      Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time / 2, pi_time),
                    #      Pulse(microwave_channel_2, section_begin_time + 5 * tau - pi_time / 2, pi_time),
                    #      Pulse(microwave_channel_1, section_begin_time + 6 * tau - pi_time / 2, pi_time),
                    #      Pulse(microwave_channel_2, section_begin_time + 7 * tau - pi_time / 2, pi_time),
                    #      Pulse(microwave_channel_1, section_begin_time + 8 * tau - pi_time / 2, pi_time)
                    #      ])
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 2 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 5 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 6 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 7 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 8 * tau - pi_time / 2, pi_time))
                    section_begin_time += 8 * tau

                # the second pi/2 pulse and readout
                pulse_sequence.append(
                    Pulse(microwave_channel_1, section_begin_time + tau / 2, three_pi_half_time))
                pulse_sequence.append(
                    Pulse('laser', section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout,
                          nv_reset_time))
                pulse_sequence.append(
                    Pulse('apd_readout',
                          section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + delay_readout,
                          meas_time))

                start_of_second_HE = section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time

                # the first pi/2 pulse
                pulse_sequence.append(
                    Pulse(microwave_channel_1, start_of_second_HE, pi_half_time))
                section_begin_time = start_of_second_HE + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                for i in range(0, number_of_pulse_blocks):
                    # pulse_sequence.extend(
                    #     # [Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time),
                    #     #  Pulse(microwave_channel_1, section_begin_time + 2 * tau - pi_time / 2, pi_time),
                    #     #  Pulse(microwave_channel_2, section_begin_time + 3 * tau - pi_time / 2, pi_time),
                    #     #  Pulse(microwave_channel_1, section_begin_time + 4 * tau - pi_time / 2, pi_time),
                    #     #  Pulse(microwave_channel_1, section_begin_time + 5 * tau - pi_time / 2, pi_time),
                    #     #  Pulse(microwave_channel_2, section_begin_time + 6 * tau - pi_time / 2, pi_time),
                    #     #  Pulse(microwave_channel_1, section_begin_time + 7 * tau - pi_time / 2, pi_time),
                    #     #  Pulse(microwave_channel_2, section_begin_time + 8 * tau - pi_time / 2, pi_time)
                    #     [Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time),
                    #      Pulse(microwave_channel_2, section_begin_time + 2 * tau - pi_time / 2, pi_time),
                    #      Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time / 2, pi_time),
                    #      Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time / 2, pi_time),
                    #      Pulse(microwave_channel_2, section_begin_time + 5 * tau - pi_time / 2, pi_time),
                    #      Pulse(microwave_channel_1, section_begin_time + 6 * tau - pi_time / 2, pi_time),
                    #      Pulse(microwave_channel_2, section_begin_time + 7 * tau - pi_time / 2, pi_time),
                    #      Pulse(microwave_channel_1, section_begin_time + 8 * tau - pi_time / 2, pi_time)
                    #      ])
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 2 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 5 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 6 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 7 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 8 * tau - pi_time / 2, pi_time))
                    section_begin_time += 8 * tau

                # the second 3*pi/2 pulse and readout
                pulse_sequence.append(
                    Pulse(microwave_channel_1, section_begin_time + tau / 2, pi_half_time))
                pulse_sequence.append(
                    Pulse('laser',
                          section_begin_time + tau / 2 + pi_half_time + delay_mw_readout,
                          nv_reset_time))
                pulse_sequence.append(
                    Pulse('apd_readout',
                          section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + delay_readout,
                          meas_time))
                # # DEER SEQUENCE
                # if self.settings['decoupling_seq']['do_deer']:
                #     start_of_DEER = section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time
                #     # the first pi/2 pulse
                #     pulse_sequence.extend([Pulse(microwave_channel_1, start_of_DEER, pi_half_time)])
                #     section_begin_time = start_of_DEER + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                #     for i in range(0, number_of_pulse_blocks):
                #         pulse_sequence.extend(
                #             [Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 1 * tau - RF_pi_time / 2, RF_pi_time),
                #              Pulse(microwave_channel_1, section_begin_time + 2 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 2 * tau - RF_pi_time / 2, RF_pi_time),
                #              Pulse(microwave_channel_2, section_begin_time + 3 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 3 * tau - RF_pi_time / 2, RF_pi_time),
                #              Pulse(microwave_channel_1, section_begin_time + 4 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 4 * tau - RF_pi_time / 2, RF_pi_time),
                #              Pulse(microwave_channel_1, section_begin_time + 5 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 5 * tau - RF_pi_time / 2, RF_pi_time),
                #              Pulse(microwave_channel_2, section_begin_time + 6 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 6 * tau - RF_pi_time / 2, RF_pi_time),
                #              Pulse(microwave_channel_1, section_begin_time + 7 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 7 * tau - RF_pi_time / 2, RF_pi_time),
                #              Pulse(microwave_channel_2, section_begin_time + 8 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 8 * tau - RF_pi_time / 2, RF_pi_time)
                #              ])
                #         section_begin_time += 8 * tau
                #
                #     # the second pi/2 pulse and readout
                #     pulse_sequence.extend([Pulse(microwave_channel_1, section_begin_time + tau / 2, pi_half_time),
                #                            Pulse('laser',
                #                                  section_begin_time + tau / 2 + pi_half_time + delay_mw_readout,
                #                                  nv_reset_time),
                #                            Pulse('apd_readout',
                #                                  section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + delay_readout,
                #                                  meas_time)
                #                            ])
                #
                #     start_of_second_DEER = section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time
                #     pulse_sequence.extend([Pulse(microwave_channel_1, start_of_second_DEER, pi_half_time)])
                #     section_begin_time = start_of_second_DEER + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                #
                #     for i in range(0, number_of_pulse_blocks):
                #         pulse_sequence.extend(
                #             [Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 1 * tau - RF_pi_time / 2, RF_pi_time),
                #              Pulse(microwave_channel_1, section_begin_time + 2 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 2 * tau - RF_pi_time / 2, RF_pi_time),
                #              Pulse(microwave_channel_2, section_begin_time + 3 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 3 * tau - RF_pi_time / 2, RF_pi_time),
                #              Pulse(microwave_channel_1, section_begin_time + 4 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 4 * tau - RF_pi_time / 2, RF_pi_time),
                #              Pulse(microwave_channel_1, section_begin_time + 5 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 5 * tau - RF_pi_time / 2, RF_pi_time),
                #              Pulse(microwave_channel_2, section_begin_time + 6 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 6 * tau - RF_pi_time / 2, RF_pi_time),
                #              Pulse(microwave_channel_1, section_begin_time + 7 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 7 * tau - RF_pi_time / 2, RF_pi_time),
                #              Pulse(microwave_channel_2, section_begin_time + 8 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 8 * tau - RF_pi_time / 2, RF_pi_time)
                #              ])
                #         section_begin_time += 8 * tau
                #
                #     # the second 3*pi/2 pulse and readout
                #     pulse_sequence.extend([Pulse(microwave_channel_1, section_begin_time + tau / 2, three_pi_half_time),
                #                            Pulse('laser',
                #                                  section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout,
                #                                  nv_reset_time),
                #                            Pulse('apd_readout',
                #                                  section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + delay_readout,
                #                                  meas_time)
                #                            ])
                # Turn on APD switch all the time
                pulse_sequence.append(Pulse('apd_switch', 0,
                                            section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + nv_reset_time))
                pulse_sequences.append(pulse_sequence)

            print('number of sequences before validation ', len(pulse_sequences))
            # return pulse_sequences, self.settings['num_averages'], tau_list, meas_time
            return pulse_sequences, tau_list, meas_time

        elif self.settings['decoupling_seq']['type'] == 'CPMG':

            for tau_total in tau_list:
                tau = tau_total / (1 * self.settings['decoupling_seq']['num_of_pulse_blocks'])
                pulse_sequence = []

                # ECHO SEQUENCE:
                start_of_first_HE = laser_off_time

                # the first pi/2 pulse
                pulse_sequence = [Pulse(microwave_channel_1, start_of_first_HE, pi_half_time)]
                section_begin_time = start_of_first_HE + pi_half_time - tau / 2  # for the first pulse, only wait tau/2

                for i in range(0, number_of_pulse_blocks):
                    pulse_sequence.append(Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                    section_begin_time += 1 * tau

                # the second pi/2 pulse and readout
                pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + tau / 2, three_pi_half_time))
                pulse_sequence.append(Pulse('laser', section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout,
                                             nv_reset_time))
                pulse_sequence.append(Pulse('apd_readout',
                                             section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + delay_readout,
                                             meas_time))

                start_of_second_HE = section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time

                # the first pi/2 pulse
                pulse_sequence.append(Pulse(microwave_channel_1, start_of_second_HE, pi_half_time))

                section_begin_time = start_of_second_HE + pi_half_time - tau / 2  # for the first pulse, only wait tau/2

                for i in range(0, number_of_pulse_blocks):
                    pulse_sequence.append(Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                    section_begin_time += 1 * tau

                # the second 3*pi/2 pulse and readout
                pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + tau / 2, pi_half_time))
                pulse_sequence.append(Pulse('laser',
                                             section_begin_time + tau / 2 + pi_half_time + delay_mw_readout,
                                             nv_reset_time))
                pulse_sequence.append(Pulse('apd_readout',
                                             section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + delay_readout,
                                             meas_time))

                # Turn on APD switch all the time
                pulse_sequence.append(Pulse('apd_switch', 0,
                                            section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + nv_reset_time))


    #             # # DEER SEQUENCE
    #             # if self.settings['decoupling_seq']['do_deer']:
    #             #     start_of_DEER = section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time
    #             #     # the first pi/2 pulse
    #             #     pulse_sequence.extend([Pulse(microwave_channel_1, start_of_DEER, pi_half_time)])
    #             #     section_begin_time = start_of_DEER + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
    #             #     for i in range(0, number_of_pulse_blocks):
    #             #         pulse_sequence.extend(
    #             #             [Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time),
    #             #              Pulse('RF_switch', section_begin_time + 1 * tau - RF_pi_time / 2, RF_pi_time)
    #             #              ])
    #             #         section_begin_time += 1 * tau
    #             #
    #             #     # the second pi/2 pulse and readout
    #             #     pulse_sequence.extend([Pulse(microwave_channel_1, section_begin_time + tau / 2, pi_half_time),
    #             #                            Pulse('laser',
    #             #                                  section_begin_time + tau / 2 + pi_half_time + delay_mw_readout,
    #             #                                  nv_reset_time),
    #             #                            Pulse('apd_readout',
    #             #                                  section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + delay_readout,
    #             #                                  meas_time)
    #             #                            ])
    #             #
    #             #     start_of_second_DEER = section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time
    #             #     pulse_sequence.extend([Pulse(microwave_channel_1, start_of_second_DEER, pi_half_time)])
    #             #     section_begin_time = start_of_second_DEER + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
    #             #
    #             #     for i in range(0, number_of_pulse_blocks):
    #             #         pulse_sequence.extend(
    #             #             [Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time),
    #             #              Pulse('RF_switch', section_begin_time + 1 * tau - RF_pi_time / 2, RF_pi_time)
    #             #              ])
    #             #         section_begin_time += 1 * tau
    #             #
    #             #     # the second 3*pi/2 pulse and readout
    #             #     pulse_sequence.extend([Pulse(microwave_channel_1, section_begin_time + tau / 2, three_pi_half_time),
    #             #                            Pulse('laser',
    #             #                                  section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout,
    #             #                                  nv_reset_time),
    #             #                            Pulse('apd_readout',
    #             #                                  section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + delay_readout,
    #             #                                  meas_time)
    #             #                            ])
    #
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

        if 'avg_block_number' in data.keys() is not None and data['avg_block_number'] is not None:
            avg_block_number = data['avg_block_number']
        else:
            avg_block_number = 0

        if 'fits' in data.keys() is not None and data['fits'] is not None:

            # counts = (-data['counts'][:,1] + data['counts'][:,0])/ (data['counts'][:,0] + data['counts'][:,1])
            tau = data['tau']
            fits = data['fits']
            axislist[0].plot(tau, data['norm_echo'])
            tauinterp = np.linspace(np.min(tau), np.max(tau), 100)
            axislist[0].plot(tauinterp, exp_offset(tauinterp, fits[0], fits[1], fits[2]))
            axislist[0].set_title(
                'Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
                    'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s)\n averaged over {:d} blocks \n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz, \n T2 (simple expo, p = 1) = {:2.1f} ns'.format(
                    self.settings['decoupling_seq']['type'],
                    self.settings['decoupling_seq']['num_of_pulse_blocks'], avg_block_number,
                    self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency'] * 1e-9,
                    fits[1]))
            axislist[0].legend(labels=('Echo', 'Exp Fit'), fontsize=10)
            axislist[0].set_ylabel('fluorescence contrast')
            axislist[0].set_xlabel('tau [ns]')

        elif 'exp_finished' in data.keys() is not None and data['exp_finished'] == 1:
            if 'norm_echo' in data.keys() is not None and data['norm_echo'] is not None:
                tau = data['tau']
                axislist[0].plot(tau, data['norm_echo'],'d')
                axislist[0].set_title(
                    '(final plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
                        'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s)\n averaged over {:d} blocks \n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz'.format(
                        self.settings['decoupling_seq']['type'],
                        self.settings['decoupling_seq']['num_of_pulse_blocks'],
                        avg_block_number,
                        self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency'] * 1e-9))
                axislist[0].legend(labels=('Echo'), fontsize=10)
                axislist[0].set_ylabel('fluorescence contrast')
                axislist[0].set_xlabel('tau [ns]')
            else:
                axislist[0].plot(data['tau'], data['counts'][:, 0])
                axislist[0].plot(data['tau'], data['counts'][:, 1])
                axislist[0].set_ylabel('fluorescence [kcps]')
                axislist[0].set_xlabel('tau [ns]')
                axislist[0].legend(labels=('Echo up {:.0f}kcps'.format(np.mean(data['counts'][:, 0])),
                                           'Echo down {:.0f}kcps'.format(np.mean(data['counts'][:, 1]))), fontsize=10)
                axislist[0].set_title(
                    '(final plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
                        'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s)\n averaged over {:d} blocks\n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz'.format(
                        self.settings['decoupling_seq']['type'], self.settings['decoupling_seq']['num_of_pulse_blocks'],
                        avg_block_number, self.settings['mw_pulses']['mw_power'],
                        self.settings['mw_pulses']['mw_frequency'] * 1e-9))

        else:
            super(PDD_N9310A, self)._plot(axislist)

            axislist[0].legend(labels=('Echo up {:.0f}kcps'.format(np.mean(data['counts'][:, 0])),
                                       'Echo down {:.0f}kcps'.format(np.mean(data['counts'][:, 1]))), fontsize=10)
            axislist[0].set_title(
                '(initial plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
                    'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s) \n averaged over {:d} blocks\n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz'.format(
                    self.settings['decoupling_seq']['type'], self.settings['decoupling_seq']['num_of_pulse_blocks'],
                    avg_block_number, self.settings['mw_pulses']['mw_power'],
                    self.settings['mw_pulses']['mw_frequency'] * 1e-9))

    def _update_plot(self, axislist):

        if len(axislist[0].lines) == 0:
            self._plot(axislist)
            return
        # if 'avg_block_number' in self.data.keys() is not None and self.data['avg_block_number'] is not None:
        #     avg_block_number =self.data['avg_block_number']
        # else:
        #     avg_block_number = 0
        super(PDD_N9310A, self)._update_plot(axislist)
        axislist[0].legend(labels=('Echo up {:.0f}kcps'.format(np.mean(self.data['counts'][:, 0])),
                                   'Echo down {:.0f}kcps'.format(np.mean(self.data['counts'][:, 1]))), fontsize=10)
        # axislist[0].set_title('(updating plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
        #     'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s) \n averaged over {:d} blocks \n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz'.format(
        #     self.settings['decoupling_seq']['type'], self.settings['decoupling_seq']['num_of_pulse_blocks'],
        #     avg_block_number, self.settings['mw_pulses']['mw_power'],
        #     self.settings['mw_pulses']['mw_frequency'] * 1e-9))
        axislist[0].set_title('(updating plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
            'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s) \n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz'.format(
            self.settings['decoupling_seq']['type'], self.settings['decoupling_seq']['num_of_pulse_blocks'],
            self.settings['mw_pulses']['mw_power'],
            self.settings['mw_pulses']['mw_frequency'] * 1e-9))

class eSensing_N9310A(PulsedExperimentBaseScript):
    """
    This script runs a PDD ( Periodic Dynamical Decoupling) sequence for different number of pi pulses.
    Applied an AC electric field to be detected, which is controled by the PB14 (aux).
    To symmetrize the sequence between the 0 and +/-1 state we reinitialize every time.
    Keysight N9310A generator is used and it has IQ modulation.
    For a single pi-pulse this is a Hahn-echo sequence. For zero pulses this is a Ramsey sequence.
    The sequence is pi/2 - tau/4 - (tau/4 - pi  - tau/4)^n - tau/4 - pi/2
    Tau/2 is the time between the center of the pulses!

    --> Ziwei Qiu 6/26/2019

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
            Parameter('time_step', 100., [2.5, 5., 10., 20., 50., 100., 200., 300., 400., 500., 600., 800., 1000.,2000.,3000.,4000.,5000., 8000., 10000.,50000., 100000., 500000.],
                      'time step increment of time between pi pulses (in ns)')
        ]),
        Parameter('decoupling_seq', [
            Parameter('type', 'spin_echo', ['spin_echo', 'CPMG', 'XY4', 'XY8'], 'type of dynamical decoupling sequences'),
            Parameter('num_of_pulse_blocks', 1, int, 'number of pulse blocks.'),
            Parameter('do_eSensing', True, bool, 'choose to perform electric sensing half of the sequence or not'),
            Parameter('eSensing_type', 'sine', ['sine', 'cosine'], 'choose to do sine or cosine electrometry'),
            Parameter('e_field_channel', 'aux', ['aux', 'aux2'], 'choose a PB channel for applying electric field. aux - fixed 5V, aux2 - tunable 0-5V')
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 200, float, 'measurement time after rabi sequence (in ns)'),
            Parameter('nv_reset_time', 2000, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 600, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 550, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('min_pulse_dur', 15, [10, 15, 50], 'Minimum allowed pulse duration (ns)'),
        Parameter('num_averages', 500000, int, 'number of averages (>100000)'),
    ]

    _INSTRUMENTS = {'daq': NI6353, 'PB': LISE607RTPulseBlaster, 'mw_gen_a': AgilentMicrowaveGenerator}

    def _function(self):
        #COMMENT_ME

        self.data['fits'] = None

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

        self.avg_block_number = super(eSensing_N9310A, self)._function()

        # this is to flag when the experiment is done or stopped and this allows proper saving even fit fails
        self.data['exp_finished'] = 1
        self.data['avg_block_number'] = self.avg_block_number

        # turn off laser, apd_switch and MW (just in case)
        self.instruments['PB']['instance'].update({'laser': {'status': False}})
        self.instruments['PB']['instance'].update({'apd_switch': {'status': False}})
        self.instruments['mw_gen_a']['instance'].update({'enable_output': False})
        self.instruments['mw_gen_a']['instance'].update({'enable_IQ': False})
        self.instruments['mw_gen_a']['instance'].update({'enable_modulation': False})

        # start fitting
        if 'counts' in self.data.keys() and 'tau' in self.data.keys() and not 0 in (self.data['counts'][:,1] + self.data['counts'][:, 0]):
            self.data['norm_echo'] = 2. * (- self.data['counts'][:, 1] + self.data['counts'][:, 0]) / (
                    self.data['counts'][:, 1] + self.data['counts'][:, 0])

            tau = self.data['tau']
            counts = self.data['norm_echo']

            fit_success = False

            try:
                fits = fit_exp_decay(tau, counts, offset = True, verbose = True)
                self.data['fits'] = fits
                fit_success = True
            except:
                self.data['fits'] = None
                self.log('T2 fit failed')
                fit_success = False

            if self.settings['decoupling_seq']['do_eSensing'] and not 0 in (self.data['counts'][:, 3] + self.data['counts'][:, 2]):
                self.data['norm_esig'] = 2. * (- self.data['counts'][:, 3] + self.data['counts'][:, 2]) / (
                            self.data['counts'][:, 3] + self.data['counts'][:, 2])

                esig_counts = self.data['norm_esig']
                try:
                    fits_esig = fit_exp_decay(tau, esig_counts, offset = True, verbose = True)
                    self.data['fits_esig'] = fits_esig
                    fit_esig_success = True
                except:
                    self.data['fits_esig'] = None
                    self.log('E-sensing T2 fit failed')
                    fit_esig_success = False

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
            microwave_channel_1 = 'microwave_i'
            microwave_channel_2 = 'microwave_q'
        else:
            microwave_channel_1 = 'microwave_q'
            microwave_channel_2 = 'microwave_i'

        if self.settings['decoupling_seq']['eSensing_type'] == 'sine':
            last_mw_channel = microwave_channel_2
        else:
            last_mw_channel = microwave_channel_1

        e_field_channel = self.settings['decoupling_seq']['e_field_channel']
        pi_time = self.settings['mw_pulses']['pi_pulse_time']
        pi_half_time = self.settings['mw_pulses']['pi_half_pulse_time']
        three_pi_half_time = self.settings['mw_pulses']['3pi_half_pulse_time']

        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']

        number_of_pulse_blocks = self.settings['decoupling_seq']['num_of_pulse_blocks']

        # for tau in tau_list:
        #     pulse_sequence = \
        #     [
        #         Pulse(microwave_channel, laser_off_time, pi_half_time),
        #         Pulse(microwave_channel, laser_off_time + pi_half_time/2. + tau - pi_time/2., pi_time),
        #         Pulse(microwave_channel, laser_off_time + pi_half_time/2. + tau + tau - pi_half_time/2., pi_half_time)
        #     ]
        #
        #     end_of_first_HE = laser_off_time + pi_half_time/2. + tau + tau - pi_half_time/2. + pi_half_time
        #
        #     pulse_sequence += [
        #          Pulse('laser', end_of_first_HE + delay_mw_readout, nv_reset_time),
        #          Pulse('apd_readout', end_of_first_HE + delay_mw_readout + delay_readout, meas_time),
        #          ]
        #
        #     start_of_second_HE = end_of_first_HE + delay_mw_readout + nv_reset_time + laser_off_time
        #
        #     pulse_sequence += \
        #     [
        #         Pulse(microwave_channel, start_of_second_HE, pi_half_time),
        #         Pulse(microwave_channel, start_of_second_HE + pi_half_time/2. + tau - pi_time/2., pi_time),
        #         Pulse(microwave_channel, start_of_second_HE + pi_half_time/2. + tau + tau - pi_half_time/2., three_pi_half_time)
        #     ]
        #
        #     end_of_second_HE = start_of_second_HE + pi_half_time/2. + tau + tau - pi_half_time/2. + pi_half_time
        #
        #     pulse_sequence += [
        #         Pulse('laser', end_of_second_HE + delay_mw_readout, nv_reset_time),
        #         Pulse('apd_readout', end_of_second_HE + delay_mw_readout + delay_readout, meas_time)
        #     ]
        #
        #     pulse_sequence.append(Pulse('apd_switch', 0, end_of_second_HE + delay_mw_readout + nv_reset_time))
        #     # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
        #     # if tau == 0 or tau>=15:
        #     pulse_sequences.append(pulse_sequence)

        if self.settings['decoupling_seq']['type'] == 'spin_echo':
            for tau_total in tau_list:
                tau = tau_total / (1 * self.settings['decoupling_seq']['num_of_pulse_blocks'])
                pulse_sequence = []

                # ECHO SEQUENCE:
                start_of_first_HE = laser_off_time
                # the first pi/2 pulse

                pulse_sequence = [Pulse(microwave_channel_1, start_of_first_HE, pi_half_time)]
                section_begin_time = start_of_first_HE + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                for i in range(0, number_of_pulse_blocks):
                    # pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                    section_begin_time += 1 * tau

                # the second pi/2 pulse and readout
                pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + tau / 2, pi_half_time))
                pulse_sequence.append(Pulse('laser', section_begin_time + tau / 2 + pi_half_time + delay_mw_readout,
                                             nv_reset_time))
                pulse_sequence.append(Pulse('apd_readout',
                                             section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + delay_readout,
                                             meas_time))
                start_of_second_HE = section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time

                # the first pi/2 pulse
                pulse_sequence.append(Pulse(microwave_channel_1, start_of_second_HE, pi_half_time))
                section_begin_time = start_of_second_HE + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                for i in range(0, number_of_pulse_blocks):
                    # pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                    section_begin_time += 1 * tau

                # the second 3*pi/2 pulse and readout
                pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + tau / 2, three_pi_half_time))
                pulse_sequence.append(Pulse('laser',
                                             section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout,
                                             nv_reset_time))
                pulse_sequence.append(Pulse('apd_readout',
                                             section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + delay_readout,
                                             meas_time))

                # E-Sensing SEQUENCE
                if self.settings['decoupling_seq']['do_eSensing']:
                    start_of_DEER = section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time
                    # the first pi/2 pulse
                    pulse_sequence.append(Pulse(microwave_channel_1, start_of_DEER, pi_half_time))
                    section_begin_time = start_of_DEER + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                    for i in range(0, number_of_pulse_blocks):
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                        if number_of_pulse_blocks == 1:
                            pulse_sequence.append(
                                Pulse(e_field_channel, section_begin_time + 1 * tau, tau / 2))
                        elif i%2 == 0:
                            pulse_sequence.append(
                                Pulse(e_field_channel, section_begin_time + 1 * tau, tau))
                        section_begin_time += 1 * tau

                    # the second pi/2 pulse and readout
                    pulse_sequence.append(Pulse(last_mw_channel, section_begin_time + tau / 2, pi_half_time))
                    pulse_sequence.append(Pulse('laser', section_begin_time + tau / 2 + pi_half_time + delay_mw_readout,
                                                nv_reset_time))
                    pulse_sequence.append(Pulse('apd_readout',
                                                section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + delay_readout,
                                                meas_time))

                    start_of_second_DEER = section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time

                    # the first pi/2 pulse
                    pulse_sequence.append(Pulse(microwave_channel_1, start_of_second_DEER, pi_half_time))
                    section_begin_time = start_of_second_DEER + pi_half_time - tau / 2  # for the first pulse, only wait tau/2

                    for i in range(0, number_of_pulse_blocks):
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                        if number_of_pulse_blocks == 1:
                            pulse_sequence.append(
                                Pulse(e_field_channel, section_begin_time + 1 * tau, tau / 2))
                        elif i % 2 == 0:
                            pulse_sequence.append(
                                Pulse(e_field_channel, section_begin_time + 1 * tau, tau))
                        section_begin_time += 1 * tau

                    # the second pi/2 pulse and readout
                    pulse_sequence.append(Pulse(last_mw_channel, section_begin_time + tau / 2, three_pi_half_time))
                    pulse_sequence.append(Pulse('laser', section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout,
                                                nv_reset_time))
                    pulse_sequence.append(Pulse('apd_readout',
                                                section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + delay_readout,
                                                meas_time))

                # Turn on APD switch all the time
                pulse_sequence.append(Pulse('apd_switch', 0,
                                            section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + nv_reset_time))

                # print(pulse_sequence)
                pulse_sequences.append(pulse_sequence)

            print('number of sequences before validation ', len(pulse_sequences))
            return pulse_sequences, tau_list, meas_time
            # return pulse_sequences, self.settings['num_averages'], tau_list, meas_time

        elif self.settings['decoupling_seq']['type'] == 'XY4':

            for tau_total in tau_list:
                tau = tau_total / (4 * self.settings['decoupling_seq']['num_of_pulse_blocks'])
                pulse_sequence = []

                # ECHO SEQUENCE:
                start_of_first_HE = laser_off_time
                # the first pi/2 pulse
                pulse_sequence = [Pulse(microwave_channel_1, start_of_first_HE, pi_half_time)]
                section_begin_time = start_of_first_HE + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                for i in range(0, number_of_pulse_blocks):
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 2 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time / 2, pi_time))

                    section_begin_time += 4 * tau

                # the second 3*pi/2 pulse and readout
                pulse_sequence.append(
                    Pulse(microwave_channel_1, section_begin_time + tau / 2, three_pi_half_time))
                pulse_sequence.append(
                    Pulse('laser', section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout,
                          nv_reset_time))
                pulse_sequence.append(
                    Pulse('apd_readout',
                          section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + delay_readout,
                          meas_time))
                start_of_second_HE = section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time

                # the first pi/2 pulse
                pulse_sequence.append(
                    Pulse(microwave_channel_1, start_of_second_HE, pi_half_time))
                section_begin_time = start_of_second_HE + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                for i in range(0, number_of_pulse_blocks):
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 2 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time / 2, pi_time))
                    section_begin_time += 4 * tau

                # the second pi/2 pulse and readout
                pulse_sequence.append(
                    Pulse(microwave_channel_1, section_begin_time + tau / 2, pi_half_time))
                pulse_sequence.append(
                    Pulse('laser',
                          section_begin_time + tau / 2 + pi_half_time + delay_mw_readout,
                          nv_reset_time))
                pulse_sequence.append(
                    Pulse('apd_readout',
                          section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + delay_readout,
                          meas_time))

                # E-Sensing SEQUENCE
                if self.settings['decoupling_seq']['do_eSensing']:
                    start_of_DEER = section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time
                    # the first pi/2 pulse
                    pulse_sequence.append(Pulse(microwave_channel_1, start_of_DEER, pi_half_time))
                    section_begin_time = start_of_DEER + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                    for i in range(0, number_of_pulse_blocks):
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                        pulse_sequence.append(
                            Pulse(e_field_channel, section_begin_time + 1 * tau, tau))
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 2 * tau - pi_time / 2, pi_time))
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time / 2, pi_time))
                        pulse_sequence.append(
                            Pulse(e_field_channel, section_begin_time + 3 * tau, tau))
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time / 2, pi_time))

                        section_begin_time += 4 * tau

                    # the second 3*pi/2 pulse and readout
                    pulse_sequence.append(
                        Pulse(last_mw_channel, section_begin_time + tau / 2, three_pi_half_time))
                    pulse_sequence.append(
                        Pulse('laser', section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout,
                              nv_reset_time))
                    pulse_sequence.append(
                        Pulse('apd_readout',
                              section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + delay_readout,
                              meas_time))

                    start_of_second_DEER = section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time

                    # the first pi/2 pulse
                    pulse_sequence.append(Pulse(microwave_channel_1, start_of_second_DEER, pi_half_time))
                    section_begin_time = start_of_second_DEER + pi_half_time - tau / 2  # for the first pulse, only wait tau/2

                    for i in range(0, number_of_pulse_blocks):
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                        pulse_sequence.append(
                            Pulse(e_field_channel, section_begin_time + 1 * tau, tau))
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 2 * tau - pi_time / 2, pi_time))
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time / 2, pi_time))
                        pulse_sequence.append(
                            Pulse(e_field_channel, section_begin_time + 3 * tau, tau))
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time / 2, pi_time))
                        section_begin_time += 4 * tau

                    # the second pi/2 pulse and readout

                    pulse_sequence.append(
                        Pulse(last_mw_channel, section_begin_time + tau / 2, pi_half_time))
                    pulse_sequence.append(
                        Pulse('laser',
                              section_begin_time + tau / 2 + pi_half_time + delay_mw_readout,
                              nv_reset_time))
                    pulse_sequence.append(
                        Pulse('apd_readout',
                              section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + delay_readout,
                              meas_time))

                pulse_sequence.append(Pulse('apd_switch', 0,
                                            section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + nv_reset_time))
                pulse_sequences.append(pulse_sequence)

            print('number of sequences before validation ', len(pulse_sequences))
            # return pulse_sequences, self.settings['num_averages'], tau_list, meas_time
            return pulse_sequences, tau_list, meas_time

        elif self.settings['decoupling_seq']['type'] == 'XY8':

            for tau_total in tau_list:
                tau = tau_total / (8 * self.settings['decoupling_seq']['num_of_pulse_blocks'])
                pulse_sequence = []

                # ECHO SEQUENCE:
                start_of_first_HE = laser_off_time

                # the first pi/2 pulse
                pulse_sequence.append(
                    Pulse(microwave_channel_1, start_of_first_HE, pi_half_time))

                section_begin_time = start_of_first_HE + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                for i in range(0, number_of_pulse_blocks):

                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 2 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 5 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 6 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 7 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 8 * tau - pi_time / 2, pi_time))
                    section_begin_time += 8 * tau

                # the second 3*pi/2 pulse and readout
                pulse_sequence.append(
                    Pulse(microwave_channel_1, section_begin_time + tau / 2, three_pi_half_time))
                pulse_sequence.append(
                    Pulse('laser', section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout,
                          nv_reset_time))
                pulse_sequence.append(
                    Pulse('apd_readout',
                          section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + delay_readout,
                          meas_time))

                start_of_second_HE = section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time

                # the first pi/2 pulse
                pulse_sequence.append(
                    Pulse(microwave_channel_1, start_of_second_HE, pi_half_time))
                section_begin_time = start_of_second_HE + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                for i in range(0, number_of_pulse_blocks):

                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 2 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 5 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 6 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 7 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 8 * tau - pi_time / 2, pi_time))
                    section_begin_time += 8 * tau

                # the second pi/2 pulse and readout
                pulse_sequence.append(
                    Pulse(microwave_channel_1, section_begin_time + tau / 2, pi_half_time))
                pulse_sequence.append(
                    Pulse('laser',
                          section_begin_time + tau / 2 + pi_half_time + delay_mw_readout,
                          nv_reset_time))
                pulse_sequence.append(
                    Pulse('apd_readout',
                          section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + delay_readout,
                          meas_time))

                # E-Sensing SEQUENCE
                if self.settings['decoupling_seq']['do_eSensing']:
                    start_of_DEER = section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time
                    # the first pi/2 pulse
                    pulse_sequence.append(Pulse(microwave_channel_1, start_of_DEER, pi_half_time))
                    section_begin_time = start_of_DEER + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                    for i in range(0, number_of_pulse_blocks):
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                        pulse_sequence.append(
                            Pulse(e_field_channel, section_begin_time + 1 * tau, tau))
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 2 * tau - pi_time / 2, pi_time))
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time / 2, pi_time))
                        pulse_sequence.append(
                            Pulse(e_field_channel, section_begin_time + 3 * tau, tau))
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time / 2, pi_time))
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 5 * tau - pi_time / 2, pi_time))
                        pulse_sequence.append(
                            Pulse(e_field_channel, section_begin_time + 5 * tau, tau))
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 6 * tau - pi_time / 2, pi_time))
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 7 * tau - pi_time / 2, pi_time))
                        pulse_sequence.append(
                            Pulse(e_field_channel, section_begin_time + 7 * tau, tau))
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 8 * tau - pi_time / 2, pi_time))
                        section_begin_time += 8 * tau

                    # the second 3*pi/2 pulse and readout
                    pulse_sequence.append(
                        Pulse(last_mw_channel, section_begin_time + tau / 2, three_pi_half_time))
                    pulse_sequence.append(
                        Pulse('laser', section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout,
                              nv_reset_time))
                    pulse_sequence.append(
                        Pulse('apd_readout',
                              section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + delay_readout,
                              meas_time))

                    start_of_second_DEER = section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time

                    # the first pi/2 pulse
                    pulse_sequence.append(Pulse(microwave_channel_1, start_of_second_DEER, pi_half_time))
                    section_begin_time = start_of_second_DEER + pi_half_time - tau / 2  # for the first pulse, only wait tau/2

                    for i in range(0, number_of_pulse_blocks):
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                        pulse_sequence.append(
                            Pulse(e_field_channel, section_begin_time + 1 * tau, tau))
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 2 * tau - pi_time / 2, pi_time))
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time / 2, pi_time))
                        pulse_sequence.append(
                            Pulse(e_field_channel, section_begin_time + 3 * tau, tau))
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time / 2, pi_time))
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 5 * tau - pi_time / 2, pi_time))
                        pulse_sequence.append(
                            Pulse(e_field_channel, section_begin_time + 5 * tau, tau))
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 6 * tau - pi_time / 2, pi_time))
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 7 * tau - pi_time / 2, pi_time))
                        pulse_sequence.append(
                            Pulse(e_field_channel, section_begin_time + 7 * tau, tau))
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 8 * tau - pi_time / 2, pi_time))
                        section_begin_time += 8 * tau

                    # the second pi/2 pulse and readout
                    pulse_sequence.append(
                        Pulse(last_mw_channel, section_begin_time + tau / 2, pi_half_time))
                    pulse_sequence.append(
                        Pulse('laser',
                              section_begin_time + tau / 2 + pi_half_time + delay_mw_readout,
                              nv_reset_time))
                    pulse_sequence.append(
                        Pulse('apd_readout',
                              section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + delay_readout,
                              meas_time))

                pulse_sequence.append(Pulse('apd_switch', 0,
                                            section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + nv_reset_time))
                pulse_sequences.append(pulse_sequence)

            print('number of sequences before validation ', len(pulse_sequences))
            # return pulse_sequences, self.settings['num_averages'], tau_list, meas_time
            return pulse_sequences, tau_list, meas_time

        elif self.settings['decoupling_seq']['type'] == 'CPMG':

            for tau_total in tau_list:
                tau = tau_total / (1 * self.settings['decoupling_seq']['num_of_pulse_blocks'])
                pulse_sequence = []

                # ECHO SEQUENCE:
                start_of_first_HE = laser_off_time

                # the first pi/2 pulse
                pulse_sequence = [Pulse(microwave_channel_1, start_of_first_HE, pi_half_time)]
                section_begin_time = start_of_first_HE + pi_half_time - tau / 2  # for the first pulse, only wait tau/2

                for i in range(0, number_of_pulse_blocks):
                    pulse_sequence.append(Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                    section_begin_time += 1 * tau

                # the second 3*pi/2 pulse and readout
                pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + tau / 2, three_pi_half_time))
                pulse_sequence.append(Pulse('laser', section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout,
                                             nv_reset_time))
                pulse_sequence.append(Pulse('apd_readout',
                                             section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + delay_readout,
                                             meas_time))

                start_of_second_HE = section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time

                # the first pi/2 pulse
                pulse_sequence.append(Pulse(microwave_channel_1, start_of_second_HE, pi_half_time))

                section_begin_time = start_of_second_HE + pi_half_time - tau / 2  # for the first pulse, only wait tau/2

                for i in range(0, number_of_pulse_blocks):
                    pulse_sequence.append(Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                    section_begin_time += 1 * tau

                # the second pi/2 pulse and readout
                pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + tau / 2, pi_half_time))
                pulse_sequence.append(Pulse('laser',
                                             section_begin_time + tau / 2 + pi_half_time + delay_mw_readout,
                                             nv_reset_time))
                pulse_sequence.append(Pulse('apd_readout',
                                             section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + delay_readout,
                                             meas_time))

                # E-Sensing SEQUENCE
                if self.settings['decoupling_seq']['do_eSensing']:
                    start_of_DEER = section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time
                    # the first pi/2 pulse
                    pulse_sequence.append(Pulse(microwave_channel_1, start_of_DEER, pi_half_time))
                    section_begin_time = start_of_DEER + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                    for i in range(0, number_of_pulse_blocks):
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                        if number_of_pulse_blocks == 1:
                            pulse_sequence.append(
                                Pulse(e_field_channel, section_begin_time + 1 * tau, tau / 2))
                        elif i % 2 == 0:
                            pulse_sequence.append(
                                Pulse(e_field_channel, section_begin_time + 1 * tau, tau))
                        section_begin_time += 1 * tau

                    # the second 3*pi/2 pulse and readout
                    pulse_sequence.append(Pulse(last_mw_channel, section_begin_time + tau / 2, three_pi_half_time))
                    pulse_sequence.append(Pulse('laser', section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout,
                                                nv_reset_time))
                    pulse_sequence.append(Pulse('apd_readout',
                                                section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + delay_readout,
                                                meas_time))

                    start_of_second_DEER = section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time

                    # the first pi/2 pulse
                    pulse_sequence.append(Pulse(microwave_channel_1, start_of_second_DEER, pi_half_time))
                    section_begin_time = start_of_second_DEER + pi_half_time - tau / 2  # for the first pulse, only wait tau/2

                    for i in range(0, number_of_pulse_blocks):
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                        if number_of_pulse_blocks == 1:
                            pulse_sequence.append(
                                Pulse(e_field_channel, section_begin_time + 1 * tau, tau / 2))
                        elif i % 2 == 0:
                            pulse_sequence.append(
                                Pulse(e_field_channel, section_begin_time + 1 * tau, tau))
                        section_begin_time += 1 * tau


                    # the second pi/2 pulse and readout
                    pulse_sequence.append(Pulse(last_mw_channel, section_begin_time + tau / 2, pi_half_time))
                    pulse_sequence.append(
                        Pulse('laser', section_begin_time + tau / 2 + pi_half_time + delay_mw_readout,
                              nv_reset_time))
                    pulse_sequence.append(Pulse('apd_readout',
                                                section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + delay_readout,
                                                meas_time))

                # Turn on APD switch all the time
                pulse_sequence.append(Pulse('apd_switch', 0,
                                            section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + nv_reset_time))


    #             # # DEER SEQUENCE
    #             # if self.settings['decoupling_seq']['do_deer']:
    #             #     start_of_DEER = section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time
    #             #     # the first pi/2 pulse
    #             #     pulse_sequence.extend([Pulse(microwave_channel_1, start_of_DEER, pi_half_time)])
    #             #     section_begin_time = start_of_DEER + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
    #             #     for i in range(0, number_of_pulse_blocks):
    #             #         pulse_sequence.extend(
    #             #             [Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time),
    #             #              Pulse('RF_switch', section_begin_time + 1 * tau - RF_pi_time / 2, RF_pi_time)
    #             #              ])
    #             #         section_begin_time += 1 * tau
    #             #
    #             #     # the second pi/2 pulse and readout
    #             #     pulse_sequence.extend([Pulse(microwave_channel_1, section_begin_time + tau / 2, pi_half_time),
    #             #                            Pulse('laser',
    #             #                                  section_begin_time + tau / 2 + pi_half_time + delay_mw_readout,
    #             #                                  nv_reset_time),
    #             #                            Pulse('apd_readout',
    #             #                                  section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + delay_readout,
    #             #                                  meas_time)
    #             #                            ])
    #             #
    #             #     start_of_second_DEER = section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time
    #             #     pulse_sequence.extend([Pulse(microwave_channel_1, start_of_second_DEER, pi_half_time)])
    #             #     section_begin_time = start_of_second_DEER + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
    #             #
    #             #     for i in range(0, number_of_pulse_blocks):
    #             #         pulse_sequence.extend(
    #             #             [Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time),
    #             #              Pulse('RF_switch', section_begin_time + 1 * tau - RF_pi_time / 2, RF_pi_time)
    #             #              ])
    #             #         section_begin_time += 1 * tau
    #             #
    #             #     # the second 3*pi/2 pulse and readout
    #             #     pulse_sequence.extend([Pulse(microwave_channel_1, section_begin_time + tau / 2, three_pi_half_time),
    #             #                            Pulse('laser',
    #             #                                  section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout,
    #             #                                  nv_reset_time),
    #             #                            Pulse('apd_readout',
    #             #                                  section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + delay_readout,
    #             #                                  meas_time)
    #             #                            ])
    #
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

        if 'avg_block_number' in data.keys() is not None and data['avg_block_number'] is not None:
            avg_block_number = data['avg_block_number']
        else:
            avg_block_number = 0

        if 'fits' in data.keys() is not None and data['fits'] is not None:

            # counts = (-data['counts'][:,1] + data['counts'][:,0])/ (data['counts'][:,0] + data['counts'][:,1])
            tau = data['tau']
            fits = data['fits']

            axislist[0].plot(tau, data['norm_echo'])
            tauinterp = np.linspace(np.min(tau), np.max(tau), 100)
            axislist[0].plot(tauinterp, exp_offset(tauinterp, fits[0], fits[1], fits[2]))


            if 'fits_esig' in data.keys() is not None and data['fits_esig'] is not None:
                fits_esig = data['fits_esig']
                axislist[0].plot(tau, data['norm_esig'])
                axislist[0].plot(tauinterp, exp_offset(tauinterp, fits_esig[0], fits_esig[1], fits_esig[2]))
                axislist[0].legend(labels=('Echo', 'Exp Fit','E-Sensing', 'Exp Fit'), fontsize=10)
                axislist[0].set_title(
                    'Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
                        'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s), {:s} electrometry\n after {:d} averages \n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz, \n T2 (p = 1) = {:2.1f} ns, T2 (E-field) = {:2.1f} ns'.format(
                        self.settings['decoupling_seq']['type'],
                        self.settings['decoupling_seq']['num_of_pulse_blocks'],
                        self.settings['decoupling_seq']['eSensing_type'], avg_block_number,
                        self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency'] * 1e-9,
                        fits[1], fits_esig[1]))
            else:
                axislist[0].legend(labels=('Echo', 'Exp Fit'), fontsize=10)
                axislist[0].set_title(
                    'Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
                        'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s)\n after {:d} averages \n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz, \n T2 (p = 1) = {:2.1f} ns'.format(
                        self.settings['decoupling_seq']['type'],
                        self.settings['decoupling_seq']['num_of_pulse_blocks'], avg_block_number,
                        self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency'] * 1e-9,
                        fits[1]))
            axislist[0].set_ylabel('fluorescence contrast')
            axislist[0].set_xlabel('tau [ns]')


        elif 'exp_finished' in data.keys() is not None and data['exp_finished'] == 1:
            if 'norm_echo' in data.keys() is not None and data['norm_echo'] is not None:
                tau = data['tau']
                axislist[0].plot(tau, data['norm_echo'],'d')
                if 'fits_esig' in data.keys() is not None and data['fits_esig'] is not None:
                    fits_esig = data['fits_esig']
                    axislist[0].plot(tau, data['norm_esig'],'o')
                    axislist[0].legend(labels=('Echo', 'E-Sensing'), fontsize=10)
                    axislist[0].set_title(
                        '(final plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
                            'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s), {:s} electrometry\n after {:d} averages \n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz'.format(
                            self.settings['decoupling_seq']['type'],
                            self.settings['decoupling_seq']['num_of_pulse_blocks'],
                            self.settings['decoupling_seq']['eSensing_type'], avg_block_number,
                            self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency'] * 1e-9))
                else:
                    axislist[0].legend(labels=('Echo'), fontsize=10)
                    axislist[0].set_title(
                        '(final plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
                            'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s)\n after {:d} averages \n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz'.format(
                            self.settings['decoupling_seq']['type'],
                            self.settings['decoupling_seq']['num_of_pulse_blocks'],
                            avg_block_number,
                            self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency'] * 1e-9))
                axislist[0].set_ylabel('fluorescence contrast')
                axislist[0].set_xlabel('tau [ns]')
            else:
                axislist[0].plot(data['tau'], data['counts'][:, 0])
                axislist[0].plot(data['tau'], data['counts'][:, 1])

                if self.settings['decoupling_seq']['do_eSensing']:
                    axislist[0].plot(data['tau'], data['counts'][:, 2])
                    axislist[0].plot(data['tau'], data['counts'][:, 3])
                    axislist[0].legend(labels=('Echo up {:.0f}kcps'.format(np.mean(data['counts'][:, 0])),
                                               'Echo down {:.0f}kcps'.format(np.mean(data['counts'][:, 1])),
                                               'E-Sensing up {:.0f}kcps'.format(np.mean(data['counts'][:, 0])),
                                               'E-Sensing down {:.0f}kcps'.format(np.mean(data['counts'][:, 1]))),
                                       fontsize=10)
                    axislist[0].set_title(
                        '(final plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
                            'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s), {:s} electrometry\n after {:d} averages\n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz'.format(
                            self.settings['decoupling_seq']['type'],
                            self.settings['decoupling_seq']['num_of_pulse_blocks'],
                            self.settings['decoupling_seq']['eSensing_type'],
                            avg_block_number, self.settings['mw_pulses']['mw_power'],
                            self.settings['mw_pulses']['mw_frequency'] * 1e-9))

                else:
                    axislist[0].legend(labels=('Echo up {:.0f}kcps'.format(np.mean(data['counts'][:, 0])),
                                           'Echo down {:.0f}kcps'.format(np.mean(data['counts'][:, 1]))), fontsize=10)
                    axislist[0].set_title(
                        '(final plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
                            'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s)\n after {:d} averages\n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz'.format(
                            self.settings['decoupling_seq']['type'], self.settings['decoupling_seq']['num_of_pulse_blocks'],
                            avg_block_number, self.settings['mw_pulses']['mw_power'],
                            self.settings['mw_pulses']['mw_frequency'] * 1e-9))
                axislist[0].set_ylabel('fluorescence [kcps]')
                axislist[0].set_xlabel('tau [ns]')

        else:
            super(eSensing_N9310A, self)._plot(axislist)

            if self.settings['decoupling_seq']['do_eSensing']:
                axislist[0].legend(labels=('Echo up {:.0f}kcps'.format(np.mean(data['counts'][:, 0])),
                                           'Echo down {:.0f}kcps'.format(np.mean(data['counts'][:, 1])),
                                           'E-Sensing up {:.0f}kcps'.format(np.mean(data['counts'][:, 0])),
                                           'E-Sensing down {:.0f}kcps'.format(np.mean(data['counts'][:, 1]))),
                                   fontsize=10)
                axislist[0].set_title(
                    '(initial plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
                        'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s), {:s} electrometry\n after {:d} averages\n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz'.format(
                        self.settings['decoupling_seq']['type'], self.settings['decoupling_seq']['num_of_pulse_blocks'],
                        self.settings['decoupling_seq']['eSensing_type'],
                        avg_block_number, self.settings['mw_pulses']['mw_power'],
                        self.settings['mw_pulses']['mw_frequency'] * 1e-9))

            else:
                axislist[0].legend(labels=('Echo up {:.0f}kcps'.format(np.mean(data['counts'][:, 0])),
                                           'Echo down {:.0f}kcps'.format(np.mean(data['counts'][:, 1]))), fontsize=10)
                axislist[0].set_title(
                    '(initial plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
                        'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s) \n after {:d} averages\n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz'.format(
                        self.settings['decoupling_seq']['type'], self.settings['decoupling_seq']['num_of_pulse_blocks'],
                        avg_block_number, self.settings['mw_pulses']['mw_power'],
                        self.settings['mw_pulses']['mw_frequency'] * 1e-9))

    def _update_plot(self, axislist):

        if len(axislist[0].lines) == 0:
            self._plot(axislist)
            return

        super(eSensing_N9310A, self)._update_plot(axislist)

        if self.settings['decoupling_seq']['do_eSensing']:
            axislist[0].legend(labels=('Echo up {:.0f}kcps'.format(np.mean(self.data['counts'][:, 0])),
                                       'Echo down {:.0f}kcps'.format(np.mean(self.data['counts'][:, 1])),
                                       'E-Sensing up {:.0f}kcps'.format(np.mean(self.data['counts'][:, 0])),
                                       'E-Sensing down {:.0f}kcps'.format(np.mean(self.data['counts'][:, 1]))),
                               fontsize=10)
            axislist[0].set_title(
                '(updating plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
                    'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s), {:s} electrometry\n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz'.format(
                    self.settings['decoupling_seq']['type'], self.settings['decoupling_seq']['num_of_pulse_blocks'],
                    self.settings['decoupling_seq']['eSensing_type'],
                    self.settings['mw_pulses']['mw_power'],
                    self.settings['mw_pulses']['mw_frequency'] * 1e-9))

        else:
            axislist[0].legend(labels=('Echo up {:.0f}kcps'.format(np.mean(self.data['counts'][:, 0])),
                                       'Echo down {:.0f}kcps'.format(np.mean(self.data['counts'][:, 1]))), fontsize=10)

            axislist[0].set_title('(updating plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
                'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s) \n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz'.format(
                self.settings['decoupling_seq']['type'], self.settings['decoupling_seq']['num_of_pulse_blocks'],
                self.settings['mw_pulses']['mw_power'],
                self.settings['mw_pulses']['mw_frequency'] * 1e-9))

class PDD_XYreadout(PulsedExperimentBaseScript):
    """
    This script runs a PDD ( Periodic Dynamical Decoupling) sequence for different number of pi pulses.
    Applied an AC electric field to be detected, which is controled by the PB14 (aux).
    To symmetrize the sequence between the 0 and +/-1 state we reinitialize every time.
    Keysight N9310A generator is used and it has IQ modulation.
    For a single pi-pulse this is a Hahn-echo sequence. For zero pulses this is a Ramsey sequence.
    The sequence is pi/2 - tau/4 - (tau/4 - pi  - tau/4)^n - tau/4 - pi/2
    Tau/2 is the time between the center of the pulses!

    At the end of the sequence, you can choose whether or not to do both X and Y readout.

    --> Ziwei Qiu 9/2/2019

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
            Parameter('time_step', 100., [2.5, 5., 10., 20., 50., 100., 200., 300., 400., 500., 600., 800., 1000.,2000.,5000., 8000., 10000.,50000., 100000., 500000.],
                      'time step increment of time between pi pulses (in ns)')
        ]),
        Parameter('do_Y_readout', True, bool, 'choose whether or not to do a Y readout'),
        Parameter('decoupling_seq', [
            Parameter('type', 'CPMG', ['spin_echo', 'CPMG', 'XY4', 'XY8'], 'type of dynamical decoupling sequences'),
            Parameter('num_of_pulse_blocks', 1, int, 'number of pulse blocks.'),
        ]),
        Parameter('do_normalization', False, bool, 'choose whether or not to do normalization on the data'),
        Parameter('read_out', [
            Parameter('meas_time', 200, float, 'measurement time after rabi sequence (in ns)'),
            Parameter('nv_reset_time', 2000, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 600, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 550, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('min_pulse_dur', 15, [10, 15, 50], 'Minimum allowed pulse duration (ns)'),
        Parameter('num_averages', 500000, int, 'number of averages (>100000)'),
    ]

    _INSTRUMENTS = {'daq': NI6353, 'PB': LISE607RTPulseBlaster, 'mw_gen_a': AgilentMicrowaveGenerator}

    def _function(self):
        #COMMENT_ME



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

        self.avg_block_number = super(PDD_XYreadout, self)._function()

        # this is to flag when the experiment is done or stopped and this allows proper saving even fit fails
        self.data['exp_finished'] = 1
        self.data['avg_block_number'] = self.avg_block_number

        # turn off laser, apd_switch and MW (just in case)
        self.instruments['PB']['instance'].update({'laser': {'status': False}})
        self.instruments['PB']['instance'].update({'apd_switch': {'status': False}})
        self.instruments['mw_gen_a']['instance'].update({'enable_output': False})
        self.instruments['mw_gen_a']['instance'].update({'enable_IQ': False})
        self.instruments['mw_gen_a']['instance'].update({'enable_modulation': False})

        # start fitting
        # self.data['fits'] = None
        self.data['norm_echo'] = None
        self.data['norm_echo_Y'] = None
        if self.settings['do_normalization'] is True:
            if 'counts' in self.data.keys() and 'tau' in self.data.keys() and not 0 in (self.data['counts'][:,1] + self.data['counts'][:, 0]):

                self.data['norm_echo'] = 2. * (- self.data['counts'][:, 1] + self.data['counts'][:, 0]) / (
                        self.data['counts'][:, 1] + self.data['counts'][:, 0])

                # tau = self.data['tau']
                # counts = self.data['norm_echo']

                # fit_success = False

                # try:
                #     fits = fit_exp_decay(tau, counts, offset = True, verbose = True)
                #     self.data['fits'] = fits
                #     fit_success = True
                # except:
                #     self.data['fits'] = None
                #     self.log('T2 fit failed')
                #     fit_success = False

                if self.settings['do_Y_readout'] and not 0 in (self.data['counts'][:, 3] + self.data['counts'][:, 2]):
                    self.data['norm_echo_Y'] = 2. * (- self.data['counts'][:, 3] + self.data['counts'][:, 2]) / (
                                self.data['counts'][:, 3] + self.data['counts'][:, 2])

                    # esig_counts = self.data['norm_esig']
                    # try:
                    #     fits_esig = fit_exp_decay(tau, esig_counts, offset = True, verbose = True)
                    #     self.data['fits_esig'] = fits_esig
                    #     fit_esig_success = True
                    # except:
                    #     self.data['fits_esig'] = None
                    #     self.log('E-sensing T2 fit failed')
                    #     fit_esig_success = False

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
            microwave_channel_1 = 'microwave_i'
            microwave_channel_2 = 'microwave_q'
        else:
            microwave_channel_1 = 'microwave_q'
            microwave_channel_2 = 'microwave_i'

        # if self.settings['decoupling_seq']['eSensing_type'] == 'sine':
        #     last_mw_channel = microwave_channel_2
        # else:
        #     last_mw_channel = microwave_channel_1
        last_mw_channel = microwave_channel_2

        # e_field_channel = self.settings['decoupling_seq']['e_field_channel']
        pi_time = self.settings['mw_pulses']['pi_pulse_time']
        pi_half_time = self.settings['mw_pulses']['pi_half_pulse_time']
        three_pi_half_time = self.settings['mw_pulses']['3pi_half_pulse_time']

        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        number_of_pulse_blocks = self.settings['decoupling_seq']['num_of_pulse_blocks']

        # for tau in tau_list:
        #     pulse_sequence = \
        #     [
        #         Pulse(microwave_channel, laser_off_time, pi_half_time),
        #         Pulse(microwave_channel, laser_off_time + pi_half_time/2. + tau - pi_time/2., pi_time),
        #         Pulse(microwave_channel, laser_off_time + pi_half_time/2. + tau + tau - pi_half_time/2., pi_half_time)
        #     ]
        #
        #     end_of_first_HE = laser_off_time + pi_half_time/2. + tau + tau - pi_half_time/2. + pi_half_time
        #
        #     pulse_sequence += [
        #          Pulse('laser', end_of_first_HE + delay_mw_readout, nv_reset_time),
        #          Pulse('apd_readout', end_of_first_HE + delay_mw_readout + delay_readout, meas_time),
        #          ]
        #
        #     start_of_second_HE = end_of_first_HE + delay_mw_readout + nv_reset_time + laser_off_time
        #
        #     pulse_sequence += \
        #     [
        #         Pulse(microwave_channel, start_of_second_HE, pi_half_time),
        #         Pulse(microwave_channel, start_of_second_HE + pi_half_time/2. + tau - pi_time/2., pi_time),
        #         Pulse(microwave_channel, start_of_second_HE + pi_half_time/2. + tau + tau - pi_half_time/2., three_pi_half_time)
        #     ]
        #
        #     end_of_second_HE = start_of_second_HE + pi_half_time/2. + tau + tau - pi_half_time/2. + pi_half_time
        #
        #     pulse_sequence += [
        #         Pulse('laser', end_of_second_HE + delay_mw_readout, nv_reset_time),
        #         Pulse('apd_readout', end_of_second_HE + delay_mw_readout + delay_readout, meas_time)
        #     ]
        #
        #     pulse_sequence.append(Pulse('apd_switch', 0, end_of_second_HE + delay_mw_readout + nv_reset_time))
        #     # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
        #     # if tau == 0 or tau>=15:
        #     pulse_sequences.append(pulse_sequence)

        if self.settings['decoupling_seq']['type'] == 'spin_echo':
            for tau_total in tau_list:
                tau = tau_total / (1 * self.settings['decoupling_seq']['num_of_pulse_blocks'])
                pulse_sequence = []

                # ECHO SEQUENCE:
                start_of_first_HE = laser_off_time
                # the first pi/2 pulse

                pulse_sequence = [Pulse(microwave_channel_1, start_of_first_HE, pi_half_time)]
                section_begin_time = start_of_first_HE + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                for i in range(0, number_of_pulse_blocks):
                    # pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                    section_begin_time += 1 * tau

                # the second pi/2 pulse and readout
                pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + tau / 2, pi_half_time))
                pulse_sequence.append(Pulse('laser', section_begin_time + tau / 2 + pi_half_time + delay_mw_readout,
                                             nv_reset_time))
                pulse_sequence.append(Pulse('apd_readout',
                                             section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + delay_readout,
                                             meas_time))
                start_of_second_HE = section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time

                # the first pi/2 pulse
                pulse_sequence.append(Pulse(microwave_channel_1, start_of_second_HE, pi_half_time))
                section_begin_time = start_of_second_HE + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                for i in range(0, number_of_pulse_blocks):
                    # pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                    section_begin_time += 1 * tau

                # the second 3*pi/2 pulse and readout
                pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + tau / 2, three_pi_half_time))
                pulse_sequence.append(Pulse('laser',
                                             section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout,
                                             nv_reset_time))
                pulse_sequence.append(Pulse('apd_readout',
                                             section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + delay_readout,
                                             meas_time))

                # Y Echo
                if self.settings['do_Y_readout']:
                    start_of_DEER = section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time
                    # the first pi/2 pulse
                    pulse_sequence.append(Pulse(microwave_channel_1, start_of_DEER, pi_half_time))
                    section_begin_time = start_of_DEER + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                    for i in range(0, number_of_pulse_blocks):
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                        # if number_of_pulse_blocks == 1:
                        #     pulse_sequence.append(
                        #         Pulse(e_field_channel, section_begin_time + 1 * tau, tau / 2))
                        # elif i%2 == 0:
                        #     pulse_sequence.append(
                        #         Pulse(e_field_channel, section_begin_time + 1 * tau, tau))
                        section_begin_time += 1 * tau

                    # the second pi/2 pulse and readout
                    pulse_sequence.append(Pulse(last_mw_channel, section_begin_time + tau / 2, pi_half_time))
                    pulse_sequence.append(Pulse('laser', section_begin_time + tau / 2 + pi_half_time + delay_mw_readout,
                                                nv_reset_time))
                    pulse_sequence.append(Pulse('apd_readout',
                                                section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + delay_readout,
                                                meas_time))

                    start_of_second_DEER = section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time

                    # the first pi/2 pulse
                    pulse_sequence.append(Pulse(microwave_channel_1, start_of_second_DEER, pi_half_time))
                    section_begin_time = start_of_second_DEER + pi_half_time - tau / 2  # for the first pulse, only wait tau/2

                    for i in range(0, number_of_pulse_blocks):
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                        # if number_of_pulse_blocks == 1:
                        #     pulse_sequence.append(
                        #         Pulse(e_field_channel, section_begin_time + 1 * tau, tau / 2))
                        # elif i % 2 == 0:
                        #     pulse_sequence.append(
                        #         Pulse(e_field_channel, section_begin_time + 1 * tau, tau))
                        section_begin_time += 1 * tau

                    # the second pi/2 pulse and readout
                    pulse_sequence.append(Pulse(last_mw_channel, section_begin_time + tau / 2, three_pi_half_time))
                    pulse_sequence.append(Pulse('laser', section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout,
                                                nv_reset_time))
                    pulse_sequence.append(Pulse('apd_readout',
                                                section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + delay_readout,
                                                meas_time))

                # Turn on APD switch all the time
                pulse_sequence.append(Pulse('apd_switch', 0,
                                            section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + nv_reset_time))

                # print(pulse_sequence)
                pulse_sequences.append(pulse_sequence)

            print('number of sequences before validation ', len(pulse_sequences))
            return pulse_sequences, tau_list, meas_time
            # return pulse_sequences, self.settings['num_averages'], tau_list, meas_time

        elif self.settings['decoupling_seq']['type'] == 'XY4':

            for tau_total in tau_list:
                tau = tau_total / (4 * self.settings['decoupling_seq']['num_of_pulse_blocks'])
                pulse_sequence = []

                # ECHO SEQUENCE:
                start_of_first_HE = laser_off_time
                # the first pi/2 pulse
                pulse_sequence = [Pulse(microwave_channel_1, start_of_first_HE, pi_half_time)]
                section_begin_time = start_of_first_HE + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                for i in range(0, number_of_pulse_blocks):
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 2 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time / 2, pi_time))

                    section_begin_time += 4 * tau

                # the second 3*pi/2 pulse and readout
                pulse_sequence.append(
                    Pulse(microwave_channel_1, section_begin_time + tau / 2, three_pi_half_time))
                pulse_sequence.append(
                    Pulse('laser', section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout,
                          nv_reset_time))
                pulse_sequence.append(
                    Pulse('apd_readout',
                          section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + delay_readout,
                          meas_time))
                start_of_second_HE = section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time

                # the first pi/2 pulse
                pulse_sequence.append(
                    Pulse(microwave_channel_1, start_of_second_HE, pi_half_time))
                section_begin_time = start_of_second_HE + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                for i in range(0, number_of_pulse_blocks):
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 2 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time / 2, pi_time))
                    section_begin_time += 4 * tau

                # the second pi/2 pulse and readout
                pulse_sequence.append(
                    Pulse(microwave_channel_1, section_begin_time + tau / 2, pi_half_time))
                pulse_sequence.append(
                    Pulse('laser',
                          section_begin_time + tau / 2 + pi_half_time + delay_mw_readout,
                          nv_reset_time))
                pulse_sequence.append(
                    Pulse('apd_readout',
                          section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + delay_readout,
                          meas_time))

                # Y Echo SEQUENCE
                if self.settings['do_Y_readout']:
                    start_of_DEER = section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time
                    # the first pi/2 pulse
                    pulse_sequence.append(Pulse(microwave_channel_1, start_of_DEER, pi_half_time))
                    section_begin_time = start_of_DEER + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                    for i in range(0, number_of_pulse_blocks):
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                        # pulse_sequence.append(
                        #     Pulse(e_field_channel, section_begin_time + 1 * tau, tau))
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 2 * tau - pi_time / 2, pi_time))
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time / 2, pi_time))
                        # pulse_sequence.append(
                        #     Pulse(e_field_channel, section_begin_time + 3 * tau, tau))
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time / 2, pi_time))

                        section_begin_time += 4 * tau

                    # the second 3*pi/2 pulse and readout
                    pulse_sequence.append(
                        Pulse(last_mw_channel, section_begin_time + tau / 2, three_pi_half_time))
                    pulse_sequence.append(
                        Pulse('laser', section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout,
                              nv_reset_time))
                    pulse_sequence.append(
                        Pulse('apd_readout',
                              section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + delay_readout,
                              meas_time))

                    start_of_second_DEER = section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time

                    # the first pi/2 pulse
                    pulse_sequence.append(Pulse(microwave_channel_1, start_of_second_DEER, pi_half_time))
                    section_begin_time = start_of_second_DEER + pi_half_time - tau / 2  # for the first pulse, only wait tau/2

                    for i in range(0, number_of_pulse_blocks):
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                        # pulse_sequence.append(
                        #     Pulse(e_field_channel, section_begin_time + 1 * tau, tau))
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 2 * tau - pi_time / 2, pi_time))
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time / 2, pi_time))
                        # pulse_sequence.append(
                        #     Pulse(e_field_channel, section_begin_time + 3 * tau, tau))
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time / 2, pi_time))
                        section_begin_time += 4 * tau

                    # the second pi/2 pulse and readout

                    pulse_sequence.append(
                        Pulse(last_mw_channel, section_begin_time + tau / 2, pi_half_time))
                    pulse_sequence.append(
                        Pulse('laser',
                              section_begin_time + tau / 2 + pi_half_time + delay_mw_readout,
                              nv_reset_time))
                    pulse_sequence.append(
                        Pulse('apd_readout',
                              section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + delay_readout,
                              meas_time))

                pulse_sequence.append(Pulse('apd_switch', 0,
                                            section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + nv_reset_time))
                pulse_sequences.append(pulse_sequence)

            print('number of sequences before validation ', len(pulse_sequences))
            # return pulse_sequences, self.settings['num_averages'], tau_list, meas_time
            return pulse_sequences, tau_list, meas_time

        elif self.settings['decoupling_seq']['type'] == 'XY8':

            for tau_total in tau_list:
                tau = tau_total / (8 * self.settings['decoupling_seq']['num_of_pulse_blocks'])
                pulse_sequence = []

                # ECHO SEQUENCE:
                start_of_first_HE = laser_off_time

                # the first pi/2 pulse
                pulse_sequence.append(
                    Pulse(microwave_channel_1, start_of_first_HE, pi_half_time))

                section_begin_time = start_of_first_HE + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                for i in range(0, number_of_pulse_blocks):

                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 2 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 5 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 6 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 7 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 8 * tau - pi_time / 2, pi_time))
                    section_begin_time += 8 * tau

                # the second 3*pi/2 pulse and readout
                pulse_sequence.append(
                    Pulse(microwave_channel_1, section_begin_time + tau / 2, three_pi_half_time))
                pulse_sequence.append(
                    Pulse('laser', section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout,
                          nv_reset_time))
                pulse_sequence.append(
                    Pulse('apd_readout',
                          section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + delay_readout,
                          meas_time))

                start_of_second_HE = section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time

                # the first pi/2 pulse
                pulse_sequence.append(
                    Pulse(microwave_channel_1, start_of_second_HE, pi_half_time))
                section_begin_time = start_of_second_HE + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                for i in range(0, number_of_pulse_blocks):

                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 2 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 5 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 6 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 7 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 8 * tau - pi_time / 2, pi_time))
                    section_begin_time += 8 * tau

                # the second pi/2 pulse and readout
                pulse_sequence.append(
                    Pulse(microwave_channel_1, section_begin_time + tau / 2, pi_half_time))
                pulse_sequence.append(
                    Pulse('laser',
                          section_begin_time + tau / 2 + pi_half_time + delay_mw_readout,
                          nv_reset_time))
                pulse_sequence.append(
                    Pulse('apd_readout',
                          section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + delay_readout,
                          meas_time))

                # Y Echo SEQUENCE
                if self.settings['do_Y_readout']:
                    start_of_DEER = section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time
                    # the first pi/2 pulse
                    pulse_sequence.append(Pulse(microwave_channel_1, start_of_DEER, pi_half_time))
                    section_begin_time = start_of_DEER + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                    for i in range(0, number_of_pulse_blocks):
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                        # pulse_sequence.append(
                        #     Pulse(e_field_channel, section_begin_time + 1 * tau, tau))
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 2 * tau - pi_time / 2, pi_time))
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time / 2, pi_time))
                        # pulse_sequence.append(
                        #     Pulse(e_field_channel, section_begin_time + 3 * tau, tau))
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time / 2, pi_time))
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 5 * tau - pi_time / 2, pi_time))
                        # pulse_sequence.append(
                        #     Pulse(e_field_channel, section_begin_time + 5 * tau, tau))
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 6 * tau - pi_time / 2, pi_time))
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 7 * tau - pi_time / 2, pi_time))
                        # pulse_sequence.append(
                        #     Pulse(e_field_channel, section_begin_time + 7 * tau, tau))
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 8 * tau - pi_time / 2, pi_time))
                        section_begin_time += 8 * tau

                    # the second 3*pi/2 pulse and readout
                    pulse_sequence.append(
                        Pulse(last_mw_channel, section_begin_time + tau / 2, three_pi_half_time))
                    pulse_sequence.append(
                        Pulse('laser', section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout,
                              nv_reset_time))
                    pulse_sequence.append(
                        Pulse('apd_readout',
                              section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + delay_readout,
                              meas_time))

                    start_of_second_DEER = section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time

                    # the first pi/2 pulse
                    pulse_sequence.append(Pulse(microwave_channel_1, start_of_second_DEER, pi_half_time))
                    section_begin_time = start_of_second_DEER + pi_half_time - tau / 2  # for the first pulse, only wait tau/2

                    for i in range(0, number_of_pulse_blocks):
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                        # pulse_sequence.append(
                        #     Pulse(e_field_channel, section_begin_time + 1 * tau, tau))
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 2 * tau - pi_time / 2, pi_time))
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time / 2, pi_time))
                        # pulse_sequence.append(
                        #     Pulse(e_field_channel, section_begin_time + 3 * tau, tau))
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time / 2, pi_time))
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 5 * tau - pi_time / 2, pi_time))
                        # pulse_sequence.append(
                        #     Pulse(e_field_channel, section_begin_time + 5 * tau, tau))
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 6 * tau - pi_time / 2, pi_time))
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 7 * tau - pi_time / 2, pi_time))
                        # pulse_sequence.append(
                        #     Pulse(e_field_channel, section_begin_time + 7 * tau, tau))
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 8 * tau - pi_time / 2, pi_time))
                        section_begin_time += 8 * tau

                    # the second pi/2 pulse and readout
                    pulse_sequence.append(
                        Pulse(last_mw_channel, section_begin_time + tau / 2, pi_half_time))
                    pulse_sequence.append(
                        Pulse('laser',
                              section_begin_time + tau / 2 + pi_half_time + delay_mw_readout,
                              nv_reset_time))
                    pulse_sequence.append(
                        Pulse('apd_readout',
                              section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + delay_readout,
                              meas_time))

                pulse_sequence.append(Pulse('apd_switch', 0,
                                            section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + nv_reset_time))
                pulse_sequences.append(pulse_sequence)

            print('number of sequences before validation ', len(pulse_sequences))
            # return pulse_sequences, self.settings['num_averages'], tau_list, meas_time
            return pulse_sequences, tau_list, meas_time

        elif self.settings['decoupling_seq']['type'] == 'CPMG':

            for tau_total in tau_list:
                tau = tau_total / (1 * self.settings['decoupling_seq']['num_of_pulse_blocks'])
                pulse_sequence = []

                # ECHO SEQUENCE:
                start_of_first_HE = laser_off_time

                # the first pi/2 pulse
                pulse_sequence = [Pulse(microwave_channel_1, start_of_first_HE, pi_half_time)]
                section_begin_time = start_of_first_HE + pi_half_time - tau / 2  # for the first pulse, only wait tau/2

                for i in range(0, number_of_pulse_blocks):
                    pulse_sequence.append(Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                    section_begin_time += 1 * tau

                # the second 3*pi/2 pulse and readout
                pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + tau / 2, three_pi_half_time))
                pulse_sequence.append(Pulse('laser', section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout,
                                             nv_reset_time))
                pulse_sequence.append(Pulse('apd_readout',
                                             section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + delay_readout,
                                             meas_time))

                start_of_second_HE = section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time

                # the first pi/2 pulse
                pulse_sequence.append(Pulse(microwave_channel_1, start_of_second_HE, pi_half_time))

                section_begin_time = start_of_second_HE + pi_half_time - tau / 2  # for the first pulse, only wait tau/2

                for i in range(0, number_of_pulse_blocks):
                    pulse_sequence.append(Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                    section_begin_time += 1 * tau

                # the second pi/2 pulse and readout
                pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + tau / 2, pi_half_time))
                pulse_sequence.append(Pulse('laser',
                                             section_begin_time + tau / 2 + pi_half_time + delay_mw_readout,
                                             nv_reset_time))
                pulse_sequence.append(Pulse('apd_readout',
                                             section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + delay_readout,
                                             meas_time))

                # Y Echo SEQUENCE
                if self.settings['do_Y_readout']:
                    start_of_DEER = section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time
                    # the first pi/2 pulse
                    pulse_sequence.append(Pulse(microwave_channel_1, start_of_DEER, pi_half_time))
                    section_begin_time = start_of_DEER + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                    for i in range(0, number_of_pulse_blocks):
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                        # if number_of_pulse_blocks == 1:
                        #     pulse_sequence.append(
                        #         Pulse(e_field_channel, section_begin_time + 1 * tau, tau / 2))
                        # elif i % 2 == 0:
                        #     pulse_sequence.append(
                        #         Pulse(e_field_channel, section_begin_time + 1 * tau, tau))
                        section_begin_time += 1 * tau

                    # the second 3*pi/2 pulse and readout
                    pulse_sequence.append(Pulse(last_mw_channel, section_begin_time + tau / 2, three_pi_half_time))
                    pulse_sequence.append(Pulse('laser', section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout,
                                                nv_reset_time))
                    pulse_sequence.append(Pulse('apd_readout',
                                                section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + delay_readout,
                                                meas_time))

                    start_of_second_DEER = section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time

                    # the first pi/2 pulse
                    pulse_sequence.append(Pulse(microwave_channel_1, start_of_second_DEER, pi_half_time))
                    section_begin_time = start_of_second_DEER + pi_half_time - tau / 2  # for the first pulse, only wait tau/2

                    for i in range(0, number_of_pulse_blocks):
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                        # if number_of_pulse_blocks == 1:
                        #     pulse_sequence.append(
                        #         Pulse(e_field_channel, section_begin_time + 1 * tau, tau / 2))
                        # elif i % 2 == 0:
                        #     pulse_sequence.append(
                        #         Pulse(e_field_channel, section_begin_time + 1 * tau, tau))
                        section_begin_time += 1 * tau


                    # the second pi/2 pulse and readout
                    pulse_sequence.append(Pulse(last_mw_channel, section_begin_time + tau / 2, pi_half_time))
                    pulse_sequence.append(
                        Pulse('laser', section_begin_time + tau / 2 + pi_half_time + delay_mw_readout,
                              nv_reset_time))
                    pulse_sequence.append(Pulse('apd_readout',
                                                section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + delay_readout,
                                                meas_time))

                # Turn on APD switch all the time
                pulse_sequence.append(Pulse('apd_switch', 0,
                                            section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + nv_reset_time))


    #             # # DEER SEQUENCE
    #             # if self.settings['decoupling_seq']['do_deer']:
    #             #     start_of_DEER = section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time
    #             #     # the first pi/2 pulse
    #             #     pulse_sequence.extend([Pulse(microwave_channel_1, start_of_DEER, pi_half_time)])
    #             #     section_begin_time = start_of_DEER + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
    #             #     for i in range(0, number_of_pulse_blocks):
    #             #         pulse_sequence.extend(
    #             #             [Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time),
    #             #              Pulse('RF_switch', section_begin_time + 1 * tau - RF_pi_time / 2, RF_pi_time)
    #             #              ])
    #             #         section_begin_time += 1 * tau
    #             #
    #             #     # the second pi/2 pulse and readout
    #             #     pulse_sequence.extend([Pulse(microwave_channel_1, section_begin_time + tau / 2, pi_half_time),
    #             #                            Pulse('laser',
    #             #                                  section_begin_time + tau / 2 + pi_half_time + delay_mw_readout,
    #             #                                  nv_reset_time),
    #             #                            Pulse('apd_readout',
    #             #                                  section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + delay_readout,
    #             #                                  meas_time)
    #             #                            ])
    #             #
    #             #     start_of_second_DEER = section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time
    #             #     pulse_sequence.extend([Pulse(microwave_channel_1, start_of_second_DEER, pi_half_time)])
    #             #     section_begin_time = start_of_second_DEER + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
    #             #
    #             #     for i in range(0, number_of_pulse_blocks):
    #             #         pulse_sequence.extend(
    #             #             [Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time),
    #             #              Pulse('RF_switch', section_begin_time + 1 * tau - RF_pi_time / 2, RF_pi_time)
    #             #              ])
    #             #         section_begin_time += 1 * tau
    #             #
    #             #     # the second 3*pi/2 pulse and readout
    #             #     pulse_sequence.extend([Pulse(microwave_channel_1, section_begin_time + tau / 2, three_pi_half_time),
    #             #                            Pulse('laser',
    #             #                                  section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout,
    #             #                                  nv_reset_time),
    #             #                            Pulse('apd_readout',
    #             #                                  section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + delay_readout,
    #             #                                  meas_time)
    #             #                            ])
    #
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

        if 'avg_block_number' in data.keys() is not None and data['avg_block_number'] is not None:
            avg_block_number = data['avg_block_number']
        else:
            avg_block_number = 0

        # if 'fits' in data.keys() is not None and data['fits'] is not None:
        #
        #     # counts = (-data['counts'][:,1] + data['counts'][:,0])/ (data['counts'][:,0] + data['counts'][:,1])
        #     tau = data['tau']
        #     fits = data['fits']
        #
        #     axislist[0].plot(tau, data['norm_echo'])
        #     tauinterp = np.linspace(np.min(tau), np.max(tau), 100)
        #     axislist[0].plot(tauinterp, exp_offset(tauinterp, fits[0], fits[1], fits[2]))
        #
        #
        #     if 'fits_esig' in data.keys() is not None and data['fits_esig'] is not None:
        #         fits_esig = data['fits_esig']
        #         axislist[0].plot(tau, data['norm_esig'])
        #         axislist[0].plot(tauinterp, exp_offset(tauinterp, fits_esig[0], fits_esig[1], fits_esig[2]))
        #         axislist[0].legend(labels=('Echo', 'Exp Fit','E-Sensing', 'Exp Fit'), fontsize=10)
        #         axislist[0].set_title(
        #             'Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
        #                 'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s), {:s} electrometry\n after {:d} averages \n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz, \n T2 (p = 1) = {:2.1f} ns, T2 (E-field) = {:2.1f} ns'.format(
        #                 self.settings['decoupling_seq']['type'],
        #                 self.settings['decoupling_seq']['num_of_pulse_blocks'],
        #                 self.settings['decoupling_seq']['eSensing_type'], avg_block_number,
        #                 self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency'] * 1e-9,
        #                 fits[1], fits_esig[1]))
        #     else:
        #         axislist[0].legend(labels=('Echo', 'Exp Fit'), fontsize=10)
        #         axislist[0].set_title(
        #             'Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
        #                 'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s)\n after {:d} averages \n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz, \n T2 (p = 1) = {:2.1f} ns'.format(
        #                 self.settings['decoupling_seq']['type'],
        #                 self.settings['decoupling_seq']['num_of_pulse_blocks'], avg_block_number,
        #                 self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency'] * 1e-9,
        #                 fits[1]))
        #     axislist[0].set_ylabel('fluorescence contrast')
        #     axislist[0].set_xlabel('tau [ns]')


        if 'exp_finished' in data.keys() is not None and data['exp_finished'] == 1:
            if 'norm_echo' in data.keys() is not None and data['norm_echo'] is not None:
                tau = data['tau']
                axislist[0].plot(tau, data['norm_echo'])
                if 'norm_echo_Y' in data.keys() is not None and data['norm_echo_Y'] is not None:
                    # fits_esig = data['fits_esig']
                    axislist[0].plot(tau, data['norm_echo_Y'])
                    axislist[0].legend(labels=('X Echo', 'Y Echo'), fontsize=10)
                    axislist[0].set_title(
                        '(final plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
                            'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s)\n after {:d} averages \n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz'.format(
                            self.settings['decoupling_seq']['type'],
                            self.settings['decoupling_seq']['num_of_pulse_blocks'], avg_block_number,
                            self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency'] * 1e-9))
                else:
                    axislist[0].legend(labels=('X Echo'), fontsize=10)
                    axislist[0].set_title(
                        '(final plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
                            'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s)\n after {:d} averages \n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz'.format(
                            self.settings['decoupling_seq']['type'],
                            self.settings['decoupling_seq']['num_of_pulse_blocks'],
                            avg_block_number,
                            self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency'] * 1e-9))
                axislist[0].set_ylabel('fluorescence contrast')
                axislist[0].set_xlabel('tau [ns]')
            else:
                axislist[0].plot(data['tau'], data['counts'][:, 0])
                axislist[0].plot(data['tau'], data['counts'][:, 1])

                if self.settings['do_Y_readout']:
                    axislist[0].plot(data['tau'], data['counts'][:, 2])
                    axislist[0].plot(data['tau'], data['counts'][:, 3])
                    axislist[0].legend(labels=('X Echo up {:.0f}kcps'.format(np.mean(data['counts'][:, 0])),
                                               'X Echo down {:.0f}kcps'.format(np.mean(data['counts'][:, 1])),
                                               'Y Echo up {:.0f}kcps'.format(np.mean(data['counts'][:, 2])),
                                               'Y Echo down {:.0f}kcps'.format(np.mean(data['counts'][:, 3]))),
                                       fontsize=10)
                    axislist[0].set_title(
                        '(final plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
                            'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s)\n after {:d} averages\n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz'.format(
                            self.settings['decoupling_seq']['type'],
                            self.settings['decoupling_seq']['num_of_pulse_blocks'],
                            avg_block_number, self.settings['mw_pulses']['mw_power'],
                            self.settings['mw_pulses']['mw_frequency'] * 1e-9))

                else:
                    axislist[0].legend(labels=('X Echo up {:.0f}kcps'.format(np.mean(data['counts'][:, 0])),
                                               'X Echo down {:.0f}kcps'.format(np.mean(data['counts'][:, 1]))),
                                       fontsize=10)
                    axislist[0].set_title(
                        '(final plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
                            'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s)\n after {:d} averages\n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz'.format(
                            self.settings['decoupling_seq']['type'], self.settings['decoupling_seq']['num_of_pulse_blocks'],
                            avg_block_number, self.settings['mw_pulses']['mw_power'],
                            self.settings['mw_pulses']['mw_frequency'] * 1e-9))
                axislist[0].set_ylabel('fluorescence [kcps]')
                axislist[0].set_xlabel('tau [ns]')

        else:
            super(PDD_XYreadout, self)._plot(axislist)

            if self.settings['do_Y_readout']:
                axislist[0].legend(labels=('X Echo up {:.0f}kcps'.format(np.mean(data['counts'][:, 0])),
                                           'X Echo down {:.0f}kcps'.format(np.mean(data['counts'][:, 1])),
                                           'Y Echo up {:.0f}kcps'.format(np.mean(data['counts'][:, 2])),
                                           'Y Echo down {:.0f}kcps'.format(np.mean(data['counts'][:, 3]))),
                                   fontsize=10)
                axislist[0].set_title(
                    '(initial plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
                        'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s)\n after {:d} averages\n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz'.format(
                        self.settings['decoupling_seq']['type'], self.settings['decoupling_seq']['num_of_pulse_blocks'],
                        avg_block_number, self.settings['mw_pulses']['mw_power'],
                        self.settings['mw_pulses']['mw_frequency'] * 1e-9))

            else:
                axislist[0].legend(labels=('X Echo up {:.0f}kcps'.format(np.mean(data['counts'][:, 0])),
                                           'X Echo down {:.0f}kcps'.format(np.mean(data['counts'][:, 1]))), fontsize=10)
                axislist[0].set_title(
                    '(initial plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
                        'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s) \n after {:d} averages\n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz'.format(
                        self.settings['decoupling_seq']['type'], self.settings['decoupling_seq']['num_of_pulse_blocks'],
                        avg_block_number, self.settings['mw_pulses']['mw_power'],
                        self.settings['mw_pulses']['mw_frequency'] * 1e-9))

    def _update_plot(self, axislist):

        if len(axislist[0].lines) == 0:
            self._plot(axislist)
            return

        super(PDD_XYreadout, self)._update_plot(axislist)

        if self.settings['do_Y_readout']:

            axislist[0].legend(labels=('X Echo up {:.0f}kcps'.format(np.mean(self.data['counts'][:, 0])),
                                       'X Echo down {:.0f}kcps'.format(np.mean(self.data['counts'][:, 1])),
                                       'Y Echo up {:.0f}kcps'.format(np.mean(self.data['counts'][:, 2])),
                                       'Y Echo down {:.0f}kcps'.format(np.mean(self.data['counts'][:, 3]))),
                               fontsize=10)

            axislist[0].set_title(
                '(updating plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
                    'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s)\n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz'.format(
                    self.settings['decoupling_seq']['type'], self.settings['decoupling_seq']['num_of_pulse_blocks'],
                    self.settings['mw_pulses']['mw_power'],
                    self.settings['mw_pulses']['mw_frequency'] * 1e-9))

        else:
            axislist[0].legend(labels=('X Echo up {:.0f}kcps'.format(np.mean(self.data['counts'][:, 0])),
                                       'X Echo down {:.0f}kcps'.format(np.mean(self.data['counts'][:, 1]))), fontsize=10)

            axislist[0].set_title('(updating plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
                'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s) \n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz'.format(
                self.settings['decoupling_seq']['type'], self.settings['decoupling_seq']['num_of_pulse_blocks'],
                self.settings['mw_pulses']['mw_power'],
                self.settings['mw_pulses']['mw_frequency'] * 1e-9))

# class eSensing_swpV(PulsedExperimentBaseScript):
#     """
#     This script does AC electrometry using NV. Tau is fixed, and voltage is swept.
#     Voltages are applied using Agilent 33120A.
#     MW pulses are applied using Keysight N9310A.
#
#     --> Ziwei Qiu 6/26/2019
#
#     """
#     _DEFAULT_SETTINGS = [
#         Parameter('mw_pulses', [
#             Parameter('mw_power', -10.0, float, 'microwave power in dB'),
#             Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
#             Parameter('microwave_channel', 'i', ['i', 'q'], 'Channel to use for the first pi/2 mw pulse'),
#             Parameter('pi_pulse_time', 50.0, float, 'time duration of a pi pulse (in ns) for channel i'),
#             Parameter('pi_half_pulse_time', 25.0, float, 'time duration of a pi/2 pulse (in ns) for channel i'),
#             Parameter('3pi_half_pulse_time', 75.0, float, 'time duration of a 3pi/2 pulse (in ns) for channel i'),
#             Parameter('Q_same_as_I', True, bool, 'check if I and Q are symmetric'),
#             Parameter('pi_pulse_time_q', 50.0, float, 'time duration of a pi pulse (in ns) for channel q'),
#             Parameter('pi_half_pulse_time_q', 25.0, float, 'time duration of a pi/2 pulse (in ns) for channel q'),
#             Parameter('3pi_half_pulse_time_q', 75.0, float, 'time duration of a 3pi/2 pulse (in ns) for channel q')
#         ]),
#         Parameter('e_fields', [
#             Parameter('min_amp', 0.1, float, '[Vpp] minimum amplitude'),
#             Parameter('max_amp', 20.0, float, '[Vpp] maximum amplitude'),
#             Parameter('amp_step', 1., [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0],
#                       '[Vpp] amplitude step increment'),
#             Parameter('offset', 0.0, float, '[Vdc] DC offset voltage'),
#             Parameter('wave_shape', 'SQUare', ['SINusoid', 'SQUare'], 'choose the wave shape'),
#             Parameter('trigger_latency', 1253.0, float, '[ns] Agilent33120A  trigger latency'),
#             Parameter('trigger_time', 500.0, float, '[ns] trigger pulse time (>100ns)'),
#         ]),
#         Parameter('tau_time', 7000, float, '[ns] choose the fixed tau time between two pi/2 pulses'),
#         Parameter('decoupling_seq', [
#             Parameter('type', 'XY4', ['spin_echo', 'CPMG', 'XY4', 'XY8'], 'type of dynamical decoupling sequences'),
#             Parameter('num_of_pulse_blocks', 1, int, 'number of pulse blocks.'),
#             Parameter('eSensing_type', 'sine', ['sine', 'cosine'], 'choose to do sine or cosine electrometry'),
#         ]),
#         # Parameter('read_out', [
#         #     Parameter('meas_time', 200, float, 'measurement time after rabi sequence (in ns)'),
#         #     Parameter('nv_reset_time', 2000, int, 'time with laser on to reset state'),
#         #     Parameter('laser_off_time', 1000, int,
#         #               'minimum laser off time before taking measurements (ns)'),
#         #     Parameter('delay_mw_readout', 600, int, 'delay between mw and readout (in ns)'),
#         #     Parameter('delay_readout', 550, int, 'delay between laser on and readout (given by spontaneous decay rate)')
#         # ]),
#         Parameter('min_pulse_dur', 15, [10, 15, 50], 'Minimum allowed pulse duration (ns)'),
#         Parameter('num_averages', 500000, int, 'number of averages (>100000)'),
#     ]
#
#     _INSTRUMENTS = {'daq': NI6353, 'PB': LISE607RTPulseBlaster, 'mw_gen_a': AgilentMicrowaveGenerator, 'fg': Agilent33120A}
#
#     def _function(self):
#
#         self.data['fits'] = None
#
#         ## Turn off green light (the pulse blaster will pulse it on when needed)
#         self.instruments['PB']['instance'].update({'laser': {'status': False}})
#         self.instruments['PB']['instance'].update({'apd_switch': {'status': False}})
#
#         # set up the RF generator
#         self.instruments['mw_gen_a']['instance'].update({'enable_IQ': True})
#         self.instruments['mw_gen_a']['instance'].update({'enable_modulation': True})
#         self.instruments['mw_gen_a']['instance'].update({'freq_mode': 'CW'})
#         self.instruments['mw_gen_a']['instance'].update({'power_mode': 'CW'})
#         self.instruments['mw_gen_a']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power']})
#         self.instruments['mw_gen_a']['instance'].update({'frequency': self.settings['mw_pulses']['mw_frequency']})
#         self.instruments['mw_gen_a']['instance'].update({'enable_output': True})
#
#         # set up the function generator
#         self.instruments['fg']['instance'].update({'display_on': False})
#         self.instruments['fg']['instance'].update({'burst_mod': True})
#         # self.instruments['fg']['instance'].update({'amplitude': self.settings['e_fields']['min_amp']})
#         self.instruments['fg']['instance'].update({'offset': self.settings['e_fields']['offset']})
#         self.instruments['fg']['instance'].update({'wave_shape': self.settings['e_fields']['wave_shape']})
#         valid = True
#         if self.settings['decoupling_seq']['type'] == 'spin_echo' or self.settings['decoupling_seq']['type'] == 'CPMG':
#             tau = self.settings['tau_time'] / (1 * self.settings['decoupling_seq']['num_of_pulse_blocks'])
#             if self.settings['decoupling_seq']['num_of_pulse_blocks'] == 1:
#                 fg_freq = 1E9 / (1 * tau)
#                 burst_cnt = 1
#                 burst_phase = 0.0
#                 if (fg_freq / 1E6) / burst_cnt > 1:
#                     self.log('Invalid setting: (fg_freq / 1E6) / burst_cnt > 1.')
#                     valid = False
#                 else:
#                     self.instruments['fg']['instance'].update({'burst_phase': burst_phase})
#                     self.instruments['fg']['instance'].update({'frequency': fg_freq})
#                     self.instruments['fg']['instance'].update({'burst_count': burst_cnt})
#             elif self.settings['decoupling_seq']['num_of_pulse_blocks']%2 == 0:
#                 fg_freq = 1E9 / (2 * tau)
#                 burst_cnt = self.settings['decoupling_seq']['num_of_pulse_blocks'] / 2
#                 burst_phase = 90.0
#                 if (fg_freq / 1E6) / burst_cnt > 1:
#                     self.log('Invalid setting: (fg_freq / 1E6) / burst_cnt > 1.')
#                     valid = False
#                 else:
#                     self.instruments['fg']['instance'].update({'burst_phase': burst_phase})
#                     self.instruments['fg']['instance'].update({'frequency': fg_freq})
#                     self.instruments['fg']['instance'].update({'burst_count': burst_cnt})
#             else:
#                 self.log('Invalid setting: num_of_pulse_blocks needs to be 1 or even number if using spin-echo or CPMG.')
#                 valid = False
#         elif self.settings['decoupling_seq']['type'] == 'XY4':
#             tau = self.settings['tau_time'] / (4 * self.settings['decoupling_seq']['num_of_pulse_blocks'])
#             fg_freq = 1E9 / (2 * tau)
#             burst_cnt = 2 * self.settings['decoupling_seq']['num_of_pulse_blocks']
#             burst_phase = 90.0
#             if (fg_freq / 1E6) / burst_cnt > 1:
#                 self.log('Invalid setting: (fg_freq / 1E6) / burst_cnt > 1.')
#                 valid = False
#             else:
#                 self.instruments['fg']['instance'].update({'burst_phase': burst_phase})
#                 self.instruments['fg']['instance'].update({'frequency': fg_freq})
#                 self.instruments['fg']['instance'].update({'burst_count': burst_cnt})
#         elif self.settings['decoupling_seq']['type'] == 'XY8':
#             tau = self.settings['tau_time'] / (8 * self.settings['decoupling_seq']['num_of_pulse_blocks'])
#             fg_freq = 1E9 / (2 * tau)
#             burst_cnt = 4 * self.settings['decoupling_seq']['num_of_pulse_blocks']
#             burst_phase = 90.0
#             if (fg_freq / 1E6) / burst_cnt > 1:
#                 self.log('Invalid setting: (fg_freq / 1E6) / burst_cnt > 1.')
#                 valid = False
#             else:
#                 self.instruments['fg']['instance'].update({'burst_phase': burst_phase})
#                 self.instruments['fg']['instance'].update({'frequency': fg_freq})
#                 self.instruments['fg']['instance'].update({'burst_count': burst_cnt})
#         else:
#             self.log('Invalid setting: choose the right decoupling_seq type')
#             valid = False
#
#         if valid:
#             self.avg_block_number = super(eSensing_swpV, self)._function()
#
#             # this is to flag when the experiment is done or stopped and this allows proper saving even fit fails
#             self.data['exp_finished'] = 1
#             self.data['avg_block_number'] = self.avg_block_number
#
#         # turn off laser, apd_switch and MW (just in case)
#         self.instruments['PB']['instance'].update({'laser': {'status': False}})
#         self.instruments['PB']['instance'].update({'apd_switch': {'status': False}})
#         self.instruments['mw_gen_a']['instance'].update({'enable_output': False})
#         self.instruments['mw_gen_a']['instance'].update({'enable_IQ': False})
#         self.instruments['mw_gen_a']['instance'].update({'enable_modulation': False})
#         self.instruments['fg']['instance'].update({'display_on': True})
#         self.instruments['fg']['instance'].update({'amplitude': 0.1})
#         self.instruments['fg']['instance'].update({'offset': 0.0})
#
#         # no fitting for now
#         if 'counts' in self.data.keys() and 'tau' in self.data.keys() and not 0 in (
#                 self.data['counts'][:, 1] + self.data['counts'][:, 0]) and not 0 in (
#                 self.data['counts'][:, 3] + self.data['counts'][:, 2]):
#
#             self.data['norm_echo'] = 2. * (- self.data['counts'][:, 1] + self.data['counts'][:, 0]) / (
#                     self.data['counts'][:, 1] + self.data['counts'][:, 0])
#             self.data['norm_esig'] = 2. * (- self.data['counts'][:, 3] + self.data['counts'][:, 2]) / (
#                     self.data['counts'][:, 3] + self.data['counts'][:, 2])
#
#             tau = self.data['tau']
#             counts = self.data['norm_echo']
#             esig_counts = self.data['norm_esig']
#             fit_success = False
#
#             # try:
#             #     fits = fit_exp_decay(tau, counts, offset = True, verbose = True)
#             #     self.data['fits'] = fits
#             #     fit_success = True
#             # except:
#             #     self.data['fits'] = None
#             #     self.log('T2 fit failed')
#             #     fit_success = False
#             #
#             # if self.settings['decoupling_seq']['do_eSensing'] and not 0 in (self.data['counts'][:, 3] + self.data['counts'][:, 2]):
#             #     self.data['norm_esig'] = 2. * (- self.data['counts'][:, 3] + self.data['counts'][:, 2]) / (
#             #                 self.data['counts'][:, 3] + self.data['counts'][:, 2])
#             #
#             #     esig_counts = self.data['norm_esig']
#                 # try:
#                 #     fits_esig = fit_exp_decay(tau, esig_counts, offset = True, verbose = True)
#                 #     self.data['fits_esig'] = fits_esig
#                 #     fit_esig_success = True
#                 # except:
#                 #     self.data['fits_esig'] = None
#                 #     self.log('E-sensing T2 fit failed')
#                 #     fit_esig_success = False
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
#         max_range = int(np.floor((self.settings['e_fields']['max_amp'] - self.settings['e_fields']['min_amp']) /
#                                  self.settings['e_fields']['amp_step']))
#         amp_list = np.array(
#             [self.settings['e_fields']['min_amp'] + i * self.settings['e_fields']['amp_step'] for i in
#              range(max_range)])
#         # ignore the sequence if the amplitude is smaller than 15Vpp
#         amp_list = [x for x in amp_list if x >= 0.1]
#
#         print('Vpp_list:', amp_list)
#
#         nv_reset_time = self.settings['read_out']['nv_reset_time']
#         delay_readout = self.settings['read_out']['delay_readout']
#
#         pi_time = self.settings['mw_pulses']['pi_pulse_time']
#         pi_half_time = self.settings['mw_pulses']['pi_half_pulse_time']
#         three_pi_half_time = self.settings['mw_pulses']['3pi_half_pulse_time']
#         if self.settings['mw_pulses']['Q_same_as_I']:
#             pi_time_q = pi_time
#             pi_half_time_q = pi_half_time
#             three_pi_half_time_q = three_pi_half_time
#         else:
#             pi_time_q = self.settings['mw_pulses']['pi_pulse_time']
#             pi_half_time_q = self.settings['mw_pulses']['pi_half_pulse_time']
#             three_pi_half_time_q = self.settings['mw_pulses']['3pi_half_pulse_time']
#
#         if self.settings['mw_pulses']['microwave_channel'] == 'i':
#             microwave_channel_1 = 'microwave_i'
#             microwave_channel_2 = 'microwave_q'
#             pi_time_1 = pi_time
#             pi_half_time_1 = pi_half_time
#             three_pi_half_time_1 = three_pi_half_time
#             pi_time_2 = pi_time_q
#             pi_half_time_2 = pi_half_time_q
#             three_pi_half_time_2 = three_pi_half_time_q
#         else:
#             microwave_channel_1 = 'microwave_q'
#             microwave_channel_2 = 'microwave_i'
#             pi_time_1 = pi_time_q
#             pi_half_time_1 = pi_half_time_q
#             three_pi_half_time_1 = three_pi_half_time_q
#             pi_time_2 = pi_time
#             pi_half_time_2 = pi_half_time
#             three_pi_half_time_2 = three_pi_half_time
#
#         if self.settings['decoupling_seq']['eSensing_type'] == 'sine':
#             last_mw_channel = microwave_channel_2
#             last_pi_time = pi_time_2
#             last_pi_half_time = pi_half_time_2
#             last_three_pi_half_time = three_pi_half_time_2
#         else:
#             last_mw_channel = microwave_channel_1
#             last_pi_time = pi_time_1
#             last_pi_half_time = pi_half_time_1
#             last_three_pi_half_time = three_pi_half_time_1
#
#
#         laser_off_time = self.settings['read_out']['laser_off_time']
#         meas_time = self.settings['read_out']['meas_time']
#         delay_mw_readout = self.settings['read_out']['delay_mw_readout']
#         trigger_latency = self.settings['e_fields']['trigger_latency']
#         trigger_time = self.settings['e_fields']['trigger_time']
#
#         number_of_pulse_blocks = self.settings['decoupling_seq']['num_of_pulse_blocks']
#         tau_total = self.settings['tau_time']
#
#         if self.settings['decoupling_seq']['type'] == 'spin_echo':
#             for amp in amp_list:
#                 tau = tau_total / (1 * self.settings['decoupling_seq']['num_of_pulse_blocks'])
#                 pulse_sequence = []
#
#                 # ECHO SEQUENCE:
#                 start_of_first_HE = laser_off_time
#
#                 # the first pi/2 pulse
#                 pulse_sequence = [Pulse(microwave_channel_1, start_of_first_HE, pi_half_time_1)]
#                 section_begin_time = start_of_first_HE + pi_half_time_1 - tau / 2  # for the first pulse, only wait tau/2
#                 for i in range(0, number_of_pulse_blocks):
#
#                     pulse_sequence.append(
#                         Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time_1 / 2, pi_time_1))
#                     section_begin_time += 1 * tau
#
#                 # the second pi/2 pulse and readout
#                 pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + tau / 2, pi_half_time_1))
#                 pulse_sequence.append(Pulse('laser', section_begin_time + tau / 2 + pi_half_time_1 + delay_mw_readout,
#                                              nv_reset_time))
#                 pulse_sequence.append(Pulse('apd_readout',
#                                              section_begin_time + tau / 2 + pi_half_time_1 + delay_mw_readout + delay_readout,
#                                              meas_time))
#                 start_of_second_HE = section_begin_time + tau / 2 + pi_half_time_1 + delay_mw_readout + nv_reset_time + laser_off_time
#
#                 # the first pi/2 pulse
#                 pulse_sequence.append(Pulse(microwave_channel_1, start_of_second_HE, pi_half_time_1))
#                 section_begin_time = start_of_second_HE + pi_half_time_1 - tau / 2  # for the first pulse, only wait tau/2
#                 for i in range(0, number_of_pulse_blocks):
#                     # pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time))
#                     pulse_sequence.append(
#                         Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time_1 / 2, pi_time_1))
#                     section_begin_time += 1 * tau
#
#                 # the second 3*pi/2 pulse and readout
#                 pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + tau / 2, three_pi_half_time_1))
#                 pulse_sequence.append(Pulse('laser',
#                                              section_begin_time + tau / 2 + three_pi_half_time_1 + delay_mw_readout,
#                                              nv_reset_time))
#                 pulse_sequence.append(Pulse('apd_readout',
#                                              section_begin_time + tau / 2 + three_pi_half_time_1 + delay_mw_readout + delay_readout,
#                                              meas_time))
#
#                 # E-Sensing SEQUENCE
#                 if True:
#                     start_of_DEER = section_begin_time + tau / 2 + three_pi_half_time_1 + delay_mw_readout + nv_reset_time + laser_off_time
#                     # the first pi/2 pulse
#                     pulse_sequence.append(Pulse(microwave_channel_1, start_of_DEER, pi_half_time_1))
#                     pulse_sequence.append(Pulse('fg_trigger', start_of_DEER - trigger_latency, trigger_time))
#                     section_begin_time = start_of_DEER + pi_half_time_1 - tau / 2  # for the first pulse, only wait tau/2
#                     for i in range(0, number_of_pulse_blocks):
#                         pulse_sequence.append(
#                             Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time_1 / 2, pi_time_1))
#                         section_begin_time += 1 * tau
#
#                     # the second pi/2 pulse and readout
#                     pulse_sequence.append(Pulse(last_mw_channel, section_begin_time + tau / 2, last_pi_half_time))
#                     pulse_sequence.append(Pulse('laser', section_begin_time + tau / 2 + last_pi_half_time + delay_mw_readout,
#                                                 nv_reset_time))
#                     pulse_sequence.append(Pulse('apd_readout',
#                                                 section_begin_time + tau / 2 + last_pi_half_time + delay_mw_readout + delay_readout,
#                                                 meas_time))
#
#                     start_of_second_DEER = section_begin_time + tau / 2 + last_pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time
#
#                     # the first pi/2 pulse
#                     pulse_sequence.append(Pulse(microwave_channel_1, start_of_second_DEER, pi_half_time_1))
#                     pulse_sequence.append(Pulse('fg_trigger', start_of_second_DEER - trigger_latency, trigger_time))
#                     section_begin_time = start_of_second_DEER + pi_half_time_1 - tau / 2  # for the first pulse, only wait tau/2
#
#                     for i in range(0, number_of_pulse_blocks):
#                         pulse_sequence.append(
#                             Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time_1 / 2, pi_time_1))
#                         section_begin_time += 1 * tau
#
#                     # the second pi/2 pulse and readout
#                     pulse_sequence.append(Pulse(last_mw_channel, section_begin_time + tau / 2, last_three_pi_half_time))
#                     pulse_sequence.append(Pulse('laser', section_begin_time + tau / 2 + last_three_pi_half_time + delay_mw_readout,
#                                                 nv_reset_time))
#                     pulse_sequence.append(Pulse('apd_readout',
#                                                 section_begin_time + tau / 2 + last_three_pi_half_time + delay_mw_readout + delay_readout,
#                                                 meas_time))
#
#                 # Turn on APD switch all the time
#                 pulse_sequence.append(Pulse('apd_switch', 0,
#                                             section_begin_time + tau / 2 + last_three_pi_half_time + delay_mw_readout + nv_reset_time))
#
#                 # print(pulse_sequence)
#                 pulse_sequences.append(pulse_sequence)
#
#             print('number of sequences before validation ', len(pulse_sequences))
#             return pulse_sequences, amp_list, meas_time
#
#         elif self.settings['decoupling_seq']['type'] == 'XY4':
#
#             for amp in amp_list:
#                 tau = tau_total / (4 * self.settings['decoupling_seq']['num_of_pulse_blocks'])
#                 pulse_sequence = []
#
#                 # ECHO SEQUENCE:
#                 start_of_first_HE = laser_off_time
#                 # the first pi/2 pulse
#                 pulse_sequence = [Pulse(microwave_channel_1, start_of_first_HE, pi_half_time_1)]
#                 section_begin_time = start_of_first_HE + pi_half_time_1 - tau / 2  # for the first pulse, only wait tau/2
#                 for i in range(0, number_of_pulse_blocks):
#                     pulse_sequence.append(
#                         Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time_2 / 2, pi_time_2))
#                     pulse_sequence.append(
#                         Pulse(microwave_channel_1, section_begin_time + 2 * tau - pi_time_1 / 2, pi_time_1))
#                     pulse_sequence.append(
#                         Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time_1 / 2, pi_time_1))
#                     pulse_sequence.append(
#                         Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time_2 / 2, pi_time_2))
#
#                     section_begin_time += 4 * tau
#
#                 # the second 3*pi/2 pulse and readout
#                 pulse_sequence.append(
#                     Pulse(microwave_channel_1, section_begin_time + tau / 2, three_pi_half_time_1))
#                 pulse_sequence.append(
#                     Pulse('laser', section_begin_time + tau / 2 + three_pi_half_time_1 + delay_mw_readout,
#                           nv_reset_time))
#                 pulse_sequence.append(
#                     Pulse('apd_readout',
#                           section_begin_time + tau / 2 + three_pi_half_time_1 + delay_mw_readout + delay_readout,
#                           meas_time))
#                 start_of_second_HE = section_begin_time + tau / 2 + three_pi_half_time_1 + delay_mw_readout + nv_reset_time + laser_off_time
#
#                 # the first pi/2 pulse
#                 pulse_sequence.append(
#                     Pulse(microwave_channel_1, start_of_second_HE, pi_half_time_1))
#                 section_begin_time = start_of_second_HE + pi_half_time_1 - tau / 2  # for the first pulse, only wait tau/2
#                 for i in range(0, number_of_pulse_blocks):
#                     pulse_sequence.append(
#                         Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time_2 / 2, pi_time_2))
#                     pulse_sequence.append(
#                         Pulse(microwave_channel_1, section_begin_time + 2 * tau - pi_time_1 / 2, pi_time_1))
#                     pulse_sequence.append(
#                         Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time_1 / 2, pi_time_1))
#                     pulse_sequence.append(
#                         Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time_2 / 2, pi_time_2))
#                     section_begin_time += 4 * tau
#
#                 # the second pi/2 pulse and readout
#                 pulse_sequence.append(
#                     Pulse(microwave_channel_1, section_begin_time + tau / 2, pi_half_time_1))
#                 pulse_sequence.append(
#                     Pulse('laser',
#                           section_begin_time + tau / 2 + pi_half_time_1 + delay_mw_readout,
#                           nv_reset_time))
#                 pulse_sequence.append(
#                     Pulse('apd_readout',
#                           section_begin_time + tau / 2 + pi_half_time_1 + delay_mw_readout + delay_readout,
#                           meas_time))
#
#                 # E-Sensing SEQUENCE
#                 if True:
#                     start_of_DEER = section_begin_time + tau / 2 + pi_half_time_1 + delay_mw_readout + nv_reset_time + laser_off_time
#                     # the first pi/2 pulse
#                     pulse_sequence.append(Pulse(microwave_channel_1, start_of_DEER, pi_half_time_1))
#                     pulse_sequence.append(Pulse('fg_trigger', start_of_DEER - trigger_latency, trigger_time))
#                     section_begin_time = start_of_DEER + pi_half_time_1 - tau / 2  # for the first pulse, only wait tau/2
#                     for i in range(0, number_of_pulse_blocks):
#                         pulse_sequence.append(
#                             Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time_2 / 2, pi_time_2))
#                         pulse_sequence.append(
#                             Pulse(microwave_channel_1, section_begin_time + 2 * tau - pi_time_1 / 2, pi_time_1))
#                         pulse_sequence.append(
#                             Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time_1 / 2, pi_time_1))
#                         pulse_sequence.append(
#                             Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time_2 / 2, pi_time_2))
#
#                         section_begin_time += 4 * tau
#
#                     # the second 3*pi/2 pulse and readout
#                     pulse_sequence.append(
#                         Pulse(last_mw_channel, section_begin_time + tau / 2, last_three_pi_half_time))
#                     pulse_sequence.append(
#                         Pulse('laser', section_begin_time + tau / 2 + last_three_pi_half_time + delay_mw_readout,
#                               nv_reset_time))
#                     pulse_sequence.append(
#                         Pulse('apd_readout',
#                               section_begin_time + tau / 2 + last_three_pi_half_time + delay_mw_readout + delay_readout,
#                               meas_time))
#
#                     start_of_second_DEER = section_begin_time + tau / 2 + last_three_pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time
#
#                     # the first pi/2 pulse
#                     pulse_sequence.append(Pulse(microwave_channel_1, start_of_second_DEER, pi_half_time_1))
#                     pulse_sequence.append(Pulse('fg_trigger', start_of_second_DEER - trigger_latency, trigger_time))
#                     section_begin_time = start_of_second_DEER + pi_half_time_1 - tau / 2  # for the first pulse, only wait tau/2
#
#                     for i in range(0, number_of_pulse_blocks):
#                         pulse_sequence.append(
#                             Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time_2 / 2, pi_time_2))
#                         pulse_sequence.append(
#                             Pulse(microwave_channel_1, section_begin_time + 2 * tau - pi_time_1 / 2, pi_time_1))
#                         pulse_sequence.append(
#                             Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time_1 / 2, pi_time_1))
#                         pulse_sequence.append(
#                             Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time_2 / 2, pi_time_2))
#                         section_begin_time += 4 * tau
#
#                     # the second pi/2 pulse and readout
#                     pulse_sequence.append(
#                         Pulse(last_mw_channel, section_begin_time + tau / 2, last_pi_half_time))
#                     pulse_sequence.append(
#                         Pulse('laser',
#                               section_begin_time + tau / 2 + last_pi_half_time + delay_mw_readout,
#                               nv_reset_time))
#                     pulse_sequence.append(
#                         Pulse('apd_readout',
#                               section_begin_time + tau / 2 + last_pi_half_time + delay_mw_readout + delay_readout,
#                               meas_time))
#
#                 pulse_sequence.append(Pulse('apd_switch', 0,
#                                             section_begin_time + tau / 2 + last_pi_half_time + delay_mw_readout + nv_reset_time))
#                 pulse_sequences.append(pulse_sequence)
#
#             print('number of sequences before validation ', len(pulse_sequences))
#             # return pulse_sequences, self.settings['num_averages'], tau_list, meas_time
#             return pulse_sequences, amp_list, meas_time
#
#         elif self.settings['decoupling_seq']['type'] == 'XY8':
#
#             for amp in amp_list:
#                 tau = tau_total / (8 * self.settings['decoupling_seq']['num_of_pulse_blocks'])
#                 pulse_sequence = []
#
#                 # ECHO SEQUENCE:
#                 start_of_first_HE = laser_off_time
#
#                 # the first pi/2 pulse
#                 pulse_sequence.append(
#                     Pulse(microwave_channel_1, start_of_first_HE, pi_half_time_1))
#
#                 section_begin_time = start_of_first_HE + pi_half_time_1 - tau / 2  # for the first pulse, only wait tau/2
#                 for i in range(0, number_of_pulse_blocks):
#
#                     pulse_sequence.append(
#                         Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time_1 / 2, pi_time_1))
#                     pulse_sequence.append(
#                         Pulse(microwave_channel_2, section_begin_time + 2 * tau - pi_time_2 / 2, pi_time_2))
#                     pulse_sequence.append(
#                         Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time_1 / 2, pi_time_1))
#                     pulse_sequence.append(
#                         Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time_2 / 2, pi_time_2))
#                     pulse_sequence.append(
#                         Pulse(microwave_channel_2, section_begin_time + 5 * tau - pi_time_2 / 2, pi_time_2))
#                     pulse_sequence.append(
#                         Pulse(microwave_channel_1, section_begin_time + 6 * tau - pi_time_1 / 2, pi_time_1))
#                     pulse_sequence.append(
#                         Pulse(microwave_channel_2, section_begin_time + 7 * tau - pi_time_2 / 2, pi_time_2))
#                     pulse_sequence.append(
#                         Pulse(microwave_channel_1, section_begin_time + 8 * tau - pi_time_1 / 2, pi_time_1))
#                     section_begin_time += 8 * tau
#
#                 # the second 3*pi/2 pulse and readout
#                 pulse_sequence.append(
#                     Pulse(microwave_channel_1, section_begin_time + tau / 2, three_pi_half_time_1))
#                 pulse_sequence.append(
#                     Pulse('laser', section_begin_time + tau / 2 + three_pi_half_time_1 + delay_mw_readout,
#                           nv_reset_time))
#                 pulse_sequence.append(
#                     Pulse('apd_readout',
#                           section_begin_time + tau / 2 + three_pi_half_time_1 + delay_mw_readout + delay_readout,
#                           meas_time))
#
#                 start_of_second_HE = section_begin_time + tau / 2 + three_pi_half_time_1 + delay_mw_readout + nv_reset_time + laser_off_time
#
#                 # the first pi/2 pulse
#                 pulse_sequence.append(
#                     Pulse(microwave_channel_1, start_of_second_HE, pi_half_time_1))
#                 section_begin_time = start_of_second_HE + pi_half_time_1 - tau / 2  # for the first pulse, only wait tau/2
#                 for i in range(0, number_of_pulse_blocks):
#
#                     pulse_sequence.append(
#                         Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time_1 / 2, pi_time_1))
#                     pulse_sequence.append(
#                         Pulse(microwave_channel_2, section_begin_time + 2 * tau - pi_time_2 / 2, pi_time_2))
#                     pulse_sequence.append(
#                         Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time_1 / 2, pi_time_1))
#                     pulse_sequence.append(
#                         Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time_2 / 2, pi_time_2))
#                     pulse_sequence.append(
#                         Pulse(microwave_channel_2, section_begin_time + 5 * tau - pi_time_2 / 2, pi_time_2))
#                     pulse_sequence.append(
#                         Pulse(microwave_channel_1, section_begin_time + 6 * tau - pi_time_1 / 2, pi_time_1))
#                     pulse_sequence.append(
#                         Pulse(microwave_channel_2, section_begin_time + 7 * tau - pi_time_2 / 2, pi_time_2))
#                     pulse_sequence.append(
#                         Pulse(microwave_channel_1, section_begin_time + 8 * tau - pi_time_1 / 2, pi_time_1))
#                     section_begin_time += 8 * tau
#
#                 # the second pi/2 pulse and readout
#                 pulse_sequence.append(
#                     Pulse(microwave_channel_1, section_begin_time + tau / 2, pi_half_time_1))
#                 pulse_sequence.append(
#                     Pulse('laser',
#                           section_begin_time + tau / 2 + pi_half_time_1 + delay_mw_readout,
#                           nv_reset_time))
#                 pulse_sequence.append(
#                     Pulse('apd_readout',
#                           section_begin_time + tau / 2 + pi_half_time_1 + delay_mw_readout + delay_readout,
#                           meas_time))
#
#                 # E-Sensing SEQUENCE
#                 if True:
#                     start_of_DEER = section_begin_time + tau / 2 + pi_half_time_1 + delay_mw_readout + nv_reset_time + laser_off_time
#                     # the first pi/2 pulse
#                     pulse_sequence.append(Pulse(microwave_channel_1, start_of_DEER, pi_half_time_1))
#                     pulse_sequence.append(Pulse('fg_trigger', start_of_DEER - trigger_latency, trigger_time))
#
#                     section_begin_time = start_of_DEER + pi_half_time_1 - tau / 2  # for the first pulse, only wait tau/2
#                     for i in range(0, number_of_pulse_blocks):
#                         pulse_sequence.append(
#                             Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time_1 / 2, pi_time_1))
#                         pulse_sequence.append(
#                             Pulse(microwave_channel_2, section_begin_time + 2 * tau - pi_time_2 / 2, pi_time_2))
#                         pulse_sequence.append(
#                             Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time_1 / 2, pi_time_1))
#                         pulse_sequence.append(
#                             Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time_2 / 2, pi_time_2))
#                         pulse_sequence.append(
#                             Pulse(microwave_channel_2, section_begin_time + 5 * tau - pi_time_2 / 2, pi_time_2))
#                         pulse_sequence.append(
#                             Pulse(microwave_channel_1, section_begin_time + 6 * tau - pi_time_1 / 2, pi_time_1))
#                         pulse_sequence.append(
#                             Pulse(microwave_channel_2, section_begin_time + 7 * tau - pi_time_2 / 2, pi_time_2))
#                         pulse_sequence.append(
#                             Pulse(microwave_channel_1, section_begin_time + 8 * tau - pi_time_1 / 2, pi_time_1))
#                         section_begin_time += 8 * tau
#
#                     # the second 3*pi/2 pulse and readout
#                     pulse_sequence.append(
#                         Pulse(last_mw_channel, section_begin_time + tau / 2, last_three_pi_half_time))
#                     pulse_sequence.append(
#                         Pulse('laser', section_begin_time + tau / 2 + last_three_pi_half_time + delay_mw_readout,
#                               nv_reset_time))
#                     pulse_sequence.append(
#                         Pulse('apd_readout',
#                               section_begin_time + tau / 2 + last_three_pi_half_time + delay_mw_readout + delay_readout,
#                               meas_time))
#
#                     start_of_second_DEER = section_begin_time + tau / 2 + last_three_pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time
#
#                     # the first pi/2 pulse
#                     pulse_sequence.append(Pulse(microwave_channel_1, start_of_second_DEER, pi_half_time_1))
#                     pulse_sequence.append(Pulse('fg_trigger', start_of_second_DEER - trigger_latency, trigger_time))
#                     section_begin_time = start_of_second_DEER + pi_half_time_1 - tau / 2  # for the first pulse, only wait tau/2
#
#                     for i in range(0, number_of_pulse_blocks):
#                         pulse_sequence.append(
#                             Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time_1 / 2, pi_time_1))
#                         pulse_sequence.append(
#                             Pulse(microwave_channel_2, section_begin_time + 2 * tau - pi_time_2 / 2, pi_time_2))
#                         pulse_sequence.append(
#                             Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time_1 / 2, pi_time_1))
#                         pulse_sequence.append(
#                             Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time_2 / 2, pi_time_2))
#                         pulse_sequence.append(
#                             Pulse(microwave_channel_2, section_begin_time + 5 * tau - pi_time_2 / 2, pi_time_2))
#                         pulse_sequence.append(
#                             Pulse(microwave_channel_1, section_begin_time + 6 * tau - pi_time_1 / 2, pi_time_1))
#                         pulse_sequence.append(
#                             Pulse(microwave_channel_2, section_begin_time + 7 * tau - pi_time_2 / 2, pi_time_2))
#                         pulse_sequence.append(
#                             Pulse(microwave_channel_1, section_begin_time + 8 * tau - pi_time_1 / 2, pi_time_1))
#                         section_begin_time += 8 * tau
#
#                     # the second pi/2 pulse and readout
#                     pulse_sequence.append(
#                         Pulse(last_mw_channel, section_begin_time + tau / 2, last_pi_half_time))
#                     pulse_sequence.append(
#                         Pulse('laser',
#                               section_begin_time + tau / 2 + last_pi_half_time + delay_mw_readout,
#                               nv_reset_time))
#                     pulse_sequence.append(
#                         Pulse('apd_readout',
#                               section_begin_time + tau / 2 + last_pi_half_time + delay_mw_readout + delay_readout,
#                               meas_time))
#
#                 pulse_sequence.append(Pulse('apd_switch', 0,
#                                             section_begin_time + tau / 2 + last_pi_half_time + delay_mw_readout + nv_reset_time))
#                 pulse_sequences.append(pulse_sequence)
#
#             print('number of sequences before validation ', len(pulse_sequences))
#             # return pulse_sequences, self.settings['num_averages'], tau_list, meas_time
#             return pulse_sequences, amp_list, meas_time
#
#         elif self.settings['decoupling_seq']['type'] == 'CPMG':
#
#             for amp in amp_list:
#                 tau = tau_total / (1 * self.settings['decoupling_seq']['num_of_pulse_blocks'])
#                 pulse_sequence = []
#
#                 # ECHO SEQUENCE:
#                 start_of_first_HE = laser_off_time
#
#                 # the first pi/2 pulse
#                 pulse_sequence = [Pulse(microwave_channel_1, start_of_first_HE, pi_half_time_1)]
#                 section_begin_time = start_of_first_HE + pi_half_time_1 - tau / 2  # for the first pulse, only wait tau/2
#
#                 for i in range(0, number_of_pulse_blocks):
#                     pulse_sequence.append(Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time_2 / 2, pi_time_2))
#                     section_begin_time += 1 * tau
#
#                 # the second 3*pi/2 pulse and readout
#                 pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + tau / 2, three_pi_half_time_1))
#                 pulse_sequence.append(Pulse('laser', section_begin_time + tau / 2 + three_pi_half_time_1 + delay_mw_readout,
#                                              nv_reset_time))
#                 pulse_sequence.append(Pulse('apd_readout',
#                                              section_begin_time + tau / 2 + three_pi_half_time_1 + delay_mw_readout + delay_readout,
#                                              meas_time))
#
#                 start_of_second_HE = section_begin_time + tau / 2 + three_pi_half_time_1 + delay_mw_readout + nv_reset_time + laser_off_time
#
#                 # the first pi/2 pulse
#                 pulse_sequence.append(Pulse(microwave_channel_1, start_of_second_HE, pi_half_time_1))
#
#                 section_begin_time = start_of_second_HE + pi_half_time_1 - tau / 2  # for the first pulse, only wait tau/2
#
#                 for i in range(0, number_of_pulse_blocks):
#                     pulse_sequence.append(Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time_2 / 2, pi_time_2))
#                     section_begin_time += 1 * tau
#
#                 # the second pi/2 pulse and readout
#                 pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + tau / 2, pi_half_time_1))
#                 pulse_sequence.append(Pulse('laser',
#                                              section_begin_time + tau / 2 + pi_half_time_1 + delay_mw_readout,
#                                              nv_reset_time))
#                 pulse_sequence.append(Pulse('apd_readout',
#                                              section_begin_time + tau / 2 + pi_half_time_1 + delay_mw_readout + delay_readout,
#                                              meas_time))
#
#                 # E-Sensing SEQUENCE
#                 if True:
#                     start_of_DEER = section_begin_time + tau / 2 + pi_half_time_1 + delay_mw_readout + nv_reset_time + laser_off_time
#                     # the first pi/2 pulse
#                     pulse_sequence.append(Pulse(microwave_channel_1, start_of_DEER, pi_half_time_1))
#                     pulse_sequence.append(Pulse('fg_trigger', start_of_DEER - trigger_latency, trigger_time))
#                     section_begin_time = start_of_DEER + pi_half_time_1 - tau / 2  # for the first pulse, only wait tau/2
#                     for i in range(0, number_of_pulse_blocks):
#                         pulse_sequence.append(
#                             Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time_2 / 2, pi_time_2))
#                         section_begin_time += 1 * tau
#
#                     # the second 3*pi/2 pulse and readout
#                     pulse_sequence.append(Pulse(last_mw_channel, section_begin_time + tau / 2, last_three_pi_half_time))
#                     pulse_sequence.append(Pulse('laser', section_begin_time + tau / 2 + last_three_pi_half_time + delay_mw_readout,
#                                                 nv_reset_time))
#                     pulse_sequence.append(Pulse('apd_readout',
#                                                 section_begin_time + tau / 2 + last_three_pi_half_time + delay_mw_readout + delay_readout,
#                                                 meas_time))
#
#                     start_of_second_DEER = section_begin_time + tau / 2 + last_three_pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time
#
#                     # the first pi/2 pulse
#                     pulse_sequence.append(Pulse(microwave_channel_1, start_of_second_DEER, pi_half_time_1))
#                     pulse_sequence.append(Pulse('fg_trigger', start_of_second_DEER - trigger_latency, trigger_time))
#                     section_begin_time = start_of_second_DEER + pi_half_time_1 - tau / 2  # for the first pulse, only wait tau/2
#
#                     for i in range(0, number_of_pulse_blocks):
#                         pulse_sequence.append(
#                             Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time_2 / 2, pi_time_2))
#                         section_begin_time += 1 * tau
#
#
#                     # the second pi/2 pulse and readout
#                     pulse_sequence.append(Pulse(last_mw_channel, section_begin_time + tau / 2, last_pi_half_time))
#                     pulse_sequence.append(
#                         Pulse('laser', section_begin_time + tau / 2 + last_pi_half_time + delay_mw_readout,
#                               nv_reset_time))
#                     pulse_sequence.append(Pulse('apd_readout',
#                                                 section_begin_time + tau / 2 + last_pi_half_time + delay_mw_readout + delay_readout,
#                                                 meas_time))
#
#                 # Turn on APD switch all the time
#                 pulse_sequence.append(Pulse('apd_switch', 0,
#                                             section_begin_time + tau / 2 + last_pi_half_time + delay_mw_readout + nv_reset_time))
#
#
#                 pulse_sequences.append(pulse_sequence)
#
#             print('number of sequences before validation ', len(pulse_sequences))
#             # return pulse_sequences, self.settings['num_averages'], tau_list, meas_time
#             return pulse_sequences, amp_list, meas_time
#
#     def _run_sweep(self, pulse_sequences, num_loops_sweep, num_daq_reads, tau_list, verbose=False):
#         """
#         Each pulse sequence specified in pulse_sequences is run num_loops_sweep consecutive times.
#
#         Args:
#             pulse_sequences: a list of pulse sequences to run, each corresponding to a different value of tau. Each
#                              sequence is a list of Pulse objects specifying a given pulse sequence
#             num_loops_sweep: number of times to repeat each sequence before moving on to the next one
#             num_daq_reads: number of times the daq must read for each sequence (generally 1, 2, or 3)
#
#         Poststate: self.data['counts'] is updated with the acquired data
#
#         """
#         # ER 20180731 set init_fluor to zero
#      #   self.data['init_fluor'] = 0.
#
#         # print('start _run_sweep now')
#
#
#         rand_indexes = list(range(len(pulse_sequences)))
#         if self.settings['randomize']:
#             random.shuffle(rand_indexes)
#
#         if verbose:
#             print(('_run_sweep number of pulse sequences', len(pulse_sequences)))
#
#         for index, sequence in enumerate(pulse_sequences):
#             # print('current index')
#             # print(index)
#             if verbose:
#                 print(('_run_sweep index', index, len(pulse_sequences)))
#
#             rand_index = rand_indexes[index]
#             self.instruments['fg']['instance'].update({'amplitude': float(tau_list[rand_index])})
#
#             # self.current_rand_index = rand_index
#             if self._abort:
#                 # print('aborted')
#                 self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})
#                 # print('mw is off')
#                 break
#
#             # print('running single sqeucnce')
#             result = self._run_single_sequence(pulse_sequences[rand_index], num_loops_sweep, num_daq_reads)  # keep entire array
#             # print('single sqeucnce done')
#             self.count_data[rand_index] = self.count_data[rand_index] + result
#             # print('before counts_to_check')
#
#             counts_to_check = self._normalize_to_kCounts(np.array(result), self.measurement_gate_width, num_loops_sweep)
#             self.data['counts'][rand_index] = self._normalize_to_kCounts(self.count_data[rand_index], self.measurement_gate_width,
#                                                                     self.current_averages)
#
#             # self.sequence_index = index
#             self.sequence_index = rand_index
#             counts_temp = counts_to_check[0]
#             # track to the NV if necessary ER 5/31/17
#             if self.settings['Tracking']['on/off']:
#                 print('doing tracking')
#                 print('counts_temp:', counts_temp)
#                 if (1+(1-self.settings['Tracking']['threshold']))*self.settings['Tracking']['init_fluor'] < counts_temp or \
#                         self.settings['Tracking']['threshold']*self.settings['Tracking']['init_fluor'] > counts_temp:
#                     print('(2-threshold)*intial_fluor:',
#                           (1 + (1 - self.settings['Tracking']['threshold'])) * self.settings['Tracking']['init_fluor'])
#                     print('threshold*intial_fluor:', self.settings['Tracking']['threshold']*self.settings['Tracking']['init_fluor'])
#                     if verbose:
#                         print('TRACKING TO NV...')
#                     self.scripts['find_nv'].run()
#                     self.scripts['find_nv'].settings['initial_point'] = self.scripts['find_nv'].data['maximum_point']
#                     print('find_nv done')
#                     if self.scripts['find_nv'].data[
#                         'fluorescence'] == 0.0:  # if it doesn't find an NV, abort the experiment
#                         print('Could not find an NV in FindNV.')
#                         self.log('Could not find an NV in FindNV.')
#                         self._abort = True
#                         return  # exit function in case no NV is found
#                     else:
#                         print('Start optimize_z.')
#                         self.scripts['optimize_z'].run()
#             self.updateProgress.emit(self._calc_progress(index))
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
#         if 'avg_block_number' in data.keys() is not None and data['avg_block_number'] is not None:
#             avg_block_number = data['avg_block_number']
#         else:
#             avg_block_number = 0
#
#         # if 'fits' in data.keys() is not None and data['fits'] is not None:
#         #
#         #     # counts = (-data['counts'][:,1] + data['counts'][:,0])/ (data['counts'][:,0] + data['counts'][:,1])
#         #     tau = data['tau']
#         #     fits = data['fits']
#         #
#         #     axislist[0].plot(tau, data['norm_echo'])
#         #     tauinterp = np.linspace(np.min(tau), np.max(tau), 100)
#         #     axislist[0].plot(tauinterp, exp_offset(tauinterp, fits[0], fits[1], fits[2]))
#         #
#         #
#         #     if 'fits_esig' in data.keys() is not None and data['fits_esig'] is not None:
#         #         fits_esig = data['fits_esig']
#         #         axislist[0].plot(tau, data['norm_esig'])
#         #         axislist[0].plot(tauinterp, exp_offset(tauinterp, fits_esig[0], fits_esig[1], fits_esig[2]))
#         #         axislist[0].legend(labels=('Echo', 'Exp Fit','E-Sensing', 'Exp Fit'), fontsize=10)
#         #         axislist[0].set_title(
#         #             'Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
#         #                 'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s), {:s} electrometry\n averaged over {:d} blocks \n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz, \n T2 (p = 1) = {:2.1f} ns, T2 (E-field) = {:2.1f} ns'.format(
#         #                 self.settings['decoupling_seq']['type'],
#         #                 self.settings['decoupling_seq']['num_of_pulse_blocks'],
#         #                 self.settings['decoupling_seq']['eSensing_type'], avg_block_number,
#         #                 self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency'] * 1e-9,
#         #                 fits[1], fits_esig[1]))
#         #     else:
#         #         axislist[0].legend(labels=('Echo', 'Exp Fit'), fontsize=10)
#         #         axislist[0].set_title(
#         #             'Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
#         #                 'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s)\n averaged over {:d} blocks \n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz, \n T2 (p = 1) = {:2.1f} ns'.format(
#         #                 self.settings['decoupling_seq']['type'],
#         #                 self.settings['decoupling_seq']['num_of_pulse_blocks'], avg_block_number,
#         #                 self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency'] * 1e-9,
#         #                 fits[1]))
#         #     axislist[0].set_ylabel('fluorescence contrast')
#         #     axislist[0].set_xlabel('tau [ns]')
#
#         if 'exp_finished' in data.keys() is not None and data['exp_finished'] == 1:
#             if 'norm_echo' in data.keys() is not None and data['norm_echo'] is not None:
#                 tau = data['tau']
#                 axislist[0].plot(tau, data['norm_echo'])
#                 if 'norm_esig' in data.keys() is not None and data['norm_esig'] is not None:
#                     axislist[0].plot(tau, data['norm_esig'])
#                     axislist[0].legend(labels=('PDD', 'E-Sensing'), fontsize=10)
#                     axislist[0].set_title(
#                         '(final plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
#                             'microwave_channel'] + '\nPDD: {:s} {:d} block(s), {:s} electrometry\n after {:d} averages \n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz \n tau = {:.2f} us'.format(
#                             self.settings['decoupling_seq']['type'],
#                             self.settings['decoupling_seq']['num_of_pulse_blocks'],
#                             self.settings['decoupling_seq']['eSensing_type'], avg_block_number,
#                             self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency'] * 1e-9, self.settings['tau_time'] * 1e-3))
#                 else:
#                     axislist[0].legend(labels=('PDD'), fontsize=10)
#                     axislist[0].set_title(
#                         '(final plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
#                             'microwave_channel'] + '\nPDD: {:s} {:d} block(s)\n after {:d} averages \n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz \n tau = {:.2f} us'.format(
#                             self.settings['decoupling_seq']['type'],
#                             self.settings['decoupling_seq']['num_of_pulse_blocks'],
#                             avg_block_number,
#                             self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency'] * 1e-9, self.settings['tau_time'] * 1e-3))
#                 axislist[0].set_ylabel('fluorescence contrast')
#                 axislist[0].set_xlabel('amplitude [Vpp]')
#             else:
#                 axislist[0].plot(data['tau'], data['counts'][:, 0])
#                 axislist[0].plot(data['tau'], data['counts'][:, 1])
#                 axislist[0].plot(data['tau'], data['counts'][:, 2])
#                 axislist[0].plot(data['tau'], data['counts'][:, 3])
#                 axislist[0].legend(labels=('PDD up {:.0f}kcps'.format(np.mean(data['counts'][:, 0])),
#                                            'PDD down {:.0f}kcps'.format(np.mean(data['counts'][:, 1])),
#                                            'E-Sensing up {:.0f}kcps'.format(np.mean(data['counts'][:, 2])),
#                                            'E-Sensing down {:.0f}kcps'.format(np.mean(data['counts'][:, 3]))),
#                                    fontsize=10)
#                 axislist[0].set_title(
#                     '(final plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
#                         'microwave_channel'] + '\nPDD: {:s} {:d} block(s), {:s} electrometry\n after {:d} averages\n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz \n tau = {:.2f} us'.format(
#                         self.settings['decoupling_seq']['type'],
#                         self.settings['decoupling_seq']['num_of_pulse_blocks'],
#                         self.settings['decoupling_seq']['eSensing_type'],
#                         avg_block_number, self.settings['mw_pulses']['mw_power'],
#                         self.settings['mw_pulses']['mw_frequency'] * 1e-9, self.settings['tau_time'] * 1e-3))
#
#
#                 axislist[0].set_ylabel('fluorescence [kcps]')
#                 axislist[0].set_xlabel('amplitude [Vpp]')
#
#         else:
#             super(eSensing_swpV, self)._plot(axislist)
#             axislist[0].set_xlabel('amplitude [Vpp]')
#             axislist[0].legend(labels=('PDD up {:.0f}kcps'.format(np.mean(data['counts'][:, 0])),
#                                        'PDD down {:.0f}kcps'.format(np.mean(data['counts'][:, 1])),
#                                        'E-Sensing up {:.0f}kcps'.format(np.mean(data['counts'][:, 2])),
#                                        'E-Sensing down {:.0f}kcps'.format(np.mean(data['counts'][:, 3]))),
#                                fontsize=10)
#             axislist[0].set_title(
#                 '(initial plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
#                     'microwave_channel'] + '\nPDD: {:s} {:d} block(s), {:s} electrometry\n after {:d} averages\n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz \n tau = {:.2f} us'.format(
#                     self.settings['decoupling_seq']['type'], self.settings['decoupling_seq']['num_of_pulse_blocks'],
#                     self.settings['decoupling_seq']['eSensing_type'],
#                     avg_block_number, self.settings['mw_pulses']['mw_power'],
#                     self.settings['mw_pulses']['mw_frequency'] * 1e-9, self.settings['tau_time'] * 1e-3))
#
#             # if self.settings['decoupling_seq']['do_eSensing']:
#             #     axislist[0].legend(labels=('Echo up {:.0f}kcps'.format(np.mean(data['counts'][:, 0])),
#             #                                'Echo down {:.0f}kcps'.format(np.mean(data['counts'][:, 1])),
#             #                                'E-Sensing up {:.0f}kcps'.format(np.mean(data['counts'][:, 0])),
#             #                                'E-Sensing down {:.0f}kcps'.format(np.mean(data['counts'][:, 1]))),
#             #                        fontsize=10)
#             #     axislist[0].set_title(
#             #         '(initial plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
#             #             'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s), {:s} electrometry\n averaged over {:d} blocks\n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz'.format(
#             #             self.settings['decoupling_seq']['type'], self.settings['decoupling_seq']['num_of_pulse_blocks'],
#             #             self.settings['decoupling_seq']['eSensing_type'],
#             #             avg_block_number, self.settings['mw_pulses']['mw_power'],
#             #             self.settings['mw_pulses']['mw_frequency'] * 1e-9))
#             #
#             # else:
#             #     axislist[0].legend(labels=('Echo up {:.0f}kcps'.format(np.mean(data['counts'][:, 0])),
#             #                                'Echo down {:.0f}kcps'.format(np.mean(data['counts'][:, 1]))), fontsize=10)
#             #     axislist[0].set_title(
#             #         '(initial plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
#             #             'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s) \n averaged over {:d} blocks\n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz'.format(
#             #             self.settings['decoupling_seq']['type'], self.settings['decoupling_seq']['num_of_pulse_blocks'],
#             #             avg_block_number, self.settings['mw_pulses']['mw_power'],
#             #             self.settings['mw_pulses']['mw_frequency'] * 1e-9))
#
#     def _update_plot(self, axislist):
#
#         if len(axislist[0].lines) == 0:
#             self._plot(axislist)
#             return
#
#         super(eSensing_swpV, self)._update_plot(axislist)
#         axislist[0].legend(labels=('PDD up {:.0f}kcps'.format(np.mean(self.data['counts'][:, 0])),
#                                    'PDD down {:.0f}kcps'.format(np.mean(self.data['counts'][:, 1])),
#                                    'E-Sensing up {:.0f}kcps'.format(np.mean(self.data['counts'][:, 2])),
#                                    'E-Sensing down {:.0f}kcps'.format(np.mean(self.data['counts'][:, 3]))),
#                            fontsize=10)
#         axislist[0].set_title(
#             '(updating plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
#                 'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s), {:s} electrometry\n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz \n tau = {:.2f} us'.format(
#                 self.settings['decoupling_seq']['type'], self.settings['decoupling_seq']['num_of_pulse_blocks'],
#                 self.settings['decoupling_seq']['eSensing_type'],
#                 self.settings['mw_pulses']['mw_power'],
#                 self.settings['mw_pulses']['mw_frequency'] * 1e-9, self.settings['tau_time'] * 1e-3))
#
#         # if self.settings['decoupling_seq']['do_eSensing']:
#         #     axislist[0].legend(labels=('Echo up {:.0f}kcps'.format(np.mean(self.data['counts'][:, 0])),
#         #                                'Echo down {:.0f}kcps'.format(np.mean(self.data['counts'][:, 1])),
#         #                                'E-Sensing up {:.0f}kcps'.format(np.mean(self.data['counts'][:, 0])),
#         #                                'E-Sensing down {:.0f}kcps'.format(np.mean(self.data['counts'][:, 1]))),
#         #                        fontsize=10)
#         #     axislist[0].set_title(
#         #         '(updating plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
#         #             'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s), {:s} electrometry\n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz'.format(
#         #             self.settings['decoupling_seq']['type'], self.settings['decoupling_seq']['num_of_pulse_blocks'],
#         #             self.settings['decoupling_seq']['eSensing_type'],
#         #             self.settings['mw_pulses']['mw_power'],
#         #             self.settings['mw_pulses']['mw_frequency'] * 1e-9))
#         #
#         # else:
#         #     axislist[0].legend(labels=('Echo up {:.0f}kcps'.format(np.mean(self.data['counts'][:, 0])),
#         #                                'Echo down {:.0f}kcps'.format(np.mean(self.data['counts'][:, 1]))), fontsize=10)
#         #
#         #     axislist[0].set_title('(updating plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
#         #         'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s) \n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz'.format(
#         #         self.settings['decoupling_seq']['type'], self.settings['decoupling_seq']['num_of_pulse_blocks'],
#         #         self.settings['mw_pulses']['mw_power'],
#         #         self.settings['mw_pulses']['mw_frequency'] * 1e-9))
#
#     def _plot_validate(self, axes_list):
#         """
#         Preview pulse sequence by plotting first and last sequence to plots 1 and 0, respectively
#         Args:
#             axes_list: List containing axes to plot to
#
#         Returns:
#
#         """
#         print('printing validated pulses')
#         axis0 = axes_list[0]
#         axis1 = axes_list[1]
#         axis0.clear()
#         axis1.clear()
#
#         pulse_sequences, tau_list, _ = self.create_pulse_sequences(logging=False)
#         # print('here')
#
#       #  if pulse_sequences[0]:
#         if pulse_sequences:
#             plot_pulses(axis0, pulse_sequences[0])
#             axis0.set_title('Pulse Visualization for Minimum Amplitude = {:f} Vpp'.format(tau_list[0]))
#             plot_pulses(axis1, pulse_sequences[-1])
#             axis1.set_title('Pulse Visualization for Maximum Amplitude = {:f} Vpp'.format(tau_list[-1]))
#         else:
#             print('no pulse sequences passed validation!!!')

class angle_sensing(PulsedExperimentBaseScript):
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_power', -15.0, float, 'microwave power in dB'),
            Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
            Parameter('microwave_channel', 'i', ['i', 'q'], 'Channel to use for the first pi/2 mw pulse'),
            Parameter('pi_pulse_time', 50.0, float, 'time duration of a pi pulse (in ns) for channel i'),
            Parameter('pi_half_pulse_time', 25.0, float, 'time duration of a pi/2 pulse (in ns) for channel i'),
            Parameter('3pi_half_pulse_time', 75.0, float, 'time duration of a 3pi/2 pulse (in ns) for channel i'),
            Parameter('Q_same_as_I', True, bool, 'check if I and Q are symmetric'),
            Parameter('pi_pulse_time_q', 50.0, float, 'time duration of a pi pulse (in ns) for channel q'),
            Parameter('pi_half_pulse_time_q', 25.0, float, 'time duration of a pi/2 pulse (in ns) for channel q'),
            Parameter('3pi_half_pulse_time_q', 75.0, float, 'time duration of a 3pi/2 pulse (in ns) for channel q')
        ]),
        Parameter('voltages', [Parameter('min_vol', 0.0, float, 'Minimum DC voltage [Vdc] on the coil'),
                               Parameter('max_vol', 1.1, float, 'Maximum DC voltage [Vdc] on the coil'),
                               Parameter('vol_step', 0.1, float, 'Voltage step [Vdc] in the sweep'),
                               ]),
        Parameter('tau_time', 1800, float, '[us] choose the fixed tau time between two pi/2 pulses'),
        Parameter('decoupling_seq', [
            Parameter('type', 'CPMG', ['spin_echo', 'CPMG'], 'type of dynamical decoupling sequences'),
            Parameter('num_of_pulse_blocks', 1, int, 'number of pulse blocks.')
        ]),
        Parameter('min_pulse_dur', 15, [10, 15, 50], 'Minimum allowed pulse duration (ns)'),
        Parameter('num_averages', 500000, int, 'number of averages (>100000)'),
    ]

    _INSTRUMENTS = {'daq': NI6353, 'PB': LISE607RTPulseBlaster, 'mw_gen_a': AgilentMicrowaveGenerator,
                    'fg': Agilent33120A}

    def _function(self):

        ## Turn off green light (the pulse blaster will pulse it on when needed)
        self.instruments['PB']['instance'].update({'laser': {'status': False}})
        self.instruments['PB']['instance'].update({'apd_switch': {'status': False}})

        # set up the RF generator
        self.instruments['mw_gen_a']['instance'].update({'enable_IQ': True})
        self.instruments['mw_gen_a']['instance'].update({'enable_modulation': True})
        self.instruments['mw_gen_a']['instance'].update({'freq_mode': 'CW'})
        self.instruments['mw_gen_a']['instance'].update({'power_mode': 'CW'})
        self.instruments['mw_gen_a']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power']})
        self.instruments['mw_gen_a']['instance'].update({'frequency': self.settings['mw_pulses']['mw_frequency']})
        self.instruments['mw_gen_a']['instance'].update({'enable_output': True})

        # set up the function generator to be in the DC mode
        self.instruments['fg']['instance'].update({'display_on': False})
        self.instruments['fg']['instance'].update({'burst_mod': False})
        self.instruments['fg']['instance'].update({'offset': 0.0})

        self.avg_block_number = super(angle_sensing, self)._function()

        # this is to flag when the experiment is done or stopped and this allows proper saving even fit fails
        self.data['exp_finished'] = 1
        self.data['avg_block_number'] = self.avg_block_number

        # turn off laser, apd_switch and MW (just in case)
        self.instruments['PB']['instance'].update({'laser': {'status': False}})
        self.instruments['PB']['instance'].update({'apd_switch': {'status': False}})
        self.instruments['mw_gen_a']['instance'].update({'enable_output': False})
        self.instruments['mw_gen_a']['instance'].update({'enable_IQ': False})
        self.instruments['mw_gen_a']['instance'].update({'enable_modulation': False})
        self.instruments['fg']['instance'].update({'display_on': True})
        self.instruments['fg']['instance'].update({'amplitude': 0.1})
        self.instruments['fg']['instance'].update({'offset': 0.0})

        # no fitting for now
        if 'counts' in self.data.keys() and 'tau' in self.data.keys() and not 0 in (
                self.data['counts'][:, 1] + self.data['counts'][:, 0]):
            self.data['norm_echo'] = 2. * (- self.data['counts'][:, 1] + self.data['counts'][:, 0]) / (
                        self.data['counts'][:, 1] + self.data['counts'][:, 0])

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
        max_range = int(np.floor((self.settings['voltages']['max_vol'] - self.settings['voltages']['min_vol']) /
                                 self.settings['voltages']['vol_step']))
        vol_list = np.array(
            [self.settings['voltages']['min_vol'] + i * self.settings['voltages']['vol_step'] for i in
             range(max_range)])
        # ignore the sequence if the input voltage is too large or too small
        vol_list = [x for x in vol_list if (x <= 15) and (x >= 0.001 or x == 0.0)]

        print('Vdc_list:', vol_list)

        nv_reset_time = self.settings['read_out']['nv_reset_time']
        delay_readout = self.settings['read_out']['delay_readout']
        pi_time = self.settings['mw_pulses']['pi_pulse_time']
        pi_half_time = self.settings['mw_pulses']['pi_half_pulse_time']
        three_pi_half_time = self.settings['mw_pulses']['3pi_half_pulse_time']
        if self.settings['mw_pulses']['Q_same_as_I']:
            pi_time_q = pi_time
            pi_half_time_q = pi_half_time
            three_pi_half_time_q = three_pi_half_time
        else:
            pi_time_q = self.settings['mw_pulses']['pi_pulse_time']
            pi_half_time_q = self.settings['mw_pulses']['pi_half_pulse_time']
            three_pi_half_time_q = self.settings['mw_pulses']['3pi_half_pulse_time']

        if self.settings['mw_pulses']['microwave_channel'] == 'i':
            microwave_channel_1 = 'microwave_i'
            microwave_channel_2 = 'microwave_q'
            pi_time_1 = pi_time
            pi_half_time_1 = pi_half_time
            three_pi_half_time_1 = three_pi_half_time
            pi_time_2 = pi_time_q
            pi_half_time_2 = pi_half_time_q
            three_pi_half_time_2 = three_pi_half_time_q
        else:
            microwave_channel_1 = 'microwave_q'
            microwave_channel_2 = 'microwave_i'
            pi_time_1 = pi_time_q
            pi_half_time_1 = pi_half_time_q
            three_pi_half_time_1 = three_pi_half_time_q
            pi_time_2 = pi_time
            pi_half_time_2 = pi_half_time
            three_pi_half_time_2 = three_pi_half_time


        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        number_of_pulse_blocks = self.settings['decoupling_seq']['num_of_pulse_blocks']
        tau_total = self.settings['tau_time']

        if self.settings['decoupling_seq']['type'] == 'spin_echo':
            for vol in vol_list:
                tau = tau_total / (1 * self.settings['decoupling_seq']['num_of_pulse_blocks'])
                pulse_sequence = []

                # ECHO SEQUENCE:
                start_of_first_HE = laser_off_time
                # the first pi/2 pulse

                pulse_sequence = [Pulse(microwave_channel_1, start_of_first_HE, pi_half_time)]
                section_begin_time = start_of_first_HE + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                for i in range(0, number_of_pulse_blocks):
                    # pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                    section_begin_time += 1 * tau

                # the second pi/2 pulse and readout
                pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + tau / 2, pi_half_time))
                pulse_sequence.append(Pulse('laser', section_begin_time + tau / 2 + pi_half_time + delay_mw_readout,
                                             nv_reset_time))
                pulse_sequence.append(Pulse('apd_readout',
                                             section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + delay_readout,
                                             meas_time))
                start_of_second_HE = section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time

                # the first pi/2 pulse
                pulse_sequence.append(Pulse(microwave_channel_1, start_of_second_HE, pi_half_time))
                section_begin_time = start_of_second_HE + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                for i in range(0, number_of_pulse_blocks):
                    # pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                    section_begin_time += 1 * tau

                # the second 3*pi/2 pulse and readout
                pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + tau / 2, three_pi_half_time))
                pulse_sequence.append(Pulse('laser',
                                             section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout,
                                             nv_reset_time))
                pulse_sequence.append(Pulse('apd_readout',
                                             section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + delay_readout,
                                             meas_time))

                # # Y Echo
                # if self.settings['do_Y_readout']:
                #     start_of_DEER = section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time
                #     # the first pi/2 pulse
                #     pulse_sequence.append(Pulse(microwave_channel_1, start_of_DEER, pi_half_time))
                #     section_begin_time = start_of_DEER + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                #     for i in range(0, number_of_pulse_blocks):
                #         pulse_sequence.append(
                #             Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                #         # if number_of_pulse_blocks == 1:
                #         #     pulse_sequence.append(
                #         #         Pulse(e_field_channel, section_begin_time + 1 * tau, tau / 2))
                #         # elif i%2 == 0:
                #         #     pulse_sequence.append(
                #         #         Pulse(e_field_channel, section_begin_time + 1 * tau, tau))
                #         section_begin_time += 1 * tau
                #
                #     # the second pi/2 pulse and readout
                #     pulse_sequence.append(Pulse(last_mw_channel, section_begin_time + tau / 2, pi_half_time))
                #     pulse_sequence.append(Pulse('laser', section_begin_time + tau / 2 + pi_half_time + delay_mw_readout,
                #                                 nv_reset_time))
                #     pulse_sequence.append(Pulse('apd_readout',
                #                                 section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + delay_readout,
                #                                 meas_time))
                #
                #     start_of_second_DEER = section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time
                #
                #     # the first pi/2 pulse
                #     pulse_sequence.append(Pulse(microwave_channel_1, start_of_second_DEER, pi_half_time))
                #     section_begin_time = start_of_second_DEER + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                #
                #     for i in range(0, number_of_pulse_blocks):
                #         pulse_sequence.append(
                #             Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                #         # if number_of_pulse_blocks == 1:
                #         #     pulse_sequence.append(
                #         #         Pulse(e_field_channel, section_begin_time + 1 * tau, tau / 2))
                #         # elif i % 2 == 0:
                #         #     pulse_sequence.append(
                #         #         Pulse(e_field_channel, section_begin_time + 1 * tau, tau))
                #         section_begin_time += 1 * tau
                #
                #     # the second pi/2 pulse and readout
                #     pulse_sequence.append(Pulse(last_mw_channel, section_begin_time + tau / 2, three_pi_half_time))
                #     pulse_sequence.append(Pulse('laser', section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout,
                #                                 nv_reset_time))
                #     pulse_sequence.append(Pulse('apd_readout',
                #                                 section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + delay_readout,
                #                                 meas_time))

                # Turn on APD switch all the time
                pulse_sequence.append(Pulse('apd_switch', 0,
                                            section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + nv_reset_time))

                # print(pulse_sequence)
                pulse_sequences.append(pulse_sequence)

        elif self.settings['decoupling_seq']['type'] == 'CPMG':

            for vol in vol_list:
                tau = tau_total / (1 * self.settings['decoupling_seq']['num_of_pulse_blocks'])
                pulse_sequence = []

                # ECHO SEQUENCE:
                start_of_first_HE = laser_off_time

                # the first pi/2 pulse
                pulse_sequence = [Pulse(microwave_channel_1, start_of_first_HE, pi_half_time)]
                section_begin_time = start_of_first_HE + pi_half_time - tau / 2  # for the first pulse, only wait tau/2

                for i in range(0, number_of_pulse_blocks):
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                    section_begin_time += 1 * tau

                # the second 3*pi/2 pulse and readout
                pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + tau / 2, three_pi_half_time))
                pulse_sequence.append(
                    Pulse('laser', section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout,
                          nv_reset_time))
                pulse_sequence.append(Pulse('apd_readout',
                                            section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + delay_readout,
                                            meas_time))

                start_of_second_HE = section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time

                # the first pi/2 pulse
                pulse_sequence.append(Pulse(microwave_channel_1, start_of_second_HE, pi_half_time))

                section_begin_time = start_of_second_HE + pi_half_time - tau / 2  # for the first pulse, only wait tau/2

                for i in range(0, number_of_pulse_blocks):
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                    section_begin_time += 1 * tau

                # the second pi/2 pulse and readout
                pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + tau / 2, pi_half_time))
                pulse_sequence.append(Pulse('laser',
                                            section_begin_time + tau / 2 + pi_half_time + delay_mw_readout,
                                            nv_reset_time))
                pulse_sequence.append(Pulse('apd_readout',
                                            section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + delay_readout,
                                            meas_time))


                # Turn on APD switch all the time
                pulse_sequence.append(Pulse('apd_switch', 0,
                                            section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + nv_reset_time))

                #             # # DEER SEQUENCE
                #             # if self.settings['decoupling_seq']['do_deer']:
                #             #     start_of_DEER = section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time
                #             #     # the first pi/2 pulse
                #             #     pulse_sequence.extend([Pulse(microwave_channel_1, start_of_DEER, pi_half_time)])
                #             #     section_begin_time = start_of_DEER + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                #             #     for i in range(0, number_of_pulse_blocks):
                #             #         pulse_sequence.extend(
                #             #             [Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time),
                #             #              Pulse('RF_switch', section_begin_time + 1 * tau - RF_pi_time / 2, RF_pi_time)
                #             #              ])
                #             #         section_begin_time += 1 * tau
                #             #
                #             #     # the second pi/2 pulse and readout
                #             #     pulse_sequence.extend([Pulse(microwave_channel_1, section_begin_time + tau / 2, pi_half_time),
                #             #                            Pulse('laser',
                #             #                                  section_begin_time + tau / 2 + pi_half_time + delay_mw_readout,
                #             #                                  nv_reset_time),
                #             #                            Pulse('apd_readout',
                #             #                                  section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + delay_readout,
                #             #                                  meas_time)
                #             #                            ])
                #             #
                #             #     start_of_second_DEER = section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time
                #             #     pulse_sequence.extend([Pulse(microwave_channel_1, start_of_second_DEER, pi_half_time)])
                #             #     section_begin_time = start_of_second_DEER + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                #             #
                #             #     for i in range(0, number_of_pulse_blocks):
                #             #         pulse_sequence.extend(
                #             #             [Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time),
                #             #              Pulse('RF_switch', section_begin_time + 1 * tau - RF_pi_time / 2, RF_pi_time)
                #             #              ])
                #             #         section_begin_time += 1 * tau
                #             #
                #             #     # the second 3*pi/2 pulse and readout
                #             #     pulse_sequence.extend([Pulse(microwave_channel_1, section_begin_time + tau / 2, three_pi_half_time),
                #             #                            Pulse('laser',
                #             #                                  section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout,
                #             #                                  nv_reset_time),
                #             #                            Pulse('apd_readout',
                #             #                                  section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + delay_readout,
                #             #                                  meas_time)
                #             #                            ])
                #
                pulse_sequences.append(pulse_sequence)

            print('number of sequences before validation ', len(pulse_sequences))
            # return pulse_sequences, self.settings['num_averages'], tau_list, meas_time
            return pulse_sequences, vol_list, meas_time

    def _run_sweep(self, pulse_sequences, num_loops_sweep, num_daq_reads, tau_list, verbose=False):
        """
        Each pulse sequence specified in pulse_sequences is run num_loops_sweep consecutive times.

        Args:
            pulse_sequences: a list of pulse sequences to run, each corresponding to a different value of tau. Each
                             sequence is a list of Pulse objects specifying a given pulse sequence
            num_loops_sweep: number of times to repeat each sequence before moving on to the next one
            num_daq_reads: number of times the daq must read for each sequence (generally 1, 2, or 3)

        Poststate: self.data['counts'] is updated with the acquired data

        """
        # ER 20180731 set init_fluor to zero
     #   self.data['init_fluor'] = 0.

        # print('start _run_sweep now')


        rand_indexes = list(range(len(pulse_sequences)))
        if self.settings['randomize']:
            random.shuffle(rand_indexes)

        if verbose:
            print(('_run_sweep number of pulse sequences', len(pulse_sequences)))

        for index, sequence in enumerate(pulse_sequences):
            # print('current index')
            # print(index)
            if verbose:
                print(('_run_sweep index', index, len(pulse_sequences)))

            rand_index = rand_indexes[index]
            self.instruments['fg']['instance'].update({'offset': float(tau_list[rand_index])})

            # self.current_rand_index = rand_index
            if self._abort:
                # print('aborted')
                self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})
                # print('mw is off')
                break

            # print('running single sqeucnce')
            result = self._run_single_sequence(pulse_sequences[rand_index], num_loops_sweep, num_daq_reads)  # keep entire array
            # print('single sqeucnce done')
            self.count_data[rand_index] = self.count_data[rand_index] + result
            # print('before counts_to_check')

            counts_to_check = self._normalize_to_kCounts(np.array(result), self.measurement_gate_width, num_loops_sweep)
            self.data['counts'][rand_index] = self._normalize_to_kCounts(self.count_data[rand_index], self.measurement_gate_width,
                                                                    self.current_averages)

            # self.sequence_index = index
            self.sequence_index = rand_index
            # counts_temp = counts_to_check[0]
            counts_temp = np.max(counts_to_check)
            # track to the NV if necessary ER 5/31/17
            if self.settings['Tracking']['on/off']:
                # print('doing tracking')
                # print('counts_temp:', counts_temp)
                print('Tracking: current counts = ', counts_temp, 'kcps')
                if (1+(1-self.settings['Tracking']['threshold']))*self.settings['Tracking']['init_fluor'] < counts_temp or \
                        self.settings['Tracking']['threshold']*self.settings['Tracking']['init_fluor'] > counts_temp:
                    print('(2-threshold)*intial_fluor:',
                          (1 + (1 - self.settings['Tracking']['threshold'])) * self.settings['Tracking']['init_fluor'])
                    print('threshold*intial_fluor:', self.settings['Tracking']['threshold']*self.settings['Tracking']['init_fluor'])
                    if verbose:
                        print('TRACKING TO NV...')
                    self.scripts['find_nv'].run()
                    self.scripts['find_nv'].settings['initial_point'] = self.scripts['find_nv'].data['maximum_point']
                    print('find_nv done')
                    if self.scripts['find_nv'].data[
                        'fluorescence'] == 0.0:  # if it doesn't find an NV, abort the experiment
                        print('Could not find an NV in FindNV.')
                        self.log('Could not find an NV in FindNV.')
                        self._abort = True
                        return  # exit function in case no NV is found
                    else:
                        print('Start optimize_z.')
                        self.scripts['optimize_z'].run()
            self.updateProgress.emit(self._calc_progress(index))

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

        if 'exp_finished' in data.keys() is not None and data['exp_finished'] == 1:
            if 'norm_echo' in data.keys() is not None and data['norm_echo'] is not None:
                tau = data['tau']
                axislist[0].plot(tau, data['norm_echo'])
                axislist[0].legend(labels=('X Echo'), fontsize=10)
                axislist[0].set_title(
                    '(final plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
                        'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s)\n after {:d} averages \n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz'.format(
                        self.settings['decoupling_seq']['type'],
                        self.settings['decoupling_seq']['num_of_pulse_blocks'],
                        avg_block_number,
                        self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency'] * 1e-9))
                axislist[0].set_ylabel('fluorescence contrast')
                axislist[0].set_xlabel('Coil Voltage [V]')
            else:
                axislist[0].plot(data['tau'], data['counts'][:, 0])
                axislist[0].plot(data['tau'], data['counts'][:, 1])

                axislist[0].legend(labels=('X Echo up {:.0f}kcps'.format(np.mean(data['counts'][:, 0])),
                                           'X Echo down {:.0f}kcps'.format(np.mean(data['counts'][:, 1]))),
                                   fontsize=10)
                axislist[0].set_title(
                    '(final plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
                        'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s)\n after {:d} averages\n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz'.format(
                        self.settings['decoupling_seq']['type'], self.settings['decoupling_seq']['num_of_pulse_blocks'],
                        avg_block_number, self.settings['mw_pulses']['mw_power'],
                        self.settings['mw_pulses']['mw_frequency'] * 1e-9))
                axislist[0].set_ylabel('fluorescence [kcps]')
                axislist[0].set_xlabel('Coil Voltage [V]')

        else:
            super(angle_sensing, self)._plot(axislist)

            axislist[0].legend(labels=('X Echo up {:.0f}kcps'.format(np.mean(data['counts'][:, 0])),
                                       'X Echo down {:.0f}kcps'.format(np.mean(data['counts'][:, 1]))), fontsize=10)
            axislist[0].set_title(
                '(initial plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
                    'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s) \n after {:d} averages\n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz'.format(
                    self.settings['decoupling_seq']['type'], self.settings['decoupling_seq']['num_of_pulse_blocks'],
                    avg_block_number, self.settings['mw_pulses']['mw_power'],
                    self.settings['mw_pulses']['mw_frequency'] * 1e-9))

    def _update_plot(self, axislist):

        if len(axislist[0].lines) == 0:
            self._plot(axislist)
            return

        super(angle_sensing, self)._update_plot(axislist)

        axislist[0].legend(labels=('X Echo up {:.0f}kcps'.format(np.mean(self.data['counts'][:, 0])),
                                   'X Echo down {:.0f}kcps'.format(np.mean(self.data['counts'][:, 1]))), fontsize=10)

        axislist[0].set_title('(updating plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
            'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s) \n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz'.format(
            self.settings['decoupling_seq']['type'], self.settings['decoupling_seq']['num_of_pulse_blocks'],
            self.settings['mw_pulses']['mw_power'],
            self.settings['mw_pulses']['mw_frequency'] * 1e-9))
        axislist[0].set_ylabel('fluorescence [kcps]')
        axislist[0].set_xlabel('Coil Voltage [V]')

class eSensing_swpV(PulsedExperimentBaseScript):
    """
    This script does AC electrometry using NV. Tau is fixed, and voltage is swept.
    Voltages are applied using Agilent 33120A.
    MW pulses are applied using Keysight N9310A.

    --> Ziwei Qiu 6/26/2019

    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_power', -10.0, float, 'microwave power in dB'),
            Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
            Parameter('microwave_channel', 'i', ['i', 'q'], 'Channel to use for the first pi/2 mw pulse'),
            Parameter('pi_pulse_time', 50.0, float, 'time duration of a pi pulse (in ns) for channel i'),
            Parameter('pi_half_pulse_time', 25.0, float, 'time duration of a pi/2 pulse (in ns) for channel i'),
            Parameter('3pi_half_pulse_time', 75.0, float, 'time duration of a 3pi/2 pulse (in ns) for channel i'),
            Parameter('Q_same_as_I', True, bool, 'check if I and Q are symmetric'),
            Parameter('pi_pulse_time_q', 50.0, float, 'time duration of a pi pulse (in ns) for channel q'),
            Parameter('pi_half_pulse_time_q', 25.0, float, 'time duration of a pi/2 pulse (in ns) for channel q'),
            Parameter('3pi_half_pulse_time_q', 75.0, float, 'time duration of a 3pi/2 pulse (in ns) for channel q')
        ]),
        Parameter('e_fields', [
            Parameter('min_amp', 0.1, float, '[Vpp] minimum amplitude'),
            Parameter('max_amp', 20.0, float, '[Vpp] maximum amplitude'),
            Parameter('amp_step', 1., [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0],
                      '[Vpp] amplitude step increment'),
            Parameter('offset', 0.0, float, '[Vdc] DC offset voltage'),
            Parameter('wave_shape', 'SQUare', ['SINusoid', 'SQUare'], 'choose the wave shape'),
            Parameter('trigger_latency', 1430.6, float, '[ns] Agilent33120A  trigger latency'),
            Parameter('trigger_time', 500.0, float, '[ns] trigger pulse time (>100ns)'),
        ]),
        Parameter('tau_time', 7000, float, '[ns] choose the fixed tau time between two pi/2 pulses'),
        Parameter('decoupling_seq', [
            Parameter('type', 'XY4', ['spin_echo', 'CPMG', 'XY4', 'XY8'], 'type of dynamical decoupling sequences'),
            Parameter('num_of_pulse_blocks', 1, int, 'number of pulse blocks.'),
            Parameter('eSensing_type', 'sine', ['sine', 'cosine'], 'choose to do sine or cosine electrometry'),
        ]),
        Parameter('min_pulse_dur', 15, [10, 15, 50], 'Minimum allowed pulse duration (ns)'),
        Parameter('num_averages', 500000, int, 'number of averages (>100000)'),
    ]

    _INSTRUMENTS = {'daq': NI6353, 'PB': LISE607RTPulseBlaster, 'mw_gen_a': AgilentMicrowaveGenerator, 'fg': Agilent33120A}

    def _function(self):

        self.data['fits'] = None

        ## Turn off green light (the pulse blaster will pulse it on when needed)
        self.instruments['PB']['instance'].update({'laser': {'status': False}})
        self.instruments['PB']['instance'].update({'apd_switch': {'status': False}})

        # set up the RF generator
        self.instruments['mw_gen_a']['instance'].update({'enable_IQ': True})
        self.instruments['mw_gen_a']['instance'].update({'enable_modulation': True})
        self.instruments['mw_gen_a']['instance'].update({'freq_mode': 'CW'})
        self.instruments['mw_gen_a']['instance'].update({'power_mode': 'CW'})
        self.instruments['mw_gen_a']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power']})
        self.instruments['mw_gen_a']['instance'].update({'frequency': self.settings['mw_pulses']['mw_frequency']})
        self.instruments['mw_gen_a']['instance'].update({'enable_output': True})

        # set up the function generator
        self.instruments['fg']['instance'].update({'display_on': False})
        self.instruments['fg']['instance'].update({'burst_mod': True})
        # self.instruments['fg']['instance'].update({'amplitude': self.settings['e_fields']['min_amp']})
        self.instruments['fg']['instance'].update({'offset': self.settings['e_fields']['offset']})
        self.instruments['fg']['instance'].update({'wave_shape': self.settings['e_fields']['wave_shape']})
        valid = True
        if self.settings['decoupling_seq']['type'] == 'spin_echo' or self.settings['decoupling_seq']['type'] == 'CPMG':
            tau = self.settings['tau_time'] / (1 * self.settings['decoupling_seq']['num_of_pulse_blocks'])
            if self.settings['decoupling_seq']['num_of_pulse_blocks'] == 1:
                fg_freq = 1E9 / (1 * tau)
                burst_cnt = 1
                burst_phase = 0.0
                if (fg_freq / 1E6) / burst_cnt > 1:
                    self.log('Invalid setting: (fg_freq / 1E6) / burst_cnt > 1.')
                    valid = False
                else:
                    self.instruments['fg']['instance'].update({'burst_phase': burst_phase})
                    self.instruments['fg']['instance'].update({'frequency': fg_freq})
                    self.instruments['fg']['instance'].update({'burst_count': burst_cnt})
            elif self.settings['decoupling_seq']['num_of_pulse_blocks']%2 == 0:
                fg_freq = 1E9 / (2 * tau)
                burst_cnt = self.settings['decoupling_seq']['num_of_pulse_blocks'] / 2
                burst_phase = 90.0
                if (fg_freq / 1E6) / burst_cnt > 1:
                    self.log('Invalid setting: (fg_freq / 1E6) / burst_cnt > 1.')
                    valid = False
                else:
                    self.instruments['fg']['instance'].update({'burst_phase': burst_phase})
                    self.instruments['fg']['instance'].update({'frequency': fg_freq})
                    self.instruments['fg']['instance'].update({'burst_count': burst_cnt})
            else:
                self.log('Invalid setting: num_of_pulse_blocks needs to be 1 or even number if using spin-echo or CPMG.')
                valid = False
        elif self.settings['decoupling_seq']['type'] == 'XY4':
            tau = self.settings['tau_time'] / (4 * self.settings['decoupling_seq']['num_of_pulse_blocks'])
            fg_freq = 1E9 / (2 * tau)
            burst_cnt = 2 * self.settings['decoupling_seq']['num_of_pulse_blocks']
            burst_phase = 90.0
            if (fg_freq / 1E6) / burst_cnt > 1:
                self.log('Invalid setting: (fg_freq / 1E6) / burst_cnt > 1.')
                valid = False
            else:
                self.instruments['fg']['instance'].update({'burst_phase': burst_phase})
                self.instruments['fg']['instance'].update({'frequency': fg_freq})
                self.instruments['fg']['instance'].update({'burst_count': burst_cnt})
        elif self.settings['decoupling_seq']['type'] == 'XY8':
            tau = self.settings['tau_time'] / (8 * self.settings['decoupling_seq']['num_of_pulse_blocks'])
            fg_freq = 1E9 / (2 * tau)
            burst_cnt = 4 * self.settings['decoupling_seq']['num_of_pulse_blocks']
            burst_phase = 90.0
            if (fg_freq / 1E6) / burst_cnt > 1:
                self.log('Invalid setting: (fg_freq / 1E6) / burst_cnt > 1.')
                valid = False
            else:
                self.instruments['fg']['instance'].update({'burst_phase': burst_phase})
                self.instruments['fg']['instance'].update({'frequency': fg_freq})
                self.instruments['fg']['instance'].update({'burst_count': burst_cnt})
        else:
            self.log('Invalid setting: choose the right decoupling_seq type')
            valid = False

        if valid:
            self.avg_block_number = super(eSensing_swpV, self)._function()

            # this is to flag when the experiment is done or stopped and this allows proper saving even fit fails
            self.data['exp_finished'] = 1
            self.data['avg_block_number'] = self.avg_block_number

        # turn off laser, apd_switch and MW (just in case)
        self.instruments['PB']['instance'].update({'laser': {'status': False}})
        self.instruments['PB']['instance'].update({'apd_switch': {'status': False}})
        self.instruments['mw_gen_a']['instance'].update({'enable_output': False})
        self.instruments['mw_gen_a']['instance'].update({'enable_IQ': False})
        self.instruments['mw_gen_a']['instance'].update({'enable_modulation': False})
        self.instruments['fg']['instance'].update({'display_on': True})
        self.instruments['fg']['instance'].update({'amplitude': 0.1})
        self.instruments['fg']['instance'].update({'offset': 0.0})

        # no fitting for now
        if 'counts' in self.data.keys() and 'tau' in self.data.keys() and not 0 in (
                self.data['counts'][:, 1] + self.data['counts'][:, 0]) and not 0 in (
                self.data['counts'][:, 3] + self.data['counts'][:, 2]):

            self.data['norm_echo'] = 2. * (- self.data['counts'][:, 1] + self.data['counts'][:, 0]) / (
                    self.data['counts'][:, 1] + self.data['counts'][:, 0])
            self.data['norm_esig'] = 2. * (- self.data['counts'][:, 3] + self.data['counts'][:, 2]) / (
                    self.data['counts'][:, 3] + self.data['counts'][:, 2])

            tau = self.data['tau']
            counts = self.data['norm_echo']
            esig_counts = self.data['norm_esig']
            fit_success = False

            # try:
            #     fits = fit_exp_decay(tau, counts, offset = True, verbose = True)
            #     self.data['fits'] = fits
            #     fit_success = True
            # except:
            #     self.data['fits'] = None
            #     self.log('T2 fit failed')
            #     fit_success = False
            #
            # if self.settings['decoupling_seq']['do_eSensing'] and not 0 in (self.data['counts'][:, 3] + self.data['counts'][:, 2]):
            #     self.data['norm_esig'] = 2. * (- self.data['counts'][:, 3] + self.data['counts'][:, 2]) / (
            #                 self.data['counts'][:, 3] + self.data['counts'][:, 2])
            #
            #     esig_counts = self.data['norm_esig']
                # try:
                #     fits_esig = fit_exp_decay(tau, esig_counts, offset = True, verbose = True)
                #     self.data['fits_esig'] = fits_esig
                #     fit_esig_success = True
                # except:
                #     self.data['fits_esig'] = None
                #     self.log('E-sensing T2 fit failed')
                #     fit_esig_success = False

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
        max_range = int(np.floor((self.settings['e_fields']['max_amp'] - self.settings['e_fields']['min_amp']) /
                                 self.settings['e_fields']['amp_step']))
        amp_list = np.array(
            [self.settings['e_fields']['min_amp'] + i * self.settings['e_fields']['amp_step'] for i in
             range(max_range)])
        # ignore the sequence if the amplitude is smaller than 15Vpp
        amp_list = [x for x in amp_list if x >= 0.1]

        print('Vpp_list:', amp_list)

        nv_reset_time = self.settings['read_out']['nv_reset_time']
        delay_readout = self.settings['read_out']['delay_readout']

        pi_time = self.settings['mw_pulses']['pi_pulse_time']
        pi_half_time = self.settings['mw_pulses']['pi_half_pulse_time']
        three_pi_half_time = self.settings['mw_pulses']['3pi_half_pulse_time']
        if self.settings['mw_pulses']['Q_same_as_I']:
            pi_time_q = pi_time
            pi_half_time_q = pi_half_time
            three_pi_half_time_q = three_pi_half_time
        else:
            pi_time_q = self.settings['mw_pulses']['pi_pulse_time']
            pi_half_time_q = self.settings['mw_pulses']['pi_half_pulse_time']
            three_pi_half_time_q = self.settings['mw_pulses']['3pi_half_pulse_time']

        if self.settings['mw_pulses']['microwave_channel'] == 'i':
            microwave_channel_1 = 'microwave_i'
            microwave_channel_2 = 'microwave_q'
            pi_time_1 = pi_time
            pi_half_time_1 = pi_half_time
            three_pi_half_time_1 = three_pi_half_time
            pi_time_2 = pi_time_q
            pi_half_time_2 = pi_half_time_q
            three_pi_half_time_2 = three_pi_half_time_q
        else:
            microwave_channel_1 = 'microwave_q'
            microwave_channel_2 = 'microwave_i'
            pi_time_1 = pi_time_q
            pi_half_time_1 = pi_half_time_q
            three_pi_half_time_1 = three_pi_half_time_q
            pi_time_2 = pi_time
            pi_half_time_2 = pi_half_time
            three_pi_half_time_2 = three_pi_half_time

        if self.settings['decoupling_seq']['eSensing_type'] == 'sine':
            last_mw_channel = microwave_channel_2
            last_pi_time = pi_time_2
            last_pi_half_time = pi_half_time_2
            last_three_pi_half_time = three_pi_half_time_2
        else:
            last_mw_channel = microwave_channel_1
            last_pi_time = pi_time_1
            last_pi_half_time = pi_half_time_1
            last_three_pi_half_time = three_pi_half_time_1


        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        trigger_latency = self.settings['e_fields']['trigger_latency']
        trigger_time = self.settings['e_fields']['trigger_time']

        number_of_pulse_blocks = self.settings['decoupling_seq']['num_of_pulse_blocks']
        tau_total = self.settings['tau_time']

        if self.settings['decoupling_seq']['type'] == 'spin_echo':
            for amp in amp_list:
                tau = tau_total / (1 * self.settings['decoupling_seq']['num_of_pulse_blocks'])
                pulse_sequence = []

                # ECHO SEQUENCE:
                start_of_first_HE = laser_off_time

                # the first pi/2 pulse
                pulse_sequence = [Pulse(microwave_channel_1, start_of_first_HE, pi_half_time_1)]
                section_begin_time = start_of_first_HE + pi_half_time_1 - tau / 2  # for the first pulse, only wait tau/2
                for i in range(0, number_of_pulse_blocks):

                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time_1 / 2, pi_time_1))
                    section_begin_time += 1 * tau

                # the second pi/2 pulse and readout
                pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + tau / 2, pi_half_time_1))
                pulse_sequence.append(Pulse('laser', section_begin_time + tau / 2 + pi_half_time_1 + delay_mw_readout,
                                             nv_reset_time))
                pulse_sequence.append(Pulse('apd_readout',
                                             section_begin_time + tau / 2 + pi_half_time_1 + delay_mw_readout + delay_readout,
                                             meas_time))
                start_of_second_HE = section_begin_time + tau / 2 + pi_half_time_1 + delay_mw_readout + nv_reset_time + laser_off_time

                # the first pi/2 pulse
                pulse_sequence.append(Pulse(microwave_channel_1, start_of_second_HE, pi_half_time_1))
                section_begin_time = start_of_second_HE + pi_half_time_1 - tau / 2  # for the first pulse, only wait tau/2
                for i in range(0, number_of_pulse_blocks):
                    # pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time_1 / 2, pi_time_1))
                    section_begin_time += 1 * tau

                # the second 3*pi/2 pulse and readout
                pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + tau / 2, three_pi_half_time_1))
                pulse_sequence.append(Pulse('laser',
                                             section_begin_time + tau / 2 + three_pi_half_time_1 + delay_mw_readout,
                                             nv_reset_time))
                pulse_sequence.append(Pulse('apd_readout',
                                             section_begin_time + tau / 2 + three_pi_half_time_1 + delay_mw_readout + delay_readout,
                                             meas_time))

                # E-Sensing SEQUENCE
                if True:
                    start_of_DEER = section_begin_time + tau / 2 + three_pi_half_time_1 + delay_mw_readout + nv_reset_time + laser_off_time
                    # the first pi/2 pulse
                    pulse_sequence.append(Pulse(microwave_channel_1, start_of_DEER, pi_half_time_1))
                    pulse_sequence.append(Pulse('fg_trigger', start_of_DEER - trigger_latency, trigger_time))
                    section_begin_time = start_of_DEER + pi_half_time_1 - tau / 2  # for the first pulse, only wait tau/2
                    for i in range(0, number_of_pulse_blocks):
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time_1 / 2, pi_time_1))
                        section_begin_time += 1 * tau

                    # the second pi/2 pulse and readout
                    pulse_sequence.append(Pulse(last_mw_channel, section_begin_time + tau / 2, last_pi_half_time))
                    pulse_sequence.append(Pulse('laser', section_begin_time + tau / 2 + last_pi_half_time + delay_mw_readout,
                                                nv_reset_time))
                    pulse_sequence.append(Pulse('apd_readout',
                                                section_begin_time + tau / 2 + last_pi_half_time + delay_mw_readout + delay_readout,
                                                meas_time))

                    start_of_second_DEER = section_begin_time + tau / 2 + last_pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time

                    # the first pi/2 pulse
                    pulse_sequence.append(Pulse(microwave_channel_1, start_of_second_DEER, pi_half_time_1))
                    pulse_sequence.append(Pulse('fg_trigger', start_of_second_DEER - trigger_latency, trigger_time))
                    section_begin_time = start_of_second_DEER + pi_half_time_1 - tau / 2  # for the first pulse, only wait tau/2

                    for i in range(0, number_of_pulse_blocks):
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time_1 / 2, pi_time_1))
                        section_begin_time += 1 * tau

                    # the second pi/2 pulse and readout
                    pulse_sequence.append(Pulse(last_mw_channel, section_begin_time + tau / 2, last_three_pi_half_time))
                    pulse_sequence.append(Pulse('laser', section_begin_time + tau / 2 + last_three_pi_half_time + delay_mw_readout,
                                                nv_reset_time))
                    pulse_sequence.append(Pulse('apd_readout',
                                                section_begin_time + tau / 2 + last_three_pi_half_time + delay_mw_readout + delay_readout,
                                                meas_time))

                # Turn on APD switch all the time
                pulse_sequence.append(Pulse('apd_switch', 0,
                                            section_begin_time + tau / 2 + last_three_pi_half_time + delay_mw_readout + nv_reset_time))

                # print(pulse_sequence)
                pulse_sequences.append(pulse_sequence)

            print('number of sequences before validation ', len(pulse_sequences))
            return pulse_sequences, amp_list, meas_time

        elif self.settings['decoupling_seq']['type'] == 'XY4':

            for amp in amp_list:
                tau = tau_total / (4 * self.settings['decoupling_seq']['num_of_pulse_blocks'])
                pulse_sequence = []

                # ECHO SEQUENCE:
                start_of_first_HE = laser_off_time
                # the first pi/2 pulse
                pulse_sequence = [Pulse(microwave_channel_1, start_of_first_HE, pi_half_time_1)]
                section_begin_time = start_of_first_HE + pi_half_time_1 - tau / 2  # for the first pulse, only wait tau/2
                for i in range(0, number_of_pulse_blocks):
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time_2 / 2, pi_time_2))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 2 * tau - pi_time_1 / 2, pi_time_1))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time_1 / 2, pi_time_1))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time_2 / 2, pi_time_2))

                    section_begin_time += 4 * tau

                # the second 3*pi/2 pulse and readout
                pulse_sequence.append(
                    Pulse(microwave_channel_1, section_begin_time + tau / 2, three_pi_half_time_1))
                pulse_sequence.append(
                    Pulse('laser', section_begin_time + tau / 2 + three_pi_half_time_1 + delay_mw_readout,
                          nv_reset_time))
                pulse_sequence.append(
                    Pulse('apd_readout',
                          section_begin_time + tau / 2 + three_pi_half_time_1 + delay_mw_readout + delay_readout,
                          meas_time))
                start_of_second_HE = section_begin_time + tau / 2 + three_pi_half_time_1 + delay_mw_readout + nv_reset_time + laser_off_time

                # the first pi/2 pulse
                pulse_sequence.append(
                    Pulse(microwave_channel_1, start_of_second_HE, pi_half_time_1))
                section_begin_time = start_of_second_HE + pi_half_time_1 - tau / 2  # for the first pulse, only wait tau/2
                for i in range(0, number_of_pulse_blocks):
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time_2 / 2, pi_time_2))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 2 * tau - pi_time_1 / 2, pi_time_1))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time_1 / 2, pi_time_1))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time_2 / 2, pi_time_2))
                    section_begin_time += 4 * tau

                # the second pi/2 pulse and readout
                pulse_sequence.append(
                    Pulse(microwave_channel_1, section_begin_time + tau / 2, pi_half_time_1))
                pulse_sequence.append(
                    Pulse('laser',
                          section_begin_time + tau / 2 + pi_half_time_1 + delay_mw_readout,
                          nv_reset_time))
                pulse_sequence.append(
                    Pulse('apd_readout',
                          section_begin_time + tau / 2 + pi_half_time_1 + delay_mw_readout + delay_readout,
                          meas_time))

                # E-Sensing SEQUENCE
                if True:
                    start_of_DEER = section_begin_time + tau / 2 + pi_half_time_1 + delay_mw_readout + nv_reset_time + laser_off_time
                    # the first pi/2 pulse
                    pulse_sequence.append(Pulse(microwave_channel_1, start_of_DEER, pi_half_time_1))
                    pulse_sequence.append(Pulse('fg_trigger', start_of_DEER - trigger_latency, trigger_time))
                    section_begin_time = start_of_DEER + pi_half_time_1 - tau / 2  # for the first pulse, only wait tau/2
                    for i in range(0, number_of_pulse_blocks):
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time_2 / 2, pi_time_2))
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 2 * tau - pi_time_1 / 2, pi_time_1))
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time_1 / 2, pi_time_1))
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time_2 / 2, pi_time_2))

                        section_begin_time += 4 * tau

                    # the second 3*pi/2 pulse and readout
                    pulse_sequence.append(
                        Pulse(last_mw_channel, section_begin_time + tau / 2, last_three_pi_half_time))
                    pulse_sequence.append(
                        Pulse('laser', section_begin_time + tau / 2 + last_three_pi_half_time + delay_mw_readout,
                              nv_reset_time))
                    pulse_sequence.append(
                        Pulse('apd_readout',
                              section_begin_time + tau / 2 + last_three_pi_half_time + delay_mw_readout + delay_readout,
                              meas_time))

                    start_of_second_DEER = section_begin_time + tau / 2 + last_three_pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time

                    # the first pi/2 pulse
                    pulse_sequence.append(Pulse(microwave_channel_1, start_of_second_DEER, pi_half_time_1))
                    pulse_sequence.append(Pulse('fg_trigger', start_of_second_DEER - trigger_latency, trigger_time))
                    section_begin_time = start_of_second_DEER + pi_half_time_1 - tau / 2  # for the first pulse, only wait tau/2

                    for i in range(0, number_of_pulse_blocks):
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time_2 / 2, pi_time_2))
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 2 * tau - pi_time_1 / 2, pi_time_1))
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time_1 / 2, pi_time_1))
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time_2 / 2, pi_time_2))
                        section_begin_time += 4 * tau

                    # the second pi/2 pulse and readout
                    pulse_sequence.append(
                        Pulse(last_mw_channel, section_begin_time + tau / 2, last_pi_half_time))
                    pulse_sequence.append(
                        Pulse('laser',
                              section_begin_time + tau / 2 + last_pi_half_time + delay_mw_readout,
                              nv_reset_time))
                    pulse_sequence.append(
                        Pulse('apd_readout',
                              section_begin_time + tau / 2 + last_pi_half_time + delay_mw_readout + delay_readout,
                              meas_time))

                pulse_sequence.append(Pulse('apd_switch', 0,
                                            section_begin_time + tau / 2 + last_pi_half_time + delay_mw_readout + nv_reset_time))
                pulse_sequences.append(pulse_sequence)

            print('number of sequences before validation ', len(pulse_sequences))
            # return pulse_sequences, self.settings['num_averages'], tau_list, meas_time
            return pulse_sequences, amp_list, meas_time

        elif self.settings['decoupling_seq']['type'] == 'XY8':

            for amp in amp_list:
                tau = tau_total / (8 * self.settings['decoupling_seq']['num_of_pulse_blocks'])
                pulse_sequence = []

                # ECHO SEQUENCE:
                start_of_first_HE = laser_off_time

                # the first pi/2 pulse
                pulse_sequence.append(
                    Pulse(microwave_channel_1, start_of_first_HE, pi_half_time_1))

                section_begin_time = start_of_first_HE + pi_half_time_1 - tau / 2  # for the first pulse, only wait tau/2
                for i in range(0, number_of_pulse_blocks):

                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time_1 / 2, pi_time_1))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 2 * tau - pi_time_2 / 2, pi_time_2))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time_1 / 2, pi_time_1))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time_2 / 2, pi_time_2))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 5 * tau - pi_time_2 / 2, pi_time_2))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 6 * tau - pi_time_1 / 2, pi_time_1))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 7 * tau - pi_time_2 / 2, pi_time_2))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 8 * tau - pi_time_1 / 2, pi_time_1))
                    section_begin_time += 8 * tau

                # the second 3*pi/2 pulse and readout
                pulse_sequence.append(
                    Pulse(microwave_channel_1, section_begin_time + tau / 2, three_pi_half_time_1))
                pulse_sequence.append(
                    Pulse('laser', section_begin_time + tau / 2 + three_pi_half_time_1 + delay_mw_readout,
                          nv_reset_time))
                pulse_sequence.append(
                    Pulse('apd_readout',
                          section_begin_time + tau / 2 + three_pi_half_time_1 + delay_mw_readout + delay_readout,
                          meas_time))

                start_of_second_HE = section_begin_time + tau / 2 + three_pi_half_time_1 + delay_mw_readout + nv_reset_time + laser_off_time

                # the first pi/2 pulse
                pulse_sequence.append(
                    Pulse(microwave_channel_1, start_of_second_HE, pi_half_time_1))
                section_begin_time = start_of_second_HE + pi_half_time_1 - tau / 2  # for the first pulse, only wait tau/2
                for i in range(0, number_of_pulse_blocks):

                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time_1 / 2, pi_time_1))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 2 * tau - pi_time_2 / 2, pi_time_2))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time_1 / 2, pi_time_1))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time_2 / 2, pi_time_2))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 5 * tau - pi_time_2 / 2, pi_time_2))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 6 * tau - pi_time_1 / 2, pi_time_1))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 7 * tau - pi_time_2 / 2, pi_time_2))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 8 * tau - pi_time_1 / 2, pi_time_1))
                    section_begin_time += 8 * tau

                # the second pi/2 pulse and readout
                pulse_sequence.append(
                    Pulse(microwave_channel_1, section_begin_time + tau / 2, pi_half_time_1))
                pulse_sequence.append(
                    Pulse('laser',
                          section_begin_time + tau / 2 + pi_half_time_1 + delay_mw_readout,
                          nv_reset_time))
                pulse_sequence.append(
                    Pulse('apd_readout',
                          section_begin_time + tau / 2 + pi_half_time_1 + delay_mw_readout + delay_readout,
                          meas_time))

                # E-Sensing SEQUENCE
                if True:
                    start_of_DEER = section_begin_time + tau / 2 + pi_half_time_1 + delay_mw_readout + nv_reset_time + laser_off_time
                    # the first pi/2 pulse
                    pulse_sequence.append(Pulse(microwave_channel_1, start_of_DEER, pi_half_time_1))
                    pulse_sequence.append(Pulse('fg_trigger', start_of_DEER - trigger_latency, trigger_time))

                    section_begin_time = start_of_DEER + pi_half_time_1 - tau / 2  # for the first pulse, only wait tau/2
                    for i in range(0, number_of_pulse_blocks):
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time_1 / 2, pi_time_1))
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 2 * tau - pi_time_2 / 2, pi_time_2))
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time_1 / 2, pi_time_1))
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time_2 / 2, pi_time_2))
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 5 * tau - pi_time_2 / 2, pi_time_2))
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 6 * tau - pi_time_1 / 2, pi_time_1))
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 7 * tau - pi_time_2 / 2, pi_time_2))
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 8 * tau - pi_time_1 / 2, pi_time_1))
                        section_begin_time += 8 * tau

                    # the second 3*pi/2 pulse and readout
                    pulse_sequence.append(
                        Pulse(last_mw_channel, section_begin_time + tau / 2, last_three_pi_half_time))
                    pulse_sequence.append(
                        Pulse('laser', section_begin_time + tau / 2 + last_three_pi_half_time + delay_mw_readout,
                              nv_reset_time))
                    pulse_sequence.append(
                        Pulse('apd_readout',
                              section_begin_time + tau / 2 + last_three_pi_half_time + delay_mw_readout + delay_readout,
                              meas_time))

                    start_of_second_DEER = section_begin_time + tau / 2 + last_three_pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time

                    # the first pi/2 pulse
                    pulse_sequence.append(Pulse(microwave_channel_1, start_of_second_DEER, pi_half_time_1))
                    pulse_sequence.append(Pulse('fg_trigger', start_of_second_DEER - trigger_latency, trigger_time))
                    section_begin_time = start_of_second_DEER + pi_half_time_1 - tau / 2  # for the first pulse, only wait tau/2

                    for i in range(0, number_of_pulse_blocks):
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time_1 / 2, pi_time_1))
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 2 * tau - pi_time_2 / 2, pi_time_2))
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time_1 / 2, pi_time_1))
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time_2 / 2, pi_time_2))
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 5 * tau - pi_time_2 / 2, pi_time_2))
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 6 * tau - pi_time_1 / 2, pi_time_1))
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 7 * tau - pi_time_2 / 2, pi_time_2))
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 8 * tau - pi_time_1 / 2, pi_time_1))
                        section_begin_time += 8 * tau

                    # the second pi/2 pulse and readout
                    pulse_sequence.append(
                        Pulse(last_mw_channel, section_begin_time + tau / 2, last_pi_half_time))
                    pulse_sequence.append(
                        Pulse('laser',
                              section_begin_time + tau / 2 + last_pi_half_time + delay_mw_readout,
                              nv_reset_time))
                    pulse_sequence.append(
                        Pulse('apd_readout',
                              section_begin_time + tau / 2 + last_pi_half_time + delay_mw_readout + delay_readout,
                              meas_time))

                pulse_sequence.append(Pulse('apd_switch', 0,
                                            section_begin_time + tau / 2 + last_pi_half_time + delay_mw_readout + nv_reset_time))
                pulse_sequences.append(pulse_sequence)

            print('number of sequences before validation ', len(pulse_sequences))
            # return pulse_sequences, self.settings['num_averages'], tau_list, meas_time
            return pulse_sequences, amp_list, meas_time

        elif self.settings['decoupling_seq']['type'] == 'CPMG':

            for amp in amp_list:
                tau = tau_total / (1 * self.settings['decoupling_seq']['num_of_pulse_blocks'])
                pulse_sequence = []

                # ECHO SEQUENCE:
                start_of_first_HE = laser_off_time

                # the first pi/2 pulse
                pulse_sequence = [Pulse(microwave_channel_1, start_of_first_HE, pi_half_time_1)]
                section_begin_time = start_of_first_HE + pi_half_time_1 - tau / 2  # for the first pulse, only wait tau/2

                for i in range(0, number_of_pulse_blocks):
                    pulse_sequence.append(Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time_2 / 2, pi_time_2))
                    section_begin_time += 1 * tau

                # the second 3*pi/2 pulse and readout
                pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + tau / 2, three_pi_half_time_1))
                pulse_sequence.append(Pulse('laser', section_begin_time + tau / 2 + three_pi_half_time_1 + delay_mw_readout,
                                             nv_reset_time))
                pulse_sequence.append(Pulse('apd_readout',
                                             section_begin_time + tau / 2 + three_pi_half_time_1 + delay_mw_readout + delay_readout,
                                             meas_time))

                start_of_second_HE = section_begin_time + tau / 2 + three_pi_half_time_1 + delay_mw_readout + nv_reset_time + laser_off_time

                # the first pi/2 pulse
                pulse_sequence.append(Pulse(microwave_channel_1, start_of_second_HE, pi_half_time_1))

                section_begin_time = start_of_second_HE + pi_half_time_1 - tau / 2  # for the first pulse, only wait tau/2

                for i in range(0, number_of_pulse_blocks):
                    pulse_sequence.append(Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time_2 / 2, pi_time_2))
                    section_begin_time += 1 * tau

                # the second pi/2 pulse and readout
                pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + tau / 2, pi_half_time_1))
                pulse_sequence.append(Pulse('laser',
                                             section_begin_time + tau / 2 + pi_half_time_1 + delay_mw_readout,
                                             nv_reset_time))
                pulse_sequence.append(Pulse('apd_readout',
                                             section_begin_time + tau / 2 + pi_half_time_1 + delay_mw_readout + delay_readout,
                                             meas_time))

                # E-Sensing SEQUENCE
                if True:
                    start_of_DEER = section_begin_time + tau / 2 + pi_half_time_1 + delay_mw_readout + nv_reset_time + laser_off_time
                    # the first pi/2 pulse
                    pulse_sequence.append(Pulse(microwave_channel_1, start_of_DEER, pi_half_time_1))
                    pulse_sequence.append(Pulse('fg_trigger', start_of_DEER - trigger_latency, trigger_time))
                    section_begin_time = start_of_DEER + pi_half_time_1 - tau / 2  # for the first pulse, only wait tau/2
                    for i in range(0, number_of_pulse_blocks):
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time_2 / 2, pi_time_2))
                        section_begin_time += 1 * tau

                    # the second 3*pi/2 pulse and readout
                    pulse_sequence.append(Pulse(last_mw_channel, section_begin_time + tau / 2, last_three_pi_half_time))
                    pulse_sequence.append(Pulse('laser', section_begin_time + tau / 2 + last_three_pi_half_time + delay_mw_readout,
                                                nv_reset_time))
                    pulse_sequence.append(Pulse('apd_readout',
                                                section_begin_time + tau / 2 + last_three_pi_half_time + delay_mw_readout + delay_readout,
                                                meas_time))

                    start_of_second_DEER = section_begin_time + tau / 2 + last_three_pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time

                    # the first pi/2 pulse
                    pulse_sequence.append(Pulse(microwave_channel_1, start_of_second_DEER, pi_half_time_1))
                    pulse_sequence.append(Pulse('fg_trigger', start_of_second_DEER - trigger_latency, trigger_time))
                    section_begin_time = start_of_second_DEER + pi_half_time_1 - tau / 2  # for the first pulse, only wait tau/2

                    for i in range(0, number_of_pulse_blocks):
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time_2 / 2, pi_time_2))
                        section_begin_time += 1 * tau


                    # the second pi/2 pulse and readout
                    pulse_sequence.append(Pulse(last_mw_channel, section_begin_time + tau / 2, last_pi_half_time))
                    pulse_sequence.append(
                        Pulse('laser', section_begin_time + tau / 2 + last_pi_half_time + delay_mw_readout,
                              nv_reset_time))
                    pulse_sequence.append(Pulse('apd_readout',
                                                section_begin_time + tau / 2 + last_pi_half_time + delay_mw_readout + delay_readout,
                                                meas_time))

                # Turn on APD switch all the time
                pulse_sequence.append(Pulse('apd_switch', 0,
                                            section_begin_time + tau / 2 + last_pi_half_time + delay_mw_readout + nv_reset_time))


                pulse_sequences.append(pulse_sequence)

            print('number of sequences before validation ', len(pulse_sequences))
            # return pulse_sequences, self.settings['num_averages'], tau_list, meas_time
            return pulse_sequences, amp_list, meas_time

    def _run_sweep(self, pulse_sequences, num_loops_sweep, num_daq_reads, tau_list,avg_num = 0, verbose=False):
        """
        Each pulse sequence specified in pulse_sequences is run num_loops_sweep consecutive times.

        Args:
            pulse_sequences: a list of pulse sequences to run, each corresponding to a different value of tau. Each
                             sequence is a list of Pulse objects specifying a given pulse sequence
            num_loops_sweep: number of times to repeat each sequence before moving on to the next one
            num_daq_reads: number of times the daq must read for each sequence (generally 1, 2, or 3)

        Poststate: self.data['counts'] is updated with the acquired data

        """
        # ER 20180731 set init_fluor to zero
     #   self.data['init_fluor'] = 0.

        # print('start _run_sweep now')


        rand_indexes = list(range(len(pulse_sequences)))
        if self.settings['randomize']:
            random.shuffle(rand_indexes)

        if verbose:
            print(('_run_sweep number of pulse sequences', len(pulse_sequences)))

        for index, sequence in enumerate(pulse_sequences):
            # print('current index')
            # print(index)
            if verbose:
                print(('_run_sweep index', index, len(pulse_sequences)))

            rand_index = rand_indexes[index]
            self.instruments['fg']['instance'].update({'amplitude': float(tau_list[rand_index])})

            # self.current_rand_index = rand_index
            if self._abort:
                # print('aborted')
                self.instruments['PB']['instance'].update({'microwave_switch': {'status': False}})
                # print('mw is off')
                break

            # print('running single sqeucnce')
            result = self._run_single_sequence(pulse_sequences[rand_index], num_loops_sweep, num_daq_reads)  # keep entire array
            # print('single sqeucnce done')
            self.count_data[rand_index] = self.count_data[rand_index] + result
            # print('before counts_to_check')

            counts_to_check = self._normalize_to_kCounts(np.array(result), self.measurement_gate_width, num_loops_sweep)
            self.data['counts'][rand_index] = self._normalize_to_kCounts(self.count_data[rand_index], self.measurement_gate_width,
                                                                    self.current_averages)

            self.count_sqrd_sum[rand_index] = self.count_sqrd_sum[rand_index] + np.array(
                counts_to_check) ** 2  # unit: kcps^2
            # counts_to_check maybe a list
            self.data['counts_err'][rand_index] = self._cal_err(self.count_sqrd_sum[rand_index],
                                                                self.data['counts'][rand_index], avg_num)

            # self.sequence_index = index
            self.sequence_index = rand_index
            counts_temp = counts_to_check[0]
            # track to the NV if necessary ER 5/31/17
            if self.settings['Tracking']['on/off']:
                print('doing tracking')
                print('counts_temp:', counts_temp)
                if (1+(1-self.settings['Tracking']['threshold']))*self.settings['Tracking']['init_fluor'] < counts_temp or \
                        self.settings['Tracking']['threshold']*self.settings['Tracking']['init_fluor'] > counts_temp:
                    print('(2-threshold)*intial_fluor:',
                          (1 + (1 - self.settings['Tracking']['threshold'])) * self.settings['Tracking']['init_fluor'])
                    print('threshold*intial_fluor:', self.settings['Tracking']['threshold']*self.settings['Tracking']['init_fluor'])
                    if verbose:
                        print('TRACKING TO NV...')
                    self.scripts['find_nv'].run()
                    self.scripts['find_nv'].settings['initial_point'] = self.scripts['find_nv'].data['maximum_point']
                    print('find_nv done')
                    if self.scripts['find_nv'].data[
                        'fluorescence'] == 0.0:  # if it doesn't find an NV, abort the experiment
                        print('Could not find an NV in FindNV.')
                        self.log('Could not find an NV in FindNV.')
                        self._abort = True
                        return  # exit function in case no NV is found
                    else:
                        print('Start optimize_z.')
                        self.scripts['optimize_z'].run()
            self.updateProgress.emit(self._calc_progress(index))

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

        # if 'fits' in data.keys() is not None and data['fits'] is not None:
        #
        #     # counts = (-data['counts'][:,1] + data['counts'][:,0])/ (data['counts'][:,0] + data['counts'][:,1])
        #     tau = data['tau']
        #     fits = data['fits']
        #
        #     axislist[0].plot(tau, data['norm_echo'])
        #     tauinterp = np.linspace(np.min(tau), np.max(tau), 100)
        #     axislist[0].plot(tauinterp, exp_offset(tauinterp, fits[0], fits[1], fits[2]))
        #
        #
        #     if 'fits_esig' in data.keys() is not None and data['fits_esig'] is not None:
        #         fits_esig = data['fits_esig']
        #         axislist[0].plot(tau, data['norm_esig'])
        #         axislist[0].plot(tauinterp, exp_offset(tauinterp, fits_esig[0], fits_esig[1], fits_esig[2]))
        #         axislist[0].legend(labels=('Echo', 'Exp Fit','E-Sensing', 'Exp Fit'), fontsize=10)
        #         axislist[0].set_title(
        #             'Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
        #                 'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s), {:s} electrometry\n averaged over {:d} blocks \n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz, \n T2 (p = 1) = {:2.1f} ns, T2 (E-field) = {:2.1f} ns'.format(
        #                 self.settings['decoupling_seq']['type'],
        #                 self.settings['decoupling_seq']['num_of_pulse_blocks'],
        #                 self.settings['decoupling_seq']['eSensing_type'], avg_block_number,
        #                 self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency'] * 1e-9,
        #                 fits[1], fits_esig[1]))
        #     else:
        #         axislist[0].legend(labels=('Echo', 'Exp Fit'), fontsize=10)
        #         axislist[0].set_title(
        #             'Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
        #                 'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s)\n averaged over {:d} blocks \n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz, \n T2 (p = 1) = {:2.1f} ns'.format(
        #                 self.settings['decoupling_seq']['type'],
        #                 self.settings['decoupling_seq']['num_of_pulse_blocks'], avg_block_number,
        #                 self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency'] * 1e-9,
        #                 fits[1]))
        #     axislist[0].set_ylabel('fluorescence contrast')
        #     axislist[0].set_xlabel('tau [ns]')

        if 'exp_finished' in data.keys() is not None and data['exp_finished'] == 1:
            if 'norm_echo' in data.keys() is not None and data['norm_echo'] is not None:
                tau = data['tau']
                axislist[0].plot(tau, data['norm_echo'])
                if 'norm_esig' in data.keys() is not None and data['norm_esig'] is not None:
                    axislist[0].plot(tau, data['norm_esig'])
                    axislist[0].legend(labels=('PDD', 'E-Sensing'), fontsize=10)
                    axislist[0].set_title(
                        '(final plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
                            'microwave_channel'] + '\nPDD: {:s} {:d} block(s), {:s} electrometry\n after {:d} averages \n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz \n tau = {:.2f} us'.format(
                            self.settings['decoupling_seq']['type'],
                            self.settings['decoupling_seq']['num_of_pulse_blocks'],
                            self.settings['decoupling_seq']['eSensing_type'], avg_block_number,
                            self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency'] * 1e-9, self.settings['tau_time'] * 1e-3))
                else:
                    axislist[0].legend(labels=('PDD'), fontsize=10)
                    axislist[0].set_title(
                        '(final plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
                            'microwave_channel'] + '\nPDD: {:s} {:d} block(s)\n after {:d} averages \n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz \n tau = {:.2f} us'.format(
                            self.settings['decoupling_seq']['type'],
                            self.settings['decoupling_seq']['num_of_pulse_blocks'],
                            avg_block_number,
                            self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency'] * 1e-9, self.settings['tau_time'] * 1e-3))
                axislist[0].set_ylabel('fluorescence contrast')
                axislist[0].set_xlabel('amplitude [Vpp]')
            else:
                axislist[0].plot(data['tau'], data['counts'][:, 0])
                axislist[0].plot(data['tau'], data['counts'][:, 1])
                axislist[0].plot(data['tau'], data['counts'][:, 2])
                axislist[0].plot(data['tau'], data['counts'][:, 3])
                axislist[0].legend(labels=('PDD up {:.0f}kcps'.format(np.mean(data['counts'][:, 0])),
                                           'PDD down {:.0f}kcps'.format(np.mean(data['counts'][:, 1])),
                                           'E-Sensing up {:.0f}kcps'.format(np.mean(data['counts'][:, 2])),
                                           'E-Sensing down {:.0f}kcps'.format(np.mean(data['counts'][:, 3]))),
                                   fontsize=10)
                axislist[0].set_title(
                    '(final plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
                        'microwave_channel'] + '\nPDD: {:s} {:d} block(s), {:s} electrometry\n after {:d} averages\n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz \n tau = {:.2f} us'.format(
                        self.settings['decoupling_seq']['type'],
                        self.settings['decoupling_seq']['num_of_pulse_blocks'],
                        self.settings['decoupling_seq']['eSensing_type'],
                        avg_block_number, self.settings['mw_pulses']['mw_power'],
                        self.settings['mw_pulses']['mw_frequency'] * 1e-9, self.settings['tau_time'] * 1e-3))


                axislist[0].set_ylabel('fluorescence [kcps]')
                axislist[0].set_xlabel('amplitude [Vpp]')

        else:
            super(eSensing_swpV, self)._plot(axislist)
            axislist[0].set_xlabel('amplitude [Vpp]')
            axislist[0].legend(labels=('PDD up {:.0f}kcps'.format(np.mean(data['counts'][:, 0])),
                                       'PDD down {:.0f}kcps'.format(np.mean(data['counts'][:, 1])),
                                       'E-Sensing up {:.0f}kcps'.format(np.mean(data['counts'][:, 2])),
                                       'E-Sensing down {:.0f}kcps'.format(np.mean(data['counts'][:, 3]))),
                               fontsize=10)
            axislist[0].set_title(
                '(initial plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
                    'microwave_channel'] + '\nPDD: {:s} {:d} block(s), {:s} electrometry\n after {:d} averages\n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz \n tau = {:.2f} us'.format(
                    self.settings['decoupling_seq']['type'], self.settings['decoupling_seq']['num_of_pulse_blocks'],
                    self.settings['decoupling_seq']['eSensing_type'],
                    avg_block_number, self.settings['mw_pulses']['mw_power'],
                    self.settings['mw_pulses']['mw_frequency'] * 1e-9, self.settings['tau_time'] * 1e-3))

            # if self.settings['decoupling_seq']['do_eSensing']:
            #     axislist[0].legend(labels=('Echo up {:.0f}kcps'.format(np.mean(data['counts'][:, 0])),
            #                                'Echo down {:.0f}kcps'.format(np.mean(data['counts'][:, 1])),
            #                                'E-Sensing up {:.0f}kcps'.format(np.mean(data['counts'][:, 0])),
            #                                'E-Sensing down {:.0f}kcps'.format(np.mean(data['counts'][:, 1]))),
            #                        fontsize=10)
            #     axislist[0].set_title(
            #         '(initial plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
            #             'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s), {:s} electrometry\n averaged over {:d} blocks\n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz'.format(
            #             self.settings['decoupling_seq']['type'], self.settings['decoupling_seq']['num_of_pulse_blocks'],
            #             self.settings['decoupling_seq']['eSensing_type'],
            #             avg_block_number, self.settings['mw_pulses']['mw_power'],
            #             self.settings['mw_pulses']['mw_frequency'] * 1e-9))
            #
            # else:
            #     axislist[0].legend(labels=('Echo up {:.0f}kcps'.format(np.mean(data['counts'][:, 0])),
            #                                'Echo down {:.0f}kcps'.format(np.mean(data['counts'][:, 1]))), fontsize=10)
            #     axislist[0].set_title(
            #         '(initial plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
            #             'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s) \n averaged over {:d} blocks\n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz'.format(
            #             self.settings['decoupling_seq']['type'], self.settings['decoupling_seq']['num_of_pulse_blocks'],
            #             avg_block_number, self.settings['mw_pulses']['mw_power'],
            #             self.settings['mw_pulses']['mw_frequency'] * 1e-9))

    def _update_plot(self, axislist):

        if len(axislist[0].lines) == 0:
            self._plot(axislist)
            return

        super(eSensing_swpV, self)._update_plot(axislist)
        axislist[0].legend(labels=('PDD up {:.0f}kcps'.format(np.mean(self.data['counts'][:, 0])),
                                   'PDD down {:.0f}kcps'.format(np.mean(self.data['counts'][:, 1])),
                                   'E-Sensing up {:.0f}kcps'.format(np.mean(self.data['counts'][:, 2])),
                                   'E-Sensing down {:.0f}kcps'.format(np.mean(self.data['counts'][:, 3]))),
                           fontsize=10)
        axislist[0].set_title(
            '(updating plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
                'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s), {:s} electrometry\n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz \n tau = {:.2f} us'.format(
                self.settings['decoupling_seq']['type'], self.settings['decoupling_seq']['num_of_pulse_blocks'],
                self.settings['decoupling_seq']['eSensing_type'],
                self.settings['mw_pulses']['mw_power'],
                self.settings['mw_pulses']['mw_frequency'] * 1e-9, self.settings['tau_time'] * 1e-3))

        # if self.settings['decoupling_seq']['do_eSensing']:
        #     axislist[0].legend(labels=('Echo up {:.0f}kcps'.format(np.mean(self.data['counts'][:, 0])),
        #                                'Echo down {:.0f}kcps'.format(np.mean(self.data['counts'][:, 1])),
        #                                'E-Sensing up {:.0f}kcps'.format(np.mean(self.data['counts'][:, 0])),
        #                                'E-Sensing down {:.0f}kcps'.format(np.mean(self.data['counts'][:, 1]))),
        #                        fontsize=10)
        #     axislist[0].set_title(
        #         '(updating plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
        #             'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s), {:s} electrometry\n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz'.format(
        #             self.settings['decoupling_seq']['type'], self.settings['decoupling_seq']['num_of_pulse_blocks'],
        #             self.settings['decoupling_seq']['eSensing_type'],
        #             self.settings['mw_pulses']['mw_power'],
        #             self.settings['mw_pulses']['mw_frequency'] * 1e-9))
        #
        # else:
        #     axislist[0].legend(labels=('Echo up {:.0f}kcps'.format(np.mean(self.data['counts'][:, 0])),
        #                                'Echo down {:.0f}kcps'.format(np.mean(self.data['counts'][:, 1]))), fontsize=10)
        #
        #     axislist[0].set_title('(updating plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
        #         'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s) \n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz'.format(
        #         self.settings['decoupling_seq']['type'], self.settings['decoupling_seq']['num_of_pulse_blocks'],
        #         self.settings['mw_pulses']['mw_power'],
        #         self.settings['mw_pulses']['mw_frequency'] * 1e-9))

    def _plot_validate(self, axes_list):
        """
        Preview pulse sequence by plotting first and last sequence to plots 1 and 0, respectively
        Args:
            axes_list: List containing axes to plot to

        Returns:

        """
        print('printing validated pulses')
        axis0 = axes_list[0]
        axis1 = axes_list[1]
        axis0.clear()
        axis1.clear()

        pulse_sequences, tau_list, _ = self.create_pulse_sequences(logging=False)
        # print('here')

      #  if pulse_sequences[0]:
        if pulse_sequences:
            plot_pulses(axis0, pulse_sequences[0])
            axis0.set_title('Pulse Visualization for Minimum Amplitude = {:f} Vpp'.format(tau_list[0]))
            plot_pulses(axis1, pulse_sequences[-1])
            axis1.set_title('Pulse Visualization for Maximum Amplitude = {:f} Vpp'.format(tau_list[-1]))
        else:
            print('no pulse sequences passed validation!!!')

# The following script is not fully working yet. One pitfall is that at short tau, the function generator might not allow applying too many pi pulses
class eSensing_swpTau(PulsedExperimentBaseScript):
    """
    This script does AC electrometry using NV. Amplitude Vpp is fixed, and tau is swept.
    Voltages are applied using Agilent 33120A.
    MW pulses are applied using Keysight N9310A.

    --> Ziwei Qiu 6/26/2019

    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_power', -10.0, float, 'microwave power in dB'),
            Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
            Parameter('microwave_channel', 'i', ['i', 'q'], 'Channel to use for the first pi/2 mw pulse'),
            Parameter('pi_pulse_time', 50.0, float, 'time duration of a pi pulse (in ns) for channel i'),
            Parameter('pi_half_pulse_time', 25.0, float, 'time duration of a pi/2 pulse (in ns) for channel i'),
            Parameter('3pi_half_pulse_time', 75.0, float, 'time duration of a 3pi/2 pulse (in ns) for channel i'),
            Parameter('Q_same_as_I', True, bool, 'check if I and Q are symmetric'),
            Parameter('pi_pulse_time_q', 50.0, float, 'time duration of a pi pulse (in ns) for channel q'),
            Parameter('pi_half_pulse_time_q', 25.0, float, 'time duration of a pi/2 pulse (in ns) for channel q'),
            Parameter('3pi_half_pulse_time_q', 75.0, float, 'time duration of a 3pi/2 pulse (in ns) for channel q')
        ]),
        Parameter('tau_times', [
            Parameter('min_time', 200, float, 'minimum time between pi pulses'),
            Parameter('max_time', 10000, float, 'maximum time between pi pulses'),
            Parameter('time_step', 100.,
                      [2.5, 5., 10., 20., 50., 100., 200., 300., 400., 500., 600., 800., 1000., 2000., 5000., 8000.,
                       10000., 50000., 100000., 500000.],
                      'time step increment of time between pi pulses (in ns)')
        ]),
        Parameter('e_fields', [
            Parameter('amp', 10, float, '[Vpp] minimum amplitude'),
            Parameter('offset', 0.0, float, '[Vdc] DC offset voltage'),
            Parameter('wave_shape', 'SQUare', ['SINusoid', 'SQUare'], 'choose the wave shape'),
            Parameter('trigger_latency', 1253.0, float, '[ns] Agilent33120A  trigger latency'),
            Parameter('trigger_time', 500.0, float, '[ns] trigger pulse time (>100ns)'),
        ]),
        Parameter('decoupling_seq', [
            Parameter('type', 'XY4', ['spin_echo', 'CPMG', 'XY4', 'XY8'], 'type of dynamical decoupling sequences'),
            Parameter('num_of_pulse_blocks', 1, int, 'number of pulse blocks.'),
            Parameter('eSensing_type', 'sine', ['sine', 'cosine'], 'choose to do sine or cosine electrometry'),
        ]),
        # Parameter('read_out', [
        #     Parameter('meas_time', 200, float, 'measurement time after rabi sequence (in ns)'),
        #     Parameter('nv_reset_time', 2000, int, 'time with laser on to reset state'),
        #     Parameter('laser_off_time', 1000, int,
        #               'minimum laser off time before taking measurements (ns)'),
        #     Parameter('delay_mw_readout', 600, int, 'delay between mw and readout (in ns)'),
        #     Parameter('delay_readout', 550, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        # ]),
        Parameter('min_pulse_dur', 15, [10, 15, 50], 'Minimum allowed pulse duration (ns)'),
        Parameter('num_averages', 500000, int, 'number of averages (>100000)'),
    ]

    _INSTRUMENTS = {'daq': NI6353, 'PB': LISE607RTPulseBlaster, 'mw_gen_a': AgilentMicrowaveGenerator, 'fg': Agilent33120A}

    def _function(self):

        self.data['fits'] = None

        ## Turn off green light (the pulse blaster will pulse it on when needed)
        self.instruments['PB']['instance'].update({'laser': {'status': False}})
        self.instruments['PB']['instance'].update({'apd_switch': {'status': False}})

        # set up the RF generator
        self.instruments['mw_gen_a']['instance'].update({'enable_IQ': True})
        self.instruments['mw_gen_a']['instance'].update({'enable_modulation': True})
        self.instruments['mw_gen_a']['instance'].update({'freq_mode': 'CW'})
        self.instruments['mw_gen_a']['instance'].update({'power_mode': 'CW'})
        self.instruments['mw_gen_a']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power']})
        self.instruments['mw_gen_a']['instance'].update({'frequency': self.settings['mw_pulses']['mw_frequency']})
        self.instruments['mw_gen_a']['instance'].update({'enable_output': True})

        # set up the function generator
        self.instruments['fg']['instance'].update({'display_on': False})
        self.instruments['fg']['instance'].update({'burst_mod': True})
        self.instruments['fg']['instance'].update({'amplitude': self.settings['e_fields']['amp']})
        self.instruments['fg']['instance'].update({'offset': self.settings['e_fields']['offset']})
        self.instruments['fg']['instance'].update({'wave_shape': self.settings['e_fields']['wave_shape']})
        valid = True
        if self.settings['decoupling_seq']['type'] == 'spin_echo' or self.settings['decoupling_seq']['type'] == 'CPMG':
            tau = self.settings['tau_time'] / (1 * self.settings['decoupling_seq']['num_of_pulse_blocks'])
            if self.settings['decoupling_seq']['num_of_pulse_blocks'] == 1:
                fg_freq = 1E9 / (1 * tau)
                burst_cnt = 1
                burst_phase = 0.0
                if (fg_freq / 1E6) / burst_cnt > 1:
                    self.log('Invalid setting: (fg_freq / 1E6) / burst_cnt > 1.')
                    valid = False
                else:
                    self.instruments['fg']['instance'].update({'burst_phase': burst_phase})
                    self.instruments['fg']['instance'].update({'frequency': fg_freq})
                    self.instruments['fg']['instance'].update({'burst_count': burst_cnt})
            elif self.settings['decoupling_seq']['num_of_pulse_blocks']%2 == 0:
                fg_freq = 1E9 / (2 * tau)
                burst_cnt = self.settings['decoupling_seq']['num_of_pulse_blocks'] / 2
                burst_phase = 90.0
                if (fg_freq / 1E6) / burst_cnt > 1:
                    self.log('Invalid setting: (fg_freq / 1E6) / burst_cnt > 1.')
                    valid = False
                else:
                    self.instruments['fg']['instance'].update({'burst_phase': burst_phase})
                    self.instruments['fg']['instance'].update({'frequency': fg_freq})
                    self.instruments['fg']['instance'].update({'burst_count': burst_cnt})
            else:
                self.log('Invalid setting: num_of_pulse_blocks needs to be 1 or even number if using spin-echo or CPMG.')
                valid = False
        elif self.settings['decoupling_seq']['type'] == 'XY4':
            tau = self.settings['tau_time'] / (4 * self.settings['decoupling_seq']['num_of_pulse_blocks'])
            fg_freq = 1E9 / (2 * tau)
            burst_cnt = 2 * self.settings['decoupling_seq']['num_of_pulse_blocks']
            burst_phase = 90.0
            if (fg_freq / 1E6) / burst_cnt > 1:
                self.log('Invalid setting: (fg_freq / 1E6) / burst_cnt > 1.')
                valid = False
            else:
                self.instruments['fg']['instance'].update({'burst_phase': burst_phase})
                self.instruments['fg']['instance'].update({'frequency': fg_freq})
                self.instruments['fg']['instance'].update({'burst_count': burst_cnt})
        elif self.settings['decoupling_seq']['type'] == 'XY8':
            tau = self.settings['tau_time'] / (8 * self.settings['decoupling_seq']['num_of_pulse_blocks'])
            fg_freq = 1E9 / (2 * tau)
            burst_cnt = 4 * self.settings['decoupling_seq']['num_of_pulse_blocks']
            burst_phase = 90.0
            if (fg_freq / 1E6) / burst_cnt > 1:
                self.log('Invalid setting: (fg_freq / 1E6) / burst_cnt > 1.')
                valid = False
            else:
                self.instruments['fg']['instance'].update({'burst_phase': burst_phase})
                self.instruments['fg']['instance'].update({'frequency': fg_freq})
                self.instruments['fg']['instance'].update({'burst_count': burst_cnt})
        else:
            self.log('Invalid setting: choose the right decoupling_seq type')
            valid = False

        if valid:
            self.avg_block_number = super(eSensing_swpTau, self)._function()

            # this is to flag when the experiment is done or stopped and this allows proper saving even fit fails
            self.data['exp_finished'] = 1
            self.data['avg_block_number'] = self.avg_block_number

        # turn off laser, apd_switch and MW (just in case)
        self.instruments['PB']['instance'].update({'laser': {'status': False}})
        self.instruments['PB']['instance'].update({'apd_switch': {'status': False}})
        self.instruments['mw_gen_a']['instance'].update({'enable_output': False})
        self.instruments['mw_gen_a']['instance'].update({'enable_IQ': False})
        self.instruments['mw_gen_a']['instance'].update({'enable_modulation': False})
        self.instruments['fg']['instance'].update({'display_on': True})

        # no fitting for now
        if 'counts' in self.data.keys() and 'tau' in self.data.keys() and not 0 in (
                self.data['counts'][:, 1] + self.data['counts'][:, 0]) and not 0 in (
                self.data['counts'][:, 3] + self.data['counts'][:, 2]):

            self.data['norm_echo'] = 2. * (- self.data['counts'][:, 1] + self.data['counts'][:, 0]) / (
                    self.data['counts'][:, 1] + self.data['counts'][:, 0])
            self.data['norm_esig'] = 2. * (- self.data['counts'][:, 3] + self.data['counts'][:, 2]) / (
                    self.data['counts'][:, 3] + self.data['counts'][:, 2])

            tau = self.data['tau']
            counts = self.data['norm_echo']
            esig_counts = self.data['norm_esig']
            fit_success = False

            # try:
            #     fits = fit_exp_decay(tau, counts, offset = True, verbose = True)
            #     self.data['fits'] = fits
            #     fit_success = True
            # except:
            #     self.data['fits'] = None
            #     self.log('T2 fit failed')
            #     fit_success = False
            #
            # if self.settings['decoupling_seq']['do_eSensing'] and not 0 in (self.data['counts'][:, 3] + self.data['counts'][:, 2]):
            #     self.data['norm_esig'] = 2. * (- self.data['counts'][:, 3] + self.data['counts'][:, 2]) / (
            #                 self.data['counts'][:, 3] + self.data['counts'][:, 2])
            #
            #     esig_counts = self.data['norm_esig']
                # try:
                #     fits_esig = fit_exp_decay(tau, esig_counts, offset = True, verbose = True)
                #     self.data['fits_esig'] = fits_esig
                #     fit_esig_success = True
                # except:
                #     self.data['fits_esig'] = None
                #     self.log('E-sensing T2 fit failed')
                #     fit_esig_success = False

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
        max_range = int(np.floor((self.settings['e_fields']['max_amp'] - self.settings['e_fields']['min_amp']) /
                                 self.settings['e_fields']['amp_step']))
        amp_list = np.array(
            [self.settings['e_fields']['min_amp'] + i * self.settings['e_fields']['amp_step'] for i in
             range(max_range)])
        # ignore the sequence if the amplitude is smaller than 15Vpp
        amp_list = [x for x in amp_list if x >= 0.1]

        print('Vpp_list:', amp_list)

        nv_reset_time = self.settings['read_out']['nv_reset_time']
        delay_readout = self.settings['read_out']['delay_readout']

        pi_time = self.settings['mw_pulses']['pi_pulse_time']
        pi_half_time = self.settings['mw_pulses']['pi_half_pulse_time']
        three_pi_half_time = self.settings['mw_pulses']['3pi_half_pulse_time']
        if self.settings['mw_pulses']['Q_same_as_I']:
            pi_time_q = pi_time
            pi_half_time_q = pi_half_time
            three_pi_half_time_q = three_pi_half_time
        else:
            pi_time_q = self.settings['mw_pulses']['pi_pulse_time']
            pi_half_time_q = self.settings['mw_pulses']['pi_half_pulse_time']
            three_pi_half_time_q = self.settings['mw_pulses']['3pi_half_pulse_time']

        if self.settings['mw_pulses']['microwave_channel'] == 'i':
            microwave_channel_1 = 'microwave_i'
            microwave_channel_2 = 'microwave_q'
            pi_time_1 = pi_time
            pi_half_time_1 = pi_half_time
            three_pi_half_time_1 = three_pi_half_time
            pi_time_2 = pi_time_q
            pi_half_time_2 = pi_half_time_q
            three_pi_half_time_2 = three_pi_half_time_q
        else:
            microwave_channel_1 = 'microwave_q'
            microwave_channel_2 = 'microwave_i'
            pi_time_1 = pi_time_q
            pi_half_time_1 = pi_half_time_q
            three_pi_half_time_1 = three_pi_half_time_q
            pi_time_2 = pi_time
            pi_half_time_2 = pi_half_time
            three_pi_half_time_2 = three_pi_half_time

        if self.settings['decoupling_seq']['eSensing_type'] == 'sine':
            last_mw_channel = microwave_channel_2
            last_pi_time = pi_time_2
            last_pi_half_time = pi_half_time_2
            last_three_pi_half_time = three_pi_half_time_2
        else:
            last_mw_channel = microwave_channel_1
            last_pi_time = pi_time_1
            last_pi_half_time = pi_half_time_1
            last_three_pi_half_time = three_pi_half_time_1


        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        trigger_latency = self.settings['e_fields']['trigger_latency']
        trigger_time = self.settings['e_fields']['trigger_time']

        number_of_pulse_blocks = self.settings['decoupling_seq']['num_of_pulse_blocks']
        tau_total = self.settings['tau_time']

        if self.settings['decoupling_seq']['type'] == 'spin_echo':
            for amp in amp_list:
                tau = tau_total / (1 * self.settings['decoupling_seq']['num_of_pulse_blocks'])
                pulse_sequence = []

                # ECHO SEQUENCE:
                start_of_first_HE = laser_off_time

                # the first pi/2 pulse
                pulse_sequence = [Pulse(microwave_channel_1, start_of_first_HE, pi_half_time_1)]
                section_begin_time = start_of_first_HE + pi_half_time_1 - tau / 2  # for the first pulse, only wait tau/2
                for i in range(0, number_of_pulse_blocks):

                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time_1 / 2, pi_time_1))
                    section_begin_time += 1 * tau

                # the second pi/2 pulse and readout
                pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + tau / 2, pi_half_time_1))
                pulse_sequence.append(Pulse('laser', section_begin_time + tau / 2 + pi_half_time_1 + delay_mw_readout,
                                             nv_reset_time))
                pulse_sequence.append(Pulse('apd_readout',
                                             section_begin_time + tau / 2 + pi_half_time_1 + delay_mw_readout + delay_readout,
                                             meas_time))
                start_of_second_HE = section_begin_time + tau / 2 + pi_half_time_1 + delay_mw_readout + nv_reset_time + laser_off_time

                # the first pi/2 pulse
                pulse_sequence.append(Pulse(microwave_channel_1, start_of_second_HE, pi_half_time_1))
                section_begin_time = start_of_second_HE + pi_half_time_1 - tau / 2  # for the first pulse, only wait tau/2
                for i in range(0, number_of_pulse_blocks):
                    # pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time_1 / 2, pi_time_1))
                    section_begin_time += 1 * tau

                # the second 3*pi/2 pulse and readout
                pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + tau / 2, three_pi_half_time_1))
                pulse_sequence.append(Pulse('laser',
                                             section_begin_time + tau / 2 + three_pi_half_time_1 + delay_mw_readout,
                                             nv_reset_time))
                pulse_sequence.append(Pulse('apd_readout',
                                             section_begin_time + tau / 2 + three_pi_half_time_1 + delay_mw_readout + delay_readout,
                                             meas_time))

                # E-Sensing SEQUENCE
                if True:
                    start_of_DEER = section_begin_time + tau / 2 + three_pi_half_time_1 + delay_mw_readout + nv_reset_time + laser_off_time
                    # the first pi/2 pulse
                    pulse_sequence.append(Pulse(microwave_channel_1, start_of_DEER, pi_half_time_1))
                    pulse_sequence.append(Pulse('fg_trigger', start_of_DEER - trigger_latency, trigger_time))
                    section_begin_time = start_of_DEER + pi_half_time_1 - tau / 2  # for the first pulse, only wait tau/2
                    for i in range(0, number_of_pulse_blocks):
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time_1 / 2, pi_time_1))
                        section_begin_time += 1 * tau

                    # the second pi/2 pulse and readout
                    pulse_sequence.append(Pulse(last_mw_channel, section_begin_time + tau / 2, last_pi_half_time))
                    pulse_sequence.append(Pulse('laser', section_begin_time + tau / 2 + last_pi_half_time + delay_mw_readout,
                                                nv_reset_time))
                    pulse_sequence.append(Pulse('apd_readout',
                                                section_begin_time + tau / 2 + last_pi_half_time + delay_mw_readout + delay_readout,
                                                meas_time))

                    start_of_second_DEER = section_begin_time + tau / 2 + last_pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time

                    # the first pi/2 pulse
                    pulse_sequence.append(Pulse(microwave_channel_1, start_of_second_DEER, pi_half_time_1))
                    pulse_sequence.append(Pulse('fg_trigger', start_of_second_DEER - trigger_latency, trigger_time))
                    section_begin_time = start_of_second_DEER + pi_half_time_1 - tau / 2  # for the first pulse, only wait tau/2

                    for i in range(0, number_of_pulse_blocks):
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time_1 / 2, pi_time_1))
                        section_begin_time += 1 * tau

                    # the second pi/2 pulse and readout
                    pulse_sequence.append(Pulse(last_mw_channel, section_begin_time + tau / 2, last_three_pi_half_time))
                    pulse_sequence.append(Pulse('laser', section_begin_time + tau / 2 + last_three_pi_half_time + delay_mw_readout,
                                                nv_reset_time))
                    pulse_sequence.append(Pulse('apd_readout',
                                                section_begin_time + tau / 2 + last_three_pi_half_time + delay_mw_readout + delay_readout,
                                                meas_time))

                # Turn on APD switch all the time
                pulse_sequence.append(Pulse('apd_switch', 0,
                                            section_begin_time + tau / 2 + last_three_pi_half_time + delay_mw_readout + nv_reset_time))

                # print(pulse_sequence)
                pulse_sequences.append(pulse_sequence)

            print('number of sequences before validation ', len(pulse_sequences))
            return pulse_sequences, amp_list, meas_time

        elif self.settings['decoupling_seq']['type'] == 'XY4':

            for amp in amp_list:
                tau = tau_total / (4 * self.settings['decoupling_seq']['num_of_pulse_blocks'])
                pulse_sequence = []

                # ECHO SEQUENCE:
                start_of_first_HE = laser_off_time
                # the first pi/2 pulse
                pulse_sequence = [Pulse(microwave_channel_1, start_of_first_HE, pi_half_time_1)]
                section_begin_time = start_of_first_HE + pi_half_time_1 - tau / 2  # for the first pulse, only wait tau/2
                for i in range(0, number_of_pulse_blocks):
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time_2 / 2, pi_time_2))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 2 * tau - pi_time_1 / 2, pi_time_1))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time_1 / 2, pi_time_1))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time_2 / 2, pi_time_2))

                    section_begin_time += 4 * tau

                # the second 3*pi/2 pulse and readout
                pulse_sequence.append(
                    Pulse(microwave_channel_1, section_begin_time + tau / 2, three_pi_half_time_1))
                pulse_sequence.append(
                    Pulse('laser', section_begin_time + tau / 2 + three_pi_half_time_1 + delay_mw_readout,
                          nv_reset_time))
                pulse_sequence.append(
                    Pulse('apd_readout',
                          section_begin_time + tau / 2 + three_pi_half_time_1 + delay_mw_readout + delay_readout,
                          meas_time))
                start_of_second_HE = section_begin_time + tau / 2 + three_pi_half_time_1 + delay_mw_readout + nv_reset_time + laser_off_time

                # the first pi/2 pulse
                pulse_sequence.append(
                    Pulse(microwave_channel_1, start_of_second_HE, pi_half_time_1))
                section_begin_time = start_of_second_HE + pi_half_time_1 - tau / 2  # for the first pulse, only wait tau/2
                for i in range(0, number_of_pulse_blocks):
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time_2 / 2, pi_time_2))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 2 * tau - pi_time_1 / 2, pi_time_1))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time_1 / 2, pi_time_1))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time_2 / 2, pi_time_2))
                    section_begin_time += 4 * tau

                # the second pi/2 pulse and readout
                pulse_sequence.append(
                    Pulse(microwave_channel_1, section_begin_time + tau / 2, pi_half_time_1))
                pulse_sequence.append(
                    Pulse('laser',
                          section_begin_time + tau / 2 + pi_half_time_1 + delay_mw_readout,
                          nv_reset_time))
                pulse_sequence.append(
                    Pulse('apd_readout',
                          section_begin_time + tau / 2 + pi_half_time_1 + delay_mw_readout + delay_readout,
                          meas_time))

                # E-Sensing SEQUENCE
                if True:
                    start_of_DEER = section_begin_time + tau / 2 + pi_half_time_1 + delay_mw_readout + nv_reset_time + laser_off_time
                    # the first pi/2 pulse
                    pulse_sequence.append(Pulse(microwave_channel_1, start_of_DEER, pi_half_time_1))
                    pulse_sequence.append(Pulse('fg_trigger', start_of_DEER - trigger_latency, trigger_time))
                    section_begin_time = start_of_DEER + pi_half_time_1 - tau / 2  # for the first pulse, only wait tau/2
                    for i in range(0, number_of_pulse_blocks):
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time_2 / 2, pi_time_2))
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 2 * tau - pi_time_1 / 2, pi_time_1))
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time_1 / 2, pi_time_1))
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time_2 / 2, pi_time_2))

                        section_begin_time += 4 * tau

                    # the second 3*pi/2 pulse and readout
                    pulse_sequence.append(
                        Pulse(last_mw_channel, section_begin_time + tau / 2, last_three_pi_half_time))
                    pulse_sequence.append(
                        Pulse('laser', section_begin_time + tau / 2 + last_three_pi_half_time + delay_mw_readout,
                              nv_reset_time))
                    pulse_sequence.append(
                        Pulse('apd_readout',
                              section_begin_time + tau / 2 + last_three_pi_half_time + delay_mw_readout + delay_readout,
                              meas_time))

                    start_of_second_DEER = section_begin_time + tau / 2 + last_three_pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time

                    # the first pi/2 pulse
                    pulse_sequence.append(Pulse(microwave_channel_1, start_of_second_DEER, pi_half_time_1))
                    pulse_sequence.append(Pulse('fg_trigger', start_of_second_DEER - trigger_latency, trigger_time))
                    section_begin_time = start_of_second_DEER + pi_half_time_1 - tau / 2  # for the first pulse, only wait tau/2

                    for i in range(0, number_of_pulse_blocks):
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time_2 / 2, pi_time_2))
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 2 * tau - pi_time_1 / 2, pi_time_1))
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time_1 / 2, pi_time_1))
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time_2 / 2, pi_time_2))
                        section_begin_time += 4 * tau

                    # the second pi/2 pulse and readout
                    pulse_sequence.append(
                        Pulse(last_mw_channel, section_begin_time + tau / 2, last_pi_half_time))
                    pulse_sequence.append(
                        Pulse('laser',
                              section_begin_time + tau / 2 + last_pi_half_time + delay_mw_readout,
                              nv_reset_time))
                    pulse_sequence.append(
                        Pulse('apd_readout',
                              section_begin_time + tau / 2 + last_pi_half_time + delay_mw_readout + delay_readout,
                              meas_time))

                pulse_sequence.append(Pulse('apd_switch', 0,
                                            section_begin_time + tau / 2 + last_pi_half_time + delay_mw_readout + nv_reset_time))
                pulse_sequences.append(pulse_sequence)

            print('number of sequences before validation ', len(pulse_sequences))
            # return pulse_sequences, self.settings['num_averages'], tau_list, meas_time
            return pulse_sequences, amp_list, meas_time

        elif self.settings['decoupling_seq']['type'] == 'XY8':

            for amp in amp_list:
                tau = tau_total / (8 * self.settings['decoupling_seq']['num_of_pulse_blocks'])
                pulse_sequence = []

                # ECHO SEQUENCE:
                start_of_first_HE = laser_off_time

                # the first pi/2 pulse
                pulse_sequence.append(
                    Pulse(microwave_channel_1, start_of_first_HE, pi_half_time_1))

                section_begin_time = start_of_first_HE + pi_half_time_1 - tau / 2  # for the first pulse, only wait tau/2
                for i in range(0, number_of_pulse_blocks):

                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time_1 / 2, pi_time_1))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 2 * tau - pi_time_2 / 2, pi_time_2))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time_1 / 2, pi_time_1))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time_2 / 2, pi_time_2))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 5 * tau - pi_time_2 / 2, pi_time_2))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 6 * tau - pi_time_1 / 2, pi_time_1))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 7 * tau - pi_time_2 / 2, pi_time_2))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 8 * tau - pi_time_1 / 2, pi_time_1))
                    section_begin_time += 8 * tau

                # the second 3*pi/2 pulse and readout
                pulse_sequence.append(
                    Pulse(microwave_channel_1, section_begin_time + tau / 2, three_pi_half_time_1))
                pulse_sequence.append(
                    Pulse('laser', section_begin_time + tau / 2 + three_pi_half_time_1 + delay_mw_readout,
                          nv_reset_time))
                pulse_sequence.append(
                    Pulse('apd_readout',
                          section_begin_time + tau / 2 + three_pi_half_time_1 + delay_mw_readout + delay_readout,
                          meas_time))

                start_of_second_HE = section_begin_time + tau / 2 + three_pi_half_time_1 + delay_mw_readout + nv_reset_time + laser_off_time

                # the first pi/2 pulse
                pulse_sequence.append(
                    Pulse(microwave_channel_1, start_of_second_HE, pi_half_time_1))
                section_begin_time = start_of_second_HE + pi_half_time_1 - tau / 2  # for the first pulse, only wait tau/2
                for i in range(0, number_of_pulse_blocks):

                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time_1 / 2, pi_time_1))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 2 * tau - pi_time_2 / 2, pi_time_2))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time_1 / 2, pi_time_1))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time_2 / 2, pi_time_2))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 5 * tau - pi_time_2 / 2, pi_time_2))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 6 * tau - pi_time_1 / 2, pi_time_1))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 7 * tau - pi_time_2 / 2, pi_time_2))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 8 * tau - pi_time_1 / 2, pi_time_1))
                    section_begin_time += 8 * tau

                # the second pi/2 pulse and readout
                pulse_sequence.append(
                    Pulse(microwave_channel_1, section_begin_time + tau / 2, pi_half_time_1))
                pulse_sequence.append(
                    Pulse('laser',
                          section_begin_time + tau / 2 + pi_half_time_1 + delay_mw_readout,
                          nv_reset_time))
                pulse_sequence.append(
                    Pulse('apd_readout',
                          section_begin_time + tau / 2 + pi_half_time_1 + delay_mw_readout + delay_readout,
                          meas_time))

                # E-Sensing SEQUENCE
                if True:
                    start_of_DEER = section_begin_time + tau / 2 + pi_half_time_1 + delay_mw_readout + nv_reset_time + laser_off_time
                    # the first pi/2 pulse
                    pulse_sequence.append(Pulse(microwave_channel_1, start_of_DEER, pi_half_time_1))
                    pulse_sequence.append(Pulse('fg_trigger', start_of_DEER - trigger_latency, trigger_time))

                    section_begin_time = start_of_DEER + pi_half_time_1 - tau / 2  # for the first pulse, only wait tau/2
                    for i in range(0, number_of_pulse_blocks):
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time_1 / 2, pi_time_1))
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 2 * tau - pi_time_2 / 2, pi_time_2))
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time_1 / 2, pi_time_1))
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time_2 / 2, pi_time_2))
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 5 * tau - pi_time_2 / 2, pi_time_2))
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 6 * tau - pi_time_1 / 2, pi_time_1))
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 7 * tau - pi_time_2 / 2, pi_time_2))
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 8 * tau - pi_time_1 / 2, pi_time_1))
                        section_begin_time += 8 * tau

                    # the second 3*pi/2 pulse and readout
                    pulse_sequence.append(
                        Pulse(last_mw_channel, section_begin_time + tau / 2, last_three_pi_half_time))
                    pulse_sequence.append(
                        Pulse('laser', section_begin_time + tau / 2 + last_three_pi_half_time + delay_mw_readout,
                              nv_reset_time))
                    pulse_sequence.append(
                        Pulse('apd_readout',
                              section_begin_time + tau / 2 + last_three_pi_half_time + delay_mw_readout + delay_readout,
                              meas_time))

                    start_of_second_DEER = section_begin_time + tau / 2 + last_three_pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time

                    # the first pi/2 pulse
                    pulse_sequence.append(Pulse(microwave_channel_1, start_of_second_DEER, pi_half_time_1))
                    pulse_sequence.append(Pulse('fg_trigger', start_of_second_DEER - trigger_latency, trigger_time))
                    section_begin_time = start_of_second_DEER + pi_half_time_1 - tau / 2  # for the first pulse, only wait tau/2

                    for i in range(0, number_of_pulse_blocks):
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time_1 / 2, pi_time_1))
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 2 * tau - pi_time_2 / 2, pi_time_2))
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time_1 / 2, pi_time_1))
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time_2 / 2, pi_time_2))
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 5 * tau - pi_time_2 / 2, pi_time_2))
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 6 * tau - pi_time_1 / 2, pi_time_1))
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 7 * tau - pi_time_2 / 2, pi_time_2))
                        pulse_sequence.append(
                            Pulse(microwave_channel_1, section_begin_time + 8 * tau - pi_time_1 / 2, pi_time_1))
                        section_begin_time += 8 * tau

                    # the second pi/2 pulse and readout
                    pulse_sequence.append(
                        Pulse(last_mw_channel, section_begin_time + tau / 2, last_pi_half_time))
                    pulse_sequence.append(
                        Pulse('laser',
                              section_begin_time + tau / 2 + last_pi_half_time + delay_mw_readout,
                              nv_reset_time))
                    pulse_sequence.append(
                        Pulse('apd_readout',
                              section_begin_time + tau / 2 + last_pi_half_time + delay_mw_readout + delay_readout,
                              meas_time))

                pulse_sequence.append(Pulse('apd_switch', 0,
                                            section_begin_time + tau / 2 + last_pi_half_time + delay_mw_readout + nv_reset_time))
                pulse_sequences.append(pulse_sequence)

            print('number of sequences before validation ', len(pulse_sequences))
            # return pulse_sequences, self.settings['num_averages'], tau_list, meas_time
            return pulse_sequences, amp_list, meas_time

        elif self.settings['decoupling_seq']['type'] == 'CPMG':

            for amp in amp_list:
                tau = tau_total / (1 * self.settings['decoupling_seq']['num_of_pulse_blocks'])
                pulse_sequence = []

                # ECHO SEQUENCE:
                start_of_first_HE = laser_off_time

                # the first pi/2 pulse
                pulse_sequence = [Pulse(microwave_channel_1, start_of_first_HE, pi_half_time_1)]
                section_begin_time = start_of_first_HE + pi_half_time_1 - tau / 2  # for the first pulse, only wait tau/2

                for i in range(0, number_of_pulse_blocks):
                    pulse_sequence.append(Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time_2 / 2, pi_time_2))
                    section_begin_time += 1 * tau

                # the second 3*pi/2 pulse and readout
                pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + tau / 2, three_pi_half_time_1))
                pulse_sequence.append(Pulse('laser', section_begin_time + tau / 2 + three_pi_half_time_1 + delay_mw_readout,
                                             nv_reset_time))
                pulse_sequence.append(Pulse('apd_readout',
                                             section_begin_time + tau / 2 + three_pi_half_time_1 + delay_mw_readout + delay_readout,
                                             meas_time))

                start_of_second_HE = section_begin_time + tau / 2 + three_pi_half_time_1 + delay_mw_readout + nv_reset_time + laser_off_time

                # the first pi/2 pulse
                pulse_sequence.append(Pulse(microwave_channel_1, start_of_second_HE, pi_half_time_1))

                section_begin_time = start_of_second_HE + pi_half_time_1 - tau / 2  # for the first pulse, only wait tau/2

                for i in range(0, number_of_pulse_blocks):
                    pulse_sequence.append(Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time_2 / 2, pi_time_2))
                    section_begin_time += 1 * tau

                # the second pi/2 pulse and readout
                pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + tau / 2, pi_half_time_1))
                pulse_sequence.append(Pulse('laser',
                                             section_begin_time + tau / 2 + pi_half_time_1 + delay_mw_readout,
                                             nv_reset_time))
                pulse_sequence.append(Pulse('apd_readout',
                                             section_begin_time + tau / 2 + pi_half_time_1 + delay_mw_readout + delay_readout,
                                             meas_time))

                # E-Sensing SEQUENCE
                if True:
                    start_of_DEER = section_begin_time + tau / 2 + pi_half_time_1 + delay_mw_readout + nv_reset_time + laser_off_time
                    # the first pi/2 pulse
                    pulse_sequence.append(Pulse(microwave_channel_1, start_of_DEER, pi_half_time_1))
                    pulse_sequence.append(Pulse('fg_trigger', start_of_DEER - trigger_latency, trigger_time))
                    section_begin_time = start_of_DEER + pi_half_time_1 - tau / 2  # for the first pulse, only wait tau/2
                    for i in range(0, number_of_pulse_blocks):
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time_2 / 2, pi_time_2))
                        section_begin_time += 1 * tau

                    # the second 3*pi/2 pulse and readout
                    pulse_sequence.append(Pulse(last_mw_channel, section_begin_time + tau / 2, last_three_pi_half_time))
                    pulse_sequence.append(Pulse('laser', section_begin_time + tau / 2 + last_three_pi_half_time + delay_mw_readout,
                                                nv_reset_time))
                    pulse_sequence.append(Pulse('apd_readout',
                                                section_begin_time + tau / 2 + last_three_pi_half_time + delay_mw_readout + delay_readout,
                                                meas_time))

                    start_of_second_DEER = section_begin_time + tau / 2 + last_three_pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time

                    # the first pi/2 pulse
                    pulse_sequence.append(Pulse(microwave_channel_1, start_of_second_DEER, pi_half_time_1))
                    pulse_sequence.append(Pulse('fg_trigger', start_of_second_DEER - trigger_latency, trigger_time))
                    section_begin_time = start_of_second_DEER + pi_half_time_1 - tau / 2  # for the first pulse, only wait tau/2

                    for i in range(0, number_of_pulse_blocks):
                        pulse_sequence.append(
                            Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time_2 / 2, pi_time_2))
                        section_begin_time += 1 * tau


                    # the second pi/2 pulse and readout
                    pulse_sequence.append(Pulse(last_mw_channel, section_begin_time + tau / 2, last_pi_half_time))
                    pulse_sequence.append(
                        Pulse('laser', section_begin_time + tau / 2 + last_pi_half_time + delay_mw_readout,
                              nv_reset_time))
                    pulse_sequence.append(Pulse('apd_readout',
                                                section_begin_time + tau / 2 + last_pi_half_time + delay_mw_readout + delay_readout,
                                                meas_time))

                # Turn on APD switch all the time
                pulse_sequence.append(Pulse('apd_switch', 0,
                                            section_begin_time + tau / 2 + last_pi_half_time + delay_mw_readout + nv_reset_time))


                pulse_sequences.append(pulse_sequence)

            print('number of sequences before validation ', len(pulse_sequences))
            # return pulse_sequences, self.settings['num_averages'], tau_list, meas_time
            return pulse_sequences, amp_list, meas_time

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

        # if 'fits' in data.keys() is not None and data['fits'] is not None:
        #
        #     # counts = (-data['counts'][:,1] + data['counts'][:,0])/ (data['counts'][:,0] + data['counts'][:,1])
        #     tau = data['tau']
        #     fits = data['fits']
        #
        #     axislist[0].plot(tau, data['norm_echo'])
        #     tauinterp = np.linspace(np.min(tau), np.max(tau), 100)
        #     axislist[0].plot(tauinterp, exp_offset(tauinterp, fits[0], fits[1], fits[2]))
        #
        #
        #     if 'fits_esig' in data.keys() is not None and data['fits_esig'] is not None:
        #         fits_esig = data['fits_esig']
        #         axislist[0].plot(tau, data['norm_esig'])
        #         axislist[0].plot(tauinterp, exp_offset(tauinterp, fits_esig[0], fits_esig[1], fits_esig[2]))
        #         axislist[0].legend(labels=('Echo', 'Exp Fit','E-Sensing', 'Exp Fit'), fontsize=10)
        #         axislist[0].set_title(
        #             'Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
        #                 'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s), {:s} electrometry\n averaged over {:d} blocks \n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz, \n T2 (p = 1) = {:2.1f} ns, T2 (E-field) = {:2.1f} ns'.format(
        #                 self.settings['decoupling_seq']['type'],
        #                 self.settings['decoupling_seq']['num_of_pulse_blocks'],
        #                 self.settings['decoupling_seq']['eSensing_type'], avg_block_number,
        #                 self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency'] * 1e-9,
        #                 fits[1], fits_esig[1]))
        #     else:
        #         axislist[0].legend(labels=('Echo', 'Exp Fit'), fontsize=10)
        #         axislist[0].set_title(
        #             'Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
        #                 'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s)\n averaged over {:d} blocks \n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz, \n T2 (p = 1) = {:2.1f} ns'.format(
        #                 self.settings['decoupling_seq']['type'],
        #                 self.settings['decoupling_seq']['num_of_pulse_blocks'], avg_block_number,
        #                 self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency'] * 1e-9,
        #                 fits[1]))
        #     axislist[0].set_ylabel('fluorescence contrast')
        #     axislist[0].set_xlabel('tau [ns]')

        if 'exp_finished' in data.keys() is not None and data['exp_finished'] == 1:
            if 'norm_echo' in data.keys() is not None and data['norm_echo'] is not None:
                tau = data['tau']
                axislist[0].plot(tau, data['norm_echo'])
                if 'norm_esig' in data.keys() is not None and data['norm_esig'] is not None:
                    axislist[0].plot(tau, data['norm_esig'])
                    axislist[0].legend(labels=('PDD', 'E-Sensing'), fontsize=10)
                    axislist[0].set_title(
                        '(final plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
                            'microwave_channel'] + '\nPDD: {:s} {:d} block(s), {:s} electrometry\n averaged over {:d} blocks \n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz \n tau = {:.2f} us'.format(
                            self.settings['decoupling_seq']['type'],
                            self.settings['decoupling_seq']['num_of_pulse_blocks'],
                            self.settings['decoupling_seq']['eSensing_type'], avg_block_number,
                            self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency'] * 1e-9, self.settings['tau_time'] * 1e-3))
                else:
                    axislist[0].legend(labels=('PDD'), fontsize=10)
                    axislist[0].set_title(
                        '(final plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
                            'microwave_channel'] + '\nPDD: {:s} {:d} block(s)\n averaged over {:d} blocks \n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz \n tau = {:.2f} us'.format(
                            self.settings['decoupling_seq']['type'],
                            self.settings['decoupling_seq']['num_of_pulse_blocks'],
                            avg_block_number,
                            self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency'] * 1e-9, self.settings['tau_time'] * 1e-3))
                axislist[0].set_ylabel('fluorescence contrast')
                axislist[0].set_xlabel('tau [ns]')
            else:
                axislist[0].plot(data['tau'], data['counts'][:, 0])
                axislist[0].plot(data['tau'], data['counts'][:, 1])
                axislist[0].plot(data['tau'], data['counts'][:, 2])
                axislist[0].plot(data['tau'], data['counts'][:, 3])
                axislist[0].legend(labels=('PDD up {:.0f}kcps'.format(np.mean(data['counts'][:, 0])),
                                           'PDD down {:.0f}kcps'.format(np.mean(data['counts'][:, 1])),
                                           'E-Sensing up {:.0f}kcps'.format(np.mean(data['counts'][:, 0])),
                                           'E-Sensing down {:.0f}kcps'.format(np.mean(data['counts'][:, 1]))),
                                   fontsize=10)
                axislist[0].set_title(
                    '(final plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
                        'microwave_channel'] + '\nPDD: {:s} {:d} block(s), {:s} electrometry\n averaged over {:d} blocks\n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz \n tau = {:.2f} us'.format(
                        self.settings['decoupling_seq']['type'],
                        self.settings['decoupling_seq']['num_of_pulse_blocks'],
                        self.settings['decoupling_seq']['eSensing_type'],
                        avg_block_number, self.settings['mw_pulses']['mw_power'],
                        self.settings['mw_pulses']['mw_frequency'] * 1e-9, self.settings['tau_time'] * 1e-3))


                axislist[0].set_ylabel('fluorescence [kcps]')
                axislist[0].set_xlabel('tau [ns]')

        else:
            super(eSensing_swpTau, self)._plot(axislist)
            axislist[0].set_xlabel('tau [ns]')
            axislist[0].legend(labels=('PDD up {:.0f}kcps'.format(np.mean(data['counts'][:, 0])),
                                       'PDD down {:.0f}kcps'.format(np.mean(data['counts'][:, 1])),
                                       'E-Sensing up {:.0f}kcps'.format(np.mean(data['counts'][:, 0])),
                                       'E-Sensing down {:.0f}kcps'.format(np.mean(data['counts'][:, 1]))),
                               fontsize=10)
            axislist[0].set_title(
                '(initial plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
                    'microwave_channel'] + '\nPDD: {:s} {:d} block(s), {:s} electrometry\n averaged over {:d} blocks\n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz \n tau = {:.2f} us'.format(
                    self.settings['decoupling_seq']['type'], self.settings['decoupling_seq']['num_of_pulse_blocks'],
                    self.settings['decoupling_seq']['eSensing_type'],
                    avg_block_number, self.settings['mw_pulses']['mw_power'],
                    self.settings['mw_pulses']['mw_frequency'] * 1e-9, self.settings['tau_time'] * 1e-3))

            # if self.settings['decoupling_seq']['do_eSensing']:
            #     axislist[0].legend(labels=('Echo up {:.0f}kcps'.format(np.mean(data['counts'][:, 0])),
            #                                'Echo down {:.0f}kcps'.format(np.mean(data['counts'][:, 1])),
            #                                'E-Sensing up {:.0f}kcps'.format(np.mean(data['counts'][:, 0])),
            #                                'E-Sensing down {:.0f}kcps'.format(np.mean(data['counts'][:, 1]))),
            #                        fontsize=10)
            #     axislist[0].set_title(
            #         '(initial plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
            #             'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s), {:s} electrometry\n averaged over {:d} blocks\n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz'.format(
            #             self.settings['decoupling_seq']['type'], self.settings['decoupling_seq']['num_of_pulse_blocks'],
            #             self.settings['decoupling_seq']['eSensing_type'],
            #             avg_block_number, self.settings['mw_pulses']['mw_power'],
            #             self.settings['mw_pulses']['mw_frequency'] * 1e-9))
            #
            # else:
            #     axislist[0].legend(labels=('Echo up {:.0f}kcps'.format(np.mean(data['counts'][:, 0])),
            #                                'Echo down {:.0f}kcps'.format(np.mean(data['counts'][:, 1]))), fontsize=10)
            #     axislist[0].set_title(
            #         '(initial plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
            #             'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s) \n averaged over {:d} blocks\n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz'.format(
            #             self.settings['decoupling_seq']['type'], self.settings['decoupling_seq']['num_of_pulse_blocks'],
            #             avg_block_number, self.settings['mw_pulses']['mw_power'],
            #             self.settings['mw_pulses']['mw_frequency'] * 1e-9))

    def _update_plot(self, axislist):

        if len(axislist[0].lines) == 0:
            self._plot(axislist)
            return

        super(eSensing_swpTau, self)._update_plot(axislist)
        axislist[0].legend(labels=('Echo up {:.0f}kcps'.format(np.mean(self.data['counts'][:, 0])),
                                   'Echo down {:.0f}kcps'.format(np.mean(self.data['counts'][:, 1])),
                                   'E-Sensing up {:.0f}kcps'.format(np.mean(self.data['counts'][:, 0])),
                                   'E-Sensing down {:.0f}kcps'.format(np.mean(self.data['counts'][:, 1]))),
                           fontsize=10)
        axislist[0].set_title(
            '(updating plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
                'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s), {:s} electrometry\n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz \n tau = {:.2f} us'.format(
                self.settings['decoupling_seq']['type'], self.settings['decoupling_seq']['num_of_pulse_blocks'],
                self.settings['decoupling_seq']['eSensing_type'],
                self.settings['mw_pulses']['mw_power'],
                self.settings['mw_pulses']['mw_frequency'] * 1e-9, self.settings['tau_time'] * 1e-3))

        # if self.settings['decoupling_seq']['do_eSensing']:
        #     axislist[0].legend(labels=('Echo up {:.0f}kcps'.format(np.mean(self.data['counts'][:, 0])),
        #                                'Echo down {:.0f}kcps'.format(np.mean(self.data['counts'][:, 1])),
        #                                'E-Sensing up {:.0f}kcps'.format(np.mean(self.data['counts'][:, 0])),
        #                                'E-Sensing down {:.0f}kcps'.format(np.mean(self.data['counts'][:, 1]))),
        #                        fontsize=10)
        #     axislist[0].set_title(
        #         '(updating plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
        #             'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s), {:s} electrometry\n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz'.format(
        #             self.settings['decoupling_seq']['type'], self.settings['decoupling_seq']['num_of_pulse_blocks'],
        #             self.settings['decoupling_seq']['eSensing_type'],
        #             self.settings['mw_pulses']['mw_power'],
        #             self.settings['mw_pulses']['mw_frequency'] * 1e-9))
        #
        # else:
        #     axislist[0].legend(labels=('Echo up {:.0f}kcps'.format(np.mean(self.data['counts'][:, 0])),
        #                                'Echo down {:.0f}kcps'.format(np.mean(self.data['counts'][:, 1]))), fontsize=10)
        #
        #     axislist[0].set_title('(updating plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
        #         'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s) \n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz'.format(
        #         self.settings['decoupling_seq']['type'], self.settings['decoupling_seq']['num_of_pulse_blocks'],
        #         self.settings['mw_pulses']['mw_power'],
        #         self.settings['mw_pulses']['mw_frequency'] * 1e-9))

class PDD_RnSIQ(PulsedExperimentBaseScript):
    """
    This script runs a PDD ( Periodic Dynamical Decoupling) sequence for different number of pi pulses.
    To symmetrize the sequence between the 0 and +/-1 state we reinitialize every time.
    The script uses Rohde&Schwarz signal generator SMB100A. The phase is controlled by an IQ mixer.
    For a single pi-pulse this is a Hahn-echo sequence. For zero pulses this is a Ramsey sequence.
    The sequence is pi/2 - tau/4 - (tau/4 - pi  - tau/4)^n - tau/4 - pi/2
    Tau/2 is the time between the center of the pulses!

    --> Ziwei Qiu 5/11/2019

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
            Parameter('time_step', 100., [2.5, 5., 10., 20., 50., 100., 200., 300., 400., 500., 600., 800., 1000.,2000.,5000., 8000., 10000.,50000., 100000., 500000.],
                      'time step increment of time between pi pulses (in ns)')
        ]),
        Parameter('decoupling_seq', [
            Parameter('type', 'spin_echo', ['spin_echo', 'CPMG', 'XY4', 'XY8'], 'type of dynamical decoupling sequences'),
            Parameter('num_of_pulse_blocks', 1, int, 'number of pulse blocks.'),
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 235, float, 'measurement time after rabi sequence (in ns)'),
            Parameter('nv_reset_time', 2000, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 600, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 550, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('min_pulse_dur', 15, [10, 15, 50], 'Minimum allowed pulse duration (ns)'),
        Parameter('num_averages', 500000, int, 'number of averages (>100000)'),
    ]

    _INSTRUMENTS = {'daq': NI6353, 'PB': LISE607RTPulseBlaster, 'mw_gen': R8SMicrowaveGenerator}

    def _function(self):
        #COMMENT_ME

        self.data['fits'] = None

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

        self.avg_block_number = super(PDD_RnSIQ, self)._function()

        # this is to flag when the experiment is done or stopped and this allows proper saving even fit fails
        self.data['exp_finished'] = 1
        self.data['avg_block_number'] = self.avg_block_number

        # turn off laser, apd_switch and MW (just in case)
        self.instruments['mw_gen']['instance'].update({'enable_output': False})
        self.instruments['PB']['instance'].update({'laser': {'status': False}})
        self.instruments['PB']['instance'].update({'apd_switch': {'status': False}})

        # start fitting
        if 'counts' in self.data.keys() and 'tau' in self.data.keys() and not 0 in (self.data['counts'][:,1] + self.data['counts'][:, 0]):
            self.data['norm_echo'] = 2. * (- self.data['counts'][:, 1] + self.data['counts'][:, 0]) / (
                    self.data['counts'][:, 1] + self.data['counts'][:, 0])
            counts = self.data['norm_echo']

            # counts = (- self.data['counts'][:, 1] + self.data['counts'][:,0]) / (self.data['counts'][:,1] + self.data['counts'][:, 0])
            tau = self.data['tau']
            fit_success = False

            try:
                fits = fit_exp_decay(tau, counts, offset = True, verbose = True)
                self.data['fits'] = fits
                fit_success = True
            except:
                self.data['fits'] = None
                self.log('T2 fit failed')
                fit_success = False

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
            microwave_channel_1 = 'microwave_i_II'
            microwave_channel_2 = 'microwave_q_II'
        else:
            microwave_channel_1 = 'microwave_q_II'
            microwave_channel_2 = 'microwave_i_II'
        pi_time = self.settings['mw_pulses']['pi_pulse_time']
        pi_half_time = self.settings['mw_pulses']['pi_half_pulse_time']
        three_pi_half_time = self.settings['mw_pulses']['3pi_half_pulse_time']

        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']

        number_of_pulse_blocks = self.settings['decoupling_seq']['num_of_pulse_blocks']

        # for tau in tau_list:
        #     pulse_sequence = \
        #     [
        #         Pulse(microwave_channel, laser_off_time, pi_half_time),
        #         Pulse(microwave_channel, laser_off_time + pi_half_time/2. + tau - pi_time/2., pi_time),
        #         Pulse(microwave_channel, laser_off_time + pi_half_time/2. + tau + tau - pi_half_time/2., pi_half_time)
        #     ]
        #
        #     end_of_first_HE = laser_off_time + pi_half_time/2. + tau + tau - pi_half_time/2. + pi_half_time
        #
        #     pulse_sequence += [
        #          Pulse('laser', end_of_first_HE + delay_mw_readout, nv_reset_time),
        #          Pulse('apd_readout', end_of_first_HE + delay_mw_readout + delay_readout, meas_time),
        #          ]
        #
        #     start_of_second_HE = end_of_first_HE + delay_mw_readout + nv_reset_time + laser_off_time
        #
        #     pulse_sequence += \
        #     [
        #         Pulse(microwave_channel, start_of_second_HE, pi_half_time),
        #         Pulse(microwave_channel, start_of_second_HE + pi_half_time/2. + tau - pi_time/2., pi_time),
        #         Pulse(microwave_channel, start_of_second_HE + pi_half_time/2. + tau + tau - pi_half_time/2., three_pi_half_time)
        #     ]
        #
        #     end_of_second_HE = start_of_second_HE + pi_half_time/2. + tau + tau - pi_half_time/2. + pi_half_time
        #
        #     pulse_sequence += [
        #         Pulse('laser', end_of_second_HE + delay_mw_readout, nv_reset_time),
        #         Pulse('apd_readout', end_of_second_HE + delay_mw_readout + delay_readout, meas_time)
        #     ]
        #
        #     pulse_sequence.append(Pulse('apd_switch', 0, end_of_second_HE + delay_mw_readout + nv_reset_time))
        #     # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
        #     # if tau == 0 or tau>=15:
        #     pulse_sequences.append(pulse_sequence)

        if self.settings['decoupling_seq']['type'] == 'spin_echo':
            for tau_total in tau_list:
                tau = tau_total / (1 * self.settings['decoupling_seq']['num_of_pulse_blocks'])
                pulse_sequence = []

                # ECHO SEQUENCE:
                start_of_first_HE = laser_off_time
                # the first pi/2 pulse

                pulse_sequence = [Pulse(microwave_channel_1, start_of_first_HE, pi_half_time)]
                section_begin_time = start_of_first_HE + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                for i in range(0, number_of_pulse_blocks):
                    pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                    section_begin_time += 1 * tau

                # the second pi/2 pulse and readout
                pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + tau / 2, pi_half_time))
                pulse_sequence.append(Pulse('laser', section_begin_time + tau / 2 + pi_half_time + delay_mw_readout,
                                             nv_reset_time))
                pulse_sequence.append(Pulse('apd_readout',
                                             section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + delay_readout,
                                             meas_time))
                start_of_second_HE = section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time

                # the first pi/2 pulse
                pulse_sequence.append(Pulse(microwave_channel_1, start_of_second_HE, pi_half_time))
                section_begin_time = start_of_second_HE + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                for i in range(0, number_of_pulse_blocks):
                    pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                    section_begin_time += 1 * tau

                # the second 3*pi/2 pulse and readout
                pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + tau / 2, three_pi_half_time))
                pulse_sequence.append(Pulse('laser',
                                             section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout,
                                             nv_reset_time))
                pulse_sequence.append(Pulse('apd_readout',
                                             section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + delay_readout,
                                             meas_time))

                # Turn on APD switch all the time
                pulse_sequence.append(Pulse('apd_switch', 0,
                                            section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + nv_reset_time))

                # # DEER SEQUENCE
                # if self.settings['decoupling_seq']['do_deer']:
                #     start_of_DEER = section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time
                #     # the first pi/2 pulse
                #     pulse_sequence.extend([Pulse(microwave_channel_1, start_of_DEER, pi_half_time)])
                #     section_begin_time = start_of_DEER + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                #     for i in range(0, number_of_pulse_blocks):
                #         pulse_sequence.extend(
                #             [Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 1 * tau - RF_pi_time / 2, RF_pi_time)
                #              ])
                #         section_begin_time += 1 * tau
                #
                #     # the second pi/2 pulse and readout
                #     pulse_sequence.extend([Pulse(microwave_channel_1, section_begin_time + tau / 2, pi_half_time),
                #                            Pulse('laser',
                #                                  section_begin_time + tau / 2 + pi_half_time + delay_mw_readout,
                #                                  nv_reset_time),
                #                            Pulse('apd_readout',
                #                                  section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + delay_readout,
                #                                  meas_time)
                #                            ])
                #
                #     start_of_second_DEER = section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time
                #     pulse_sequence.extend([Pulse(microwave_channel_1, start_of_second_DEER, pi_half_time)])
                #     section_begin_time = start_of_second_DEER + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                #
                #     for i in range(0, number_of_pulse_blocks):
                #         pulse_sequence.extend(
                #             [Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 1 * tau - RF_pi_time / 2, RF_pi_time)
                #              ])
                #         section_begin_time += 1 * tau
                #
                #     # the second 3*pi/2 pulse and readout
                #     pulse_sequence.extend([Pulse(microwave_channel_1, section_begin_time + tau / 2, three_pi_half_time),
                #                            Pulse('laser',
                #                                  section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout,
                #                                  nv_reset_time),
                #                            Pulse('apd_readout',
                #                                  section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + delay_readout,
                #                                  meas_time)
                #                            ])

                # print(pulse_sequence)
                pulse_sequences.append(pulse_sequence)

            print('number of sequences before validation ', len(pulse_sequences))
            return pulse_sequences, tau_list, meas_time
            # return pulse_sequences, self.settings['num_averages'], tau_list, meas_time

        elif self.settings['decoupling_seq']['type'] == 'XY4':

            for tau_total in tau_list:
                tau = tau_total / (4 * self.settings['decoupling_seq']['num_of_pulse_blocks'])
                pulse_sequence = []

                # ECHO SEQUENCE:
                start_of_first_HE = laser_off_time
                # the first pi/2 pulse
                pulse_sequence = [Pulse(microwave_channel_1, start_of_first_HE, pi_half_time)]
                section_begin_time = start_of_first_HE + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                for i in range(0, number_of_pulse_blocks):
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 2 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time / 2, pi_time))

                    section_begin_time += 4 * tau

                # the second pi/2 pulse and readout
                pulse_sequence.append(
                    Pulse(microwave_channel_1, section_begin_time + tau / 2, three_pi_half_time))
                pulse_sequence.append(
                    Pulse('laser', section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout,
                          nv_reset_time))
                pulse_sequence.append(
                    Pulse('apd_readout',
                          section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + delay_readout,
                          meas_time))
                start_of_second_HE = section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time

                # the first pi/2 pulse
                pulse_sequence.append(
                    Pulse(microwave_channel_1, start_of_second_HE, pi_half_time))
                section_begin_time = start_of_second_HE + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                for i in range(0, number_of_pulse_blocks):
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 2 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time / 2, pi_time))
                    section_begin_time += 4 * tau

                # the second 3*pi/2 pulse and readout
                pulse_sequence.append(
                    Pulse(microwave_channel_1, section_begin_time + tau / 2, pi_half_time))
                pulse_sequence.append(
                    Pulse('laser',
                          section_begin_time + tau / 2 + pi_half_time + delay_mw_readout,
                          nv_reset_time))
                pulse_sequence.append(
                    Pulse('apd_readout',
                          section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + delay_readout,
                          meas_time))

                # # DEER SEQUENCE
                # if self.settings['decoupling_seq']['do_deer']:
                #     start_of_DEER = section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time
                #     # the first pi/2 pulse
                #     pulse_sequence.extend([Pulse(microwave_channel_1, start_of_DEER, pi_half_time)])
                #     section_begin_time = start_of_DEER + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                #     for i in range(0, number_of_pulse_blocks):
                #         pulse_sequence.extend(
                #             [Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 1 * tau - RF_pi_time / 2, RF_pi_time),
                #              Pulse(microwave_channel_1, section_begin_time + 2 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 2 * tau - RF_pi_time / 2, RF_pi_time),
                #              Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 3 * tau - RF_pi_time / 2, RF_pi_time),
                #              Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 4 * tau - RF_pi_time / 2, RF_pi_time)
                #              ])
                #         section_begin_time += 4 * tau
                #
                #     # the second pi/2 pulse and readout
                #     pulse_sequence.extend([Pulse(microwave_channel_1, section_begin_time + tau / 2, pi_half_time),
                #                            Pulse('laser',
                #                                  section_begin_time + tau / 2 + pi_half_time + delay_mw_readout,
                #                                  nv_reset_time),
                #                            Pulse('apd_readout',
                #                                  section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + delay_readout,
                #                                  meas_time)
                #                            ])
                #
                #     start_of_second_DEER = section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time
                #     pulse_sequence.extend([Pulse(microwave_channel_1, start_of_second_DEER, pi_half_time)])
                #     section_begin_time = start_of_second_DEER + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                #
                #     for i in range(0, number_of_pulse_blocks):
                #         pulse_sequence.extend(
                #             [Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 1 * tau - RF_pi_time / 2, RF_pi_time),
                #              Pulse(microwave_channel_1, section_begin_time + 2 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 2 * tau - RF_pi_time / 2, RF_pi_time),
                #              Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 3 * tau - RF_pi_time / 2, RF_pi_time),
                #              Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 4 * tau - RF_pi_time / 2, RF_pi_time)
                #              ])
                #         section_begin_time += 4 * tau
                #
                #     # the second 3*pi/2 pulse and readout
                #     pulse_sequence.extend([Pulse(microwave_channel_1, section_begin_time + tau / 2, three_pi_half_time),
                #                            Pulse('laser',
                #                                  section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout,
                #                                  nv_reset_time),
                #                            Pulse('apd_readout',
                #                                  section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + delay_readout,
                #                                  meas_time)
                #                            ])
                # Turn on APD switch all the time
                pulse_sequence.append(Pulse('apd_switch', 0,
                                            section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + nv_reset_time))
                pulse_sequences.append(pulse_sequence)

            print('number of sequences before validation ', len(pulse_sequences))
            # return pulse_sequences, self.settings['num_averages'], tau_list, meas_time
            return pulse_sequences, tau_list, meas_time

        elif self.settings['decoupling_seq']['type'] == 'XY8':

            for tau_total in tau_list:
                tau = tau_total / (8 * self.settings['decoupling_seq']['num_of_pulse_blocks'])
                pulse_sequence = []

                # ECHO SEQUENCE:
                start_of_first_HE = laser_off_time

                # the first pi/2 pulse
                pulse_sequence.append(
                    Pulse(microwave_channel_1, start_of_first_HE, pi_half_time))

                section_begin_time = start_of_first_HE + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                for i in range(0, number_of_pulse_blocks):
                    # pulse_sequence.extend(
                    #     # [Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time),
                    #     #  Pulse(microwave_channel_1, section_begin_time + 2 * tau - pi_time / 2, pi_time),
                    #     #  Pulse(microwave_channel_2, section_begin_time + 3 * tau - pi_time / 2, pi_time),
                    #     #  Pulse(microwave_channel_1, section_begin_time + 4 * tau - pi_time / 2, pi_time),
                    #     #  Pulse(microwave_channel_1, section_begin_time + 5 * tau - pi_time / 2, pi_time),
                    #     #  Pulse(microwave_channel_2, section_begin_time + 6 * tau - pi_time / 2, pi_time),
                    #     #  Pulse(microwave_channel_1, section_begin_time + 7 * tau - pi_time / 2, pi_time),
                    #     #  Pulse(microwave_channel_2, section_begin_time + 8 * tau - pi_time / 2, pi_time)
                    #     [Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time),
                    #      Pulse(microwave_channel_2, section_begin_time + 2 * tau - pi_time / 2, pi_time),
                    #      Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time / 2, pi_time),
                    #      Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time / 2, pi_time),
                    #      Pulse(microwave_channel_2, section_begin_time + 5 * tau - pi_time / 2, pi_time),
                    #      Pulse(microwave_channel_1, section_begin_time + 6 * tau - pi_time / 2, pi_time),
                    #      Pulse(microwave_channel_2, section_begin_time + 7 * tau - pi_time / 2, pi_time),
                    #      Pulse(microwave_channel_1, section_begin_time + 8 * tau - pi_time / 2, pi_time)
                    #      ])
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 2 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 5 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 6 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 7 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 8 * tau - pi_time / 2, pi_time))
                    section_begin_time += 8 * tau

                # the second pi/2 pulse and readout
                pulse_sequence.append(
                    Pulse(microwave_channel_1, section_begin_time + tau / 2, three_pi_half_time))
                pulse_sequence.append(
                    Pulse('laser', section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout,
                          nv_reset_time))
                pulse_sequence.append(
                    Pulse('apd_readout',
                          section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + delay_readout,
                          meas_time))

                start_of_second_HE = section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time

                # the first pi/2 pulse
                pulse_sequence.append(
                    Pulse(microwave_channel_1, start_of_second_HE, pi_half_time))
                section_begin_time = start_of_second_HE + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                for i in range(0, number_of_pulse_blocks):
                    # pulse_sequence.extend(
                    #     # [Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time),
                    #     #  Pulse(microwave_channel_1, section_begin_time + 2 * tau - pi_time / 2, pi_time),
                    #     #  Pulse(microwave_channel_2, section_begin_time + 3 * tau - pi_time / 2, pi_time),
                    #     #  Pulse(microwave_channel_1, section_begin_time + 4 * tau - pi_time / 2, pi_time),
                    #     #  Pulse(microwave_channel_1, section_begin_time + 5 * tau - pi_time / 2, pi_time),
                    #     #  Pulse(microwave_channel_2, section_begin_time + 6 * tau - pi_time / 2, pi_time),
                    #     #  Pulse(microwave_channel_1, section_begin_time + 7 * tau - pi_time / 2, pi_time),
                    #     #  Pulse(microwave_channel_2, section_begin_time + 8 * tau - pi_time / 2, pi_time)
                    #     [Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time),
                    #      Pulse(microwave_channel_2, section_begin_time + 2 * tau - pi_time / 2, pi_time),
                    #      Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time / 2, pi_time),
                    #      Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time / 2, pi_time),
                    #      Pulse(microwave_channel_2, section_begin_time + 5 * tau - pi_time / 2, pi_time),
                    #      Pulse(microwave_channel_1, section_begin_time + 6 * tau - pi_time / 2, pi_time),
                    #      Pulse(microwave_channel_2, section_begin_time + 7 * tau - pi_time / 2, pi_time),
                    #      Pulse(microwave_channel_1, section_begin_time + 8 * tau - pi_time / 2, pi_time)
                    #      ])
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 2 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 5 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 6 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_2, section_begin_time + 7 * tau - pi_time / 2, pi_time))
                    pulse_sequence.append(
                        Pulse(microwave_channel_1, section_begin_time + 8 * tau - pi_time / 2, pi_time))
                    section_begin_time += 8 * tau

                # the second 3*pi/2 pulse and readout
                pulse_sequence.append(
                    Pulse(microwave_channel_1, section_begin_time + tau / 2, pi_half_time))
                pulse_sequence.append(
                    Pulse('laser',
                          section_begin_time + tau / 2 + pi_half_time + delay_mw_readout,
                          nv_reset_time))
                pulse_sequence.append(
                    Pulse('apd_readout',
                          section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + delay_readout,
                          meas_time))
                # # DEER SEQUENCE
                # if self.settings['decoupling_seq']['do_deer']:
                #     start_of_DEER = section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time
                #     # the first pi/2 pulse
                #     pulse_sequence.extend([Pulse(microwave_channel_1, start_of_DEER, pi_half_time)])
                #     section_begin_time = start_of_DEER + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                #     for i in range(0, number_of_pulse_blocks):
                #         pulse_sequence.extend(
                #             [Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 1 * tau - RF_pi_time / 2, RF_pi_time),
                #              Pulse(microwave_channel_1, section_begin_time + 2 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 2 * tau - RF_pi_time / 2, RF_pi_time),
                #              Pulse(microwave_channel_2, section_begin_time + 3 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 3 * tau - RF_pi_time / 2, RF_pi_time),
                #              Pulse(microwave_channel_1, section_begin_time + 4 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 4 * tau - RF_pi_time / 2, RF_pi_time),
                #              Pulse(microwave_channel_1, section_begin_time + 5 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 5 * tau - RF_pi_time / 2, RF_pi_time),
                #              Pulse(microwave_channel_2, section_begin_time + 6 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 6 * tau - RF_pi_time / 2, RF_pi_time),
                #              Pulse(microwave_channel_1, section_begin_time + 7 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 7 * tau - RF_pi_time / 2, RF_pi_time),
                #              Pulse(microwave_channel_2, section_begin_time + 8 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 8 * tau - RF_pi_time / 2, RF_pi_time)
                #              ])
                #         section_begin_time += 8 * tau
                #
                #     # the second pi/2 pulse and readout
                #     pulse_sequence.extend([Pulse(microwave_channel_1, section_begin_time + tau / 2, pi_half_time),
                #                            Pulse('laser',
                #                                  section_begin_time + tau / 2 + pi_half_time + delay_mw_readout,
                #                                  nv_reset_time),
                #                            Pulse('apd_readout',
                #                                  section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + delay_readout,
                #                                  meas_time)
                #                            ])
                #
                #     start_of_second_DEER = section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time
                #     pulse_sequence.extend([Pulse(microwave_channel_1, start_of_second_DEER, pi_half_time)])
                #     section_begin_time = start_of_second_DEER + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                #
                #     for i in range(0, number_of_pulse_blocks):
                #         pulse_sequence.extend(
                #             [Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 1 * tau - RF_pi_time / 2, RF_pi_time),
                #              Pulse(microwave_channel_1, section_begin_time + 2 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 2 * tau - RF_pi_time / 2, RF_pi_time),
                #              Pulse(microwave_channel_2, section_begin_time + 3 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 3 * tau - RF_pi_time / 2, RF_pi_time),
                #              Pulse(microwave_channel_1, section_begin_time + 4 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 4 * tau - RF_pi_time / 2, RF_pi_time),
                #              Pulse(microwave_channel_1, section_begin_time + 5 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 5 * tau - RF_pi_time / 2, RF_pi_time),
                #              Pulse(microwave_channel_2, section_begin_time + 6 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 6 * tau - RF_pi_time / 2, RF_pi_time),
                #              Pulse(microwave_channel_1, section_begin_time + 7 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 7 * tau - RF_pi_time / 2, RF_pi_time),
                #              Pulse(microwave_channel_2, section_begin_time + 8 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 8 * tau - RF_pi_time / 2, RF_pi_time)
                #              ])
                #         section_begin_time += 8 * tau
                #
                #     # the second 3*pi/2 pulse and readout
                #     pulse_sequence.extend([Pulse(microwave_channel_1, section_begin_time + tau / 2, three_pi_half_time),
                #                            Pulse('laser',
                #                                  section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout,
                #                                  nv_reset_time),
                #                            Pulse('apd_readout',
                #                                  section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + delay_readout,
                #                                  meas_time)
                #                            ])
                # Turn on APD switch all the time
                pulse_sequence.append(Pulse('apd_switch', 0,
                                            section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + nv_reset_time))
                pulse_sequences.append(pulse_sequence)

            print('number of sequences before validation ', len(pulse_sequences))
            # return pulse_sequences, self.settings['num_averages'], tau_list, meas_time
            return pulse_sequences, tau_list, meas_time

        elif self.settings['decoupling_seq']['type'] == 'CPMG':

            for tau_total in tau_list:
                tau = tau_total / (1 * self.settings['decoupling_seq']['num_of_pulse_blocks'])
                pulse_sequence = []

                # ECHO SEQUENCE:
                start_of_first_HE = laser_off_time

                # the first pi/2 pulse
                pulse_sequence = [Pulse(microwave_channel_1, start_of_first_HE, pi_half_time)]
                section_begin_time = start_of_first_HE + pi_half_time - tau / 2  # for the first pulse, only wait tau/2

                for i in range(0, number_of_pulse_blocks):
                    pulse_sequence.append(Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                    section_begin_time += 1 * tau

                # the second pi/2 pulse and readout
                pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + tau / 2, three_pi_half_time))
                pulse_sequence.append(Pulse('laser', section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout,
                                             nv_reset_time))
                pulse_sequence.append(Pulse('apd_readout',
                                             section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + delay_readout,
                                             meas_time))

                start_of_second_HE = section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time

                # the first pi/2 pulse
                pulse_sequence.append(Pulse(microwave_channel_1, start_of_second_HE, pi_half_time))

                section_begin_time = start_of_second_HE + pi_half_time - tau / 2  # for the first pulse, only wait tau/2

                for i in range(0, number_of_pulse_blocks):
                    pulse_sequence.append(Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                    section_begin_time += 1 * tau

                # the second 3*pi/2 pulse and readout
                pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + tau / 2, pi_half_time))
                pulse_sequence.append(Pulse('laser',
                                             section_begin_time + tau / 2 + pi_half_time + delay_mw_readout,
                                             nv_reset_time))
                pulse_sequence.append(Pulse('apd_readout',
                                             section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + delay_readout,
                                             meas_time))

                # Turn on APD switch all the time
                pulse_sequence.append(Pulse('apd_switch', 0,
                                            section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + nv_reset_time))


    #             # # DEER SEQUENCE
    #             # if self.settings['decoupling_seq']['do_deer']:
    #             #     start_of_DEER = section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time
    #             #     # the first pi/2 pulse
    #             #     pulse_sequence.extend([Pulse(microwave_channel_1, start_of_DEER, pi_half_time)])
    #             #     section_begin_time = start_of_DEER + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
    #             #     for i in range(0, number_of_pulse_blocks):
    #             #         pulse_sequence.extend(
    #             #             [Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time),
    #             #              Pulse('RF_switch', section_begin_time + 1 * tau - RF_pi_time / 2, RF_pi_time)
    #             #              ])
    #             #         section_begin_time += 1 * tau
    #             #
    #             #     # the second pi/2 pulse and readout
    #             #     pulse_sequence.extend([Pulse(microwave_channel_1, section_begin_time + tau / 2, pi_half_time),
    #             #                            Pulse('laser',
    #             #                                  section_begin_time + tau / 2 + pi_half_time + delay_mw_readout,
    #             #                                  nv_reset_time),
    #             #                            Pulse('apd_readout',
    #             #                                  section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + delay_readout,
    #             #                                  meas_time)
    #             #                            ])
    #             #
    #             #     start_of_second_DEER = section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time
    #             #     pulse_sequence.extend([Pulse(microwave_channel_1, start_of_second_DEER, pi_half_time)])
    #             #     section_begin_time = start_of_second_DEER + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
    #             #
    #             #     for i in range(0, number_of_pulse_blocks):
    #             #         pulse_sequence.extend(
    #             #             [Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time),
    #             #              Pulse('RF_switch', section_begin_time + 1 * tau - RF_pi_time / 2, RF_pi_time)
    #             #              ])
    #             #         section_begin_time += 1 * tau
    #             #
    #             #     # the second 3*pi/2 pulse and readout
    #             #     pulse_sequence.extend([Pulse(microwave_channel_1, section_begin_time + tau / 2, three_pi_half_time),
    #             #                            Pulse('laser',
    #             #                                  section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout,
    #             #                                  nv_reset_time),
    #             #                            Pulse('apd_readout',
    #             #                                  section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + delay_readout,
    #             #                                  meas_time)
    #             #                            ])
    #
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

        if 'avg_block_number' in data.keys() is not None and data['avg_block_number'] is not None:
            avg_block_number = data['avg_block_number']
        else:
            avg_block_number = 0

        if 'fits' in data.keys() is not None and data['fits'] is not None:

            # counts = (-data['counts'][:,1] + data['counts'][:,0])/ (data['counts'][:,0] + data['counts'][:,1])
            tau = data['tau']
            fits = data['fits']
            axislist[0].plot(tau, data['norm_echo'])
            tauinterp = np.linspace(np.min(tau), np.max(tau), 100)
            axislist[0].plot(tauinterp, exp_offset(tauinterp, fits[0], fits[1], fits[2]))
            axislist[0].set_title(
                'R&S SMB100A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
                    'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s)\n averaged over {:d} blocks \n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz, \n T2 (simple expo, p = 1) = {:2.1f} ns'.format(
                    self.settings['decoupling_seq']['type'],
                    self.settings['decoupling_seq']['num_of_pulse_blocks'], avg_block_number,
                    self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency'] * 1e-9,
                    fits[1]))
            axislist[0].legend(labels=('Echo', 'Exp Fit'), fontsize=10)
            axislist[0].set_ylabel('fluorescence contrast')
            axislist[0].set_xlabel('tau [ns]')

        elif 'exp_finished' in data.keys() is not None and data['exp_finished'] == 1:
            if 'norm_echo' in data.keys() is not None and data['norm_echo'] is not None:
                tau = data['tau']
                axislist[0].plot(tau, data['norm_echo'],'d')
                axislist[0].set_title(
                    '(final plot) R&S SMB100A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
                        'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s)\n averaged over {:d} blocks \n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz'.format(
                        self.settings['decoupling_seq']['type'],
                        self.settings['decoupling_seq']['num_of_pulse_blocks'],
                        avg_block_number,
                        self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency'] * 1e-9))
                axislist[0].legend(labels=('Echo'), fontsize=10)
                axislist[0].set_ylabel('fluorescence contrast')
                axislist[0].set_xlabel('tau [ns]')
            else:
                axislist[0].plot(data['tau'], data['counts'][:, 0])
                axislist[0].plot(data['tau'], data['counts'][:, 1])
                axislist[0].set_ylabel('fluorescence [kcps]')
                axislist[0].set_xlabel('tau [ns]')
                axislist[0].legend(labels=('Echo up {:.0f}kcps'.format(np.mean(data['counts'][:, 0])),
                                           'Echo down {:.0f}kcps'.format(np.mean(data['counts'][:, 1]))), fontsize=10)
                axislist[0].set_title(
                    '(final plot) R&S SMB100A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
                        'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s)\n averaged over {:d} blocks\n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz'.format(
                        self.settings['decoupling_seq']['type'], self.settings['decoupling_seq']['num_of_pulse_blocks'],
                        avg_block_number, self.settings['mw_pulses']['mw_power'],
                        self.settings['mw_pulses']['mw_frequency'] * 1e-9))

        else:
            super(PDD_RnSIQ, self)._plot(axislist)

            axislist[0].legend(labels=('Echo up {:.0f}kcps'.format(np.mean(data['counts'][:, 0])),
                                       'Echo down {:.0f}kcps'.format(np.mean(data['counts'][:, 1]))), fontsize=10)
            axislist[0].set_title(
                '(initial plot) R&S SMB100A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
                    'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s) \n averaged over {:d} blocks\n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz'.format(
                    self.settings['decoupling_seq']['type'], self.settings['decoupling_seq']['num_of_pulse_blocks'],
                    avg_block_number, self.settings['mw_pulses']['mw_power'],
                    self.settings['mw_pulses']['mw_frequency'] * 1e-9))

    def _update_plot(self, axislist):

        if len(axislist[0].lines) == 0:
            self._plot(axislist)
            return
        # if 'avg_block_number' in self.data.keys() is not None and self.data['avg_block_number'] is not None:
        #     avg_block_number =self.data['avg_block_number']
        # else:
        #     avg_block_number = 0
        super(PDD_RnSIQ, self)._update_plot(axislist)
        axislist[0].legend(labels=('Echo up {:.0f}kcps'.format(np.mean(self.data['counts'][:, 0])),
                                   'Echo down {:.0f}kcps'.format(np.mean(self.data['counts'][:, 1]))), fontsize=10)
        # axislist[0].set_title('(updating plot) Keysight N9310A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
        #     'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s) \n averaged over {:d} blocks \n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz'.format(
        #     self.settings['decoupling_seq']['type'], self.settings['decoupling_seq']['num_of_pulse_blocks'],
        #     avg_block_number, self.settings['mw_pulses']['mw_power'],
        #     self.settings['mw_pulses']['mw_frequency'] * 1e-9))
        axislist[0].set_title('(updating plot) R&S SMB100A, $\pi/2$ pulse carved by ' + self.settings['mw_pulses'][
            'microwave_channel'] + '\nT2 measurement: {:s} {:d} block(s) \n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz'.format(
            self.settings['decoupling_seq']['type'], self.settings['decoupling_seq']['num_of_pulse_blocks'],
            self.settings['mw_pulses']['mw_power'],
            self.settings['mw_pulses']['mw_frequency'] * 1e-9))

class PDD_RnS(PulsedExperimentBaseScript):
    """
    ATTENTION: The microwave generator used in this script doesn't have IQ option...
    This script runs a PDD ( Periodic Dynamical Decoupling) sequence for different number of pi pulses.
    To symmetrize the sequence between the 0 and +/-1 state we reinitialize every time.
    For a single pi-pulse this is a Hahn-echo sequence. For zero pulses this is a Ramsey sequence.
    The sequence is pi/2 - tau/4 - (tau/4 - pi  - tau/4)^n - tau/4 - pi/2
    Tau/2 is the time between the center of the pulses!

    --> Ziwei Qiu 3/26/2019

    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_pulses', [
            Parameter('mw_power', -10.0, float, 'microwave power in dB'),
            Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
            # Parameter('microwave_channel', 'i', ['i', 'q'], 'Channel to use for mw pulses'),
            Parameter('pi_pulse_time', 50.0, float, 'time duration of a pi pulse (in ns)'),
            Parameter('pi_half_pulse_time', 25.0, float, 'time duration of a pi/2 pulse (in ns)'),
            Parameter('3pi_half_pulse_time', 75.0, float, 'time duration of a 3pi/2 pulse (in ns)')
        ]),
        Parameter('tau_times', [
            Parameter('min_time', 200, float, 'minimum time between pi pulses'),
            Parameter('max_time', 10000, float, 'maximum time between pi pulses'),
            Parameter('time_step', 100., [2.5, 5., 10., 20., 50., 100., 200., 500., 1000., 10000., 100000., 500000.],
                      'time step increment of time between pi pulses (in ns)')
        ]),
        Parameter('decoupling_seq', [
            # Parameter('type', 'spin_echo', ['spin_echo', 'XY4', 'XY8', 'CPMG'], 'type of dynamical decoupling sequences'),
            Parameter('type', 'spin_echo', ['spin_echo'],
                      'type of dynamical decoupling sequences (only spin-echo is available since no IQ option is available)'),
            Parameter('num_of_pulse_blocks', 1, int, 'number of pulse blocks.'),
        ]),
        Parameter('read_out', [
            Parameter('meas_time', 235, float, 'measurement time after rabi sequence (in ns)'),
            Parameter('nv_reset_time', 2000, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 600, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 550, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('min_pulse_dur', 15, [10, 15, 50], 'Minimum allowed pulse duration (ns)'),
        Parameter('num_averages', 100000, int, 'number of averages'),
    ]

    # _INSTRUMENTS = {'daq': NI6259, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator}
    _INSTRUMENTS = {'daq': NI6353, 'PB': LISE607RTPulseBlaster, 'mw_gen': R8SMicrowaveGenerator}

    def _function(self):
        #COMMENT_ME

        ## Turn off green light (the pulse blaster will pulse it on when needed)
        self.instruments['PB']['instance'].update({'laser': {'status': False}})

        self.data['fits'] = None
        # self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
        # self.instruments['mw_gen']['instance'].update({'enable_modulation': True}) # ER 20181018
        # self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power']})
        # self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulses']['mw_frequency']})
        self.instruments['mw_gen']['instance'].update({'enable_output': True})
        self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_pulses']['mw_power']})
        self.instruments['mw_gen']['instance'].update({'frequency': self.settings['mw_pulses']['mw_frequency']})
        self.instruments['mw_gen']['instance'].update({'freq_mode': 'CW'})
        self.instruments['mw_gen']['instance'].update({'power_mode': 'CW'})

        # self.avg_block_number = super(PDD_RnS, self)._function(self.data)
        self.avg_block_number = super(PDD_RnS, self)._function()

        # this is to flag when the experiment is done or stopped and this allows proper saving even fit fails
        self.data['exp_finished'] = 1

        # turn off laser, apd_switch and MW (just in case)
        self.instruments['PB']['instance'].update({'laser': {'status': False}})
        self.instruments['PB']['instance'].update({'apd_switch': {'status': False}})
        self.instruments['mw_gen']['instance'].update({'enable_output': False})


        # start fitting
        if 'counts' in self.data.keys() and 'tau' in self.data.keys() and not 0 in (self.data['counts'][:,1] + self.data['counts'][:, 0]):
            self.data['norm_echo'] = 2. * (- self.data['counts'][:, 1] + self.data['counts'][:, 0]) / (
                    self.data['counts'][:, 1] + self.data['counts'][:, 0])
            counts = self.data['norm_echo']

            # counts = (- self.data['counts'][:, 1] + self.data['counts'][:,0]) / (self.data['counts'][:,1] + self.data['counts'][:, 0])
            tau = self.data['tau']
            fit_success = False

            try:
                fits = fit_exp_decay(tau, counts, offset = True, verbose = True)
                self.data['fits'] = fits
                fit_success = True
            except:
                self.data['fits'] = None
                self.log('T2 fit failed')
                fit_success = False

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
        # tau_list = range(int(max(15, self.settings['tau_times']['time_step'])), int(self.settings['tau_times']['max_time'] + 15),
        #                  self.settings['tau_times']['time_step'])
        # JG 16-08-25 changed (15ns min spacing is taken care of later):
        # tau_list = np.arange(self.settings['tau_times']['min_time'], self.settings['tau_times']['max_time'],self.settings['tau_times']['time_step'])
        # tau_list = np.ndarray.tolist(tau_list) # 20180731 ER convert to list

        max_range = int(np.floor((self.settings['tau_times']['max_time'] - self.settings['tau_times']['min_time']) /
                                 self.settings['tau_times']['time_step']))
        tau_list = np.array(
            [self.settings['tau_times']['min_time'] + i * self.settings['tau_times']['time_step'] for i in
             range(max_range)])

        # ignore the sequence if the mw-pulse is shorter than 15ns (0 is ok because there is no mw pulse!)
        tau_list = [x for x in tau_list if x == 0 or x >= MIN_DURATION]

        nv_reset_time = self.settings['read_out']['nv_reset_time']
        delay_readout = self.settings['read_out']['delay_readout']
        # microwave_channel = 'microwave_' + self.settings['mw_pulses']['microwave_channel']
        microwave_channel_1 = 'microwave_switch_II'
        microwave_channel_2 = 'microwave_switch' # Should use I/Q!!!
        pi_time = self.settings['mw_pulses']['pi_pulse_time']
        pi_half_time = self.settings['mw_pulses']['pi_half_pulse_time']
        three_pi_half_time = self.settings['mw_pulses']['3pi_half_pulse_time']

        laser_off_time = self.settings['read_out']['laser_off_time']
        meas_time = self.settings['read_out']['meas_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']

        number_of_pulse_blocks = self.settings['decoupling_seq']['num_of_pulse_blocks']

        # for tau in tau_list:
        #     pulse_sequence = \
        #     [
        #         Pulse(microwave_channel, laser_off_time, pi_half_time),
        #         Pulse(microwave_channel, laser_off_time + pi_half_time/2. + tau - pi_time/2., pi_time),
        #         Pulse(microwave_channel, laser_off_time + pi_half_time/2. + tau + tau - pi_half_time/2., pi_half_time)
        #     ]
        #
        #     end_of_first_HE = laser_off_time + pi_half_time/2. + tau + tau - pi_half_time/2. + pi_half_time
        #
        #     pulse_sequence += [
        #          Pulse('laser', end_of_first_HE + delay_mw_readout, nv_reset_time),
        #          Pulse('apd_readout', end_of_first_HE + delay_mw_readout + delay_readout, meas_time),
        #          ]
        #
        #     start_of_second_HE = end_of_first_HE + delay_mw_readout + nv_reset_time + laser_off_time
        #
        #     pulse_sequence += \
        #     [
        #         Pulse(microwave_channel, start_of_second_HE, pi_half_time),
        #         Pulse(microwave_channel, start_of_second_HE + pi_half_time/2. + tau - pi_time/2., pi_time),
        #         Pulse(microwave_channel, start_of_second_HE + pi_half_time/2. + tau + tau - pi_half_time/2., three_pi_half_time)
        #     ]
        #
        #     end_of_second_HE = start_of_second_HE + pi_half_time/2. + tau + tau - pi_half_time/2. + pi_half_time
        #
        #     pulse_sequence += [
        #         Pulse('laser', end_of_second_HE + delay_mw_readout, nv_reset_time),
        #         Pulse('apd_readout', end_of_second_HE + delay_mw_readout + delay_readout, meas_time)
        #     ]
        #
        #     pulse_sequence.append(Pulse('apd_switch', 0, end_of_second_HE + delay_mw_readout + nv_reset_time))
        #     # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
        #     # if tau == 0 or tau>=15:
        #     pulse_sequences.append(pulse_sequence)

        if self.settings['decoupling_seq']['type'] == 'spin_echo':
            for tau_total in tau_list:
                tau = tau_total / (1 * self.settings['decoupling_seq']['num_of_pulse_blocks'])
                pulse_sequence = []
                # ECHO SEQUENCE:

                start_of_first_HE = laser_off_time
                # the first pi/2 pulse

                pulse_sequence = [Pulse(microwave_channel_1, start_of_first_HE, pi_half_time)]


                section_begin_time = start_of_first_HE + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                for i in range(0, number_of_pulse_blocks):

                    pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                    section_begin_time += 1 * tau

                # the second pi/2 pulse and readout
                pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + tau / 2, pi_half_time))
                pulse_sequence.append(Pulse('laser', section_begin_time + tau / 2 + pi_half_time + delay_mw_readout,
                                             nv_reset_time))
                pulse_sequence.append(Pulse('apd_readout',
                                             section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + delay_readout,
                                             meas_time))
                start_of_second_HE = section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time

                # the first pi/2 pulse
                pulse_sequence.append(Pulse(microwave_channel_1, start_of_second_HE, pi_half_time))
                section_begin_time = start_of_second_HE + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                for i in range(0, number_of_pulse_blocks):
                    pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time))
                    section_begin_time += 1 * tau

                # the second 3*pi/2 pulse and readout
                pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + tau / 2, three_pi_half_time))
                pulse_sequence.append(Pulse('laser',
                                             section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout,
                                             nv_reset_time))
                pulse_sequence.append(Pulse('apd_readout',
                                             section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + delay_readout,
                                             meas_time))

                # Turn on APD switch all the time
                pulse_sequence.append(Pulse('apd_switch', 0,
                                            section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + nv_reset_time))

                # # DEER SEQUENCE
                # if self.settings['decoupling_seq']['do_deer']:
                #     start_of_DEER = section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time
                #     # the first pi/2 pulse
                #     pulse_sequence.extend([Pulse(microwave_channel_1, start_of_DEER, pi_half_time)])
                #     section_begin_time = start_of_DEER + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                #     for i in range(0, number_of_pulse_blocks):
                #         pulse_sequence.extend(
                #             [Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 1 * tau - RF_pi_time / 2, RF_pi_time)
                #              ])
                #         section_begin_time += 1 * tau
                #
                #     # the second pi/2 pulse and readout
                #     pulse_sequence.extend([Pulse(microwave_channel_1, section_begin_time + tau / 2, pi_half_time),
                #                            Pulse('laser',
                #                                  section_begin_time + tau / 2 + pi_half_time + delay_mw_readout,
                #                                  nv_reset_time),
                #                            Pulse('apd_readout',
                #                                  section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + delay_readout,
                #                                  meas_time)
                #                            ])
                #
                #     start_of_second_DEER = section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time
                #     pulse_sequence.extend([Pulse(microwave_channel_1, start_of_second_DEER, pi_half_time)])
                #     section_begin_time = start_of_second_DEER + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
                #
                #     for i in range(0, number_of_pulse_blocks):
                #         pulse_sequence.extend(
                #             [Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time),
                #              Pulse('RF_switch', section_begin_time + 1 * tau - RF_pi_time / 2, RF_pi_time)
                #              ])
                #         section_begin_time += 1 * tau
                #
                #     # the second 3*pi/2 pulse and readout
                #     pulse_sequence.extend([Pulse(microwave_channel_1, section_begin_time + tau / 2, three_pi_half_time),
                #                            Pulse('laser',
                #                                  section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout,
                #                                  nv_reset_time),
                #                            Pulse('apd_readout',
                #                                  section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + delay_readout,
                #                                  meas_time)
                #                            ])

                # print(pulse_sequence)
                pulse_sequences.append(pulse_sequence)

            print('number of sequences before validation ', len(pulse_sequences))
            return pulse_sequences, tau_list, meas_time
            # return pulse_sequences, self.settings['num_averages'], tau_list, meas_time

        # elif self.settings['decoupling_seq']['type'] == 'XY4':
        #     for tau_total in tau_list:
        #         tau = tau_total / (4 * self.settings['decoupling_seq']['num_of_pulse_blocks'])
        #         pulse_sequence = []
        #         # ECHO SEQUENCE:
        #
        #         start_of_first_HE = laser_off_time
        #         # the first pi/2 pulse
        #         pulse_sequence.extend([Pulse(microwave_channel_1, start_of_first_HE, pi_half_time)])
        #         section_begin_time = start_of_first_HE + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
        #         for i in range(0, number_of_pulse_blocks):
        #             pulse_sequence.extend(
        #                 [Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time),
        #                  Pulse(microwave_channel_1, section_begin_time + 2 * tau - pi_time / 2, pi_time),
        #                  Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time / 2, pi_time),
        #                  Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time / 2, pi_time)
        #                  ])
        #             section_begin_time += 4 * tau
        #
        #         # the second pi/2 pulse and readout
        #         pulse_sequence.extend([Pulse(microwave_channel_1, section_begin_time + tau / 2, pi_half_time),
        #                                Pulse('laser', section_begin_time + tau / 2 + pi_half_time + delay_mw_readout,
        #                                      nv_reset_time),
        #                                Pulse('apd_readout',
        #                                      section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + delay_readout,
        #                                      meas_time)
        #                                ])
        #
        #         start_of_second_HE = section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time
        #         # the first pi/2 pulse
        #         pulse_sequence.extend([Pulse(microwave_channel_1, start_of_second_HE, pi_half_time)])
        #         section_begin_time = start_of_second_HE + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
        #         for i in range(0, number_of_pulse_blocks):
        #             pulse_sequence.extend(
        #                 [Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time),
        #                  Pulse(microwave_channel_1, section_begin_time + 2 * tau - pi_time / 2, pi_time),
        #                  Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time / 2, pi_time),
        #                  Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time / 2, pi_time)
        #                  ])
        #             section_begin_time += 4 * tau
        #
        #         # the second 3*pi/2 pulse and readout
        #         pulse_sequence.extend([Pulse(microwave_channel_1, section_begin_time + tau / 2, three_pi_half_time),
        #                                Pulse('laser',
        #                                      section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout,
        #                                      nv_reset_time),
        #                                Pulse('apd_readout',
        #                                      section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + delay_readout,
        #                                      meas_time)
        #                                ])
        #
        #         # DEER SEQUENCE
        #         if self.settings['decoupling_seq']['do_deer']:
        #             start_of_DEER = section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time
        #             # the first pi/2 pulse
        #             pulse_sequence.extend([Pulse(microwave_channel_1, start_of_DEER, pi_half_time)])
        #             section_begin_time = start_of_DEER + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
        #             for i in range(0, number_of_pulse_blocks):
        #                 pulse_sequence.extend(
        #                     [Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time),
        #                      Pulse('RF_switch', section_begin_time + 1 * tau - RF_pi_time / 2, RF_pi_time),
        #                      Pulse(microwave_channel_1, section_begin_time + 2 * tau - pi_time / 2, pi_time),
        #                      Pulse('RF_switch', section_begin_time + 2 * tau - RF_pi_time / 2, RF_pi_time),
        #                      Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time / 2, pi_time),
        #                      Pulse('RF_switch', section_begin_time + 3 * tau - RF_pi_time / 2, RF_pi_time),
        #                      Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time / 2, pi_time),
        #                      Pulse('RF_switch', section_begin_time + 4 * tau - RF_pi_time / 2, RF_pi_time)
        #                      ])
        #                 section_begin_time += 4 * tau
        #
        #             # the second pi/2 pulse and readout
        #             pulse_sequence.extend([Pulse(microwave_channel_1, section_begin_time + tau / 2, pi_half_time),
        #                                    Pulse('laser',
        #                                          section_begin_time + tau / 2 + pi_half_time + delay_mw_readout,
        #                                          nv_reset_time),
        #                                    Pulse('apd_readout',
        #                                          section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + delay_readout,
        #                                          meas_time)
        #                                    ])
        #
        #             start_of_second_DEER = section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time
        #             pulse_sequence.extend([Pulse(microwave_channel_1, start_of_second_DEER, pi_half_time)])
        #             section_begin_time = start_of_second_DEER + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
        #
        #             for i in range(0, number_of_pulse_blocks):
        #                 pulse_sequence.extend(
        #                     [Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time),
        #                      Pulse('RF_switch', section_begin_time + 1 * tau - RF_pi_time / 2, RF_pi_time),
        #                      Pulse(microwave_channel_1, section_begin_time + 2 * tau - pi_time / 2, pi_time),
        #                      Pulse('RF_switch', section_begin_time + 2 * tau - RF_pi_time / 2, RF_pi_time),
        #                      Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time / 2, pi_time),
        #                      Pulse('RF_switch', section_begin_time + 3 * tau - RF_pi_time / 2, RF_pi_time),
        #                      Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time / 2, pi_time),
        #                      Pulse('RF_switch', section_begin_time + 4 * tau - RF_pi_time / 2, RF_pi_time)
        #                      ])
        #                 section_begin_time += 4 * tau
        #
        #             # the second 3*pi/2 pulse and readout
        #             pulse_sequence.extend([Pulse(microwave_channel_1, section_begin_time + tau / 2, three_pi_half_time),
        #                                    Pulse('laser',
        #                                          section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout,
        #                                          nv_reset_time),
        #                                    Pulse('apd_readout',
        #                                          section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + delay_readout,
        #                                          meas_time)
        #                                    ])
        #
        #         pulse_sequences.append(pulse_sequence)
        #
        #     print('number of sequences before validation ', len(pulse_sequences))
        #     return pulse_sequences, self.settings['num_averages'], tau_list, meas_time
        #
        # elif self.settings['decoupling_seq']['type'] == 'XY8':
        #     for tau_total in tau_list:
        #         tau = tau_total / (8 * self.settings['decoupling_seq']['num_of_pulse_blocks'])
        #         pulse_sequence = []
        #         # ECHO SEQUENCE:
        #
        #         start_of_first_HE = laser_off_time
        #         # the first pi/2 pulse
        #         pulse_sequence.extend([Pulse(microwave_channel_1, start_of_first_HE, pi_half_time)])
        #         section_begin_time = start_of_first_HE + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
        #         for i in range(0, number_of_pulse_blocks):
        #             pulse_sequence.extend(
        #                 # [Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time),
        #                 #  Pulse(microwave_channel_1, section_begin_time + 2 * tau - pi_time / 2, pi_time),
        #                 #  Pulse(microwave_channel_2, section_begin_time + 3 * tau - pi_time / 2, pi_time),
        #                 #  Pulse(microwave_channel_1, section_begin_time + 4 * tau - pi_time / 2, pi_time),
        #                 #  Pulse(microwave_channel_1, section_begin_time + 5 * tau - pi_time / 2, pi_time),
        #                 #  Pulse(microwave_channel_2, section_begin_time + 6 * tau - pi_time / 2, pi_time),
        #                 #  Pulse(microwave_channel_1, section_begin_time + 7 * tau - pi_time / 2, pi_time),
        #                 #  Pulse(microwave_channel_2, section_begin_time + 8 * tau - pi_time / 2, pi_time)
        #                 [Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time),
        #                  Pulse(microwave_channel_2, section_begin_time + 2 * tau - pi_time / 2, pi_time),
        #                  Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time / 2, pi_time),
        #                  Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time / 2, pi_time),
        #                  Pulse(microwave_channel_2, section_begin_time + 5 * tau - pi_time / 2, pi_time),
        #                  Pulse(microwave_channel_1, section_begin_time + 6 * tau - pi_time / 2, pi_time),
        #                  Pulse(microwave_channel_2, section_begin_time + 7 * tau - pi_time / 2, pi_time),
        #                  Pulse(microwave_channel_1, section_begin_time + 8 * tau - pi_time / 2, pi_time)
        #                  ])
        #             section_begin_time += 8 * tau
        #
        #         # the second pi/2 pulse and readout
        #         pulse_sequence.extend([Pulse(microwave_channel_1, section_begin_time + tau / 2, pi_half_time),
        #                                Pulse('laser', section_begin_time + tau / 2 + pi_half_time + delay_mw_readout,
        #                                      nv_reset_time),
        #                                Pulse('apd_readout',
        #                                      section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + delay_readout,
        #                                      meas_time)
        #                                ])
        #
        #         start_of_second_HE = section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time
        #         # the first pi/2 pulse
        #         pulse_sequence.extend([Pulse(microwave_channel_1, start_of_second_HE, pi_half_time)])
        #         section_begin_time = start_of_second_HE + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
        #         for i in range(0, number_of_pulse_blocks):
        #             pulse_sequence.extend(
        #                 # [Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time),
        #                 #  Pulse(microwave_channel_1, section_begin_time + 2 * tau - pi_time / 2, pi_time),
        #                 #  Pulse(microwave_channel_2, section_begin_time + 3 * tau - pi_time / 2, pi_time),
        #                 #  Pulse(microwave_channel_1, section_begin_time + 4 * tau - pi_time / 2, pi_time),
        #                 #  Pulse(microwave_channel_1, section_begin_time + 5 * tau - pi_time / 2, pi_time),
        #                 #  Pulse(microwave_channel_2, section_begin_time + 6 * tau - pi_time / 2, pi_time),
        #                 #  Pulse(microwave_channel_1, section_begin_time + 7 * tau - pi_time / 2, pi_time),
        #                 #  Pulse(microwave_channel_2, section_begin_time + 8 * tau - pi_time / 2, pi_time)
        #                 [Pulse(microwave_channel_1, section_begin_time + 1 * tau - pi_time / 2, pi_time),
        #                  Pulse(microwave_channel_2, section_begin_time + 2 * tau - pi_time / 2, pi_time),
        #                  Pulse(microwave_channel_1, section_begin_time + 3 * tau - pi_time / 2, pi_time),
        #                  Pulse(microwave_channel_2, section_begin_time + 4 * tau - pi_time / 2, pi_time),
        #                  Pulse(microwave_channel_2, section_begin_time + 5 * tau - pi_time / 2, pi_time),
        #                  Pulse(microwave_channel_1, section_begin_time + 6 * tau - pi_time / 2, pi_time),
        #                  Pulse(microwave_channel_2, section_begin_time + 7 * tau - pi_time / 2, pi_time),
        #                  Pulse(microwave_channel_1, section_begin_time + 8 * tau - pi_time / 2, pi_time)
        #                  ])
        #             section_begin_time += 8 * tau
        #
        #         # the second 3*pi/2 pulse and readout
        #         pulse_sequence.extend([Pulse(microwave_channel_1, section_begin_time + tau / 2, three_pi_half_time),
        #                                Pulse('laser',
        #                                      section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout,
        #                                      nv_reset_time),
        #                                Pulse('apd_readout',
        #                                      section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + delay_readout,
        #                                      meas_time)
        #                                ])
        #
        #         # DEER SEQUENCE
        #         if self.settings['decoupling_seq']['do_deer']:
        #             start_of_DEER = section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time
        #             # the first pi/2 pulse
        #             pulse_sequence.extend([Pulse(microwave_channel_1, start_of_DEER, pi_half_time)])
        #             section_begin_time = start_of_DEER + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
        #             for i in range(0, number_of_pulse_blocks):
        #                 pulse_sequence.extend(
        #                     [Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time),
        #                      Pulse('RF_switch', section_begin_time + 1 * tau - RF_pi_time / 2, RF_pi_time),
        #                      Pulse(microwave_channel_1, section_begin_time + 2 * tau - pi_time / 2, pi_time),
        #                      Pulse('RF_switch', section_begin_time + 2 * tau - RF_pi_time / 2, RF_pi_time),
        #                      Pulse(microwave_channel_2, section_begin_time + 3 * tau - pi_time / 2, pi_time),
        #                      Pulse('RF_switch', section_begin_time + 3 * tau - RF_pi_time / 2, RF_pi_time),
        #                      Pulse(microwave_channel_1, section_begin_time + 4 * tau - pi_time / 2, pi_time),
        #                      Pulse('RF_switch', section_begin_time + 4 * tau - RF_pi_time / 2, RF_pi_time),
        #                      Pulse(microwave_channel_1, section_begin_time + 5 * tau - pi_time / 2, pi_time),
        #                      Pulse('RF_switch', section_begin_time + 5 * tau - RF_pi_time / 2, RF_pi_time),
        #                      Pulse(microwave_channel_2, section_begin_time + 6 * tau - pi_time / 2, pi_time),
        #                      Pulse('RF_switch', section_begin_time + 6 * tau - RF_pi_time / 2, RF_pi_time),
        #                      Pulse(microwave_channel_1, section_begin_time + 7 * tau - pi_time / 2, pi_time),
        #                      Pulse('RF_switch', section_begin_time + 7 * tau - RF_pi_time / 2, RF_pi_time),
        #                      Pulse(microwave_channel_2, section_begin_time + 8 * tau - pi_time / 2, pi_time),
        #                      Pulse('RF_switch', section_begin_time + 8 * tau - RF_pi_time / 2, RF_pi_time)
        #                      ])
        #                 section_begin_time += 8 * tau
        #
        #             # the second pi/2 pulse and readout
        #             pulse_sequence.extend([Pulse(microwave_channel_1, section_begin_time + tau / 2, pi_half_time),
        #                                    Pulse('laser',
        #                                          section_begin_time + tau / 2 + pi_half_time + delay_mw_readout,
        #                                          nv_reset_time),
        #                                    Pulse('apd_readout',
        #                                          section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + delay_readout,
        #                                          meas_time)
        #                                    ])
        #
        #             start_of_second_DEER = section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time
        #             pulse_sequence.extend([Pulse(microwave_channel_1, start_of_second_DEER, pi_half_time)])
        #             section_begin_time = start_of_second_DEER + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
        #
        #             for i in range(0, number_of_pulse_blocks):
        #                 pulse_sequence.extend(
        #                     [Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time),
        #                      Pulse('RF_switch', section_begin_time + 1 * tau - RF_pi_time / 2, RF_pi_time),
        #                      Pulse(microwave_channel_1, section_begin_time + 2 * tau - pi_time / 2, pi_time),
        #                      Pulse('RF_switch', section_begin_time + 2 * tau - RF_pi_time / 2, RF_pi_time),
        #                      Pulse(microwave_channel_2, section_begin_time + 3 * tau - pi_time / 2, pi_time),
        #                      Pulse('RF_switch', section_begin_time + 3 * tau - RF_pi_time / 2, RF_pi_time),
        #                      Pulse(microwave_channel_1, section_begin_time + 4 * tau - pi_time / 2, pi_time),
        #                      Pulse('RF_switch', section_begin_time + 4 * tau - RF_pi_time / 2, RF_pi_time),
        #                      Pulse(microwave_channel_1, section_begin_time + 5 * tau - pi_time / 2, pi_time),
        #                      Pulse('RF_switch', section_begin_time + 5 * tau - RF_pi_time / 2, RF_pi_time),
        #                      Pulse(microwave_channel_2, section_begin_time + 6 * tau - pi_time / 2, pi_time),
        #                      Pulse('RF_switch', section_begin_time + 6 * tau - RF_pi_time / 2, RF_pi_time),
        #                      Pulse(microwave_channel_1, section_begin_time + 7 * tau - pi_time / 2, pi_time),
        #                      Pulse('RF_switch', section_begin_time + 7 * tau - RF_pi_time / 2, RF_pi_time),
        #                      Pulse(microwave_channel_2, section_begin_time + 8 * tau - pi_time / 2, pi_time),
        #                      Pulse('RF_switch', section_begin_time + 8 * tau - RF_pi_time / 2, RF_pi_time)
        #                      ])
        #                 section_begin_time += 8 * tau
        #
        #             # the second 3*pi/2 pulse and readout
        #             pulse_sequence.extend([Pulse(microwave_channel_1, section_begin_time + tau / 2, three_pi_half_time),
        #                                    Pulse('laser',
        #                                          section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout,
        #                                          nv_reset_time),
        #                                    Pulse('apd_readout',
        #                                          section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + delay_readout,
        #                                          meas_time)
        #                                    ])
        #
        #         pulse_sequences.append(pulse_sequence)
        #
        #     print('number of sequences before validation ', len(pulse_sequences))
        #     return pulse_sequences, self.settings['num_averages'], tau_list, meas_time

    # This is not a real CPMG seq because there is not I/Q option...
    #     elif self.settings['decoupling_seq']['type'] == 'CPMG':
    #         for tau_total in tau_list:
    #             tau = tau_total / (1 * self.settings['decoupling_seq']['num_of_pulse_blocks'])
    #             pulse_sequence = []
    #
    #             # ECHO SEQUENCE:
    #             start_of_first_HE = laser_off_time
    #
    #             # the first pi/2 pulse
    #             pulse_sequence = [Pulse(microwave_channel_1, start_of_first_HE, pi_half_time)]
    #
    #             section_begin_time = start_of_first_HE + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
    #
    #             for i in range(0, number_of_pulse_blocks):
    #                 pulse_sequence.append(Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time))
    #                 section_begin_time += 1 * tau
    #
    #             # the second pi/2 pulse and readout
    #             pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + tau / 2, pi_half_time))
    #             pulse_sequence.append(Pulse('laser', section_begin_time + tau / 2 + pi_half_time + delay_mw_readout,
    #                                          nv_reset_time))
    #             pulse_sequence.append(Pulse('apd_readout',
    #                                          section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + delay_readout,
    #                                          meas_time))
    #
    #             start_of_second_HE = section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time
    #
    #             # the first pi/2 pulse
    #             pulse_sequence.append(Pulse(microwave_channel_1, start_of_second_HE, pi_half_time))
    #
    #             section_begin_time = start_of_second_HE + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
    #
    #             for i in range(0, number_of_pulse_blocks):
    #                 pulse_sequence.append(Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time))
    #                 section_begin_time += 1 * tau
    #
    #             # the second 3*pi/2 pulse and readout
    #             pulse_sequence.append(Pulse(microwave_channel_1, section_begin_time + tau / 2, three_pi_half_time))
    #             pulse_sequence.append(Pulse('laser',
    #                                          section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout,
    #                                          nv_reset_time))
    #             pulse_sequence.append(Pulse('apd_readout',
    #                                          section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + delay_readout,
    #                                          meas_time))
    #
    #             # Turn on APD switch all the time
    #             pulse_sequence.append(Pulse('apd_switch', 0,
    #                                         section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + nv_reset_time))
    #
    #
    #             # # DEER SEQUENCE
    #             # if self.settings['decoupling_seq']['do_deer']:
    #             #     start_of_DEER = section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time
    #             #     # the first pi/2 pulse
    #             #     pulse_sequence.extend([Pulse(microwave_channel_1, start_of_DEER, pi_half_time)])
    #             #     section_begin_time = start_of_DEER + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
    #             #     for i in range(0, number_of_pulse_blocks):
    #             #         pulse_sequence.extend(
    #             #             [Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time),
    #             #              Pulse('RF_switch', section_begin_time + 1 * tau - RF_pi_time / 2, RF_pi_time)
    #             #              ])
    #             #         section_begin_time += 1 * tau
    #             #
    #             #     # the second pi/2 pulse and readout
    #             #     pulse_sequence.extend([Pulse(microwave_channel_1, section_begin_time + tau / 2, pi_half_time),
    #             #                            Pulse('laser',
    #             #                                  section_begin_time + tau / 2 + pi_half_time + delay_mw_readout,
    #             #                                  nv_reset_time),
    #             #                            Pulse('apd_readout',
    #             #                                  section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + delay_readout,
    #             #                                  meas_time)
    #             #                            ])
    #             #
    #             #     start_of_second_DEER = section_begin_time + tau / 2 + pi_half_time + delay_mw_readout + nv_reset_time + laser_off_time
    #             #     pulse_sequence.extend([Pulse(microwave_channel_1, start_of_second_DEER, pi_half_time)])
    #             #     section_begin_time = start_of_second_DEER + pi_half_time - tau / 2  # for the first pulse, only wait tau/2
    #             #
    #             #     for i in range(0, number_of_pulse_blocks):
    #             #         pulse_sequence.extend(
    #             #             [Pulse(microwave_channel_2, section_begin_time + 1 * tau - pi_time / 2, pi_time),
    #             #              Pulse('RF_switch', section_begin_time + 1 * tau - RF_pi_time / 2, RF_pi_time)
    #             #              ])
    #             #         section_begin_time += 1 * tau
    #             #
    #             #     # the second 3*pi/2 pulse and readout
    #             #     pulse_sequence.extend([Pulse(microwave_channel_1, section_begin_time + tau / 2, three_pi_half_time),
    #             #                            Pulse('laser',
    #             #                                  section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout,
    #             #                                  nv_reset_time),
    #             #                            Pulse('apd_readout',
    #             #                                  section_begin_time + tau / 2 + three_pi_half_time + delay_mw_readout + delay_readout,
    #             #                                  meas_time)
    #             #                            ])
    #
    #             pulse_sequences.append(pulse_sequence)
    #
    #         print('number of sequences before validation ', len(pulse_sequences))
    #         # return pulse_sequences, self.settings['num_averages'], tau_list, meas_time
    #         return pulse_sequences, tau_list, meas_time

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

        if 'fits' in data.keys() is not None and data['fits'] is not None:

            # counts = (-data['counts'][:,1] + data['counts'][:,0])/ (data['counts'][:,0] + data['counts'][:,1])
            tau = data['tau']
            fits = data['fits']

            # axislist[0].plot(tau, counts)
            axislist[0].plot(tau, data['norm_echo'])
            tauinterp = np.linspace(np.min(tau), np.max(tau), 100)
            axislist[0].plot(tauinterp, exp_offset(tauinterp, fits[0], fits[1], fits[2]))
            axislist[0].set_title('T2 decay time (simple exponential, p = 1): {:2.1f} ns'.format(fits[1]))

            axislist[0].set_title(
                'T2 measurement \n {:s} {:d} block(s)\n averaged over {:d} blocks \n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz, \n T2 (simple expo, p = 1) = {:2.1f} ns'.format(
                    self.settings['decoupling_seq']['type'],
                    self.settings['decoupling_seq']['num_of_pulse_blocks'], self.avg_block_number,
                    self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency'] * 1e-9,
                    fits[1]))
            axislist[0].legend(labels=('Echo', 'Exp Fit'), fontsize=10)
            axislist[0].set_ylabel('fluorescence contrast')
            axislist[0].set_xlabel('tau [ns]')

        elif 'exp_finished' in data.keys() is not None and data['exp_finished'] == 1:
            if 'norm_echo' in data.keys() is not None and data['norm_echo'] is not None:
                tau = data['tau']
                axislist[0].plot(tau, data['norm_echo'],'d')
                axislist[0].set_title(
                    '(final plot) T2 measurement \n {:s} {:d} block(s)\n averaged over {:d} blocks \n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz'.format(
                        self.settings['decoupling_seq']['type'],
                        self.settings['decoupling_seq']['num_of_pulse_blocks'],
                        self.avg_block_number,
                        self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency'] * 1e-9))
                axislist[0].legend(labels=('Echo'), fontsize=10)
                axislist[0].set_ylabel('fluorescence contrast')
                axislist[0].set_xlabel('tau [ns]')
            else:
                super(PDD_RnS, self)._plot(axislist)
                axislist[0].legend(labels=('Echo up {:.0f}kcps'.format(np.mean(data['counts'][:, 0])),
                                           'Echo down {:.0f}kcps'.format(np.mean(data['counts'][:, 1]))), fontsize=10)
                axislist[0].set_title(
                    '(final plot) T2 measurement \n {:s} {:d} block(s) \n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz'.format(
                        self.settings['decoupling_seq']['type'], self.settings['decoupling_seq']['num_of_pulse_blocks'],
                        self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency'] * 1e-9))

        else:
            super(PDD_RnS, self)._plot(axislist)
            axislist[0].legend(labels=('Echo up {:.0f}kcps'.format(np.mean(data['counts'][:, 0])),
                                       'Echo down {:.0f}kcps'.format(np.mean(data['counts'][:, 1]))), fontsize=10)
            axislist[0].set_title(
                '(initial plot) T2 measurement \n {:s} {:d} block(s) \n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz'.format(
                    self.settings['decoupling_seq']['type'], self.settings['decoupling_seq']['num_of_pulse_blocks'],
                    self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency'] * 1e-9))

    def _update_plot(self, axislist):

        if len(axislist[0].lines) == 0:
            self._plot(axislist)
            return
        super(PDD_RnS, self)._update_plot(axislist)
        axislist[0].legend(labels=('Echo up {:.0f}kcps'.format(np.mean(self.data['counts'][:, 0])),
                                   'Echo down {:.0f}kcps'.format(np.mean(self.data['counts'][:, 1]))), fontsize=10)
        axislist[0].set_title('(updating plot) T2 measurement \n {:s} {:d} block(s) \n mw-power:{:.0f}dBm, mw_freq:{:.4f} GHz'.format(
            self.settings['decoupling_seq']['type'], self.settings['decoupling_seq']['num_of_pulse_blocks'],
            self.settings['mw_pulses']['mw_power'], self.settings['mw_pulses']['mw_frequency'] * 1e-9))


# class HahnEcho_bothIQ(PulsedExperimentBaseScript): # ER 20181013
#     """
# This script runs a Hahn echo on the NV to find the Hahn echo T2. To symmetrize the sequence between the 0 and +/-1 state we reinitialize every time.
# The difference between this script and HahnEcho is that the phases of the pulses should be xyx, instead of xxx; i.e., we use BOTH IQ channels now.
# The channel 'microwave_channel' in the settings corresponds to the pi/2 and 3pi/2 pulses. The other one the user chooses will be for the pi pulse.
# We do this to check if we can extend T2 with something like this, which may help mitigate pulse errors.
#
#     """
#     _DEFAULT_SETTINGS = [
#         Parameter('mw_pulses', [
#             Parameter('mw_power', -45.0, float, 'microwave power in dB'),
#             Parameter('mw_frequency', 2.87e9, float, 'microwave frequency in Hz'),
#             Parameter('microwave_channel', 'i', ['i', 'q'], 'Channel to use for mw pi/2 pulses (other one will be for the pi pulses)'),
#             Parameter('pi_pulse_time', 50.0, float, 'time duration of a pi pulse (in ns)'),
#             Parameter('pi_half_pulse_time', 25.0, float, 'time duration of a pi/2 pulse (in ns)'),
#             Parameter('3pi_half_pulse_time', 75.0, float, 'time duration of a 3pi/2 pulse (in ns)')
#         ]),
#         Parameter('tau_times', [
#             Parameter('min_time', 500, float, 'minimum time between pi pulses'),
#             Parameter('max_time', 10000, float, 'maximum time between pi pulses'),
#             Parameter('time_step', 5, [2.5, 5, 10, 20, 50, 100, 200, 500, 1000, 10000, 100000, 500000],
#                   'time step increment of time between pi pulses (in ns)')
#         ]),
#         Parameter('read_out', [
#             Parameter('meas_time', 250, float, 'measurement time after rabi sequence (in ns)'),
#             Parameter('nv_reset_time', 1750, int, 'time with laser on to reset state'),
#             Parameter('laser_off_time', 1000, int,
#                       'minimum laser off time before taking measurements (ns)'),
#             Parameter('delay_mw_readout', 1000, int, 'delay between mw and readout (in ns)'),
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
#         super(HahnEcho_bothIQ, self)._function(self.data)
#
#         counts = (- self.data['counts'][:, 1] + self.data['counts'][:,0]) / (self.data['counts'][:,1] + self.data['counts'][:, 0])
#         tau = self.data['tau']
#
#         try:
#             fits = fit_exp_decay(tau, counts, offset = True, verbose = True)
#             self.data['fits'] = fits
#         except:
#             self.data['fits'] = None
#             self.log('t2 fit failed')
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
#         tau_list = np.arange(self.settings['tau_times']['min_time'], self.settings['tau_times']['max_time'],self.settings['tau_times']['time_step'])
#         tau_list = np.ndarray.tolist(tau_list) # 20180731 ER convert to list
#
#         # ignore the sequence if the mw-pulse is shorter than 15ns (0 is ok because there is no mw pulse!)
#         tau_list = [x for x in tau_list if x == 0 or x >= 15]
#
#         nv_reset_time = self.settings['read_out']['nv_reset_time']
#         delay_readout = self.settings['read_out']['delay_readout']
#         microwave_channel = 'microwave_' + self.settings['mw_pulses']['microwave_channel']
#         if self.settings['mw_pulses']['microwave_channel'] == 'i':
#             mw_chan_pi = 'q'
#         else:
#             mw_chan_pi = 'i'
#         microwave_channel_pi = 'microwave_' + mw_chan_pi
#         pi_time = self.settings['mw_pulses']['pi_pulse_time']
#         pi_half_time = self.settings['mw_pulses']['pi_half_pulse_time']
#         three_pi_half_time = self.settings['mw_pulses']['3pi_half_pulse_time']
#
#         laser_off_time = self.settings['read_out']['laser_off_time']
#         meas_time = self.settings['read_out']['meas_time']
#         delay_mw_readout = self.settings['read_out']['delay_mw_readout']
#
#         for tau in tau_list:
#             pulse_sequence = \
#             [
#                 Pulse(microwave_channel, laser_off_time, pi_half_time),
#                 Pulse(microwave_channel_pi, laser_off_time + pi_half_time + tau - pi_time/2., pi_time),
#                 Pulse(microwave_channel, laser_off_time + pi_half_time + tau + tau, pi_half_time)
#             ]
#
#             end_of_first_HE = laser_off_time + pi_half_time + tau + tau + pi_half_time
#
#             pulse_sequence += [
#                  Pulse('laser', end_of_first_HE + delay_mw_readout, nv_reset_time),
#                  Pulse('apd_readout', end_of_first_HE + delay_mw_readout + delay_readout, meas_time),
#                  ]
#
#             start_of_second_HE = end_of_first_HE + delay_mw_readout + nv_reset_time + laser_off_time
#
#             pulse_sequence += \
#             [
#                 Pulse(microwave_channel, start_of_second_HE, pi_half_time),
#                 Pulse(microwave_channel_pi, start_of_second_HE + pi_half_time + tau - pi_time/2., pi_time),
#                 Pulse(microwave_channel, start_of_second_HE + pi_half_time + tau + tau, three_pi_half_time)
#             ]
#
#             end_of_second_HE = start_of_second_HE + pi_half_time + tau + tau + three_pi_half_time
#
#             pulse_sequence += [
#                 Pulse('laser', end_of_second_HE + delay_mw_readout, nv_reset_time),
#                 Pulse('apd_readout', end_of_second_HE + delay_mw_readout + delay_readout, meas_time)
#             ]
#             # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
#             # if tau == 0 or tau>=15:
#             pulse_sequences.append(pulse_sequence)
#
#         return pulse_sequences, tau_list, meas_time
#
# class HahnEchoManyNVs(Script):
#     _DEFAULT_SETTINGS = [
#         Parameter('esr_peak', 'upper', ['upper', 'lower', 'both'], 'if ESR fits two peaks, defines which one to use')
#     ]
#     _INSTRUMENTS = {}
#     _SCRIPTS = {'select_NVs': SelectPoints, 'ESR': ESR, 'Rabi': Rabi, 'HE': HahnEcho}
#
#     def _function(self):
#         for num, nv_loc in enumerate(self.scripts['select_NVs'].data['nv_locations']):
#             if self._abort:
#                 break
#             find_NV_rabi = self.scripts['Rabi'].scripts['find_nv']
#             find_NV_rabi.settings['initial_point']['x'] = nv_loc[0]
#             find_NV_rabi.settings['initial_point']['y'] = nv_loc[1]
#             find_NV_rabi.run()
#             self.scripts['ESR'].settings['tag'] = 'esr_NV' + str(num)
#             self.scripts['ESR'].run()
#             fit_params = self.scripts['ESR'].data['fit_params']
#             if fit_params is None:
#                 continue
#             if len(fit_params) == 4:
#                 freqs = [fit_params[2]]
#             elif len(fit_params == 6):
#                 if self.settings['esr_peak'] == 'lower':
#                     freqs = [fit_params[4]]
#                 elif self.settings['esr_peak'] == 'upper':
#                     freqs = [fit_params[5]]
#                 elif self.settings['esr_peak'] == 'both':
#                     freqs = [fit_params[4], fit_params[5]]
#             for freq in freqs:
#                 if self._abort:
#                     break
#                 print('running rabi')
#                 rabi = self.scripts['Rabi']
#                 rabi.settings['tag'] = 'rabi_NV' + str(num)
#                 rabi.settings['mw_pulses']['mw_frequency'] = float(freq)
#                 print('about to run rabi')
#                 rabi.run()
#                 rabi_fit = rabi.data['fits']
#                 if rabi_fit is None:
#                     continue
#                 pi_time = abs((np.pi - rabi_fit[2])/rabi_fit[1])
#                 pi_time = min(max(np.round(pi_time / 2.5) * 2.5, 15.), 300.) #round to nearest 2.5 ns
#                 pi_half_time = min(max((np.pi / 2 - rabi_fit[2]) / rabi_fit[1], 15.), 300)
#                 three_pi_half_time = min(max((3 * np.pi / 2 - rabi_fit[2]) / rabi_fit[1], 15.), 300)
#                 find_NV_HE = self.scripts['HE'].scripts['find_nv']
#                 find_NV_HE.settings['initial_point']['x'] = find_NV_rabi.data['maximum_point']['x']
#                 find_NV_HE.settings['initial_point']['y'] = find_NV_rabi.data['maximum_point']['y']
#                 HE = self.scripts['HE']
#                 HE.settings['mw_pulses']['mw_frequency'] = float(freq)
#                 HE.settings['mw_pulses']['pi_time'] = float(pi_time)
#                 HE.settings['mw_pulses']['pi_half_time'] = float(pi_half_time)
#                 HE.settings['mw_pulses']['3pi_half_time'] = float(three_pi_half_time)
#                 HE.settings['tag'] = 'HE' + '_NV' + str(num)
#                 HE.run()
#
#     def plot(self, figure_list):
#         if self._current_subscript_stage is not None:
#             if self._current_subscript_stage['current_subscript'] is not None:
#                 self._current_subscript_stage['current_subscript'].plot(figure_list)
#
#
#     def skip_next(self):
#         for script in self.scripts.values():
#             script.stop()
