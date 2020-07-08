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
from b26_toolkit.instruments import NI6353, LISE607RTPulseBlaster, R8SMicrowaveGenerator, Pulse, AgilentMicrowaveGenerator, AgilentMicrowaveGeneratorII
from b26_toolkit.plotting.plots_1d import plot_esr, plot_pulses
from pylabcontrol.core import Parameter
import random

class PulsedESR_RnS(PulsedExperimentBaseScript): # ER 20170616 - wrote for symmetry between 0 and -1 state
    """
    This script applies a microwave pulse at fixed power and durations for varying frequencies.
    The Keysight (Agilent) signal generator N9310A is used as MW source.
    Uses double_init scheme.

    Modified from B26 script. Do full averages on one frequency before going to another frequency.

    -- Edited by Ziwei Qiu 8/31/2019
    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_power', -45.0, float, 'microwave power in dB'),
        Parameter('tau_mw', 80, float, 'the time duration of the microwaves (in ns)'),
        Parameter('freq_start', 2.82e9, float, 'start frequency of scan in Hz'),
        Parameter('freq_stop', 2.92e9, float, 'end frequency of scan in Hz'),
        Parameter('freq_points', 50, int, 'number of frequencies in scan in Hz'),
        Parameter('microwave_channel', 'i', ['i', 'q'], 'Channel to use for mw pulses'),
        Parameter('read_out', [
            Parameter('meas_time', 200, float, 'measurement time after rabi sequence (in ns)'),
            Parameter('nv_reset_time', 2000, int, 'time with laser on to reset state'),
            Parameter('laser_off_time', 1000, int,
                      'minimum laser off time before taking measurements (ns)'),
            Parameter('delay_mw_readout', 600, int, 'delay between mw and readout (in ns)'),
            Parameter('delay_readout', 550, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        ]),
        Parameter('num_averages', 1000000, int, 'number of averages'),

    ]

    # _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator}
    _INSTRUMENTS = {'daq': NI6353, 'PB': LISE607RTPulseBlaster, 'mw_gen_a': AgilentMicrowaveGenerator}


    def _function(self):
        #COMMENT_ME

        ## Turn off green light (the pulse blaster will pulse it on when needed)
        self.instruments['PB']['instance'].update({'laser': {'status': False}})
        self.instruments['PB']['instance'].update({'apd_switch': {'status': False}})

        # self.instruments['mw_gen_a']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen_a']['instance'].update({'enable_IQ': True})
        self.instruments['mw_gen_a']['instance'].update({'enable_modulation': True})
        self.instruments['mw_gen_a']['instance'].update({'power_mode': 'CW'})
        self.instruments['mw_gen_a']['instance'].update({'amplitude': self.settings['freq_start']})
        self.instruments['mw_gen_a']['instance'].update({'freq_mode': 'CW'})
        self.instruments['mw_gen_a']['instance'].update({'enable_output': True})

        # assert self.settings['freq_start'] < self.settings['freq_stop']
        if self.settings['freq_start'] >= self.settings['freq_stop']:
            self.log('ATTENTION: freq_start should be smaller than freq_stop. Experiment stopped.')
            self._abort = True

        self.data = {'mw_frequencies': np.linspace(self.settings['freq_start'], self.settings['freq_stop'],
                                                   self.settings['freq_points']), 'esr_counts': []}

        for i, mw_frequency in enumerate(self.data['mw_frequencies']):
            self._loop_count = i
            self.instruments['mw_gen_a']['instance'].update({'frequency': float(mw_frequency)})
            super(PulsedESR_RnS, self)._function(self.data) # this might mean all averages are done on one freq before going to another freq
            self.data['esr_counts'].append(self.data['counts'])

    # def _calc_progress(self):
    #     #COMMENT_ME
    #     # todo: change to _calc_progress(self, index):
    #     progress = int(100. * (self._loop_count) / self.settings['freq_points'])
    #     return progress

    def _create_pulse_sequences(self):

        '''

        Returns: pulse_sequences, num_averages, tau_list
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
            meas_time: the width (in ns) of the daq measurement

        '''

        tau = self.settings['tau_mw']
        pulse_sequences = []
        tau_list = [tau]

        meas_time = self.settings['read_out']['meas_time']
        nv_reset_time = self.settings['read_out']['nv_reset_time']
        laser_off_time = self.settings['read_out']['laser_off_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        delay_readout = self.settings['read_out']['delay_readout']
        microwave_channel = 'microwave_' + self.settings['microwave_channel']

        for tau in tau_list:
            pulse_sequence = \
                [Pulse('laser', laser_off_time + tau + delay_mw_readout, nv_reset_time),
                 Pulse('apd_readout', laser_off_time + tau + delay_mw_readout + delay_readout, meas_time),
                 ]
            # if tau is 0 there is actually no mw pulse

            if tau > 0:
                pulse_sequence.append(
                    Pulse(microwave_channel, laser_off_time + tau + delay_mw_readout + nv_reset_time + laser_off_time,
                          tau))

            pulse_sequence.append(Pulse('laser',
                                        laser_off_time + tau + delay_mw_readout + nv_reset_time + laser_off_time + tau + delay_mw_readout,
                                        nv_reset_time))
            pulse_sequence.append(Pulse('apd_readout',
                                        laser_off_time + tau + delay_mw_readout + nv_reset_time + laser_off_time + tau + delay_mw_readout + delay_readout,
                                        meas_time))
            pulse_sequence.append(Pulse('apd_switch', 0,
                                        laser_off_time + tau + delay_mw_readout + nv_reset_time + laser_off_time + tau + delay_mw_readout + nv_reset_time))

            # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
            # if tau == 0 or tau>=15:
            pulse_sequences.append(pulse_sequence)

        return pulse_sequences, tau_list, meas_time

    def _plot(self, axes_list, data = None):
        '''
        Plot 1: self.data['tau'], the list of times specified for a given experiment, verses self.data['counts'], the data
        received for each time
        Plot 2: the pulse sequence performed at the current time (or if plotted statically, the last pulse sequence
        performed

        Args:
            axes_list: list of axes to write plots to (uses first 2)
            data (optional) dataset to plot, if not provided use self.data
        '''
        if data is None:
            data = self.data

        mw_frequencies = data['mw_frequencies']
        esr_counts = data['esr_counts']
        axis1 = axes_list[0]
        if not esr_counts == []:
            counts = esr_counts
            plot_esr(axis1, mw_frequencies[0:len(counts)], counts)
            axis1.hold(False)
        axis2 = axes_list[1]
        plot_pulses(axis2, self.pulse_sequences[0])

    def _update_plot(self, axes_list):
        mw_frequencies = self.data['mw_frequencies']
        esr_counts = self.data['esr_counts']
        axis1 = axes_list[0]
        if not esr_counts == []:
            counts = esr_counts
            plot_esr(axis1, mw_frequencies[0:len(counts)], counts)
            axis1.hold(False)
            # axis2 = axes_list[1]
            # update_pulse_plot(axis2, self.pulse_sequences[0])


class PulsedESR(PulsedExperimentBaseScript):  # ER 20170616 - wrote for symmetry between 0 and -1 state
    """
    This script applies a microwave pulse at fixed power and durations for varying frequencies.
    The Keysight (Agilent) signal generator N9310A is used as MW source.
    Uses double_init scheme.
    This script first goes through all frequencies before doing another averaging block.

    -- Edited by Ziwei Qiu 8/31/2019
    """
    _DEFAULT_SETTINGS = [
        Parameter('mw_power', -45.0, float, 'microwave power in dB'),
        Parameter('tau_mw', 80, float, 'the time duration of the microwaves (in ns)'),
        Parameter('freq_start', 2.87e9, float, 'start frequency of scan'),
        Parameter('freq_stop', 3e7, float, 'end frequency of scan'),
        Parameter('range_type', 'center_range', ['start_stop', 'center_range'],
                  'start_stop: freq. range from freq_start to freq_stop. center_range: centered at freq_start and width freq_stop'),
        Parameter('freq_points', 50, int, 'number of frequencies in scan'),

        Parameter('microwave_channel', 'i', ['i', 'q'], 'Channel to use for mw pulses'),
        # Parameter('read_out', [
        #     Parameter('meas_time', 200, float, 'measurement time after rabi sequence (in ns)'),
        #     Parameter('nv_reset_time', 2000, int, 'time with laser on to reset state'),
        #     Parameter('laser_off_time', 1000, int,
        #               'minimum laser off time before taking measurements (ns)'),
        #     Parameter('delay_mw_readout', 600, int, 'delay between mw and readout (in ns)'),
        #     Parameter('delay_readout', 550, int, 'delay between laser on and readout (given by spontaneous decay rate)')
        # ]),
        Parameter('min_pulse_dur', 15, [10, 15, 50], 'Minimum allowed pulse duration (ns)'),
        Parameter('num_averages', 1000000, int, 'number of averages'),

    ]
    _INSTRUMENTS = {'daq': NI6353, 'PB': LISE607RTPulseBlaster, 'mw_gen_a': AgilentMicrowaveGenerator}

    def _function(self):

        ## Turn off green light (the pulse blaster will pulse it on when needed)
        self.instruments['PB']['instance'].update({'laser': {'status': False}})
        self.instruments['PB']['instance'].update({'apd_switch': {'status': False}})

        # self.instruments['mw_gen_a']['instance'].update({'modulation_type': 'IQ'})
        self.instruments['mw_gen_a']['instance'].update({'enable_IQ': True})
        self.instruments['mw_gen_a']['instance'].update({'enable_modulation': True})
        self.instruments['mw_gen_a']['instance'].update({'power_mode': 'CW'})
        self.instruments['mw_gen_a']['instance'].update({'amplitude': self.settings['mw_power']})
        self.instruments['mw_gen_a']['instance'].update({'freq_mode': 'CW'})
        self.instruments['mw_gen_a']['instance'].update({'amplitude': self.settings['freq_start']})
        self.instruments['mw_gen_a']['instance'].update({'enable_output': True})

        valid = True

        # contruct the frequency array
        if self.settings['range_type'] == 'start_stop':
            if self.settings['freq_start'] > self.settings['freq_stop']:
                self.log('end freq. must be larger than start freq when range_type is start_stop. Abort script')
                valid = False
                self._abort = True

            if self.settings['freq_start'] < 9.0E3 or self.settings['freq_stop'] > 3.00E9:  # freq range of the Keysight N9310A
                self.log('start or stop frequency out of bounds')
                valid = False
                self._abort = True

            freq_values = np.linspace(self.settings['freq_start'], self.settings['freq_stop'],
                                      self.settings['freq_points'])
            self.data = {'mw_frequencies': freq_values}
            freq_range = max(freq_values) - min(freq_values)

        elif self.settings['range_type'] == 'center_range':
            if self.settings['freq_start'] < 2 * self.settings['freq_stop']:
                self.log(
                    'end freq. (range) must be smaller than 2x start freq (center) when range_type is center_range. Abort script')
                valid = False
                self._abort = True
            freq_values = np.linspace(self.settings['freq_start'] - self.settings['freq_stop'] / 2,
                                      self.settings['freq_start'] + self.settings['freq_stop'] / 2,
                                      self.settings['freq_points'])
            self.data = {'mw_frequencies': freq_values}
            freq_range = max(freq_values) - min(freq_values)

            if self.settings['freq_stop'] > 1e9:
                self.log('freq_stop (range) is quite large --- did you mean to set \'range_type\' to \'start_stop\'? ')
        else:
            self.log('unknown range parameter. Abort script')
            valid = False
            self._abort = True


        if valid:
            self.avg_block_number = super(PulsedESR, self)._function()

            # this is to flag when the experiment is done or stopped and this allows proper saving even fit fails
            self.data['exp_finished'] = 1
            self.data['avg_block_number'] = self.avg_block_number

        # turn off laser, apd_switch and MW (just in case)
        self.instruments['PB']['instance'].update({'laser': {'status': False}})
        self.instruments['PB']['instance'].update({'apd_switch': {'status': False}})
        self.instruments['mw_gen_a']['instance'].update({'enable_output': False})
        self.instruments['mw_gen_a']['instance'].update({'enable_IQ': False})
        self.instruments['mw_gen_a']['instance'].update({'enable_modulation': False})

        # no fitting for now, just normalize the data
        if 'counts' in self.data.keys() and 'tau' in self.data.keys() and not 0 in (self.data['counts'][:, 0]) :

            self.data['norm_esr'] = self.data['counts'][:, 1] / self.data['counts'][:, 0]
            tau = self.data['tau']
            counts = self.data['norm_esr']
            fit_success = False


    def _create_pulse_sequences(self):

        '''

        Returns: pulse_sequences, num_averages, tau_list
            pulse_sequences: a list of pulse sequences, each corresponding to a different time 'tau' that is to be
            scanned over. Each pulse sequence is a list of pulse objects containing the desired pulses. Each pulse
            sequence must have the same number of daq read pulses
            num_averages: the number of times to repeat each pulse sequence
            tau_list: the list of times tau, with each value corresponding to a pulse sequence in pulse_sequences
            meas_time: the width (in ns) of the daq measurement

        '''
        pulse_sequences = []
        freq_list = self.data['mw_frequencies']
        print('freq_list:', freq_list)

        tau = self.settings['tau_mw']
        meas_time = self.settings['read_out']['meas_time']
        nv_reset_time = self.settings['read_out']['nv_reset_time']
        laser_off_time = self.settings['read_out']['laser_off_time']
        delay_mw_readout = self.settings['read_out']['delay_mw_readout']
        delay_readout = self.settings['read_out']['delay_readout']
        microwave_channel = 'microwave_' + self.settings['microwave_channel']

        for freq in freq_list:
            pulse_sequence = \
                [Pulse('laser', laser_off_time + tau + delay_mw_readout, nv_reset_time),
                 Pulse('apd_readout', laser_off_time + tau + delay_mw_readout + delay_readout, meas_time),
                 ]
            # if tau is 0 there is actually no mw pulse

            if tau > 0:
                pulse_sequence.append(
                    Pulse(microwave_channel, laser_off_time + tau + delay_mw_readout + nv_reset_time + laser_off_time,
                          tau))

            pulse_sequence.append(Pulse('laser',
                                        laser_off_time + tau + delay_mw_readout + nv_reset_time + laser_off_time + tau + delay_mw_readout,
                                        nv_reset_time))
            pulse_sequence.append(Pulse('apd_readout',
                                        laser_off_time + tau + delay_mw_readout + nv_reset_time + laser_off_time + tau + delay_mw_readout + delay_readout,
                                        meas_time))
            pulse_sequence.append(Pulse('apd_switch', 0,
                                        laser_off_time + tau + delay_mw_readout + nv_reset_time + laser_off_time + tau + delay_mw_readout + nv_reset_time))

            # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
            # if tau == 0 or tau>=15:
            pulse_sequences.append(pulse_sequence)

        print('number of sequences before validation ', len(pulse_sequences))
        # return pulse_sequences, self.settings['num_averages'], tau_list, meas_time
        return pulse_sequences, freq_list, meas_time

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
            self.instruments['mw_gen_a']['instance'].update({'frequency': float(tau_list[rand_index])})


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
            if 'norm_esr' in data.keys() is not None and data['norm_esr'] is not None:
                mean_fluor = np.mean(data['counts'][:, 0])
                axislist[0].plot(data['tau'] * 1e-6, data['norm_esr']) # convert to MHz
                axislist[0].legend(labels=('normalized ESR data'), fontsize=10)
                axislist[0].set_title('(final plot) Keysight N9310A, MW channel is ' + self.settings[
                    'microwave_channel'] + '\n after {:d} averages \n Ref fluor: {:0.1f}kcps \n mw-power:{:.0f} dBm, tau:{:.2f} ns'.format(
                    avg_block_number, mean_fluor, self.settings['mw_power'], self.settings['tau_mw']))
                axislist[0].set_ylabel('fluorescence contrast')
                axislist[0].set_xlabel('MW frequency [MHz]')

            else:
                axislist[0].plot(data['tau'] * 1e-6, data['counts'][:, 0])
                axislist[0].plot(data['tau'] * 1e-6, data['counts'][:, 1])

                axislist[0].legend(labels=('Reference {:.0f}kcps'.format(np.mean(data['counts'][:, 0])),
                                           'ESR {:.0f}kcps'.format(np.mean(data['counts'][:, 1]))), fontsize=10)
                axislist[0].set_title('(final plot) Keysight N9310A, MW channel is ' + self.settings[
                    'microwave_channel'] + '\n after {:d} averages \nmw-power:{:.0f} dBm, tau:{:.2f} ns'.format(
                    avg_block_number, self.settings['mw_power'], self.settings['tau_mw']))
                axislist[0].set_ylabel('fluorescence contrast')
                axislist[0].set_xlabel('MW frequency [MHz]')

        else:
            super(PulsedESR, self)._plot(axislist)
            axislist[0].set_xlabel('MW frequency [GHz]')
            axislist[0].legend(labels=('Reference {:.0f}kcps'.format(np.mean(data['counts'][:, 0])),
                                       'ESR {:.0f}kcps'.format(np.mean(data['counts'][:, 1]))), fontsize=10)
            axislist[0].set_title('(initial plot) Keysight N9310A, MW channel is ' + self.settings[
                'microwave_channel'] + '\n after {:d} averages \nmw-power:{:.0f} dBm, tau:{:.2f} ns'.format(
                avg_block_number, self.settings['mw_power'], self.settings['tau_mw']))

    def _update_plot(self, axislist):

        if len(axislist[0].lines) == 0:
            self._plot(axislist)
            return

        super(PulsedESR, self)._update_plot(axislist)
        axislist[0].legend(labels=('Reference {:.0f}kcps'.format(np.mean(self.data['counts'][:, 0])),
                                   'ESR {:.0f}kcps'.format(np.mean(self.data['counts'][:, 1]))), fontsize=10)

        axislist[0].set_title('(updating plot) Keysight N9310A, MW channel is ' + self.settings[
            'microwave_channel'] + '\nmw-power:{:.0f} dBm, tau:{:.2f} ns'.format(
            self.settings['mw_power'], self.settings['tau_mw']))

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
            axis0.set_title('Pulse Visualization for Minimum Frequency = {:.2f} MHz'.format(tau_list[0]*1e-6))
            plot_pulses(axis1, pulse_sequences[-1])
            axis1.set_title('Pulse Visualization for Maximum Frequency = {:.2f} MHz'.format(tau_list[-1]*1e-6))
        else:
            print('no pulse sequences passed validation!!!')











# The following is from old B26 script
# class PulsedESR(PulsedExperimentBaseScript): # ER 20170616 - wrote for symmetry between 0 and -1 state
#     """
# This script applies a microwave pulse at fixed power and durations for varying frequencies.
# Uses double_init scheme.
#
#     """
#     _DEFAULT_SETTINGS = [
#         Parameter('mw_power', -45.0, float, 'microwave power in dB'),
#         Parameter('tau_mw', 80, float, 'the time duration of the microwaves (in ns)'),
#         Parameter('num_averages', 1000000, int, 'number of averages'),
#         Parameter('freq_start', 2.82e9, float, 'start frequency of scan in Hz'),
#         Parameter('freq_stop', 2.92e9, float, 'end frequency of scan in Hz'),
#         Parameter('freq_points', 100, int, 'number of frequencies in scan in Hz'),
#         Parameter('read_out', [
#             Parameter('meas_time', 250, float, 'measurement time after rabi sequence (in ns)'),
#             Parameter('nv_reset_time', 1750, int, 'time with laser on to reset state'),
#             Parameter('laser_off_time', 1000, int,
#                       'minimum laser off time before taking measurements (ns)'),
#             Parameter('delay_mw_readout', 100, int, 'delay between mw and readout (in ns)'),
#             Parameter('delay_readout', 30, int, 'delay between laser on and readout (given by spontaneous decay rate)')
#         ]),
#     ]
#
#     _INSTRUMENTS = {'NI6259': NI6259, 'NI9402': NI9402, 'PB': B26PulseBlaster, 'mw_gen': MicrowaveGenerator}
#
#     def _function(self):
#         #COMMENT_ME
#         self.instruments['mw_gen']['instance'].update({'modulation_type': 'IQ'})
#         self.instruments['mw_gen']['instance'].update({'amplitude': self.settings['mw_power']})
#         self.instruments['mw_gen']['instance'].update({'enable_output': True})
#
#         assert self.settings['freq_start'] < self.settings['freq_stop']
#
#         self.data = {'mw_frequencies': np.linspace(self.settings['freq_start'], self.settings['freq_stop'],
#                                                    self.settings['freq_points']), 'esr_counts': []}
#
#         for i, mw_frequency in enumerate(self.data['mw_frequencies']):
#             self._loop_count = i
#             self.instruments['mw_gen']['instance'].update({'frequency': float(mw_frequency)})
#             super(PulsedESR, self)._function(self.data)
#             self.data['esr_counts'].append(self.data['counts'])
#
#     # def _calc_progress(self):
#     #     #COMMENT_ME
#     #     # todo: change to _calc_progress(self, index):
#     #     progress = int(100. * (self._loop_count) / self.settings['freq_points'])
#     #     return progress
#
#     def _plot(self, axes_list, data = None):
#         '''
#         Plot 1: self.data['tau'], the list of times specified for a given experiment, verses self.data['counts'], the data
#         received for each time
#         Plot 2: the pulse sequence performed at the current time (or if plotted statically, the last pulse sequence
#         performed
#
#         Args:
#             axes_list: list of axes to write plots to (uses first 2)
#             data (optional) dataset to plot, if not provided use self.data
#         '''
#         if data is None:
#             data = self.data
#
#         mw_frequencies = data['mw_frequencies']
#         esr_counts = data['esr_counts']
#         axis1 = axes_list[0]
#         if not esr_counts == []:
#             counts = esr_counts
#             plot_esr(axis1, mw_frequencies[0:len(counts)], counts)
#             axis1.hold(False)
#         axis2 = axes_list[1]
#         plot_pulses(axis2, self.pulse_sequences[0])
#
#     def _update_plot(self, axes_list):
#         mw_frequencies = self.data['mw_frequencies']
#         esr_counts = self.data['esr_counts']
#         axis1 = axes_list[0]
#         if not esr_counts == []:
#             counts = esr_counts
#             plot_esr(axis1, mw_frequencies[0:len(counts)], counts)
#             axis1.hold(False)
#             # axis2 = axes_list[1]
#             # update_pulse_plot(axis2, self.pulse_sequences[0])
#
#     def _create_pulse_sequences(self):
#
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
#
#         tau = self.settings['tau_mw']
#         pulse_sequences = []
#         tau_list = [tau]
#
#         nv_reset_time = self.settings['read_out']['nv_reset_time']
#         delay_readout = self.settings['read_out']['delay_readout']
#         microwave_channel = 'microwave_i'
#
#         laser_off_time = self.settings['read_out']['laser_off_time']
#         meas_time = self.settings['read_out']['meas_time']
#         delay_mw_readout = self.settings['read_out']['delay_mw_readout']
#
#         for tau in tau_list:
#             pulse_sequence = \
#                 [Pulse('laser', laser_off_time + tau + 2 * 40, nv_reset_time)]
#                #  Pulse('apd_readout', laser_off_time + tau + 2 * 40 + delay_readout, meas_time),
#                #  ]
#             # if tau is 0 there is actually no mw pulse
#             if tau > 0:
#                 pulse_sequence += [
#                     Pulse(microwave_channel, laser_off_time + tau + 2 * 40 + nv_reset_time + laser_off_time, tau)]
#
#             pulse_sequence += [
#                 Pulse('laser',
#                       laser_off_time + tau + 2 * 40 + nv_reset_time + laser_off_time + tau + 2 * 40 + delay_mw_readout,
#                       nv_reset_time),
#                 Pulse('apd_readout',
#                       laser_off_time + tau + 2 * 40 + nv_reset_time + laser_off_time + tau + 2 * 40 + delay_mw_readout + delay_readout,
#                       meas_time)
#             ]
#             # ignore the sequence is the mw is shorter than 15ns (0 is ok because there is no mw pulse!)
#             # if tau == 0 or tau>=15:
#             pulse_sequences.append(pulse_sequence)
#
#         return pulse_sequences, tau_list, self.settings['read_out']['meas_time']