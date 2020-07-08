"""
    This file is part of b26_toolkit, a PyLabControl add-on for experiments in Harvard LISE B26.
    Copyright (C) <2016>  Arthur Safira, Jan Gieseler, Aaron Kabcenell

    Foobar is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Foobar is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
"""

import numpy as np
from copy import deepcopy
import time
from scipy.signal import savgol_filter


from pylabcontrol.core import Script, Parameter
from b26_toolkit.plotting.plots_1d import plot_counts, plot_counts_vs_pos
from b26_toolkit.scripts.galvo_scan.confocal_scan_LISE607RT import ConfocalScan
from b26_toolkit.scripts.set_laser import SetConfocal
from b26_toolkit.instruments import LISE607RTPulseBlaster, NI6353

class optimize(Script):
    """
                optimize x, y and z
                written by ZQ 2/12/2019 3:29 pm

    """
    _DEFAULT_SETTINGS = [
            Parameter('optimizing_x', False, bool, 'if true, opmitimize x'),
            Parameter('optimizing_y', False, bool, 'if true, opmitimize y'),
            Parameter('optimizing_z', True, bool, 'if true, opmitimize z'),

            Parameter('sweep_range',
                      [Parameter('x', 0.3, float, 'x voltage range for optimizing scan'),
                       Parameter('y', 0.3, float, 'y voltage range for optimizing scan'),
                       Parameter('z', 0.6, float, 'z voltage range for optimizing scan')
                       ]),
            Parameter('num_points',
                      [Parameter('x', 125, int, 'number of x points to scan'),
                       Parameter('y', 125, int, 'number of y points to scan'),
                       Parameter('z', 41, int, 'number of z points to scan')
                       ]),
            # Parameter('smoothing',
            #           [Parameter('window_size'),
            #            [Parameter('x', 21, int, 'window size for savitzky_golay filtering of data x, must be odd'),
            #            Parameter('y', 21, int, 'window size for savitzky_golay filtering of data y, must be odd'),
            #            Parameter('z', 11, int, 'window size for savitzky_golay filtering of data z, must be odd')
            #            ],
            #           Parameter('polynomial_order'),
            #           [Parameter('x', 3, int, 'polynomial order for savitzky_golay filtering of data x, must be odd'),
            #            Parameter('y', 3, int, 'polynomial order for savitzky_golay filtering of data x, must be odd'),
            #            Parameter('z', 3, int, 'polynomial order for savitzky_golay filtering of data x, must be odd')
            #            ]]
            #           ),
            Parameter('smoothing_window_size',
                      [Parameter('x', 21, int, 'window size for savitzky_golay filtering of data x, must be odd and smaller than data size'),
                       Parameter('y', 21, int, 'window size for savitzky_golay filtering of data y, must be odd and smaller than data size'),
                       Parameter('z', 5, int, 'window size for savitzky_golay filtering of data z, must be odd and smaller than data size')
                       ]),
            Parameter('smoothing_polynomial_order',
                      [Parameter('x', 3, int, 'polynomial order for savitzky_golay filtering of data x, must be smaller than window size'),
                       Parameter('y', 3, int, 'polynomial order for savitzky_golay filtering of data y, must be smaller than window size'),
                       Parameter('z', 3, int, 'polynomial order for savitzky_golay filtering of data z, must be smaller than window size')
                       ]),

            Parameter('time_per_pt',
                      [Parameter('galvo', .015, [.0005, .001, .002, .005, .01, .015, .02, 0.05, 0.1, 0.2],
                                 'time in s to measure at each point'),
                       Parameter('z-piezo', .1, [.1, .25, .5, 1.], 'time in s to measure at each point for 1D z-scans only'),
                       ]),
            Parameter('settle_time',
                      [Parameter('galvo', .0005, [.0002, .0005, .001], 'wait time between points to allow galvo to settle'),
                       # Parameter('galvo', .0005, [.0005, .001, .002, .005, 0.01, 0.02, 0.05, 0.1],
                       #           'wait time between points to allow galvo to settle'),
                       Parameter('z-piezo', .25, [.25],
                                 'settle time for objective z-motion (measured for oil objective to be ~10ms, in reality appears to be much longer)'),
                       ]),
            Parameter('DAQ_channels',
                      [Parameter('x_ao_channel', 'ao0', ['ao0', 'ao1', 'ao2', 'ao3'],
                                 'Daq channel used for x voltage analog output'),
                       Parameter('y_ao_channel', 'ao1', ['ao0', 'ao1', 'ao2', 'ao3'],
                                 'Daq channel used for y voltage analog output'),
                       Parameter('z_ao_channel', 'ao2', ['ao0', 'ao1', 'ao2', 'ao3'],
                                 'Daq channel used for z voltage analog output'),
                       Parameter('counter_channel', 'ctr0', ['ctr0', 'ctr1', 'ctr2', 'ctr3'],
                                 'Daq channel used for counter')
                       ]),
        ]

    _INSTRUMENTS = {'PB': LISE607RTPulseBlaster, 'daq': NI6353}
    _SCRIPTS = {'set_focus': SetConfocal, '1d_scan': ConfocalScan}

    def __init__(self, scripts, name = None, settings = None, instruments=None, log_function = None, timeout = 1000000000, data_path = None):

        Script.__init__(self, name, scripts = scripts, settings=settings, instruments = instruments, log_function=log_function, data_path = data_path)
        self.daq_in = self.instruments['daq']['instance']
        self.daq_out = self.instruments['daq']['instance']

    def _function(self):

        # self.data = {'final_point': [],
        #                      'initial_point': [],
        #                      'max_fluor': [],
        #                      'fluor_vector': [[],[],[]],
        #                      'extent': [],
        #
        #                      }
        # turn on laser and apd_switch
        print('turn on laser and APD readout channel.')
        self.instruments['PB']['instance'].update({'laser': {'status': True}})
        self.instruments['PB']['instance'].update({'apd_switch': {'status': True}})

        initial_position = self.instruments['daq']['instance'].get_analog_voltages(
            [self.settings['DAQ_channels']['x_ao_channel'], self.settings['DAQ_channels']['y_ao_channel'],
             self.settings['DAQ_channels']['z_ao_channel']])

        self.data['initial_position'] = initial_position
        print('initial positions are:')
        print(initial_position)

        self.scripts['set_focus'].settings['point'].update({'x': initial_position[0]})
        self.scripts['set_focus'].settings['point'].update({'y': initial_position[1]})
        self.scripts['set_focus'].settings['point'].update({'z': initial_position[2]})

        self.scripts['1d_scan'].settings['point_a'].update(
            {'x': initial_position[0], 'y': initial_position[1], 'z': initial_position[2]})
        self.scripts['1d_scan'].settings['point_b'].update(
            {'x': self.settings['sweep_range']['x'], 'y': self.settings['sweep_range']['y'],
             'z': self.settings['sweep_range']['z']})
        self.scripts['1d_scan'].update({'RoI_mode': 'center'})
        self.scripts['1d_scan'].settings['num_points'].update(
            {'x': self.settings['num_points']['x'], 'y': self.settings['num_points']['y'],
             'z': self.settings['num_points']['z']})
        self.scripts['1d_scan'].settings['time_per_pt'].update(
            {'galvo': self.settings['time_per_pt']['galvo'], 'z-piezo': self.settings['time_per_pt']['z-piezo']})
        self.scripts['1d_scan'].settings['settle_time'].update(
            {'galvo': self.settings['settle_time']['galvo'], 'z-piezo': self.settings['settle_time']['z-piezo']})
        self.scripts['1d_scan'].settings['DAQ_channels'].update(
            {'x_ao_channel': self.settings['DAQ_channels']['x_ao_channel'],
             'y_ao_channel': self.settings['DAQ_channels']['y_ao_channel'],
             'z_ao_channel': self.settings['DAQ_channels']['z_ao_channel'],
             'counter_channel': self.settings['DAQ_channels']['counter_channel']})
        self.data['flag'] = 'flag'

        if self.settings['optimizing_x']:
            self.data['flag'] = 'x'
            print('now optimizing x...')
            self.scripts['1d_scan'].update({'scan_axes': 'x'})
            self.scripts['1d_scan'].run()
            self.data['original_image_data_x'] = deepcopy(self.scripts['1d_scan'].data['image_data'])
            # self.data['image_data_x'] = self.data['original_image_data_x']
            self.data['image_data_x'] = savgol_filter(self.data['original_image_data_x'],
                                                      self.settings['smoothing_window_size']['x'],
                                                      self.settings['smoothing_polynomial_order']['x'])
            self.data['extent_x'] = deepcopy(self.scripts['1d_scan'].data['bounds'])
            self.data['max_fluor_x'] = np.amax(self.data['image_data_x'])
            self.data['maximum_point_x'] = self.data['extent_x'][0] + (self.data['extent_x'][1] - self.data['extent_x'][0]) / (
                        len(self.data['image_data_x']) - 1) * float(np.argmax(self.data['image_data_x']))
            self.log('optimal x is Vx = {:.4}'.format(self.data['maximum_point_x']))
            print('optimal x is Vx = {:.4}'.format(self.data['maximum_point_x']))
            self.scripts['set_focus'].settings['point'].update({'x': self.data['maximum_point_x']})
            self.scripts['set_focus'].run()
        if self.settings['optimizing_y']:
            self.data['flag'] = 'y'
            print('now optimizing y...')
            self.scripts['1d_scan'].update({'scan_axes': 'y'})
            self.scripts['1d_scan'].run()
            self.data['original_image_data_y'] = deepcopy(self.scripts['1d_scan'].data['image_data'])
            # self.data['image_data_y'] = self.data['original_image_data_y']
            self.data['image_data_y'] = savgol_filter(self.data['original_image_data_y'],
                                                      self.settings['smoothing_window_size']['y'],
                                                      self.settings['smoothing_polynomial_order']['y'])
            self.data['extent_y'] = deepcopy(self.scripts['1d_scan'].data['bounds'])
            self.data['max_fluor_y'] = np.amax(self.data['image_data_y'])
            self.data['maximum_point_y'] = self.data['extent_y'][0] + (self.data['extent_y'][1] - self.data['extent_y'][0]) / (
                        len(self.data['image_data_y']) - 1) * float(np.argmax(self.data['image_data_y']))
            self.log('optimal y is Vy = {:.4}'.format(self.data['maximum_point_y']))
            print('optimal y is Vy = {:.4}'.format(self.data['maximum_point_y']))
            self.scripts['set_focus'].settings['point'].update({'y': self.data['maximum_point_y']})
            self.scripts['set_focus'].run()
        if self.settings['optimizing_z']:
            self.data['flag'] = 'z'
            self.scripts['1d_scan'].update({'scan_axes': 'z'})
            self.scripts['1d_scan'].run()
            self.data['original_image_data_z'] = deepcopy(self.scripts['1d_scan'].data['image_data'])
            # self.data['image_data_z'] = self.data['original_image_data_z']
            self.data['image_data_z'] = savgol_filter(self.data['original_image_data_z'],
                                                      self.settings['smoothing_window_size']['z'],
                                                      self.settings['smoothing_polynomial_order']['z'])
            self.data['extent_z'] = deepcopy(self.scripts['1d_scan'].data['bounds'])
            self.data['max_fluor_z'] = np.amax(self.data['image_data_z'])
            self.data['maximum_point_z'] = self.data['extent_z'][0] + (self.data['extent_z'][1] - self.data['extent_z'][0]) / (
                        len(self.data['image_data_z']) - 1) * float(np.argmax(self.data['image_data_z']))
            self.log('optimal z is Vz = {:.4}'.format(self.data['maximum_point_z']))
            print('optimal z is Vz = {:.4}'.format(self.data['maximum_point_z']))
            self.scripts['set_focus'].settings['point'].update({'z': self.data['maximum_point_z']})
            self.scripts['set_focus'].run()


        final_position = self.instruments['daq']['instance'].get_analog_voltages(
            [self.settings['DAQ_channels']['x_ao_channel'], self.settings['DAQ_channels']['y_ao_channel'],
             self.settings['DAQ_channels']['z_ao_channel']])

        self.data['final_position'] = final_position
        print('current positions are Vx={:.4}, Vy={:.4}, Vz={:.4}'.format(final_position[0], final_position[1],                                                          final_position[2]))
        # self.log('current positions are Vx={:.4}, Vy={:.4}, Vz={:.4}'.format(final_position[0], final_position[1],
        #                                                                final_position[2]))
        # turn off laser and apd_switch
        print('turn off laser and APD readout channel.')
        self.instruments['PB']['instance'].update({'laser': {'status': False}})
        self.instruments['PB']['instance'].update({'apd_switch': {'status': False}})

    def plot_data(self, axes_list, data):

        if self.data['flag'] == 'x':
            data['original_fluor_vector'] = self.data['original_image_data_x']
            data['fluor_vector'] = self.data['image_data_x']
            data['extent'] = self.data['extent_x']
            data['maximum_point'] = self.data['maximum_point_x']
            data['max_fluor'] = self.data['max_fluor_x']
            # plot_counts(axes_list[0], data['fluor_vector'],
            #             np.linspace(data['extent'][0], data['extent'][1], len(data['fluor_vector'])), 'x [V]')
            plot_counts_vs_pos(axes_list[0], data['original_fluor_vector'],
                               np.linspace(data['extent'][0], data['extent'][1], len(data['fluor_vector'])),
                               x_label='x [V]')
            plot_counts_vs_pos(axes_list[0], data['fluor_vector'],
                               np.linspace(data['extent'][0], data['extent'][1], len(data['fluor_vector'])),
                               x_label='x [V]')
            if data['maximum_point'] and data['maximum_point_x']:
                axes_list[0].plot(data['maximum_point'], data['max_fluor'], 'ro')
                axes_list[0].text(data['maximum_point'], data['max_fluor'] - .002, 'optimal x', color='r', fontsize=8)
                print('Optimal x: {:0.2f} V, Max fluor: {:0.2f} kcps'.format(data['maximum_point'], data['max_fluor']))

        elif self.data['flag'] == 'y':
            data['original_fluor_vector'] = self.data['original_image_data_y']
            data['fluor_vector'] = self.data['image_data_y']
            data['extent'] = self.data['extent_y']
            data['maximum_point'] = self.data['maximum_point_y']
            data['max_fluor'] = self.data['max_fluor_y']
            # plot_counts(axes_list[0], data['fluor_vector'],
            #             np.linspace(data['extent'][0], data['extent'][1], len(data['fluor_vector'])), 'y [V]')
            plot_counts_vs_pos(axes_list[0], data['original_fluor_vector'],
                               np.linspace(data['extent'][0], data['extent'][1], len(data['fluor_vector'])),
                               x_label='y [V]')
            plot_counts_vs_pos(axes_list[0], data['fluor_vector'],
                               np.linspace(data['extent'][0], data['extent'][1], len(data['fluor_vector'])),
                               x_label='y [V]')
            if data['maximum_point'] and data['max_fluor']:
                # axes_list[0].hold(True)
                axes_list[0].plot(data['maximum_point'], data['max_fluor'], 'go')
                axes_list[0].text(data['maximum_point'], data['max_fluor'] - .002, 'optimal y', color='r', fontsize=8)
                print('Optimal y: {:0.2f} V, Max fluor: {:0.2f} kcps'.format(data['maximum_point'], data['max_fluor']))
        elif self.data['flag'] == 'z':

            data['original_fluor_vector'] = self.data['original_image_data_z']
            data['fluor_vector'] = self.data['image_data_z']
            data['extent'] = self.data['extent_z']
            data['maximum_point'] = self.data['maximum_point_z']
            data['max_fluor'] = self.data['max_fluor_z']

            # plot_counts(axes_list[0], data['fluor_vector'],
            #             np.linspace(data['extent'][0], data['extent'][1], len(data['fluor_vector'])), 'z [V]')
            plot_counts_vs_pos(axes_list[0], data['original_fluor_vector'],
                               np.linspace(data['extent'][0], data['extent'][1], len(data['fluor_vector'])),
                               x_label='z [V]')
            plot_counts_vs_pos(axes_list[0], data['fluor_vector'],
                               np.linspace(data['extent'][0], data['extent'][1], len(data['fluor_vector'])),
                               x_label='z [V]')


            if data['maximum_point'] and data['max_fluor']:
                # axes_list[0].hold(True)
                axes_list[0].plot(data['maximum_point'], data['max_fluor'], 'bo')
                axes_list[0].text(data['maximum_point'], data['max_fluor'] - .002, 'optimal z', color='r', fontsize=8)
                print('Optimal z: {:0.2f} V, Max fluor: {:0.2f} kcps'.format(data['maximum_point'], data['max_fluor']))
                # axes_list[0].hold(False)

    def _plot(self, axes_list, data=None):
        """
        Plots the confocal scan image
        Args:
            axes_list: list of axes objects on which to plot the galvo scan on the first axes object
            data: data (dictionary that contains keys image_data, extent) if not provided use self.data
        """

        if data is None:
            data = self.data

        if self._current_subscript_stage['current_subscript'] == self.scripts['1d_scan']:

            self.scripts['1d_scan']._plot(axes_list)

        else:

            self.plot_data(axes_list, data)


        # plot_counts(axes_list[0], data['fluor_vector'],np.linspace(data['extent'][0],data['extent'][1],len(data['fluor_vector'])),'z [V]')
        # axes_list[0].hold(True)
        # axes_list[0].plot(data['maximum_point'], data['max_fluor'], 'ro')
        # axes_list[0].hold(False)


    def _update_plot(self, axes_list):
        """
        updates the galvo scan image
        Args:
            axes_list: list of axes objects on which to plot plots the esr on the first axes object
        """

        if self._current_subscript_stage['current_subscript'] == self.scripts['1d_scan']:
            self.scripts['1d_scan']._update_plot(axes_list)

        # if self.data['maximum_point'] and self.data['max_fluor']:
        #     axes_list[0].hold(True)
        #     axes_list[0].plot(self.data['maximum_point'], self.data['max_fluor'], 'ro')
        #     axes_list[0].hold(False)


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

        # only pick the first figure from the figure list, this avoids that get_axes_layout clears all the figures
        return super(optimize, self).get_axes_layout([figure_list[1]])


if __name__ == '__main__':
    script, failed, instruments = Script.load_and_append(script_dict={'optimize': 'optimize'})

    print(script)
    print(failed)
    # print(instruments)

