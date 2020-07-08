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

import time
from collections import deque
import numpy as np

from b26_toolkit.instruments import NI6353
from b26_toolkit.plotting.plots_1d import plot_counts, plot_counts_vs_pos, update_1d_simple, update_counts_vs_pos
from pylabcontrol.core import Parameter, Script
from b26_toolkit.instruments import LISE607RTPulseBlaster


import numpy as np
from matplotlib import patches

from b26_toolkit.instruments import PM100D, IntensityWheel
from pylabcontrol.core import Script, Parameter
from b26_toolkit.scripts.daq_read_counter import Daq_Read_Counter

from collections import deque
import time

class IntensityWheel_Calibration(Script):
    """
       This script reads power meter as a function of intensity wheel
       -- Ziwei Qiu 5/16/2019
    """
    _DEFAULT_SETTINGS = [
        Parameter('to-do', 'move', ['move', 'sweep', 'read'], 'Choose to move to a point, do a sweep or just read the current position'),
        Parameter('move_to', 255.0, float, 'move to the position in deg'),
        Parameter('start', 0.0, float, 'start position of the intensity wheel in deg'),
        Parameter('stop', 360.0, float, 'stop position of the intensity wheel in deg'),
        Parameter('num_of_points', 20, int, 'number of points in sweeping the intensity wheel position'),
        Parameter('gain_factor', 14, float, 'set the gain factor of the power meter'),
        Parameter('ini_stab_time', 8, int, 'laser power stabilization time (sec)'),
    ]

    _INSTRUMENTS = {'power_meter':PM100D, 'intensity_wheel': IntensityWheel,'PB': LISE607RTPulseBlaster}
    _SCRIPTS = {}

    def __init__(self, instruments, scripts=None, name=None, settings=None, log_function=None, data_path=None):
        Script.__init__(self, name, settings=settings, scripts=scripts, instruments=instruments,
                        log_function=log_function, data_path=data_path)

        self.data = None

    def _function(self):

        self.instruments['power_meter']['instance'].update({'gain_factor': self.settings['gain_factor']})

        if self.settings['to-do'] == 'move':
            self.instruments['intensity_wheel']['instance'].update({'angle': self.settings['move_to']})
            current_position = self.instruments['intensity_wheel']['instance']._get_position()
            current_velocity = self.instruments['intensity_wheel']['instance']._get_velocity()
            self.log('IntensityWheel moved to: ' + str(current_position) + 'deg, velocity = ' + str(current_velocity) + 'deg/s')
        elif self.settings['to-do'] == 'read':
            current_position = self.instruments['intensity_wheel']['instance']._get_position()
            current_velocity = self.instruments['intensity_wheel']['instance']._get_velocity()
            self.log('IntensityWheel position: ' + str(current_position) + 'deg, velocity = ' + str(current_velocity) + 'deg/s')

        else:
            original_position = self.instruments['intensity_wheel']['instance']._get_position()

            self.data = {'pos': [], 'power': []}
            pos_values = np.linspace(self.settings['start'], self.settings['stop'], self.settings['num_of_points'])
            power_data = np.zeros(self.settings['num_of_points'])  # for the raw esr data
            self.data.update({'pos': pos_values})

            # turn on laser and (optionally) stabilize the power for a few seconds
            self.instruments['PB']['instance'].update({'laser': {'status': True}})
            if self.settings['ini_stab_time'] > 0:
                self.log('laser power stabilization for ' + str(self.settings['ini_stab_time']) + ' sec')
                time.sleep(self.settings['ini_stab_time'])

            pos_index = 0
            for pos in pos_values:
                if self._abort:
                    break
                self.instruments['intensity_wheel']['instance'].update({'angle': float(pos)})
                power_data[pos_index] = self.instruments['power_meter']['instance'].get_power()
                self.data.update({'power': power_data})
                progress = self._calc_progress(pos_index)
                self.updateProgress.emit(progress)
                pos_index += 1

            self.data['exp_finished'] = 1  # this is used to flag when experiment is done and it allows proper saving to datasets

            # go back to the original position and turn off the laser
            self.instruments['intensity_wheel']['instance'].update({'angle': original_position})
            current_power = self.instruments['power_meter']['instance'].get_power()
            self.log('Intensity wheel moved back to the original position = {:0.2f} deg'.format(original_position))
            self.log('Optical power is {:0.4f} mW'.format(current_power))
            self.instruments['PB']['instance'].update({'laser': {'status': False}})

    def _calc_progress(self, pos_index):
        #COMMENT_ME

        progress = float(pos_index) / self.settings['num_of_points'] * 100.
        self.progress = progress
        return int(progress)

    def _plot(self, axes_list, data = None):
        """
        plotting function for esr
        Args:
            axes_list: list of axes objects on which to plot plots the esr on the first axes object
            data: data (dictionary that contains keys frequency, data and fit_params) if not provided use self.data
        Returns:

        """

        if data is None:
            data = self.data
        if len(data['power']) > 0:

            plot_counts_vs_pos(axes_list[0], data['power'], data['pos'], x_label='Intensity wheel [deg]', y_label='Power [mW]', title = 'Intensity Wheel Calibration')

    def _update_plot(self, axes_list, data = None):
        if data is None:
            data = self.data
        if data:

            update_counts_vs_pos(axes_list[0], data['power'], data['pos'])

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
        return super(IntensityWheel_Calibration, self).get_axes_layout([figure_list[1]])

class Saturation(Script):
    """
        This script measures the saturation curve of NV.
        -- Ziwei Qiu 5/17/2019
    """

    _DEFAULT_SETTINGS = [
        Parameter('start', 270.0, float, 'start position of the intensity wheel in deg'),
        Parameter('stop', 340.0, float, 'stop position of the intensity wheel in deg'),
        Parameter('num_of_points', 20, int, 'number of points in sweeping the intensity wheel position'),
        Parameter('time_per_point', 2.0, float, 'measurement time for each point in sec'),
        Parameter('gain_factor', 14.0, float, 'set the gain factor of the power meter'),
        Parameter('ini_stab_time', 8, int, 'laser power stabilization time (sec)'),

    ]
    _INSTRUMENTS = {'power_meter': PM100D, 'intensity_wheel': IntensityWheel, 'PB': LISE607RTPulseBlaster}
    _SCRIPTS = {'counter': Daq_Read_Counter}

    def __init__(self, instruments, scripts=None, name=None, settings=None, log_function=None, data_path=None):
        Script.__init__(self, name, settings=settings, scripts=scripts, instruments=instruments,
                        log_function=log_function, data_path=data_path)

        self.data = None

    def _function(self):

        self.instruments['power_meter']['instance'].update({'gain_factor': self.settings['gain_factor']})
        original_position = self.instruments['intensity_wheel']['instance']._get_position()
        self.scripts['counter'].update({'total_int_time': self.settings['time_per_point']})
        self.scripts['counter'].update({'integration_time': 0.25})
        self.scripts['counter'].update({'laser_on_before': False})
        self.scripts['counter'].update({'laser_off_after': False})
        self.data = {'pos': [], 'power': [], 'fluor': []}
        pos_values = np.linspace(self.settings['start'], self.settings['stop'], self.settings['num_of_points'])
        power_data = np.zeros(self.settings['num_of_points'])
        fluor_data = np.zeros(self.settings['num_of_points'])
        self.data.update({'power': power_data})
        self.data.update({'fluor': fluor_data})
        self.data.update({'pos': pos_values})

        # turn on laser and (optionally) stabilize the power for a few seconds
        self.instruments['PB']['instance'].update({'laser': {'status': True}})
        if self.settings['ini_stab_time'] > 0:
            self.log('laser power stabilization for ' + str(self.settings['ini_stab_time']) + ' sec')
            time.sleep(self.settings['ini_stab_time'])

        pos_index = 0

        for pos in pos_values:
            if self._abort:
                break
            self.instruments['intensity_wheel']['instance'].update({'angle': float(pos)})
            power_data[pos_index] = self.instruments['power_meter']['instance'].get_power()
            self.data.update({'power': power_data})

            time.sleep(0.1)


            self.scripts['counter'].run()
            fluor_data[pos_index] = np.mean(self.scripts['counter'].data['counts'])
            self.data.update({'fluor': fluor_data})

            print('fluor:', self.data['fluor'])
            print('power:', self.data['power'])

            progress = self._calc_progress(pos_index)
            self.updateProgress.emit(progress)
            pos_index += 1

        # this is used to flag when experiment is done and it allows proper saving to datasets
        self.data['exp_finished'] = 1

        # go back to the original position and turn off the laser
        self.instruments['intensity_wheel']['instance'].update({'angle': original_position})
        current_power = self.instruments['power_meter']['instance'].get_power()
        self.log('Intensity wheel moved back to the original position = {:0.2f} deg'.format(original_position))
        self.log('Optical power is {:0.4f} mW'.format(current_power))
        self.instruments['PB']['instance'].update({'laser': {'status': False}})

    def _calc_progress(self, pos_index):
        #COMMENT_ME

        progress = float(pos_index) / self.settings['num_of_points'] * 100.
        self.progress = progress
        return int(progress)

    def _plot(self, axes_list, data = None):
        """
        plotting function for esr
        Args:
            axes_list: list of axes objects on which to plot plots the esr on the first axes object
            data: data (dictionary that contains keys frequency, data and fit_params) if not provided use self.data
        Returns:

        """

        if data is None:
            data = self.data
            # print('fluor:', data['fluor'])
            # print('power:', data['power'])
        if len(data['power']) > 0 and len(data['fluor']) > 0:

            plot_counts_vs_pos(axes_list[0], data['fluor'], data['power'], x_label='Power [mW]', y_label='kCounts/sec', title = 'Saturation Curve', marker = 'd')

    def _update_plot(self, axes_list, data = None):
        if data is None:
            data = self.data
        if data:

            update_counts_vs_pos(axes_list[0], data['fluor'], data['power'])

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
        return super(Saturation, self).get_axes_layout([figure_list[1]])

if __name__ == '__main__':
    script = {}
    instr = {}
    # IntensityWheel_Calibration
    # print('line 513')
    # script, failed, instr = Script.load_and_append({'Daq_Read_Cntr': 'Daq_Read_Cntr'}, script, instr)
    script, failed, instr = Script.load_and_append({'IntensityWheel_Calibration': 'IntensityWheel_Calibration'}, script, instr)

    print(script)
    print(failed)
    print(instr)
