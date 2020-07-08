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
import time

# from b26_toolkit.src.instruments import NI6259, NI9263, NI9402
from b26_toolkit.instruments import NI6353
from b26_toolkit.plotting.plots_2d import plot_fluorescence_new, update_fluorescence
from pylabcontrol.core import Script, Parameter
from b26_toolkit.plotting.plots_1d import plot_counts, update_counts_vs_pos, plot_counts_vs_pos
from b26_toolkit.instruments import LISE607RTPulseBlaster


class ConfocalScan(Script):
    """
            confocal scan x, y and z
            updated by ZQ 1/4/2019 5:44 pm

    """

    _DEFAULT_SETTINGS = [
        Parameter('scan_axes', 'xy', ['xy','xz','yz','x','y','z'],'Choose 2D or 1D confocal scan to perform'),
        Parameter('point_a',
                  [Parameter('x', 0, float, 'x-coordinate [V]'),
                   Parameter('y', 0, float, 'y-coordinate [V]'),
                   Parameter('z', 5, float, 'z-coordinate [V]')
                   ]),
        Parameter('point_b',
                  [Parameter('x', 1.0, float, 'x-coordinate [V]'),
                   Parameter('y', 1.0, float, 'y-coordinate [V]'),
                   Parameter('z', 10.0, float, 'z-coordinate [V]')
                   ]),
        Parameter('RoI_mode', 'center', ['corner', 'center'], 'mode to calculate region of interest.\n \
                                                           corner: pta and ptb are diagonal corners of rectangle.\n \
                                                           center: pta is center and pta is extend or rectangle'),
        Parameter('num_points',
                  [Parameter('x', 125, int, 'number of x points to scan'),
                   Parameter('y', 125, int, 'number of y points to scan'),
                   Parameter('z', 51, int, 'number of z points to scan')
                   ]),
        Parameter('time_per_pt',
                  [Parameter('galvo', .002, [.0005, .001, .002, .005, .01, .015, .02, 0.05, 0.1, 0.2], 'time in s to measure at each point'),
                   Parameter('z-piezo', .5, [.25, .5, 1.], 'time in s to measure at each point for 1D z-scans only'),
                   ]),
        Parameter('settle_time',
                  [Parameter('galvo', .0005, [.0002,.0005, .001], 'wait time between points to allow galvo to settle'),
                   # Parameter('galvo', .0005, [.0005, .001, .002, .005, 0.01, 0.02, 0.05, 0.1],
                   #           'wait time between points to allow galvo to settle'),
                   Parameter('z-piezo', .25, [.25], 'settle time for objective z-motion (measured for oil objective to be ~10ms, in reality appears to be much longer)'),
                   ]),
        # Parameter('min_counts_plot', -1, int, 'Rescales colorbar with this as the minimum counts on replotting'),
        Parameter('max_counts_plot', -1, int, 'Rescales colorbar with this as the maximum counts on replotting'),
        Parameter('min_counts_plot', -1, int, 'Rescales colorbar with this as the maximum counts on replotting'),
        Parameter('DAQ_channels',
                   [Parameter('x_ao_channel', 'ao0', ['ao0', 'ao1', 'ao2', 'ao3'], 'Daq channel used for x voltage analog output'),
                    Parameter('y_ao_channel', 'ao1', ['ao0', 'ao1', 'ao2', 'ao3'], 'Daq channel used for y voltage analog output'),
                    Parameter('z_ao_channel', 'ao2', ['ao0', 'ao1', 'ao2', 'ao3'], 'Daq channel used for z voltage analog output'),
                    Parameter('counter_channel', 'ctr0', ['ctr0', 'ctr1', 'ctr2', 'ctr3'], 'Daq channel used for counter')
                  ]),
        # Parameter('ending_behavior', 'return_to_start', ['return_to_start', 'return_to_origin', 'leave_at_corner'], 'return to the corn'),
        # Parameter('daq_type', 'PCI', ['PCI', 'cDAQ'], 'Type of daq to use for scan')
    ]

    # _INSTRUMENTS = {'NI6259':  NI6259, 'NI9263': NI9263, 'NI9402': NI9402}
    _INSTRUMENTS = {'daq':  NI6353, 'PB': LISE607RTPulseBlaster}
    _SCRIPTS = {}

    def __init__(self, instruments, name=None, settings=None, log_function=None, data_path=None):
        '''
        Initializes ConfocalScan script for use in gui

        Args:
            instruments: list of instrument objects
            name: name to give to instantiated script object
            settings: dictionary of new settings to pass in to override defaults
            log_function: log function passed from the gui to direct log calls to the gui log
            data_path: path to save data

        '''
        Script.__init__(self, name, settings=settings, instruments=instruments, log_function=log_function,
                        data_path=data_path)
        # defines which daqs contain the input and output based on user selection of daq interface
        # if self.settings['daq_type'] == 'PCI':
        self.daq_in = self.instruments['daq']['instance']
        self.daq_out = self.instruments['daq']['instance']
        # elif self.settings['daq_type'] == 'cDAQ':
        #     self.daq_in = self.instruments['NI9402']['instance']
        #     self.daq_out = self.instruments['NI9263']['instance']

    def _function(self):
        """
        Executes threaded galvo scan
        """

        # update_time = datetime.datetime.now()

        # self._plot_refresh = True

        # self._plotting = True

        # turn on laser and apd_switch
        print('turn on laser and APD readout channel.')
        self.instruments['PB']['instance'].update({'laser': {'status': True}})
        self.instruments['PB']['instance'].update({'apd_switch': {'status': True}})

        def scan2D():
            self._recording = False

            self.clockAdjust = int(
                (self.settings['time_per_pt']['galvo'] + self.settings['settle_time']['galvo']) / self.settings['settle_time']['galvo'])

            self.var1_array = np.repeat(
                np.linspace(self.var1range[0], self.var1range[1], self.settings['num_points'][self.var1], endpoint=True),
                self.clockAdjust)
            self.var2_array = np.linspace(self.var2range[0], self.var2range[1], self.settings['num_points'][self.var2], endpoint=True)

            self.data = {'image_data': np.zeros((self.settings['num_points'][self.var2], self.settings['num_points'][self.var1])),
                         'bounds': [self.var1range[0], self.var1range[1], self.var2range[0], self.var2range[1]],
                         'varinitialpos': self.varinitialpos}
            self.data['extent'] = [self.var1range[0], self.var1range[1], self.var2range[1], self.var2range[0]]
            # self.data['varcalib'] = [self.settings['um_per_V'][self.var1],self.settings['um_per_V'][self.var1]]
            # self.data['varlbls'] = [self.var1 + ' [$\mu$m]',self.var2 + ' [$\mu$m]']
            self.data['varlbls'] = [self.var1 + ' [V]', self.var2 + ' [V]']

            # objective takes longer to settle after big jump, so give it time before starting scan:
            if self.settings['scan_axes'] == 'xz' or self.settings['scan_axes'] == 'yz':
                self.daq_out.set_analog_voltages(
                    {self.settings['DAQ_channels'][self.var1channel]: self.var1_array[0],
                     self.settings['DAQ_channels'][self.var2channel]: self.var2_array[0]})
                time.sleep(1)

            for var2Num in range(0, len(self.var2_array)):

                if self._abort:
                    break

                # set galvo to initial point of next line
                self.initPt = [self.var1_array[0], self.var2_array[var2Num]]
                self.daq_out.set_analog_voltages(
                    {self.settings['DAQ_channels'][self.var1channel]: self.initPt[0],
                     self.settings['DAQ_channels'][self.var2channel]: self.initPt[1]})

                # initialize APD thread
                ctrtask = self.daq_in.setup_counter(
                    self.settings['DAQ_channels']['counter_channel'],
                    len(self.var1_array) + 1)
                aotask = self.daq_out.setup_AO([self.settings['DAQ_channels'][self.var1channel]],
                                               self.var1_array, ctrtask)

                # start counter and scanning sequence
                self.daq_out.run(aotask)
                self.daq_in.run(ctrtask)
                self.daq_out.waitToFinish(aotask)
                self.daq_out.stop(aotask)
                var1LineData, _ = self.daq_in.read(ctrtask)
                self.daq_in.stop(ctrtask)
                diffData = np.diff(var1LineData)
                summedData = np.zeros(int(len(self.var1_array) / self.clockAdjust))

                for i in range(0, int((len(self.var1_array) / self.clockAdjust))):
                    pxarray = diffData[(i * self.clockAdjust + 1):(i * self.clockAdjust + self.clockAdjust - 1)]
                    normalization = len(pxarray) / self.sample_rate / 0.001
                    summedData[i] = np.sum(pxarray)/normalization

                #summedData = np.flipud(summedData)

                # also normalizing to kcounts/sec
                # self.data['image_data'][var2Num] = summedData * (.001 / self.settings['time_per_pt']['galvo'])
                self.data['image_data'][var2Num] = summedData
                self.progress = float(var2Num + 1) / len(self.var2_array) * 100
                self.updateProgress.emit(int(self.progress))

        def scan1D():
            self._recording = False
            self.clockAdjust = int(
                (self.settings['time_per_pt']['galvo'] + self.settings['settle_time']['galvo']) / self.settings['settle_time']['galvo'])

            self.var1_array = np.repeat(
                np.linspace(self.var1range[0], self.var1range[1], self.settings['num_points'][self.var1],endpoint=True),self.clockAdjust)
            # self.var2_array = np.linspace(self.var2range[0], self.var2range[1], self.settings['num_points'][self.var2], endpoint=True)
            self.data = {'image_data': np.zeros(self.settings['num_points'][self.var1]),
                         'bounds': [self.var1range[0], self.var1range[1]],
                         'varinitialpos': self.varinitialpos}
            self.data['varlbls'] = self.var1 + ' [V]'
            # print('right before while loop')

            # while True:
            for x in range(1):
                if self._abort:
                    break
                # print('set galvo to initial point of 1D scan')
                # set galvo to initial point of 1D scan
                self.initPt = self.var1_array[0]
                self.daq_out.set_analog_voltages({self.settings['DAQ_channels'][self.var1channel]: self.initPt})

                # print('initialize APD thread')
                # initialize APD thread
                ctrtask = self.daq_in.setup_counter(
                    self.settings['DAQ_channels']['counter_channel'],len(self.var1_array) + 1)
                aotask = self.daq_out.setup_AO([self.settings['DAQ_channels'][self.var1channel]],self.var1_array, ctrtask)

                # start counter and scanning sequence
                # print('start counter and scanning sequence')
                self.daq_out.run(aotask)
                self.daq_in.run(ctrtask)
                self.daq_out.waitToFinish(aotask)
                self.daq_out.stop(aotask)

                var1LineData, _ = self.daq_in.read(ctrtask)
                self.daq_in.stop(ctrtask)
                diffData = np.diff(var1LineData)
                # print('summedData = np.zeros(len(self.var1_array) / self.clockAdjust)')
                # print('self.clockAdjust)=')
                # print(self.clockAdjust)
                # print('len(self.var1_array)=')
                # print(len(self.var1_array))

                summedData = np.zeros(int(len(self.var1_array) / self.clockAdjust))
                # print('start iteration')
                for i in range(0, int((len(self.var1_array) / self.clockAdjust))):
                    # print('i=')
                    # print(i)
                    pxarray = diffData[(i * self.clockAdjust + 1):(i * self.clockAdjust + self.clockAdjust - 1)]
                    normalization = len(pxarray) / self.sample_rate / 0.001
                    summedData[i] = np.sum(pxarray)/normalization
                # also normalizing to kcounts/sec
                # self.data['image_data'] = summedData * (.001 / self.settings['time_per_pt']['galvo'])

                # print('iteration done')
                self.data['image_data'] = summedData

                self.progress = 50.
                self.updateProgress.emit(int(self.progress))

        def scanZ():
            self._recording = False

            nsamples = int((self.settings['time_per_pt']['z-piezo'] + self.settings['settle_time']['z-piezo']) * self.sample_rate)
            self.var1_array = np.linspace(self.var1range[0], self.var1range[1], self.settings['num_points'][self.var1], endpoint=True)

            self.data = {'image_data': np.zeros((self.settings['num_points'][self.var1])),
                         'bounds': [self.var1range[0], self.var1range[1]],
                         'varinitialpos': self.varinitialpos}
            self.data['varlbls'] = self.var1 + ' [V]'

            # objective takes longer to settle after a big jump, so give it time before starting scan:
            self.daq_out.set_analog_voltages({self.settings['DAQ_channels'][self.var1channel]: self.var1_array[0]})
            time.sleep(1)

            for var1Num in range(0, len(self.var1_array)):
                if self._abort:
                    break
                # initialize APD thread
                ctrtask = self.daq_in.setup_counter(self.settings['DAQ_channels']['counter_channel'],nsamples)
                aotask = self.daq_out.setup_AO([self.settings['DAQ_channels'][self.var1channel]], self.var1_array[var1Num]*np.ones(nsamples), ctrtask)

                # start counter and scanning sequence
                self.daq_out.run(aotask)
                self.daq_in.run(ctrtask)
                self.daq_out.waitToFinish(aotask)
                self.daq_out.stop(aotask)
                samparray, _ = self.daq_in.read(ctrtask)
                self.daq_in.stop(ctrtask)
                diffData = np.diff(samparray)

                # sum and normalize to kcounts/sec
                # self.data['image_data'][var1Num] = np.sum(diffData) * (.001 / self.settings['time_per_pt']['z-piezo'])
                normalization = len(diffData) / self.sample_rate / 0.001
                self.data['image_data'][var1Num] = np.sum(diffData) / normalization

                self.progress = float(var1Num + 1) / len(self.var1_array) * 100
                self.updateProgress.emit(int(self.progress))

        # if self.settings['daq_type'] == 'PCI':
        initial_position = self.daq_out.get_analog_voltages(
            [self.settings['DAQ_channels']['x_ao_channel'], self.settings['DAQ_channels']['y_ao_channel'], self.settings['DAQ_channels']['z_ao_channel']])

        print('initial positions are:')
        print(initial_position)
        [self.xVmin, self.xVmax, self.yVmin, self.yVmax, self.zVmin, self.zVmax] = self.pts_to_extent(
            self.settings['point_a'],
            self.settings['point_b'],
            self.settings['RoI_mode'])

        self.sample_rate = float(1) / self.settings['settle_time']['galvo']

        self.daq_out.settings['analog_output'][
            self.settings['DAQ_channels']['x_ao_channel']]['sample_rate'] = self.sample_rate
        self.daq_out.settings['analog_output'][
            self.settings['DAQ_channels']['y_ao_channel']]['sample_rate'] = self.sample_rate
        self.daq_out.settings['analog_output'][
            self.settings['DAQ_channels']['z_ao_channel']]['sample_rate'] = self.sample_rate
        self.daq_in.settings['digital_input'][
            self.settings['DAQ_channels']['counter_channel']]['sample_rate'] = self.sample_rate

        # print('ready for scanning')

        # depending to axes to be scanned, assigns the correct channels to be scanned and scan ranges, then starts the 2D or 1D scan
        self.var1range = 0
        self.var2range = 0
        if self.settings['scan_axes'] == 'xy':
            self.var1 = 'x'
            self.var2 = 'y'
            self.var1channel = 'x_ao_channel'
            self.var2channel = 'y_ao_channel'
            self.var1range = [self.xVmin,self.xVmax]
            self.var2range = [self.yVmin,self.yVmax]
            self.varinitialpos = [initial_position[0],initial_position[1]]
            scan2D()
        elif self.settings['scan_axes'] == 'xz':
            self.var1 = 'x'
            self.var2 = 'z'
            self.var1channel = 'x_ao_channel'
            self.var2channel = 'z_ao_channel'
            self.var1range = [self.xVmin,self.xVmax]
            self.var2range = [self.zVmin,self.zVmax]
            self.varinitialpos = [initial_position[0],initial_position[2]]
            scan2D()
        elif self.settings['scan_axes'] == 'yz':
            self.var1 = 'y'
            self.var2 = 'z'
            self.var1channel = 'y_ao_channel'
            self.var2channel = 'z_ao_channel'
            self.var1range = [self.yVmin,self.yVmax]
            self.var2range = [self.zVmin,self.zVmax]
            self.varinitialpos = [initial_position[1],initial_position[2]]
            scan2D()
        elif self.settings['scan_axes'] == 'x':
            self.var1 = 'x'
            self.var1channel = 'x_ao_channel'
            self.var1range = [self.xVmin,self.xVmax]
            self.varinitialpos = initial_position[0]
            scan1D()
        elif self.settings['scan_axes'] == 'y':
            self.var1 = 'y'
            self.var1channel = 'y_ao_channel'
            self.var1range = [self.yVmin,self.yVmax]
            self.varinitialpos = initial_position[1]
            scan1D()
        elif self.settings['scan_axes'] == 'z':
            self.var1 = 'z'
            self.var1channel = 'z_ao_channel'
            self.var1range = [self.zVmin,self.zVmax]
            self.varinitialpos = initial_position[2]
            scanZ()


        self.daq_out.set_analog_voltages(
            {self.settings['DAQ_channels']['x_ao_channel']: initial_position[0],
            self.settings['DAQ_channels']['y_ao_channel']: initial_position[1],
            self.settings['DAQ_channels']['z_ao_channel']: initial_position[2]})
        print('voltage returned to initial values')

        # turn off laser and apd_switch
        self.instruments['PB']['instance'].update({'laser': {'status': False}})
        self.instruments['PB']['instance'].update({'apd_switch': {'status': False}})
        print('laser and APD readout is off.')


    def get_confocal_location(self):
        """
        Returns the current position of the galvo. Requires a daq with analog inputs internally routed to the analog
        outputs (ex. NI6353. Note that the cDAQ does not have this capability).
        Returns: list with two floats, which give the x and y position of the galvo mirror
        """
        confocal_position = self.daq_out.get_analog_voltages([
            self.settings['DAQ_channels']['x_ao_channel'],
            self.settings['DAQ_channels']['y_ao_channel'],
            self.settings['DAQ_channels']['z_ao_channel']]
        )
        return confocal_position

    def set_confocal_location(self, confocal_position):
        """
        sets the current position of the confocal
        confocal_position: list with three floats, which give the x, y, z positions of the confocal (galvo mirrors and objective)
        """
        print('\t'.join(map(str, confocal_position)))
        if confocal_position[0] > 10 or confocal_position[0] < -10 or confocal_position[1] > 10 or confocal_position[1] < -10:
            raise ValueError('The script attempted to set the galvo position to an illegal position outside of +- 10 V')
        if confocal_position[2] > 10 or confocal_position[2] < 0:
            raise ValueError('The script attempted to set the objective position to an illegal position outside of 0-10 V')

        pt = confocal_position
        # daq API only accepts either one point and one channel or multiple points and multiple channels
        pt = np.transpose(np.column_stack((pt[0], pt[1], pt[2])))
        pt = (np.repeat(pt, 3, axis=1))

        task = self.daq_out.setup_AO(
            [self.settings['DAQ_channels']['x_ao_channel'], self.settings['DAQ_channels']['y_ao_channel'], self.settings['DAQ_channels']['z_ao_channel']], pt)
        self.daq_out.run(task)
        self.daq_out.waitToFinish(task)
        self.daq_out.stop(task)

    @staticmethod
    def pts_to_extent(pta, ptb, roi_mode):
        """

        Args:
            pta: point a
            ptb: point b
            roi_mode:   mode how to calculate region of interest
                        corner: pta and ptb are diagonal corners of rectangle.
                        center: pta is center and ptb is extend or rectangle

        Returns: extend of region of interest [xVmin, xVmax, yVmax, yVmin]

        """
        if roi_mode == 'corner':
            xVmin = min(pta['x'], ptb['x'])
            xVmax = max(pta['x'], ptb['x'])
            yVmin = min(pta['y'], ptb['y'])
            yVmax = max(pta['y'], ptb['y'])
            zVmin = min(pta['z'], ptb['z'])
            zVmax = max(pta['z'], ptb['z'])
        elif roi_mode == 'center':
            xVmin = pta['x'] - float(ptb['x']) / 2.
            xVmax = pta['x'] + float(ptb['x']) / 2.
            yVmin = pta['y'] - float(ptb['y']) / 2.
            yVmax = pta['y'] + float(ptb['y']) / 2.
            zVmin = pta['z'] - float(ptb['z']) / 2.
            zVmax = pta['z'] + float(ptb['z']) / 2.
        return [xVmin, xVmax, yVmin, yVmax, zVmin, zVmax]

    def plot(self, figure_list):
        # Choose whether to plot results in top or bottom figure
        # print('plot')
        if 'image_data' in self.data.keys() is not None:

            if np.ndim(self.data['image_data'])==2:
                super(ConfocalScan, self).plot([figure_list[0]])
            elif np.ndim(self.data['image_data'])==1:
                super(ConfocalScan, self).plot([figure_list[1]])

    def _plot(self, axes_list, data=None):
        """
        Plots the confocal scan image
        Args:
            axes_list: list of axes objects on which to plot the galvo scan on the first axes object
            data: data (dictionary that contains keys image_data, extent) if not provided use self.data
        """
        # print('_plot')
        # print('np.ndim(data')

        if data is None:
            data = self.data
        # plot_fluorescence_new(data['image_data'], data['extent'], self.data['varcalib'], self.data['varlbls'], self.data['varinitialpos'], axes_list[0], min_counts=self.settings['min_counts_plot'], max_counts=self.settings['max_counts_plot'])
        # print(np.ndim(data['image_data']))
        if np.ndim(data['image_data'])==2:
            # plot_fluorescence_new(data['image_data'], data['extent'], self.data['varlbls'], self.data['varinitialpos'], axes_list[0], min_counts=self.settings['min_counts_plot'], max_counts=self.settings['max_counts_plot'])
            plot_fluorescence_new(data['image_data'], data['extent'],axes_list[0], max_counts=self.settings['max_counts_plot'], min_counts=self.settings['min_counts_plot'], axes_labels = self.settings['scan_axes'])
        elif np.ndim(data['image_data'])==1:
            # plot_counts(axes_list[0], data['image_data'], axes_labels=self.data['varlbls'])
            plot_counts_vs_pos(axes_list[0], data['image_data'], np.linspace(data['bounds'][0],data['bounds'][1],len(data['image_data'])), x_label = data['varlbls'])



    def _update_plot(self, axes_list):
        """
        updates the galvo scan image
        Args:
            axes_list: list of axes objects on which to plot plots the esr on the first axes object
        """
        # print('_update_plot')
        if np.ndim(self.data['image_data'])==2:
            # update_fluorescence(self.data['image_data'], axes_list[0], self.settings['min_counts_plot'], self.settings['max_counts_plot'])
            update_fluorescence(self.data['image_data'], axes_list[0], max_counts= self.settings['max_counts_plot'], min_counts=self.settings['min_counts_plot'])
        elif np.ndim(self.data['image_data']) == 1:
            # plot_counts(axes_list[0], self.data['image_data'], np.linspace(self.data['bounds'][0],self.data['bounds'][1],len(self.data['image_data'])), self.data['varlbls'])
            update_counts_vs_pos(axes_list[0], self.data['image_data'], np.linspace(self.data['bounds'][0],self.data['bounds'][1],len(self.data['image_data'])))

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
        return super(ConfocalScan, self).get_axes_layout([figure_list[0]])


if __name__ == '__main__':
    script, failed, instruments = Script.load_and_append(script_dict={'ConfocalScan': 'ConfocalScan'})

    print(script)
    print(failed)
    # print(instruments)

