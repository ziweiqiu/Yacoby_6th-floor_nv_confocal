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
from matplotlib import patches
import matplotlib.collections
import scipy.spatial
import time

# from b26_toolkit.instruments import NI6259, NI9263
from b26_toolkit.instruments import NI6353
from pylabcontrol.core import Script, Parameter

class SetLaser(Script):
    """
This script points the laser to a point
updated by ZQ 1/2/2019 6:09 pm
    """

    _DEFAULT_SETTINGS = [
        Parameter('point',
                  [Parameter('x', 0.0, float, 'x-coordinate'),
                   Parameter('y', 0.0, float, 'y-coordinate')
                   ]),
        Parameter('patch_size', 0.005, [0.0005, 0.005, 0.05, 0.5], 'size of the red circle'),
        Parameter('DAQ_channels',
            [Parameter('x_ao_channel', 'ao0', ['ao0', 'ao1', 'ao2', 'ao3'], 'Daq channel used for x voltage analog output'),
            Parameter('y_ao_channel', 'ao1', ['ao0', 'ao1', 'ao2', 'ao3'], 'Daq channel used for y voltage analog output')
            ]),
        Parameter('daq_type', 'PCI', ['PCI'], 'Type of daq to use for scan')
    ]

    # _INSTRUMENTS = {'NI6259':  NI6259, 'NI9263': NI9263}
    _INSTRUMENTS = {'NI6353': NI6353}

    _SCRIPTS = {}


    def __init__(self, instruments = None, scripts = None, name = None, settings = None, log_function = None, data_path = None):
        """
        Example of a script that emits a QT signal for the gui
        Args:
            name (optional): name of script, if empty same as class name
            settings (optional): settings for this script, if empty same as default settings
        """
        Script.__init__(self, name, settings = settings, instruments = instruments, scripts = scripts, log_function= log_function, data_path = data_path)


    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """
        pt = (self.settings['point']['x'], self.settings['point']['y'])

        # daq API only accepts either one point and one channel or multiple points and multiple channels
        pt = np.transpose(np.column_stack((pt[0],pt[1])))
        pt = (np.repeat(pt, 2, axis=1))
        # print(pt)

        self._setup_daq()

        task = self.daq_out.setup_AO([self.settings['DAQ_channels']['x_ao_channel'], self.settings['DAQ_channels']['y_ao_channel']], pt)

        self.daq_out.run(task)
        self.daq_out.waitToFinish(task)
        self.daq_out.stop(task)
        self.log('laser set to Vx={:.4}, Vy={:.4}'.format(self.settings['point']['x'], self.settings['point']['y']))

        # if self.settings['daq_type'] == 'PCI':
        #     self.daq_out = self.instruments['NI6259']['instance']
        # elif self.settings['daq_type'] == 'cDAQ':
        #     self.daq_out = self.instruments['NI9263']['instance']
        self.daq_out = self.instruments['NI6353']['instance']

    def _setup_daq(self):
        # defines which daqs contain the input and output based on user selection of daq interface
        # if self.settings['daq_type'] == 'PCI':
        #     self.daq_out = self.instruments['NI6259']['instance']
        # elif self.settings['daq_type'] == 'cDAQ':
        #     self.daq_out = self.instruments['NI9263']['instance']
        self.daq_out = self.instruments['NI6353']['instance']

    def get_galvo_position(self):
        """
        reads the current position from the x and y channels and returns it
        Returns:

        """
        if self.settings['daq_type'] == 'PCI':
            galvo_position = self.daq_out.get_analog_voltages([
                self.settings['DAQ_channels']['x_ao_channel'],
                self.settings['DAQ_channels']['y_ao_channel']]
            )
        elif self.settings['daq_type'] == 'cDAQ':
            print("WARNING cDAQ doesn't allow to read values")
            galvo_position = []

        return galvo_position
    #must be passed figure with galvo plot on first axis

    def plot(self, figure_list):
        axes_Image = figure_list[0].axes[0]

        # removes patches
        [child.remove() for child in axes_Image.get_children() if isinstance(child, patches.Circle)]

        patch = patches.Circle((self.settings['point']['x'], self.settings['point']['y']), self.settings['patch_size'], fc='r')
        axes_Image.add_patch(patch)

class SetConfocal(Script):
    """
This script points the laser to a point
updated by ZQ 1/2/2019 6:09 pm
    """

    _DEFAULT_SETTINGS = [
        Parameter('point',
                  [Parameter('x', 0.0, float, 'x-coordinate'),
                   Parameter('y', 0.0, float, 'y-coordinate'),
                   Parameter('z', 5.0, float, 'z-coordinate')
                   ]),
        Parameter('patch_size', 0.005, [0.0005, 0.005, 0.05, 0.5], 'size of the red circle'),
        Parameter('DAQ_channels',
            [Parameter('x_ao_channel', 'ao0', ['ao0', 'ao1', 'ao2', 'ao3'], 'Daq channel used for x voltage analog output'),
            Parameter('y_ao_channel', 'ao1', ['ao0', 'ao1', 'ao2', 'ao3'], 'Daq channel used for y voltage analog output'),
             Parameter('z_ao_channel', 'ao2', ['ao0', 'ao1', 'ao2', 'ao3'],
                       'Daq channel used for z voltage analog output'),
             Parameter('sample_rates',5000,[1000,2000,5000],'set same sample rates for all three ao channels [Hz]')
            ]),
        Parameter('daq_type', 'PCI', ['PCI'], 'Type of daq to use for scan')
    ]

    # _INSTRUMENTS = {'NI6259':  NI6259, 'NI9263': NI9263}
    _INSTRUMENTS = {'NI6353': NI6353}
    _SCRIPTS = {}


    def __init__(self, instruments = None, scripts = None, name = None, settings = None, log_function = None, data_path = None):
        """
        Example of a script that emits a QT signal for the gui
        Args:
            name (optional): name of script, if empty same as class name
            settings (optional): settings for this script, if empty same as default settings
        """
        Script.__init__(self, name, settings = settings, instruments = instruments, scripts = scripts, log_function= log_function, data_path = data_path)


    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """
        pt = (self.settings['point']['x'], self.settings['point']['y'], self.settings['point']['z'])

        # daq API only accepts either one point and one channel or multiple points and multiple channels

        pt = np.transpose(np.column_stack((pt[0], pt[1], pt[2])))
        pt = (np.repeat(pt,3, axis=1))

        # print(pt)

        self._setup_daq()

        # force all channels to have the same sample rates (ZQ 1/7/2019 7:37pm)
        self.daq_out.settings['analog_output']['ao0']['sample_rate'] = self.settings['DAQ_channels']['sample_rates']
        self.daq_out.settings['analog_output']['ao1']['sample_rate'] = self.settings['DAQ_channels']['sample_rates']
        self.daq_out.settings['analog_output']['ao2']['sample_rate'] = self.settings['DAQ_channels']['sample_rates']
        self.daq_out.settings['analog_output']['ao3']['sample_rate'] = self.settings['DAQ_channels']['sample_rates']



        task = self.daq_out.setup_AO(
            [self.settings['DAQ_channels']['x_ao_channel'], self.settings['DAQ_channels']['y_ao_channel'],self.settings['DAQ_channels']['z_ao_channel']], pt)


        self.daq_out.run(task)
        self.daq_out.waitToFinish(task)
        self.daq_out.stop(task)

        confocal_position = self.daq_out.get_analog_voltages(
            [self.settings['DAQ_channels']['x_ao_channel'], self.settings['DAQ_channels']['y_ao_channel'],
             self.settings['DAQ_channels']['z_ao_channel']])


        print('laser is set to Vx={:.4}, Vy={:.4}, Vz={:.4}'.format(confocal_position[0], confocal_position[1], confocal_position[2]))
        self.log('laser is set to Vx={:.4}, Vy={:.4}, Vz={:.4}'.format(confocal_position[0], confocal_position[1],
                                                                 confocal_position[2]))

        # if self.settings['daq_type'] == 'PCI':
        #     self.daq_out = self.instruments['NI6259']['instance']
        # elif self.settings['daq_type'] == 'cDAQ':
        #     self.daq_out = self.instruments['NI9263']['instance']
        self.daq_out = self.instruments['NI6353']['instance']

    def _setup_daq(self):
        # defines which daqs contain the input and output based on user selection of daq interface
        # if self.settings['daq_type'] == 'PCI':
        #     self.daq_out = self.instruments['NI6259']['instance']
        # elif self.settings['daq_type'] == 'cDAQ':
        #     self.daq_out = self.instruments['NI9263']['instance']
        self.daq_out = self.instruments['NI6353']['instance']

    def get_galvo_position(self):
        """
        reads the current position from the x and y channels and returns it
        Returns:

        """
        if self.settings['daq_type'] == 'PCI':
            galvo_position = self.daq_out.get_analog_voltages([
                self.settings['DAQ_channels']['x_ao_channel'],
                self.settings['DAQ_channels']['y_ao_channel'],
                self.settings['DAQ_channels']['z_ao_channel']]
            )
        elif self.settings['daq_type'] == 'cDAQ':
            print("WARNING cDAQ doesn't allow to read values")
            galvo_position = []

        return galvo_position
    #must be passed figure with galvo plot on first axis

    def plot(self, figure_list):
        axes_Image = figure_list[0].axes[0]



        # removes patches
        [child.remove() for child in axes_Image.get_children() if isinstance(child, patches.Circle)]

        patch = patches.Circle((self.settings['point']['x'], self.settings['point']['y']), self.settings['patch_size'], fc='r')
        axes_Image.add_patch(patch)

# class ClickSetLaser(Script):
#     """
#     This script selects a point on an image and points the laser to that point.
#     This script is not working. Mouse click event only works when SelectPoints script is coming from PyLabControl, not from B26toolkit.
#     -- Ziwei Qiu 7/6/2019
#     """
#
#     _DEFAULT_SETTINGS = [
#         # Parameter('point',
#         #           [Parameter('x', 0.0, float, 'x-coordinate'),
#         #            Parameter('y', 0.0, float, 'y-coordinate')
#         #            ]),
#         Parameter('patch_size', 0.005, [0.0005, 0.005, 0.05, 0.5], 'size of the red circle'),
#         Parameter('DAQ_channels',
#             [Parameter('x_ao_channel', 'ao0', ['ao0', 'ao1', 'ao2', 'ao3'], 'Daq channel used for x voltage analog output'),
#             Parameter('y_ao_channel', 'ao1', ['ao0', 'ao1', 'ao2', 'ao3'], 'Daq channel used for y voltage analog output')
#             ]),
#         Parameter('daq_type', 'PCI', ['PCI'], 'Type of daq to use for scan')
#     ]
#
#     _INSTRUMENTS = {'NI6353': NI6353}
#     _SCRIPTS = {}
#
#
#     def __init__(self, instruments = None, scripts = None, name = None, settings = None, log_function = None, data_path = None):
#         """
#         Example of a script that emits a QT signal for the gui
#         Args:
#             name (optional): name of script, if empty same as class name
#             settings (optional): settings for this script, if empty same as default settings
#         """
#         Script.__init__(self, name, settings = settings, instruments = instruments, scripts = scripts, log_function= log_function, data_path = data_path)
#
#         self.text = []
#         self.patch_collection = None
#         self.plot_settings = {}
#
#     def _function(self):
#         """
#         This is the actual function that will be executed. It uses only information that is provided in the settings property
#         will be overwritten in the __init__
#         """
#
#
#
#         self.data = {'nv_locations': [], 'image_data': None, 'extent': None}
#         self.progress = 50
#         self.updateProgress.emit(self.progress)
#         # keep script alive while NVs are selected
#         while not self.data['nv_locations'] and not self._abort:
#             time.sleep(1)
#
#         if self.data['nv_locations']:
#             self.data['nv_locations'] = np.asarray(self.data['nv_locations'])
#             pt = (self.data['nv_locations'][0], self.data['nv_locations'][1])
#             # pt = (self.settings['point']['x'], self.settings['point']['y'])
#             # daq API only accepts either one point and one channel or multiple points and multiple channels
#             pt = np.transpose(np.column_stack((pt[0],pt[1])))
#             pt = (np.repeat(pt, 2, axis=1))
#             # print(pt)
#             self._setup_daq()
#             task = self.daq_out.setup_AO([self.settings['DAQ_channels']['x_ao_channel'], self.settings['DAQ_channels']['y_ao_channel']], pt)
#             self.daq_out.run(task)
#             self.daq_out.waitToFinish(task)
#             self.daq_out.stop(task)
#             self.log('laser set to Vx={:.4}, Vy={:.4}'.format(self.settings['point']['x'], self.settings['point']['y']))
#             self.daq_out = self.instruments['NI6353']['instance']
#
#     def _setup_daq(self):
#         # defines which daqs contain the input and output based on user selection of daq interface
#         # if self.settings['daq_type'] == 'PCI':
#         #     self.daq_out = self.instruments['NI6259']['instance']
#         # elif self.settings['daq_type'] == 'cDAQ':
#         #     self.daq_out = self.instruments['NI9263']['instance']
#         self.daq_out = self.instruments['NI6353']['instance']
#
#     def get_galvo_position(self):
#         """
#         reads the current position from the x and y channels and returns it
#         Returns:
#
#         """
#         galvo_position = self.daq_out.get_analog_voltages([
#             self.settings['DAQ_channels']['x_ao_channel'],
#             self.settings['DAQ_channels']['y_ao_channel']]
#         )
#
#         # if self.settings['daq_type'] == 'PCI':
#         #     galvo_position = self.daq_out.get_analog_voltages([
#         #         self.settings['DAQ_channels']['x_ao_channel'],
#         #         self.settings['DAQ_channels']['y_ao_channel']]
#         #     )
#         # elif self.settings['daq_type'] == 'cDAQ':
#         #     print("WARNING cDAQ doesn't allow to read values")
#         #     galvo_position = []
#
#         return galvo_position
#
#     #must be passed figure with galvo plot on first axis
#     def plot(self, figure_list):
#         '''
#         Plots a dot on top of each selected NV, with a corresponding number denoting the order in which the NVs are
#         listed.
#         Precondition: must have an existing image in figure_list[0] to plot over
#         Args:
#             figure_list:
#         '''
#         # if there is not image data get it from the current plot
#         if not self.data == {} and self.data['image_data'] is None:
#             axes = figure_list[0].axes[0]
#             if len(axes.images) > 0:
#                 self.data['image_data'] = np.array(axes.images[0].get_array())
#                 self.data['extent'] = np.array(axes.images[0].get_extent())
#                 self.plot_settings['cmap'] = axes.images[0].get_cmap().name
#                 self.plot_settings['xlabel'] = axes.get_xlabel()
#                 self.plot_settings['ylabel'] = axes.get_ylabel()
#                 self.plot_settings['title'] = axes.get_title()
#                 self.plot_settings['interpol'] = axes.images[0].get_interpolation()
#         Script.plot(self, figure_list)
#
#     # must be passed figure with galvo plot on first axis
#     def _plot(self, axes_list):
#         '''
#         Plots a dot on top of each selected NV, with a corresponding number denoting the order in which the NVs are
#         listed.
#         Precondition: must have an existing image in figure_list[0] to plot over
#         Args:
#             figure_list:
#         '''
#         axes = axes_list[0]
#         if self.plot_settings:
#             axes.imshow(self.data['image_data'], cmap=self.plot_settings['cmap'],
#                         interpolation=self.plot_settings['interpol'], extent=self.data['extent'])
#             axes.set_xlabel(self.plot_settings['xlabel'])
#             axes.set_ylabel(self.plot_settings['ylabel'])
#             axes.set_title(self.plot_settings['title'])
#         self._update(axes_list)
#
#     def _update(self, axes_list):
#         # note: may be able to use blit to make things faster
#         axes = axes_list[0]
#         patch_size = self.settings['patch_size']
#         # first clear all old patches (circles and numbers), then redraw all
#         if self.patch_collection:
#             try:
#                 self.patch_collection.remove()
#                 for text in self.text:
#                     text.remove()
#             except ValueError:
#                 pass
#         patch_list = []
#         if (self.data['nv_locations'] is not None):
#             for index, pt in enumerate(self.data['nv_locations']):
#                 circ = patches.Circle((pt[0], pt[1]), patch_size, fc='b')
#                 patch_list.append(circ)
#                 # cap number of drawn numbers at 400 since drawing text is extremely slow and they're all so close together
#                 # as to be unreadable anyways
#                 if len(self.data['nv_locations']) <= 400:
#                     text = axes.text(pt[0], pt[1], '{:d}'.format(index),
#                                      horizontalalignment='center',
#                                      verticalalignment='center',
#                                      color='white'
#                                      )
#                     self.text.append(text)
#             # patch collection used here instead of adding individual patches for speed
#             self.patch_collection = matplotlib.collections.PatchCollection(patch_list)
#             axes.add_collection(self.patch_collection)
#
#
#     def toggle_NV(self, pt):
#         '''
#         If there is not currently a selected NV within self.settings[patch_size] of pt, adds it to the selected list. If
#         there is, removes that point from the selected list.
#         Args:
#             pt: the point to add or remove from the selected list
#         Poststate: updates selected list
#         '''
#         print('toggle_NV is called ZQ')
#         if not self.data['nv_locations']:  # if self.data is empty so this is the first point
#             self.data['nv_locations'].append(pt)
#             self.data['image_data'] = None  # clear image data
#         else:
#             # use KDTree to find NV closest to mouse click
#             tree = scipy.spatial.KDTree(self.data['nv_locations'])
#             # does a search with k=1, that is a search for the nearest neighbor, within distance_upper_bound
#             d, i = tree.query(pt, k=1, distance_upper_bound=self.settings['patch_size'])
#
#             # removes NV if previously selected
#             if d is not np.inf:
#                 self.data['nv_locations'].pop(i)
#             # adds NV if not previously selected
#             else:
#                 self.data['nv_locations'].append(pt)


# class SetLaser_cDAQ(SetLaser):
#     """
# This script points the laser to a point
#     """
#
#     _DEFAULT_SETTINGS = [
#         Parameter('point',
#                   [Parameter('x', -0.4, float, 'x-coordinate'),
#                    Parameter('y', -0.4, float, 'y-coordinate')
#                    ]),
#         Parameter('DAQ_channels',
#             [Parameter('x_ao_channel', 'ao0', ['ao0', 'ao1', 'ao2', 'ao3'], 'Daq channel used for x voltage analog output'),
#             Parameter('y_ao_channel', 'ao3', ['ao0', 'ao1', 'ao2', 'ao3'], 'Daq channel used for y voltage analog output')
#             ])
#     ]
#
#     _INSTRUMENTS = {'daq_out':  NI9263}
#
#     def _function(self):
#         """
#         This is the actual function that will be executed. It uses only information that is provided in the settings property
#         will be overwritten in the __init__
#         """
#         pt = (self.settings['point']['x'], self.settings['point']['y'])
#
#         # daq API only accepts either one point and one channel or multiple points and multiple channels
#         pt = np.transpose(np.column_stack((pt[0],pt[1])))
#         pt = (np.repeat(pt, 2, axis=1))
#
#         task = self.instruments['daq_out']['instance'].setup_AO([self.settings['DAQ_channels']['x_ao_channel'], self.settings['DAQ_channels']['y_ao_channel']], pt)
#         self.instruments['daq_out']['instance'].run(task)
#         self.instruments['daq_out']['instance'].waitToFinish(task)
#         self.instruments['daq_out']['instance'].stop(task)
#         self.log('laser set to Vx={:.4}, Vy={:.4}'.format(self.settings['point']['x'], self.settings['point']['y']))

if __name__ == '__main__':
    from pylabcontrol.core import Instrument

    # instruments, instruments_failed = Instrument.load_and_append({'daq':  'NI6259'})

    # script, failed, instruments = Script.load_and_append(script_dict={'SetLaser_cDAQ': 'SetLaser_cDAQ'})
    script, failed, instruments = Script.load_and_append(script_dict={'SetLaser': 'SetLaser'})
    script, failed, instruments = Script.load_and_append(script_dict={'SetConfocal': 'SetConfocal'})

    print(script)
    print(failed)
    # print(instruments)