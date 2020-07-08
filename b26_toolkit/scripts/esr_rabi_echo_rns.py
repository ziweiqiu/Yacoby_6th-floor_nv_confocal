import time
from pylabcontrol.core import Script, Parameter
from b26_toolkit.scripts.esr_RnS import ESR_FastSwp_RnS_FitGuaranteed
from b26_toolkit.scripts.pulse_sequences.rabi import Rabi_RnS
# from b26_toolkit.instruments import R8SMicrowaveGenerator, NI6353, LISE607RTPulseBlaster
from b26_toolkit.scripts.daq_read_counter import Daq_Read_Counter

# class ESRandRabi_RnS(Script):
#     """
#         MagnetSweep2D sweeps the position of the automatic translation stages, in 1D or 2D scans, and records NV fluorescence using Daq_Read_Counter and/or measure the ESR of NV, at each point in the scan.
#         --> Last edited by ZQ 3/15/2019
#     """
#
#     _DEFAULT_SETTINGS = [
#
#         Parameter('esr_settings', [
#             Parameter('esr_mw_pwr', -10, float, 'microwave power for ESR scan'),
#             Parameter('esr_avg_min', 12, int, 'minimum number of esr averages'),
#             Parameter('esr_avg_max', 50, int, 'maximum number of esr averages'),
#             Parameter('esr_cntr_freq', 2.82e9, float, 'center frequency for ESR scan'),
#             Parameter('esr_freq_range', 8.5e7, float, 'frequency range for ESR scan (suggest 6e7 - 9e7)'),
#             Parameter('esr_num_of_pts', 65, int, 'number of frequency points for ESR scan'),
#             Parameter('esr_time_per_pt', 0.02, [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.1],
#                       'integration time for each point in the fast ESR scan (suggest 0.02-0.04)'),
#             Parameter('esrfit_num_of_peaks', 1, [-1, 1, 2],
#                       'specify number of peaks for fitting. if not specifying the number of peaks, choose -1'),
#             Parameter('esrfit_minimum_counts', 0.9, float,
#                       'minumum counts for an ESR to not be considered noise (suggest 0.8 - 1.01 if esr is normalized)'),
#             Parameter('esrfit_contrast_factor', 3.0, float,
#                       'minimum contrast for an ESR to not be considered noise (suggest 3.0-4.0)')
#         ]),
#         Parameter('rabi_settings', [
#             Parameter('tau_min', 15, float, 'minimum time for rabi oscillations (in ns)'),
#             Parameter('tau_max', 1000, float, 'total time of rabi oscillations (in ns)'),
#             Parameter('tau_step', 20., [2.5, 5., 10., 20., 50., 100., 200., 500., 1000., 10000., 100000., 500000.],
#                       'time step increment of rabi pulse duration (in ns)'),
#             Parameter('num_averages', 500000, int, 'number of averages (>100000)'),
#         ]),
#         Parameter('read_out', [
#             Parameter('meas_time', 100, float, 'apd readout duration (in ns)'),
#             Parameter('nv_reset_time', 2000, int, 'time with laser on to reset state'),
#             Parameter('laser_off_time', 1000, int,
#                       'minimum laser off time before taking measurements (ns)'),
#             Parameter('delay_mw_readout', 600, int, 'delay between mw and readout (in ns)'),
#             Parameter('delay_readout', 530, int,
#                       '[ns] delay between laser on and APD readout (given by spontaneous decay rate)')
#         ])
#     ]
#
#     _INSTRUMENTS = {}
#
#     _SCRIPTS = {'esr': ESR_FastSwp_RnS_FitGuaranteed, 'rabi': Rabi_RnS, 'daq_read_counter': Daq_Read_Counter}
#
#     def __init__(self, instruments = None, scripts = None, name = None, settings = None, log_function = None, data_path = None):
#         """
#         Example of a script that makes use of an instrument
#         Args:
#             instruments: instruments the script will make use of
#             name (optional): name of script, if empty same as class name
#             settings (optional): settings for this script, if empty same as default settings
#         """
#
#         # call init of superclass
#         Script.__init__(self, name, settings=settings, instruments=instruments, scripts=scripts,
#                         log_function=log_function, data_path=data_path)
#
#
#     def _function(self):
#         # Set the right parameters for the ESR scan
#         self.scripts['esr'].settings['power_out'] = self.settings['esr_settings']['esr_mw_pwr']
#         self.scripts['esr'].settings['esr_avg_min'] = self.settings['esr_settings']['esr_avg_min']
#         self.scripts['esr'].settings['esr_avg_max'] = self.settings['esr_settings']['esr_avg_max']
#         self.scripts['esr'].settings['freq_start'] = self.settings['esr_settings']['esr_cntr_freq']
#         self.scripts['esr'].settings['freq_stop'] = self.settings['esr_settings']['esr_freq_range']
#         self.scripts['esr'].settings['range_type'] = 'center_range'
#         self.scripts['esr'].settings['freq_points'] = self.settings['esr_settings']['esr_num_of_pts']
#         self.scripts['esr'].settings['time_per_pt'] = self.settings['esr_settings']['esr_time_per_pt']
#         self.scripts['esr'].settings['fit_constants']['num_of_peaks'] = self.settings['esr_settings'][
#             'esrfit_num_of_peaks']
#         self.scripts['esr'].settings['fit_constants']['minimum_counts'] = self.settings['esr_settings'][
#             'esrfit_minimum_counts']
#         self.scripts['esr'].settings['fit_constants']['contrast_factor'] = self.settings['esr_settings'][
#             'esrfit_contrast_factor']
#
#         # Start ESR
#         print('==> Start measuring ESR...')
#         print('self._plot_refresh', self._plot_refresh)
#         self.scripts['esr'].run()
#         esr_fit_data = self.scripts['esr'].data['fit_params']
#         print('     len(esr_fit_data) =  ', esr_fit_data)
#
#         if self.scripts['esr'].data['fit_params'] is not None:
#             if len(self.scripts['esr'].data['fit_params']) == 4:
#                 self.rabi_frequency = self.scripts['esr'].data['fit_params'][2]
#             elif len(self.scripts['esr'].data['fit_params']) == 6:
#                 self.rabi_frequency = self.scripts['esr'].data['fit_params'][4]
#             else:
#                 self.log('Could not get fit parameters from esr script. Experiment aborted.')
#                 print('Could not get fit parameters from esr script. Experiment aborted.')
#                 self._abort = True
#                 # raise RuntimeError('Could not get fit parameters from esr script')
#
#             # if self.rabi_frequency < self.scripts['esr'].settings['freq_start']:
#             if self.rabi_frequency < (
#                     self.scripts['esr'].settings['freq_start'] - 0.5 * self.scripts['esr'].settings['freq_stop']):
#                 self.log(
#                     'Resonance frequency found ({:0.2e}) Hz was below esr sweep range, aborting rabi attempt'.format(
#                         self.rabi_frequency))
#             # elif self.rabi_frequency > self.scripts['esr'].settings['freq_stop']:
#             elif self.rabi_frequency > (
#                     self.scripts['esr'].settings['freq_start'] + 0.5 * self.scripts['esr'].settings['freq_stop']):
#                 self.log(
#                     'Resonance frequency found ({:0.2e}) Hz was above esr sweep range, aborting rabi attempt'.format(
#                         self.rabi_frequency))
#             else:
#                 self.log(
#                     'Starting rabi with frequency {:.4f} GHz and power {:.2f} dBm'.format(self.rabi_frequency / 1E9,
#                                                                                           self.settings['esr_settings'][
#                                                                                               'esr_mw_pwr']))
#                 print('Starting rabi with frequency {:.4f} GHz and power {:.2f} dBm'.format(self.rabi_frequency / 1E9,
#                                                                                             self.settings[
#                                                                                                 'esr_settings'][
#                                                                                                 'esr_mw_pwr']))
#                 self.scripts['rabi'].settings['mw_pulses']['mw_power'] = self.settings['esr_settings']['esr_mw_pwr']
#                 self.scripts['rabi'].settings['mw_pulses']['mw_frequency'] = float(self.rabi_frequency)
#                 self.scripts['rabi'].run()
#         else:
#             self.log('No resonance frequency found skipping rabi attempt')
#
#     # def _plot(self, axes_list, data=None):
#     #     # COMMENT_ME
#     #
#     #     if data is None:
#     #         data = self.data
#     #     if self.settings['scan_axis'] in ['x', 'y', 'z'] and self.settings['to-do'] == 'sweep':
#     #         print('(1D plot)')
#     #         if self.data['counts'] is not None:
#     #             lbls1 = ['magnet position ' + self.settings['scan_axis'] + ' [mm]',
#     #                      'counts [kcps]', 'Fluorescence']
#     #             plot_magnet_sweep1D_Fluor([axes_list[0]], self.data['positions'], np.array(self.data['counts']),
#     #                                       lbls1, x_r=self.data['positions_r'], y1_r=np.array(self.data['counts_r']))
#     #
#     #         if self.settings['exp_settings']['to_plot'] ==  'contrast': # to plot contrast
#     #             if self.data['esr_fo'] is not None and self.data['esr_ctrst'] is not None:
#     #                 lbls2 = ['magnet position ' + self.settings['scan_axis'] + ' [mm]', 'f0 [Hz]',
#     #                          'contrast', 'ESR']
#     #                 plot_magnet_sweep1D_ESR([axes_list[2], axes_list[3]], self.data['positions'],
#     #                                         np.array(self.data['esr_fo']),
#     #                                         np.array(self.data['esr_ctrst']), lbls2, x_r=self.data['positions_r'],
#     #                                         y1_r=np.array(self.data['esr_fo_r']), y2_r=np.array(self.data['esr_ctrst_r']))
#     #             if self.data['esr2_fo'] is not None and self.data['esr2_ctrst'] is not None:
#     #                 lbls3 = ['magnet position ' + self.settings['scan_axis'] + ' [mm]', 'f0 [Hz]',
#     #                          'contrast', 'ESR 2']
#     #                 plot_magnet_sweep1D_ESR([axes_list[4], axes_list[5]], self.data['positions'],
#     #                                         np.array(self.data['esr2_fo']),
#     #                                         np.array(self.data['esr2_ctrst']), lbls3, x_r=self.data['positions_r'],
#     #                                         y1_r=np.array(self.data['esr2_fo_r']), y2_r=np.array(self.data['esr2_ctrst_r']))
#     #         else:
#     #             if self.data['esr_fo'] is not None and self.data['esr_wo'] is not None:
#     #                 lbls2 = ['magnet position ' + self.settings['scan_axis'] + ' [mm]', 'f0 [Hz]',
#     #                          'wo[Hz]', 'ESR']
#     #                 plot_magnet_sweep1D_ESR([axes_list[2], axes_list[3]], self.data['positions'],
#     #                                         np.array(self.data['esr_fo']),
#     #                                         np.array(self.data['esr_wo']), lbls2, x_r=self.data['positions_r'],
#     #                                         y1_r=np.array(self.data['esr_fo_r']), y2_r=np.array(self.data['esr_wo_r']))
#     #             if self.data['esr2_fo'] is not None and self.data['esr2_wo'] is not None:
#     #                 lbls3 = ['magnet position ' + self.settings['scan_axis'] + ' [mm]', 'f0 [Hz]',
#     #                          'wo[Hz]', 'ESR 2']
#     #                 plot_magnet_sweep1D_ESR([axes_list[4], axes_list[5]], self.data['positions'],
#     #                                         np.array(self.data['esr2_fo']),
#     #                                         np.array(self.data['esr2_wo']), lbls3, x_r=self.data['positions_r'],
#     #                                         y1_r=np.array(self.data['esr2_fo_r']), y2_r=np.array(self.data['esr2_wo_r']))
#     #
#     #     elif self.settings['scan_axis'] in ['xy', 'yx', 'yz', 'zy', 'zx', 'xz'] and self.settings['to-do'] == 'sweep':
#     #         print('(2D plot)')
#     #         if self.data['counts'] is not None:
#     #             plot_magnet_sweep2D_Fluor(data['counts'], data['bounds'], axes_list[0], axes_labels=self.settings['scan_axis'], title = 'Fluor Forward', colorbar_name = 'kcps')
#     #         if self.data['counts_r'] is not None:
#     #             plot_magnet_sweep2D_Fluor(data['counts_r'], data['bounds'], axes_list[2], axes_labels=self.settings['scan_axis'], title = 'Fluor Backward', colorbar_name = 'kcps')
#     #         if self.data['esr_fo'] is not None:
#     #             plot_magnet_sweep2D_Fluor(data['esr_fo']/1e6, data['bounds'], axes_list[3], axes_labels=self.settings['scan_axis'], title = 'ESR Forward', colorbar_name = 'MHz')
#     #         if self.data['esr_fo_r'] is not None:
#     #             plot_magnet_sweep2D_Fluor(data['esr_fo_r']/1e6, data['bounds'], axes_list[4], axes_labels=self.settings['scan_axis'], title = 'ESR Backward', colorbar_name = 'MHz')
#     #         # elif self.data['esr_wo'] is not None:
#     #         #     plot_magnet_sweep2D_Fluor(data['esr_wo']/1e6, data['bounds'], axes_list[4], axes_labels=self.settings['scan_axis'], title = 'ESR width Forward', colorbar_name = 'MHz')
#     #         if self.data['esr2_fo'] is not None:
#     #             plot_magnet_sweep2D_Fluor(data['esr2_fo']/1e6, data['bounds'], axes_list[5], axes_labels=self.settings['scan_axis'], title = 'ESR2 Forward', colorbar_name = 'MHz')
#     #         if self.data['esr2_fo_r'] is not None:
#     #             plot_magnet_sweep2D_Fluor(data['esr2_fo_r']/1e6, data['bounds'], axes_list[6], axes_labels=self.settings['scan_axis'], title = 'ESR2 Backward', colorbar_name = 'MHz')
#     #         # elif self.data['esr2_wo'] is not None:
#     #         #     plot_magnet_sweep2D_Fluor(data['esr2_wo']/1e6, data['bounds'], axes_list[6], axes_labels=self.settings['scan_axis'], title = 'ESR2 width Forward', colorbar_name = 'MHz')
#
#     def _update_plot(self, axes_list):
#
#         print('_update_plot')
#
#
#         if self._current_subscript_stage['current_subscript'] is self.scripts['esr'] and self.scripts['esr'].is_running:
#             self.scripts['esr']._update_plot([axes_list[1]])
#
#
#     def get_axes_layout(self, figure_list):
#         """
#         returns the axes objects the script needs to plot its data
#         this overwrites the default get_axis_layout in PyLabControl.src.core.scripts
#         Args:
#             figure_list: a list of figure objects
#         Returns:
#             axes_list: a list of axes objects
#
#         """
#         axes_list = []
#
#         if self._plot_refresh is True:
#             for fig in figure_list:
#                 fig.clf()
#             # 5 subplots in total
#             axes_list.append(figure_list[0].add_subplot(111))  # axes_list[0]
#             axes_list.append(figure_list[1].add_subplot(111))  # axes_list[1]
#
#
#         else:
#             axes_list.append(figure_list[0].axes[0])
#             axes_list.append(figure_list[1].axes[0])
#
#         return axes_list

class ESRandRabi_RnS(Script):
    """
        MagnetSweep2D sweeps the position of the automatic translation stages, in 1D or 2D scans, and records NV fluorescence using Daq_Read_Counter and/or measure the ESR of NV, at each point in the scan.
        --> Last edited by ZQ 3/15/2019
    """

    _DEFAULT_SETTINGS = [
        Parameter('to-do', 'sweep', ['move', 'sweep'], 'Choose to move to a point or do a magnet sweep'),
        Parameter('servo_initial',
                  [Parameter('initialize', True, bool,
                             'whether or not to intialize the servo position before sweeping? (highly recommended)'),
                   Parameter('Xservo', 9.0, float, 'initial position of Xservo'),
                   Parameter('Yservo', 4.0, float, 'initial position of Yservo'),
                   Parameter('Zservo', 5.0, float, 'initial position of Zservo'),
                   Parameter('moving_velocity', 0.5, float, 'servo moving velocity (mm/s)'),
                   Parameter('Xservo_min', 0.0, float, 'minimum allowed position of Xservo'),
                   Parameter('Xservo_max', 23.0, float, 'maximum allowed position of Xservo'),
                   Parameter('Yservo_min', 0.0, float, 'minimum allowed position of Yservo'),
                   Parameter('Yservo_max', 13.0, float, 'maximum allowed position of Yservo'),
                   Parameter('Zservo_min', 0.0, float, 'minimum allowed position of Zservo'),
                   Parameter('Zservo_max', 25.0, float, 'maximum allowed position of Zservo'),
                   ]),

        Parameter('scan_axis', 'xy', ['xy', 'yx', 'xz', 'zx', 'yz', 'zy', 'x', 'y', 'z'],
                  'Choose 2D or 1D magnet sweep to perform'),
        Parameter('move_to',
                  [Parameter('x', 15.0, float, 'move to x-coordinate [mm]'),
                   Parameter('y', 10.0, float, 'move to y-coordinate [mm]'),
                   Parameter('z', 10.0, float, 'move to z-coordinate [mm]')
                   ]),
        Parameter('sweep_center',
                  [Parameter('x', 15.0, float, 'x-coordinate [mm] of the sweep center'),
                   Parameter('y', 10.0, float, 'y-coordinate [mm] of the sweep center'),
                   Parameter('z', 5.0, float, 'z-coordinate [mm] of the sweep center')
                   ]),
        Parameter('sweep_span',
                  [Parameter('x', 10.0, float, 'x-coordinate [mm]'),
                   Parameter('y', 10.0, float, 'y-coordinate [mm]'),
                   Parameter('z', 5.0, float, 'z-coordinate [mm]')
                   ]),
        Parameter('num_points',
                  [Parameter('x', 11, int, 'number of x points to scan'),
                   Parameter('y', 11, int, 'number of y points to scan'),
                   Parameter('z', 11, int, 'number of z points to scan')
                   ]),
        Parameter('exp_to_do', [Parameter('fluorescence', True, bool, 'measure the NV fluorescence'),
                                Parameter('esr', True, bool, 'measure the ESR of NV'),
                                Parameter('esr2', True, bool,
                                          'measure the ESR of NV at two different frequencies')]),
        Parameter('exp_settings', [
            Parameter('fluorescence_time_per_pt', 0.4, float, 'time for fluorescence measurement at each point (s)'),
            Parameter('esr_mw_pwr', -10, float, 'microwave power for ESR scan'),
            Parameter('esr_avg_min', 12, int, 'minimum number of esr averages'),
            Parameter('esr_avg_max', 50, int, 'maximum number of esr averages'),
            Parameter('esr_cntr_freq', 2.82e9, float, 'center frequency for ESR scan'),
            Parameter('esr_freq_range', 8.5e7, float, 'frequency range for ESR scan (suggest 6e7 - 9e7)'),
            Parameter('esr2_cntr_freq', 2.99e9, float, 'center frequency for the second ESR scan'),
            Parameter('esr2_freq_range', 8.5e7, float, 'frequency range for the second ESR scan(suggest 6e7 - 9e7)'),
            Parameter('esr_num_of_pts', 65, int, 'number of frequency points for ESR scan'),
            Parameter('esr_time_per_pt', 0.02, [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.1],
                      'integration time for each point in the fast ESR scan (suggest 0.02-0.04)'),
            Parameter('esrfit_num_of_peaks', 1, [-1, 1, 2],
                      'specify number of peaks for fitting. if not specifying the number of peaks, choose -1'),
            Parameter('esrfit_minimum_counts', 0.9, float,
                      'minumum counts for an ESR to not be considered noise (suggest 0.8 - 1.01 if esr is normalized)'),
            Parameter('esrfit_contrast_factor', 3.0, float,
                      'minimum contrast for an ESR to not be considered noise (suggest 3.0-4.0)'),
            Parameter('to_plot', 'contrast', ['fwhm', 'contrast'], 'choose to plot fwhm or contrast in 1D sweep')
        ]),
        Parameter('tracking_settings', [Parameter('track_focus', 'autofocus', ['optimize_z', 'autofocus', 'None'],
                                                  'choose the method for tracking (optimize_z is recommended)'),
                                        Parameter('track_focus_every_N', 18, int, 'track every N points'),
                                        Parameter('track_to_nv', True, bool,
                                                  'check to use find_nv to track to the NV'),
                                        Parameter('track_to_nv_every_N', 5, int, 'track every N points'),
                                        Parameter('track_frequency', True, bool,
                                                  'keep track of the frequency and set it to the central frequency of the next ESR scan (recommended)'),
                                        Parameter('track_frequency_every_N', 1, int, 'track every N points')]),
        Parameter('optimize_z_settings', [Parameter('sweep_range', 0.6, float, 'z voltage range for optimizing scan (suggest 0.6)'),
                                          Parameter('num_points', 41, float, 'number of z points to scan (suggest 41)'),
                                          ]),
        Parameter('autofocus_settings', [Parameter('scan_width', 0.6, float, 'z voltage range for optimizing scan (suggest 0.6-0.9)'),
                                         Parameter('num_sweep_points', 6, int,
                                                   'number of values to sweep between min and max voltage (suggest 6-10)'),
                                         ]),
        Parameter('find_nv_settings',
                  [Parameter('sweep_range', 0.35, float, 'voltage range to sweep over to find a max (suggest 0.4)'),
                   Parameter('num_points', 61, int, 'number of points to sweep in the sweep range'),
                   Parameter('nv_size', 21, int,
                             'size of nv in pixels - need to be refined!! needs to be odd number!!!'),
                   Parameter('min_mass', 80, int, 'TEMP: brightness of nv - need to be refined!! (suggest 60-100)'),
                   ])
    ]

    _INSTRUMENTS = {}

    _SCRIPTS = {'esr': ESR_FastSwp_RnS_FitGuaranteed}

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

    def _get_instr(self):
        """
        Assigns an instrument relevant to the 1D scan axis.
        """
        if self.settings['scan_axis'] == 'x':
            return self.instruments['XServo']['instance']
        elif self.settings['scan_axis'] == 'y':
            return self.instruments['YServo']['instance']
        elif self.settings['scan_axis'] == 'z':
            return self.instruments['ZServo']['instance']

    def _get_instr_2D(self):
        """
        Assigns an instrument relevant to the 2D scan axis.
        """
        if self.settings['scan_axis'] == 'xy':
            return self.instruments['XServo']['instance'], self.instruments['YServo']['instance']
        elif self.settings['scan_axis'] == 'yx':
            return self.instruments['YServo']['instance'], self.instruments['XServo']['instance']
        elif self.settings['scan_axis'] == 'xz':
            return self.instruments['XServo']['instance'], self.instruments['ZServo']['instance']
        elif self.settings['scan_axis'] == 'zx':
            return self.instruments['ZServo']['instance'], self.instruments['XServo']['instance']
        elif self.settings['scan_axis'] == 'yz':
            return self.instruments['YServo']['instance'], self.instruments['ZServo']['instance']
        elif self.settings['scan_axis'] == 'zy':
            return self.instruments['ZServo']['instance'], self.instruments['YServo']['instance']


    @staticmethod
    def pts_to_extent(pta, ptb):
        """

        Args:
            pta: point a
            ptb: point b
            roi_mode:   mode how to calculate region of interest
                        corner: pta and ptb are diagonal corners of rectangle.
                        center: pta is center and ptb is extend or rectangle

        Returns: extend of region of interest [xVmin, xVmax, yVmax, yVmin]

        """
        xVmin = pta['x'] - float(ptb['x']) / 2.
        xVmax = pta['x'] + float(ptb['x']) / 2.
        yVmin = pta['y'] - float(ptb['y']) / 2.
        yVmax = pta['y'] + float(ptb['y']) / 2.
        zVmin = pta['z'] - float(ptb['z']) / 2.
        zVmax = pta['z'] + float(ptb['z']) / 2.

        return [xVmin, xVmax, yVmin, yVmax, zVmin, zVmax]



    def do_esr(self, esr_cntr_freq, esr_freq_range, label = None, index = -1, verbose=False):

        # update the tag of the esr script
        if label is not None and index >= 0:
            self.scripts['esr'].settings['tag'] = label + '_ind' + str(index)
        elif label is not None:
            self.scripts['esr'].settings['tag'] = label
        elif index >= 0:
            self.scripts['esr'].settings['tag'] = 'esr_ind' + str(index)

        # set the right parameters for the ESR scan
        self.scripts['esr'].settings['power_out'] = self.settings['exp_settings']['esr_mw_pwr']
        self.scripts['esr'].settings['esr_avg_min'] = self.settings['exp_settings']['esr_avg_min']
        self.scripts['esr'].settings['esr_avg_max'] = self.settings['exp_settings']['esr_avg_max']
        self.scripts['esr'].settings['freq_start'] = float(esr_cntr_freq)
        self.scripts['esr'].settings['freq_stop'] = float(esr_freq_range)
        self.scripts['esr'].settings['range_type'] = 'center_range'
        self.scripts['esr'].settings['freq_points'] = self.settings['exp_settings']['esr_num_of_pts']
        self.scripts['esr'].settings['time_per_pt'] = self.settings['exp_settings']['esr_time_per_pt']
        self.scripts['esr'].settings['fit_constants']['num_of_peaks'] = self.settings['exp_settings'][
            'esrfit_num_of_peaks']
        self.scripts['esr'].settings['fit_constants']['minimum_counts'] = self.settings['exp_settings'][
            'esrfit_minimum_counts']
        self.scripts['esr'].settings['fit_constants']['contrast_factor'] = self.settings['exp_settings'][
            'esrfit_contrast_factor']
        print('==> Start measuring ESR...')
        self.scripts['esr'].run()
        esr_fit_data = self.scripts['esr'].data['fit_params']
        if verbose:
            print('len(esr_fit_data) =  ', esr_fit_data)

        return esr_fit_data

    def _function(self):
        esr_fit_data = self.do_esr(2.87E9,
                                   self.settings['exp_settings']['esr_freq_range'], label='esr1', index=2)

    # def _plot(self, axes_list, data=None):
        # COMMENT_ME

        # if data is None:
        #     data = self.data
        # if self.settings['scan_axis'] in ['x', 'y', 'z'] and self.settings['to-do'] == 'sweep':
        #     print('(1D plot)')
        #     if self.data['counts'] is not None:
        #         lbls1 = ['magnet position ' + self.settings['scan_axis'] + ' [mm]',
        #                  'counts [kcps]', 'Fluorescence']
        #         plot_magnet_sweep1D_Fluor([axes_list[0]], self.data['positions'], np.array(self.data['counts']),
        #                                   lbls1, x_r=self.data['positions_r'], y1_r=np.array(self.data['counts_r']))
        #
        #     if self.settings['exp_settings']['to_plot'] ==  'contrast': # to plot contrast
        #         if self.data['esr_fo'] is not None and self.data['esr_ctrst'] is not None:
        #             lbls2 = ['magnet position ' + self.settings['scan_axis'] + ' [mm]', 'f0 [Hz]',
        #                      'contrast', 'ESR']
        #             plot_magnet_sweep1D_ESR([axes_list[2], axes_list[3]], self.data['positions'],
        #                                     np.array(self.data['esr_fo']),
        #                                     np.array(self.data['esr_ctrst']), lbls2, x_r=self.data['positions_r'],
        #                                     y1_r=np.array(self.data['esr_fo_r']), y2_r=np.array(self.data['esr_ctrst_r']))
        #         if self.data['esr2_fo'] is not None and self.data['esr2_ctrst'] is not None:
        #             lbls3 = ['magnet position ' + self.settings['scan_axis'] + ' [mm]', 'f0 [Hz]',
        #                      'contrast', 'ESR 2']
        #             plot_magnet_sweep1D_ESR([axes_list[4], axes_list[5]], self.data['positions'],
        #                                     np.array(self.data['esr2_fo']),
        #                                     np.array(self.data['esr2_ctrst']), lbls3, x_r=self.data['positions_r'],
        #                                     y1_r=np.array(self.data['esr2_fo_r']), y2_r=np.array(self.data['esr2_ctrst_r']))
        #     else:
        #         if self.data['esr_fo'] is not None and self.data['esr_wo'] is not None:
        #             lbls2 = ['magnet position ' + self.settings['scan_axis'] + ' [mm]', 'f0 [Hz]',
        #                      'wo[Hz]', 'ESR']
        #             plot_magnet_sweep1D_ESR([axes_list[2], axes_list[3]], self.data['positions'],
        #                                     np.array(self.data['esr_fo']),
        #                                     np.array(self.data['esr_wo']), lbls2, x_r=self.data['positions_r'],
        #                                     y1_r=np.array(self.data['esr_fo_r']), y2_r=np.array(self.data['esr_wo_r']))
        #         if self.data['esr2_fo'] is not None and self.data['esr2_wo'] is not None:
        #             lbls3 = ['magnet position ' + self.settings['scan_axis'] + ' [mm]', 'f0 [Hz]',
        #                      'wo[Hz]', 'ESR 2']
        #             plot_magnet_sweep1D_ESR([axes_list[4], axes_list[5]], self.data['positions'],
        #                                     np.array(self.data['esr2_fo']),
        #                                     np.array(self.data['esr2_wo']), lbls3, x_r=self.data['positions_r'],
        #                                     y1_r=np.array(self.data['esr2_fo_r']), y2_r=np.array(self.data['esr2_wo_r']))
        #
        # elif self.settings['scan_axis'] in ['xy', 'yx', 'yz', 'zy', 'zx', 'xz'] and self.settings['to-do'] == 'sweep':
        #     print('(2D plot)')
        #     if self.data['counts'] is not None:
        #         plot_magnet_sweep2D_Fluor(data['counts'], data['bounds'], axes_list[0], axes_labels=self.settings['scan_axis'], title = 'Fluor Forward', colorbar_name = 'kcps')
        #     if self.data['counts_r'] is not None:
        #         plot_magnet_sweep2D_Fluor(data['counts_r'], data['bounds'], axes_list[2], axes_labels=self.settings['scan_axis'], title = 'Fluor Backward', colorbar_name = 'kcps')
        #     if self.data['esr_fo'] is not None:
        #         plot_magnet_sweep2D_Fluor(data['esr_fo']/1e6, data['bounds'], axes_list[3], axes_labels=self.settings['scan_axis'], title = 'ESR Forward', colorbar_name = 'MHz')
        #     if self.data['esr_fo_r'] is not None:
        #         plot_magnet_sweep2D_Fluor(data['esr_fo_r']/1e6, data['bounds'], axes_list[4], axes_labels=self.settings['scan_axis'], title = 'ESR Backward', colorbar_name = 'MHz')
        #     # elif self.data['esr_wo'] is not None:
        #     #     plot_magnet_sweep2D_Fluor(data['esr_wo']/1e6, data['bounds'], axes_list[4], axes_labels=self.settings['scan_axis'], title = 'ESR width Forward', colorbar_name = 'MHz')
        #     if self.data['esr2_fo'] is not None:
        #         plot_magnet_sweep2D_Fluor(data['esr2_fo']/1e6, data['bounds'], axes_list[5], axes_labels=self.settings['scan_axis'], title = 'ESR2 Forward', colorbar_name = 'MHz')
        #     if self.data['esr2_fo_r'] is not None:
        #         plot_magnet_sweep2D_Fluor(data['esr2_fo_r']/1e6, data['bounds'], axes_list[6], axes_labels=self.settings['scan_axis'], title = 'ESR2 Backward', colorbar_name = 'MHz')
        #     # elif self.data['esr2_wo'] is not None:
        #     #     plot_magnet_sweep2D_Fluor(data['esr2_wo']/1e6, data['bounds'], axes_list[6], axes_labels=self.settings['scan_axis'], title = 'ESR2 width Forward', colorbar_name = 'MHz')

    def _update_plot(self, axes_list):

        print('_update_plot')
        if self._current_subscript_stage['current_subscript'] is self.scripts['esr'] and self.scripts['esr'].is_running:
            print('updating esr plot now')

        # if self._current_subscript_stage['current_subscript'] is self.scripts['daq_read_counter'] and self.scripts['daq_read_counter'].is_running:
        #     self.scripts['daq_read_counter']._plot([axes_list[1]])
        # elif self._current_subscript_stage['current_subscript'] is self.scripts['esr'] and self.scripts['esr'].is_running:
        #     print('updating esr plot now')
        #     self.scripts['esr']._update_plot([axes_list[1]])
        # elif self._current_subscript_stage['current_subscript'] is self.scripts['find_nv'] and self.scripts['find_nv'].is_running:
        #     if self.flag_find_nv_plot:
        #         # print('self.flag_find_nv_plot is', self.flag_find_nv_plot)
        #         self.scripts['find_nv']._plot([axes_list[1]], colorbar=0) # this is to remove colorbar
        #         self.flag_find_nv_plot = False
        #     else:
        #         self.scripts['find_nv']._update_plot([axes_list[1]])
        # elif self._current_subscript_stage['current_subscript'] is self.scripts['optimize_z'] and self.scripts['optimize_z'].is_running:
        #     # print('optimize_z is running, update plot')
        #     if self.flag_optimize_z_plot:
        #         self.scripts['optimize_z']._plot([axes_list[1]])
        #         self.flag_optimize_z_plot = False
        #     else:
        #         self.scripts['optimize_z']._update_plot ([axes_list[1]])
        # elif self._current_subscript_stage['current_subscript'] is self.scripts['autofocus'] and self.scripts['autofocus'].is_running:
        #     # print('autofocus is running, update plot')
        #     if self.flag_autofocus_plot:
        #         # print('self.flag_autofocus_plot is', self.flag_autofocus_plot)
        #         self.scripts['autofocus']._plot([axes_list[0], axes_list[1]], colorbar=0)
        #         self.flag_autofocus_plot = False
        #         self.flag_image0_update_plot = False
        #     else:
        #         self.scripts['autofocus']._update_plot([axes_list[0], axes_list[1]])
        #         self.flag_image0_update_plot = False
        # else:
        #     if self.settings['scan_axis'] in ['x', 'y', 'z'] and self.settings['to-do'] == 'sweep':
        #         print('(updating 1D plot)')
        #
        #         if self.data['counts'] is not None:
        #             lbls1 = ['magnet position ' + self.settings['scan_axis'] + ' [mm]',
        #                      'counts [kcps]', 'Fluorescence']
        #             plot_magnet_sweep1D_Fluor([axes_list[0]], self.data['positions'], np.array(self.data['counts']),
        #                                       lbls1, x_r=self.data['positions_r'], y1_r=np.array(self.data['counts_r']))
        #         if self.settings['exp_settings']['to_plot'] == 'contrast':  # to plot contrast
        #             if self.data['esr_fo'] is not None and self.data['esr_ctrst'] is not None:
        #                 lbls2 = ['magnet position ' + self.settings['scan_axis'] + ' [mm]', 'f0 [Hz]',
        #                          'contrast', 'ESR']
        #                 plot_magnet_sweep1D_ESR([axes_list[2], axes_list[3]], self.data['positions'],
        #                                         np.array(self.data['esr_fo']),
        #                                         np.array(self.data['esr_ctrst']), lbls2, x_r=self.data['positions_r'],
        #                                         y1_r=np.array(self.data['esr_fo_r']),
        #                                         y2_r=np.array(self.data['esr_ctrst_r']))
        #             if self.data['esr2_fo'] is not None and self.data['esr2_ctrst'] is not None:
        #                 lbls3 = ['magnet position ' + self.settings['scan_axis'] + ' [mm]', 'f0 [Hz]',
        #                          'contrast', 'ESR 2']
        #                 plot_magnet_sweep1D_ESR([axes_list[4], axes_list[5]], self.data['positions'],
        #                                         np.array(self.data['esr2_fo']),
        #                                         np.array(self.data['esr2_ctrst']), lbls3, x_r=self.data['positions_r'],
        #                                         y1_r=np.array(self.data['esr2_fo_r']),
        #                                         y2_r=np.array(self.data['esr2_ctrst_r']))
        #         else: # to plot width
        #             if self.data['esr_fo'] is not None and self.data['esr_wo'] is not None:
        #                 lbls2 = ['magnet position ' + self.settings['scan_axis'] + ' [mm]', 'f0 [Hz]',
        #                          'wo[Hz]', 'ESR']
        #                 plot_magnet_sweep1D_ESR([axes_list[2], axes_list[3]], self.data['positions'],
        #                                         np.array(self.data['esr_fo']),
        #                                         np.array(self.data['esr_wo']), lbls2, x_r=self.data['positions_r'],
        #                                         y1_r=np.array(self.data['esr_fo_r']),
        #                                         y2_r=np.array(self.data['esr_wo_r']))
        #             if self.data['esr2_fo'] is not None and self.data['esr2_wo'] is not None:
        #                 lbls3 = ['magnet position ' + self.settings['scan_axis'] + ' [mm]', 'f0 [Hz]',
        #                          'wo[Hz]', 'ESR 2']
        #                 plot_magnet_sweep1D_ESR([axes_list[4], axes_list[5]], self.data['positions'],
        #                                         np.array(self.data['esr2_fo']),
        #                                         np.array(self.data['esr2_wo']), lbls3, x_r=self.data['positions_r'],
        #                                         y1_r=np.array(self.data['esr2_fo_r']),
        #                                         y2_r=np.array(self.data['esr2_wo_r']))
        #
        #     elif self.settings['scan_axis'] in ['xy', 'yx', 'yz', 'zy', 'zx', 'xz'] and self.settings['to-do'] == 'sweep':
        #         print('(updating 2D plot)')
        #         if self.data['counts'] is not None:
        #             if self.flag_image0_update_plot:
        #                 update_magnet_sweep2D_Fluor(self.data['counts'], axes_list[0])
        #             else:
        #                 plot_magnet_sweep2D_Fluor(self.data['counts'], self.data['bounds'], axes_list[0],
        #                                           axes_labels=self.settings['scan_axis'], title='Fluor Forward',
        #                                           colorbar_name='kcps')
        #                 self.flag_image0_update_plot = True
        #         if self.data['counts_r'] is not None:
        #             update_magnet_sweep2D_Fluor(self.data['counts_r'],axes_list[2])
        #         if self.data['esr_fo'] is not None:
        #             update_magnet_sweep2D_Fluor(self.data['esr_fo']/1e6, axes_list[3])
        #         if self.data['esr_fo_r'] is not None:
        #             update_magnet_sweep2D_Fluor(self.data['esr_fo_r']/1e6, axes_list[4])
        #         if self.data['esr2_fo'] is not None:
        #             update_magnet_sweep2D_Fluor(self.data['esr2_fo']/1e6, axes_list[5])
        #         if self.data['esr2_fo_r'] is not None:
        #             update_magnet_sweep2D_Fluor(self.data['esr2_fo_r']/1e6, axes_list[6])

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
        if self.settings['scan_axis'] in ['x', 'y', 'z']:
            if self._plot_refresh is True:
                for fig in figure_list:
                    fig.clf()
                # 5 subplots in total
                axes_list.append(figure_list[0].add_subplot(131))  # axes_list[0]
                axes_list.append(figure_list[1].add_subplot(111))  # axes_list[1]
                axes_list.append(figure_list[0].add_subplot(232))  # axes_list[2]
                axes_list.append(figure_list[0].add_subplot(235))  # axes_list[3]
                axes_list.append(figure_list[0].add_subplot(233))  # axes_list[4]
                axes_list.append(figure_list[0].add_subplot(236))  # axes_list[5]

            else:
                axes_list.append(figure_list[0].axes[0])
                axes_list.append(figure_list[1].axes[0])
                axes_list.append(figure_list[0].axes[1])
                axes_list.append(figure_list[0].axes[2])
                axes_list.append(figure_list[0].axes[3])
                axes_list.append(figure_list[0].axes[4])

            return axes_list

        elif self.settings['scan_axis'] in ['xy', 'yx', 'yz', 'zy', 'zx', 'xz']:
            if self._plot_refresh is True:
                for fig in figure_list:
                    fig.clf()
                # 6 subplots in total
                axes_list.append(figure_list[0].add_subplot(231))  # axes_list[0]
                axes_list.append(figure_list[1].add_subplot(111))  # axes_list[1]
                axes_list.append(figure_list[0].add_subplot(234))  # axes_list[2]
                axes_list.append(figure_list[0].add_subplot(232))  # axes_list[3]
                axes_list.append(figure_list[0].add_subplot(235))  # axes_list[4]
                axes_list.append(figure_list[0].add_subplot(233))  # axes_list[5]
                axes_list.append(figure_list[0].add_subplot(236))  # axes_list[6]

            else:
                axes_list.append(figure_list[0].axes[0])
                axes_list.append(figure_list[1].axes[0])
                axes_list.append(figure_list[0].axes[1])
                axes_list.append(figure_list[0].axes[2])
                axes_list.append(figure_list[0].axes[3])
                axes_list.append(figure_list[0].axes[4])
                axes_list.append(figure_list[0].axes[5])

            return axes_list











