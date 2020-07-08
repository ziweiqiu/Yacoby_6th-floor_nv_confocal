"""
    This file is part of b26_toolkit, a PyLabControl add-on for experiments in Harvard LISE B26.
    Copyright (C) <2016>  Arthur Safira, Jan Gieseler, Aaron Kabcenell

    PyLabControl is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    PyLabControl is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with PyLabControl.  If not, see <http://www.gnu.org/licenses/>.
"""

# clr is python for .net
import clr # run pip install pythonnet
import sys, os
from pylabcontrol.core.read_write_functions import get_config_value


dll_path = get_config_value('KINESIS_DLL_PATH', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'config.txt'))
# print('dll_path:',dll_path)
sys.path.insert(0,dll_path)
# JG: July 27 2016 uncommented folowing line: don't use import *!
# from PyLabControl.src.core.instruments import *
from pylabcontrol.core import Parameter, Instrument

# ctypes DLL load failed: Probably a C++ dll was provided, which is incompatable with ctypes, possibly due to name
# mangling. Instead, we use the .net framework with python for .net to interface with the dll
#ctypes.cdll.LoadLibrary("Thorlabs.MotionControl.TCube.DCServo.dll")

# makes each dll, corresponding to a namespace, avaliable to python at runtime
clr.AddReference('ThorLabs.MotionControl.DeviceManagerCLI')
clr.AddReference('Thorlabs.MotionControl.TCube.DCServoCLI')
clr.AddReference('System')

# imports classes from the namespaces. All read as unresolved references because python doesn't know about the dlls
# until runtime
from Thorlabs.MotionControl.DeviceManagerCLI import DeviceManagerCLI
from Thorlabs.MotionControl.TCube.DCServoCLI import TCubeDCServo
from Thorlabs.MotionControl.KCube.DCServoCLI import KCubeDCServo
# adds .NET stuctures corresponding to primitives
from System import Decimal, Double

class TDC001_II(Instrument):
    '''
    Class to control the thorlabs TDC001 servo. Note that ALL DLL FUNCTIONS TAKING NUMERIC INPUT REQUIRE A SYSTEM.DECIMAL
    VALUE. Check help doc at C:\Program Files\Thorlabs\Kinesis\Thorlabs.MotionControl.DotNet_API for the DLL api.
    The class communicates with the device over USB.
    This is the working version. 5/16/2019

    --> last editted by ZQ for use in LISE 607 NV RT confocal setup on 3/12/2019
    '''

    _DEFAULT_SETTINGS = Parameter([
        Parameter('serial_number', 83860291, int, 'serial number written on device'),
        Parameter('position', 0.0, float, 'servo position (from 0 to 25 in mm)'),
        Parameter('velocity', 0.5, float, 'servo maximum velocity in mm/s'),
        Parameter('max_allowed_velocity', 2.5, float, 'servo maximum allowed velocity in mm/s'),
        Parameter('lower_limit', 0.0, float, 'servo maximum position in mm'),
        Parameter('upper_limit', 25.0, float, 'servo maximum position in mm'),

    ])

    startup_magnet_Z = float(0)
    # TDC001_lower_lim = float(0)
    # TDC001_upper_lim = float(25)  # for 25mm translation stage Z825B

    def __init__(self, name = None, settings = None):

        super(TDC001_II, self).__init__(name, settings)
        try:
            DeviceManagerCLI.BuildDeviceList()
            serial_number_list = DeviceManagerCLI.GetDeviceList(TCubeDCServo.DevicePrefix)
            print('serial_number_list:', str(serial_number_list))
        except (Exception):
            print("Exception raised by BuildDeviceList")
        if not (str(self.settings['serial_number']) in serial_number_list):
            print(str(self.settings['serial_number']) + " is NOT a valid serial number")
            raise
        else:
            print(str(self.settings['serial_number']) + " is a valid serial number")

        self.device = TCubeDCServo.CreateTCubeDCServo(str(self.settings['serial_number']))

        if(self.device == None):
            print(str(self.settings['serial_number']) + ' is not a TCubeDCServo')
            raise

        try:
            # self.device.Connect(str(self.settings['serial_number']))
            self._connect()
            # print(self.is_connected())
            # if self.is_connected():
            #     print('connection is successful')
            # else:
            #     print('failed to connect')
        except Exception:
            print('Failed to open device ' + str(self.settings['serial_number']))
            raise

        if not self.device.IsSettingsInitialized():
            try:
                self.device.WaitForSettingsInitialized(5000)
            except Exception:
                print("Settings failed to initialize")
                raise

        self.device.StartPolling(250)

        motorSettings = self.device.GetMotorConfiguration(str(self.settings['serial_number']))
        currentDeviceSettings = self.device.MotorDeviceSettings
        self.startup_position = self._get_position()
        self.startup_velocity = self._get_velocity()
        # self.TDC001_lower_lim = self.settings['lower_limit']
        # self.TDC001_upper_lim = self.settings['upper_limit']
        print('Connection to ' + str(self.settings['serial_number']) + ' is successful!')
        print('TDC001 startup position = ' + str(self.startup_position) + 'mm')
        print('TDC001 startup velocity = ' + str(self.startup_velocity) + 'mm/s')

    def _connect(self, verbose=False):

        print('Now connecting to ' + str(self.settings['serial_number']) )
        self.device.Connect(str(self.settings['serial_number']))

    def update(self, settings):
        '''
        Updates internal settings, as well as the position and velocity set on the physical device
        Args:
            settings: A dictionary in the form of settings as seen in default settings
        '''
        super(TDC001_II, self).update(settings)
        if self._settings_initialized: # this is super super important!!
            for key, value in settings.items():
            # for key, value in settings.iteritems():
                if key == 'position':
                    self._move_servo(value)
                elif key == 'velocity':
                    self._set_velocity(value)

    @property
    def _PROBES(self):
        return{
            'position': 'servo position in mm',
            'velocity': 'servo velocity in mm/s',
            'device': 'TDC001 device',
            'serial_number': 'serial number of device'

        }

    def read_probes(self, key):
        assert key in self._PROBES.keys()
        assert isinstance(key, str)

        #query always returns string, need to cast to proper return type
        if key in ['position']:
            return self._get_position()
        elif key in ['velocity']:
            return self._get_velocity()

    @property
    def is_connected(self):
        DeviceManagerCLI.BuildDeviceList()
        return(str(self.settings['serial_number']) in DeviceManagerCLI.GetDeviceList(TCubeDCServo.DevicePrefix))

    def __del__(self):
        '''
        Cleans up TDC001 connection
        :PostState: TDC001 is disconnected
        '''
        print('TDC001 is closing...')
        self.device.StopPolling()
        self.device.Disconnect()
        print('TDC001 is now closed')

    def goto_home(self):
        '''
        Recenters device at the home position. Raises an exception on failure.
        '''
        try:
            self.device.Home(60000)
        except Exception:
            print("Failed to move to position")
            raise

    def _move_servo(self, position, velocity = 0, verbose=True):
        '''
        Move servo to given position with given maximum velocity. Raises an exception on failure.
        :param position: position in mm, ranges from 0-25
        :param velocity: maximum velocity in mm/s, ranges from 0-2.5
        :PostState: servo has moved
        '''
        if (position >= self.settings['lower_limit']) & (position <= self.settings['upper_limit']):
            try:
                if(velocity != 0):
                    self._set_velocity(velocity)
                current_velocity = self.device.GetVelocityParams()
                old_position = self._get_position()
                if verbose:
                    print('---> Current position is {0} mm. Current velocity is {1} mm/s.'.format(old_position, current_velocity.MaxVelocity))
                    print('     Device is moving to indicated position = {0} mm... (be patient)'.format(position))

                self.device.MoveTo(self._Py_Decimal(position), 60000)

                # print(dir(self.device))
                current_position = self._get_position()
                print('     ---> Done: device current position is {0} mm.'.format(current_position))
                return True

            except Exception:
                print("FAILED to move to the position " + str(position) + " mm")
                return False
                # raise
        else:
            print(str(position) + " is out of allowed range [" + str(self.settings['lower_limit'])+","+str(self.settings['upper_limit'])+"]; TDC001 did not move!")
            return False




        # if verbose:
        #     current_position = self._get_position()
        #     print('TDC001 current position = ' + str(current_position))

        # if verbose:
        #     self._get_position()
            # print('Device now at position {0} mm'.format(position_now))

    def _get_position(self):
        '''
        :return: position of servo
        '''
        # print('Device now at position {0} mm'.format(self.device.Position()))
        # print(self.device.Position)
        return self._Undo_Decimal(self.device.Position)

    def _set_velocity(self, velocity, verbose = True):
        '''
        :param maximum velocity in mm/s, ranges from 0-2.5
        :PostState: velocity changed in hardware
        '''
        if velocity > self.settings['max_allowed_velocity']:
            print('Maximum allowed velocity is ' + str(self.settings['max_allowed_velocity']) + 'mm/s; No action!')
        elif(velocity != 0):
            # print(dir(self.device))
            velPars = self.device.GetVelocityParams()
            if verbose:
                print('---> Current velocity = {0} mm/s.'.format(velPars.MaxVelocity))
                print('     Setting to indicated velocity = {0} mm/s...'.format(velocity))
            velPars.MaxVelocity = self._Py_Decimal(velocity)
            self.device.SetVelocityParams(velPars)

                # print('---> Setting to indicated velocity = {:0.3f} mm/s.'.format(float(velocity)))
            velPars_new = self.device.GetVelocityParams()
            if verbose:
                print('     ---> Done: device current velocity is {0} mm/s.'.format(velPars_new.MaxVelocity))
                # print('     ---> Device current velocity is {:0.3f} mm/s.'.format(float(velPars_new.MaxVelocity)))

    def _get_velocity(self):
        '''
        :return: maximum velocity setting
        '''
        return self._Undo_Decimal(self.device.GetVelocityParams().MaxVelocity)

    def _Py_Decimal(self, value):
        '''
        Casting a python double to System.Decimal results in the Decimal having only integer values, likely due to an
        improper selection of the overloaded Decimal function. Casting it first to System.Double, which always maintains
        precision, then from Double to Decimal, where the proper overloaded function is clear, bypasses this issue
        :param value: a python double
        :return: the input as a System.Decimal
        '''
        return Decimal(Double(value))

    def _Undo_Decimal(self, value):
        '''
        Casting back from System.Decimal to a python float fails due to overloading issues, but one can successfully
        cast back to a string. Thus, we use a two-part cast to return to python numeric types
        :param value: a System.Decimal
        :return: the input as a python float
        '''
        return float(str(value))

class KDC101_II(Instrument):
    '''
    Class to control the thorlabs KDC101 servo. Note that ALL DLL FUNCTIONS TAKING NUMERIC INPUT REQUIRE A SYSTEM.DECIMAL
    VALUE. Check help doc at C:\Program Files\Thorlabs\Kinesis\Thorlabs.MotionControl.DotNet_API for the DLL api.
    The class communicates with the device over USB.

    -- Ziwei Qiu 5/16/2019
    '''

    # startup_magnet_azimuth = 0

    _DEFAULT_SETTINGS = Parameter([
        Parameter('serial_number', 27500496, int, 'serial number written on device'),
        Parameter('angle', 0.0, float, 'servo angle (from 0 to 360 deg)'),
        Parameter('angular_velocity', 5.0, float, 'servo angular velocity in deg/s'),
        Parameter('max_allowed_velocity', 5.0, float, 'maximum allowed angular velocity in deg/s'),
        # Parameter('angular_acceleration', 10, float, 'servo angular acceleration in deg/s/s'),
        Parameter('lower_limit', 0.0, float, 'servo maximum velocity in mm/s'),
        Parameter('upper_limit', 360.0, float, 'servo maximum velocity in mm/s'),
    ])

    def __init__(self, name = None, settings = None):
        super(KDC101_II, self).__init__(name, settings)
        try:
            DeviceManagerCLI.BuildDeviceList()
            serial_number_list = DeviceManagerCLI.GetDeviceList(KCubeDCServo.DevicePrefix)
            # print('serial_number_list:', serial_number_list)
        except (Exception):
            print("Exception raised by BuildDeviceList")
        if not (str(self.settings['serial_number']) in serial_number_list):
            print(str(self.settings['serial_number']) + " is not a valid serial number")
            raise
        else:
            print(str(self.settings['serial_number']) + " is a valid serial number")

        self.device = KCubeDCServo.CreateKCubeDCServo(str(self.settings['serial_number']))

        if(self.device == None):
            print(str(self.settings['serial_number']) + " is NOT a KCubeDCServo")
            raise
        else:
            print(str(self.settings['serial_number']) + " is a KCubeDCServo")


        try:
            self._connect()
            #
            #
            # if self.is_connected():
            #     print('connection to ' + self.settings['serial_number'] + ' is successful')
            # else:
            #     print('failed to connect to '+ self.settings['serial_number'])
        except Exception:
            print('Failed to open device: ' + str(self.settings['serial_number']))
            raise

        if not self.device.IsSettingsInitialized():
            try:
                self.device.WaitForSettingsInitialized(5000)
            except Exception:
                print("Settings failed to initialize")
                raise

        self.device.StartPolling(250)

        motorSettings = self.device.GetMotorConfiguration(str(self.settings['serial_number']))
        currentDeviceSettings = self.device.MotorDeviceSettings
        self.startup_position = self._get_position()
        self.startup_velocity = self._get_velocity()
        print('Connection to ' + str(self.settings['serial_number']) + ' is successful!')
        print('KDC001 startup position = ' + str(self.startup_position) + 'deg')
        print('KDC001 startup velocity = ' + str(self.startup_velocity) + 'deg/s')

    def _connect(self, verbose=False):

        print('Now connecting to ' + str(self.settings['serial_number']))
        self.device.Connect(str(self.settings['serial_number']))

    def update(self, settings):
        '''
        Updates internal settings, as well as the position and velocity set on the physical device
        Args:
            settings: A dictionary in the form of settings as seen in default settings
        '''
        super(KDC101_II, self).update(settings)
        if self._settings_initialized:  # this is super super important!!
            for key, value in settings.items():
            # for key, value in settings.iteritems():
                if key == 'angle':
                    self._move_servo(value)
                elif key == 'angular_velocity':
                    self._set_velocity(value)

    @property
    def _PROBES(self):
        return{
            'position': 'servo position in deg',
            'velocity': 'servo velocity in deg/s',
            'device': 'KDC101 device',
            'serial_number': 'serial number of device'

        }

    def read_probes(self, key):
        assert key in self._PROBES.keys()
        assert isinstance(key, str)

        #query always returns string, need to cast to proper return type
        if key in ['position']:
            return self._get_position()
        elif key in ['velocity']:
            return self._get_velocity()

    @property
    def is_connected(self):
        DeviceManagerCLI.BuildDeviceList()
        return(str(self.settings['serial_number']) in DeviceManagerCLI.GetDeviceList(KCubeDCServo.DevicePrefix))

    def __del__(self):
        '''
        Cleans up TDC001 connection
        :PostState: TDC001 is disconnected
        '''

        print('KDC101 is closing...')
        self.device.StopPolling()
        self.device.Disconnect()
        print('KDC101 is now closed')

    def goto_home(self):
        '''
        Recenters device at the home position. Raises an exception on failure.
        '''
        try:
            self.device.Home(60000)
        except Exception:
            print("Failed to move to position")
            raise

    def _move_servo(self, position, velocity = 0, verbose=True):
        '''
        Move servo to given position with given maximum velocity. Raises an exception on failure.
        :param position: position in deg, ranges from 0-360
        :param velocity: maximum velocity in deg/s, ranges from 0-2.5
        :PostState: servo has moved
        '''
        # if True:
        #     try:
        #         if(velocity != 0):
        #             self._set_velocity(velocity)
        #         self.device.MoveTo(self._Py_Decimal(position), 60000)
        #         print('Moved KDC101 to ' + str(position) + 'deg')
        #     except Exception:
        #         print('FAILED to move magnet azimuth to ' + str(position) + 'deg')
        #         # raise
        # else:
        #     print(str(position) + ' is out of allowed range [' + str(self.magnetZ_lower_lim) + ',' + str(
        #         self.magnetZ_upper_lim) + ']; no magnet azimuth move initiated')



        if (position >= self.settings['lower_limit']) & (position <= self.settings['upper_limit']):
            try:
                if(velocity != 0):
                    self._set_velocity(velocity)
                current_velocity = self.device.GetVelocityParams()
                old_position = self._get_position()
                if verbose:
                    print('---> Current position is {0} deg. Current velocity is {1} deg/s.'.format(old_position, current_velocity.MaxVelocity))
                    print('     Device is moving to indicated position = {0} deg... (be patient)'.format(position))

                self.device.MoveTo(self._Py_Decimal(position), 60000)

                # print(dir(self.device))
                current_position = self._get_position()
                print('     ---> Done: device current position is {0} deg.'.format(current_position))
                return True

            except Exception:
                print("FAILED to move to the position " + str(position) + " mm")
                return False
                # raise
        else:
            print(str(position) + " is out of allowed range [" + str(self.settings['lower_limit'])+","+str(self.settings['upper_limit'])+"]; TDC001 did not move!")
            return False

    def _get_position(self):
        '''
        :return: position of servo
        '''
        # print('_get_position')
        return self._Undo_Decimal(self.device.Position)

    def _set_velocity(self, velocity, verbose = True):
        '''
        :param maximum velocity in mm/s, ranges from 0-2.5
        :PostState: velocity changed in hardware
        '''
        # if(velocity != 0):
        #     velPars = self.device.GetVelocityParams()
        #     velPars.MaxVelocity = self._Py_Decimal(velocity)
        #     self.device.SetVelocityParams(velPars)

        if velocity > self.settings['max_allowed_velocity']:
            print('Maximum allowed velocity is ' + str(self.settings['max_allowed_velocity']) + 'deg/s; No action!')
        elif(velocity != 0):
            # print(dir(self.device))
            velPars = self.device.GetVelocityParams()
            if verbose:
                print('---> Current velocity = {0} deg/s.'.format(velPars.MaxVelocity))
                print('     Setting to indicated velocity = {0} deg/s...'.format(velocity))
            velPars.MaxVelocity = self._Py_Decimal(velocity)
            self.device.SetVelocityParams(velPars)

                # print('---> Setting to indicated velocity = {:0.3f} mm/s.'.format(float(velocity)))
            velPars_new = self.device.GetVelocityParams()
            if verbose:
                print('     ---> Done: device current velocity is {0} deg/s.'.format(velPars_new.MaxVelocity))
                # print('     ---> Device current velocity is {:0.3f} mm/s.'.format(float(velPars_new.MaxVelocity)))

    def _get_velocity(self):
        '''
        :return: maximum velocity setting
        '''
        # print('_get_velocity')
        return self._Undo_Decimal(self.device.GetVelocityParams().MaxVelocity)

    def _Py_Decimal(self, value):
        '''
        Casting a python double to System.Decimal results in the Decimal having only integer values, likely due to an
        improper selection of the overloaded Decimal function. Casting it first to System.Double, which always maintains
        precision, then from Double to Decimal, where the proper overloaded function is clear, bypasses this issue
        :param value: a python double
        :return: the input as a System.Decimal
        '''
        return Decimal(Double(value))

    def _Undo_Decimal(self, value):
        '''
        Casting back from System.Decimal to a python float fails due to overloading issues, but one can successfully
        cast back to a string. Thus, we use a two-part cast to return to python numeric types
        :param value: a System.Decimal
        :return: the input as a python float
        '''
        return float(str(value))

class MagnetX(TDC001_II):
    '''

    Same as TDC001 except the serial number is fixed
    --> last editted by ZQ for use in LISE 607 NV RT confocal setup on 3/12/2019

    '''
    _DEFAULT_SETTINGS = Parameter([
        Parameter('serial_number', 83860546, [83860546], 'serial number written on device'),
        Parameter('position', 0.0, float, 'servo position (from 0 to 25 in mm)'),
        Parameter('velocity', 0.5, float, 'servo maximum velocity in mm/s'),
        Parameter('max_allowed_velocity', 2.5, float, 'servo maximum allowed velocity in mm/s'),
        Parameter('lower_limit', 0, float, 'servo maximum velocity in mm/s'),
        Parameter('upper_limit', 25, float, 'servo maximum velocity in mm/s'),

    ])

    _PROBES = {'position': 'current position of stage', 'velocity':'current velocity of stage', 'device': 'TDC001 device', 'serial_number': 'serial number of device'}

    def __init__(self, name=None, settings=None):
        super(MagnetX, self).__init__()

class MagnetY(TDC001_II):
    '''

    Same as TDC001 except the serial number is fixed
    --> last editted by ZQ for use in LISE 607 NV RT confocal setup on 3/12/2019

    '''
    _DEFAULT_SETTINGS = Parameter([
        Parameter('serial_number', 83860305, [83860305], 'serial number written on device'),
        Parameter('position', 0.0, float, 'servo position (from 0 to 25 in mm)'),
        Parameter('velocity', 0.5, float, 'servo maximum velocity in mm/s'),
        Parameter('max_allowed_velocity', 2.5, float, 'servo maximum allowed velocity in mm/s'),
        Parameter('lower_limit', 0, float, 'servo maximum velocity in mm/s'),
        Parameter('upper_limit', 25, float, 'servo maximum velocity in mm/s'),

    ])

    _PROBES = {'position': 'current position of stage', 'velocity':'current velocity of stage', 'device': 'TDC001 device', 'serial_number': 'serial number of device'}

    def __init__(self, name=None, settings=None):
        super(MagnetY, self).__init__()

class MagnetZ(TDC001_II):
    '''

    Same as TDC001 except the serial number is fixed
    --> last editted by ZQ for use in LISE 607 NV RT confocal setup on 3/12/2019

    '''
    _DEFAULT_SETTINGS = Parameter([
        Parameter('serial_number', 83860291, [83860291], 'serial number written on device'),
        Parameter('position', 0.0, float, 'servo position (from 0 to 25 in mm)'),
        Parameter('velocity', 0.5, float, 'servo maximum velocity in mm/s'),
        Parameter('max_allowed_velocity', 2.5, float, 'servo maximum allowed velocity in mm/s'),
        Parameter('lower_limit', 0, float, 'servo maximum velocity in mm/s'),
        Parameter('upper_limit', 25, float, 'servo maximum velocity in mm/s'),

    ])

    _PROBES = {'position': 'current position of stage', 'velocity':'current velocity of stage', 'device': 'TDC001 device', 'serial_number': 'serial number of device'}

    def __init__(self, name=None, settings=None):
        super(MagnetZ, self).__init__()

class SampleX(TDC001_II):
    '''

    Same as TDC001 except the serial number is fixed
    --> last editted by ZQ for use in LISE 607 NV RT confocal setup on 3/12/2019

    '''
    _DEFAULT_SETTINGS = Parameter([
        Parameter('serial_number', 83860519, [83860519], 'serial number written on device'),
        Parameter('position', 0.0, float, 'servo position (from 0 to 25 in mm)'),
        Parameter('velocity', 0.5, float, 'servo maximum velocity in mm/s'),
        Parameter('max_allowed_velocity', 2.5, float, 'servo maximum allowed velocity in mm/s'),
        Parameter('lower_limit', 0, float, 'servo maximum velocity in mm/s'),
        Parameter('upper_limit', 12, float, 'servo maximum velocity in mm/s'),

    ])

    _PROBES = {'position': 'current position of stage', 'velocity':'current velocity of stage', 'device': 'TDC001 device', 'serial_number': 'serial number of device'}

    def __init__(self, name=None, settings=None):
        super(SampleX, self).__init__()

class SampleY(TDC001_II):
    '''

    Same as TDC001 except the serial number is fixed
    --> last editted by ZQ for use in LISE 607 NV RT confocal setup on 3/12/2019

    '''
    _DEFAULT_SETTINGS = Parameter([
        Parameter('serial_number', 83860538, [83860538], 'serial number written on device'),
        Parameter('position', 0.0, float, 'servo position (from 0 to 25 in mm)'),
        Parameter('velocity', 0.5, float, 'servo maximum velocity in mm/s'),
        Parameter('max_allowed_velocity', 2.5, float, 'servo maximum allowed velocity in mm/s'),
        Parameter('lower_limit', 0, float, 'servo maximum velocity in mm/s'),
        Parameter('upper_limit', 12, float, 'servo maximum velocity in mm/s'),

    ])

    _PROBES = {'position': 'current position of stage', 'velocity':'current velocity of stage', 'device': 'TDC001 device', 'serial_number': 'serial number of device'}

    def __init__(self, name=None, settings=None):
        super(SampleY, self).__init__()

class SampleZ(TDC001_II):
    '''

    Same as TDC001 except the serial number is fixed
    --> last editted by ZQ for use in LISE 607 NV RT confocal setup on 3/12/2019

    '''
    _DEFAULT_SETTINGS = Parameter([
        Parameter('serial_number', 83854376, [83854376], 'serial number written on device'),
        Parameter('position', 0.0, float, 'servo position (from 0 to 25 in mm)'),
        Parameter('velocity', 0.5, float, 'servo maximum velocity in mm/s'),
        Parameter('max_allowed_velocity', 2.5, float, 'servo maximum allowed velocity in mm/s'),
        Parameter('lower_limit', 0, float, 'servo maximum velocity in mm/s'),
        Parameter('upper_limit', 12, float, 'servo maximum velocity in mm/s'),

    ])

    _PROBES = {'position': 'current position of stage', 'velocity':'current velocity of stage', 'device': 'TDC001 device', 'serial_number': 'serial number of device'}

    def __init__(self, name=None, settings=None):
        super(SampleZ, self).__init__()

class IntensityWheel(KDC101_II):
    '''

    Same as KDC001 except the serial number is fixed
    --> last editted by ZQ for use in LISE 607 NV RT confocal setup on 5/16/2019

    '''
    _DEFAULT_SETTINGS = Parameter([
        Parameter('serial_number', 27500496, [27500496], 'serial number written on device'),
        Parameter('angle', 100.0, float, 'servo angle (from 0 to 360 deg)'),
        Parameter('angular_velocity', 5.0, float, 'servo angular velocity in deg/s'),
        Parameter('max_allowed_velocity', 10.0, float, 'maximum allowed angular velocity in deg/s'),
        Parameter('lower_limit', 0.0, float, 'servo maximum velocity in mm/s'),
        Parameter('upper_limit', 360.0, float, 'servo maximum velocity in mm/s'),
    ])

    _PROBES = {'position': 'current position of stage', 'velocity':'current velocity of stage', 'device': 'KDC001 device', 'serial_number': 'serial number of device'}

    def __init__(self, name=None, settings=None):
        super(IntensityWheel, self).__init__()


if __name__ == '__main__':
    #A test function for the device. Tries to connect to the
    a = TDC001_II()
    a.is_connected