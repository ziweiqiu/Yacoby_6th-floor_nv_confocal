from pylabcontrol.core import Parameter, Instrument
import ctypes
import time

class TLI_DeviceInfo(ctypes.Structure):
    _fields_ = [("typeID", ctypes.c_ulong),
                ("description", (65 * ctypes.c_char)),
                ("serialNo", (9 * ctypes.c_char)),
                ("PID", ctypes.c_ulong),
                ("isKnownType", ctypes.c_bool),
                ("motorType", ctypes.c_int),
                ("isPiezoDevice", ctypes.c_bool),
                ("isLaser", ctypes.c_bool),
                ("isCustomType", ctypes.c_bool),
                ("isRack", ctypes.c_bool),
                ("maxChannels", ctypes.c_short)]

class TDC001(Instrument):
    """
    Class to control the thorlabs TDC001 servo. Note that ALL DLL FUNCTIONS TAKING NUMERIC INPUT REQUIRE A SYSTEM.DECIMAL
    VALUE. Check help doc at C:\Program Files\Thorlabs\Kinesis\Thorlabs.MotionControl.DotNet_API for the DLL api.
    The class communicates with the device over USB.

    --> last editted by ZQ for use in LISE 607 NV RT confocal setup on 2/18/2019 3:37pm
    --> modified from B26 class
    """

    _DEFAULT_SETTINGS = Parameter([
        Parameter('serial_number', 83860546, [83860546], 'serial number written on device'),
        Parameter('position', 3.0, float, 'servo position (0 to 25) [mm]'),
        Parameter('velocity', 0.01, float, 'servo velocity (0 to 2.6) [mm/s]. If set to zero, instrument default will be used')
    ])

    _PROBES = {}

    def __init__(self, name=None, settings=None):

        super(TDC001, self).__init__(name, settings)
        # self._servo_library = ctypes.cdll.LoadLibrary('C:\\Program Files\\Thorlabs\\Kinesis\\Thorlabs.MotionControl.KCube.DCServo.dll')
        self._servo_library = ctypes.cdll.LoadLibrary(
            'C:\\Program Files\\Thorlabs\\Kinesis\\Thorlabs.MotionControl.TCube.DCServo.dll')

        # conversion from mm to device units at
        # https://www.thorlabs.com/software/apt/APT_Communications_Protocol_Rev_15.pdf, page 16
        self._position_encoder2mm_conversion_factor = 34304
        self._velocity_encoder2mm_conversion_factor = 767367.49
        self._acceleration_encoder2mm_conversion_factor = 261.93
        self._acceleration = 1 * self._acceleration_encoder2mm_conversion_factor # use hard-coded acceleration of 1 mm/s^2 (can be changed to parameter if we want)
        self._manually_set_library_inputs_and_outputs()
        self._connect()

    def _manually_set_library_inputs_and_outputs(self):
        """
        Sets the input and output types for each servo library call we make.

        """
        self._servo_library.TLI_BuildDeviceList.restypes = ctypes.c_short

        self._servo_library.TLI_GetDeviceListSize.restypes = ctypes.c_short

        self._servo_library.TLI_GetDeviceListByTypeExt.argtypes = [ctypes.POINTER(ctypes.c_char), ctypes.c_ulong,
                                                                   ctypes.c_int]
        self._servo_library.TLI_GetDeviceListByTypeExt.restypes = ctypes.c_short

        self._servo_library.TLI_GetDeviceInfo.argtypes = [ctypes.POINTER(ctypes.c_char), ctypes.POINTER(TLI_DeviceInfo)]
        self._servo_library.TLI_GetDeviceInfo.restypes = ctypes.c_short

        self._servo_library.CC_StartPolling.argtypes = [ctypes.POINTER(ctypes.c_char), ctypes.c_int]
        self._servo_library.CC_StartPolling.restypes = ctypes.c_bool

        self._servo_library.CC_Open.argtypes = [ctypes.POINTER(ctypes.c_char)]
        self._servo_library.CC_Open.restypes = ctypes.c_short

        self._servo_library.CC_Close.argtypes = [ctypes.POINTER(ctypes.c_char)]
        self._servo_library.CC_Close.restypes = ctypes.c_short

        self._servo_library.CC_ClearMessageQueue.argtypes = [ctypes.POINTER(ctypes.c_char)]
        self._servo_library.CC_ClearMessageQueue.restypes = ctypes.c_short

        self._servo_library.CC_WaitForMessage.argtypes = [ctypes.POINTER(ctypes.c_char),
                                                          ctypes.POINTER(ctypes.c_ushort),
                                                          ctypes.POINTER(ctypes.c_ushort),
                                                          ctypes.POINTER(ctypes.c_ulong)]
        self._servo_library.CC_WaitForMessage.restypes = ctypes.c_bool

        self._servo_library.CC_GetVelParams.argtypes = [ctypes.POINTER(ctypes.c_char),
                                                        ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long)]

        self._servo_library.CC_GetVelParams.restypes = ctypes.c_short

        self._servo_library.CC_MoveToPosition.argtypes = [ctypes.POINTER(ctypes.c_char), ctypes.c_int]
        self._servo_library.CC_MoveToPosition.restypes = ctypes.c_short

        self._servo_library.CC_SetVelParams.argtypes = [ctypes.POINTER(ctypes.c_char), ctypes.c_int, ctypes.c_int]
        self._servo_library.CC_SetVelParams.restypes = ctypes.c_short

        self._servo_library.CC_StopPolling.argtypes = [ctypes.POINTER(ctypes.c_char)]

        self._servo_library.CC_Home.argtypes = [ctypes.POINTER(ctypes.c_char)]
        self._servo_library.CC_Home.restypes = ctypes.c_short

        self._servo_library.CC_GetPosition.argtypes = [ctypes.POINTER(ctypes.c_char)]
        self._servo_library.CC_GetPosition.restypes = ctypes.c_int

    def _connect(self, verbose=True):

        """
        makes connection to the device through the serial port. Also opens connection to the device.
        """

        # this version of the serial number is useful
        self._serial_num_int = self.settings['serial_number']

        self._serial_num = ctypes.c_char_p(bytes(str(self.settings['serial_number']), "utf-8"))

        if self._servo_library.TLI_BuildDeviceList() == 0:
            num_devices = self._servo_library.TLI_GetDeviceListSize()
            if verbose:
                print('Number of devices detected: {0}'.format(num_devices))

            # The servo library fills in a byte string of connected device serial numbers, separated by a comma.
            # The length of this string will be the length of the serial number (8 bytes), plus a byte for a comma,
            # for each connected device.
            # Here, we first pre-allocate the string, then send it to the library function, and then examine it.

            connected_devices_string_length = num_devices*9
            string_of_serial_nums = ctypes.create_string_buffer(connected_devices_string_length)

            # self._servo_library.TLI_GetDeviceListByTypeExt(string_of_serial_nums,
            #                                                ctypes.c_ulong(connected_devices_string_length+1),
            #                                                ctypes.c_int(27)) # KDC101 serial number starts with 27
            self._servo_library.TLI_GetDeviceListByTypeExt(string_of_serial_nums,
                                                           ctypes.c_ulong(connected_devices_string_length + 1),
                                                           ctypes.c_int(83)) # TDC001 serial number starts with 83



            list_of_connected_serial_nums = string_of_serial_nums.raw.decode("utf-8")
            # list_of_connected_serial_nums = string_of_serial_nums.raw.decode('ascii')
            print('print(list_of_connected_serial_nums)')
            print(list_of_connected_serial_nums)

            if str(self.settings['serial_number']) not in list_of_connected_serial_nums:
                error_msg = 'No servo with the given serial number was detected.\n'
                error_msg += 'Given Serial Number: {0}\n'.format(self.settings['serial_number'])
                if list_of_connected_serial_nums:
                    error_msg += 'Connected Devices: {0}\n'.format(list_of_connected_serial_nums)
                else:
                    error_msg += 'Connected Devices: None\n'
                raise AttributeError(error_msg)

            elif verbose:
                print('Found device with matching serial number.')
                device_info = TLI_DeviceInfo()
                self._servo_library.TLI_GetDeviceInfo(self._serial_num, ctypes.byref(device_info))
                print("Description: ", device_info.description)
                print("Serial No: ", device_info.serialNo)
                print("USB PID: ", device_info.PID)

        if self.settings['velocity'] > 0.0: # update the velocity if not set to default
            self.set_velocity()
            # current_velocity = self.get_velocity()
            # print('velocity is set to ({0} mm/s)'.format(current_velocity))

        # open connection to the device. The other option is to NOT call _open_device, and require the script to handle device opening and closing
        self._open_device()

        # get the current position of the device and set the position setting to this value. That makes sure that the device does not move when you create the instrument.
        self.settings['position'] = self.get_position()
        # todo(emma): implement self.get_velocity()
        #self.velocity = self.get_velocity()

    def set_position(self, verbose=True):

        """
        sets position of device in mm
        """
        position = ctypes.c_int(int(self._position_encoder2mm_conversion_factor * self.settings['position']))
        self._servo_library.CC_MoveToPosition(self._serial_num, position) # command sent to device
        message_type = ctypes.c_ushort(0)
        message_id = ctypes.c_ushort(0)
        message_data = ctypes.c_ulong(0)
        if verbose:
            print('Device is moving to indicated position ({0} mm)'.format(self.settings['position']))

        # wait for it to actually move
        self._servo_library.CC_WaitForMessage(self._serial_num, message_type, message_id, message_data)
        while message_type.value != 2 or message_id.value != 1:
            self._servo_library.CC_WaitForMessage(self._serial_num, message_type, message_id, message_data)

        if verbose:
            position = self.get_position()
            print('Device now at position {0} mm'.format(position))

    def home(self, verbose=True):

        """
        Calibrates origin: scans full range and sets stage position to zero
        """

        self._servo_library.CC_Home(self._serial_num) # command sent to device!
        if verbose:
            print('Device is homing.')
        # wait for it to actually home
        message_type = ctypes.c_ushort(0)
        message_id = ctypes.c_ushort(0)
        message_data = ctypes.c_ulong(0)
        self._servo_library.CC_WaitForMessage(self._serial_num, message_type, message_id, message_data)
        while message_type.value != 2 or message_id.value != 0:
            self._servo_library.CC_WaitForMessage(self._serial_num, message_type, message_id, message_data)

        if verbose:
            print('Device finished homing')

    def _open_device(self, verbose=True):

        """
        Opens the serial connection to the device, which lets you send commands to and receive data from the device.
        Also starts polling all of the devices connected to the computer.
        """
        # print('_open_device now')

        if not self._servo_library.CC_Open(
                self._serial_num):
            if verbose:
                print('Starting to poll')
            milliseconds = ctypes.c_int(200)
            self._servo_library.CC_StartPolling(self._serial_num,
                                                milliseconds)  # continuous polling allows you to keep track of what the instrument is doing in real time
            time.sleep(3)
            self._servo_library.CC_ClearMessageQueue(self._serial_num)

    def _close_device(self, verbose=True):

        """
        Stops polling and closes the serial connection.
        """
        if verbose:
            print('TDC001 is closing')
        self._servo_library.CC_StopPolling(self._serial_num)
        self._servo_library.CC_Close(self._serial_num)
        time.sleep(3)
        self._servo_library.CC_ClearMessageQueue(self._serial_num)

    def get_position(self, verbose=True):

        """
        returns position of stage in mm
        """
        print('getting current position now')
        position = self._servo_library.CC_GetPosition(
            self._serial_num) / self._position_encoder2mm_conversion_factor
        if verbose:
            print('Position of device is currently {0} mm'.format(position))
        return position

    def get_velocity(self, verbose=True):

        max_vel = 2.6

        """
        returns velocity (when moving) of stage in mm/s
        """
       # raise NotImplementedError

        ## todo: fix the input arguments of get_velocity ER 20180519
        print('getting current velocity now')
        velocity_pointer_type = ctypes.POINTER(ctypes.c_long) #ctypes.POINTER(ctypes.c_long)#ctypes.LP_c_long()#ctypes.POINTER(ctypes.c_int)
        acceleration_pointer_type = ctypes.POINTER(ctypes.c_long) #ctypes.POINTER(ctypes.c_long)#ctypes.LP_c_long()#ctypes.POINTER(ctypes.c_int)
        velocity_pointer = ctypes.byref(ctypes.c_long())#velocity_pointer_type()
        # velocity_pointer = ctypes.byref(ctypes.c_long(max_vel*self._velocity_encoder2mm_conversion_factor))  # velocity_pointer_type()

        acceleration_pointer = ctypes.byref(ctypes.c_long())#acceleration_pointer_type()
        # vel = self._servo_library.CC_GetVelParams(self._serial_num, velocity_pointer, acceleration_pointer)
        vel = self._servo_library.CC_GetVelParams(self._serial_num, acceleration_pointer, velocity_pointer)

        if verbose:
            print(vel)
            print(type(vel))
            print('Velocity of device is currently {0} mm/s'.format(vel))
        return vel

    def set_velocity(self, verbose=True):

        """
        sets velocity (when moving) of stage in mm/s
        """
        # maximum velocity of instrument
        max_vel = 2.6
        if self.settings['velocity'] > max_vel:
            self.settings['velocity'] = max_vel

        velocity = ctypes.c_int(int(self._velocity_encoder2mm_conversion_factor * self.settings['velocity']))
        acceleration = ctypes.c_int(int(self._acceleration))
        self._servo_library.CC_SetVelParams(self._serial_num, acceleration, velocity)
        print(velocity)
        print('Setting to indicated velocity: {0} mm/s'.format( velocity.value / self._velocity_encoder2mm_conversion_factor))
        if verbose:
            pass
            # current_vel = self.get_velocity()
            # print('Device velocity was set to {0} mm/s'.format(current_vel))
          # todo(emma): implement get_velocity()

    def update(self, settings, verbose=True):
        super(TDC001, self).update(settings)
        if self._settings_initialized:
            for key, value in settings.items():
                if key == 'position':
                    self.set_position()
                elif key == 'velocity':
                    self.set_velocity()

    def read_probes(self, key = None):
        assert key in list(self._PROBES.keys())
        assert isinstance(key, str)

        if key in ['position']:
            return (self.get_position(True))
        if key in ['serial_number']:
            return (self._serial_num_int)
        elif key in ['velocity']:
            # todo(emma): implement get_velocity
            pass

    def __del__(self):
        # when done with device we have to close the connections
        # called when GUI is closed
        self._close_device()
        print('thorlab TDC001 is closed.')

# the following instrument is NOT WORKING!!!
class KDC101(Instrument):
    """
    Class to control the thorlabs TDC001 servo. Note that ALL DLL FUNCTIONS TAKING NUMERIC INPUT REQUIRE A SYSTEM.DECIMAL
    VALUE. Check help doc at C:\Program Files\Thorlabs\Kinesis\Thorlabs.MotionControl.DotNet_API for the DLL api.
    The class communicates with the device over USB.

    NOT WORKING!!!

    --> last editted by ZQ for use in LISE 607 NV RT confocal setup on 2/18/2019 3:37pm
    --> modified from B26 class
    """

    _DEFAULT_SETTINGS = Parameter([
        Parameter('serial_number', 27500496, [27500496], 'serial number written on device'),
        Parameter('position', 3.0, float, 'servo position (0 to 360) [deg]'),
        Parameter('velocity', 1.0, float, 'servo velocity (0 to 2.5) [deg/s]. If set to zero, instrument default will be used')
    ])

    _PROBES = {}

    def __init__(self, name=None, settings=None):

        super(KDC101, self).__init__(name, settings)
        # self._servo_library = ctypes.cdll.LoadLibrary('C:\\Program Files\\Thorlabs\\Kinesis\\Thorlabs.MotionControl.KCube.DCServo.dll')
        self._servo_library = ctypes.cdll.LoadLibrary(
            'C:\\Program Files\\Thorlabs\\Kinesis\\Thorlabs.MotionControl.TCube.DCServo.dll')

        # conversion from mm to device units at
        # https://www.thorlabs.com/software/apt/APT_Communications_Protocol_Rev_15.pdf, page 16
        # self._position_encoder2mm_conversion_factor = 34304
        # self._velocity_encoder2mm_conversion_factor = 767367.49
        # self._acceleration_encoder2mm_conversion_factor = 261.93
        self._position_encoder2deg_conversion_factor = 1919.6418578623391
        self._velocity_encoder2deg_conversion_factor = 42941.66
        self._acceleration_encoder2deg_conversion_factor = 14.66
        self._acceleration = 1 * self._acceleration_encoder2deg_conversion_factor # use hard-coded acceleration of 1 mm/s^2 (can be changed to parameter if we want)
        self._manually_set_library_inputs_and_outputs()
        self._connect()

    def _manually_set_library_inputs_and_outputs(self):
        """
        Sets the input and output types for each servo library call we make.

        """
        self._servo_library.TLI_BuildDeviceList.restypes = ctypes.c_short

        self._servo_library.TLI_GetDeviceListSize.restypes = ctypes.c_short

        self._servo_library.TLI_GetDeviceListByTypeExt.argtypes = [ctypes.POINTER(ctypes.c_char), ctypes.c_ulong,
                                                                   ctypes.c_int]
        self._servo_library.TLI_GetDeviceListByTypeExt.restypes = ctypes.c_short

        self._servo_library.TLI_GetDeviceInfo.argtypes = [ctypes.POINTER(ctypes.c_char), ctypes.POINTER(TLI_DeviceInfo)]
        self._servo_library.TLI_GetDeviceInfo.restypes = ctypes.c_short

        self._servo_library.CC_StartPolling.argtypes = [ctypes.POINTER(ctypes.c_char), ctypes.c_int]
        self._servo_library.CC_StartPolling.restypes = ctypes.c_bool

        self._servo_library.CC_Open.argtypes = [ctypes.POINTER(ctypes.c_char)]
        self._servo_library.CC_Open.restypes = ctypes.c_short

        self._servo_library.CC_Close.argtypes = [ctypes.POINTER(ctypes.c_char)]
        self._servo_library.CC_Close.restypes = ctypes.c_short

        self._servo_library.CC_ClearMessageQueue.argtypes = [ctypes.POINTER(ctypes.c_char)]
        self._servo_library.CC_ClearMessageQueue.restypes = ctypes.c_short

        self._servo_library.CC_WaitForMessage.argtypes = [ctypes.POINTER(ctypes.c_char),
                                                          ctypes.POINTER(ctypes.c_ushort),
                                                          ctypes.POINTER(ctypes.c_ushort),
                                                          ctypes.POINTER(ctypes.c_ulong)]
        self._servo_library.CC_WaitForMessage.restypes = ctypes.c_bool

        self._servo_library.CC_GetVelParams.argtypes = [ctypes.POINTER(ctypes.c_char),
                                                        ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long)]

        self._servo_library.CC_GetVelParams.restypes = ctypes.c_short

        self._servo_library.CC_MoveToPosition.argtypes = [ctypes.POINTER(ctypes.c_char), ctypes.c_int]
        self._servo_library.CC_MoveToPosition.restypes = ctypes.c_short

        self._servo_library.CC_SetVelParams.argtypes = [ctypes.POINTER(ctypes.c_char), ctypes.c_int, ctypes.c_int]
        self._servo_library.CC_SetVelParams.restypes = ctypes.c_short

        self._servo_library.CC_StopPolling.argtypes = [ctypes.POINTER(ctypes.c_char)]

        self._servo_library.CC_Home.argtypes = [ctypes.POINTER(ctypes.c_char)]
        self._servo_library.CC_Home.restypes = ctypes.c_short

        self._servo_library.CC_GetPosition.argtypes = [ctypes.POINTER(ctypes.c_char)]
        self._servo_library.CC_GetPosition.restypes = ctypes.c_int

    def _connect(self, verbose=True):

        """
        makes connection to the device through the serial port. Also opens connection to the device.
        """

        # this version of the serial number is useful
        self._serial_num_int = self.settings['serial_number']

        self._serial_num = ctypes.c_char_p(bytes(str(self.settings['serial_number']), "utf-8"))

        if self._servo_library.TLI_BuildDeviceList() == 0:
            num_devices = self._servo_library.TLI_GetDeviceListSize()
            if verbose:
                print('Number of devices detected: {0}'.format(num_devices))

            # The servo library fills in a byte string of connected device serial numbers, separated by a comma.
            # The length of this string will be the length of the serial number (8 bytes), plus a byte for a comma,
            # for each connected device.
            # Here, we first pre-allocate the string, then send it to the library function, and then examine it.

            connected_devices_string_length = num_devices*9
            string_of_serial_nums = ctypes.create_string_buffer(connected_devices_string_length)

            self._servo_library.TLI_GetDeviceListByTypeExt(string_of_serial_nums,
                                                           ctypes.c_ulong(connected_devices_string_length+1),
                                                           ctypes.c_int(27)) # KDC101 serial number starts with 27
            # self._servo_library.TLI_GetDeviceListByTypeExt(string_of_serial_nums,
            #                                                ctypes.c_ulong(connected_devices_string_length + 1),
            #                                                ctypes.c_int(83)) # TDC001 serial number starts with 83



            list_of_connected_serial_nums = string_of_serial_nums.raw.decode("utf-8")
            # list_of_connected_serial_nums = string_of_serial_nums.raw.decode('ascii')
            print('print(list_of_connected_serial_nums)')
            print(list_of_connected_serial_nums)

            if str(self.settings['serial_number']) not in list_of_connected_serial_nums:
                error_msg = 'No servo with the given serial number was detected.\n'
                error_msg += 'Given Serial Number: {0}\n'.format(self.settings['serial_number'])
                if list_of_connected_serial_nums:
                    error_msg += 'Connected Devices: {0}\n'.format(list_of_connected_serial_nums)
                else:
                    error_msg += 'Connected Devices: None\n'
                raise AttributeError(error_msg)

            elif verbose:
                print('Found device with matching serial number.')
                device_info = TLI_DeviceInfo()
                self._servo_library.TLI_GetDeviceInfo(self._serial_num, ctypes.byref(device_info))
                print("Description: ", device_info.description)
                print("Serial No: ", device_info.serialNo)
                print("USB PID: ", device_info.PID)

        if self.settings['velocity'] > 0.0: # update the velocity if not set to default
            self.set_velocity()
            # current_velocity = self.get_velocity()
            # print('velocity is set to ({0} mm/s)'.format(current_velocity))

        # open connection to the device. The other option is to NOT call _open_device, and require the script to handle device opening and closing
        self._open_device()

        # get the current position of the device and set the position setting to this value. That makes sure that the device does not move when you create the instrument.
        self.settings['position'] = self.get_position()
        # todo(emma): implement self.get_velocity()
        #self.velocity = self.get_velocity()

    def set_position(self, verbose=True):

        """
        sets position of device in mm
        """
        position = ctypes.c_int(int(self._position_encoder2deg_conversion_factor * self.settings['position']))
        self._servo_library.CC_MoveToPosition(self._serial_num, position) # command sent to device
        message_type = ctypes.c_ushort(0)
        message_id = ctypes.c_ushort(0)
        message_data = ctypes.c_ulong(0)
        if verbose:
            print('Device is moving to indicated position ({0} deg)'.format(self.settings['position']))

        # wait for it to actually move
        self._servo_library.CC_WaitForMessage(self._serial_num, message_type, message_id, message_data)
        while message_type.value != 2 or message_id.value != 1:
            self._servo_library.CC_WaitForMessage(self._serial_num, message_type, message_id, message_data)

        if verbose:
            position = self.get_position()
            print('Device now at position {0} deg'.format(position))

    def home(self, verbose=True):

        """
        Calibrates origin: scans full range and sets stage position to zero
        """

        self._servo_library.CC_Home(self._serial_num) # command sent to device!
        if verbose:
            print('Device is homing.')
        # wait for it to actually home
        message_type = ctypes.c_ushort(0)
        message_id = ctypes.c_ushort(0)
        message_data = ctypes.c_ulong(0)
        self._servo_library.CC_WaitForMessage(self._serial_num, message_type, message_id, message_data)
        while message_type.value != 2 or message_id.value != 0:
            self._servo_library.CC_WaitForMessage(self._serial_num, message_type, message_id, message_data)

        if verbose:
            print('Device finished homing')

    def _open_device(self, verbose=True):

        """
        Opens the serial connection to the device, which lets you send commands to and receive data from the device.
        Also starts polling all of the devices connected to the computer.
        """
        # print('_open_device now')

        if not self._servo_library.CC_Open(
                self._serial_num):
            if verbose:
                print('Starting to poll')
            milliseconds = ctypes.c_int(200)
            self._servo_library.CC_StartPolling(self._serial_num,
                                                milliseconds)  # continuous polling allows you to keep track of what the instrument is doing in real time
            time.sleep(3)
            self._servo_library.CC_ClearMessageQueue(self._serial_num)

    def _close_device(self, verbose=True):

        """
        Stops polling and closes the serial connection.
        """

        self._servo_library.CC_StopPolling(self._serial_num)
        self._servo_library.CC_Close(self._serial_num)
        if verbose:
            print('Device closed')

    def get_position(self, verbose=True):

        """
        returns position of stage in mm
        """
        print('getting current position now')
        position = self._servo_library.CC_GetPosition(
            self._serial_num) / self._position_encoder2deg_conversion_factor
        if verbose:
            print('Position of device is currently {0} deg'.format(position))
        return position

    def get_velocity(self, verbose=True):

        max_vel = 2.5

        """
        returns velocity (when moving) of stage in deg/s
        """
       # raise NotImplementedError

        ## todo: fix the input arguments of get_velocity ER 20180519
        print('getting current velocity now')
        velocity_pointer_type = ctypes.POINTER(ctypes.c_long) #ctypes.POINTER(ctypes.c_long)#ctypes.LP_c_long()#ctypes.POINTER(ctypes.c_int)
        acceleration_pointer_type = ctypes.POINTER(ctypes.c_long) #ctypes.POINTER(ctypes.c_long)#ctypes.LP_c_long()#ctypes.POINTER(ctypes.c_int)
        velocity_pointer = ctypes.byref(ctypes.c_long())#velocity_pointer_type()
        # velocity_pointer = ctypes.byref(ctypes.c_long(max_vel*self._velocity_encoder2mm_conversion_factor))  # velocity_pointer_type()

        acceleration_pointer = ctypes.byref(ctypes.c_long())#acceleration_pointer_type()
        # vel = self._servo_library.CC_GetVelParams(self._serial_num, velocity_pointer, acceleration_pointer)
        vel = self._servo_library.CC_GetVelParams(self._serial_num, acceleration_pointer, velocity_pointer)

        if verbose:
            print(vel)
            print(type(vel))
            print('Velocity of device is currently {0} deg/s'.format(vel))
        return vel

    def set_velocity(self, verbose=True):

        """
        sets velocity (when moving) of stage in mm/s
        """
        # maximum velocity of instrument
        max_vel = 2.5
        if self.settings['velocity'] > max_vel:
            self.settings['velocity'] = max_vel

        velocity = ctypes.c_int(int(self._velocity_encoder2deg_conversion_factor * self.settings['velocity']))
        acceleration = ctypes.c_int(int(self._acceleration))
        self._servo_library.CC_SetVelParams(self._serial_num, acceleration, velocity)
        print(velocity)
        print('Setting to indicated velocity: {0} deg/s'.format( velocity.value / self._velocity_encoder2deg_conversion_factor))
        if verbose:
            pass
            # current_vel = self.get_velocity()
            # print('Device velocity was set to {0} mm/s'.format(current_vel))
          # todo(emma): implement get_velocity()


    def update(self, settings, verbose=True):
        super(KDC101, self).update(settings)
        if self._settings_initialized:
            for key, value in settings.items():
                if key == 'position':
                    self.set_position()
                elif key == 'velocity':
                    self.set_velocity()

    def read_probes(self, key = None):
        assert key in list(self._PROBES.keys())
        assert isinstance(key, str)

        if key in ['position']:
            return (self.get_position(True))
        if key in ['serial_number']:
            return (self._serial_num_int)
        elif key in ['velocity']:
            # todo(emma): implement get_velocity
            pass

    def __del__(self):
        # when done with device we have to close the connections
        # called when GUI is closed
        self._close_device()

# class KDC101(Instrument):
#     '''
#     Class to control the thorlabs KDC101 servo. Note that ALL DLL FUNCTIONS TAKING NUMERIC INPUT REQUIRE A SYSTEM.DECIMAL
#     VALUE. Check help doc at C:\Program Files\Thorlabs\Kinesis\Thorlabs.MotionControl.DotNet_API for the DLL api.
#     The class communicates with the device over USB.
#     '''
#
#     startup_magnet_azimuth = 0
#
#     _DEFAULT_SETTINGS = Parameter([
#         Parameter('serial_number', 27001615, int, 'serial number written on device'),
#         Parameter('angle', 0, float, 'servo angle (from 0 to 360 deg)'),
#         Parameter('angular_velocity', 0, float, 'servo angular velocity in deg/s'),
#         Parameter('angular_acceleration', 0, float, 'servo angular acceleration in deg/s/s')
#     ])
#
#     def __init__(self, name = None, settings = None):
#         super(KDC101, self).__init__(name, settings)
#         try:
#             DeviceManagerCLI.BuildDeviceList()
#             serial_number_list = DeviceManagerCLI.GetDeviceList(KCubeDCServo.DevicePrefix)
#         except (Exception):
#             print("Exception raised by BuildDeviceList")
#         if not (str(self.settings['serial_number']) in serial_number_list):
#             print(str(self.settings['serial_number']) + " is not a valid serial number")
#             raise
#
#         self.device = KCubeDCServo.CreateKCubeDCServo(str(self.settings['serial_number']))
#         if(self.device == None):
#             print(self.settings['serial_number'] + " is not a KCubeDCServo")
#             raise
#
#         try:
#             self.device.Connect(str(self.settings['serial_number']))
#         except Exception:
#             print('Failed to open device ' + str(self.settings['serial_number']))
#             raise
#
#         if not self.device.IsSettingsInitialized():
#             try:
#                 self.device.WaitForSettingsInitialized(5000)
#             except Exception:
#                 print("Settings failed to initialize")
#                 raise
#
#         self.device.StartPolling(250)
#
#         motorSettings = self.device.GetMotorConfiguration(str(self.settings['serial_number']))
#         currentDeviceSettings = self.device.MotorDeviceSettings
#         self.startup_magnet_azimuth = self._get_position()
#         print('magnet azimuthal startup angle = ' + str(self.startup_magnet_azimuth) + "deg")
#
#     def update(self, settings):
#         '''
#         Updates internal settings, as well as the position and velocity set on the physical device
#         Args:
#             settings: A dictionary in the form of settings as seen in default settings
#         '''
#         super(KDC101, self).update(settings)
#         for key, value in settings.iteritems():
#             if key == 'position':
#                 self._move_servo(value)
#             elif key == 'velocity':
#                 self._set_velocity(value)
#
#     @property
#     def _PROBES(self):
#         return{
#             'position': 'servo position in mm',
#             'velocity': 'servo velocity in mm/s'
#         }
#
#     def read_probes(self, key):
#         assert key in self._PROBES.keys()
#         assert isinstance(key, str)
#
#         #query always returns string, need to cast to proper return type
#         if key in ['position']:
#             return self._get_position()
#         elif key in ['velocity']:
#             return self._get_velocity()
#
#     @property
#     def is_connected(self):
#         DeviceManagerCLI.BuildDeviceList()
#         return(str(self.settings['serial_number']) in DeviceManagerCLI.GetDeviceList(KCubeDCServo.DevicePrefix))
#
#     def __del__(self):
#         '''
#         Cleans up TDC001 connection
#         :PostState: TDC001 is disconnected
#         '''
#         self.device.StopPolling()
#         self.device.Disconnect()
#
#     def goto_home(self):
#         '''
#         Recenters device at the home position. Raises an exception on failure.
#         '''
#         try:
#             self.device.Home(60000)
#         except Exception:
#             print("Failed to move to position")
#             raise
#
#     def _move_servo(self, position, velocity = 0):
#         '''
#         Move servo to given position with given maximum velocity. Raises an exception on failure.
#         :param position: position in deg, ranges from 0-360
#         :param velocity: maximum velocity in deg/s, ranges from 0-2.5
#         :PostState: servo has moved
#         '''
#         if True:
#             try:
#                 if(velocity != 0):
#                     self._set_velocity(velocity)
#                 self.device.MoveTo(self._Py_Decimal(position), 60000)
#                 print("Moved magnet azimuth to " + str(position) + "deg")
#             except Exception:
#                 print("FAILED to move magnet azimuth to " + str(position) + "deg")
#                 # raise
#         else:
#             print(str(position) + " is out of allowed range [" + str(self.magnetZ_lower_lim)+","+str(self.magnetZ_upper_lim)+"]; no magnet azimuth move initiated")
#
#     def _get_position(self):
#         '''
#         :return: position of servo
#         '''
#         return self._Undo_Decimal(self.device.Position)
#
#     def _set_velocity(self, velocity):
#         '''
#         :param maximum velocity in mm/s, ranges from 0-2.5
#         :PostState: velocity changed in hardware
#         '''
#         if(velocity != 0):
#             velPars = self.device.GetVelocityParams()
#             velPars.MaxVelocity = self._Py_Decimal(velocity)
#             self.device.SetVelocityParams(velPars)
#
#     def _get_velocity(self):
#         '''
#         :return: maximum velocity setting
#         '''
#         return self._Undo_Decimal(self.device.GetVelocityParams().MaxVelocity)
#
#
#
#     def _Py_Decimal(self, value):
#         '''
#         Casting a python double to System.Decimal results in the Decimal having only integer values, likely due to an
#         improper selection of the overloaded Decimal function. Casting it first to System.Double, which always maintains
#         precision, then from Double to Decimal, where the proper overloaded function is clear, bypasses this issue
#         :param value: a python double
#         :return: the input as a System.Decimal
#         '''
#         return Decimal(Double(value))
#
#     def _Undo_Decimal(self, value):
#         '''
#         Casting back from System.Decimal to a python float fails due to overloading issues, but one can successfully
#         cast back to a string. Thus, we use a two-part cast to return to python numeric types
#         :param value: a System.Decimal
#         :return: the input as a python float
#         '''
#         return float(str(value))
class L607RT_KDC101_intensity_wheel(KDC101):
    '''

    Same as TDC001 except adds safety limits and specifies serial number for the x axis
    NOT WORKING!!!
    --> last editted by ZQ for use in LISE 607 NV RT confocal setup on 2/18/2019 6:51pm

    '''
    _DEFAULT_SETTINGS = Parameter([ # NB the serial number and min_, max_pos values will be overwritten
        Parameter('serial_number', 27500496, [27500496], 'serial number written on device'),
        Parameter('position', 20, float, 'servo position (0 to 360) [deg]'),
        Parameter('velocity', 1, float,
                  'servo velocity (0 to 2.5) [deg/s]. If set to zero, instrument default will be used')
    ])

    _PROBES = {'position': 'current position of stage', 'velocity':'current velocity of stage', 'serial_number': 'serial number of device'}

    def __init__(self, name=None, settings=None):
        self.max_pos = 360.
        self.min_pos = 0.
        super(L607RT_KDC101_intensity_wheel, self).__init__()

    def set_position(self, verbose=True):
        '''

        sets position of the stage x if safety limits are met

        '''

        # check if safety limits are met
        if self.settings['position'] < self.max_pos and self.settings['position'] > self.min_pos:
            super(L607RT_KDC101_intensity_wheel, self).set_position(self)
        else:
            print('didnt make the safety cut! doing nothing')
          #  raise AttributeError('position is outside safety limits!! Doing nothing.')

class L607RT_TDC001_sampleX(TDC001):
    '''

    Same as TDC001 except adds safety limits and specifies serial number for the x axis
    --> last editted by ZQ for use in LISE 607 NV RT confocal setup on 2/18/2019 6:51pm

    '''
    _DEFAULT_SETTINGS = Parameter([ # NB the serial number and min_, max_pos values will be overwritten
        Parameter('serial_number', 83860519, [83860519], 'serial number written on device'),
        Parameter('position', 5.9645, float, 'servo position (0 to 25) [mm]'),
        Parameter('velocity', 0.005, float,
                  'servo velocity (0 to 2.6) [mm/s]. If set to zero, instrument default will be used')
    ])

    _PROBES = {'position': 'current position of stage', 'velocity':'current velocity of stage', 'serial_number': 'serial number of device'}

    def __init__(self, name=None, settings=None):
        self.max_pos = 25.
        self.min_pos = 0.
        super(L607RT_TDC001_sampleX, self).__init__()

    def set_position(self, verbose=True):
        '''

        sets position of the stage x if safety limits are met

        '''

        # check if safety limits are met
        if self.settings['position'] < self.max_pos and self.settings['position'] > self.min_pos:
            super(L607RT_TDC001_sampleX, self).set_position(self)
        else:
            print('didnt make the safety cut! doing nothing')
          #  raise AttributeError('position is outside safety limits!! Doing nothing.')

class L607RT_TDC001_sampleY(TDC001):
    '''

    Same as TDC001 except adds safety limits and specifies serial number for the x axis
    --> last editted by ZQ for use in LISE 607 NV RT confocal setup on 2/18/2019 6:51pm

    '''
    _DEFAULT_SETTINGS = Parameter([ # NB the serial number and min_, max_pos values will be overwritten
        Parameter('serial_number', 83860538, [83860538], 'serial number written on device'),
        Parameter('position', 11.1209, float, 'servo position (0 to 25) [mm]'),
        Parameter('velocity', 0.005, float,
                  'servo velocity (0 to 2.6) [mm/s]. If set to zero, instrument default will be used')
    ])

    _PROBES = {'position': 'current position of stage', 'velocity':'current velocity of stage', 'serial_number': 'serial number of device'}

    def __init__(self, name=None, settings=None):
        self.max_pos = 25.
        self.min_pos = 0.
        super(L607RT_TDC001_sampleY, self).__init__()

    def set_position(self, verbose=True):
        '''

        sets position of the stage x if safety limits are met

        '''

        # check if safety limits are met
        if self.settings['position'] < self.max_pos and self.settings['position'] > self.min_pos:
            super(L607RT_TDC001_sampleY, self).set_position(self)
        else:
            print('didnt make the safety cut! doing nothing')
          #  raise AttributeError('position is outside safety limits!! Doing nothing.')

class L607RT_TDC001_sampleZ(TDC001):
    '''

    Same as TDC001 except adds safety limits and specifies serial number for the x axis
    --> last editted by ZQ for use in LISE 607 NV RT confocal setup on 2/18/2019 6:51pm

    '''
    _DEFAULT_SETTINGS = Parameter([ # NB the serial number and min_, max_pos values will be overwritten
        Parameter('serial_number', 83854376, [83854376], 'serial number written on device'),
        Parameter('position', 4.8155, float, 'servo position (0 to 25) [mm]'),
        Parameter('velocity', 0.005, float,
                  'servo velocity (0 to 2.6) [mm/s]. If set to zero, instrument default will be used')
    ])

    _PROBES = {'position': 'current position of stage', 'velocity':'current velocity of stage', 'serial_number': 'serial number of device'}

    def __init__(self, name=None, settings=None):
        self.max_pos = 25.
        self.min_pos = 0.
        super(L607RT_TDC001_sampleZ, self).__init__()

    def set_position(self, verbose=True):
        '''

        sets position of the stage x if safety limits are met

        '''

        # check if safety limits are met
        if self.settings['position'] < self.max_pos and self.settings['position'] > self.min_pos:
            super(L607RT_TDC001_sampleZ, self).set_position(self)
        else:
            print('didnt make the safety cut! doing nothing')
          #  raise AttributeError('position is outside safety limits!! Doing nothing.')

class L607RT_TDC001_magnetX(TDC001):
    '''

    Same as TDC001 except adds safety limits and specifies serial number for the x axis
    --> last editted by ZQ for use in LISE 607 NV RT confocal setup on 2/18/2019 6:51pm

    '''
    _DEFAULT_SETTINGS = Parameter([ # NB the serial number and min_, max_pos values will be overwritten
        Parameter('serial_number', 83860546, [83860546], 'serial number written on device'),
        Parameter('position', 12.9, float, 'servo position (0 to 25) [mm]'),
        Parameter('velocity', 0.005, float,
                  'servo velocity (0 to 2.6) [mm/s]. If set to zero, instrument default will be used')
    ])

    _PROBES = {'position': 'current position of stage', 'velocity':'current velocity of stage', 'serial_number': 'serial number of device'}

    def __init__(self, name=None, settings=None):
        self.max_pos = 25.
        self.min_pos = 0.
        super(L607RT_TDC001_magnetX, self).__init__()

    def set_position(self, verbose=True):
        '''

        sets position of the stage x if safety limits are met

        '''

        # check if safety limits are met
        if self.settings['position'] <= self.max_pos and self.settings['position'] >= self.min_pos:
            super(L607RT_TDC001_magnetX, self).set_position(self)
        else:
            print('desired position:', self.settings['position'])
            print('max position:', self.max_pos)
            print('min position:', self.min_pos)
            print('did not make the safety cut! doing nothing')
            # self.log('did not make the safety cut! doing nothing')
          #  raise AttributeError('position is outside safety limits!! Doing nothing.')

class L607RT_TDC001_magnetY(TDC001):
    '''

    Same as TDC001 except adds safety limits and specifies serial number for the x axis
    --> last editted by ZQ for use in LISE 607 NV RT confocal setup on 2/18/2019 6:51pm

    '''
    _DEFAULT_SETTINGS = Parameter([ # NB the serial number and min_, max_pos values will be overwritten
        Parameter('serial_number', 83860305, [83860305], 'serial number written on device'),
        Parameter('position', 5.7942, float, 'servo position (0 to 25) [mm]'),
        Parameter('velocity', 0.005, float,
                  'servo velocity (0 to 2.6) [mm/s]. If set to zero, instrument default will be used')
    ])

    _PROBES = {'position': 'current position of stage', 'velocity':'current velocity of stage', 'serial_number': 'serial number of device'}

    def __init__(self, name=None, settings=None):
        self.max_pos = 25.
        self.min_pos = 0.
        super(L607RT_TDC001_magnetY, self).__init__()

    def set_position(self, verbose=True):
        '''

        sets position of the stage x if safety limits are met

        '''

        # check if safety limits are met
        if self.settings['position'] < self.max_pos and self.settings['position'] > self.min_pos:
            super(L607RT_TDC001_magnetY, self).set_position(self)
        else:
            print('didnt make the safety cut! doing nothing')
          #  raise AttributeError('position is outside safety limits!! Doing nothing.')

class L607RT_TDC001_magnetZ(TDC001):
    '''

    Same as TDC001 except adds safety limits and specifies serial number for the x axis
    --> last editted by ZQ for use in LISE 607 NV RT confocal setup on 2/18/2019 6:51pm

    '''
    _DEFAULT_SETTINGS = Parameter([ # NB the serial number and min_, max_pos values will be overwritten
        Parameter('serial_number', 83860291, [83860291], 'serial number written on device'),
        Parameter('position', 19.6110, float, 'servo position (0 to 25) [mm]'),
        Parameter('velocity', 0.0, float,
                  'servo velocity (0 to 2.6) [mm/s]. If set to zero, instrument default will be used')
    ])

    _PROBES = {'position': 'current position of stage', 'velocity':'current velocity of stage', 'serial_number': 'serial number of device'}

    def __init__(self, name=None, settings=None):
        self.max_pos = 25.
        self.min_pos = 0.
        super(L607RT_TDC001_magnetZ, self).__init__()

    def set_position(self, verbose=True):
        '''

        sets position of the stage x if safety limits are met

        '''

        # check if safety limits are met
        if self.settings['position'] < self.max_pos and self.settings['position'] > self.min_pos:
            super(L607RT_TDC001_magnetZ, self).set_position(self)
        else:
            print('didnt make the safety cut! doing nothing')
          #  raise AttributeError('position is outside safety limits!! Doing nothing.')

# class KDC001(Instrument): Is this a typo? should be TDC001
#     """
#     Class to control the thorlabs KDC001 servo. Note that ALL DLL FUNCTIONS TAKING NUMERIC INPUT REQUIRE A SYSTEM.DECIMAL
#     VALUE. Check help doc at C:\Program Files\Thorlabs\Kinesis\Thorlabs.MotionControl.DotNet_API for the DLL api.
#     The class communicates with the device over USB.
#
#     --> last editted by Ziwei for use in LISE 607 NV RT confocal setup on 2/10/2019 3:37pm
#     --> modified from B26 class
#     """
#
#     _DEFAULT_SETTINGS = Parameter([
#         Parameter('serial_number', 27501971, int, 'serial number written on device'),
#         Parameter('position', 3.0, float, 'servo position (0 to 25) [mm]'),
#         Parameter('velocity', 0.01, float, 'servo velocity (0 to 2.6) [mm/s]. If set to zero, instrument default will be used')
#     ])
#
#     _PROBES = {}
#
#     def __init__(self, name=None, settings=None):
#
#         super(KDC001, self).__init__(name, settings)
#         self._servo_library = ctypes.cdll.LoadLibrary('C:\\Program Files\\Thorlabs\\Kinesis\\Thorlabs.MotionControl.KCube.DCServo.dll')
#
#         # conversion from mm to device units at
#         # https://www.thorlabs.com/software/apt/APT_Communications_Protocol_Rev_15.pdf, page 16
#         self._position_encoder2mm_conversion_factor = 34304
#         self._velocity_encoder2mm_conversion_factor = 767367.49
#         self._acceleration_encoder2mm_conversion_factor = 261.93
#         self._acceleration = 1 * self._acceleration_encoder2mm_conversion_factor # use hard-coded acceleration of 1 mm/s^2 (can be changed to parameter if we want)
#         self._manually_set_library_inputs_and_outputs()
#         self._connect()
#
#     def _manually_set_library_inputs_and_outputs(self):
#         """
#         Sets the input and output types for each servo library call we make.
#
#         """
#         self._servo_library.TLI_BuildDeviceList.restypes = ctypes.c_short
#
#         self._servo_library.TLI_GetDeviceListSize.restypes = ctypes.c_short
#
#         self._servo_library.TLI_GetDeviceListByTypeExt.argtypes = [ctypes.POINTER(ctypes.c_char), ctypes.c_ulong,
#                                                                    ctypes.c_int]
#         self._servo_library.TLI_GetDeviceListByTypeExt.restypes = ctypes.c_short
#
#         self._servo_library.TLI_GetDeviceInfo.argtypes = [ctypes.POINTER(ctypes.c_char), ctypes.POINTER(TLI_DeviceInfo)]
#         self._servo_library.TLI_GetDeviceInfo.restypes = ctypes.c_short
#
#         self._servo_library.CC_StartPolling.argtypes = [ctypes.POINTER(ctypes.c_char), ctypes.c_int]
#         self._servo_library.CC_StartPolling.restypes = ctypes.c_bool
#
#         self._servo_library.CC_Open.argtypes = [ctypes.POINTER(ctypes.c_char)]
#         self._servo_library.CC_Open.restypes = ctypes.c_short
#
#         self._servo_library.CC_Close.argtypes = [ctypes.POINTER(ctypes.c_char)]
#         self._servo_library.CC_Close.restypes = ctypes.c_short
#
#         self._servo_library.CC_ClearMessageQueue.argtypes = [ctypes.POINTER(ctypes.c_char)]
#         self._servo_library.CC_ClearMessageQueue.restypes = ctypes.c_short
#
#         self._servo_library.CC_WaitForMessage.argtypes = [ctypes.POINTER(ctypes.c_char),
#                                                           ctypes.POINTER(ctypes.c_ushort),
#                                                           ctypes.POINTER(ctypes.c_ushort),
#                                                           ctypes.POINTER(ctypes.c_ulong)]
#         self._servo_library.CC_WaitForMessage.restypes = ctypes.c_bool
#
#         self._servo_library.CC_GetVelParams.argtypes = [ctypes.POINTER(ctypes.c_char),
#                                                         ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long)]
#
#         self._servo_library.CC_GetVelParams.restypes = ctypes.c_short
#
#         self._servo_library.CC_MoveToPosition.argtypes = [ctypes.POINTER(ctypes.c_char), ctypes.c_int]
#         self._servo_library.CC_MoveToPosition.restypes = ctypes.c_short
#
#         self._servo_library.CC_SetVelParams.argtypes = [ctypes.POINTER(ctypes.c_char), ctypes.c_int, ctypes.c_int]
#         self._servo_library.CC_SetVelParams.restypes = ctypes.c_short
#
#         self._servo_library.CC_StopPolling.argtypes = [ctypes.POINTER(ctypes.c_char)]
#
#         self._servo_library.CC_Home.argtypes = [ctypes.POINTER(ctypes.c_char)]
#         self._servo_library.CC_Home.restypes = ctypes.c_short
#
#         self._servo_library.CC_GetPosition.argtypes = [ctypes.POINTER(ctypes.c_char)]
#         self._servo_library.CC_GetPosition.restypes = ctypes.c_int
#
#     def _connect(self, verbose=True):
#
#         """
#         makes connection to the device through the serial port. Also opens connection to the device.
#         """
#
#         # this version of the serial number is useful
#         self._serial_num_int = self.settings['serial_number']
#         self._serial_num = ctypes.c_char_p(bytes(str(self.settings['serial_number']), "utf-8"))
#
#         if self._servo_library.TLI_BuildDeviceList() == 0:
#             num_devices = self._servo_library.TLI_GetDeviceListSize()
#             if verbose:
#                 print('Number of devices detected: {0}'.format(num_devices))
#
#             # The servo library fills in a byte string of connected device serial numbers, separated by a comma.
#             # The length of this string will be the length of the serial number (8 bytes), plus a byte for a comma,
#             # for each connected device.
#             # Here, we first pre-allocate the string, then send it to the library function, and then examine it.
#
#             connected_devices_string_length = num_devices*9
#             string_of_serial_nums = ctypes.create_string_buffer(connected_devices_string_length)
#             self._servo_library.TLI_GetDeviceListByTypeExt(string_of_serial_nums,
#                                                            ctypes.c_ulong(connected_devices_string_length+1),
#                                                            ctypes.c_int(27))
#
#             list_of_connected_serial_nums = string_of_serial_nums.raw.decode("utf-8")
#
#             if str(self.settings['serial_number']) not in list_of_connected_serial_nums:
#                 error_msg = 'No servo with the given serial number was detected.\n'
#                 error_msg += 'Given Serial Number: {0}\n'.format(self.settings['serial_number'])
#                 if list_of_connected_serial_nums:
#                     error_msg += 'Connected Devices: {0}\n'.format(list_of_connected_serial_nums)
#                 else:
#                     error_msg += 'Connected Devices: None\n'
#                 raise AttributeError(error_msg)
#
#             elif verbose:
#                 print('Found device with matching serial number.')
#                 device_info = TLI_DeviceInfo()
#                 self._servo_library.TLI_GetDeviceInfo(self._serial_num, ctypes.byref(device_info))
#                 print("Description: ", device_info.description)
#                 print("Serial No: ", device_info.serialNo)
#                 print("USB PID: ", device_info.PID)
#
#         if self.settings['velocity'] > 0.0: # update the velocity if not set to default
#             self.set_velocity()
#
#         # open connection to the device. The other option is to NOT call _open_device, and require the script to handle device opening and closing
#         self._open_device()
#
#         # get the current position of the device and set the position setting to this value. That makes sure that the device does not move when you create the instrument.
#         self.settings['position'] = self.get_position()
#         # todo(emma): implement self.get_velocity()
#         #self.velocity = self.get_velocity()
#
#     def set_position(self, verbose=True):
#
#         """
#         sets position of device in mm
#         """
#         position = ctypes.c_int(int(self._position_encoder2mm_conversion_factor * self.settings['position']))
#         self._servo_library.CC_MoveToPosition(self._serial_num, position) # command sent to device
#         message_type = ctypes.c_ushort(0)
#         message_id = ctypes.c_ushort(0)
#         message_data = ctypes.c_ulong(0)
#         if verbose:
#             print('Device is moving to indicated position ({0} mm)'.format(self.settings['position']))
#
#         # wait for it to actually move
#         self._servo_library.CC_WaitForMessage(self._serial_num, message_type, message_id, message_data)
#         while message_type.value != 2 or message_id.value != 1:
#             self._servo_library.CC_WaitForMessage(self._serial_num, message_type, message_id, message_data)
#
#         if verbose:
#             position = self.get_position()
#             print('Device now at position {0} mm'.format(position))
#
#     def home(self, verbose=True):
#
#         """
#         Calibrates origin: scans full range and sets stage position to zero
#         """
#
#         self._servo_library.CC_Home(self._serial_num) # command sent to device!
#         if verbose:
#             print('Device is homing.')
#         # wait for it to actually home
#         message_type = ctypes.c_ushort(0)
#         message_id = ctypes.c_ushort(0)
#         message_data = ctypes.c_ulong(0)
#         self._servo_library.CC_WaitForMessage(self._serial_num, message_type, message_id, message_data)
#         while message_type.value != 2 or message_id.value != 0:
#             self._servo_library.CC_WaitForMessage(self._serial_num, message_type, message_id, message_data)
#
#         if verbose:
#             print('Device finished homing')
#
#     def _open_device(self, verbose=True):
#
#         """
#         Opens the serial connection to the device, which lets you send commands to and receive data from the device.
#         Also starts polling all of the devices connected to the computer.
#         """
#
#         if not self._servo_library.CC_Open(
#                 self._serial_num):
#             if verbose:
#                 print('Starting to poll')
#             milliseconds = ctypes.c_int(200)
#             self._servo_library.CC_StartPolling(self._serial_num,
#                                                 milliseconds)  # continuous polling allows you to keep track of what the instrument is doing in real time
#             time.sleep(3)
#             self._servo_library.CC_ClearMessageQueue(self._serial_num)
#
#     def _close_device(self, verbose=True):
#
#         """
#         Stops polling and closes the serial connection.
#         """
#
#         self._servo_library.CC_StopPolling(self._serial_num)
#         self._servo_library.CC_Close(self._serial_num)
#         if verbose:
#             print('Device closed')
#
#     def get_position(self, verbose=True):
#
#         """
#         returns position of stage in mm
#         """
#
#         position = self._servo_library.CC_GetPosition(
#             self._serial_num) / self._position_encoder2mm_conversion_factor
#         if verbose:
#             print('Position of device is currently {0} mm'.format(position))
#         return position
#
#     def get_velocity(self, verbose=True):
#
#         """
#         returns velocity (when moving) of stage in mm/s
#         """
#        # raise NotImplementedError
#
#         ## todo: fix the input arguments of get_velocity ER 20180519
#         velocity_pointer_type = ctypes.POINTER(ctypes.c_long) #ctypes.POINTER(ctypes.c_long)#ctypes.LP_c_long()#ctypes.POINTER(ctypes.c_int)
#         acceleration_pointer_type = ctypes.POINTER(ctypes.c_long) #ctypes.POINTER(ctypes.c_long)#ctypes.LP_c_long()#ctypes.POINTER(ctypes.c_int)
#         velocity_pointer = ctypes.byref(ctypes.c_long())#velocity_pointer_type()
#         acceleration_pointer = ctypes.byref(ctypes.c_long())#acceleration_pointer_type()
#         vel = self._servo_library.CC_GetVelParams(self.serial_num, velocity_pointer, acceleration_pointer)
#         if verbose:
#             print('Velocity of device is currently {0} mm/s'.format(vel))
#         return vel
#
#     def set_velocity(self, verbose=True):
#
#         """
#         sets velocity (when moving) of stage in mm/s
#         """
#         # maximum velocity of instrument
#         max_vel = 2.6
#         if self.settings['velocity'] > max_vel:
#             self.settings['velocity'] = max_vel
#
#         velocity = ctypes.c_int(int(self._velocity_encoder2mm_conversion_factor * self.settings['velocity']))
#         acceleration = ctypes.c_int(int(self._acceleration))
#         self._servo_library.CC_SetVelParams(self._serial_num, acceleration, velocity)
#         if verbose:
#             pass
#           #  set_vel = self.get_velocity()
#           # todo(emma): implement get_velocity()
#            # print('Device velocity was set to {0} mm/s'.format(set_vel))
#
#     def update(self, settings, verbose=True):
#         super(KDC001, self).update(settings)
#         if self._settings_initialized:
#             for key, value in settings.items():
#                 if key == 'position':
#                     self.set_position()
#                 elif key == 'velocity':
#                     self.set_velocity()
#
#     def read_probes(self, key = None):
#         assert key in list(self._PROBES.keys())
#         assert isinstance(key, str)
#
#         if key in ['position']:
#             return (self.get_position(True))
#         if key in ['serial_number']:
#             return (self._serial_num_int)
#         elif key in ['velocity']:
#             # todo(emma): implement get_velocity
#             pass
#
#     def __del__(self):
#         # when done with device we have to close the connections
#         # called when GUI is closed
#         self._close_device()
# class KDC001_old(Instrument):
#     """
#     Class to control the thorlabs KDC001 servo. Note that ALL DLL FUNCTIONS TAKING NUMERIC INPUT REQUIRE A SYSTEM.DECIMAL
#     VALUE. Check help doc at C:\Program Files\Thorlabs\Kinesis\Thorlabs.MotionControl.DotNet_API for the DLL api.
#     The class communicates with the device over USB.
#     """
#
#     _DEFAULT_SETTINGS = Parameter([
#         Parameter('serial_number', 27501971, int, 'serial number written on device'),
#         Parameter('position', 3.0, float, 'servo position (0 to 25) [mm]'),
#         Parameter('velocity', 1.0, float, 'servo velocity (0 to 2.6) [mm/s]. If set to zero, instrument default will be used')
#     ])
#
#     _PROBES = {}
#
#     def __init__(self, name=None, settings=None):
#
#         super(KDC001, self).__init__(name, settings)
#         self._servo_library = ctypes.cdll.LoadLibrary('C:\\Program Files\\Thorlabs\\Kinesis\\Thorlabs.MotionControl.KCube.DCServo.dll')
#
#         # conversion from mm to device units at
#         # https://www.thorlabs.com/software/apt/APT_Communications_Protocol_Rev_15.pdf, page 16
#         self._position_encoder2mm_conversion_factor = 34304
#         self._velocity_encoder2mm_conversion_factor = 767367.49
#         self._acceleration_encoder2mm_conversion_factor = 261.93
#         self._acceleration = 1 * self._acceleration_encoder2mm_conversion_factor # use hard-coded acceleration of 1 mm/s^2 (can be changed to parameter if we want)
#         self._manually_set_library_inputs_and_outputs()
#         self._connect()
#
#     def _manually_set_library_inputs_and_outputs(self):
#         """
#         Sets the input and output types for each servo library call we make.
#
#         """
#         self._servo_library.TLI_BuildDeviceList.restypes = ctypes.c_short
#
#         self._servo_library.TLI_GetDeviceListSize.restypes = ctypes.c_short
#
#         self._servo_library.TLI_GetDeviceListByTypeExt.argtypes = [ctypes.POINTER(ctypes.c_char), ctypes.c_ulong,
#                                                                    ctypes.c_int]
#         self._servo_library.TLI_GetDeviceListByTypeExt.restypes = ctypes.c_short
#
#         self._servo_library.TLI_GetDeviceInfo.argtypes = [ctypes.POINTER(ctypes.c_char), ctypes.POINTER(TLI_DeviceInfo)]
#         self._servo_library.TLI_GetDeviceInfo.restypes = ctypes.c_short
#
#         self._servo_library.CC_StartPolling.argtypes = [ctypes.POINTER(ctypes.c_char), ctypes.c_int]
#         self._servo_library.CC_StartPolling.restypes = ctypes.c_bool
#
#         self._servo_library.CC_Open.argtypes = [ctypes.POINTER(ctypes.c_char)]
#         self._servo_library.CC_Open.restypes = ctypes.c_short
#
#         self._servo_library.CC_Close.argtypes = [ctypes.POINTER(ctypes.c_char)]
#         self._servo_library.CC_Close.restypes = ctypes.c_short
#
#         self._servo_library.CC_ClearMessageQueue.argtypes = [ctypes.POINTER(ctypes.c_char)]
#         self._servo_library.CC_ClearMessageQueue.restypes = ctypes.c_short
#
#         self._servo_library.CC_WaitForMessage.argtypes = [ctypes.POINTER(ctypes.c_char),
#                                                           ctypes.POINTER(ctypes.c_ushort),
#                                                           ctypes.POINTER(ctypes.c_ushort),
#                                                           ctypes.POINTER(ctypes.c_ulong)]
#         self._servo_library.CC_WaitForMessage.restypes = ctypes.c_bool
#
#         self._servo_library.CC_GetVelParams.argtypes = [ctypes.POINTER(ctypes.c_char),
#                                                         ctypes.POINTER(ctypes.c_long), ctypes.POINTER(ctypes.c_long)]
#
#         self._servo_library.CC_GetVelParams.restypes = ctypes.c_short
#
#         self._servo_library.CC_MoveToPosition.argtypes = [ctypes.POINTER(ctypes.c_char), ctypes.c_int]
#         self._servo_library.CC_MoveToPosition.restypes = ctypes.c_short
#
#         self._servo_library.CC_SetVelParams.argtypes = [ctypes.POINTER(ctypes.c_char), ctypes.c_int, ctypes.c_int]
#         self._servo_library.CC_SetVelParams.restypes = ctypes.c_short
#
#         self._servo_library.CC_StopPolling.argtypes = [ctypes.POINTER(ctypes.c_char)]
#
#         self._servo_library.CC_Home.argtypes = [ctypes.POINTER(ctypes.c_char)]
#         self._servo_library.CC_Home.restypes = ctypes.c_short
#
#         self._servo_library.CC_GetPosition.argtypes = [ctypes.POINTER(ctypes.c_char)]
#         self._servo_library.CC_GetPosition.restypes = ctypes.c_int
#
#     def _connect(self, verbose=True):
#
#         """
#         makes connection to the device through the serial port. Also opens connection to the device.
#         """
#
#         # this version of the serial number is useful
#         self._serial_num_int = self.settings['serial_number']
#         self._serial_num = ctypes.c_char_p(bytes(str(self.settings['serial_number']), "utf-8"))
#
#         if self._servo_library.TLI_BuildDeviceList() == 0:
#             num_devices = self._servo_library.TLI_GetDeviceListSize()
#             if verbose:
#                 print('Number of devices detected: {0}'.format(num_devices))
#
#             # The servo library fills in a byte string of connected device serial numbers, separated by a comma.
#             # The length of this string will be the length of the serial number (8 bytes), plus a byte for a comma,
#             # for each connected device.
#             # Here, we first pre-allocate the string, then send it to the library function, and then examine it.
#
#             connected_devices_string_length = num_devices*9
#             string_of_serial_nums = ctypes.create_string_buffer(connected_devices_string_length)
#             self._servo_library.TLI_GetDeviceListByTypeExt(string_of_serial_nums,
#                                                            ctypes.c_ulong(connected_devices_string_length+1),
#                                                            ctypes.c_int(27))
#
#             list_of_connected_serial_nums = string_of_serial_nums.raw.decode("utf-8")
#
#             if str(self.settings['serial_number']) not in list_of_connected_serial_nums:
#                 error_msg = 'No servo with the given serial number was detected.\n'
#                 error_msg += 'Given Serial Number: {0}\n'.format(self.settings['serial_number'])
#                 if list_of_connected_serial_nums:
#                     error_msg += 'Connected Devices: {0}\n'.format(list_of_connected_serial_nums)
#                 else:
#                     error_msg += 'Connected Devices: None\n'
#                 raise AttributeError(error_msg)
#
#             elif verbose:
#                 print('Found device with matching serial number.')
#                 device_info = TLI_DeviceInfo()
#                 self._servo_library.TLI_GetDeviceInfo(self._serial_num, ctypes.byref(device_info))
#                 print("Description: ", device_info.description)
#                 print("Serial No: ", device_info.serialNo)
#                 print("USB PID: ", device_info.PID)
#
#         if self.settings['velocity'] > 0.0: # update the velocity if not set to default
#             self.set_velocity()
#
#         # open connection to the device. The other option is to NOT call _open_device, and require the script to handle device opening and closing
#         self._open_device()
#
#         # get the current position of the device and set the position setting to this value. That makes sure that the device does not move when you create the instrument.
#         self.settings['position'] = self.get_position()
#         # todo(emma): implement self.get_velocity()
#         #self.velocity = self.get_velocity()
#
#     def set_position(self, verbose=True):
#
#         """
#         sets position of device in mm
#         """
#         position = ctypes.c_int(int(self._position_encoder2mm_conversion_factor * self.settings['position']))
#         self._servo_library.CC_MoveToPosition(self._serial_num, position) # command sent to device
#         message_type = ctypes.c_ushort(0)
#         message_id = ctypes.c_ushort(0)
#         message_data = ctypes.c_ulong(0)
#         if verbose:
#             print('Device is moving to indicated position ({0} mm)'.format(self.settings['position']))
#
#         # wait for it to actually move
#         self._servo_library.CC_WaitForMessage(self._serial_num, message_type, message_id, message_data)
#         while message_type.value != 2 or message_id.value != 1:
#             self._servo_library.CC_WaitForMessage(self._serial_num, message_type, message_id, message_data)
#
#         if verbose:
#             position = self.get_position()
#             print('Device now at position {0} mm'.format(position))
#
#     def home(self, verbose=True):
#
#         """
#         Calibrates origin: scans full range and sets stage position to zero
#         """
#
#         self._servo_library.CC_Home(self._serial_num) # command sent to device!
#         if verbose:
#             print('Device is homing.')
#         # wait for it to actually home
#         message_type = ctypes.c_ushort(0)
#         message_id = ctypes.c_ushort(0)
#         message_data = ctypes.c_ulong(0)
#         self._servo_library.CC_WaitForMessage(self._serial_num, message_type, message_id, message_data)
#         while message_type.value != 2 or message_id.value != 0:
#             self._servo_library.CC_WaitForMessage(self._serial_num, message_type, message_id, message_data)
#
#         if verbose:
#             print('Device finished homing')
#
#     def _open_device(self, verbose=True):
#
#         """
#         Opens the serial connection to the device, which lets you send commands to and receive data from the device.
#         Also starts polling all of the devices connected to the computer.
#         """
#
#         if not self._servo_library.CC_Open(
#                 self._serial_num):
#             if verbose:
#                 print('Starting to poll')
#             milliseconds = ctypes.c_int(200)
#             self._servo_library.CC_StartPolling(self._serial_num,
#                                                 milliseconds)  # continuous polling allows you to keep track of what the instrument is doing in real time
#             time.sleep(3)
#             self._servo_library.CC_ClearMessageQueue(self._serial_num)
#
#     def _close_device(self, verbose=True):
#
#         """
#         Stops polling and closes the serial connection.
#         """
#
#         self._servo_library.CC_StopPolling(self._serial_num)
#         self._servo_library.CC_Close(self._serial_num)
#         if verbose:
#             print('Device closed')
#
#     def get_position(self, verbose=True):
#
#         """
#         returns position of stage in mm
#         """
#
#         position = self._servo_library.CC_GetPosition(
#             self._serial_num) / self._position_encoder2mm_conversion_factor
#         if verbose:
#             print('Position of device is currently {0} mm'.format(position))
#         return position
#
#     def get_velocity(self, verbose=True):
#
#         """
#         returns velocity (when moving) of stage in mm/s
#         """
#        # raise NotImplementedError
#
#         ## todo: fix the input arguments of get_velocity ER 20180519
#         velocity_pointer_type = ctypes.POINTER(ctypes.c_long) #ctypes.POINTER(ctypes.c_long)#ctypes.LP_c_long()#ctypes.POINTER(ctypes.c_int)
#         acceleration_pointer_type = ctypes.POINTER(ctypes.c_long) #ctypes.POINTER(ctypes.c_long)#ctypes.LP_c_long()#ctypes.POINTER(ctypes.c_int)
#         velocity_pointer = ctypes.byref(ctypes.c_long())#velocity_pointer_type()
#         acceleration_pointer = ctypes.byref(ctypes.c_long())#acceleration_pointer_type()
#         vel = self._servo_library.CC_GetVelParams(self.serial_num, velocity_pointer, acceleration_pointer)
#         if verbose:
#             print('Velocity of device is currently {0} mm/s'.format(vel))
#         return vel
#
#     def set_velocity(self, verbose=True):
#
#         """
#         sets velocity (when moving) of stage in mm/s
#         """
#         # maximum velocity of instrument
#         max_vel = 2.6
#         if self.settings['velocity'] > max_vel:
#             self.settings['velocity'] = max_vel
#
#         velocity = ctypes.c_int(int(self._velocity_encoder2mm_conversion_factor * self.settings['velocity']))
#         acceleration = ctypes.c_int(int(self._acceleration))
#         self._servo_library.CC_SetVelParams(self._serial_num, acceleration, velocity)
#         if verbose:
#             pass
#           #  set_vel = self.get_velocity()
#           # todo(emma): implement get_velocity()
#            # print('Device velocity was set to {0} mm/s'.format(set_vel))
#
#     def update(self, settings, verbose=True):
#         super(KDC001, self).update(settings)
#         if self._settings_initialized:
#             for key, value in settings.items():
#                 if key == 'position':
#                     self.set_position()
#                 elif key == 'velocity':
#                     self.set_velocity()
#
#     def read_probes(self, key = None):
#         assert key in list(self._PROBES.keys())
#         assert isinstance(key, str)
#
#         if key in ['position']:
#             return (self.get_position(True))
#         if key in ['serial_number']:
#             return (self._serial_num_int)
#         elif key in ['velocity']:
#             # todo(emma): implement get_velocity
#             pass
#
#     def __del__(self):
#         # when done with device we have to close the connections
#         # called when GUI is closed
#         self._close_device()
#
# class B26KDC001x_old(KDC001):
#     '''
#
#     Same as KDC001 except adds safety limits and specifies serial number for the x axis
#
#     '''
#     _DEFAULT_SETTINGS = Parameter([ # NB the serial number and min_, max_pos values will be overwritten
#         Parameter('serial_number', 27501971, int, 'serial number written on device'),
#         Parameter('position', 3.0, float, 'servo position (0 to 25) [mm]'),
#         Parameter('velocity', 0.0, float,
#                   'servo velocity (0 to 2.6) [mm/s]. If set to zero, instrument default will be used')
#     ])
#
#     _PROBES = {'position': 'current position of stage', 'velocity':'current velocity of stage', 'serial_number': 'serial number of device'}
#
#     def __init__(self, name=None, settings=None):
#         self.max_pos = 25.
#         self.min_pos = 0.
#         super(B26KDC001x, self).__init__()
#
#     def set_position(self, verbose=True):
#         '''
#
#         sets position of the stage x if safety limits are met
#
#         '''
#
#         # check if safety limits are met
#         if self.settings['position'] < self.max_pos and self.settings['position'] > self.min_pos:
#             super(B26KDC001x, self).set_position(self)
#         else:
#             print('didnt make the safety cut! doing nothing')
#           #  raise AttributeError('position is outside safety limits!! Doing nothing.')
#
# class B26KDC001y_old(KDC001):
#     '''
#
#     Same as KDC001 except adds safety limits and specifies serial number for the y axis
#
#     '''
#     _DEFAULT_SETTINGS = Parameter([ # NB the serial number and min_, max_pos values will be overwritten
#         Parameter('serial_number', 27501986, int, 'serial number written on device'),
#         Parameter('position', 3.0, float, 'servo position (0 to 25) [mm]'),
#         Parameter('velocity', 0.0, float,
#                   'servo velocity (0 to 2.6) [mm/s]. If set to zero, instrument default will be used')
#     ])
#
#     _PROBES = {'position': 'current position of stage', 'velocity':'current velocity of stage', 'serial_number': 'serial number of device'}
#
#     def __init__(self, name=None, settings=None):
#         self.max_pos = 19.7
#         self.min_pos = 0.
#         super(B26KDC001y, self).__init__()
#
#     def set_position(self, verbose=True):
#         '''
#
#         sets position of the stage y if safety limits are met
#
#         '''
#
#         # check if safety limits are met
#         if self.settings['position'] < self.max_pos and self.settings['position'] > self.min_pos:
#             super(B26KDC001y, self).set_position(self)
#         else:
#             print('didnt make the safety cut! doing nothing')
#           #  raise AttributeError('position is outside safety limits!! Doing nothing.')
#
# class B26KDC001z_old(KDC001):
#     '''
#
#     Same as KDC001 except adds safety limits and specifies serial number for the y axis
#
#     '''
#     _DEFAULT_SETTINGS = Parameter([ # NB the serial number and min_, max_pos values will be overwritten
#         Parameter('serial_number', 27001862, int, 'serial number written on device'),
#         Parameter('position', 3.0, float, 'servo position (0 to 25) [mm]'),
#         Parameter('velocity', 0.0, float,
#                   'servo velocity (0 to 2.6) [mm/s]. If set to zero, instrument default will be used')
#     ])
#
#     _PROBES = {'position': 'current position of stage', 'velocity':'current velocity of stage', 'serial_number': 'serial number of device'}
#
#     def __init__(self, name=None, settings=None):
#         self.max_pos = 16.7 #16.0
#         self.min_pos = 0.
#         super(B26KDC001z, self).__init__()
#
#     def set_position(self, verbose=True):
#         '''
#
#         sets position of the stage z if safety limits are met
#
#         '''
#
#         # check if safety limits are met
#         if self.settings['position'] < self.max_pos and self.settings['position'] > self.min_pos:
#             super(B26KDC001z, self).set_position(self)
#         else:
#             print('didnt make the safety cut! doing nothing')
#           #  raise AttributeError('position is outside safety limits!! Doing nothing.')

if __name__ == '__main__':
    a = L607RT_TDC001_sampleX()
    b = L607RT_TDC001_sampleY()
    c = L607RT_TDC001_sampleZ()
    d = L607RT_TDC001_magnetX()
    e = L607RT_TDC001_magnetY()
    f = L607RT_TDC001_magnetZ()
  #  a.set_velocity()
 #   a.set_position()