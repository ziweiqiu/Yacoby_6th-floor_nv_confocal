import serial
from pylabcontrol.core import Instrument, Parameter

class WhiteLEDController(Instrument):
    """
    This class controls the white light intensity with a current source. The current source is recognized as a serial port.
    The class communicates with the device over RS232 using pyserial.

    --> written by ZQ 2/21/2019 2:43pm

    """

    # ASCII Characters used for controller communication
    # ETX = chr(3)
    # CR = chr(13)
    # LF = chr(10)
    # ENQ = chr(5)
    # ACK = chr(6)
    # NAK = chr(21)

    # _possible_com_ports = ['COM' + str(i) for i in range(0, 256)]

    _DEFAULT_SETTINGS = Parameter([
        # Parameter('port', 'COM5', _possible_com_ports, 'com port to which the white light controller is connected'),
        Parameter('port', 'COM5', ['COM5'], 'com port to which the white light controller is connected'),
        Parameter('intensity', 20.0, float, 'light intensity between 0 and 100'),
        Parameter('timeout', 1.0, float, 'amount of time to wait for a response from the white led controller for each query'),
        Parameter('baudrate', 9600, int, 'baudrate of serial communication with gauge')
    ])
    # _DEFAULT_SETTINGS = Parameter([
    #     Parameter('serial_number', 83860546, [83860546], 'serial number written on device'),
    #     Parameter('position', 3.0, float, 'servo position (0 to 25) [mm]'),
    #     Parameter('velocity', 0.01, float,
    #               'servo velocity (0 to 2.6) [mm/s]. If set to zero, instrument default will be used')
    # ])

    _PROBES = {} # this is necessary otherwise the software will crash upon expanding the settings....

    def __init__(self, name=None, settings=None):
        """
        The serial connection should be setup with the following parameters:
        1 start bit, 8 data bits, 1 stop bit, No parity bit, no hardware
        handshake. These are all default for Serial and therefore not input
        below
        """
        super(WhiteLEDController, self).__init__(name, settings)

        self._connect()

    def _connect(self, verbose = False):

        try:
            self.serial_connection = serial.Serial(port=self.settings['port'], baudrate=self.settings['baudrate'],
                                                   timeout=self.settings['timeout'], parity=serial.PARITY_NONE,
                                                   stopbits=serial.STOPBITS_ONE, bytesize=serial.EIGHTBITS)
            if verbose:
                if self.is_connected():
                    print('port is connected.')
                else:
                    print('fails to connect the white LED light controller')

        except IOError:  # if port is already opened, close it and open it again and print message
            self.serial_connection.close()
            self.serial_connection.open()
            if verbose:
                print("port was already open, was closed and opened again!")

    def update(self, settings):
        super(WhiteLEDController, self).update(settings)
        for key, value in settings.items():
            if key == 'port':
                self._connect()
            elif key == 'intensity':
                if self.is_connected():
                    # self.serial_connection.write('KRDG? A' + self.CR + self.LF)
                    MAX_number = 65535
                    MIN_number = 0
                    intensity = round(MAX_number*self.settings['intensity']/100.0)
                    intensity_str ='B' + str(intensity) +  ';'
                    # intensity_str = 'B' + str(intensity)
                    print('Setting intensity = {0} %'.format(self.settings['intensity']))

                    # print(self.settings['intensity'])
                    # print('Writing to the serial port:')
                    # print(intensity_str)
                    self.serial_connection.write(intensity_str.encode('utf-8'))
                else:
                    print('Fails to write because the white LED light controller is not connected!')

    def is_connected(self):
        """
        checks if serial connection is still open with instrument.

        :return: boolean connection status
        """
        return self.serial_connection.isOpen()

    def __del__(self):
        """
        Destructor, to close the serial connection when the instance is this class is garbage collected
        """
        self.serial_connection.close()
        if not self.is_connected():
            print('port is closed.')
        else:
            print('port is still open.')


if __name__ == '__main__':
    # instruments, failed = Instrument.load_and_append(
    #     instrument_dict={'WhiteLEDController': WhiteLEDController})
    a = WhiteLEDController()


        # print((instruments['WhiteLEDController']._get_temperature()))



