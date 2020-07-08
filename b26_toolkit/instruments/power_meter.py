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

import visa
import pyvisa.errors
from pylabcontrol.core import Parameter, Instrument
import time


class PM100D(Instrument):
    """
        This class implements the Thorlabs power meter PM100D. The class commuicates with the device over USB.
        manual: https://www.physics.utoronto.ca/~phy326/opt/PM100D-Manual.pdf

        -- Ziwei Qiu 5/16/2019
    """

    _DEFAULT_SETTINGS = Parameter([
        Parameter('VISA_address', 'USB0::0x1313::0x8078::P0010992::0::INSTR',
                  ['USB0::0x1313::0x8078::P0010992::0::INSTR'], 'VISA address of the instrument'),
        Parameter('gain_factor', 5.2, float, 'Type in the calibrated gain factor'),
        # Parameter('parameter', 'Power', ['Power', 'Power_dBm', 'Current', 'Irradiance'],
        #           'select the parameter to display'),
        # Parameter('resolution', 'Medium', ['Low', 'Medium', 'High'],'select the display resolution'),
        # Parameter('wavelength', '535nm', ['535nm', '633nm', '1064nm'], 'select the light wavelength'),
        # Parameter('bandwidth', 'High', ['Low', 'High'], 'select the detection bandwidth'),
        # Parameter('ranges', '1.4mW', ['1.4uW', '14uW', '140uW', '1.4mW', '14mW', '140mW'], 'select the detection bandwidth'),
        # Parameter('auto_range', False, bool, 'enable or disable the autorange'),
        # Parameter('attenuation', 0, int, 'select the attenutation in dB'),
    ])

    def __init__(self, name=None, settings=None):

        super(PM100D, self).__init__(name, settings)

        # XXXXX MW ISSUE = START
        # ===========================================
        # Issue where visa.ResourceManager() takes 4 minutes no longer happens after using pdb to debug (??? not sure why???)
        try:
            self._connect()
        except pyvisa.errors.VisaIOError:
            print('No power meter PM100D Detected!. Check that you are using the correct communication type')
            raise
        except Exception as e:
            print('error in __init__')
            raise (e)
        # Configure for power measurement
        self.srs.write('CONF:POW')
        self.get_power()

    def _connect(self, verbose = False):
        # print('_connect')
        rm = visa.ResourceManager()
        self.srs = rm.open_resource(self.settings['VISA_address'])
        # self.srs.query('*IDN?')
        if verbose:
            print('Power meter PM100D is connected: ' + self.srs.query('*IDN?'))

    def update(self, settings):
        """
        Updates the internal settings of the MicrowaveGenerator, and then also updates physical parameters such as
        frequency, amplitude, modulation type, etc in the hardware
        Args:
            settings: a dictionary in the standard settings format

        """

        super(PM100D, self).update(settings)

    def get_power(self):
        # power = self.srs.query('READ?')
        power = float(self.srs.query('READ?')) * self.settings['gain_factor'] * 1000 # convert to mW
        print('Gain factor is {:0.4f}'.format(self.settings['gain_factor']))
        print('Optical power is {:0.4f} mW'.format(power))

        return power

    @property

    def _PROBES(self):
        return None

    def read_probes(self, key):
        pass

    @property
    def is_connected(self):
        try:
            self.srs.query('*IDN?')  # arbitrary call to check connection, throws exception on failure to get response
            return True
        except pyvisa.errors.VisaIOError:
            print('error in is_connected')
            return False


    def __del__(self, verbose = True):
        # when done with device we have to close the connections
        # called when GUI is closed
        if verbose:
            print('PM100D is closing..')

        self.srs.close()


if __name__ == '__main__':

    power_meter = PM100D()
