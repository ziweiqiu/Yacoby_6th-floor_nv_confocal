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
import numpy
import time

from pylabcontrol.core import Parameter, Instrument

# RANGE_MIN = 2025000000 #2.025 GHz
# RANGE_MIN = -0.500 # V, minimum voltage for the SRS IQ
# RANGE_MAX = 0.500 #V, maximum power for the SRS IQ

class Agilent33120A(Instrument):
    """
    This class implements the Keysight AFG/AWG 33120A. The class commuicates with the
    device over GPIB using pyvisa.
    Attention that this instrument has trigger latency = 1.144us.

    -- Ziwei Qiu 6/25/2019
    """

    _DEFAULT_SETTINGS = Parameter([

        Parameter('VISA_address', 'GPIB1::10::INSTR', ['GPIB1::10::INSTR'], 'VISA address of the instrument'),
        Parameter('display_on', False, bool, 'Switches the update of the display on/off for fast performance'),
        Parameter('output_load', 'INFinity', ['INFinity'], 'open-circuit termination'),
        Parameter('trigger_latency', 1144.0, float,'trigger latency in ns'),
        Parameter('frequency', 4e6, float, 'Frequency in Hz. Do not exceed 5MHz, highest limit for the burst mode'),
        Parameter('amplitude', 5.0, float, 'Amplitude in Vp-p. This is the number assuming 50ohm termination'),
        Parameter('offset', 0.0, float, 'DC offset in Volt'),
        Parameter('wave_shape', 'SINusoid', ['SINusoid', 'SQUare', 'DC'], 'choose the wave shape'),
        Parameter('burst_mod', False, bool,'Turn on/off the burst modulation. External trigger is always used.'),
        Parameter('burst_count', 2, int, 'Sets the number of cycles per burst (1 to 50,000 cycles)'),
        Parameter('burst_phase', 2.0, float, 'Unit: degree. From -360 deg to 360 deg.'),

    ])

    def __init__(self, name=None, settings=None):

        super(Agilent33120A, self).__init__(name, settings)
        try:
            self._connect()
        except pyvisa.errors.VisaIOError:
            print('No AFG33120A Detected!. Check that you are using the correct communication type')
            raise
        except Exception as e:
            raise (e)

    def _connect(self, verbose = False):

        rm = visa.ResourceManager()
        self.afg = rm.open_resource(self.settings['VISA_address'])
        self.afg.query('*IDN?')
        if verbose:
            print('Agilent 33120A is connected: ' + self.afg.query('*IDN?'))
        time.sleep(0.1)
        self.afg.write('OUTPut:LOAD INF') # open-circuit termination

    def update(self, settings, verbose = False):
        """
        Updates the internal settings of the MicrowaveGenerator, and then also updates physical parameters such as
        frequency, amplitude, modulation type, etc in the hardware
        Args:
            settings: a dictionary in the standard settings format
        """
        super(Agilent33120A, self).update(settings)
        # ===========================================
        FREQ_MIN = 0.01 # 10mHz
        FREQ_MAX = 5000000  # 5MHz for burst
        AMP_MIN = 0  # 0Vpp
        AMP_MAX = 20  # 20Vpp
        CNT_MIN = 1 # minimum burst counts = 1
        CNT_MAX = 50000 # maximum burst counts = 50000
        PHASE_MIN = -360.0  # minimum burst phase = -360
        PHASE_MAX = 360.0  # maximum burst phase = 360

        for key, value in settings.items():
            if verbose:
                print('KEY: ',key)
                print('VALUE: ', value)

            if key == 'VISA_address':
                time.sleep(0.1)
                self._connect()

            elif not key == 'trigger_latency':
                if key == 'burst_mod':
                    if value:
                        time.sleep(0.1)
                        self.afg.write('TRIG:SOUR EXT')
                        time.sleep(0.1)
                        self.afg.write('BM:SOUR INT')
                    value = self._output_to_internal(value)
                if key == 'display_on':
                    value = self._output_to_internal(value)
                elif key == 'frequency':
                    if value > FREQ_MAX:
                        print('Invalid frequency. All frequencies must be between 10mHz and 5MHz. Set to 5MHz instead')
                        value = 5000000.0
                    elif value < FREQ_MIN:
                        print('Invalid frequency. All frequencies must be between 10mHz and 5MHz. Set to 10mHz instead')
                        value = 0.01
                elif key == 'amplitude' or key == 'offset':
                    if value > AMP_MAX:
                        print('Invalid amplitude/offset. All amplitudes/offsets must be between 0 and 20Vpp. Set to 20Vpp instead')
                        value = 20.0
                    elif value < AMP_MIN:
                        print('Invalid amplitude/offset. All amplitudes/offsets must be between 0 and 20Vpp. Set to 0Vpp instead')
                        value = 0.0
                elif key == 'burst_count':
                    if value > CNT_MAX:
                        print(
                            'Invalid burst counts. Burst counts must be between 1 and 50000 cycles. Set to 50000 instead')
                        value = 50000
                    if value < CNT_MIN:
                        print(
                            'Invalid burst counts. Burst counts must be between 1 and 50000 cycles. Set to 1 instead')
                        value = 1
                elif key == 'burst_phase':
                    if value > PHASE_MAX:
                        print('Invalid burst phase. Burst phases must be between -360 and 360 deg. Set to 360 instead')
                        value = 360.0
                    if value < PHASE_MIN:
                        print('Invalid burst phase. Burst phases must be between -360 and 360 deg. Set to -360 instead')
                        value = -360.0
                key = self._param_to_internal(key)
                if self._settings_initialized:
                    time.sleep(0.1)
                    self.afg.write(key + ' ' + str(value))

    @property
    def _PROBES(self):
        return{
        }

    def read_probes(self, key):

        # assert hasattr(self, 'awg') #will cause read_probes to fail if connection not yet established, such as when called in init
        assert(self._settings_initialized) #will cause read_probes to fail if settings (and thus also connection) not yet initialized
        assert key in list(self._PROBES.keys())

        #query always returns string, need to cast to proper return type

        if key == 'trigger_latency':
            value = True
        elif key == 'display_on' or key == 'burst_mod':
            key_internal = self._param_to_internal(key)
            value = int(self.srs.query(key_internal + '?'))
            if value == 1:
                value = True
            elif value == 0:
                value = False
        else:
            key_internal = self._param_to_internal(key)
            value = float(self.srs.query(key_internal + '?'))

        return value

    @property
    def is_connected(self):
        try:
            self.awg.query('*IDN?') # arbitrary call to check connection, throws exception on failure to get response
            return True
        except pyvisa.errors.VisaIOError:
            return False

    def _param_to_internal(self, param, key0=0):
        """
        Converts settings parameters to the corresponding key used for GPIB commands in the SRS.
        Args:
            param: settings parameter, ex. enable_output

        Returns: GPIB command, ex. ENBR

        """
        if param == 'display_on':
            return 'DISPlay'
        elif param == 'output_load':
            return 'OUTPut:LOAD'
        elif param == 'frequency':
            return 'FREQ'
        elif param == 'amplitude':
            return 'VOLT'
        elif param == 'offset':
            return 'VOLT:OFFS'
        elif param == 'wave_shape':
            return 'FUNC:SHAP'
        elif param == 'burst_mod':
            return 'BM:STAT'
        elif param == 'burst_count':
            return 'BM:NCYC'
        elif param == 'burst_phase':
            return 'BM:PHAS'
        else:
            print('KeyError with param = ', param)
            raise KeyError

    def _output_to_internal(self, value):
        if value == True:
            return 'ON'
        elif value == False:
            return 'OFF'
        else:
            print('Agilent33120A _output_to_internal KeyError')
            raise KeyError

if __name__ == '__main__':
    arb = Agilent33120A()
