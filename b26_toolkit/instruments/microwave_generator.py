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


class AgilentMicrowaveGenerator(Instrument):
    """
        This class implements the Agilent microwave generator N9310A. The class commuicates with the device over USB.
        manual: https://www.keysight.com/upload/cmc_upload/All/N9310AUsersGuide.pdf

        -Ziwei (1/31/2019 4:05pm)
    """

    _DEFAULT_SETTINGS = Parameter([
        Parameter('VISA_address', 'USB0::0x0957::0x2018::01152687::0::INSTR',
                  ['USB0::0x0957::0x2018::01152687::0::INSTR'], 'VISA address of the instrument'),
        Parameter('reference_oscillator', 'EXT10MHZ', ['INT10MHZ', 'EXT2MHZ', 'EXT5MHZ', 'EXT10MHZ'],
                  'choose the reference oscillator'),
        Parameter('enable_output', False, bool, 'Type-N RF output enabled'),
        Parameter('enable_modulation', False, bool, 'enable modulation'),
        Parameter('enable_IQ', False, bool, 'enable IQ modulation (only IQ is implemented for now)'),
        # Parameter('modulation_type', 'IQ', ['IQ'],'choose the modulation type (only IQ is implemented for now)'),
        Parameter('freq_mode', 'CW', ['CW', 'Sweep'], 'select the frequency mode'),
        Parameter('power_mode', 'CW', ['CW', 'Sweep'], 'select the power mode'),
        Parameter('edge', 'Positive', ['Positive', 'Negative'], 'select the external triggerring edge'),
        Parameter('swp_direction', 'UP', ['UP', 'DOWN'], 'select the sweep direction'),
        Parameter('frequency', 2.87e9, float, 'RF frequency in Hz, range: 9kHz to 3 GHz'),
        Parameter('freq_start', 100e6, float, 'start frequency in Hz in sweep mode'),
        Parameter('freq_stop', 400e6, float, 'stop frequency in Hz in sweep mode'),
        Parameter('freq_pts', 100, float, 'number of sweep steps in freq sweep mode'),
        Parameter('amplitude', -50, float, 'RF Type-N power in dBm, range: -127 to +20dBm'),
        Parameter('pwr_start', -20, float, 'start power in dBm in sweep mode'),
        Parameter('pwr_stop', 0, float, 'stop power in dBm in sweep mode'),
        Parameter('pwr_pts', 20, float, 'number of sweep steps in power sweep mode'),
        Parameter('LF_enable_output', False, bool, 'BNC LF output enabled'),
        Parameter('LF_frequency', 100, float, 'LF frequency in Hz, range: 20Hz - 80kHz'),
        Parameter('LF_amplitude', 0, float, 'LF output amplitude in V, range: 0 to 3Vp')
    ])

    def __init__(self, name=None, settings=None):

        super(AgilentMicrowaveGenerator, self).__init__(name, settings)

        # XXXXX MW ISSUE = START
        # ===========================================
        # Issue where visa.ResourceManager() takes 4 minutes no longer happens after using pdb to debug (??? not sure why???)
        try:
            self._connect()
        except pyvisa.errors.VisaIOError:
            print('No Agilent Microwave Generator Detected!. Check that you are using the correct communication type')
            raise
        except Exception as e:
            print('error in __init__')
            raise (e)

    def _connect(self, verbose = False):
        # print('_connect')
        rm = visa.ResourceManager()
        self.srs = rm.open_resource(self.settings['VISA_address'])
        # self.srs.query('*IDN?')
        if verbose:
            print('Agilent N9310A RF generator connected: ' + self.srs.query('*IDN?'))

    def update(self, settings):
        """
        Updates the internal settings of the MicrowaveGenerator, and then also updates physical parameters such as
        frequency, amplitude, modulation type, etc in the hardware
        Args:
            settings: a dictionary in the standard settings format

        """

        super(AgilentMicrowaveGenerator, self).update(settings)
        # XXXXX MW ISSUE = START
        # ===========================================
        RF_RANGE_MIN = 9000  # 9 kHz
        RF_RANGE_MAX = 3000000000  # 3 GHZ
        LF_RANGE_MIN = 20  # 20Hz
        LF_RANGE_MAX = 80000  # 80 kHz

        RF_PWR_MIN = -127  # -127dBm
        RF_PWR_MAX = 20  # 20dBm
        LF_AMP_MIN = 0  # 0Vp
        LF_AMP_MAX = 3  # 3Vp

        STEP_PTS_MIN = 2
        STEP_PTS_MAX = 1001

        for key, value in settings.items():
            # print (key + ' is updated')
            # print (settings)
            # print(' ')
            error_flag = False
            # print('_update')
            if key == 'VISA_address':
                time.sleep(0.1)
                self._connect()
            elif key == 'reference_oscillator':
                time.sleep(0.1)
                self.srs.write(':SYSTem:REFerence:FREQuency' + ' ' + str(value))
                print('Agilent MW generator current reference oscillator :' + self.srs.query(':SYSTem:REFerence:FREQuency?'))
            else:
                if key in ['enable_output', 'enable_modulation', 'enable_IQ','LF_enable_output', 'edge']:
                    value = self._output_to_internal(value)
                elif key in ['frequency', 'freq_start', 'freq_stop']:
                    if value > RF_RANGE_MAX or value < RF_RANGE_MIN:
                        error_flag = True
                        print('Invalid RF frequency. RF frequencies must be between 9kHz and 3GHz.')
                    else:
                        unit = 'Hz'
                elif key in ['amplitude', 'pwr_start', 'pwr_stop']:
                    if value > RF_PWR_MAX or value < RF_PWR_MIN:
                        error_flag = True
                        print('Invalid RF power. RF power must be between -127 and 20dBm.')
                    else:
                        unit = 'dBm'
                elif key == 'LF_frequency':
                    if value > LF_RANGE_MAX or value < LF_RANGE_MIN:
                        error_flag = True
                        print('Invalid LF frequency. LF frequencies must be between 20Hz and 80kHz.')
                        # raise ValueError("Invalid frequency. RF frequencies must be between 20Hz and 80kHz.")
                    else:
                        unit = 'Hz'
                elif key == 'LF_amplitude':
                    if value > LF_AMP_MAX or value < LF_AMP_MIN:
                        error_flag = True
                        print('Invalid LF amplitude. LF amplitude must be between 0 and 3Vp.')
                        # raise ValueError("Invalid frequency. RF frequencies must be between 20Hz and 80kHz.")
                    else:
                        unit = 'V'
                elif key == 'freq_pts' or key == 'pwr_pts':
                    if value > STEP_PTS_MAX or value < STEP_PTS_MIN:
                        error_flag = True
                        print('Invalid step points. RF frequencies must be between' + str(STEP_PTS_MIN) + 'and' + str(STEP_PTS_MAX))
                elif key == 'freq_mode':
                    time.sleep(0.1)
                    self.srs.write(':SWEep:REPeat CONTinuous')
                    time.sleep(0.1)
                    self.srs.write(':SWEep:TYPE STEP')
                    time.sleep(0.1)
                    self.srs.write(':FREQuency:RF:SCALe LIN')
                    time.sleep(0.1)
                    self.srs.write(':SWEep:STRG IMMediate')
                    time.sleep(0.1)
                    self.srs.write(':SWEep:PTRG EXT')
                    value = self._output_to_internal(value)
                elif key == 'power_mode':
                    # The LF sweep and amplitude sweep scale is default to linear scale only
                    time.sleep(0.1)
                    self.srs.write(':SWEep:REPeat CONTinuous')
                    time.sleep(0.1)
                    self.srs.write(':SWEep:TYPE STEP')
                    time.sleep(0.1)
                    self.srs.write(':SWEep:STRG IMMediate')
                    time.sleep(0.1)
                    self.srs.write(':SWEep:PTRG EXT')
                    value = self._output_to_internal(value)

                if not error_flag:
                    internal_key = self._param_to_internal(key)
                    if self._settings_initialized:
                        if key in ['frequency', 'freq_start', 'freq_stop', 'amplitude', 'pwr_start', 'pwr_stop','LF_frequency', 'LF_amplitude']:
                            command = internal_key + ' ' + str(value) + ' ' + str(unit)
                            time.sleep(0.1)
                            self.srs.write(command)
                            # print('Command: ', command)
                        else:
                            command = internal_key + ' ' + str(value)
                            time.sleep(0.1)
                            self.srs.write(command)
                            # print('Command: ', command)
                else:
                    print('Invalid input - no updating on Agilent N9310A MW generator')

    @property
    def _PROBES(self):
        return {
            'enable_output': 'if type-N output is enabled',
            'frequency': 'frequency of output in Hz',
            'amplitude': 'RF type-N amplitude in dBm',
            'LF_enable_output': 'if type-N output is enabled',
            'LF_frequency': 'frequency of output in Hz',
            'LF_amplitude': 'LF amplitude in V'
        }

    def read_probes(self, key):
        # assert hasattr(self, 'srs') #will cause read_probes to fail if connection not yet established, such as when called in init
        assert (
            self._settings_initialized)  # will cause read_probes to fail if settings (and thus also connection) not yet initialized
        assert key in self._PROBES.keys()

        # query always returns string, need to cast to proper return type
        # if key in ['enable_output', 'enable_modulation']:
        if key == 'enable_output' or key == 'LF_enable_output':
            key_internal = self._param_to_internal(key)
            value = int(self.srs.query(key_internal + '?'))
            if value == 1:
                value = True
            elif value == 0:
                value = False
        # elif key in ['modulation_type','modulation_function','pulse_modulation_function']:
        #     key_internal = self._param_to_internal(key)
        #     value = int(self.srs.query(key_internal + '?'))
        #     if key == 'modulation_type':
        #         value = self._internal_to_mod_type(value)
        #     elif key == 'modulation_function':
        #         value = self._internal_to_mod_func(value)
        #     elif key == 'pulse_modulation_function':
        #         value = self._internal_to_pulse_mod_func(value)
        else:
            key_internal = self._param_to_internal(key)
            value = float(self.srs.query(key_internal + '?'))

        return value

    @property
    def is_connected(self):
        try:
            self.srs.query('*IDN?')  # arbitrary call to check connection, throws exception on failure to get response
            return True
        except pyvisa.errors.VisaIOError:
            print('error in is_connected')
            return False

    def reset_sweep(self):

        # the following command effectively resets the frequency sweep from the initial value
        self.srs.write(':SWEep:RF:STATe ON')
        self.srs.write(':SWEep:STRG IMMediate')
        self.srs.write(':SWEep:PTRG EXT')

    def _param_to_internal(self, param):
        """
        Converts settings parameters to the corresponding key used for GPIB commands in the SRS.
        Args:
            param: settings parameter, ex. enable_output

        Returns: GPIB command, ex. ENBR

        """
        if param == 'enable_output':
            return ':RFOutput:STATE'
        elif param == 'enable_modulation':
            return ':MOD:STATe'
        elif param == 'enable_IQ':
            return ':IQ:STATe'
        elif param == 'freq_mode':
            return ':SWEep:RF:STATe'
        elif param == 'power_mode':
            return ':SWEep:AMPLitude:STATe'
        elif param == 'edge':
            return ':SWEep:STRG:SLOPe'
        elif param == 'swp_direction':
            return ':SWEep:DIRection'
        elif param == 'frequency':
            return ':FREQ:CW'
        elif param == 'freq_start':
            return ':SWEep:RF:STARt'
        elif param == 'freq_stop':
            return ':SWEep:RF:STOP'
        elif param == 'freq_pts':
            return ':SWEep:STEP:POINts'
        elif param == 'amplitude':
            return ':AMPLITUDE:CW'
        elif param == 'pwr_start':
            return ':SWEep:AMPLitude:STARt'
        elif param == 'pwr_stop':
            return ':SWEep:AMPLitude:STOP'
        elif param == 'pwr_pts':
            return ':SWEep:STEP:POINts' # same command as freq points
        elif param == 'LF_enable_output':
            return ':LFOutput:STATE'
        elif param == 'LF_frequency':
            return ':LFOutput:FREQ'
        elif param == 'LF_amplitude':
            return ':LFOutput:AMPL'
        else:
            print('Error in _param_to_internal')
            raise KeyError

    def _output_to_internal(self, value):
        if value == True or value == 'Sweep':
            return 'ON'
        elif value == False or value == 'CW':
            return 'OFF'
        elif value == 'Positive':
            return 'EXTP'
        elif value == 'Negative':
            return 'EXTN'
        else:
            print('error in _output_to_internal')
            raise KeyError

    def __del__(self, verbose = True):
        # when done with device we have to close the connections
        # called when GUI is closed
        if verbose:
            print('N9310A is closing..')

        self.srs.close()

class AgilentMicrowaveGeneratorII(AgilentMicrowaveGenerator):
    """
        This class implements the Agilent microwave generator N9310A. The class commuicates with the device over USB.
        manual: https://www.keysight.com/upload/cmc_upload/All/N9310AUsersGuide.pdf

        -Ziwei (2/21/2019 4:05pm)
    """

    _DEFAULT_SETTINGS = Parameter([
        Parameter('VISA_address', 'USB0::0x0957::0x2018::01152673::0::INSTR',
                  ['USB0::0x0957::0x2018::01152673::0::INSTR'], 'VISA address of the instrument'),
        Parameter('reference_oscillator', 'EXT10MHZ', ['INT10MHZ', 'EXT2MHZ', 'EXT5MHZ', 'EXT10MHZ'],
                  'choose the reference oscillator'),
        Parameter('enable_output', False, bool, 'Type-N RF output enabled'),
        # Parameter('enable_modulation', False, bool, 'enable modulation'),
        # Parameter('enable_IQ', False, bool, 'enable IQ modulation (only IQ is implemented for now)'),
        # Parameter('modulation_type', 'IQ', ['IQ'],'choose the modulation type (only IQ is implemented for now)'),
        Parameter('freq_mode', 'CW', ['CW', 'Sweep'], 'select the frequency mode'),
        Parameter('power_mode', 'CW', ['CW', 'Sweep'], 'select the power mode'),
        Parameter('edge', 'Positive', ['Positive', 'Negative'], 'select the external triggerring edge'),
        Parameter('swp_direction', 'UP', ['UP', 'DOWN'], 'select the sweep direction'),
        Parameter('frequency', 2.87e9, float, 'RF frequency in Hz, range: 9kHz to 3 GHz'),
        Parameter('freq_start', 100e6, float, 'start frequency in Hz in sweep mode'),
        Parameter('freq_stop', 400e6, float, 'stop frequency in Hz in sweep mode'),
        Parameter('freq_pts', 100, float, 'number of sweep steps in freq sweep mode'),
        Parameter('amplitude', -50, float, 'RF Type-N power in dBm, range: -127 to +20dBm'),
        Parameter('pwr_start', -20, float, 'start power in dBm in sweep mode'),
        Parameter('pwr_stop', 0, float, 'stop power in dBm in sweep mode'),
        Parameter('pwr_pts', 20, float, 'number of sweep steps in power sweep mode'),
        Parameter('LF_enable_output', False, bool, 'BNC LF output enabled'),
        Parameter('LF_frequency', 100, float, 'LF frequency in Hz, range: 20Hz - 80kHz'),
        Parameter('LF_amplitude', 0, float, 'LF output amplitude in V, range: 0 to 3Vp')

    ])

class R8SMicrowaveGenerator(Instrument):
    """
    This class implements the ROHDE & SCHWARZ microwave generator SMB100A. The class commuicates with the
    device over USB or GPIB using pyvisa.

    -Ziwei (12/20/2018 2:30pm)
    """

    _DEFAULT_SETTINGS = Parameter([
        Parameter('VISA_address','GPIB0::29::INSTR',['GPIB0::29::INSTR'],'VISA address of the instrument'),
        Parameter('reference_oscillator', 'INT10MHZ', ['INT10MHZ', 'EXT5MHZ', 'EXT10MHZ'],
                  'choose the reference oscillator'),
        Parameter('enable_output', False, bool, 'Type-N output enabled'),
        Parameter('display_on', True, bool,'Switches the update of the display on/off.'),
        Parameter('freq_mode','CW',['CW','Sweep'],'select the frequency mode'),
        Parameter('power_mode', 'CW', ['CW','Sweep'], 'select the power mode'),
        Parameter('edge', 'Positive', ['Positive', 'Negative'], 'select the triggerring edge'),
        Parameter('frequency', 252e6, float, 'frequency in Hz, or with label in other units ex 300 MHz'),
        Parameter('freq_start', 100e6, float, 'start frequency in Hz in sweep mode'),
        Parameter('freq_stop', 400e6, float, 'stop frequency in Hz in sweep mode'),
        Parameter('freq_pts', 100, float, 'number of sweep steps in freq sweep mode'),
        Parameter('amplitude', -45, float, 'Type-N amplitude in dBm'),
        Parameter('pwr_start', -20, float, 'start power in dBm in sweep mode'),
        Parameter('pwr_stop', 0, float, 'stop power in dBm in sweep mode'),
        Parameter('pwr_pts', 20, float, 'number of sweep steps in power sweep mode')
    ])

    def __init__(self, name=None, settings=None):

        super(R8SMicrowaveGenerator, self).__init__(name, settings)

        # XXXXX MW ISSUE = START
        #===========================================
        # Issue where visa.ResourceManager() takes 4 minutes no longer happens after using pdb to debug (??? not sure why???)
        try:
            self._connect()
        except pyvisa.errors.VisaIOError:
            print('No Microwave Controller Detected!. Check that you are using the correct communication type')
            raise
        except Exception as e:
            raise(e)
        #XXXXX MW ISSUE = END
        #===========================================

    def _connect(self, verbose = False):
        rm = visa.ResourceManager()
        # self.srs = rm.open_resource('USB0::0x0AAD::0x006E::102953::INSTR')
        self.srs = rm.open_resource(self.settings['VISA_address'])
        self.srs.query('*IDN?')
        if verbose:
            print('Rohde&Schwarz SMB100A connected: '+  self.srs.query('*IDN?'))
    #Doesn't appear to be necessary, can't manually make two sessions conflict, rms may share well
    # def __del__(self):
    #     self.srs.close()
    def update(self, settings):
        """
        Updates the internal settings of the MicrowaveGenerator, and then also updates physical parameters such as
        frequency, amplitude, modulation type, etc in the hardware
        Args:
            settings: a dictionary in the standard settings format

        """
        # print('Rhode Schwartz start updating')
        super(R8SMicrowaveGenerator, self).update(settings)
        # XXXXX MW ISSUE = START
        # ===========================================
        RANGE_MIN = 9000
        RANGE_MAX = 3200000000  # 3.2 GHZ
        if self._settings_initialized:
            for key, value in settings.items():
                # print (key + ' is updated')
                # print (settings)
                # print(' ')
                if key == 'VISA_address':

                    self._connect()
                elif key == 'reference_oscillator':
                    # print('updating reference_oscillator')
                    if value == 'INT10MHZ':
                        self.srs.write('ROSC:SOUR INT')
                    elif value == 'EXT10MHZ':
                        self.srs.write('ROSC:SOUR EXT')
                        self.srs.write('ROSC:EXT:FREQ 10 MHz')
                    else:
                        self.srs.write('ROSC:SOUR EXT')
                        self.srs.write('ROSC:EXT:FREQ 5 MHz')
                    # print('R&S MW generator current reference oscillator :' + self.srs.query(
                    #     'ROSC:SOUR?' + ' ' + self.srs.query('ROSC:EXT:FREQ?')))
                    print('R&S MW generator current reference oscillator :' + self.srs.query(
                        'ROSC:SOUR?') + self.srs.query('ROSC:EXT:FREQ?'))


                else:
                    # print('current key', key)
                    if key == 'enable_output' or key == 'display_on':
                        value = self._output_to_internal(value)
                    elif key == 'frequency':
                        if value > RANGE_MAX or value < RANGE_MIN:
                            print("Invalid frequency. All frequencies must be between 9kHz and 3.2 GHz.")

                    elif key == 'freq_mode':
                        self.srs.write('SOUR:SWE:FREQ:MODE STEP')
                        self.srs.write('SOUR:SWE:FREQ:SPAC LIN')
                        self.srs.write('TRIG:FSW:SOUR EXT')
                    elif key == 'power_mode':
                        self.srs.write('SOUR:SWE:POW:MODE STEP')
                        self.srs.write('TRIG:PSW:SOUR EXT')
                    key = self._param_to_internal(key)
                    # if self._settings_initialized:
                    self.srs.write(key + ' ' + str(value))
                    # print(key + ' ' + str(value))

        # print('Rhode Schwartz update done')




                # only send update to instrument if connection to instrument has been established
                 # frequency change operation timed using timeit.timeit and
                                                           # completion confirmed by query('*OPC?'), found delay of <10ms
                    # print(self.srs.query('*OPC?'))

        # XXXXX MW ISSUE = END
        # ===========================================

    @property
    def _PROBES(self):
        # return{
        #     'enable_output': 'if type-N output is enabled',
        #     'frequency': 'frequency of output in Hz',
        #     'amplitude': 'type-N amplitude in dBm',
        #     'phase': 'phase',
        #     'enable_modulation': 'is modulation enabled',
        #     'modulation_type': 'Modulation Type: 0= AM, 1=FM, 2= PhaseM, 3= Freq sweep, 4= Pulse, 5 = Blank, 6=IQ',
        #     'modulation_function': 'Modulation Function: 0=Sine, 1=Ramp, 2=Triangle, 3=Square, 4=Noise, 5=External',
        #     'pulse_modulation_function': 'Pulse Modulation Function: 3=Square, 4=Noise(PRBS), 5=External',
        #     'dev_width': 'Width of deviation from center frequency in FM'
        # }

        return {
            'enable_output': 'if type-N output is enabled',
            'frequency': 'frequency of output in Hz',
            'amplitude': 'type-N amplitude in dBm'
        }

    def read_probes(self, key):
        # assert hasattr(self, 'srs') #will cause read_probes to fail if connection not yet established, such as when called in init
        assert(self._settings_initialized) #will cause read_probes to fail if settings (and thus also connection) not yet initialized
        assert key in self._PROBES.keys()

        #query always returns string, need to cast to proper return type
        # if key in ['enable_output', 'enable_modulation']:
        if key == 'enable_output':
            key_internal = self._param_to_internal(key)
            value = int(self.srs.query(key_internal + '?'))
            if value == 1:
                value = True
            elif value == 0:
                value = False
        # elif key in ['modulation_type','modulation_function','pulse_modulation_function']:
        #     key_internal = self._param_to_internal(key)
        #     value = int(self.srs.query(key_internal + '?'))
        #     if key == 'modulation_type':
        #         value = self._internal_to_mod_type(value)
        #     elif key == 'modulation_function':
        #         value = self._internal_to_mod_func(value)
        #     elif key == 'pulse_modulation_function':
        #         value = self._internal_to_pulse_mod_func(value)
        else:
            key_internal = self._param_to_internal(key)
            value = float(self.srs.query(key_internal + '?'))

        return value

    @property
    def is_connected(self):
        try:
            self.srs.query('*IDN?') # arbitrary call to check connection, throws exception on failure to get response
            return True
        except pyvisa.errors.VisaIOError:
            return False

    def reset_sweep(self):
        # print('resetting sweep')
        self.srs.write('SWE:RES')

    def _param_to_internal(self, param):
        """
        Converts settings parameters to the corresponding key used for GPIB commands in the SRS.
        Args:
            param: settings parameter, ex. enable_output

        Returns: GPIB command, ex. ENBR

        """

        if param == 'enable_output':
            return 'OUTP'
        if param == 'display_on':
            return 'SYST:DISP:UPD'
        elif param == 'freq_mode':
            return 'SOUR:FREQ:MODE'
        elif param == 'power_mode':
            return 'SOUR:POW:MODE'
        elif param == 'edge':
            return 'INP:TRIG:SLOP'
        elif param == 'frequency':
            return 'FREQ'
        elif param == 'freq_start':
            return 'FREQ:STAR'
        elif param == 'freq_stop':
            return 'FREQ:STOP'
        elif param == 'freq_pts':
            return 'SWE:POIN'
        elif param == 'amplitude':
            return 'SOUR:POW:POW'
        elif param == 'pwr_start':
            return 'POW:STAR'
        elif param == 'pwr_stop':
            return 'POW:STOP'
        elif param == 'pwr_pts':
            return 'SWE:POW:POIN'
        # elif param == 'phase':
        #     return 'PHAS'
        # elif param == 'enable_modulation':
        #     return 'MODL'
        # elif param == 'modulation_type':
        #     return 'TYPE'
        # elif param == 'modulation_function':
        #     return 'MFNC'
        # elif param == 'pulse_modulation_function':
        #     return 'PFNC'
        # elif param == 'dev_width':
        #     return 'FDEV'
        else:

            raise KeyError

    def _output_to_internal(self, value):
        if value == True:
            return 'ON'
        elif value == False:
            return 'OFF'
        else:
            print('Rhode&Schwartz MW generator _output_to_internal KeyError')
            raise KeyError

# the following is not working don't use it -ZQ 3/2/2019
class RnSMicrowaveGenerator(Instrument):
    """
    This class implements the ROHDE & SCHWARZ microwave generator SMB100A. The class commuicates with the
    device over USB or GPIB using pyvisa.
    This class has more settings on the frequency sweep.

    -Ziwei (12/20/2018 2:30pm)
    """

    _DEFAULT_SETTINGS = Parameter([
        Parameter('VISA_address','GPIB0::29::INSTR',['GPIB0::29::INSTR'],'VISA address of the instrument'),
        Parameter('enable_output', False, bool, 'Type-N output enabled'),
        Parameter('freq_mode','CW',['CW','Sweep'],'select the frequency mode'),
        Parameter('power_mode', 'CW', ['CW','Sweep'], 'select the power mode'),
        Parameter('edge', 'Positive', ['Positive', 'Negative'], 'select the triggerring edge'),
        Parameter('frequency', 252e6, float, 'frequency in Hz, or with label in other units ex 300 MHz'),
        Parameter('amplitude', -45, float, 'Type-N amplitude in dBm'),
        Parameter('freq_sweep', [Parameter('freq_sweep_state', False, bool, 'turn on frequency sweep'),
                                 Parameter('freq_start', 100e6, float, 'start frequency in Hz in sweep mode'),
                                 Parameter('freq_stop', 400e6, float, 'stop frequency in Hz in sweep mode'),
                                 Parameter('freq_pts', 100, float, 'number of sweep steps in freq sweep mode')]),
        Parameter('pwr_sweep', [Parameter('pwr_sweep_state', False, bool, 'turn on power sweep'),
                                Parameter('pwr_start', -20, float, 'start power in dBm in sweep mode'),
                                Parameter('pwr_stop', 0, float, 'stop power in dBm in sweep mode'),
                                Parameter('pwr_pts', 20, float, 'number of sweep steps in power sweep mode')]),
    ])

    def __init__(self, name=None, settings=None):
        super(RnSMicrowaveGenerator, self).__init__(name, settings)
        # XXXXX MW ISSUE = START
        #===========================================
        # Issue where visa.ResourceManager() takes 4 minutes no longer happens after using pdb to debug (??? not sure why???)
        try:
            self._connect()
        except pyvisa.errors.VisaIOError:
            print('No Microwave Controller Detected!. Check that you are using the correct communication type')
            raise
        except Exception as e:
            raise(e)
        #XXXXX MW ISSUE = END
        #===========================================

    def _connect(self):
        rm = visa.ResourceManager()
        # self.srs = rm.open_resource('USB0::0x0AAD::0x006E::102953::INSTR')

        self.srs = rm.open_resource(self.settings['VISA_address'])
        self.srs.query('*IDN?')
        print('Rohde&Schwarz SMB100A connected: '+  self.srs.query('*IDN?'))

    #Doesn't appear to be necessary, can't manually make two sessions conflict, rms may share well
    # def __del__(self):
    #     self.srs.close()

    def update(self, settings):
        """
        Updates the internal settings of the MicrowaveGenerator, and then also updates physical parameters such as
        frequency, amplitude, modulation type, etc in the hardware
        Args:
            settings: a dictionary in the standard settings format

        """
        # print('Rhode Schwartz start updating')
        super(RnSMicrowaveGenerator, self).update(settings)
        # XXXXX MW ISSUE = START
        # ===========================================
        RANGE_MIN = 9000
        RANGE_MAX = 3200000000  # 3.2 GHZ
        for key, value in settings.items():
            # print (key + ' is updated')
            # print (settings)
            # print(' ')
            if key == 'VISA_address':
                self._connect()
            else:
                if key == 'enable_output':
                    value = self._output_to_internal(value)
                elif key == 'frequency':
                    if value > RANGE_MAX or value < RANGE_MIN:
                        print('Invalid frequency. All frequencies must be between 9kHz and 3.2 GHz.')
                        # raise ValueError("Invalid frequency. All frequencies must be between 9kHz and 3.2 GHz.")
                elif key == 'freq_mode':
                    self.srs.write('SOUR:SWE:FREQ:MODE STEP')
                    self.srs.write('SOUR:SWE:FREQ:SPAC LIN')
                    self.srs.write('TRIG:FSW:SOUR EXT')
                elif key == 'power_mode':
                    self.srs.write('SOUR:SWE:POW:MODE STEP')
                    self.srs.write('TRIG:PSW:SOUR EXT')
                key = self._param_to_internal(key)
                if self._settings_initialized:
                    self.srs.write(key + ' ' + str(value))

                # only send update to instrument if connection to instrument has been established
                 # frequency change operation timed using timeit.timeit and
                                                           # completion confirmed by query('*OPC?'), found delay of <10ms
                    # print(self.srs.query('*OPC?'))



    @property
    def _PROBES(self):


        return {
            'enable_output': 'if type-N output is enabled',
            'frequency': 'frequency of output in Hz',
            'amplitude': 'type-N amplitude in dBm'
        }

    def read_probes(self, key):
        # assert hasattr(self, 'srs') #will cause read_probes to fail if connection not yet established, such as when called in init
        assert(self._settings_initialized) #will cause read_probes to fail if settings (and thus also connection) not yet initialized
        assert key in self._PROBES.keys()

        #query always returns string, need to cast to proper return type
        # if key in ['enable_output', 'enable_modulation']:
        if key == 'enable_output':
            key_internal = self._param_to_internal(key)
            value = int(self.srs.query(key_internal + '?'))
            if value == 1:
                value = True
            elif value == 0:
                value = False
        # elif key in ['modulation_type','modulation_function','pulse_modulation_function']:
        #     key_internal = self._param_to_internal(key)
        #     value = int(self.srs.query(key_internal + '?'))
        #     if key == 'modulation_type':
        #         value = self._internal_to_mod_type(value)
        #     elif key == 'modulation_function':
        #         value = self._internal_to_mod_func(value)
        #     elif key == 'pulse_modulation_function':
        #         value = self._internal_to_pulse_mod_func(value)
        else:
            key_internal = self._param_to_internal(key)
            value = float(self.srs.query(key_internal + '?'))

        return value

    @property
    def is_connected(self):
        try:
            self.srs.query('*IDN?') # arbitrary call to check connection, throws exception on failure to get response
            return True
        except pyvisa.errors.VisaIOError:
            return False

    def _param_to_internal(self, param):
        """
        Converts settings parameters to the corresponding key used for GPIB commands in the SRS.
        Args:
            param: settings parameter, ex. enable_output

        Returns: GPIB command, ex. ENBR

        """
        if param == 'enable_output':
            return 'OUTP'
        elif param == 'freq_mode':
            return 'SOUR:FREQ:MODE'
        elif param == 'power_mode':
            return 'SOUR:POW:MODE'
        elif param == 'edge':
            return 'INP:TRIG:SLOP'
        elif param == 'frequency':
            return 'FREQ'
        elif param == 'freq_start':
            return 'FREQ:STAR'
        elif param == 'freq_stop':
            return 'FREQ:STOP'
        elif param == 'freq_pts':
            return 'SWE:POIN'
        elif param == 'amplitude':
            return 'SOUR:POW:POW'
        elif param == 'pwr_start':
            return 'POW:STAR'
        elif param == 'pwr_stop':
            return 'POW:STOP'
        elif param == 'pwr_pts':
            return 'SWE:POW:POIN'
        # elif param == 'phase':
        #     return 'PHAS'
        # elif param == 'enable_modulation':
        #     return 'MODL'
        # elif param == 'modulation_type':
        #     return 'TYPE'
        # elif param == 'modulation_function':
        #     return 'MFNC'
        # elif param == 'pulse_modulation_function':
        #     return 'PFNC'
        # elif param == 'dev_width':
        #     return 'FDEV'
        else:

            raise KeyError

    def _output_to_internal(self, value):
        if value == True:
            return 'ON'
        elif value == False:
            return 'OFF'
        else:

            raise KeyError

# the following one might not be working
class MicrowaveGenerator(Instrument):
    """
    This class implements the Stanford Research Systems SG384 microwave generator. The class commuicates with the
    device over GPIB using pyvisa.
    """

    _DEFAULT_SETTINGS = Parameter([
        Parameter('connection_type', 'RS232', ['GPIB', 'RS232'], 'type of connection to open to controller'),
        Parameter('port', 8, list(range(0, 31)), 'GPIB or COM port on which to connect'), ## JG: what out for the ports this might be different on each computer and might cause issues when running export default
        Parameter('GPIB_num', 0, int, 'GPIB device on which to connect'),
        Parameter('enable_output', False, bool, 'Type-N output enabled'),
        Parameter('frequency', 3e9, float, 'frequency in Hz, or with label in other units ex 300 MHz'),
        Parameter('amplitude', -60, float, 'Type-N amplitude in dBm'),
        Parameter('phase', 0, float, 'output phase'),
        Parameter('enable_modulation', True, bool, 'enable modulation'),
        Parameter('modulation_type', 'FM', ['AM', 'FM', 'PhaseM', 'Freq sweep', 'Pulse', 'Blank', 'IQ'],
                  'Modulation Type: 0= AM, 1=FM, 2= PhaseM, 3= Freq sweep, 4= Pulse, 5 = Blank, 6=IQ'),
        Parameter('modulation_function', 'External', ['Sine', 'Ramp', 'Triangle', 'Square', 'Noise', 'External'],
                  'Modulation Function: 0=Sine, 1=Ramp, 2=Triangle, 3=Square, 4=Noise, 5=External'),
        Parameter('pulse_modulation_function', 'External', ['Square', 'Noise(PRBS)', 'External'], 'Pulse Modulation Function: 3=Square, 4=Noise(PRBS), 5=External'),
        Parameter('dev_width', 32e6, float, 'Width of deviation from center frequency in FM')
    ])

    def __init__(self, name=None, settings=None):

        super(MicrowaveGenerator, self).__init__(name, settings)

        # XXXXX MW ISSUE = START
        #===========================================
        # Issue where visa.ResourceManager() takes 4 minutes no longer happens after using pdb to debug (??? not sure why???)
        try:
            self._connect()
        except pyvisa.errors.VisaIOError:
            print('No Microwave Controller Detected!. Check that you are using the correct communication type')
            raise
        except Exception as e:
            raise(e)
        #XXXXX MW ISSUE = END
        #===========================================

    def _connect(self):
        rm = visa.ResourceManager()
        if self.settings['connection_type'] == 'GPIB':
            self.srs = rm.open_resource(
                'GPIB' + str(self.settings['GPIB_num']) + '::' + str(self.settings['port']) + '::INSTR')
        elif self.settings['connection_type'] == 'RS232':
            self.srs = rm.open_resource('COM' + str(self.settings['port']))
            self.srs.baud_rate = 115200
        self.srs.query('*IDN?')

    #Doesn't appear to be necessary, can't manually make two sessions conflict, rms may share well
    # def __del__(self):
    #     self.srs.close()

    def update(self, settings):
        """
        Updates the internal settings of the MicrowaveGenerator, and then also updates physical parameters such as
        frequency, amplitude, modulation type, etc in the hardware
        Args:
            settings: a dictionary in the standard settings format
        """
        super(MicrowaveGenerator, self).update(settings)
        # XXXXX MW ISSUE = START
        # ===========================================
        RANGE_MIN = 1012500000
        RANGE_MAX = 4050000000  # 4.050 GHZ
        for key, value in settings.items():
            if key == 'connection_type':
                self._connect()
            elif not (key == 'port' or key == 'GPIB_num'):
                if self.settings.valid_values[key] == bool: #converts booleans, which are more natural to store for on/off, to
                    value = int(value)                #the integers used internally in the SRS
                elif key == 'modulation_type':
                    value = self._mod_type_to_internal(value)
                elif key == 'modulation_function':
                    value = self._mod_func_to_internal(value)
                elif key == 'pulse_modulation_function':
                    value = self._pulse_mod_func_to_internal
                # elif key == 'frequency':
                #     if value > RANGE_MAX or value < RANGE_MIN:
                #         raise ValueError("Invalid frequency. All frequencies must be between 2.025 GHz and 4.050 GHz.")
                key = self._param_to_internal(key)

                # only send update to instrument if connection to instrument has been established
                if self._settings_initialized:
                    self.srs.write(key + ' ' + str(value)) # frequency change operation timed using timeit.timeit and
                                                           # completion confirmed by query('*OPC?'), found delay of <10ms
                    # ER 20180904
                   # if key == 'FREQ':
                   #     print('frequency set to: ', float(self.srs.query('FREQ?')))
                    # print(self.srs.query('*OPC?'))

        # XXXXX MW ISSUE = END
        # ===========================================

    @property
    def _PROBES(self):
        return{
            'enable_output': 'if type-N output is enabled',
            'frequency': 'frequency of output in Hz',
            'amplitude': 'type-N amplitude in dBm',
            'phase': 'phase',
            'enable_modulation': 'is modulation enabled',
            'modulation_type': 'Modulation Type: 0= AM, 1=FM, 2= PhaseM, 3= Freq sweep, 4= Pulse, 5 = Blank, 6=IQ',
            'modulation_function': 'Modulation Function: 0=Sine, 1=Ramp, 2=Triangle, 3=Square, 4=Noise, 5=External',
            'pulse_modulation_function': 'Pulse Modulation Function: 3=Square, 4=Noise(PRBS), 5=External',
            'dev_width': 'Width of deviation from center frequency in FM'
        }

    def read_probes(self, key):
        # assert hasattr(self, 'srs') #will cause read_probes to fail if connection not yet established, such as when called in init
        assert(self._settings_initialized) #will cause read_probes to fail if settings (and thus also connection) not yet initialized
        assert key in list(self._PROBES.keys())

        #query always returns string, need to cast to proper return type
        if key in ['enable_output', 'enable_modulation']:
            key_internal = self._param_to_internal(key)
            value = int(self.srs.query(key_internal + '?'))
            if value == 1:
                value = True
            elif value == 0:
                value = False
        elif key in ['modulation_type','modulation_function','pulse_modulation_function']:
            key_internal = self._param_to_internal(key)
            value = int(self.srs.query(key_internal + '?'))
            if key == 'modulation_type':
                value = self._internal_to_mod_type(value)
            elif key == 'modulation_function':
                value = self._internal_to_mod_func(value)
            elif key == 'pulse_modulation_function':
                value = self._internal_to_pulse_mod_func(value)
        else:
            key_internal = self._param_to_internal(key)
            value = float(self.srs.query(key_internal + '?'))

        return value

    @property
    def is_connected(self):
        try:
            self.srs.query('*IDN?') # arbitrary call to check connection, throws exception on failure to get response
            return True
        except pyvisa.errors.VisaIOError:
            return False

    def _param_to_internal(self, param):
        """
        Converts settings parameters to the corresponding key used for GPIB commands in the SRS.
        Args:
            param: settings parameter, ex. enable_output

        Returns: GPIB command, ex. ENBR

        """
        if param == 'enable_output':
            return 'ENBR'
        elif param == 'frequency':
            return 'FREQ'
        elif param == 'amplitude':
            return 'AMPR'
        elif param == 'phase':
            return 'PHAS'
        elif param == 'enable_modulation':
            return 'MODL'
        elif param == 'modulation_type':
            return 'TYPE'
        elif param == 'modulation_function':
            return 'MFNC'
        elif param == 'pulse_modulation_function':
            return 'PFNC'
        elif param == 'dev_width':
            return 'FDEV'
        else:
            raise KeyError

    def _mod_type_to_internal(self, value):
        #COMMENT_ME
        if value == 'AM':
            return 0
        elif value == 'FM':
            return 1
        elif value == 'PhaseM':
            return 2
        elif value == 'Freq sweep':
            return 3
        elif value == 'Pulse':
            return 4
        elif value == 'Blank':
            return 5
        elif value == 'IQ':
            return 6
        else:
            raise KeyError

    def _internal_to_mod_type(self, value):
        #COMMENT_ME
        if value == 0:
            return 'AM'
        elif value == 1:
            return 'FM'
        elif value == 2:
            return 'PhaseM'
        elif value == 3:
            return 'Freq sweep'
        elif value == 4:
            return 'Pulse'
        elif value == 5:
            return 'Blank'
        elif value == 6:
            return 'IQ'
        else:
            raise KeyError

    def _mod_func_to_internal(self, value):
        #COMMENT_ME
        if value == 'Sine':
            return 0
        elif value == 'Ramp':
            return 1
        elif value == 'Triangle':
            return 2
        elif value == 'Square':
            return 3
        elif value == 'Noise':
            return 4
        elif value == 'External':
            return 5
        else:
            raise KeyError

    def _internal_to_mod_func(self, value):
        #COMMENT_ME
        if value == 0:
            return 'Sine'
        elif value == 1:
            return 'Ramp'
        elif value == 2:
            return 'Triangle'
        elif value == 3:
            return 'Square'
        elif value == 4:
            return 'Noise'
        elif value == 5:
            return 'External'
        else:
            raise KeyError

    def _pulse_mod_func_to_internal(self, value):
        #COMMENT_ME
        if value == 'Square':
            return 3
        elif value == 'Noise(PRBS)':
            return 4
        elif value == 'External':
            return 5
        else:
            raise KeyError

    def _internal_to_pulse_mod_func(self, value):
        #COMMENT_ME
        if value == 3:
            return 'Square'
        elif value == 4:
            return 'Noise(PRBS)'
        elif value == 5:
            return 'External'
        else:
            raise KeyError

if __name__ == '__main__':
    # from pylabcontrol.core import Instrument
    #
    # instruments =        {"MicrowaveGenerator": {
    #         "class": "MicrowaveGenerator",
    #         "settings": {
    #             "enable_output": False,
    #             "enable_modulation": True,
    #             "amplitude": -60,
    #             "GPIB_num": 0,
    #             "frequency": 3000000000.0,
    #             "dev_width": 32000000.0,
    #             "pulse_modulation_function": "External",
    #             "phase": 0,
    #             "modulation_function": "External",
    #             "port": 27,
    #             "modulation_type": "FM"
    #         }
    #     }}
    #
    # instr = {}
    #
    # instr, loaded_failed = Instrument.load_and_append(
    #     {name: instruments[name] for name in instruments.keys()}, instr)
    #
    # print(instr)
    # print(loaded_failed)
    # import pylabcontrol.instruments.microwave_generator MicrowaveGenerator

    mw = MicrowaveGenerator()
    # mw = MicrowadveGenerator(settings={'enable_modulation': True, 'frequency': 3000000000.0, 'dev_width': 32000000.0, 'pulse_modulation_function': 'External', 'phase': 0, 'port': 27, 'modulation_type': 'FM', 'enable_output': False, 'GPIB_num': 0, 'amplitude': -60, 'modulation_function': 'External'})

    print((mw.srs))

    # instrument_name= 'MicrowaveGenerator'
    # instrument_settings = {'enable_modulation': True, 'frequency': 3000000000.0, 'dev_width': 32000000.0, 'pulse_modulation_function': 'External', 'phase': 0, 'port': 27, 'modulation_type': 'FM', 'enable_output': False, 'GPIB_num': 0, 'amplitude': -60, 'modulation_function': 'External'}
    # class_of_instrument = MicrowaveGenerator
    # instrument_instance = class_of_instrument(name=instrument_name, settings=instrument_settings)
