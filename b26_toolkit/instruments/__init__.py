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

# from .gauge_controller import PressureGauge, PumpLinePressureGauge, ChamberPressureGauge
# from .spectrum_analyzer import SpectrumAnalyzer
# from .ni_daq import NI6259, NI9263, NI9402, NI9219
from .ni_daq import NI6353
# from .ni_daq_v2 import NIPCIe6353
# from .piezo_controller import PiezoController
# from .zurich_instruments import ZIHF2
from .pulse_blaster import LISE607RTPulseBlaster, Pulse
# from .maestro import MaestroLightControl
# from .attocube import Attocube
from .microwave_generator import R8SMicrowaveGenerator, AgilentMicrowaveGenerator, AgilentMicrowaveGeneratorII

# from .microwave_generator import MicrowaveGenerator
# from .magnet_coils import MagnetCoils
# from .temperature_controller import TemperatureController
from .whiteLED_controller import WhiteLEDController
from .power_meter import PM100D

# from .montana import CryoStation
# from .newport_smc100 import SMC100
from .awg import Agilent33120A
# from .keysight_oscilloscope import Oscilloscope
# from .thorlabs_kcube import KDC001, TLI_DeviceInfo, B26KDC001x, B26KDC001y, B26KDC001z
# from .thorlabs_kcube import TDC001, TLI_DeviceInfo, L607RT_TDC001_sampleX, L607RT_TDC001_sampleY, L607RT_TDC001_sampleZ, L607RT_TDC001_magnetX, L607RT_TDC001_magnetY, L607RT_TDC001_magnetZ
from .dcservo_kinesis_dll import TDC001_II, KDC101_II, MagnetX, MagnetY, MagnetZ, IntensityWheel
# from .magnet_coils import MagnetCoils
# from .ueye_camera import UEyeCamera
# from .optotune_lens import OptotuneLens
# from .R8S_microwave_generator import R8SMicrowaveGenerator