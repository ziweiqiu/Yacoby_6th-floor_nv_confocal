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

# from test_script import ScriptTest
from .galvo_scan.galvo_scan import GalvoScan
from .galvo_scan.confocal_scan_LISE607RT import ConfocalScan
from .set_laser import SetLaser, SetConfocal
from .daq_read_counter import Daq_Read_Counter
from .flip_mount_controls import FlipMountControls
from .magnet_sweep import MagnetSweep, MagnetSweep2D, MagnetSweep1DAdv, MagnetCoilSweep, ElectricScan1D, ElectricScan2D
# from .optimize import optimizeZ
# from .take_image_camera import TakeImage
# from .esr import ESR
from .esr_RnS import ESR_RnS, ESR_FastSwp_RnS, ESR_FastSwp_RnS_FitGuaranteed, CalMagnetAngle
from .esr_N9310A import ESR_N9310A
from .esr_rabi_echo import ESR_Rabi_RnS, ESR_Rabi_Echo_RnS
# from .esr_rabi_echo_rns import ESRandRabi_RnS
# from .esr_dithering import ESR_FM_Dither
# from .esr_two_freq_continuous import ESRTwoFreqContinuous
# from .spec_analyzer_get_spectrum import SpecAnalyzerGetSpectrum
# from .zi_sweeper import ZISweeper
# from .zi_high_res_sweep import ZISweeperHighResolution
from .find_nv import FindNV
from .read_power import IntensityWheel_Calibration, Saturation
# from .atto_scan import AttoStep
from .pulse_sequences.rabi import Rabi_RnS, Rabi_N9310A, Rabi_RnSIQ, Ramsey_RnS, Ramsey_N9310A, Ramsey_CoilVol, Ramsey_RnSIQ, TwoPhotonRabi
from .pulse_sequences.hahn_echo import PDD_N9310A, eSensing_N9310A, eSensing_swpV, PDD_RnSIQ, PDD_XYreadout, angle_sensing
from .pulse_sequences.pulsed_esr import PulsedESR
from .pulse_sequences.t1 import T1SingleInit
from .pulse_sequences.calibration import GrDelayMeas, IQCalibration_N9310A
# # from .pulse_sequences import XY8_k, T1, Rabi, PDD, XY4, T1SingleInit, PulsedESR, \
# #     HahnEcho, XY4, XYXY, ReadoutStartTimeWithoutMW, ReadoutStartTime, ReadoutDuration, CPMG, \
# #     HahnEchoManyNVs, RabiPowerSweepSingleTau
# from .pulse_sequences.rabi import Rabi
# from .esr_and_rabi import ESRAndRabi
# # from .spec_analyzer_get_spectrum import KeysightGetSpectrum
# from .light_control import ApplyLightControlSettings, CameraOn
from .correlate_images import Track_Correlate_Images, Take_And_Correlate_Images
# from .autofocus import AutoFocusDAQ
# from .autofocus import AutoFocusDAQ, AutoFocusTwoPoints, AutoFocusTwoPointsFR, AutoFocusDaqSMC, AutoFocusCameraSMC
# from .record_pressures import RecordPressures
# from .set_magnetic_coils import SetMagneticCoils
# from .align_magnetic_field_to_NV import AlignFieldToNV
# from .Ni_9263_polarization_controller import Ni9263_BalancePolarization
# from .stability_with_microwaves import Stability_With_Microwaves
# from .read_temperature_lakeshore import ReadTemperatureLakeshore
