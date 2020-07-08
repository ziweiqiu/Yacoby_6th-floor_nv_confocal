import time
import numpy as np

from b26_toolkit.instruments import NI6353
from pylabcontrol.core import Parameter, Script



class FlipMountControls(Script):
    """
    This script flips on and off the mirrors by sending NI digital output pulses.
    written by ZQ 2/12/2019 6:30pm
    """
    _DEFAULT_SETTINGS = [
        Parameter('white_light_camera_flip', False, bool, 'Flip the white light/camera pellicle'),
        Parameter('APD_lightblock_flip', False, bool, 'Flip the APD light block'),
        Parameter('power_meter_flip', False, bool, 'Flip the power meter'),
        Parameter('G2_flip', False, bool, 'Flip the G2 mirror'),
        Parameter('do_channels',
                  [
                      Parameter('white_light_camera_flip', 'do1', ['do0', 'do1', 'do3', 'do4'],
                                'channel controlling white_light_camera_flip'),
                      Parameter('APD_lightblock_flip', 'do4', ['do0', 'do1', 'do3', 'do4'],
                                'channel controlling white_light_camera_flip'),
                      Parameter('power_meter_flip', 'do3', ['do0', 'do1', 'do3', 'do4'],
                                'channel controlling white_light_camera_flip'),
                      Parameter('G2_flip', 'do0', ['do0', 'do1', 'do3', 'do4'],
                                'channel controlling white_light_camera_flip'),
                  ])
    ]

    _INSTRUMENTS = {'daq': NI6353}

    _SCRIPTS = {}

    def __init__(self, instruments, scripts=None, name=None, settings=None, log_function=None, data_path=None):
        """
        Example of a script that emits a QT signal for the gui
        Args:
            name (optional): name of script, if empty same as class name
            settings (optional): settings for this script, if empty same as default settings
        """
        Script.__init__(self, name, settings=settings, scripts=scripts, instruments=instruments,
                        log_function=log_function, data_path=data_path)


    def _function(self):
        """
        This is the actual function that will be executed. It uses only information that is provided in the settings property
        will be overwritten in the __init__
        """
        self._setup_daq()

        for key, value in self.settings.items():
            if key is not 'path' and key is not 'tag' and key is not 'save' and key is not 'do_channels' and value is True:
                print(key)
                self.daq_out.set_digital_output(output_dict={self.settings['do_channels'][key]: True})
                time.sleep(0.1)
                self.daq_out.set_digital_output(output_dict={self.settings['do_channels'][key]: False})
                time.sleep(0.1)
                # ZQ 2/13/2019: sending pulses on different lines simultaneously doesn't work... need to do it sequentially



        # self.daq_out.run(task)
        # self.daq_out.waitToFinish(task)
        # self.daq_out.stop(task)
        self.log('flip mounts are updated')


    def _setup_daq(self):
        # defines which daqs contain the input and output based on user selection of daq interface

        self.daq_out = self.instruments['daq']['instance']


if __name__ == '__main__':
    script = {}
    instr = {}
    # print('line 511')
    FlipMountControls
    # print('line 513')
    # script, failed, instr = Script.load_and_append({'Daq_Read_Cntr': 'Daq_Read_Cntr'}, script, instr)
    script, failed, instr = Script.load_and_append({'FlipMountControls': 'FlipMountControls'}, script, instr)

    print(script)
    print(failed)
    print(instr)
