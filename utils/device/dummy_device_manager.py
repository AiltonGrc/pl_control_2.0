from device.device_manager import DeviceManager
import time
import typing
import numpy as np


class DummyDeviceManager(DeviceManager):

    def connect(self) -> typing.Tuple[bool, bool, bool]:
        print("connect to dummy device manager")
        time.sleep(0.2)
        return True, True, True

    def disconnect(self) -> None:
        print("disconnect from dummy device manager")
        return

    def get_linear_pos(self) -> float:
        return np.random.rand() * 300

    def get_piezo_pos(self) -> float:
        return np.random.rand() * 50

    def switch_servo(self, state: bool):
        pass

    def set_linear_velocity(self, value: float) -> None:
        time.sleep(0.2)

    def set_linear_acceleration(self, value: float) -> None:
        time.sleep(0.2)

    def set_linear_deceleration(self, value: float) -> None:
        time.sleep(0.2)

    def goto_linear(self, position: float) -> None:
        time.sleep(0.2)

    def goto_piezo(self, position: float) -> None:
        time.sleep(0.2)

    def trig(self, low_handshake=True) -> None:
        time.sleep(0.2)
