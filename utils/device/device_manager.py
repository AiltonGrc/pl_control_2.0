from abc import ABC, abstractmethod
import typing


class DeviceManager(ABC):

    def __init__(self):
        self.linear = None
        self.piezo = None
        self.trigger = None

    @abstractmethod
    def connect(self) -> typing.Tuple[bool, bool, bool]:
        pass

    @abstractmethod
    def disconnect(self) -> None:
        pass

    @abstractmethod
    def get_linear_pos(self) -> float:
        pass

    @abstractmethod
    def get_piezo_pos(self) -> float:
        pass

    @abstractmethod
    def switch_servo(self, state: bool):
        pass

    @abstractmethod
    def set_linear_velocity(self, value: float) -> None:
        pass

    @abstractmethod
    def set_linear_acceleration(self, value: float) -> None:
        pass

    @abstractmethod
    def set_linear_deceleration(self, value: float) -> None:
        pass

    @abstractmethod
    def goto_linear(self, position: float) -> None:
        pass

    @abstractmethod
    def goto_piezo(self, position: float) -> None:
        pass

    @abstractmethod
    def trig(self, low_handshake=True) -> None:
        pass
