from device.device_manager import DeviceManager
from pipython import GCSDevice # Physikalische Instrumente Python Interface
import time


class PIDeviceManager(DeviceManager):

    def __init__(self):
        super(PIDeviceManager, self).__init__()

        # ID's for moving stages
        self.linear_id = 1
        self.piezo_id = 'A'
        # DIO's for trigger
        self.in_id = 4
        self.out_id = 1

    def connect(self):
        try:
            self.disconnect()
        except:
            print("cannot close connection because they are not there yet")
            pass

        linear_connected = True
        try:
            self.linear = GCSDevice()
            self.linear.ConnectUSB(serialnum='0145500570')
            print(self.linear.qIDN().strip())
        except:
            linear_connected = False

        piezo_connected = True
        try:
            self.piezo = GCSDevice('E-816')
            self.piezo.ConnectUSB(serialnum='20296')
            print(self.piezo.qIDN().strip())
        except:
            piezo_connected = False

        return linear_connected, piezo_connected, linear_connected

    def disconnect(self) -> None:
        for device in [self.linear, self.piezo, self.trigger]:
            if device is not None:
                device.CloseConnection()

    def get_linear_pos(self) -> float:
        return self.linear.qPOS(f'{self.linear_id}')[f'{self.linear_id}']

    def get_piezo_pos(self) -> float:
        return self.piezo.qPOS(self.piezo_id)[self.piezo_id]

    def set_linear_velocity(self, value: float) -> None:
        self.linear.VEL(f'{self.linear_id}', value)

    def set_linear_acceleration(self, value: float) -> None:
        self.linear.ACC(f'{self.linear_id}', value)

    def set_linear_deceleration(self, value: float) -> None:
        self.linear.DEC(f'{self.linear_id}', value)

    def switch_servo(self, state: bool):
        self.linear.SVO(f'{self.linear_id}', state)

    def goto_linear(self, position: float) -> None:
        # switch servo on to be able to move
        self.switch_servo(True)
        self.linear.MOV(f'{self.linear_id}', position)
        # wait until on target
        on_target = 0
        while on_target < 10:
            if self.linear.qONT()[f'{self.linear_id}']:
                on_target += 1
            else:
                on_target = 0
            time.sleep(0.05)
        # switch servo off to prevent additional motion
        self.switch_servo(False)


    def goto_piezo(self, position: float) -> None:
        self.piezo.MOV(self.piezo_id, position)
        on_target = 0
        while on_target < 3:
            if self.piezo.qONT(self.piezo_id)[self.piezo_id]:
                on_target += 1
            else:
                on_target = 0
            time.sleep(0.05)

    def wait_input(self, pi_device: GCSDevice, in_id: int, state: bool, max_waittime: float = 10.) -> bool:
        start = time.time()

        while True:
            if pi_device.qDIO(f"{in_id}")[f"{in_id}"] == state:
                return True
            elif time.time() - start > max_waittime:
                return False

    def trig(self, low_handshake=True):
        # for now only print internal message if handshake fails
        self.linear.DIO(f"{self.out_id}", True)  # rising edge triggering spectrometer rec.
        if self.wait_input(self.linear, self.in_id, True) is not True:
            print("Error in HIGH handshake")

        if low_handshake:
            self.linear.DIO(f"{self.out_id}", False)  # falling edge resets output channel for next rec.
            if self.wait_input(self.linear, self.in_id, False) is not True:
                print("Error in LOW handshake")