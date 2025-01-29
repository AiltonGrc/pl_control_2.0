import os
from datetime import datetime
from pathlib import Path
from PyQt5.QtCore import QObject, pyqtSignal, QUrl


class ErrorLogger(QObject):
    logging_updated = pyqtSignal()
    error_signal = pyqtSignal(tuple)
    device_busy_signal = pyqtSignal(bool)

    def __init__(self, filename):
        super().__init__()

        # create logging file (.txt)
        marple_path = Path(__file__).parents[2]
        self.logger_path = marple_path / 'logfiles'

        os.makedirs(self.logger_path, exist_ok=True)

        file_suffix = datetime.now().strftime("%Y%m%d_%H-%M-%S")

        self.filename = filename + file_suffix + ".txt"

        self.file = self.logger_path / self.filename

        current_path = os.getcwd()
        os.chdir(self.logger_path)
        with open(self.filename, 'w') as f:
            f.write("Error log: \n")
        os.chdir(current_path)
        self.source_path = QUrl.fromLocalFile(str(self.file))

    def log_error_msg(self, error_msg: tuple):
        print("Kenny logds in")
        device, state, info = error_msg
        if state:
            message_type = "Warning"
        else:
            message_type = "Error"

        log_message = f"{message_type}\t-\t{device}: {info}"

        current_path = os.getcwd()
        os.chdir(self.logger_path)
        with open(self.filename, 'a') as f:
            f.write(log_message)
            f.write("\n\r")

        os.chdir(current_path)
        self.logging_updated.emit()


if __name__ == "__main__":
    logger = ErrorLogger("test_log")
    logger.log_error_msg(("Test", False, "Error: XY-Stage failed"))
    print(logger.file)