from PyQt5.QtGui import QDoubleValidator, QValidator
import typing

class QDoubleRangeValidator(QDoubleValidator):

    def __init__(self, bottom, top, decimals):
        super(QDoubleRangeValidator, self).__init__(bottom, top, decimals, notation=QDoubleValidator.StandardNotation)

    def validate(self, a0: str, a1: int) -> typing.Tuple[QValidator.State, str, int]:
        state, a0, a1 = super(QDoubleRangeValidator, self).validate(a0, a1)
        if state == QValidator.State.Intermediate and a0 != "":
            return QValidator.State.Invalid, a0, a1
        return state, a0, a1

