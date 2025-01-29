The michelson_interface.ui file can be edited by the QTDesigner. Upon changes in there, the michelson_interface.py
file needs to be created anew based on the updated ui file. This can be done by executing the command

    python -m PyQt5.uic.pyuic -x michelson_interface.ui -o michelson_interface.py

from a console.