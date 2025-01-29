import time

import pygame
from PyQt5.QtCore import QThread, pyqtSignal,QObject
import json
import os
from pathlib import Path

# This is a simple class that will help us print to the screen.
# It has nothing to do with the joysticks, just outputting the
# information.
class PS4_controller(QThread):
    button = pyqtSignal(str,float)

    def __init__(self):
        super(PS4_controller,self).__init__()

        pygame.init()

        self.joysticks = []
        self.terminate = False

        cwd = Path.cwd()

        with open(os.path.join(cwd,'utils','device','ps4_keys.json'), 'r+') as file:
            self.button_keys = json.load(file)

        self.analog_keys = {0:0, 1:0, 2:0, 3:0, 4:-1, 5: -1 }

        self.joysticks.append(pygame.joystick.Joystick(0))
        for joystick in self.joysticks:
            try:
                joystick.init()
                print(joystick)
            except Exception as e:
                print("Error in init joystic", e)

        print("__init__ DONE")



    def run(self) -> None:

        while not self.terminate:
            button = None
            value = 0

            ################################# CHECK PLAYER INPUT #################################
            try:
                event = pygame.event.poll()
            except:
                self.terminate = True
            if event.type == pygame.QUIT:
                self.terminate = True
            if event.type == pygame.KEYDOWN:
                ############### UPDATE SPRITE IF SPACE IS PRESSED #################################
                pass

            if event.type != 0:
                button, value = "event: {} \n Event Type: {} \n".format(event, event.type), 0

                self.button.emit(button, 10 * float(value))

            if event.type == 0: continue

            # HANDLES BUTTON PRESSES
            elif event.type == 1539:
                if event.button == self.button_keys['left_arrow']:
                    button,value = 'left',1
                elif event.button == self.button_keys['right_arrow']:
                    button,value = 'right',1
                elif event.button == self.button_keys['down_arrow']:
                    button,value = 'down',1
                elif event.button == self.button_keys['up_arrow']:
                    button,value = 'up',1
                elif event.button == self.button_keys['L1']:
                    button,value = 'L1',0.1
                elif event.button == self.button_keys['R1']:
                    button,value = 'R1',0.1
                else:
                    button,value = event.type,0
            # HANDLES BUTTON RELEASES
            # if event.type == pygame.JOYBUTTONUP:
            #     if event.button == self.button_keys['left_arrow']:
            #         self.button.emit('LEFT = False')
            #     if event.button == self.button_keys['right_arrow']:
            #         self.button.emit('RIGHT = False')
            #     if event.button == self.button_keys['down_arrow']:
            #         self.button.emit('DOWN = False')
            #     if event.button == self.button_keys['up_arrow']:
            #         self.button.emit('UP = False')

            # HANDLES ANALOG INPUTS
            elif event.type == 1536:
                value = event.value
                # print(analog_keys)
                # Horizontal Analog - left joyaxis
                if event.axis == 0:
                    if value < 0:
                        button,value = 'left',abs(value/2)
                    elif value > 0:
                        button,value = 'right',value/2
                    if value>0.4:value=0.6
                # Vertical Analog
                elif event.axis == 1:
                    if value < 0:
                        button,value = 'up',abs(value/2)
                    elif value > 0:
                        button,value = 'down',value/2
                    if value>0.1:value=0.1
                #Vertical analgue - right
                elif event.axis == 2:
                    if value < 0:
                        button,value = 'z-up',.2
                    elif value > .0:
                        button,value = 'z-down',.2
                # Horizontal Analog - right
                elif event.axis == 3:
                    if value < 0:
                        continue
                    elif value > 0:
                        continue
                    # Triggers
                # Triggers
                elif event.axis == 4:  # Left trigger
                    if value > .2:
                        button,value = 'LT',.1
                    else: continue
                elif event.axis == 5:  # Left trigger
                    if value > .2:
                        button,value = 'RT',.1

            else:
                button, value = "event: {} \n Event Type: {} \n".format(event, event.type), 0
            try:
                self.button.emit(button, 10 * float(value))
            except:
                pass
            time.sleep(.1)

            if len(pygame.event.get()) > 250:
                pygame.event.clear()


        if self.terminate:
            print("End controller thread.")