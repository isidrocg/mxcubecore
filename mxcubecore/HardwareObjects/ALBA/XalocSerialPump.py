#  Project: MXCuBE
#  https://github.com/mxcube.
#
#  This file is part of MXCuBE software.
#
#  MXCuBE is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  MXCuBE is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with MXCuBE.  If not, see <http://www.gnu.org/licenses/>.

"""
[Name] ALBATestPump
[Description]
HwObj used to control Shimadzu pump for SSX
[Signals]
- valueChanged
"""

from __future__ import print_function

import logging
import math

# from mxcubecore.BaseHardwareObjects import Device
from mxcubecore.BaseHardwareObjects import HardwareObject

__credits__ = ["ALBA Synchrotron"]
__version__ = "2.3"
__category__ = "Test"


class XalocSerialPump(HardwareObject):

    def __init__(self, *args):
        self.logger = logging.getLogger("HWR.SerialPump")
        HardwareObject.__init__(self, *args)
        
        # Initialize properties
        self.max_flow = None
        self.min_flow = None
        self.max_pressure = None
        self.set_point = None
        self.flow = None
        self.sample_vol = None
        
        self.capillary_id = None
        self.reservoir = None

        # Initialize channels
        self.channelflow = None
        self.channelpressure = None
        self.channelpumping = None
        self.channelunit = None
        self.phase_channel = None

        # Initialize commands
        self._cmdStartPumping = None
        self._cmdStopPumping = None


    def init(self):
        self.logger.debug("Initializing {0}".format(self.__class__.__name__))

        # Properties for PID control. Instead of properties in xml, read them
        # from the ui for live control
        # self.max_flow = self.get_property("max_flow")
        self.min_flow = self.get_property("min_flow")
        self.max_pressure = self.get_property("max_press")
        # self.set_point = self.get_property("setpoint")
        


        # Properties for sample consumption and jet speed calculation
        # self.capillary_id = self.get_property("capillary_id")
        # self.reservoir = self.get_property("reservoir")
        # if self.reservoir == 40:
        #     self.amp_factor = 14
        # elif self.reservoir == 20:
        #     self.amp_factor = 34

        

        ## Initialize pumped_volume and previous time
        #self.pumped_volume = 0
        #self.previous_time = time.time()
        #self.jet_speed = 0


        # Channels
        self.channelflow = self.get_channel_object("PumpFlow")
        if self.channelflow is not None:
            self.connect(self.channelflow, "update", self.flow_changed)
            # self.channelflow.connectSignal('update', self.flow_changed)

        self.channelunit = self.get_channel_object("PumpPressureUnit")
        if self.channelunit is not None:
            self.connect(self.channelunit, "update", self.unit_changed)
            # self.pressure_unit = self.channelunit.get_value()
            
        self.channelpressure = self.get_channel_object("PumpPressure")
        if self.channelpressure is not None:
            self.connect(self.channelpressure, "update", self.pressure_changed)
            # self.channelpressure.connectSignal('update', self.pressure_changed)
            self.logger.debug("Pressure: {0}".format(self.channelpressure.get_value()))

        self.channelpumping = self.get_channel_object("PumpPumping")
        if self.channelpumping is not None:
            self.connect(self.channelpumping, "update", self.pumping_changed)
            # self.channelpumping.connectSignal('update', self.pumping_changed)
            self.logger.debug("Pumping: {0}".format(self.channelpumping.get_value()))
                
        # Commands
        self._cmdStartPumping = self.get_command_object("_cmdStartPumping")
        self._cmdStopPumping = self.get_command_object("_cmdStopPumping")
        
        # Command to go to supervisor transfer and serial. Implement serial 
        # phase in supervisor. This is for button to prepare for serial or 
        # restore to monocrystal. Implement warning to accept, long time to 
        # change from/to serial

        #self.go_transfer_cmd = self.getCommandObject("go_transfer")
        #self.go_serial_cmd = self.getCommandObject("go_serial")
        
    def getFlow(self):
        """
        Get the pump flow in microliters per minute.

        Returns: float

        """
        return self.channelflow.get_value() * 1000

    def get_jet_speed(self, pump_flow, reservoir, capillary_id):
        """
        Given the pump_flow in ul/min, the reservoir used (40 or 20 ul) and the
        capillary id (in um), get the jet speed in micrometers per second.

        Returns: float

        """
        # Pump amplification factor:
        if reservoir == 40:
            amp_factor = 14
        elif reservoir == 20:
            amp_factor = 34

        # Convert pump speed (ul/min) to extruder speed (ul/s) 
        extruder_speed = pump_flow / (amp_factor * 60)

        # Capillary area in um2
        capillary_area = math.pi * (capillary_id / 2) ** 2

        # Jet speed:
        # extruder_speed (ES) is volume per unit of time,
        # Volume is Area * Height, jet speed will be the height per unit of time
        # so ES = V/s, V = A*h -> h = V/A. Jet speed (JS) is then 
        # JS = h/s = ES/A. Apply factor 1E9 to get um/s
        jet_speed = extruder_speed * 1E9 / capillary_area

        return round(jet_speed, 2)

    
    #def getPumping(self):
        #"""
        #Get the pumping state of the serial pump.

        #Returns: str

        #"""
        #pumping = self.channelpumping.getValue()
        #if pumping:
            #return "Pumping"
        #else:
            #return "Paused"
        
    def getPressure(self, pressure):
        '''Calculate pressure in bar in all the pump units available:
        Deduced pressure units: {0: 'bar', 1: 'psi', 2: 'Mpa', 3: 'kgf/cm2'}
        Unit conversions: 10 bar = 1 MPa / 145 psi = 1 MPa / 1 MPa = 1 MPa / 10,197 kgf/cm2 = 1 MPa
        
        Returns pump pressure in bar
        '''
        pump_units = {0: ['bar', 1.0], 1: ['psi', 14.5], 2: ['Mpa', 0.1], 3: ['kgf/cm2', 1.0197]}
        
        # Best setting is to retrieve pump pressure in psi, as the device server
        # gives only integers. Then convert to bar
        pressure_unit = self.channelunit.get_value()
        # pressure_unit = self.pressure_unit

        # pressure = self.channelpressure.getValue()
        
        pressure_bar = pressure / pump_units[pressure_unit][1]
        pressure_bar = round(pressure_bar, 1)
        
        return pressure_bar
        
    def flow_changed(self, value):
        # self.flow = value --> Signal emits flow in ml/min, maybe better to use value instead getting flow again. 
        # self.emit_flow_changed()
        
        # Convert the received flow value in ml/min to ul/min
        flow = value * 1000
        
        if flow != self.flow:
            self.flow = flow
            self.emit('flowChanged', flow)
            # self.emit('test', flow)
            self.logger.debug("flow_changed: flowChanged emitted ({})".format(flow))
            # self.logger.debug("flow_changed: value ({})".format(value))
        # else:
            # self.logger.debug("Flow not changed")
        
    # def emit_flow_changed(self):
        # self.logger.debug(" emitting pump info")
        # if self.flow != 1000:
            # self.emit("flowChanged", (self.flow, ))
    
    
    # def flowChanged(self):
    #     """
    #     Slot receiving the zoom position changed and emitting the signal.
    #     Additionally, we could also actuate the backlight intensity here.

    #     Args:
    #         position: Zoom int position

    #     Returns: None

    #     """
    #     # TODO: Arrives an int, emits string
    #     # We should use the position received.
    #     #self.logger.debug("positionChange event received: {}".format(position))
    #     #position = self.getCurrentPositionName()
    #     self.current_flow = self.getFlow 

    #     if current_flow != self.flow:
    #         self.logger.debug(
    #             "Flow changed: {} -> {}".format(self.flow, current_flow)
    #         )
    #         self.flow = current_flow
    #         self.emit("flowChanged", self.flow)
    #         self.logger.debug("Emitting Flow: {0}".format(self.current_flow))
    #         #self.emit("predefinedPositionChanged", (position, 0))
    #         # TODO: blight range is 0-30
    #         # Links the zoom position with the backlight value
    #         # self.chan_blight.setValue(int(self.current_position.split()[0]))
            
    def pressure_changed(self, value):
        
        #Convert input value to bar, pump unit may be changed in DS
        pressure_bar = self.getPressure(value)
        
        #self.logger.debug("pressure_changed: value ({})".format(value))
        self.emit('pressureChanged', pressure_bar) # Check here the units and emit always in bar
        # self.logger.debug("pressure_changed: pressureChanged emitted ({})".format(pressure_bar))

        # Stop pump before hardware alarm/stop, prevent manual intervention
        if pressure_bar > 28:
            self._cmdStopPumping()
            self.logger.debug("High pressure ({}), stopping pump".format(pressure_bar))
        
            
    def pumping_changed(self, value):
        self.logger.debug("pumping_changed: value ({})".format(value))
        self.emit('pumpingChanged', value)
        self.logger.debug("pumping_changed: pumpingChanged emitted ({})".format(value))

    def unit_changed(self, value):
        self.logger.debug("unit_changed: value ({})".format(value))
        self.emit('unitChanged', value)
        self.logger.debug("unit_changed: unitChanged emitted ({})".format(value))
        self.pressure_unit = value

    def update_values(self):
        # We re-emit the channel value received in the event to any Qt client
        flow = self.getFlow()
        self.emit('flowChanged', flow)
        
        pressure = self.channelpressure.get_value()
        pressure_bar = self.getPressure(pressure)
        self.emit('pressureChanged', pressure_bar)

        pumping = self.channelpumping.get_value()
        self.emit('pumpingChanged', pumping)
        

    def flow_control(self, pressure, target_pressure, max_flow, kp):
        '''
        '''
        self.logger.debug("Flow control activated")
        # Kp = 0.01
        Kp = kp
        MV_bar = self.getFlow() # in ul/min
        # SP = self.set_point
        SP = target_pressure
        PV = pressure
        
        e = SP - PV
        P = Kp*e
        MV = (MV_bar + P)

        self.logger.debug("MV_bar = %s" % str(MV_bar))
        self.logger.debug("SP = %s" % str(SP))
        self.logger.debug("PV = %s" % str(PV))
        self.logger.debug("Kp = %s" % str(Kp))
        self.logger.debug("P = %s" % str(P))
        self.logger.debug("MV = %s" % str(MV))
        
        # if MV > self.max_flow:
        #     MV = self.max_flow
        if MV > max_flow:
            MV = max_flow
        elif MV < self.min_flow:
            MV = self.min_flow
        
        self.logger.debug("Flow control: new_flow = %s" % str(MV))
        #self.channelflow.setValue(0.001)
        # self.channelflow.set_value(MV)
        self.set_flow(MV)

    def set_flow(self, ul_min):

        self.channelflow.set_value(ul_min / 1000)

    def consumed_sample(self):
        '''
        '''
        if self.sample_vol:
            flow = self.getFlow()
            previous_time = self.time
            time_now = time.perf_counter()
            delta_minutes = (time_now - previous_time) / 60

            volume = flow * delta_minutes

            self.time = time_now
            self.sample_vol -= volume

        

    #def _pid(self, Kp, Ki, Kd, MV_bar=0):
        #Connect this function to pressure signal, then update the flow according
        #to the read pressure value
        
        #"""Creates PID controllers with specified gain for proportional (Kp),
        #integral (Ki) and derivative (Kd) components"""
        ## initialize stored data
        #e_prev = 0
        #t_prev = time.time() - 100
        #I = 0
        
        ## initial control
        #MV = MV_bar
        
        #while True:
            ## yield MV, wait for new t, PV, SP
            #t, PV, SP = yield MV
            
            ## PID calculations
            #e = SP - PV
            
            #P = Kp*e
            #I = I + Ki*e*(t - t_prev)
            #D = Kd*(e - e_prev)/(t - t_prev)
            
            ## For PID tunning:
            #MV = MV_bar + P

            ## Uncomment when PID is tuned
            ## MV = MV_bar + P + I + D

            ## Control flow for max/min desired pump flow (Could be DS attributes)
            ## Check this, we are changing MV and maintaining error, 
            ## this affects the next cycle calculation.
            #if MV > self.max_flow:
                #MV = self.max_flow
                #I = 0 # reset integral solves?
            #elif MV < self.min_flow:
                #MV = self.min_flow
                #I = 0 # reset integral solves?
            
            ## update stored data for next iteration
            #e_prev = e
            #t_prev = t



    def start_stop_pump(self):
        # Start pump if not running
        if not self.channelpumping.get_value():
            self.logger.debug("Starting pump...")
            self._cmdStartPumping()
            #self.start_time = time.time()
            return True
        # Stop pump if running
        else:
            self.logger.debug("Stopping pump...")
            self._cmdStopPumping()
            #self.stop_time = time.time()
            return False
            
        #in all cases, return the final status of the pump to change the label of the button
        

    #def read_super_phase(self):
        #"""
        #Returns supervisor phase (CurrentPhase attribute from Beamline 
        #Supervisor Tango DS)
        #@return: str (TRANSFER, SAMPLE, ...)
        #"""
        #return self.phase_channel.getValue().upper()

    #def diff_send_transfer(self):
        #"""
        #Checks if beamline supervisor is in TRANSFER phase.
        #If is not the case, It sends the diff to TRANFER phase.
        #Returns a boolean value indication if the diff is in TRANSFER phase.
        #TODO: Incorporate cstage to transfer phase logic in supervisor
        #@return: boolean
        #"""
        #if self.read_super_phase() == "TRANSFER":
            #logging.getLogger("user_level_log").error(
                #"Supervisor is already in transfer phase"
            #)
            #return True

        #self.go_transfer_cmd()
        #ret = self._wait_phase_done("TRANSFER")
        #return ret

    #def diff_send_serial(self):
        #"""
        #Checks if beamline supervisor is in SAMPLE phase (i.e. sample changer in SAMPLE phase too).
        #If is not the case, It sends the sample changer to SAMPLE phase.
        #Returns a boolean value indication if the sample changer is in SAMPLE phase.
        #@return: boolean
        #"""
        #if self.read_super_phase() == "SERIAL":
            #logging.getLogger("user_level_log").error(
                #"Supervisor is already in serial phase"
            #)
            #return True

        #self.go_serial_cmd()
        #ret = self._wait_phase_done("SERIAL")
        #return ret

    #def _wait_phase_done(self, final_phase):
        #"""
        #Method to wait a phase change. When supervisor reaches the final phase, the diffracometer
        #returns True.
        #@final_phase: target phase
        #@return: boolean
        #"""
        #while True:
            #state = str(self.super_state_channel.getValue())
            #if state == "ON":
                #logging.getLogger("user_level_log").error(
                    #"Supervisor is in ON state. Returning"
                #)
                #break
            #elif str(state) != "MOVING":
                #logging.getLogger("user_level_log").error(
                    #"Supervisor is in a funny state %s" % str(state)
                #)
                #return False

            #logging.getLogger("HWR").debug("Supervisor waiting to finish phase change")
            #time.sleep(0.2)

        #time.sleep(0.1)

        #if self.read_super_phase().upper() != final_phase:
            #logging.getLogger("user_level_log").error(
                #"Supervisor is not yet in %s phase. Aborting load" % final_phase
            #)
            #return False
        #else:
            #return True
            
