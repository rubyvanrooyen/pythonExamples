#!/usr/bin/python

import serial, socket, string, time

class SCPI:
  PORT = 5025
  BAUDRATE = 9600

  ## Connect to the R&S signal generator
  def __init__(self,
               host=None, port=PORT,           # set up socket connection
               device=None, baudrate=BAUDRATE, # set up serial port
               timeout=1):
    if host and device:
      raise RuntimeError('Only one connection can be initaited at a time.\nSelect socket or serial connection.\n')

    # Ethernet socket connection
    self.connection = None
    if host:
      self.connection = 'socket'
      self.s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
      self.s.connect((host, port))
      self.s.settimeout(1)
    elif device:
      self.connection = 'serial'
      self.s = serial.Serial(device, baudrate, timeout, rtscts=0)
    else:
      raise RuntimeError('No connections specified.\n')

  # Querie instrument identificaton
  def deviceID(self):
    self.write("*IDN?")
    print "DEVICE: " + self.read()

  # send query / command via relevant port comm
  def write(self, command):
    if self.connection == 'serial':
      self.s.write(command + '\r\n')
    else:
      self.s.send(command+ '\n')
    time.sleep(1)
  def read(self, bufsize=128):
    if self.connection == 'serial':
      return self.s.readline()
    else:
#       return self.s.recv(128)
      return self.s.recv(bufsize)

  # reset
  def reset(self):
    self.write("*RST")
    self.write("*CLS")

  def fsuSetup(self):
    self.write("SYST:DISP:UPD ON")
    # set the ampl ref level 5dBm
    self.write("DISP:WIND:TRAC:Y:RLEV 5dBm")
    # set the rf attenuation to 0dB
    self.write("INP:ATT %.3f dB" % (0.0,))
    # set the resolution bandwidth to 300KHz
    self.write("BAND:AUTO OFF")
    self.write("BAND %.3fKHz" % (300,))
    # set the video bandwidth to 300KHz
    self.write("BAND:VID:AUTO OFF")
    self.write("BAND:VID %.3fKHz" % (300,))
    # set frequency range, start at 0MHz and stop at 400MHz
    self.write("FREQ:STAR %.2fMHz"%(0,))
    self.write("FREQ:STOP %.2fMHz"%(400,))
    # set sweep time to 15ms
    self.write("SWE:TIME:AUTO OFF")
    self.write("SWE:TIME %.3f s" % (0.015,))

  def noisePwr(self):
    self.write("INIT:CONT OFF") # single sweep
    # select absolute power measurement
    self.write("SENS:POW:ACH:MODE ABS")
    # set channel bandwidth of the main transmission channel
    self.write("SENS:POW:ACH:ACP 0")
    self.write("SENS:POW:ACH:BWID %.3fMHz" %(256,))
    # switch on power measurements
    self.write("CALC:MARK:FUNC:POW ON")
    # calculation of channel power
    self.write("CALC:MARK:FUNC:POW:SEL CPOW")
    self.write("INIT;*WAI") # start and wait to complete the sweep
    self.write("CALC:MARK:FUNC:POW:RES? CPOW")
    return float(self.read(bufsize=32768))


  def cwPeak(self):
    # determine the max peak
    self.write("INIT:CONT OFF") # single sweep
    self.write("CALC:MARK ON")  # switch on marker 1
    self.write("CALC:MARK:COUN ON") # switch on frequency counter for marker 1
    self.write("INIT:IMM;*WAI") # start and wait to complete the sweep
    self.write("CALC:MARK:MAX:AUTO ON") # switch on automatic peak search for marker 1
    self.write("CALC:MARK:COUN:FREQ?") # outputs the measured value
    freq_hz = float(self.read(bufsize=32768))
    self.write("CALC:MARK:Y?")# outputs the measured value for marker 1
    amp_dbm = float(self.read(bufsize=65536))
    return [freq_hz, amp_dbm]
  # close the comms port to the R&S signal generator
  def __close__(self):
    self.s.close()

if __name__ == '__main__':

  # SMB100A R&S Signal Generator IP address
  fsu_ip='192.168.1.123'
  socket_port=5025
## Using SCPI class for comms to FSU R&S Spectrum Analyser
  print 'Connecting to device IP:', fsu_ip
  fsu_obj=SCPI(fsu_ip)

  try:
    fsu_obj.deviceID()
    fsu_obj.reset()
    # setup FSU SA
    fsu_obj.fsuSetup()
    # measure channel power
    channel_power = fsu_obj.noisePwr()
    print "Channel power: %f dBm" % channel_power
    # measure the CW ampl
    [freq_hz, amp_dbm] = fsu_obj.cwPeak()
    print "CW at %f MHz has amplitude %f dBm" % (freq_hz/1e6, amp_dbm)
  except: # Exception as e:
    print e
    pass # need to close port

  print 'Closing all ports...'
  try:
    fsu_obj.__close__()
  except:
    pass # socket already closed

# -fin-

