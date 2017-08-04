from optparse import OptionParser
import sys
import numpy

##Some general parameters and functions for Analysis
class Parameters():
  def __init__(self, Vpp, TdegC, Fs):
    # Assume sinusoidal input signal
    self.Vpp = Vpp
    self.InputImpedance = 50. # Ohm
    self.N = 8. # bits
    self.SNR = 6.02 * self.N + 1.76 # dB
    self.AmbientTemp = TdegC + 273.15 # K
    self.Boltzmann = 1.3806504e-23 # J/K
    self.Fs = Fs # Hz

  # convert peak voltage to rms volt
  def Vpp2Vrms(self, Vpp):
    return float(Vpp) / (2.*numpy.sqrt(2.)) # V

  # convert degrees cen to kelvin
  def deg2kel(self, TdegC):
    return float(TdegC) + 273.15 # K

  # convert linear power to effective temp in kelvin
  def watt2kel(self, Pwatt, BW):
    return Pwatt / (self.Boltzmann*BW)

  # convert dBm power to linear power
  def dbm2watt(self, Pdbm):
    return 0.001 * 10.**(Pdbm/10.)

  # convert linear power to dBm power
  def watt2dbm(self, Pw):
    return 10. * numpy.log10(Pw/0.001)

  # calculate integrated power
  def rms(self, data):
    return numpy.sqrt(numpy.sum(data**2)/len(data))

##Computing the Effective ADC Noise Figure
#  Reference discussion on --  https://katfs.kat.ac.za/pmwiki/DBE/A2D
def DBENoiseFigure( Bandwidth, Vpp, TdegC, Fs, SNR = -1.):
  sys.stdout.write('The ADC noise floor (NF) is the difference between the thermal and common noise floors expressed relatively to the ADC full-scale range.\n')

  # Create ADC object with default theoretical values
  ADC = Parameters(Vpp, TdegC, Fs)
  # Set input parameters values corretly
  AmbientTemp = ADC.AmbientTemp
  Vrms = ADC.Vpp2Vrms(Vpp) 
  AmbientTemp = ADC.deg2kel(TdegC)
  if SNR < 0: SNR = ADC.SNR

  sys.stdout.write('\t Equation: NFeff [dB] = PFS [dBm] - SNR [dB] - (10*log10(kT) + 30) [dBm] - 10*log10(Fs/2)\n')
  sys.stdout.write('\t\t Full-Scale Power Level (PFS) = ')
  PFS = (10. * numpy.log10(Vrms**2. / ADC.InputImpedance)) + 30. # dBm
  sys.stdout.write('%4.2f dBm \n' %PFS)
  sys.stdout.write('\t\t ADC SNR in first Nyquist bandwidth (SNR) = ')
  sys.stdout.write('%4.2f dB \n' %SNR)
  ##The effect of clock jitter is not implemented here and will be evaluated when the measured SNR values are known.
  # For full-scale input signals the ADC noise floor is mainly increased due to clock jitter. The effects of aperture and sampling clock jitter on an ideal ADCs SNR can be predicted. The rms value for a full-scale input sine wave is:
  # SNRjitter [dB] = -20log10(2pi*fin*tj), where fin is the ADC sine wave input frequency in Hz and tj is the rms jitter in seconds.
  # Thus the effective signal-to-noise ratio of the ADC (SNRADC) is computed:
  # SNR [dB] = - 20log10sqrt[10(-SNR/10) + 10(-SNRjitter/10)), where SNR is due to quantization noise.
  sys.stdout.write('\t\t Calculate thermal noise power at room temperature (10*log10(kT) + 30) = ')
  ThermalNP = 10. * numpy.log10(ADC.Boltzmann*AmbientTemp) + 30 # dBm
  sys.stdout.write('%4.2f dBm \n' %ThermalNP)
  sys.stdout.write('\t\t Processing gain (10*log10(Fs/2))= ')
  ProcessGain = 10. * numpy.log10(Bandwidth) # dB
  sys.stdout.write('%4.2f dB \n' %ProcessGain)

  sys.stdout.write('\t ADC Effective Noise Figure: NFeff = ')
  NFeff = PFS - SNR - ThermalNP - ProcessGain # dBm
  sys.stdout.write('%4.2f dBm \n' %NFeff)


# All default values assigned will be theoretical approximations
if __name__ == '__main__':

  usage = "\n python %prog [options]"

  parser = OptionParser(usage=usage, version="%prog 1.0")
  parser.add_option('--enbw',
		    action='store',
		    dest='enbw',
		    type='float',
		    default=1.,
		    help='Equivalent noise bandwidth.')
  parser.add_option('--channels',
		    action='store',
		    dest='M',
		    type='float',
		    default=0.,
		    help='Number of FFT frequency channels (default=%default).')
  parser.add_option('--Fs',
		    action='store',
		    dest='Fs',
		    type='float',
		    default=800e6,
		    help='ADC sampling frequency (default=%default Hz).')
  parser.add_option('--snr',
		    action='store',
		    dest='snr',
		    type='float',
		    default=-1,
		    help='Signal-to-noise ratio (dB).')
  parser.add_option('--temp',
		    action='store',
		    dest='ambient',
		    type='float',
		    default=25.,
		    help='Room temperature (default = %default degC).')
  parser.add_option('--Vpp',
		    action='store',
		    dest='Vpp',
		    type='float',
		    default=2.,
		    help='Peak-to-peak voltage (default=%default V).')
  (opts, args) = parser.parse_args()

  # Noise Figure Computation
  B = 1 # Hz
  if opts.M > 0:
    B = opts.enbw*(opts.Fs/opts.M)
  Bandwidth = opts.Fs/(2.*B)
  DBENoiseFigure(Bandwidth,          # Hz 
                 opts.Vpp,           # V
                 opts.ambient,       # degC
                 opts.Fs,            # Hz
                 SNR=opts.snr)       # dB

# -fin- 
