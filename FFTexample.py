from optparse import OptionParser
import sys
import heapq, numpy, pylab, time


def generateSine(fs,             # sampling frequency [Hz]
                 fin,            # signal frequency [Hz]
                 amplitude=[1.], # amplitude
                 time=0.,        # signal generation time [sec]
                 N=0,            # number of signal samples
                 noise=0.):      # amplitude of additive white noise
  # simulate time vector
  t = []
  if time > 0:
    # t = arange(start, end, step)
    t = numpy.arange(0, time, 1./fs) # [s]
  elif N > 0:
    t = numpy.arange(0, N) / fs

  fin = numpy.array(fin)
  a = numpy.array(amplitude)
  if len(amplitude) != len(fin):
    if len(amplitude) == 1:
      a = numpy.ones(fin.shape)
    if len(amplitude) > len(fin):
      noise_amp = a[-1]
      a = a[:len(fin)]

  # generate signal
  signal = numpy.zeros(t.shape) # generate 0 vector of lenght 't'
  for i in range(0,len(fin)):
    signal = signal + amplitude[i]*numpy.sin(2*numpy.pi*fin[i]*t)

  # add white noise to the signal
  if noise > 0:
    thenoise = numpy.random.randn(len(t),1)
    signal = signal + thenoise.transpose().flatten()

  return [t, signal]

def FFT(fs,                 # sampling frequency [Hz]
        data,               # sequence to process
        withpadding=False): # add zero padding to FFT (generally only for noisy signals)
  # Remove DC components
  data = data - numpy.mean(data)
# Compute the FFT of the captured signal using an algorithm for the MATLAB helpdesk
  L = len(data)
  NFFT = L
  if withpadding:
    # Get the next power of 2 from the length of the data series
    NFFT=int(2**numpy.ceil(numpy.log2(L)))

# Now use the fft function to compute the DFT of the sequence.
  # Division by L because we want to use the fft to approximate the continuous Fourier transform under the assumption
  # that the signal samples in y were taken at intervals of dt=1/L
  Fourier = numpy.fft.fft(data, NFFT)/L

  # All the following computation assumes that Fourier is a complex array -- ensure this assumption
  if not numpy.iscomplexobj(Fourier):
    # Make real array as complex with zero imaginary components
    Fourier = numpy.array(Fourier, dtype=numpy.complex)

# Use the abs function to obtain the magnitude of the data, the angle
# function to obtain the phase information, and unwrap to remove phase
# jumps greater than pi to their 2*pi complement
  # Magnitude
  magnitude = 2*numpy.abs(Fourier[1:NFFT/2])
# The two-sided amplitude spectrum shows half the peak amplitude at the
# positive and half at the negative frequencies. To convert to the
# single-sided form, multiply each frequency other than DC by two, and
# discard the second half o the array
  # Phase
  phase = numpy.unwrap(numpy.angle(Fourier[1:NFFT/2])) # unwrapping needs to happen in radians!!

  # Frequencies corresponding to FFT
  frequencies = numpy.fft.fftfreq(NFFT, 1./fs)[1:NFFT/2]

  return [frequencies, magnitude, phase]

def generateDisplay(sampletime, # time axis for generated signal [sec]
                    sinewave,   # generated CW sequence
                    spectrum,   # frequency bins of FFT [Hz]
                    magnitude,  # FFT amplitudes
                    phase):     # FFT phases
  pylab.ion()
  pylab.figure()
  pylab.subplot(3,1,1)
  pylab.plot(sampletime, sinewave,'b')
  pylab.ylabel('Sine Wave')
  pylab.subplot(3,1,2)
  # amplitude spectrum = amplitude vs frequency
  # power spectrum = power vs frequency
  # power = 20*numpy.log10(amplitude)
  pylab.plot(spectrum, magnitude, 'b')
  pylab.ylabel('Magnitude [V]')
  pylab.subplot(3,1,3)
  # phase spectrum = phase vs frequency
  pylab.plot(spectrum, phase, 'b')
  pylab.ylabel('Phase [Deg]')
  pylab.xlabel('Frequency [Hertz]')
  pylab.draw()
  # force redraw to get pylab to show data plot
  pylab.draw()
  time.sleep(5)
  pylab.close('all')
  pylab.ioff()

if __name__ == '__main__':

  usage = "\n python %prog [options]"
  parser = OptionParser(usage=usage)
  parser.add_option('-d', '--display',
                    action='store_true',
                    dest='display',
                    default=False,
                    help='Display spectrum components of generated test signals')
  (opts, args) = parser.parse_args()

  # Sample frequency
  fs = 1000. # [Hz]
  # Signal (observation) length
  signaltime = 0.3 # [sec]
  # number of samples
  N = 300 # [fs*time]

## Generate a single frequency sinewave
  # Input signal frequency (or frequency range)
  fin = [50.] # [Hz]
  sys.stdout.write('\nSignal = generateSine(%5.2f, [%4.2f], time=0.3)\n' % (fs, fin[0]))
  [sampletime, sinewave] = generateSine(fs, fin, time=signaltime)
  [spectrum, magnitude, phase] = FFT(fs,sinewave)
  if opts.display:
    generateDisplay(sampletime, sinewave, spectrum, magnitude, phase)
  sys.stdout.write('\tAmplitude spectrum = %4.2f [V]\n' % numpy.max(magnitude))
  sys.stdout.write('\tFundamental frequency from FFT = %4.2f [Hz]\n' % spectrum[numpy.argmax(magnitude)])

## Generate a single frequency sinewave with noise
  fin = [50.] # [Hz]
  sys.stdout.write('\nNoisy signal = generateSine(%5.2f, [%4.2f], time=0.3, noise=0.2)\n' % (fs, fin[0]))
  # [sampletime, sinewave] = generateSine(fs, fin, time=signaltime, noise=0.2)
  [sampletime, sinewave] = generateSine(fs, fin, N=300, noise=0.2)
  [spectrum, magnitude, phase] = FFT(fs,sinewave, withpadding=True) # Padding the FFT of a noisy signal improves the phase information
  if opts.display:
    generateDisplay(sampletime, sinewave, spectrum, magnitude, phase)
  sys.stdout.write('\tNoisy amplitude spectrum = %4.2f [V]\n' % numpy.max(magnitude))
  sys.stdout.write('\tFundamental frequency from FFT = %4.2f [Hz]\n' % spectrum[numpy.argmax(magnitude)])

# Generate a multiple frequency sinewave
  fin = [50., 70.] # [Hz]
  sys.stdout.write('\nSignal multiple frequencies = generateSine(%5.2f, [%4.2f, %4.2f], amplitude=[1., 0.3], time=0.3)\n' % (fs, fin[0], fin[1]))
  [sampletime, sinewave] = generateSine(fs, fin, amplitude=[1., 0.3], time=0.3)
  [spectrum, magnitude, phase] = FFT(fs,sinewave)
  if opts.display:
    generateDisplay(sampletime, sinewave, spectrum, magnitude, phase)
  largest2_value = heapq.nlargest(2, magnitude)
  largest2_index = map(list(magnitude).index, largest2_value)
  sys.stdout.write('\tAmplitude spectrum = [%4.2f, %4.2f] [V]\n' % (largest2_value[0], largest2_value[1]))
  sys.stdout.write('\tFundamental frequency from FFT = [%4.2f, %4.2f] [Hz]\n' % (spectrum[largest2_index[0]], spectrum[largest2_index[1]]))

# Generate a multiple frequency sinewave with noise
  fin = [50., 70.] # [Hz]
  sys.stdout.write('\nSignal multiple frequencies = generateSine(%5.2f, [%4.2f, %4.2f], amplitude=[1., 0.3], time=0.3, noise=0.2)\n' % (fs, fin[0], fin[1]))
  [sampletime, sinewave] = generateSine(fs, fin, amplitude=[1., 0.3], time=0.3, noise=0.2)
  [spectrum, magnitude, phase] = FFT(fs,sinewave, withpadding=True)
  if opts.display:
    generateDisplay(sampletime, sinewave, spectrum, magnitude, phase)
  largest2_value = heapq.nlargest(2, magnitude)
  largest2_index = map(list(magnitude).index, largest2_value)
  sys.stdout.write('\tAmplitude spectrum = [%4.2f, %4.2f] [V]\n' % (largest2_value[0], largest2_value[1]))
  sys.stdout.write('\tFundamental frequency from FFT = [%4.2f, %4.2f] [Hz]\n' % (spectrum[largest2_index[0]], spectrum[largest2_index[1]]))


# -fin-
