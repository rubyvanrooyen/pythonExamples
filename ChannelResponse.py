import os, re, sys
import numpy, pylab, string, time

import scipy.signal

## Any finite impulse response (FIR) polyphase-DFT filterbank can be modelled by an equivalent windowed DFT
## Display filter frequency response for both overlapping segment and non overlapping segment DFTs
def windowDFT_sinesweep():
  # (0) Initialization and input signals
  sys.stdout.write('Simulate a linear sine sweep from 0 Hz to F1 Hz...\n')
  dt = 0.00001 # Time step (dt = 1/fs the sampling frequency)
  D  = 12. # Signal duration in seconds
  # T1, F1 = The instantaneous frequency F1 is achieved at time T1
  T1 = 1.
  F1 = 30.
  fs = 1./dt # sampling frequency
  L = numpy.ceil(fs*D)+1 # signal duration (samples)
  t = numpy.arange(0, D, dt) # discrete-time axis (seconds)
  # t = numpy.linspace(0, D, D/dt) # discrete-time axis (seconds)
  input_signal = scipy.signal.chirp(t, f0=0, f1=F1, t1=T1, method='linear') # sine sweep from 0 Hz to F1 Hz (<=fs/2 Nyquist)

  # (1) Take the samples of the impulse response of the prototype lowpass filter as the window coefficients
  TotalTaps=4.
  PFBSize=10.
  Channels=float(2**PFBSize) # radix-2 for FFT size
  fwidth=1.
  alltaps = 2**PFBSize * TotalTaps # length of finite-impulse response prototype filter = KNt = (# Channels)(# Taps)
  WindowType='hamming'
  windowval = numpy.hamming(alltaps)
  window_fir_coeffs = windowval * numpy.sinc(fwidth*(numpy.arange(0,alltaps)/float(2**PFBSize)-TotalTaps/2.))

  # (2) Evaluation over the input signal in steps of K samples
  fft_out = numpy.array([])
  NrSegments = numpy.floor((len(input_signal) + (Channels+1))/Channels) - TotalTaps # overlapping segments with windowsize = alltaps
  for i in range(1, int(NrSegments)+1):
    sys.stdout.write('\tprocessing overlapping window %i of %i\n' %(i, NrSegments))
    segment_start = (i-1)*Channels
    segment_end = (TotalTaps + (i-1))*Channels-1
    # (2a) Perform an KNt-point windowed DFT
    windowed_segment = input_signal[segment_start:segment_end+1] * window_fir_coeffs
    fft_column = numpy.abs(numpy.fft.fft(windowed_segment))
    fft_column = fft_column.reshape(1,len(fft_column))/numpy.max(fft_column)
    if len(fft_out) < 1:
      fft_out = fft_column
    fft_out = numpy.vstack((fft_out, fft_column))

  # ABS(X) is the complex modulus (magnitude)
  # |x+iy] = sqrt(x^2+y^2)
  pylab.ion()
  pylab.figure(1)
  pylab.subplot(211)
  pylab.hold(True)
  pylab.plot(10.*numpy.log10(fft_out[:,9:21]))
  pylab.plot(-3.*numpy.ones(fft_out[:,7:11].shape))
  pylab.hold(False)
  pylab.title('Overlapping segments filtered with DFT')
  pylab.ylabel('Power [dBm]')

  # non overlapping segments
  fft_out = numpy.array([])
  for i in range(1, numpy.int(len(input_signal)/alltaps)):
    sys.stdout.write('\tprocessing non-overlapping window %i...\n'% i)
    # windowed_segment = input_signal[i*alltaps : (i+1)*alltaps-1]) * window_fir_coeffs
    windowed_segment = input_signal[i*alltaps : (i+1)*alltaps] * window_fir_coeffs
    fft_column = numpy.abs(numpy.fft.fft(windowed_segment))
    fft_column = fft_column.reshape(1,len(fft_column))/numpy.max(fft_column)
    if len(fft_out) < 1:
      fft_out = fft_column
    fft_out = numpy.vstack((fft_out, fft_column))

  pylab.subplot(212)
  pylab.hold(True)
  pylab.plot(10.*numpy.log10(fft_out[:,9:21]))
  pylab.plot(-3.*numpy.ones(fft_out[:,7:11].shape))
  pylab.hold(False)
  pylab.title('Non-Overlapping segments filtered with DFT')
  pylab.ylabel('Power [dBm]')
  pylab.xlabel('Windowed Channels')
  pylab.savefig('SpectralResponseHamming.png')
  pylab.draw()
  # force redraw to get pylab to show data plot
  pylab.draw()
  time.sleep(5)
  pylab.close('all')
  pylab.ioff()

# All default values assigned will be theoretical approximations
if __name__ == '__main__':

  # Model Channel Response
  windowDFT_sinesweep()

# -fin- 


