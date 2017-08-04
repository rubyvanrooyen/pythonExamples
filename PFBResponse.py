from optparse import OptionParser
import scipy.signal
import os, sys, time
import numpy, pylab

## Any finite impulse response (FIR) polyphase-DFT filterbank can be modelled by an equivalent windowed DFT
## Display filter frequency response for both overlapping segment and non overlapping segment DFTs
def PFB_sinesweep(nChannels, nTaps, nFreq, nPoints, ImName=None, Noise=False):

  # Calculate the coefficients for the prototype lowpass filter. Take the samples of the impulse response of the prototype lowpass filter as
  # the winow coefficients
  fwidth=1.
  alltaps = nChannels * nTaps# length of finite-impulse response prototype filter = KNt = (# Channels)(# Taps)
  windowval = numpy.hamming(alltaps) # Hamming window = Better DFT lowpass
  window_fir_coeffs = windowval * numpy.sinc(fwidth*(numpy.arange(0,alltaps)/nChannels-nTaps/2.))
  # must be the same length as the prototype low-pass FIR filter = nChannels*nTaps
  if (nChannels*nTaps - len(window_fir_coeffs)) > 0:
    window_fir_coeffs = numpy.hstack([window_fir_coeffs, numpy.zeros((1, nChannels*nTaps - len(window_fir_coeffs)))])

  # Reshape the filter coefficients into a matrix whos rows represent the individual polyphase filters to be distributed among the filter bank
  filter_bank = numpy.transpose(numpy.fliplr(numpy.reshape(window_fir_coeffs, (nTaps, nChannels))))

  # Filter a sinusoid that is stepped in frequency from 0 to pi radians
  w = numpy.linspace(0, numpy.pi, nFreq) # frequency vector from 0 to pi
  t = numpy.arange(0, nChannels*nTaps*nPoints) # time vector
  rms_out = numpy.array([]) # store the power of the filtered signal
  for j in range(0,len(w)):
    print j, 'out of', len(w)
    input_signal = numpy.sin(w[j]*t) # signal to filter
    if Noise:
      input_signal = input_signal + numpy.random.rand(len(input_signal))
    # reshape the input so that it represents parallel channels of data going into the filter bank
    if (nChannels*numpy.ceil(len(input_signal)/nChannels) - len(input_signal)) > 0:
      input_signal = numpy.hstack((input_signal, numpy.zeros((nChannels*numpy.ceil(len(input_signal)/nChannels) - len(input_signal), 1))[0]))
    parallel_signals = numpy.transpose(numpy.reshape(input_signal, (len(input_signal)/nChannels, nChannels)))

    # make the output the same size as the input
    fft_out = numpy.zeros(parallel_signals.shape)
    # Apply FIR filter bank
    for channel in range(int(nChannels)):
      channel_filter = filter_bank[channel,:]
      channel_signal = parallel_signals[channel,:]
      fft_out[channel,:] = scipy.signal.lfilter(channel_filter, 1, channel_signal)

    # FFT
    FFT_out = numpy.fft.fft(fft_out, axis=0)
    FFT_nw = numpy.fft.fft(parallel_signals, axis=0)
    # Store the output power
#     FFT_rms = numpy.sqrt(numpy.mean(numpy.power(numpy.abs(FFT_out),2), axis=1)) # RMS voltage
#     FFT_rms_nw = numpy.sqrt(numpy.mean(numpy.power(numpy.abs(FFT_nw),2), axis=1)) # RMS voltage
    FFT_rms = numpy.mean(numpy.abs(FFT_out), axis=1) # Average amplitude
    FFT_rms_nw = numpy.mean(numpy.abs(FFT_nw), axis=1) # Average amplitude
    if len(rms_out) < 1:
      rms_out = FFT_rms.reshape((len(FFT_rms),1))
      rms_nw = FFT_rms_nw.reshape((len(FFT_rms_nw),1))
    else:
      rms_out = numpy.hstack((rms_out, FFT_rms.reshape((len(FFT_rms),1))))
      rms_nw = numpy.hstack((rms_nw, FFT_rms_nw.reshape((len(FFT_rms_nw),1))))

  # Normalize values to get 0dB values at max
  for k in range(int(nChannels)):
    rms_out[k,:] = rms_out[k,:]/numpy.max(rms_out[k,:])
    rms_nw[k,:] = rms_nw[k,:]/numpy.max(rms_nw[k,:])

  start_idx = 1
  end_idx = -1
  rms_out_dB = 20.*numpy.log10(rms_out[:,start_idx:end_idx])
  rms_nw_dB = 20.*numpy.log10(rms_nw[:,start_idx:end_idx])

  pylab.ion()
  pylab.figure(1)
  channel_num = int(nChannels/4)
  pylab.plot(w[start_idx:end_idx],rms_out_dB[channel_num,:], 'r', w[start_idx:end_idx],rms_nw_dB[channel_num,:], 'm')
  pylab.legend(('PFB','FFT'))
#   pylab.title('Filter Bank Frequency Response - RMS power')
  pylab.title('Filter Bank Frequency Response - Average Amplitude')
  pylab.xlabel('Frequency (normalized to channel center)')
  pylab.ylabel('Magnitude Response (dB)')
  pylab.xticks(numpy.arange((nChannels+1)/2)*w[-1]/nChannels*2, numpy.arange((nChannels+1)/2)-nChannels/4)
  pylab.axis('tight')
  pylab.draw()
  # force redraw to get pylab to show data plot
  pylab.draw()
  time.sleep(5)
  if ImName:
    pylab.savefig('FilterResponse_'+ImName, dpi=300)
    print 'Output image as:', 'FilterResponse_'+ImName
  pylab.close()

  pylab.ion()
  pylab.figure(1)
  pylab.hold(True)
  if Noise:
    for k in range(2,int(nChannels/2)):
      pylab.plot(w[start_idx:end_idx],rms_out_dB[k,:])
  else:
    for k in range(int(nChannels)):
      pylab.plot(w[start_idx:end_idx],rms_out_dB[k,:])
  pylab.axhline(y=-6, color='r')
  pylab.hold(False)
#   pylab.title('Filter Bank Frequency Response - RMS power = sqrt(mean(abs(FFT).^2))')
  pylab.title('Filter Bank Frequency Response - Average Amplitude')
  pylab.xlabel('Frequency (normalized to channel center)')
  pylab.ylabel('Magnitude Response (dB)')
  pylab.xticks(numpy.arange((nChannels+1)/2)*w[-1]/nChannels*2, numpy.arange((nChannels+1)/2)-nChannels/4)
  pylab.axis('tight')
  pylab.draw()
  # force redraw to get pylab to show data plot
  pylab.draw()
  time.sleep(5)
  if ImName:
    pylab.savefig('WideBand_'+ImName, dpi=300)
    print 'Output image as:', 'WideBand'+ImName
  pylab.close()


# All default values assigned will be theoretical approximations
if __name__ == '__main__':

  usage = "\n python %prog [options] \
\nExamples: \
\n\tDesign your own filterbank:\
\n\t\tpython %prog --channels=32 --taps=8 --resolution=500 --freq=1600 \
\n\tFilter response of noisy sine wave and output graphs as images:\
\n\t\tpython %prog -n -o output"

  parser = OptionParser(usage=usage)
  parser.add_option('--pfb',
                    action='store',
                    dest='pfb_size',
                    type='float',
                    default=None,
                    help='Specify the PFB size. Not to be confused with the number of channels, since it does not have to be radix 2')
  parser.add_option('--channels',
                    action='store',
                    dest='n_channels',
                    type='float',
                    default=None,
                    help='Number of channels in the filter bank (radix-2 for FFT size)')
  parser.add_option('--taps',
                    action='store',
                    dest='n_taps',
                    type='string',
                    default='4',
                    help='Number of taps in each FIR filter (default=%default). More than one number can be specified: 4,8')
  parser.add_option('--freqs',
                    action='store',
                    dest='n_freqs',
                    type='float',
                    default=200,
                    help='Number of frequencies to sweep (default=%default)')
  parser.add_option('--resolution',
                    action='store',
                    dest='n_points',
                    type='float',
                    default=100,
                    help='Number of output points from each channel (default=%default)')
  parser.add_option('-o', '--output',
                    action='store',
                    dest='imagename',
                    type='string',
                    default=None,
                    help='Output channel response displays to .png images')
  parser.add_option('-n', '--noise',
                    action='store_true',
                    dest='noise',
                    default=False,
                    help='Add Gaussian noise to test Sine input signal')
  (opts, args) = parser.parse_args()

  # Initialize the variables to define the filters and the filter bank
  PFBSize = 4.
  nChannels = float(2**PFBSize)
  if opts.n_channels:
    nChannels = opts.n_channels
  elif opts.pfb_size:
    nChannels = float(2**opts.pfb_size)
  taps = []
  try:
    taps = numpy.array(opts.n_taps.split(','), dtype=numpy.float)
  except ValueError:
    print "String of taps should be without spaces: e.g. --taps=4,8"
    sys.exit(1)
  nFreq = opts.n_freqs
  nPoints = opts.n_points
  ImName = None
  if opts.imagename:
    ImName = os.path.splitext(opts.imagename)[0]+'.png'
  for nTaps in taps:
    # Model Channel Response
    PFB_sinesweep(nChannels, nTaps, nFreq, nPoints, ImName = ImName, Noise = opts.noise)

# -fin-


