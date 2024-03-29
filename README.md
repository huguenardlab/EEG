# EEG
## openephysreader is code for marking seizures so that they can later be analyzed, e.g., power spectra and coherence.
## OpenEphysLoader is a wrapper function, by JMS, which in turn calls files available from https://github.com/open-ephys/analysis-tools
## harmonic_filter2.m (uploaded 3/28/2024) will attempt to remove line-frequency noise (i.e. 60 Hz hum) from a signal by fitting a sine wave to the signal to estimate amplitude and phase of the line noise, and then subtract this from the original signal. When this works well, it does not distort the signal like a narrow band bandpass filter (notch filter) does.
