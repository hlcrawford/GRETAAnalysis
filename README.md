# GRETAAnalysis

Code compiles to 4 executables:

1) Analyze -- analyzes waveform data, either written via the CS forward buffers or directly from the DM (compile flag chooses which)
2) getTau -- fits waveforms to extract a single tau decay constant
3) getRate -- counts LEDs to extract a rate for a given file (faster than full analysis)
4) readGreta -- analyzes SFB processed events, i.e. energy + time + waveform snippet for all segments in a crystal in one event
