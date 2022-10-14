function [ confidenceMetric ] = computeConfidenceMetric2( data, Fs, fft_size )
%% data是相位差分数据
Fs = 50;
fs = Fs;
FFT_size = fft_size;
freqBreathTop = 0.6;   % 先验知识
freqBreathLow = 0.1;
spectrumIndexStart = floor(freqBreathLow*FFT_size/fs) + 1;
spectrumIndexEnd = ceil(freqBreathTop*FFT_size/fs) + 1;

BreathSpectrum = fft(data,FFT_size);
BreathAbsSpectrum = abs(BreathSpectrum);
BreathAbsSpectrum = BreathAbsSpectrum.^2;

sumSignal = sum(BreathAbsSpectrum(1:FFT_size/2+1));
sumPeak = sum(BreathAbsSpectrum(spectrumIndexStart:spectrumIndexEnd));

confidenceMetric = sumPeak/sumSignal;

end

