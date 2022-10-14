function [ result1, result2 ] = old_optAlgorithm( GUIdata, name, selected_index )
%% 承接旧日的两个算法
%% 解析数据
rangeProf = reshape(GUIdata.rangeProf,22,4,size(GUIdata.numFrame,2));
numFrame = GUIdata.numFrame;
subName = name;
selectFrame = 1:length(numFrame);
numFrames = length(selectFrame);

rangeBinStartIndex = 5;
rangeBinEndIndex = 22;

focus_theta = [0,20,40];
rangeBinIndexPhase = 10;
numSelectBins = 16;

rangeBin_val = zeros(rangeBinEndIndex,numFrames);
pTemp_Prev = zeros(4,22);
pRangeProfileClutterRemoved = zeros(22,4,numFrames);
rangeBinPhase = zeros(1,numFrames);

dataPrev1 = zeros(1,numSelectBins);
dataPrev2 = zeros(1,numSelectBins);
dataCurr = zeros(1,numSelectBins);

phasePrevFrame = zeros(1,numSelectBins);
delayBreath_HF_Phase = zeros(numSelectBins,6*3+2);
tempPhasePre = zeros(1,numSelectBins);
outPhasePerBin = zeros(1,numSelectBins);
phaseChaPerBin = zeros(2,numSelectBins);

bin_val = zeros(16,numFrames);

outPhase = zeros(1,numFrames);
binSel_record = zeros(1,numFrames);

guiFlag_Clutteremoval = 0;

%breath filter parameters
delayBreath = zeros(numSelectBins,6*2+2);
filtercoefsBreath = [1.0000, 0, -1.0000, 1.0000, -1.9835, 0.9837;
    1.0000, 0, -1.0000, 1.0000, -1.8898, 0.8977];
scalevalsBreath = [0.0427, 0.0427, 1.0000];
numstagesBreath = 2;

%breath high filter parameters
delayBreath_HF = zeros(numSelectBins,6*3+2);
filtercoefsBreath_HF = [1.0000, 0, -1.0000, 1.0000, -1.988050715235064, 0.988214663473026;
    1.0000, 0, -1.0000, 1.0000, -1.743792905978169, 0.797913309708027;
    1.0000, 0, -1.0000, 1.0000, -1.782959838965117, 0.785792576664767];
scalevalsBreath_HF = [0.112624285734591, 0.112624285734591, 0.107103711667616, 1.0000];
numstagesBreath_HF = 3;

%heart filter parameters
delayHeart = zeros(1,6*4+2);
filtercoefsHeart = [1.0000, 0, -1.0000, 1.0000, -0.5306, 0.5888;
    1.0000, 0, -1.0000, 1.0000, -1.8069, 0.8689;
    1.0000, 0, -1.0000, 1.0000, -1.4991, 0.5887;
    1.0000, 0, -1.0000, 1.0000, -0.6654, 0.2099];
scalevalsHeart = [0.4188, 0.4188, 0.3611, 0.3611, 1];
numstagesHeart = 4;

% 差分数据的毛刺滤波
Num = [0.0853386857109242,0.0886216297925902,0.0912405722095520,0.0931469621823841,...
    0.0943052579023246,0.0946937844044499,0.0943052579023246,0.0931469621823841,...
    0.0912405722095520,0.0886216297925902,0.0853386857109242];
numFirCoef = length(Num);
xValueBuffer = zeros(numSelectBins,numFirCoef);

range_angle_binVal = zeros(1,16);
bufferSizeForCenter = 100;


%Breath_CircularBuffer
idxBinSelected = 0;
circularBufferSizeBreath = 512 + numSelectBins;
breathWaveDiff_CircularBuffer = zeros(numSelectBins,circularBufferSizeBreath);
breathWave_CircularBuffer =zeros(numSelectBins,circularBufferSizeBreath);
breathWaveDiffAfterFilter = zeros(numSelectBins,circularBufferSizeBreath);

BreathWaveDiffConfidenceMetric = zeros(1,numSelectBins);
BreathWaveDiffConfidenceMetric_Record = zeros(numSelectBins,numFrames);
breathWaveDiffAfterFilter_record = zeros(numSelectBins,numFrames);
confidenceMetricAveDiffReord = zeros(numSelectBins,8);


%% circleFitting
circleFittingSize = 200;
circleFittingCenter = zeros(1, numSelectBins);
circleFittingCenter_log = zeros(numSelectBins, numFrames);

%% 毛刺噪声能量统计
medianSize = 7;
ImpluseNoiseSize = 256;
phaseCircluarBuffer = zeros(numSelectBins, ImpluseNoiseSize);
ImpluseNoiseBufferBuffer = zeros(numSelectBins, ImpluseNoiseSize);
ImpluseNoiseEnergy = zeros(numSelectBins, 1);
phaseTemp = zeros(numSelectBins, 1);
phaseTemp_log = zeros(numSelectBins, numFrames);
instantImpluseNoiseEnergy = zeros(numSelectBins, 1);
ImpluseNoiseEnergy_Log = zeros(numSelectBins, numFrames);
instantImpluseNoiseEnergy_Log = zeros(numSelectBins, numFrames);

% cValue = zeros(numSelectBins,1);
% cValueRecord = zeros(numSelectBins,numFrames);
% maxPeakVal = zeros(numSelectBins,1);
% maxPeakValRecord = zeros(numSelectBins,numFrames);

FLAG_REMOVE_IMPULSE_NOISE = 1;

%% 使用一阶RC去除原数据中的DC成分
fs = 50;
fc = 0.1;
fc2 = 0.01;
Ts = 1/fs;
alphaClutter = Ts/(Ts+1/(2*pi*fc));
breathFilterRecord = zeros(numSelectBins,numFrames);
baselineLength = 200;
outPhasePre = 0;
idxTemp = 0;
flag_bad_frame_record = zeros(1,numFrames);

%%
% 处理主体
for gFrameCount = selectFrame(1):selectFrame(end)
    % for gFrameCount = 1:1:500
    idxTemp = idxTemp + 1;
    rangeBinMax = 0;
    maxVal = 0;
    rangeBinMaxClutter = 0;
    maxValClutter = 0;
    
    % find rangeBinIndexPhase
    for rangeBinIndex = rangeBinStartIndex:rangeBinEndIndex
        if guiFlag_Clutteremoval == 1
            temp_Curr = squeeze(rangeProf(rangeBinIndex,:,gFrameCount)).';
            currRangeIndex = rangeBinIndex - rangeBinStartIndex+1;
            pTemp_Prev(:,currRangeIndex) = alphaClutter*temp_Curr + (1-alphaClutter)*pTemp_Prev(currRangeIndex);
            currVal = abs(sum(temp_Curr - pTemp_Prev(:,currRangeIndex)));
            pRangeProfileClutterRemoved(rangeBinIndex,:,idxTemp) = temp_Curr - pTemp_Prev(:,currRangeIndex);
            
            if currVal > maxValClutter
                maxValClutter = currVal;
                rangeBinMaxClutter = rangeBinIndex;
            end
        else
            temp_Curr = rangeProf(rangeBinIndex,1,gFrameCount);
            rangeBin_val(rangeBinIndex, idxTemp) = temp_Curr;
            absVal = abs(temp_Curr);
            if absVal > maxVal
                maxVal = absVal;
                rangeBinMax = rangeBinIndex;
            end
        end
    end
    
    if mod(idxTemp,128) == 1
        if guiFlag_Clutteremoval == 1
            rangeBinIndexMaxPhase = rangeBinMaxClutter;
        else
            rangeBinIndexMaxPhase = rangeBinMax;
        end
    end
    
    for rangeBinIndex = 1:numSelectBins
        idx = floor((rangeBinIndex-1)/3);
        theta = focus_theta(mod(rangeBinIndex-1, 3)+1);
        
        if rangeBinIndex == numSelectBins
            if rangeBinIndexMaxPhase == 0
                select_bin = 10;
            else
                select_bin = rangeBinIndexMaxPhase;
            end
            temp_bin_val = rangeProf(select_bin,1,gFrameCount);
        else
            select_bin = rangeBinIndexPhase + idx;
            temp_bin_val = 0;
            w = single(hanning(4));
            for idx = 1:4
                A = exp(-1j*pi*(idx-1)*sind(theta))*w(idx);
                temp_bin_val = temp_bin_val + (A*rangeProf(select_bin,idx,gFrameCount));
            end
        end
        
        if isnan(temp_bin_val)
            flag_bad_frame_record(idxTemp) = 1;
            temp_bin_val = bin_val(rangeBinIndex,idxTemp-1);
        end
        
        %% 对相位的差分数据进行去毛刺滤波
        range_angle_binVal(rangeBinIndex) = temp_bin_val;
        bin_val(rangeBinIndex,idxTemp) = temp_bin_val;
        
        temp_bin_val = temp_bin_val - circleFittingCenter(rangeBinIndex);
        rangeBinPhase = atan2(real(temp_bin_val), imag(temp_bin_val));
        phaseUsedComputation = wrapToPi(rangeBinPhase - phasePrevFrame(rangeBinIndex));
        phasePrevFrame(rangeBinIndex) = rangeBinPhase;
        
        %% 统计毛刺噪声能量
        for idxloop = 1:ImpluseNoiseSize-1
            phaseCircluarBuffer(rangeBinIndex,idxloop) = phaseCircluarBuffer(rangeBinIndex,idxloop+1);
        end
        phaseCircluarBuffer(rangeBinIndex, ImpluseNoiseSize) = phaseUsedComputation;
        
        tempWave = phaseCircluarBuffer(rangeBinIndex, ImpluseNoiseSize-6:ImpluseNoiseSize);
        phaseDifftemp = median(tempWave);
        phaseTemp(rangeBinIndex) = sum(phaseCircluarBuffer(rangeBinIndex,:));
        phaseTemp_log(rangeBinIndex,idxTemp) =  abs(phaseTemp(rangeBinIndex));
        outDiffval = tempWave(floor(medianSize/2)+1) - phaseDifftemp;
        
        idxforBuffer = mod(idxTemp-1, ImpluseNoiseSize)+1;
        ImpluseNoiseBufferBuffer(rangeBinIndex,idxforBuffer) = outDiffval;
        ImpluseNoiseEnergy(rangeBinIndex) = sum(ImpluseNoiseBufferBuffer(rangeBinIndex,:).^2);
        instantImpluseNoiseEnergy(rangeBinIndex) = outDiffval*outDiffval;
        
        ImpluseNoiseEnergy_Log(rangeBinIndex,idxTemp) = ImpluseNoiseEnergy(rangeBinIndex);
        instantImpluseNoiseEnergy_Log(rangeBinIndex,idxTemp) = instantImpluseNoiseEnergy(rangeBinIndex);
        
        %% 毛刺滤波
        if FLAG_REMOVE_IMPULSE_NOISE
            dataPrev2(rangeBinIndex) = dataPrev1(rangeBinIndex);
            dataPrev1(rangeBinIndex) = dataCurr(rangeBinIndex);
            dataCurr(rangeBinIndex) = phaseUsedComputation;
            phaseUsedComputation = Filter_RemoveImpluseNoise(dataPrev2(rangeBinIndex), dataPrev1(rangeBinIndex), dataCurr(rangeBinIndex), 1.5);
        end
        
        phaseChaPerBin(1,rangeBinIndex) = phaseUsedComputation;
        outPhasePerBin(rangeBinIndex) = phaseUsedComputation + outPhasePerBin(rangeBinIndex);
        
        %% IIR+FIR滤波
        tempPhase = alphaClutter*outPhasePerBin(rangeBinIndex) + (1-alphaClutter)*tempPhasePre(rangeBinIndex);
        tempPhasePre(rangeBinIndex) = tempPhase;
        % breathWaveDiffAfterFilter_record(rangeBinIndex,idxTemp) = tempPhase;
        tempPhase = outPhasePerBin(rangeBinIndex) - tempPhase;
        for idxBuffer = 1:numFirCoef-1
            xValueBuffer(rangeBinIndex,idxBuffer) =  xValueBuffer(rangeBinIndex,idxBuffer+1);
        end
        xValueBuffer(rangeBinIndex,numFirCoef) = tempPhase;
        yout = xValueBuffer(rangeBinIndex,:)*Num';
        for i = 2:circularBufferSizeBreath
            breathWaveDiffAfterFilter(rangeBinIndex,i-1) = breathWaveDiffAfterFilter(rangeBinIndex, i);
        end
        
        breathWaveDiffAfterFilter(rangeBinIndex, circularBufferSizeBreath) = yout;
        breathWaveDiffAfterFilter_record(rangeBinIndex,idxTemp) = yout;
        
        %% 对波形数据的IIR滤波（0.1-0.6Hz）
        delayBreathTempForPhase = delayBreath_HF_Phase(rangeBinIndex,:);
        [outputFilterBreathOutForPhase, delayBreathTempForPhase] = Filter_IIR_BiquadCascade(outPhasePerBin(rangeBinIndex), filtercoefsBreath, scalevalsBreath, delayBreathTempForPhase, numstagesBreath);
        delayBreath_HF_Phase(rangeBinIndex,:) = delayBreathTempForPhase;
        
        for i = 2:circularBufferSizeBreath
            breathWave_CircularBuffer(rangeBinIndex,i-1) = breathWave_CircularBuffer(rangeBinIndex, i);
        end
        breathWave_CircularBuffer(rangeBinIndex, circularBufferSizeBreath) = outputFilterBreathOutForPhase;
        
        %% 对波形差分数据的IIR滤波（0.1-2Hz）
        delayBreathTemp = delayBreath_HF(rangeBinIndex,:);
        [outputFilterBreathOut, delayBreathTemp] = Filter_IIR_BiquadCascade(phaseUsedComputation, filtercoefsBreath_HF, scalevalsBreath_HF, delayBreathTemp, numstagesBreath_HF);
        delayBreath_HF(rangeBinIndex,:) = delayBreathTemp;
        
        for i = 2:circularBufferSizeBreath
            breathWaveDiff_CircularBuffer(rangeBinIndex,i-1) = breathWaveDiff_CircularBuffer(rangeBinIndex, i);
        end
        breathWaveDiff_CircularBuffer(rangeBinIndex, circularBufferSizeBreath) = int32(100000*outputFilterBreathOut);
        breathFilterRecord(rangeBinIndex, idxTemp) = outputFilterBreathOut;
        
    end
    
    % Call FFT in a loop
    idxBinSelected  = idxBinSelected + 1;
    if idxBinSelected > numSelectBins
        idxBinSelected = 1;
    end
    
    %% 对每个bin做圆拟合 并判断是否更新中心点
    if(idxTemp > 512)
        tempBinVal = bin_val(idxBinSelected,(idxTemp-511):idxTemp);
        circleFittingDataBuffer = tempBinVal(1:2:end);
        [cenTemp, radius] = circleFit2(circleFittingDataBuffer);
        
        tempWave = wrapToPi(diff(angle(tempBinVal - cenTemp)));
        confidenceMetricTemp = computeConfidenceMetric2(tempWave,fs,512);
        
        cenPrev = circleFittingCenter(idxBinSelected);
        prevWave = wrapToPi(diff(angle(tempBinVal - cenPrev)));
        confidenceMetricPrev = computeConfidenceMetric2(prevWave,fs,512);
        
        if(confidenceMetricPrev < confidenceMetricTemp || confidenceMetricTemp > 0.35)
            circleFittingCenter(idxBinSelected) = cenTemp;
        end
        
        circleFittingCenter_log(:,idxTemp) = circleFittingCenter.';
    end
    
    numIndexAroundPeak = 3;
    idxStr = 17;
    idxEnd = circularBufferSizeBreath;
    
    tempDataBreathWave = breathWaveDiff_CircularBuffer(idxBinSelected,idxStr:idxEnd);
    BreathWaveDiffConfidenceMetric(idxBinSelected) = computeConfidenceMetric(tempDataBreathWave,numIndexAroundPeak);
    
    for idx = 1:7
        confidenceMetricAveDiffReord(idxBinSelected,idx) = confidenceMetricAveDiffReord(idxBinSelected,idx+1);
    end
    confidenceMetricAveDiffReord(idxBinSelected,8) = BreathWaveDiffConfidenceMetric(idxBinSelected);
    BreathWaveDiffConfidenceMetric_Record(:,idxTemp) = mean(confidenceMetricAveDiffReord,2);
    
    confidenceMetric = mean(confidenceMetricAveDiffReord,2);
    
    if idxTemp <= 256
        binSel = floor(GUIdata.outIndex(gFrameCount)/10) + 1;
        %         binSel = numSelectBins;
    else
        
        %% 计算每个bin的呼吸置信系数——以自身及相邻bin的权重 毛刺噪声能量
        if mod(numFrame(gFrameCount), 128) == 1
            idxBinSelected = 0;
            BreathCMRecord = zeros(1,numSelectBins);
            
            for idxX = 1:4
                BreathCMRecord((idxX-1)*3+1) = 0.75*confidenceMetric((idxX-1)*3+1)+0.35*confidenceMetric((idxX-1)*3+2)+0.25*confidenceMetric((idxX-1)*3+4)+0.15*confidenceMetric((idxX-1)*3+5);
                BreathCMRecord((idxX-1)*3+2) = 0.75*confidenceMetric((idxX-1)*3+2)+0.35*confidenceMetric((idxX-1)*3+3)+0.15*confidenceMetric((idxX-1)*3+4)+0.25*confidenceMetric((idxX-1)*3+5);
                BreathCMRecord((idxX-1)*3+3) = 0.75*confidenceMetric((idxX-1)*3+3)+0.35*confidenceMetric((idxX-1)*3+2)+0.25*confidenceMetric((idxX-1)*3+5)+0.15*confidenceMetric((idxX-1)*3+6);
            end
            BreathCMRecord(13) = 0.9*confidenceMetric(13)+0.6*confidenceMetric(14);
            BreathCMRecord(14) = 0.9*confidenceMetric(14)+0.6*confidenceMetric(15);
            BreathCMRecord(15) = 0.9*confidenceMetric(15)+0.6*confidenceMetric(14);
            
            [maxCM,idxMaxCM] = max(BreathCMRecord);
            
            if ((maxCM > 0.65)  && ((BreathCMRecord(binSel) < 0.45 ) || confidenceMetric(binSel) < 0.3 ))
                [~,idxCMSort] = sort(BreathCMRecord,'descend');
                numCandidata = 4;
                binSel = idxCMSort(1);
                tempVal = ImpluseNoiseEnergy(idxCMSort(1));
                for idx = 2:numCandidata
                    if confidenceMetric(idxCMSort(idx)) > maxCM/1.5-0.1
                        if ImpluseNoiseEnergy(idxCMSort(idx)) < tempVal
                            binSel = idxCMSort(idx);
                            tempVal = ImpluseNoiseEnergy(idxCMSort(idx));
                        end
                    end
                end
            end
        end
    end
    
    binSel_record(1,idxTemp) = binSel;
    
    %% 以速度来计算增量
    velocityfftSize = 512;
    lambda = 5;
    dataVelocitySize = 16;
    vel = (-velocityfftSize/2:velocityfftSize/2-1)/velocityfftSize*lambda/2*fs;
    if idxTemp > 16
        dataIn = bin_val(binSel, idxTemp-dataVelocitySize+1:idxTemp) - circleFittingCenter(binSel);
        pSpectrum = fftshift(fft(dataIn,velocityfftSize));
        if(sum(abs(pSpectrum).^2) < 512*512*500 && confidenceMetric(binSel) < 0.3)
            velCompute = 0;
        else
            velCompute = vel(2:end)*abs(pSpectrum(2:end))'/sum(abs(pSpectrum(2:end)))/50/5*4*pi;
        end
    else
        velCompute = 0;
    end
    outPhase(idxTemp) = outPhasePre - velCompute;
    outPhasePre = outPhase(idxTemp);
    %
    disp([subName,',老算法处理帧数: ',num2str(gFrameCount)]);
end
result1 = GUIdata.outHeart / 8000;
result2 = outPhase;

end


function output = Filter_RemoveImpluseNoise(dataPrev2, dataPrev1, dataCurr, thresh)
backwardDiff = dataPrev1 - dataPrev2;
forwardDiff = dataPrev1 - dataCurr;
x1 = 0;
x2 = 2;
y1 = dataPrev2;
y2 = dataCurr;
x = 1;
if(((forwardDiff>thresh)&&(backwardDiff>thresh))||((forwardDiff<-thresh)&&(backwardDiff<-thresh)))
    y = y1 + (((x-x1)*(y2-y1))/(x2-x1));
else
    y = dataPrev1;
end
output = y;
end

function [output, delayout] = Filter_IIR_BiquadCascade(DataIn, filtercoefs, scalevals, delay, numStages)

input = DataIn;
for indexstage = 1:numStages
    indextemp = 3*(indexstage-1)+1;
    delay(indextemp) = scalevals(indexstage)*input-filtercoefs(indexstage,5)*delay(indextemp+1)-filtercoefs(indexstage,6)*delay(indextemp+2);
    y = filtercoefs(indexstage,1)*delay(indextemp)+filtercoefs(indexstage,2)*delay(indextemp+1)+filtercoefs(indexstage,3)*delay(indextemp+2);
    delay(indextemp+2) = delay(indextemp+1);
    delay(indextemp+1) = delay(indextemp);
    input = y;
end
delayout = delay;
output = y;
end

function confidenceMetric = computeConfidenceMetric(data, numIndexAroundPeak)

FFT_size = 512;
fs = 50;
freqBreathTop = 0.6;
freqBreathLow = 0.1;
spectrumIndexStart = floor(freqBreathLow*FFT_size/fs) + 1;
spectrumIndexEnd = ceil(freqBreathTop*FFT_size/fs) + 1;

BreathSpectrum = fft(data,FFT_size);
BreathAbsSpectrum = abs(BreathSpectrum);
BreathAbsSpectrum = BreathAbsSpectrum.^2;
[pks,locs] = findpeaks(BreathAbsSpectrum(1:spectrumIndexEnd));
if isempty(pks)
    [~,peakIndex] = max(BreathAbsSpectrum(spectrumIndexStart:spectrumIndexEnd));
    peakIndex = peakIndex + spectrumIndexStart - 1;
    startInd = peakIndex - numIndexAroundPeak;
    if startInd <= 1
        startInd = 2;
    end
    sumBreath = sum(BreathAbsSpectrum(spectrumIndexStart:spectrumIndexEnd));
    sumSignal = sum(BreathAbsSpectrum(1:FFT_size/2+1));
    confidenceMetric = sumBreath/sumSignal;
else
    [~,maxIndex] = max(pks);
    peakIndex = locs(maxIndex);
    peakIndex = peakIndex + spectrumIndexStart - 1;
    startInd = peakIndex - numIndexAroundPeak;
    endInd = peakIndex + numIndexAroundPeak;
    if startInd <= 1
        startInd = 2;
    end
    for idx = 1:numIndexAroundPeak
        tempIdx = peakIndex+idx;
        if BreathAbsSpectrum(tempIdx) > BreathAbsSpectrum(peakIndex)
            endInd = tempIdx-1;
            break;
        end
    end
    sumSignal = sum(BreathAbsSpectrum(1:FFT_size/4+1));
    sumPeak = sum(BreathAbsSpectrum(startInd:endInd));
    if abs(sumSignal - sumPeak) < 0.0001
        confidenceMetric = 0;
    else
        confidenceMetric = sumPeak/sumSignal;
    end
end
end
