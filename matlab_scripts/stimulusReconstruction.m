clear variables
close all

addpath ('~/CARLsim5/tools/offline_analysis_toolbox/')
addpath ('./util')
%% Analysis setup
indiForAnalysis = '4';

plotWeights = 0;
plotResults = 0;

MSTpixel = 8;
nMST = MSTpixel * MSTpixel;
numSpeeds = 5;
vecSpeeds = 2.^(linspace(log2(0.5), log2(16), numSpeeds));
numDirs = 8;
vecDirs = (0:numDirs-1)/numDirs*2*pi;
nMT = prod(dimPx)*numSpeeds*numDirs;

dimPx=[15 15];

dir = '';
VFile = 'V-8dir-5speed.csv';
weightFile = [dir, 'conn_MT_MST_', indiForAnalysis, '.dat'];
SRFile = [dir 'MST-fr.csv'];
%% Read input
V = csvread(VFile);
%% Read test trial indices
numTest = 160;

trialsAll = dlmread([dir, 'trials.csv'], ',');
if (size(trialsAll,1) > 1)
    trials = trialsAll(str2num(indiForAnalysis)+1, 1:numTest);
else
    trials = trialsAll(1,1:numTest);
end
trials = trials + 1; % matlab index starts from 1
%% Read connection weights
CR = ConnectionReader(weightFile);
[allTimestamps, allWeights] = CR.readWeights();

weightData = reshape(allWeights(end,:), CR.getNumNeuronsPost(), CR.getNumNeuronsPre());
weightData(isnan(weightData)) = 0;
sortedWts = weightData;
%% Read spiking activity
spkDataRaw = load(SRFile);
if size(spkDataRaw,1) > numTest
    spkDataRaw = spkDataRaw(str2num(indiForAnalysis)*numTest+1:(str2num(indiForAnalysis)+1)*numTest,:);
end
%% Plot weight vector
wts = sortedWts';
idxWts = 1:size(wts,2);

skipSmall = false;
% 
if plotWeights
    figure(1);
    [qxWts,qyWts]=plotWts(wts, idxWts, dimPx, vecDirs, ...
	vecSpeeds);
end
%% Plot reconstruction result for each trial
recMTTest = zeros(numTest,nMT);
for trialInd = 1:numTest
    trial = trials(trialInd);  
    spkData = reshape(spkDataRaw(trialInd,:), [MSTpixel MSTpixel]);
    sortedSpkData = spkData;
    %% plot MST activity
    if plotResults
        maxD = max(max(spkData));
        h2 = figure(2);
        c = gray;
        c = flipud(c);
        colormap(c);
        imagesc(sortedSpkData, [0 maxD])
        axis image equal
        title(['rate = [0, ' num2str(maxD) '] Hz'])
        xlabel('nrX')
        ylabel('nrY')
        colorbar
    end
    %% reconstruct MT
    recMT = transpose(sortedWts) * sortedSpkData(:);
    idxMT = 1; % plot one stimulus at a time
    recMTTest(trialInd,:) = recMT;

    if plotResults
        figure(3);
        subplot(1,2,2)
        [qxMT,qyMT] = plotPopulationVector(recMT, idxMT, dimPx, vecDirs, ...
        vecSpeeds, skipSmall);
    end

    %% plot original stimulus
    if plotResults
        subplot(1,2,1);
        [qxOrig, qyOrig] = plotPopulationVector(V, trial, dimPx, vecDirs, ...
            vecSpeeds, skipSmall);
        % display correlation between original and reconstructed input
        disp(['Trial ', num2str(trial), ': ', num2str(corr(recMT, V(:,trial)))])
        pause();
    end
   
end