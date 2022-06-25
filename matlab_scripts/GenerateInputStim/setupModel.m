function [flow,MT] = setupModel(vecSpeedsDeg, vecDirs)
% script to set up FlowRaudies2012 and ModelMT

if nargin<2
    % dir tuning is uniformly distributed, 45deg increments
    numDirs = 8;
    vecDirs = (0:numDirs-1)/numDirs*2*pi;
else
    numDirs = numel(vecDirs);
end

if nargin<1
	% speed tuning is log-Gaussian, log-uniform distribution between 0.5
	% and 16 deg/s (octave-spaced bins)
	numSpeeds = 5;
	vecSpeedsDeg = 2.^(linspace(log2(0.5), log2(16), numSpeeds));
else
	numSpeeds = numel(vecSpeedsDeg);
end


%% set up optic flow geometry

dimPx = [15 15]; % pixel dimensions of the flow field (columns, rows)
dimMeters = [0.01 0.01]; % abs col/row dimensions of image plane
f = 0.01; % focal length
flow = FlowRaudies2012(dimPx, f, dimMeters);
%% set up MT model

% set angular velocity
vecSpeeds = deg2rad(vecSpeedsDeg);

% params from Mineault et al. (2012), Perrone & Stone (1998)
sigmaD = 3; % aka \kappa of von Mises, corresponds to 90deg

% params from Nover et al. (2005)
sigmaS = 1.16;
s0 = deg2rad(0.33);

%... but disable nonlinearity
gainIO = 1;
betaIO = 1;
% betaIO = 0.4;
MT = ModelMT(vecDirs, sigmaD, vecSpeeds, sigmaS, s0, gainIO, betaIO);

end