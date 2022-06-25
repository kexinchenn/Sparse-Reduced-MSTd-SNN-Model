%% 
clear variables
close all
%% Initialize optic flow model and MT model
data_dir = 'input_data/';

% MT neurons are tuned to 8 directions and 5 speeds
numSpeeds = 5;
vecSpeeds = 2.^(linspace(log2(0.5), log2(16), numSpeeds));
numDirs = 8;
vecDirs = (0:numDirs-1)/numDirs*2*pi;
nMT = prod(dimPx) * numSpeeds * numDirs;

% initialize models
[flow, MT] = setupModel(vecSpeeds, vecDirs);

%% Set up parameters for generating the training and validation dataset
origFlowTemplate = csvread('16-flow-field-type-param.csv');
numOrigPattern = size(origFlowTemplate, 1);

% Indices of each motion type
rotation = [1,5];      % clockwise and counter clockwise rotation
translation = [9:16];  % translation along eight directions (45 deg intervals in 360 deg)
spiral = [2,4,6,8];    % intermediate spiral motions
expanContract = [3,7]; % expansion and contraction

skipSmall=0;
%% Generate optic flow patterns for the training and validation dataset
numFlows = 1280;
numRep = numFlows/16;
V = zeros(nMT, numFlows);     % MT response 
specs = zeros(numFlows, 12);  % parameters of the optic flow pattern

flow.setDepth('backplane', struct('depth', 1));

tr = 1; % starting ind
for r = 1:numRep
    for i = 1:numOrigPattern
        p = origFlowTemplate(i,:);

        p(1) = normrnd(p(1), 10);
        p(2) = normrnd(p(2), 10);
        
        p(4) = normrnd(p(4), 10);
        p(5) = normrnd(p(5), 10);

        [x,y,z] = sph2cartGu(p(1), p(2), p(3));
        T = [x y z];
        TT = [p(1), p(2), p(3)];

        [x,y,z] = sph2cartGu(p(4), p(5), p(6));
        R = [x y z];
        RR = [p(4), p(5), p(6)];

        [vx,vy] = flow.getFlow(T, R);    

        V(:, tr) = MT.getResponseFromFlow(vx, vy);
        specs(tr, :) = [T R TT RR];
        
        tr = tr + 1;
    end
end
% csvwrite([data_dir 'V-8dir-5speed.csv'], V)
%% Generate dataset for analysis [Gaussian tuning in spiral space]
numFlows = 16;
V = zeros(nMT, numFlows);

flow.setDepth('backplane', struct('depth', 1));
for i = 1:numFlows
    p = origFlowTemplate(i,:);

    [x,y,z] = sph2cartGu(p(1), p(2), p(3));
    T = [x y z];

    [x,y,z] = sph2cartGu(p(4), p(5), p(6));
    R = [x y z];

    [vx,vy] = flow.getFlow(T, R);

    V(:, i) = MT.getResponseFromFlow(vx, vy);
end
% csvwrite([data_dir 'graziano-V-8dir-5speed.csv'], V)
%% Generate dataset for analysis [3D translation and rotation heading selectivity]
flow.setDepth('dotcloud', struct('depthMin',0.05,'depthMax',0.4));

numSamples = 26;

modes = {'trans','rot'};

% elevation: five possible angles
ele = [zeros(1,8) 45*ones(1,8) -45*ones(1,8) -90 90]; % [0 45 -45 -90 90]

% 8 directions for ele=[0 45 -45], then up, then down
azi = [repmat(0:45:359, 1, 3) 0 0];

for m=1:length(modes)
    mode = modes{m};
    V = zeros(prod(flow.dimPx)*MT.numDirs*MT.numSpeeds, numSamples);
    specs = zeros(numSamples, 6);
    switch lower(mode)
        case 'trans'
            % translational stimuli had peak velocity v=0.3m/s (Gu'06)
            rad = ones(size(azi))*0.3;
        case 'rot'
            % rotational stimuli had peak velocity v=20deg/s (Takahashi'07)
            rad = deg2rad(ones(size(azi))*20);
    end
    numStim = numel(ele);
 
    for s=1:numStim
        switch lower(mode)
            case 'trans'
                % translational stimulus has zero rotation
                [x,y,z] = sph2cartGu(azi(s), ele(s), rad(s));
                T = [x y z];
                R = [0 0 0];
            case 'rot'
                % rotational stimulus has zero translation
                [x,y,z] = sph2cartGu(azi(s), ele(s), rad(s));
                T = [0 0 0];
                R = [x y z];
        end

        [vx,vy] = flow.getFlow(T, R);

        V(:,s) = MT.getResponseFromFlow(vx, vy);
        specs(s, 1:3) = T;
        specs(s, 4:6) = R;
    end
%     csvwrite([data_dir, mode, '-V-8dir-5speed-dotcloud.csv'], V)
end
%% Generate dataset for analysis [heading discrimination]
flow.setDepth('dotcloud', struct('depthMin',0.05, 'depthMax',0.75));

azi = 0:45:359;
ele = zeros(1,numel(azi));
rad = ones(1,numel(azi));

numReps = 150;
numSamples = numel(azi) * numReps;

V = zeros(prod(flow.dimPx)*MT.numDirs*MT.numSpeeds, numSamples);
specs = zeros(numSamples, 6);

sampleCount = 1;
for r = 1:numReps
    for s = 1:numel(azi)
        % translational stimulus has zero rotation
        [x,y,z] = sph2cartGu(azi(s), ele(s), rad(s));
        T = [x y z];
        R = [0 0 0];

        [vx,vy] = flow.getFlow(T, R);

        V(:,sampleCount) = MT.getResponseFromFlow(vx, vy);
        specs(sampleCount, 1:3) = T;
        specs(sampleCount, 4:6) = R;
        sampleCount = sampleCount + 1;
    end
end
csvwrite([data_dir,'headingGu2010-V-8dir-5speed.csv'], V)
%% Generate dataset for analysis [Population encoding of heading]

numStim = 10000;                 
                                                                      
depths = 2.^(1:5);                       % distance to backplane, meters

% linear movements
aziLinRange = [60 120]; % 90+-30
eleLinRange = [-30 30]; % 0+-30
linSpeeds = linspace(0.5, 1.5, 5);

% no rotational movements
aziAngRange = [0 0];
eleAngRange = [0 0];
angSpeeds = 0;

V_heading = [];
specs_heading = [];

for d=depths
    % Make a total of `numStim` flow fields
	numStimPerD = ceil(numStim/numel(depths));

	flow.setDepth('backplane', struct('depth', d));
    
    % Generate flow fields. Save both the flow fields (`V`) as well as the
    % heading parameters that generated the flow fields (`specs`). The
    % latter are the six variables of linear and angular velocities: These
    % are the variables we are trying to predict.
	[V_heading(:,end+1:end+numStimPerD), specs_heading(end+1:end+numStimPerD,:)] = ...
		generateInputMatrixRand(flow, MT, numStimPerD, linSpeeds, angSpeeds, ...
		aziLinRange, eleLinRange, aziAngRange, eleAngRange);
end

% Correct for rounding errors
V_heading = V_heading(:,1:numStim);
specs_heading = specs_heading(1:numStim,:);

% csvwrite([data_dir 'heading-V-8dir-5speed.csv'], V_heading)
% csvwrite([data_dir 'heading-specs-8dir-5speed.csv'], specs_heading)