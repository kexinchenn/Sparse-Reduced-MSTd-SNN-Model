function [V,specs] = generateInputMatrixRand(flow, MT, numSamples, ...
    linSpeeds, angSpeeds, aziLinRange, eleLinRange, ...
	aziAngRange, eleAngRange, shuffle)
if nargin<10,shuffle=false;end
if nargin<9 || isempty(eleAngRange),eleAngRange=[-90 90];end
if nargin<8 || isempty(aziAngRange),aziAngRange=[0 359];end
if nargin<7 || isempty(eleLinRange),eleLinRange=[-90 90];end
if nargin<6 || isempty(aziLinRange),aziLinRange=[0 359];end
if nargin<5 || isempty(angSpeeds),angSpeeds=[0 5 10];end
if nargin<4 || isempty(linSpeeds),linSpeeds=[0.25 0.5 1];end
if nargin<3,numSamples=1000;end

if ~isscalar(numSamples)
	error('numSamples must be a scalar')
end
if ~isnumeric(linSpeeds) || ~isvector(linSpeeds)
	error('linSpeeds must be a numeric vector')
end
if ~isnumeric(angSpeeds) || ~isvector(angSpeeds)
	error('angSpeeds must be a numeric vector')
end
if ~islogical(shuffle)
	error('shuffle must be true/false')
end
%% Initialize

% pick random elements from the list of specified speeds
ls = linSpeeds(round(rand(1, numSamples) * (numel(linSpeeds)-1)) + 1);
as = angSpeeds(round(rand(1, numSamples) * (numel(angSpeeds)-1)) + 1);

% azimuth in [0,2*pi)
linAzi = mod(rand(1, numSamples) * diff(aziLinRange) + aziLinRange(1), 360);
angAzi = mod(rand(1, numSamples) * diff(aziAngRange) + aziAngRange(1), 360);

% elevation in [-pi/2, pi/2]
linEle = rand(1, numSamples) * diff(eleLinRange) + eleLinRange(1);
angEle = rand(1, numSamples) * diff(eleAngRange) + eleAngRange(1);
%% Generate

if MT.numSpeeds > 0
	V = zeros(prod(flow.dimPx)*MT.numDirs*MT.numSpeeds, numSamples);
else
	V = zeros(prod(flow.dimPx)*MT.numDirs, numSamples);
end
specs = zeros(numSamples, 14);
for s=1:numSamples
    if ls(s) > 0
		[x,y,z] = sph2cartGu(linAzi(s), linEle(s), ls(s));
        T = [x y z];
        TT = [linAzi(s) linEle(s) ls(s)];
    else
        T = [0 0 0];
        TT = [0 0 0];
    end

    if as(s) > 0
    	[x,y,z] = sph2cartGu(angAzi(s), angEle(s), as(s));
        R = [x y z];
        RR = [angAzi(s) angEle(s) as(s)];
    else
        R = [0 0 0];
        RR = [0 0 0];
    end
    
    heading = flow.f*[T(1) T(2)]/T(3);
	
	[vx,vy] = flow.getFlow(T, R);
	V(:,s) = MT.getResponseFromFlow(vx, vy);
	specs(s,:) = [T R TT RR heading];
end


%% Shuffle

if shuffle
	p = randperm(size(V,2));
	V = V(:, p);
	specs = specs(p, :);
end

end