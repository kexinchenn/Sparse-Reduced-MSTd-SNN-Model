classdef ModelMT < handle
	%% Public Methods
	methods
		function obj = ModelMT(vecDirs, sigmaD, vecSpeeds, sigmaS, ...
				s0, gainIO, betaIO)
			if nargin<1,vecDirs=(0:7)/4*pi;end
			if nargin<2,sigmaD=3;end
			if nargin<3,vecSpeeds=1;end
			if nargin<4,sigmaS=1.16;end
            if nargin<5,s0=0.33;end
			if nargin<6,gainIO=1;end
			if nargin<7,betaIO=1;end
			
			obj.numDirs = numel(vecDirs);
            obj.vecDirs = vecDirs;
			obj.numSpeeds = numel(vecSpeeds);
            obj.vecSpeeds = vecSpeeds;
			obj.sigmaS = sigmaS;
            obj.s0 = s0;
			obj.sigmaD = sigmaD;
			obj.gainIO = gainIO;
			obj.betaIO = betaIO;
			
			obj.fncDir = @(th,th_pref) (exp(obj.sigmaD.*(cos(th-th_pref)-1)));
            obj.fncSpeed = @(sp,sp_pref) (exp(-log((sp+obj.s0) ./ ...
                (sp_pref+obj.s0)).^2 / (2*obj.sigmaS^2)));
			obj.fncIO = @(x) (obj.gainIO * max(x, 0) .^ obj.betaIO);
		end
		
		function respMT = getResponseFromFlow(obj, vx, vy)
			% respMT = MT.getResponseFromFlow(vx, vy) computes MT 
			% population response from a flow field (vx,vy). Local
			% direction of flow will be given by the angle between vx and
			% vy. Local speed will be given by the vector length.
            vx = vx(:); vy = vy(:);
            numPixels = length(vx);
            noFlowPixels = isnan(vx);
            pixelsWithFlow = setdiff(1:numPixels,noFlowPixels);
            
            dir(noFlowPixels) = nan;
			dir(pixelsWithFlow) = atan2(vy(pixelsWithFlow), vx(pixelsWithFlow));
			
            speed(noFlowPixels) = nan;
            speed(setdiff(1:numPixels,noFlowPixels)) = sqrt(vx(pixelsWithFlow).^2 + vy(pixelsWithFlow).^2);
			respMT = obj.calcResponse(dir, speed);
		end
		
		function respMT = getResponseFromDirSpeed(obj, dir, speed)
			% respMT = MT.getResponseFromDirSpeed(dir, speed) computes MT
			% population response from a vector/matrix of direction and
			% (optional) speed values.
			% If no speed values are given, a speed of 1/ms is assumed. If
			% speed is a scalar, that speed is assumed to appear everywhere
			% in the stimulus.
			if nargin<3
				% no speed given: assume 1m/s
				speed=ones(size(dir));
			else
				if isscalar(speed)
					% assume scalar speed, reshape to match dir size
					speed = speed * ones(size(dir));
				end
			end
			respMT = obj.calcResponse(dir(:), speed(:));
		end
		
		function plotResponseFromFlow(obj, vx, vy, mode)
			resp = obj.getResponseFromFlow(vx, vy);
			if ~ischar(mode)
				error('mode must be a string')
			end
			
			switch lower(mode)
				case 'popvec'
					obj.plotResponsePopVec(resp);
				case 'ml'
					obj.plotResponseML(resp);
				case 'heat'
					obj.plotResponseHeat(resp);
				otherwise
					error(['Unrecognized plotting mode "' mode '"'])
			end
		end
		
		function new = clone(obj)
			% make a deep copy
			save('tmp.mat', 'obj');
			Foo = load('tmp.mat');
			new = Foo.obj;
			delete('tmp.mat');
		end
	end
	
	%% Private Methods
	methods (Hidden, Access = private)
		function respMT = calcResponse(obj, dir, speed)
			% Private method to calculate the MT population response based
			% on direction and speed input.
			% DIR and SPEED must be vectors, flattened from 2D flow field.
			assert(isvector(dir))
% 			assert(isvector(speed) && all(speed>=0))
% 			vecSpeeds = 2.^linspace(log(1), log(40), obj.numSpeeds);
% 			vecDirections = (0:obj.numDirs-1)*2*pi/obj.numDirs;
            noFlowPixels = isnan(dir);
            pixelsWithFlow = setdiff(1:length(dir),noFlowPixels);
            
			respMT = [];
			if obj.numSpeeds > 0
				% include speed tuning
				for s=obj.vecSpeeds
					for d=obj.vecDirs
                        resp = zeros(length(dir),1);
                        resp(noFlowPixels) = 0;
                        resp(pixelsWithFlow) = obj.fncIO(...
							obj.fncDir(dir(pixelsWithFlow),d) .* obj.fncSpeed(speed(pixelsWithFlow),s));
                        respMT = [respMT; resp];
					end
				end
			else
				% direction tuning only
				for d=obj.vecDirs
					respMT = [respMT; obj.fncIO(obj.fncDir(dir, d))];
				end
			end
		end
		
		function plotResponsePopVec(obj, resp)
			error('Not yet implemented')
		end
		
		function plotResponseML(obj, resp)
			error('Not yet implemented')
		end
		
		function plotResponseHeat(obj, resp)
			numPx = numel(resp) / obj.numDirs / obj.numSpeeds;
			for s=1:numel(obj.vecSpeeds)
				subplot(numel(obj.vecSpeeds),1,s)
				idxStart = (s-1)*obj.numDirs*numPx + 1;
				idxEnd = s*obj.numDirs*numPx;
                r = reshape(resp(idxStart:idxEnd), sqrt(numPx), ...
                    sqrt(numPx), obj.numDirs);
                r = reshape(permute(r,[2 1 3]), sqrt(numPx), []);
                imagesc(flipud(r), [0 1])
				colorbar
				title(['speed = ' num2str(obj.vecSpeeds(s))])
			end
		end
	end
	
	%% Properties
	properties %(SetAccess = private)
		numDirs;	  % number of directions
        vecDirs;
		numSpeeds;	  % number of speeds
        vecSpeeds;
		
		sigmaS;		  % speed tuning width
        s0;           % speed offset param
		sigmaD;		  % direction tuning width
		
		betaIO;       % exponent of input-output function
		gainIO;       % gain factor of input output function
		
		fncDir;
		fncSpeed;
		fncIO;
	end
end
