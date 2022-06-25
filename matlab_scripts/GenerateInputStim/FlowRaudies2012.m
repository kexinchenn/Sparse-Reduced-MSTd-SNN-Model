classdef FlowRaudies2012 < handle
	%% Public Methods
	methods
		function obj = FlowRaudies2012(dimPx, f, dimMeters)
			% flow = FlowRaudies2012(dim, f, imgDim) creates a new optic
			% flow geometry according to Longuet-Higgins & Prazdny (1980)
			% and Raudies et al. (2012).
			%
			% This model uses a pinhole camera with focal length F and with
			% the image plane col \in [-dimPx(1) dimPx(1)] and row \in
			% [-dimPx(2) dimPx(2)].
			%
			% DIM          - pixel dimensions: [numPxCol numPxRow]
			% F            - focal length (meters)
			% DIMMETERS    - image plane dimensions (meters): [mCol mRow]
			if nargin<3,dimMeters=[0.01 0.01];end
			if nargin<2,f=-0.01;end
			if nargin<1,dimPx=[15 15];end
			
			if ~isvector(dimMeters) || ~isnumeric(dimMeters) || ...
					any(dimMeters <= 0) || numel(dimMeters) ~= 2
				error('dimMeters must be [mCol mRow]')
			end
			if ~isvector(dimPx) || ~isnumeric(dimPx) || ...
					any(dimPx <= 0) || numel(dimPx) ~= 2
				error('dimPx must be [dimPxCol dimPxRow]')
			end
			if ~isscalar(f)
				error('Focal length f must be a scalar')
			end
			
			obj.dimPx = dimPx;
			obj.f = f;
			obj.dimMeters = dimMeters;
			
			obj.initDefaultParams()
			obj.initPinholeModel()
		end
		
		function setDepth(obj, mode, options)
			% flow.setDepth(mode, options) sets the depth value of every
			% pixel according to MODE.
			% MODE can be either 'backplane' or 'dotcloud'. OPTIONS is a
			% struct that provides the required arguments for the depth
			% mode.
			%
			% MODE          - 'backplane': Pixels lie on a plane. Their
			%                 depth is given according to the distance and
			%                 slant of the backplane.
			%                 'dotcloud': Pixels are in a dot cloud of
			%                 given depth interval. Their depths are drawn
			%                 randomly (uniformly) from a depth interval.
			% OPTIONS       - struct specifying required parameters.
			%                 For 'backplane':
			%                   - 'depth': depth of the plane (meters)
			%                   - 'alpha': optional slant of plane (deg)
			%                 For 'dotcloud':
			%                   - 'depthMin': min depth (meters)
			%                   - 'depthMax': max depth (meters)
			if ~ischar(mode)
				error('Mode must be a string')
			end
			if ~isstruct(options)
				error('Options must be a struct')
			end
			
			obj.needToCalcZ = true;
			obj.modeZ = mode;
			obj.optionsZ = options;
			obj.calcZ();
		end
		
		function setNoise(obj, mode, options)
			% flow.setNoise(mode, options) sets the noise level of every
			% pixel according to MODE.
			% MODE can be either 'gaussian' or 'outlier'. OPTIONS is a
			% struct that provides the required arguments for the noise
			% mode.
			%
			% MODE          - 'gaussian': Adds Gaussian noise to every
			%                 pixel. Requires OPTIONS fields mu and sigma.
			%               - 'outlier': TODO
			% OPTIONS       - struct specifying required parameters.
			%                 For 'gaussian':
			%                  - 'mu': mean of Gaussian
			%                  - 'sigma': standard deviation
			if ~ischar(mode)
				error('Mode must be a string')
			end
			
			if isempty(mode)
				obj.applyNoise = false;
			else
				if ~isstruct(options)
					error('Options must be a struct')
				end
				
				if strcmpi(mode,'gaussian')
					if ~isfield(options, 'mu')
						error('Options must include a field "mu"')
					end
					if ~isfield(options, 'sigma')
						error('Options must include a field "sigma"')
					end
				elseif strcmpi(mode,'outlier')
					error('Not yet implemented')
				else
					error(['Unknown noise mode "' mode '"'])
				end
				
				obj.applyNoise = true;
				obj.noiseMode = mode;
				obj.noiseOptions = options;
			end
		end
		
		
		
		function [vCol, vRow, col, row] = getFlow(obj, T, R)
			% [vCol,vRow,col,row] = flow.getFlow(T, R) returns the optic
			% flow for translational velocity T and angular velocity R
			% (both in cartesian coordinates); that is, T = [Tx Ty Tz], and
			% R = [Rx Ry Rz].
			%
			% Units of translational velocity are meters per sec (or per
			% frame, to be more exact). Units of angular velocity are
			% radians per sec (or per frame).
			%
			% Returns x and y components of flow (vCol and vRow,
			% respectively), as well as the x and y coordinates (col and
			% row, respectively).
			%
			% Use flow.plotFlow(T, R) to quiver-plot the result instead.
			%
			% Must set up a backplane first via flow.setBackPlane or
			% flow.setRandomDepthCloud.
			%
			% T         - translational velocity (m/s). T = [Tx Ty Tz]
			% R         - angular velocity (m/s). R = [Rx Ry Rz]
			col = obj.col;
			row = obj.row;
			[vCol,vRow] = obj.calcFlow(T, R);
		end
		
		function [vColT,vRowT,vColR,vRowR,col,row] = getFlowComp(obj,T,R)
			% [vColT,vRowT,vColR,vRowR,col,row] = flow.getFlowComp(T,R)
			% returns the translational and rotational flow components in
			% separate matrices.
			% T and R are the same as in flow.getFlow.
			%
			% Units of translational velocity are meters per sec (or per
			% frame, to be more exact). Units of angular velocity are
			% radians per sec (or per frame).
			%
			% Returns x and y components of translational flow (vColT and
			% vRowT, respectively), x and y components of rotational flow
			% (vColR and vRowR, respectively) as well as the x and y
			% coordinates (col and row, respectively).
			%
			% Use flow.plotFlow(T, R) to quiver-plot the result instead.
			%
			% Must set up a backplane first via flow.setBackPlane or
			% flow.setRandomDepthCloud.
			%
			% T         - translational velocity (m/s). T = [Tx Ty Tz]
			% R         - angular velocity (m/s). R = [Rx Ry Rz]
			col = obj.col;
			row = obj.row;
			[~,~,vColT,vRowT,vColR,vRowR] = obj.calcFlow(T, R);
		end
		
		
		function plotFlow(obj, T, R)
			% flow.plotFlow(T, R) plots the optic flow for translational
			% velocity T and angular velocity R (both in Cartesian
			% coordinates); that is, T = [Tx Ty Tz], and R = [Rx Ry Rz].
			%
			% Use flow.plotFlowComp(T, R) to plot the translational and
			% rotational components individually.
			% Use flow.getFlow(T, R) to access [vCol,vRow] at [col,row]
			% instead.
			%
			% Must set up a backplane first via flow.setBackPlane or
			% flow.setRandomDepthCloud.
			%
			% T         - translational velocity (m/s). T = [Tx Ty Tz]
			% R         - angular velocity (m/s). R = [Rx Ry Rz]
			[vCol,vRow] = obj.calcFlow(T, R);
% 			quiver(obj.col, obj.row, vCol, vRow, 'linewidth', 1.5)
            px = ceil(15/10);
            quiver(obj.col(1:px:end,1:px:end), obj.row(1:px:end,1:px:end), ...
                vCol(1:px:end,1:px:end), vRow(1:px:end,1:px:end),...
        		'linewidth', 1.5)
%                     title(['Basis Vector ' num2str(idx(i))])
%             axis off
            set(gca,'XTick',[], 'YTick', [])
            
			minMaxCol = [min(obj.col(:)) max(obj.col(:))];
			minMaxRow = [min(obj.row(:)) max(obj.row(:))];
			heading = obj.f*[T(1) T(2)]/T(3);
			if heading(1) >= minMaxCol(1) && ...
					heading(1) <= minMaxCol(2) && ...
					heading(2) >= minMaxRow(1) && ...
					heading(2) <= minMaxRow(2)
				% display FOE only if within quiver boundaries
				hold on
				plot(obj.f*T(1)/T(3), obj.f*T(2)/T(3), '+', ...
					'MarkerSize', 10, 'LineWidth', 2)
				hold off
			end
			axis equal
			axis(1.2*[min(obj.col(:)) max(obj.col(:)) ...
				min(obj.row(:)) max(obj.row(:))])
% 			xlabel('x (m)')
% 			ylabel('y (m)')
		end
		
		function plotFlowComp(obj, T, R, subplotArgs)
			% flow.plotFlowComp(T, R) plots the optic flow components
			% separately for translational velocity T and angular velocity
			% R (both in cartesian coordinates); that is, T = [Tx Ty Tz],
			% and R = [Rx Ry Rz].
			% Optionally, subplots can be specified for each component.
			% Each row of SUBPLOTARGS should specify the three arguments
			% fed to MATLAB's SUBPLOT function. If the variable contains
			% three rows, the combined flow field is plotted in a third
			% panel. Else only translational and rotational components are
			% plotted.
			%
			% Use flow.plotFlow(T, R) to plot the full flow field instead.
			% Use flow.getFlow(T, R) to access [vCol,vRow] at [col,row]
			% instead.
			%
			% Must set up a backplane first via flow.setBackPlane or
			% flow.setRandomDepthCloud.
			%
			% T         - translational velocity (m/s). T = [Tx Ty Tz]
			% R         - angular velocity (m/s). R = [Rx Ry Rz]
			% SUBPLOTS  - subplot specifier. First row: translational,
			%             second row: rotational, third row (optional):
			%             combined. Default: [1 3 1;1 3 2;1 3 3]
			if nargin<4,subplotArgs=[1 3 1;1 3 2;1 3 3];end
			[vCol,vRow,vColT,vRowT,vColR,vRowR] = obj.calcFlow(T, R);
			
			if ~ismatrix(subplotArgs) || ~isnumeric(subplotArgs) || ...
					~any(size(subplotArgs,1)==[2 3]) || ...
					size(subplotArgs,2) ~= 3
				error('subplotArgs must be 2x3 or 3x3')
			end
			
			% plot translational flow and FOE marker
			subplot(subplotArgs(1,1), subplotArgs(1,2), subplotArgs(1,3))
			quiver(obj.col, obj.row, vColT, vRowT)
			hold on
			plot(obj.f*T(1)/T(3), obj.f*T(2)/T(3), '+', ...
				'MarkerSize', 10, 'LineWidth', 2)
			hold off
			axis equal
			axis(1.2*[min(obj.col(:)) max(obj.col(:)) ...
				min(obj.row(:)) max(obj.row(:))])
			xlabel('x (m)')
			ylabel('y (m)')
			title('Translational Flow')
			
			% plot rotational flow
			subplot(subplotArgs(2,1), subplotArgs(2,2), subplotArgs(2,3))
			quiver(obj.col, obj.row, vColR, vRowR)
			axis equal
			axis(1.2*[min(obj.col(:)) max(obj.col(:)) ...
				min(obj.row(:)) max(obj.row(:))])
			xlabel('x (m)')
			ylabel('y (m)')
			title('Rotational Flow')
			
			% plot combined flow and FOE marker
			if size(subplotArgs,1) == 3
				subplot(subplotArgs(3,1), subplotArgs(3,2), ...
					subplotArgs(3,3))
				quiver(obj.col, obj.row, vCol, vRow)
				hold on
				plot(obj.f*T(1)/T(3), obj.f*T(2)/T(3), '+', ...
					'MarkerSize', 10, 'LineWidth', 2)
				hold off
				axis equal
				axis(1.2*[min(obj.col(:)) max(obj.col(:)) ...
					min(obj.row(:)) max(obj.row(:))])
				xlabel('x (m)')
				ylabel('y (m)')
				title('Combined Flow')
			end
        end
		
		function new = clone(obj)
			% make a deep copy
			save('tmp.mat', 'obj');
			Foo = load('tmp.mat');
			new = Foo.obj;
			delete('tmp.mat');
		end
		
		function print(obj)
			fprintf('FlowRaudies2012 object\n')
			fprintf('----------------------\n')
			fprintf(' - "%s"\n', obj.modeZ)
			if isfield(obj.optionsZ, 'depth')
				fprintf(' - depth: %.2f\n', obj.optionsZ.depth)
			end
			if isfield(obj.optionsZ, 'depthMin')
				fprintf(' - depthMin: %.2f\n', obj.optionsZ.depthMin)
			end
			if isfield(obj.optionsZ, 'depthMax')
				fprintf(' - depthMax: %.2f\n', obj.optionsZ.depthMax)
			end
			if isfield(obj.optionsZ, 'alpha')
				fprintf(' - alpha: %.2f\n', obj.optionsZ.alpha)
			end
			fprintf(' - [%d %d] projected onto [%.2f %.2f]\n', ...
				obj.dimPx(1), obj.dimPx(2), ...
				obj.dimMeters(1), obj.dimMeters(2))
			fprintf(' - f = %.2f\n', obj.f)
        end
        function fovDeg = getFieldOfViewDeg(obj)
			% fovDeg = flow.getFieldOfViewDeg() returns the field of view
			% (FoV) in degrees. FoV depends on camera geometry, namely the
			% focal length and the image plane dimensions. FoV is important
			% to determine angular velocities in px/s.
			imgColRg = range(obj.col(:));
			fovDeg = 2*rad2deg(atan2(imgColRg/2, abs(obj.f)));
		end
	end
	
	
	%% Private Methods
	methods (Hidden, Access = private)
		function initDefaultParams(obj)
			obj.Z = [];
			obj.modeZ = [];
			obj.optionsZ = [];
			obj.needToCalcZ = true;
			obj.setNoise('');
		end
		
		function initPinholeModel(obj)
			[obj.row, obj.col] = meshgrid(linspace(-obj.dimMeters(2), ...
				obj.dimMeters(2), obj.dimPx(2)), ...
				linspace(-obj.dimMeters(1), obj.dimMeters(1), ...
				obj.dimPx(1)));
		end
		
		function [vCol,vRow,vColT,vRowT,vColR,vRowR] = calcFlow(obj,T,R)
			% calculates optic flow field
			if ~isvector(T) || ~isnumeric(T) || numel(T) ~= 3
				error('flow.calcFlow: T must be [Tx Ty Tz]')
			end
			if ~isvector(R) || ~isnumeric(R) || numel(R) ~= 3
				error('flow.calcFlow: R must be [Rx Ry Rz]')
			end
			if isempty(obj.modeZ)
				error('flow.calcFlow: Must call flow.setDepth first')
			end
			
			% translational motion (Tx, Ty, Tz)
			Tx = T(1);
			Ty = T(2);
			Tz = T(3);
			
			% rotational motion in (Rx, Ry, Rz)
			Rx = R(1);
			Ry = R(2);
			Rz = R(3);
			
			vColT = (obj.col*Tz - obj.f*Tx) ./ obj.Z;
			vColR = (obj.col.*obj.row*Rx - (obj.f^2+obj.col.^2)*Ry + ...
				obj.f*obj.row*Rz) / obj.f;
			vCol = vColT + vColR;
			
			vRowT = (obj.row*Tz - obj.f*Ty) ./ obj.Z;
			vRowR = ((obj.f^2+obj.row.^2)*Rx - obj.row.*obj.col*Ry - ...
				obj.f*obj.col*Rz) / obj.f;
			vRow = vRowT + vRowR;
			
			if obj.applyNoise
				if strcmpi(obj.noiseMode,'gaussian')
					noiseCol = normrnd(obj.noiseOptions.mu, ...
						obj.noiseOptions.sigma * ones(size(vCol)));
					vCol = vCol + noiseCol;
					noiseRow = normrnd(obj.noiseOptions.mu, ...
						obj.noiseOptions.sigma * ones(size(vRow)));
					vRow = vRow + noiseRow;
				elseif strcmpi(obj.noiseMode,'outlier')
					error('Not yet implemented')
				else
					error(['Unknown noise mode "' obj.noiseMode '"'])
				end
			end
		end
		
		function calcZ(obj)
			if ~obj.needToCalcZ
				% don't need to recalculate Z
				return
			end
			
			switch lower(obj.modeZ)
				case 'backplane'
					if ~isfield(obj.optionsZ, 'depth')
						error(['In ''backplane'' mode, options ' ...
							'struct must have a field ''depth'''])
					end
					if ~isscalar(obj.optionsZ.depth)
						error('Field ''depth'' must be a scalar')
					end
					if isfield(obj.optionsZ, 'alpha') && ...
							~isscalar(obj.optionsZ.alpha)
						error('Field ''alpha'' must be a scalar')
					end
					
					if isfield(obj.optionsZ, 'depthMax') && ...
							~isscalar(obj.optionsZ.depthMax)
						error('Field ''depthMax'' must be a scalar')
					end
					if ~isfield(obj.optionsZ, 'alpha') || ...
							isnan(obj.optionsZ.alpha)
						% single depth value
						obj.Z = obj.optionsZ.depth * ...
							ones(size(obj.col));
					else
						alphaRad = deg2rad(obj.optionsZ.alpha);
						obj.Z = obj.optionsZ.depth * obj.f ./ ...
							(obj.row.*cos(alphaRad) + obj.f * ...
							sin(alphaRad));
					end
					
					if isfield(obj.optionsZ, 'depthMax')
						% clip depth
						obj.Z = min(obj.optionsZ.depthMax, obj.Z);
					end
					
					% no randomness: need to calc Z only once
					obj.needToCalcZ = false;
				case 'dotcloud'
					if ~isfield(obj.optionsZ, 'depthMin') || ...
							~isfield(obj.optionsZ, 'depthMax')
						error(['In ''dotcloud'' mode, options ' ...
							'struct must have fields ''depthMin'' ' ...
							'and ''depthMax'''])
					end
					if ~isscalar(obj.optionsZ.depthMin) || ...
							~isscalar(obj.optionsZ.depthMax)
						error(['''depthMin'' and ''depthMax'' ' ...
							'must be scalars'])
					end
					
					rangeZ = diff([obj.optionsZ.depthMin ...
						obj.optionsZ.depthMax]);
					obj.Z = rand(size(obj.row)) * rangeZ + ...
						obj.optionsZ.depthMin;
					
					% make sure flag is still set
					obj.needToCalcZ = true;
				otherwise
					error(['Unknown depth mode "' obj.modeZ '"'])
			end
        end
	end
	
	
	%% Properties
	properties %(SetAccess = private)
		dimPx;        % pixel dimensions of flow field: [numPxCol numPxRow]
		dimMeters;    % image plane dimensions (meters): [mCol mRow]
		f;            % focal length (meters)
		
		needToCalcZ;  % flag whether we have to recalculate Z
		Z;            % depth of each pixel (meters)
		modeZ;		  % depth mode: 'backplane' or 'dotcloud'
		optionsZ;	  % optional argument struct for depth mode
		
		applyNoise;   % flag whether to apply noise
		noiseMode;    % noise mode: '' (none), 'gaussian', 'outlier'
		noiseOptions; % optional struct for noise
		
		row;          % meshgrid for y-dimension (rows; meters)
		col;          % meshgrid for x-dimension (columns; meters)
	end
end