function [qx,qy]=plotWts(wts, idx, dimPx, vecDirs, ...
	vecSpeeds, skipSmall)
    if nargin<2 || isempty(idx),idx=1:size(wts,2);end
    if nargin<3 || isempty(dimPx),dimPx=[15 15];end
    if nargin<4 || isempty(vecDirs),vecDirs=(0:7)*pi/4;end
    if nargin<5 || isempty(vecSpeeds),vecSpeeds=2.^(linspace(log2(0.5), log2(32), 5));end
    if nargin<6,skipSmall=false;end

    numDirs = numel(vecDirs);
    numSpeeds = numel(vecSpeeds);
    px = ceil(min(dimPx)/10);

    assert(size(wts,1) == prod(dimPx)*numel(vecDirs)*numel(vecSpeeds))
    %% Compute population vector

    qx = zeros(dimPx(1), dimPx(2), numel(idx));
    qy = zeros(dimPx(1), dimPx(2), numel(idx));

    for i=1:numel(idx)
        % process each pixel, estimate direction and speed
        W = reshape(wts(:,idx(i)),dimPx(1),dimPx(2),numDirs,numSpeeds);
        len = zeros(dimPx);
        for r=1:dimPx(1)
            for c=1:dimPx(2)
                dirNumX = 0;
                dirNumY = 0;
                dirDenom = 0;
                speedNum = 0;
                speedDenom = 0;
                for s=1:numSpeeds
                    ss = vecSpeeds(s);
                    for d=1:numDirs
                        % population vector for direction
                        dd = vecDirs(d);
                        dirNumX = dirNumX + W(r,c,d,s)*cos(dd);
                        dirNumY = dirNumY + W(r,c,d,s)*sin(dd);
                        dirDenom = dirDenom + W(r,c,d,s);

                        % population vector for speed
                        speedNum = speedNum + W(r,c,d,s)*ss;
                        speedDenom = speedDenom + W(r,c,d,s);
                    end
                end
                len(r,c) = dirNumX^2 + dirNumY^2;

                % make sure direction vector is normalized
                % this will determine the orientation of the flow vector
                inferredDir = [dirNumX dirNumY]/dirDenom;
                inferredDir = inferredDir / norm(inferredDir);

                % speed will determine the magnitude of the flow vector
                inferredSpeed = speedNum / speedDenom;

                qx(r,c,i) = inferredDir(1)*inferredSpeed;
                qy(r,c,i) = inferredDir(2)*inferredSpeed;
            end
        end

        if skipSmall
            qqx = qx(:,:,i);
            qqy = qy(:,:,i);
            qqx(len < 0.2*max(len(:))) = 0;
            qqy(len < 0.2*max(len(:))) = 0;
            qx(:,:,i) = qqx;
            qy(:,:,i) = qqy;
        end
    end


    %% Plot population vector

    nr=floor(sqrt(numel(idx)));
    nc=ceil(numel(idx)/nr);
    [yMesh,xMesh] = meshgrid(1:dimPx(2),1:dimPx(1));
    i = 1;

    for col=1:nc
        for row=1:nr
            if numel(idx) > 1
                subplot('Position',[(col-1)*1/nc (nr-row)*1/nr 1/nc 1/nr])
                hold off
            end

            quiver(xMesh(1:px:end,1:px:end), yMesh(1:px:end,1:px:end), ...
                qx(1:px:end,1:px:end,i), qy(1:px:end,1:px:end,i),...
                'linewidth', 1.5)
            axis equal
            set(gca,'XTick',[], 'YTick', [])
            hold off
            i = i + 1;
        end
    end

end