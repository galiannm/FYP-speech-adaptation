function [roiIdx, scalpXY] = selectScalpROI(chanlocs, nChannels, centre, roiRadius)
%SELECTSCALPROI Select a circle in the EEGLAB 2-D scalp projection.
% Nose up; head outline radius 0.5; negative y is posterior.
% Uses the same projection as FYP_forwardTRF_singleFeature.m.
validateattributes(centre,{'numeric'},{'real','finite','size',[1,2]});
validateattributes(roiRadius,{'numeric'},{'real','finite','scalar','positive'});
if ~isstruct(chanlocs) || numel(chanlocs) ~= nChannels || ...
        ~all(isfield(chanlocs,{'theta','radius'}))
    error('selectScalpROI:MissingCoordinates', ...
        'Spatial ROI requires theta/radius coordinates for every model channel.');
end

scalpXY = nan(nChannels,2);
for ch = 1:nChannels
    theta = chanlocs(ch).theta;
    radius = chanlocs(ch).radius;
    if isnumeric(theta) && isreal(theta) && isscalar(theta) && isfinite(theta) && ...
            isnumeric(radius) && isreal(radius) && isscalar(radius) && ...
            isfinite(radius) && radius >= 0
        scalpXY(ch,:) = [-radius*sind(theta),radius*cosd(theta)];
    end
end
validCoordinates = all(isfinite(scalpXY),2);
distance = sqrt(sum((scalpXY-centre).^2,2));
roiIdx = find(validCoordinates & distance <= roiRadius);
if isempty(roiIdx)
    error('selectScalpROI:EmptyROI', ...
        'No electrodes fall within the spatial ROI. Check coordinates, centre and radius.');
end
if any(~validCoordinates)
    warning('selectScalpROI:InvalidCoordinates', ...
        'Excluded %d channels with missing scalp coordinates.',sum(~validCoordinates));
end
end
