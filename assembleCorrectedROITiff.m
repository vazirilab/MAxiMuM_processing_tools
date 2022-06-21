function [imageData,frameRate,pixelResolution] = assembleCorrectedROITiff(filename)

%% Load in data and determine size of dimensions

[roiData, roiGroup, header, ~] = scanimage.util.getMroiDataFromTiff(filename); % load in data throuhg scanimage utility
numROIs = numel(roiData); % number of ROIs (ASSUMES THEY ARE ORDERED LEFT TO RIGHT)
totalFrame = length(roiData{1}.imageData{1}); % total number of frames in data set
totalChannel = length(roiData{1}.imageData); % number of channels
frameRate = header.SI.hRoiManager.scanVolumeRate;
sizeXY = roiGroup.rois(1,1).scanfields.sizeXY;
FOV = 157.5.*sizeXY;
numPX = roiGroup.rois(1,1).scanfields.pixelResolutionXY;
pixelResolution = mean(FOV./numPX);

% find size of final assembled FOV
strip=1; channel=1; frame=1;
stripTemp = roiData{1,strip}.imageData{1,channel}{1,frame}{1,1};
roiSize = size(stripTemp);
corr = returnScanOffset2(stripTemp,1); % find offset correction
stripTemp = fixScanPhase(stripTemp,corr,1); % fix scan phase
val = round(size(stripTemp,1)*0.03); % trim excess
stripTemp = stripTemp(val:end,7:138,:,:);

fullPXx = size(stripTemp,1);            % final FOV x size (px)
fullPXy = size(stripTemp,2) * numROIs;  % final FOV y size (px)

%% Assemble frames
imageData = zeros(fullPXx, fullPXy, totalChannel,totalFrame, 'single');
for channel = 1:totalChannel
    disp(['Assembling channel ' num2str(channel) ' of ' num2str(totalChannel) '...'])
    frameTemp = zeros(fullPXx, fullPXy / numROIs, 1, totalFrame, 'single');
    for strip = 1:numROIs

        % Generate the time series of each ROI in the data
        
        stripTemp = zeros(roiSize(1), roiSize(2), 1, totalFrame);
        for frame = 1:totalFrame
            stripTemp(:,:,1,frame) = single(roiData{1,strip}.imageData{1,channel}{1,frame}{1,1});
        end
        corr = returnScanOffset2(stripTemp,1); % find offset correction
        stripTemp = fixScanPhase(stripTemp,corr,1); % fix scan phase
        val = round(size(stripTemp,1)*0.03); % trim excess
        stripTemp = stripTemp(val:end,7:138,:,:);
    
        yy = (strip-1)*fullPXy / numROIs+1:strip*fullPXy / numROIs;
        frameTemp(:,yy,1,:) = stripTemp;

    end

    imageData(:,:,channel,:) = frameTemp;
end