function stats = extractCellStatistics(I, settings)
basal = median(I);
dFOverF = (I - basal) / basal;
dFOverF = imgaussfilt(dFOverF, 5) - imgaussfilt(dFOverF, 150);

[~,~,w,p] = findpeaks(dFOverF, settings.timestep,'MinPeakWidth',5,'MinPeakProminence',0.05, 'Annotate', 'extents');

delete = w.*p < 3;

w(delete) = [];
p(delete) = [];

nPeaks = length(w);
WHM = mean(w);
amp = mean(p);
dFOverF(dFOverF<0) = 0;
integratedIntensity = trapz(dFOverF); %mean(p.*w/2);

if nPeaks == 0
    WHM = 0;
    amp = 0;
end

stats.basal = basal;
stats.WHM = WHM;
stats.amp = amp;
stats.nPeaks = nPeaks;
stats.integratedIntensity = integratedIntensity;
end

