function stats = extractWellStatistics(data, dye, settings)
    [X, Y, T] = size(data);

    %% Segment cells
    pxlList = segmentCl8GCaMP(data, dye, settings);
    cellNumber = length(pxlList);
    
    %% Extract stats from each cell
    for i = cellNumber:-1:1
        pxlListtmp = pxlList{i};
        profile = nan(1, T);
        for t = T:-1:1
            I_Mat = data(:,:,t);
            I_Array = I_Mat(:);
            
            profile(t) = mean(I_Array(pxlListtmp));
        end
        cellStats(i) = extractCellStatistics(profile, settings);
    end
    
    for i = cellNumber:-1:1
        cellStats(i).pxlSize = length(pxlList{i});
    end
    
    %% Package well-level stats
    if cellNumber > 1
        stats.meanAmp       = mean([cellStats.amp]);
        stats.meanBasal     = mean([cellStats.basal]);
        stats.meanWHM       = mean([cellStats.WHM]) * settings.timestep;
        stats.meanfreq      = mean([cellStats.nPeaks]) / (T * settings.timestep);
        stats.meanInteg     = mean([cellStats.integratedIntensity]) / settings.timestep;
    
        stats.cellStats     = cellStats;
    else
        stats.meanAmp       = NaN;
        stats.meanBasal     = NaN;
        stats.meanWHM       = NaN;
        stats.meanfreq      = NaN;        
        stats.meanInteg     = NaN;
    end
    
    stats.cellNumber        = cellNumber;
end

