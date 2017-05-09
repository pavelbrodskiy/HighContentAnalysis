function processConfocalData(labels, settings)

%% Process data
mkdir(settings.thruData);

for i = 1:length(labels)
    if exist([settings.thruData 'LowFrequency_f' num2str(i) '.mat'], 'file')
       continue 
    end
    
    disp(['Processing file: ' labels{i}]);
    
    r = bfGetReader([settings.inData 'LowFrequency' filesep labels{i} '.nd'], 0);
     
    numSeries = r.getSeriesCount();
    
    globalMetadata = r.getGlobalMetadata();
    
    for s = 1:numSeries
        seriesMatOut = [settings.thruData 'LowFrequency_s' num2str(s) '_f' num2str(i) '.mat'];
        
        r.setSeries(s - 1);
        stackSize = [r.getSizeX(), r.getSizeY(), r.getSizeZ(), r.getSizeT(), r.getSizeC()];
        
        hyperstack = zeros(stackSize, 'uint16');
        
        for f = 1:r.getImageCount()
            % Determine image coordinates
            zct = r.getZCTCoords(f - 1);
            qz = zct(1) + 1;
            qc = zct(2) + 1;
            qt = zct(3) + 1;
            
            % Insert frame into hyperstack
            hyperstack(:,:,qz,qt,qc) = bfGetPlane(r, f);
        end
        
        zProj = permute(max(hyperstack, [], 3),[1,2,4,5,3]);
        
        % Extract metadata table for this series
        seriesMetadata = r.getSeriesMetadata();
        javaMethod('merge', 'loci.formats.MetadataTools', ...
            globalMetadata, seriesMetadata, 'Global ');
        metadata.series(s) = seriesMetadata;
        metadata.series1(s) = r.getMetadataStore();
        
        % Save this series to a matlab file
        save(seriesMatOut, ...
            'hyperstack', 'zProj')
    end
    metadata.global = globalMetadata;
    directory = [settings.inData 'LowFrequency' filesep labels{i} '.nd'];
    save([settings.thruData 'LowFrequency_f' num2str(i) '.mat'], ...
            'metadata', 'directory')
    r.close();
end

%% Get metadata
% reader = bfGetReader([PathName FileName]);
% omeMeta = reader.getMetadataStore();



% seriesNumber = reader.getSeriesCount;
%
% for series = 0 %:seriesNumber
%     reader.setSeries(series)
%
%     r = bfopen(reader)
%
% end

% reader.getSizeX()
% reader.getSizeY()
% reader.getSizeZ()
% reader.getSizeC()
% reader.getSizeT()
% reader.getImageCount()










