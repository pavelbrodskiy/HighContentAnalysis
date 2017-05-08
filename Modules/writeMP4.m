function writeMP4(Img, name, frameRate)
[xMon, yMon, tMon, cMon] = size(Img);
if cMon ~= 1 && cMon ~= 3
    error('Video wrong size')
end

for c = 1:cMon
    
    if nargin < 3
        frameRate = 12;
    end
    
    I = double(Img(:,:,:,c));
    
    minWindow = prctile(I(:), 1);
    maxWindow = prctile(I(:), 99.7);
    
    I = double(I - minWindow);
    I = I / (maxWindow - minWindow);
    I = I * 255;
    I = uint8(I);
    
    rescale = max([xMon, yMon]./[1088, 1920]);
    rescale = max(rescale, 1);
    new = round([xMon, yMon] / rescale);
    
    imgfinal(:,:,c,:) = I;
end

v = VideoWriter([name '.mp4'],'MPEG-4');
v.FrameRate = frameRate;
open(v)

for t = 1:tMon
    compressedFrame = imresize(imgfinal(:,:,:,t),new);
    writeVideo(v, compressedFrame);
end
close(v)


end

