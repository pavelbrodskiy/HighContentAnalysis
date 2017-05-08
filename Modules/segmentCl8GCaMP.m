function pxlList = segmentCl8GCaMP(data, dye, settings)
I = dye; %median(data, 3) - settings.background;
I = imclearborder(I);
I = wiener2(I, [6 6]);
bw = im2bw(I, graythresh(I)*0.75);
bw2 = imfill(bw,'holes');
bw3 = imopen(bw2, strel('disk',0));
bw4 = bwareaopen(bw3, 10);
CC = bwconncomp(bw4, 4);
pxlList = CC.PixelIdxList;
end