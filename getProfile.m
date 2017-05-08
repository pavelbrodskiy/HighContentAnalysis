imshow(a{1}(:,:,1),[])
[y,x] = ginput(1);
profile = squeeze(mean(mean(a{1}(uint8(x)-4:uint8(x)+4,uint8(y)-4:uint8(y)+4,:),2),1));
t = linspace(0, 400*3/60, 401);
plot(t, profile);
xlabel('Time (min)')
ylabel('Fluorescent Intensity')