clearvars
close all

% No user defined functions needed unless you don't want to use the 'dir' command.
% Visualize the light sheet along the z-axis. The camera, of course, takes
% images in the xy-plane, so in order to see the third dimension, must
% interpolate the intensity between the pictures. 
% Must be in the folder with only the .bmp pics of the LS.

pixel = 1.67; % [um] pitch of camera pixels
useDir = 1;   % Use the 'dir' function when you are in the correct directory
sheetHorizontal = 1; % 1 if sheet is oriented horizontally. 0 otherwise.

% the horizontal region containing non-zero values.
% sheetRegion = 1900:3500; %04NOV_1pixel
% sheetRegion = 393:1812; %06NOV
% sheetRegion = 1000:1900; %08NOV vertical region containing sheet
 sheetRegion = 701:2050; %20NOV vertical region w/ sheet. Sheet is HORIZONTAL.
L = length(sheetRegion);

if useDir==1
    f = dir;
    for i=1:length(f)-2
        filename(i,:) = f(i+2).name;
    end
    path = f(1).folder;
else
    [path,filename] = LoadFilenames04NOV_2pixel;
end

N=length(filename);

if sheetHorizontal==1
    A = imread([path,filename(1,:)]).';
else
    A = imread([path,filename(1,:)]);
end
sz = size(A);
B = zeros(sz(1),sz(2),N);

%--------------------------------------------------------------------------
% Find the brightest 1.5 mm of the sheet and charaterize that. Discard the
% other 3 mm.
%--------------------------------------------------------------------------
if sheetHorizontal==1
    for i=1:N
        B(:,:,i) = imread([path,filename(i,:)]).';
    end
else
    for i=1:N
        B(:,:,i) = imread([path,filename(i,:)]);
    end
end

pixels = length(B(:,1,1));
numPixelsAvg = 25; % The number of neighboring pixels to average
Bavg = zeros(floor(pixels/numPixelsAvg),length(B(1,:,1)),length(B(1,1,:)));
Imax=zeros(length(B(1,1,:)),1); indMax=zeros(length(B(1,1,:)),1);

[h,index]=max(mean(B(:,:,:)));
h=squeeze(h); index=squeeze(index);
for i = 1:floor(pixels/numPixelsAvg)
    begin=numPixelsAvg*(i-1)+1;
    stop = numPixelsAvg*i;
    Bavg(i,:,:) = mean(B(begin:stop,:,:));
end
for i=1:N
    [Imax(i),indMax(i)] = max(Bavg(:,index(i),i));
end
[ImaxSorted,iSorted]=sort(Imax,'ascend'); % Sort by brightest intensity, which corrseponds to most in focus
indMaxSorted = indMax(iSorted);
center = round(mean(indMaxSorted(1:5))*numPixelsAvg); % Find the center of the bright region 
                                                      % of the sheet by averaging the 5 most in focus images.

s = round(100/pixel);  % The number of camera pixels separating each image in the folder. Each image is taken 100um apart.

I = zeros(N*(s-1) + 1,length(sheetRegion)+1000);                      

for i=1:N
    I(60*(i-1)+1,501:(end-500)) = B(center,sheetRegion,i);
end

[r,c] = find(I);

[maxT,colT]= max(I(r(1),:));
[maxB,colB]= max(I(r(end),:));

slopeCols = (colB - colT)/(r(end)-r(1));


h=length(r); w=length(c);
delta=zeros(h-1,1);
slope=zeros(h-1,1);

for j=1:(h-1)
    if r(j)==r(1)
       col_o = c(j);
    else
        col_o = col_f;
    end
    
    if r((j+1))>r(j) && col_o+(colB-colT) <= (length(sheetRegion)+1000) 
        col_f = col_o + round(slopeCols*(r(j+1)-r(j)));
        delta(j) = I(r((j+1)),col_f) - I(r(j),col_o);
        slope(j) = delta(j)/(r((j+1))-r(j));
        
        for k = 1:(r((j+1))-r((j))-1)
            I(k+r(j), round(k*slopeCols) + col_o) = I(r(j),col_o) + slope(j)*k;
        end
    end
end
%
% 06NOV_1pixel
%{
I2=I';
I2(552:610, 315:321) = 220;
position = [80,550];
I2(800:805,340:939) = 220;
position2 = [560,815];
%}
%04NOV
%{
I2=imrotate(I,90);
I2(462:520, 315:321) = 220;
position = [80,460];
I2(705:710,390:989) = 220;
position2 = [610,705];
%}
%20NOV
%
I2=I';
I2(1162:1220, 515:521) = 220;
position = [220,1160];
I2(1105:1110,590:1189) = 220;
position2 = [810,1105];
%}
text = '100 microns';
text2 = '1 mm';
I3= insertText(I2,position,text,'BoxOpacity',0,'TextColor',[220 220 220],'FontSize',42);
I4= insertText(I3,position2,text2,'BoxOpacity',0,'TextColor',[220 220 220],'FontSize',42);
imshow(mat2gray(I4))
set(gca,'DataAspectRatio',[1 1 1])
xlabel('z-direction (optical axis)')
ylabel('y-direction')
%}
