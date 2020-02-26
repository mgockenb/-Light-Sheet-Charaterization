clearvars
close all

% The only files in the folder should be .bmp files which are images of the light sheet taken at different positions. Each image should
% be separated by the same amount; input that amount as the variable "imgSep".
f = dir;
for i=1:length(f)-2
    filename{i} = f(i+2).name;
end
path = [f(3).folder,'\'];

sheetHorizontal = 0; % 1 if sheet is oriented horizontally. 0 otherwise.

yStart = 0; % The y-axis reading of the first image minus 0.1
imgSep = 0.1; % The z distance in mm between ea img.

% The region containing non-zero values. The ones already below are for images taken on the dates listed.
% sheetRegion=701:2050; % 20NOV, vertical region containing sheet. Sheet is HORIZONTAL.
% sheetRegion=524:2024; % 25NOV, vertical region containing sheet. Sheet is HORIZONTAL.
% sheetRegion=501:1700; % 29JAN, vertical region containing sheet. Sheet is HORIZONTAL.
sheetRegion = 330:730; % 21FEB, horizontal region containing sheet. Sheet is VERTICAL.

L = length(sheetRegion);

N=length(filename);

if sheetHorizontal==1
    A = imread([path,filename{1}]).';
else
    A = imread([path,filename{1}]);
end
sz = size(A);
B = zeros(sz(1),sz(2),N);

%--------------------------------------------------------------------------
% Find the brightest 1.5 mm of the sheet and charaterize that. Discard the
% other 3 mm.
%--------------------------------------------------------------------------
if sheetHorizontal==1
    for i=1:N
        B(:,:,i) = imread([path,filename{i}]).';
    end
else
    for i=1:N
        B(:,:,i) = imread([path,filename{i}]);
    end
end

pixels = length(B(:,1,1));
numPixelsAvg = 25; % The number of neighboring pixels to average
Bavg = zeros(floor(pixels/numPixelsAvg),length(B(1,:,1)),length(B(1,1,:)));
Imax=zeros(length(B(1,1,:)),1);indMax=zeros(length(B(1,1,:)),1);

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
pixel = 2.4; %Second CMOS:2.4;%First CMOS:1.67; % [um] pixel size
x = (1:L)'*pixel; % length scale from pixels to micrometers

sheetWidth = round(1500/pixel); % 900 pixels is approx 1.5mm for the 1st CMOS
startSheet = center - sheetWidth/2; % The vertical height to start considering the sheet

stp = 30; % How big of a step to take between measuring widths of intensity profiles

I = zeros(sheetWidth,L,N); % hold all intensity maps
I_avg = zeros(L,N); % average lateral intensity profile
param = zeros(floor(sheetWidth/stp),14,N); %hold fit parameters
P = zeros(floor(sheetWidth/stp),2,N); % peak values and location of peaks
s = zeros(floor(sheetWidth/stp),3,N);
y = zeros(floor(sheetWidth/stp),L,N);

% fit 4 Gaussians, a line, and a constant to the intensity profiles
G = @(a1,a2,a3,b1,b2,b3,c1,c2,c3,d1,d2,d3,e,f,x) a1*exp(-(x-a2).^2/(2*a3)^2) +...
    b1*exp(-(x-b2).^2/(2*b3)^2) + c1*exp(-(x-c2).^2/(2*c3)^2) +...
    d1*exp(-(x-d2).^2/(2*d3)^2) + e*x + f;
G2=@(k,x) G(k(1),k(2),k(3),k(4),k(5),k(6),k(7),k(8),k(9),k(10),k(11),k(12),k(13),k(14),x);

for i=1:N
    I (:,:,i) = B(startSheet:(startSheet+sheetWidth-1),sheetRegion,i);
    I_avg(:,i) = mean(I(:,:,i));
    
    for j=1:floor(sheetWidth/stp)
        m=j*stp; % Choose the right row to analyze
        [M,c] = max(I(m,:,i));
        a = mean(I(m,:,i));
        c = c*pixel;
        
        %04NOV
        % guess = [M-a,c,17,6,c,72,26,c,8,2,M-250,19,0.006,a];
        %06NOV
         guess = [M-a,c,60,13,c,12,55,c,3,1,290,33,0.002,a];
        optns=optimset('tolX',1e-14, 'Tolfun',1e-30); %'Display','iter',
        param(j,:,i) = lsqcurvefit(G2,guess,x,I(m,:,i)',zeros(14,1),...
            [200,x(end),250,200,x(end),250,200,x(end),250,200,x(end),250,10,100],optns);
    
        y(j,:,i) = G2(param(j,:,i),x)-param(j,13,i)*x-param(j,14,i);
        [P(j,1,i),P(j,2,i)] = max(y(j,:,i));
        P(j,2,i) = P(j,2,i)*pixel;
        
        syms w
        s(j,1,i)=vpasolve(G2(param(j,:,i),w)-param(j,13,i)*w-param(j,14,i)==P(j,1,i)/2,...
            w,[1,P(j,2,i)]);
        s(j,2,i)=vpasolve(G2(param(j,:,i),w)-param(j,13,i)*w-param(j,14,i)==P(j,1,i)/2,...
            w,[P(j,2,i),x(end)]);
        s(j,3,i)=(s(j,2,i)-s(j,1,i)); %in micrometers
    end
    %plot(x,y(j,:,i))
    %hold on
end

s_avg = squeeze(mean(s(:,3,:)));
s_std = squeeze(std(s(:,3,:)));

figure;
scatter((1:N)*imgSep+yStart,s(1,3,:),'.')
hold on
for i=2:floor(sheetWidth/stp)
    scatter((1:N)*imgSep+yStart,s(i,3,:),'.')
end
errorbar((1:N)*imgSep+yStart,s_avg,s_std)
xlabel('Length [mm]')
ylabel('FWHM of sheet [\mum]')
hold off
