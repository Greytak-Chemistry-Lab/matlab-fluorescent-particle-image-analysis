%% Particle image analysis template updated 2025-12-18
% This is a template representative of our approach to particle image
% analysis, exemplified by analysis of true-color micrographs of
% fluorescent beads, for which we are interested in comparison of the
% intensity observed in the different color channels for each bead.
% Greytak Chemistry Laboratory -- Example image from Lavigne Laboratory
% University of South Carolina

% This script allows you to:
% - load a single image from a file (demo is for 8-bit per channel)
% - display auto-scaled, separated color channels with a colorbar
% - select, clear, save, and load a polygon "region of interest" based on 
%   the true color image or any of the channels 
% - select a crop area for display (regions outside the crop area ARE
%   included in calculations: if this is not desired, adjust the ROI, 
%   or crop the image file elsewere and then load it here
% - identify and label particles within the ROI based on a threshold criterion
% - collect a variety of statistics for each of the particles in the ROI
% - produce histograms and scatter plots showing the distribution of
%   statistics, or correlations between statistics, among particles within
%   the ROI

clear
%% setup

% Define colormaps (all length 256)
redmap = zeros(256,3);
redmap(:,1)=linspace(0,1,256);
greenmap = zeros(256,3);
greenmap(:,2)=linspace(0,1,256);
bluemap = zeros(256,3);
bluemap(:,3)=linspace(0,1,256);

%%%%%%% Specify "magic numbers" guiding display and thresholding
% sets the bins used to track pixel or average brightness
levels = 0:255;
% levels=(levels(1:end-1)+levels(2:end))./2; % gets midpoints between level centers

max_level=255; % for brightness scaling of RGB image
% tolerance for stretchlim in display of separated channels, set to 0 for
% full range display 
channel_tol=0.0001; 
% threshold_channel='green';
% threshold_value=50; % 0 for auto, enter a non-zero level for manual

%%%%%%% Image analysis flags
max_scaling_flag=1; % enable brightness scaling of RGB image
sRGB_flag=0; % apply sRGB to linear correction (gamma correction) for levels
       
%% select an image file 
% In this demo, data is loaded from a RGB (color) image file. This cell can
% be modified to load channels from separate files for each channel if
% needed.

% The demo furthermore uses images with 8 bits per channel (max level 255).
% Modifications for 16 bit images and/or sRGB gamma correction could be
% used.

% Choose file with dialog box:
if(~exist('rgbpath'))
    rgbpath='jpg Images';
end
[rgbfile,rgbpath]=uigetfile('*.*','Select an image file',rgbpath);
% rgbpath='..\example_data';  % full or relative path to folder where image is
% rgbfile = 'ImageA.jpg'; % name of RGB image file

titlestr = ['Input file ' rgbfile];

%% (re-)load and display image from the selected RGB file
clear red green blue merged
merged = imread(fullfile(rgbpath,rgbfile));
if(sRGB_flag)
    merged=rgb2lin(merged);
end
red = merged(:,:,1);        
green = merged(:,:,2); 
blue = merged(:,:,3); 

if(~exist('myroi','var'))
    myroi = true(size(red)); % initialize all-on, pick subregion only if we want
end

% merged
figure(1);clf;
if(max_scaling_flag==1)
    imshow(merged(:,:,:).*(max_level/(double(max(merged(myroi))))));
else
    imshow(merged)
end
title(titlestr);
text(0.02*max(xlim),0.02*max(ylim),titlestr,'color','w')

% % red slice
figure(2);clf;
imshow(red,[0 max(stretchlim(red(myroi),channel_tol))]*max_level)
% imshow(red,[0 max(red(myroi))])
colormap(redmap); colorbar;
title(titlestr);
 
% % green slice
figure(3);clf;
imshow(green,[0 max(stretchlim(green(myroi),channel_tol))]*max_level)
% imshow(green,[0 max(green(myroi))])
colormap(greenmap); colorbar;
title(titlestr) 

% % blue slice
figure(4);clf;
imshow(blue,[0 max(stretchlim(blue(myroi),channel_tol))]*max_level)
% imshow(blue,[0 max(blue(myroi))])
colormap(bluemap); colorbar;
title(titlestr)


% cut to crop area
if(exist('cropx','var') && length(cropx)>1)
    for(i=find(ishandle(1:4)))
        figure(i);
        xlim(cropx);
        ylim(cropy); % assume this exists too
    end
end

% display ROI outline
if(exist('roix','var') && length(roix)>2)
    for(i=find(ishandle(1:4)))
        figure(i);
        line(roix,roiy,'color','w','linestyle','-','linewidth',2);
        line(roix,roiy,'color','k','linestyle',':','linewidth',2);
    end
end

%% Pick crop area
figure(1) % brings figure 1 to front
cropx = [];
cropy = [];
croprect=imrect; % NOTE: imrect is now discouraged, Rectangle() is preferred
croppos=ceil(croprect.getPosition); % [ xmin ymin width height]
cropx = [ croppos(1) ( croppos(1) + croppos(3) ) ];
cropy = [ croppos(2) ( croppos(2) + croppos(4) ) ];
croprect.delete;
% apply to current figure
xlim(cropx);
ylim(cropy); % assume this exists too
% apply to all figures
if(exist('cropx','var') && length(cropx)>1)
    for(i=find(ishandle(1:4)))
        figure(i);
        xlim(cropx);
        ylim(cropy); % assume this exists too
    end
end
%% Clear crop area
cropx = [];
cropy = [];
clear croppos
for(i=find(ishandle(1:4)))
    figure(i);
    xlim('auto');
    ylim('auto'); % assume this exists too
end
%% Save crop area
save mycrop cropx cropy
%% Load crop area
load mycrop 
% apply to all figures
if(exist('cropx','var') && length(cropx)>1)
    for(i=find(ishandle(1:4)))
        figure(i);
        xlim(cropx);
        ylim(cropy); % assume this exists too
    end
end
%% Pick ROI
figure(1) % brings figure 1 to front
myroi = true(size(red)); % initialize all-on, pick subregion only if we want
roix = [];
roiy = [];
% uncomment if you want to pick a subregion 
delete(findall(gcf,'type','line'));
[myroi roix roiy]=roipoly;
% display ROI outline
if(length(roix)>2)
    line(roix,roiy,'color','w','linestyle','-','linewidth',2);
    line(roix,roiy,'color','k','linestyle',':','linewidth',2);
end
%% Reset ROI 
myroi = true(size(red)); % initialize all-on, pick subregion only if we want
roix = [];
roiy = [];
%% Save an ROI
save myroi myroi roix roiy
%% Load an ROI
load myroi
% display ROI outline
if(exist('roix','var') && length(roix)>2)
    for(i=find(ishandle(1:4)))
        figure(i);
        line(roix,roiy,'color','w','linestyle','-','linewidth',2);
        line(roix,roiy,'color','k','linestyle',':','linewidth',2);
    end
end

%% find particles
% the settings below are now defined at the top
% threshold_channel='green';
% threshold_value=50; % 0 for auto, enter a non-zero level for manual

switch threshold_channel
    case 'red'
        if(threshold_value==0)
            particle_mask=imbinarize(red);
        else
            particle_mask=red>threshold_value;
        end
    case 'green'
        if(threshold_value==0)
            particle_mask=imbinarize(green);
        else
            particle_mask=green>threshold_value;
        end
    case 'blue'
        if(threshold_value==0)
            particle_mask=imbinarize(blue);
        else
            particle_mask=blue>threshold_value;
        end
end
particle_mask=particle_mask & myroi;

% display mask of all particles, or label image (indexed color image with
% each particle marked by a different color index level)
figure(5)
% imshow(particle_mask)
% imshow(label2rgb(bwlabel(particle_mask)))
imshow(bwlabel(particle_mask),[])
colormap('colorcube')

% cut to crop area
if(exist('cropx','var') && length(cropx)>1)
    for(i=find(ishandle(1:5)))
        figure(i);
        xlim(cropx);
        ylim(cropy); % assume this exists too
    end
end

% display ROI outline
if(exist('roix','var') && length(roix)>2)
    for(i=find(ishandle(1:5)))
        figure(i);
        line(roix,roiy,'color','w','linestyle','-','linewidth',2);
        line(roix,roiy,'color','k','linestyle',':','linewidth',2);
    end
end

%% get particle data
particle_label=bwlabel(particle_mask);
particles=struct('area',0,'perimeter',0,'red_total',0,'green_total',0,'blue_total',0);
for(i=1:max(particle_label(:)))
    particles(i).area=sum(particle_label(:)==i);
    particles(i).perimeter=bwperim(particle_label==i);
    particles(i).red_total=sum(red(particle_label==i));
    particles(i).green_total=sum(green(particle_label==i));
    particles(i).blue_total=sum(blue(particle_label==i));
end

%% How to analyze your results now that you have separately labeled each particle
% plot any property versus any other to look for correlations
figure(6)
plot([particles.red_total],[particles.green_total],'ks')
xlim([0 max(xlim)])
ylim([0 max(ylim)])
xlabel('Red total')
ylabel('Green total')

% histogram of particle properties learn more about Matlab's 'histogram'
% and related functions from the documentation. Enclosing [particles.area]
% in [] makes it pass a vector of the values of the area field in every
% element of the 'particles' structure array to the histogram function.
figure(7) 
histogram([particles.area],[0:1000:11000])
xlabel('Particle area')
ylabel('Number found')
