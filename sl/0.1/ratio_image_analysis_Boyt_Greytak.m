% For analysis images of colored or fluorescent particles
% 12/22/2025


%% setup

% Define colormaps (all length 256)
redmap=zeros(256,3);
redmap(:,1)=linspace(0,1,256);
greenmap=zeros(256,3);
greenmap(:,2)=linspace(0,1,256);
bluemap=zeros(256,3);
bluemap(:,3)=linspace(0,1,256);
ratiomap=zeros(256,3);
ratiomap(1:128,1)=linspace(0,1,128);
ratiomap(129:256,1)=1;
ratiomap(1:128,2)=1;
ratiomap(129:256,2)=linspace(1,0,128);
ratiomap(1,:)=0;

ratio_min=0.0;
ratio_max=1.0;

ratiolevels=linspace(ratio_min,ratio_max,101);
ratiolevels=(ratiolevels(2:end) + ratiolevels(1:end-1))./2;

% sets the bins used to track pixel or average brightness
levels=0:255;

% set thresholding values to exclude dark and/or saturated regions
 red_squelch=0;
 red_limit=+Inf;
 green_squelch=0;
 green_limit=+Inf; 
 blue_squelch=0;
 blue_limit=+Inf;  
 
 labelstyle='particle'; % 'grid', 'particle', or 'single' 

%% define dataset

% path for input file location [folder of jpg(s)]
indir='D:\jpg Images\';
clear outfile 
% path for output file location
outdir=['D:\Output Text File\'];
outfile='test'; 


% set up a list of files corresponding to a single (or each) experiment

files=dir([indir '*.jpg']);
for(i=1:length(files))
    inputfile(i).rgbfile=files(i).name;
end

% initialize the (manual) region of interest
clear red merged
roix=[];
roiy=[];
merged=imread([indir inputfile(1).rgbfile]);
red=merged(:,:,1);        
myroi=logical(ones(size(red)));

clear red merged files

%% index slices required

% initialize -- set histograms for pixel intensities and label (particle
% average) intensities to zero
    for(theinputfile=1:length(inputfile))
        inputfile(theinputfile).redhist ...
                =zeros(size(levels));
        inputfile(theinputfile).greenhist ...
                =zeros(size(levels));
        inputfile(theinputfile).bluehist ...
                =zeros(size(levels));
        inputfile(theinputfile).brighthist ...
                =zeros(size(levels)); 
        inputfile(theinputfile).validlabels ...
                =[]; 
    end
       
clear  newslice thelocation theslice thefile
%% show an image
theinputfile=1;

titlestr=['Input file ' num2str(theinputfile)];

% aggregate images
clear red green blue merged
merged=imread([indir inputfile(theinputfile).rgbfile]);
red=merged(:,:,1);        
green=merged(:,:,2); 
blue=merged(:,:,3); 

% merged
figure(1);clf;
imshow(merged)
title(titlestr);
text(0.02*max(xlim),0.02*max(ylim),titlestr,'color','w')

% % red slice
figure(2);clf;
imshow(red,[0 max(red(myroi))])
colormap(redmap); colorbar;
title(titlestr);
% 
% % green slice
figure(3);clf;
imshow(green,[0 max(green(myroi))])
colormap(greenmap); colorbar;
title(titlestr)
% 
% % blue slice
figure(4);clf;
    imshow(blue,[0 max(blue(myroi))])
    colormap(bluemap); colorbar;
    title(titlestr)


% display ROI outline
if(exist('roix','var') && length(roix)>2)
    for(i=find(ishandle(1:4)))
        figure(i);
        line(roix,roiy,'color','w','linestyle','-','linewidth',2);
        line(roix,roiy,'color','k','linestyle',':','linewidth',2);
    end
end

%% Reset ROI 
myroi=true(size(red)); % initialize all-on, pick subregion only if we want
roix=[];
roiy=[];

%% Load an ROI

% load myroi matching desired experimental/ instrumental parameters
load 'D:\myroi.mat'

% display ROI outline
if(exist('roix','var') && length(roix)>2)
    for(i=find(ishandle(1:4)))
        figure(i);
        line(roix,roiy,'color','w','linestyle','-','linewidth',2);
        line(roix,roiy,'color','k','linestyle',':','linewidth',2);
    end
end
%% compute statistics for each slice

clear red green blue displaymask

% set values to 1:[# of jpg files in folder] to analyze all files
myinputfiles=1:5;

showimages=0; % set to 1 if you want it to display each time
saveimages=0; % write labeled images out to output directory ... doesn't work yet

if(exist('outfile','var'))
    if(~exist(outdir,'dir'))
        mkdir(outdir);
    end
end
particletable=fopen([outdir outfile '_particles.txt'],'w');

% now gather values for selected slices
for(theinputfile=myinputfiles)
    % print to command line the index of the file we are currently studying    
    fprintf(1,'\n Inputfile %d',theinputfile);
       
    % read data from the current input image of interest
    clear red green blue merged
    merged=imread([indir inputfile(theinputfile).rgbfile]);
    red=merged(:,:,1);        
    green=merged(:,:,2); 
    blue=merged(:,:,3); 

    % Apply correction factors
    if(exist('ff_red','var'))
        red=uint16(double(red).*ff_red); 
    end
    if(exist('ff_green','var'))
        green=uint16(double(green).*ff_green); 
    end

    
    % Thresholding
    % full auto
    mymask = im2bw(merged,graythresh(merged));

    % manual
    mymask = mymask & ...
        (red <= red_limit) & (red >= red_squelch) & ...
        (green <= green_limit) & (green >= green_squelch) & ...
        (blue <= blue_limit) & (blue >= blue_squelch);
    
    % include a manually-chosen region of interest if desired ... if it
    % doesn't exist, we create one that includes all of the image
    if(~exist('myroi','var'))
        myroi=true(size(red));  roix=[]; roiy=[];
    end
    
    if(showimages == 1) % display threshold mask and ROI
        figure(2);
        imshow(mymask & myroi,[]);
        titlestr=['mymask & myroi: Inputfile ' num2str(theinputfile)  ': ' inputfile(theinputfile).rgbfile];
        text(0.02*max(xlim),0.02*max(ylim),titlestr,'color','w')
        % draw ROI box
        if(exist('roix','var') && length(roix)>2)
            line(roix,roiy,'color','w','linestyle','-','linewidth',2);
            line(roix,roiy,'color','k','linestyle',':','linewidth',2);
        end
    end
                      
    % Break image into labeled regions
    mylabel=uint16(zeros(size(red)));
    switch labelstyle
        case 'grid'
            % impose 4x4 grid at image center
            clear x y
            x=floor((0:size(red,2)-1)/(size(red,2)/6));
            y=floor((0:size(red,1)-1)/(size(red,1)/6));
            x(x>4)=0;
            y(y>4)=0;
            [x y]=meshgrid(x,y);
            x(y==0)=0;
            y(x==0)=0;
            y=(y-1)*4;
            y(y<0)=0;
            mylabel=x+y;
            clear x y                 
            mylabel(~mymask)=0;
            mylabel(~myroi)=0;
        case 'particle'
            % break image into labeled thresholded islands (e.g., cells)
            mylabel=bwlabel(mymask);
        otherwise
            % default (trivial case): single tile for entire image
            mylabel(:)=1;
            mylabel(~mymask)=0;
            mylabel(~myroi)=0;
    end
       
    % Display all labeled regions for current slice
    if(showimages==1)
        figure(3);
        titlestr=['All labels: Inputfile ' num2str(theinputfile)  ': ' inputfile(theinputfile).rgbfile];
        imshow(mylabel.*mymask,[]);colormap('colorcube')
        text(0.02*max(xlim),0.02*max(ylim),titlestr,'color','w')
        % draw ROI box
        if(exist('roix','var') && length(roix)>2)
            line(roix,roiy,'color','w','linestyle','-','linewidth',2);
            line(roix,roiy,'color','k','linestyle',':','linewidth',2);
        end
    end
    
    % Compute statistics by pixel (actually just histograms here;
    % averages and variances can be determined later for display)
    redhist=hist(red(mymask & myroi),levels);
    greenhist=hist(green(mymask & myroi),levels);
    bluehist=hist(blue(mymask & myroi),levels);
    brightimg=rgb2gray(merged);
    brighthist=hist(brightimg(mymask & myroi),levels);

    % Compute statistics for all color channels, by label
    validlabels=false(1,max(mylabel(:)));
    for(j=1:max(mylabel(:)))
        if(~(mylabel==j)) % mostly for grid-based labels: skip if no pixels have this label
            continue
        end
        % get particle properties -- characterize entire particle even if
        % part is OUTSIDE of ROI - these will be labeled not "valid" post-process
        % mean intensity
        label(j).red=mean(red(mylabel==j));
        label(j).green=mean(green(mylabel==j));
        label(j).blue=mean(blue(mylabel==j));
        label(j).bright=mean(brightimg(mylabel==j));
        % std dev within each label
        label(j).redstd=std(double(red(mylabel==j)));
        label(j).greenstd=std(double(green(mylabel==j)));
        label(j).bluestd=std(double(blue(mylabel==j)));
        label(j).brightstd=std(double(brightimg(mylabel==j)));
        % get percentile values within each label
        tolerance=[0.05 0.95];
        label(j).redhi=255*max(stretchlim(red(mylabel==j),tolerance));
        label(j).greenhi=255*max(stretchlim(green(mylabel==j),tolerance));
        label(j).bluehi=255*max(stretchlim(blue(mylabel==j),tolerance));
        label(j).brighthi=255*max(stretchlim(brightimg(mylabel==j),tolerance));
        
        labelprops=regionprops((mylabel==j & mymask),'Area','Perimeter');

        label(j).area=sum(mylabel(:)==j); % more reliable: works if disconnected
        label(j).perimeter=labelprops.Perimeter;
        label(j).circ=4*pi()*(label(j).area/(label(j).perimeter^2));   
              
        validlabels(j)=true; % all valid unless disqualifeid below
        % disqualify labels that extend outside ROI borders
        if(find(mylabel==j & ~myroi)) 
            validlabels(j)=false;
        end
        % disqualify labels that touch ROI edge (including edge of
        % image)
        if(find(mylabel==j & bwperim(myroi))) 
            validlabels(j)=false;
        end   
     
        % disqualify particles that we don't like (extremely small
        % particles)
        % add other parameters here to match experimental design
        if(label(j).area < 1000)
            validlabels(j)=false;
        end
    
    end
    
    if(showimages==1 || saveimages==1) % shows all valid particles on screen once process is over
        myvalidlabel=mylabel;
        for(j=1:length(validlabels))
            if(validlabels(j) == 0)
                myvalidlabel(myvalidlabel==j)=0;
            end
        end
    end
    if(showimages==1)
        % make a copy of labeled indexed color image, and delete invalid
        % particles
        figure(4);
        imshow(myvalidlabel,[]);colormap('colorcube')
        titlestr=['Valid labels: Inputfile ' num2str(theinputfile)  ': ' inputfile(theinputfile).rgbfile];
        text(0.02*max(xlim),0.02*max(ylim),titlestr,'color','w')
        % draw ROI box
        if(exist('roix','var') && length(roix)>2)
            line(roix,roiy,'color','w','linestyle','-','linewidth',2);
            line(roix,roiy,'color','k','linestyle',':','linewidth',2);
        end
    end
   
    % Copy values for current image file to inputfile structure
    inputfile(theinputfile).redhist ...
        =redhist;
    inputfile(theinputfile).greenhist ...
        =greenhist;
    inputfile(theinputfile).bluehist ...
        =bluehist;
    inputfile(theinputfile).brighthist ...
        =brighthist;
    inputfile(theinputfile).validlabels ...
        =validlabels;
    inputfile(theinputfile).label ...
        =label;

    % print particle statistics to logfile
    labelfields=fieldnames(label);
    
    % print header row, if desired, if it hasn't been printed already
    if(ftell(particletable)==0)   % tell us if we're at the beginning of the file
        % file name
        fprintf(particletable,'Inputfile ',inputfile(theinputfile).rgbfile);
        % print file number, particle number, valid yes/no
        fprintf(particletable,'FileNo. ParticleNo. Valid ');
        % print label data fields
        for(k=1:length(labelfields)) 
            fprintf(particletable,'%s ',labelfields{k});
        end
        % end of line
        fprintf(particletable,'\r\n');
    end
    
    % generate a line for each particle
    for(j=1:length(validlabels))
        % print input filename, if desired
        fprintf(particletable,'%s ',inputfile(theinputfile).rgbfile);
        % print file number, particle number, valid yes/no
        fprintf(particletable,'%d %d %d ', ...
            theinputfile,j,validlabels(j));
        % print label data fields
        for(k=1:length(labelfields)) 
            fprintf(particletable,'%d ',label(j).(labelfields{k}));
        end
        % end of line
        fprintf(particletable,'\r\n');
    end
    
    % Build up the union of the masks of (valid) labels for all images, for
    % diagnostic purposes
    if( ~exist('displaymask','var') )
        displaymask=zeros(size(red));
    end
    if(exist('displaymask','var') )
        for(j=find(validlabels))
            displaymask=displaymask | (mylabel==j);
        end
    end

end

fclose(particletable);
save([outdir outfile '_data.mat'],'inputfile');

if(exist('displaymask','var'))
    figure(5);clf
    imshow(displaymask)
    titlestr=['All labels, all files'];
    text(0.02*max(xlim),0.02*max(ylim),titlestr,'color','w')
    % draw ROI box
    if(exist('roix','var') && length(roix)>2)
        line(roix,roiy,'color','w','linestyle','-','linewidth',2);
        line(roix,roiy,'color','k','linestyle',':','linewidth',2);
    end
end

clear red green blue redhist greenhist bluehist ratiohist
clear validlabels labelratio labelred labelblue labelgreen 
clear labelredhist labelgreenhist labelratiohist
clear i j
clear mymask particletable
