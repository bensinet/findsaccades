clear all; close all;
% find saccades in TSLO extraction
%in order to run this code you will need extractsaccadesnewv1,
%mediansmooth, and ginputc in your matlab path
%this code calls extractsaccadesnewv1 where you will use a gui to select
%manually confirm the selection of saccades. To add a saccade simply left
%click the point where there is a saccade, to remove a saccade right click.
%To clear the everything press 'e' then select the saccades you want.
%When you are done hit enter
% This program assumes that a dataset has been created by GKR_SBS_AOSLO 
% software from an SLO movie
% portions (c)2006 SBStevenson@uh.edu
% edited by Ethan Bensinger

[matname, matpath] = uigetfile('*.mat','Find a MAT file with frameshifts data in it','MultiSelect','On');
load([matpath matname]); 
%load in the matfiles generated from stablized from raw

    %%% this next section of code is from Girish
    degrees_per_pixel = [1/100,1/100]; %pixels per degree x,y  %%% this value must be checked each time!
    %%% he converts to degrees
    arcminutes_per_pixel = [60/100,60/100];
    %%% he converts to arcmin
    
    xy = frameshifts_strips_spline.*repmat(degrees_per_pixel,length(frameshifts_strips_spline),1);
    xy = xy - ones(length(xy),1) * mean(xy);
    %center it around 0
    figure(541);plot(timeaxis_secs, xy, ':');
    %create the first plot with the raw eye trace
% The TSLO and AOSLO have an artifact that shows up at the frame rate of 30 Hz.
% The artifact occurs because of distortion in the reference frame, and also because of torsion
% The reference distorsion is removed by calculating the frame averaged eye motion, but the torsion
% artifact changes over time. For this I have used a running average frame, using a median smooth
% We should be able to extract torsion and remove it that way, but torsion signals are much smaller
% than x and y in the SLO
    xwrapped = reshape(xy(:,1), 16, length(xy)/16); % assumes 16 strips per frame 480 hz sampling 30 hz frame rate
    ywrapped = reshape(xy(:,2), 16, length(xy)/16);
    xwrapped = xwrapped - ones(16,1) * xwrapped(1,:);
    ywrapped = ywrapped - ones(16,1) * ywrapped(1,:);
    xwrapped_smoothed = mediansmooth(xwrapped', 20)'; % transposed to smooth on the 2nd dimension, second argument is the range to smooth. Since we have wrapped it now we have one row per frame so 15 is 1/2 second
    ywrapped_smoothed = mediansmooth(ywrapped', 20)';
    xyfixed = xy - [xwrapped_smoothed(:), ywrapped_smoothed(:)];
    figure(541);hold on; plot(timeaxis_secs, xyfixed);hold off
    %create the plot with the filered traces
% using a median smooth we have removed the 30 Hz so now we extract the saccades
sacstatsmoothed = extractsaccades326(xyfixed, timeaxis_secs, 10, 2);
%do another basic smooth to try and remove some of the jittering artifact
%we found in a few videos
         savefig(figure(541),strcat(matpath,matname,'xypositions.fig'));
         savefig(figure(301),strcat(matpath,matname,'xysacpositionsmoothed.fig'));
         savefig(figure(302),strcat(matpath,matname,'velocity.fig'));
         savefig(figure(303),strcat(matpath,matname,'direction&amplitude.fig'));
         savefig(figure(304),strcat(matpath,matname,'xypositionsmoothedclean.fig'));
         %save all the figures made within extractsaccadesnewv1
%      end
sacdirection=sacstatsmoothed(9);
sacamp=sacstatsmoothed(10);
sacvel=sacstatsmoothed(12);
save(strcat(matpath,matname,'sacdirectionampsandvel.mat'),'sacamp','sacdirection','sacvel')
%get the data on every saccade into a mat file
sacnum=sacstatsmoothed{1};
avesacvel=sacstatsmoothed{2};
maxsacvel=sacstatsmoothed{3};
avesacamps=sacstatsmoothed{4};
sacfreq=sacstatsmoothed{5};
driftamp=sacstatsmoothed{6};
driftdirect=sacstatsmoothed{7};
driftvel=sacstatsmoothed{8};
T=table(sacnum,avesacvel,maxsacvel,avesacamps,sacfreq,driftamp,driftdirect,driftvel);
writetable(T,strcat(matpath,matname,'sacstats.csv'),'Delimiter',',','QuoteStrings',true)
%write all the averaged variables of interest into an excel file labelled
%sacstats.csv into your matpath