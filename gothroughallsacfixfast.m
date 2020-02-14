%this code will run through all the files that have matvelocity figures and
%also have corresponding mat files from the stablized from raw output
%pulling the figures may return an error the first time you use it simply
%run it again and it should work
% Contact Ethan Bensinger (bensinet@gmail.com) with questions
clearvars;
%select all the matvelocity figures
[fnames,pname] = uigetfile('*matvelocity.fig','Open AVI files','MultiSelect','On');
%set up a new figure to generate some new plots
amps=figure;
%lets define all the thresholds and variables we will use below 
degrees_per_pixel = [1/100,1/100]; %pixels per degree x,y  % this value must be checked each time!
mediansmoothwindow=20; %argument is the range to the mediansmooth function. Since we have wrapped the input now we have one row per frame so 15 is 1/2 second
gausssmoothsd=3; %arguement is the sd to the smoothgaussr function
sgolaywindow=15;
sgolaypoly=4;
ampthresh=degrees_per_pixel; %amplitude threshold in degrees (this value should scale with the degrees_per_pixel since that effects the noise of the registeration). deg/pixel/2 is slightly above Sheehy et al. noise levels
windowbefore=20; %number of samples to look before the sac mark
lenswobblelen=8; %number of samples that lens wobble takes
windowafter=windowbefore+lenswobblelen; %number of samples to look after the sac mark
tooclosethresh=4; %checks if saccades are too close using this threshold

startthresh=3; %velocity threshold for the start of the saccade
endthresh=4; %velocity threshold for the end of the saccade (must be higher to deal with high velocity lens wobble
defaultstart=15; %if the velocity of the saccade never reaches the start threshold how many samples back should we set the start
defaultend=defaultstart+lenswobblelen;
agecollumn=6; %just defining collumns to get stats we need
mscollumn=2; %just defining collumns to get stats we need

%this just starts the counter for how many saccades didn't pass
tooclosestart=0;
toocloseend=0;
toosmall=0;
tooslow=0;

%edit the file below that has the corresponding excel sheet. This excel
%sheet needs to have the first collumn be the corresponding ID numbers and
%make sure that it has 2 sheets to it since it will write the outputs to
%the second sheet.
filenameforexcel='C:\Users\Lenovo\Downloads\PHI_TSLORPT.xlsx';
[num,text,raw]=xlsread(filenameforexcel,1);
%idall grabs all the id's that are in the sheet. If we have files from id's
%not in the sheet these are thrown out.
idall=num(:,1);
msall=num(:,mscollumn);
ageall=num(:,agecollumn);
% here we start the loop going through all the data
for i=1:size(fnames,2)
    clear sacstarts sacends refractoryind
    close
    %tracenum(i) grabs the corresponding tracenumber V0__ from the filename
    tracenum(i)=str2double(fnames{i}(8))*1000+str2double(fnames{i}(13));
    idnumall{i}=fnames{i}(1:6);
    %idnum(i) grabs the corresponding idnumber from the filename
    idnum(i)=str2double(fnames{i}(1:5));
    
    %grab the eye from the file name and check if it is correct
    eyel=fnames{i}(6);
    if eyel=='L' || eyel=='l'
        eye(i)=1;
    elseif eyel=='R' || eyel=='r'
        eye(i)=2;
    else
        eye(i)=0;
        incorrecteye{i}=fnames{i};
    end
    ind=find(idall==idnum(i));
    %we check if the id filename matches the id in the sheet
    if isempty(ind)
        idnotinset{i}=fnames{i};
        continue
    end
    
    %grab the ms vs control and patient age from excel sheet
    ms(i)=msall(ind);
    age(i)=ageall(ind);
%     idnum(i)=idnumb(i);
%     eye(i)=eyeb(i);
%     tracenum(i)=tracenumb(i);
    %open the velocity figure grab all of it's objects, the only 
    %one we are interested in is xdata which contains where we marked the 
    %index of the microsaccades. It is a cell array so call it with xdata{1,1}
    f=openfig([pname, fnames{i}],'invisible');
    dataObjs=findobj(f,'Type','line');
    xdata = get(dataObjs, 'XData');
       
    % get the name of the file output from stablized from raw and check
    % that it is in our directory
    name=[pname, fnames{i}];
    if exist(name(1:end-12),'file')==2
        load(name(1:end-12),'frameshifts_strips_spline','timeaxis_secs','blinkframes')
    else
        missingmat{i}=name(1:end-12);
        continue
    end
    
    %below are to put it to the correct format. This is directly copied from 
    %find saccades 
    xy = frameshifts_strips_spline.*repmat(degrees_per_pixel,length(frameshifts_strips_spline),1);
    xy = xy - ones(length(xy),1) * mean(xy);
    xwrapped = reshape(xy(:,1), 16, length(xy)/16); % assumes 16 strips per frame 480 hz sampling 30 hz frame rate
    ywrapped = reshape(xy(:,2), 16, length(xy)/16);
    xwrapped = xwrapped - ones(16,1) * xwrapped(1,:);
    ywrapped = ywrapped - ones(16,1) * ywrapped(1,:);
    xwrapped_smoothed = mediansmooth(xwrapped', mediansmoothwindow)'; % transposed to smooth on the 2nd dimension, 
    ywrapped_smoothed = mediansmooth(ywrapped', mediansmoothwindow)';
    xyfixed = xy - [xwrapped_smoothed(:), ywrapped_smoothed(:)];
    
    % % %     we fix the end points to they don't go crazy (doesn't really work)
    xyfixed(2) = xyfixed(3) + 1;
    xyfixed(1) = xyfixed(2) - 1;
    xyfixed(end-1) = xyfixed(end-2) + 1;
    xyfixed(end) = xyfixed(end-1) - 1;
    xyfixed(2,2) = xyfixed(3,2) + 1;
    xyfixed(1,2) = xyfixed(2,2) - 1;
    xyfixed(end-1,2) = xyfixed(end-2,2) + 1;
    xyfixed(end,2) = xyfixed(end-1,2) - 1;
% % %     leyepositiondata is where I think there was a mix up before causing the velocity 
% % %     mismatch before 
% % %     Here we use a gaussian smooth filter
% % %     the second input to this function is the SD of the filter
%     xyfixed(outlierind,:)=[]; timeaxis_secs(outlierind)=[];
    eyevelxy = [0,0; diff(xyfixed) ./ repmat(diff(timeaxis_secs),1,2)];
    eyeacc= [0; (sqrt(diff(eyevelxy(:,1)).^2+diff(eyevelxy(:,2)).^2)) ./ diff(timeaxis_secs)];
    badacc=find(eyeacc>50000);
    xyfixed(badacc,:)=[];timeaxis_secs(badacc)=[];eyeacc(badacc)=[];
    eyevelxy = [0,0; diff(xyfixed) ./ repmat(diff(timeaxis_secs),1,2)];
    eyeaccxy = [0,0; diff(eyevelxy) ./ repmat(diff(timeaxis_secs),1,2)];
    eyepositiondata=smoothgaussr(xyfixed,gausssmoothsd);
    
% % %     in the future we could try using an sgolay filter or rloess to preserve
% % %     the velocity better like the commented line below
% % %     eyepositiondata=smoothdata(xyfixed,'sgolay');
    
% % %     we create the figure with the xy motion data and later we 
% % %     will add the starts and ends of saccades later on
%     figure(amps)
%     hold off
%     plot(timeaxis_secs,eyepositiondata(:,1),'r.-')
%     hold on
%     plot(timeaxis_secs,eyepositiondata(:,2),'b.-')
%     
    %calculation of eye velocity sqrt(diff(x^2)+diff(y^2))/diff(time)
    sampletimedata=timeaxis_secs;
    eyevelmedsmooth= [0; sqrt(diff(xyfixed(:,1)).^2+diff(xyfixed(:,2)).^2)./diff(sampletimedata(:))];
    %accmedsmooth= [0; abs(diff(eyevelmedsmooth)./diff(sampletimedata))];
    eyevelmedgausssmooth = [0; sqrt(diff(eyepositiondata(:,1)).^2+diff(eyepositiondata(:,2)).^2)./diff(sampletimedata(:))];
    
    %are just a check to make sure we are using the right
    %saccade times from the matvelocity figure
    if isempty(xdata)
        continue
    elseif size(xdata{1,1},2)>size(xdata{2,1},2)
        sac_click=sort(xdata{2,1});
    else
        sac_click=sort(xdata{1,1});
    end
    
    %loop through all the saccade times and turn them into time indexes
    sacfailvelthresh=zeros(size(sac_click));
    sacfailampthresh=zeros(size(sac_click));
    sactooclose=zeros(size(sac_click));
    for j=1:size(sac_click,2)
         saccadesample=find(timeaxis_secs==sac_click(j));
         if isempty(saccadesample)
             [d,saccadesample]=min(abs(timeaxis_secs-sac_click(j)));
         end
         
         %is the search window over which we will look for the start
         %and end of the saccade this one is really important!!!!
         samplesaroundsac=saccadesample-windowbefore:saccadesample+windowafter;
         
         %check to make sure that the period of the saccade doesn't
         %extend past the time axis
         if samplesaroundsac(end)>size(timeaxis_secs,1)
             indend=find(samplesaroundsac==size(timeaxis_secs,1)-2);
             samplesaroundsac(indend:end)=[];
         elseif samplesaroundsac(1)<1
             samplesaroundsac=3:samplesaroundsac(end);
         end
         
         %only gets called if the whole saccade is outside the timeaxis
         if isempty(samplesaroundsac)
             badmatsorvel{i}=fnames{i};
             peakvelateach(j)=[];
             ampsateach(j)=[];
             dirateach(j)=[];
             sacstarts(j)=[];
             sacends(j)=[];
         else
             
            %grab all the smoothed velocities and find the ones over a set 
            % threshold for starts and ends
            allvelsacs=abs(eyevelmedgausssmooth(samplesaroundsac(1):samplesaroundsac(end)));
            veloverthreshstart=find(allvelsacs>startthresh);%change this to change start thresh
            veloverthreshend=find(allvelsacs>endthresh); %change this to change end thresh
            
            %check if it never reach or only reached once the threshold for
            %the start of the microsaccade. then we simply set the beginning to 
            %be ~ 20 ms before that selection ***put flag below christy
            if isempty(veloverthreshstart) || veloverthreshstart(end)==1
                sacfailvelthresh(j)=saccadesample;
                samplesaroundsac(1:windowbefore-defaultstart)=[];
            %otherwise we want to go 10ms before where it reaches the
            %velocity threshold over our window
            else
                threshvelstart=allvelsacs(veloverthreshstart);
                samplesaroundsac(1:veloverthreshstart(1)-8)=[];
                
                %Below is something I was trying but didn't work better for
                %now if two saccades are in our starting window then we
                %only use one of them.
%                 if isempty(find(diff(find(diff(threshvelstart)>startthresh))>1))~=1 %weird way of checking whether the window is too big
%                     if exist('sacstarts', 'var') ~=0
%                         if min(abs(sacstarts-samplesaroundsac(1)))<tooclosethresh %check if there is another saccade marked in the window
%                         %if there is we want to choose the other part of the saccade
%                             figure
%                             plot(allvelsacs)
%                             [val,startsind]=min(abs(sacstarts-samplesaroundsac(1)));
%                             inflectionpoints=find(diff(veloverthreshstart)>1);
%                             sacstarts(startsind)=samplesaroundsac(veloverthreshstart(inflectionpoints(1)+1));
%                             samplesaroundsac(1:veloverthreshstart(1)-8)=[];
%                         end
%                     end
%                 end
            end
                
            
            %check if it never reached or only reached once the
            %threshold for the end of the microsacade then we simply set 
            %the end to be ~ 30 ms after the selection. ***Look at this and
            %put flag here Christy
            if  isempty(veloverthreshend) || veloverthreshend(end)==1
                sacfailvelthresh(j)=saccadesample;
                samplesaroundsac(windowbefore+defaultend:end)=[];
                tooslow=tooslow+1
                
            %if not we check if we check the times that it reached that
            %threshold then went back below the threshold. Here the first
            %statement checks if it never stayed below that threshold for
            %more than 1 sample (needed incase there was some interpolation
            %during the saccade). We then go ~15ms after that point to deal
            %with lens wobble
            elseif isempty(find(diff(veloverthreshend)>2))
                if  veloverthreshend(end)+lenswobblelen<samplesaroundsac(end)-samplesaroundsac(1)
                    samplesaroundsac(veloverthreshend(end)+lenswobblelen:end)=[];
                end
                
            %if not we check if we check the times that it reached that
            %threshold then went below the threshold. Here the first
            %statement checks that it did stay below that threshold for
            %more than 1 sample which is indicative of a second saccade
            %within the search window. We then go ~15ms after that point to 
            %deal with lens wobble 
            elseif isempty(find(diff(veloverthreshend)>2))~=1
                drop=find(diff(veloverthreshend)>2);
                samplesaroundsac(veloverthreshend(drop(1))+lenswobblelen:end)=[];
               
            %if it only reached that threshold velocity toward the end of
            %the 100msec interval then we leave it alone and define starts
            %and ends
            end
            
            %let's check if this saccade is too close to another saccade
            %first 2 are to check if starts to starts are too close. This is probably
            %overkill but I like to be safe. The 23 here comes from 30ms
            %for the saccade to finish and 15ms for the lens wobble 
            if exist('sacstarts', 'var') ~=0
            if min(abs(sacstarts-samplesaroundsac(1)))<tooclosethresh
                sactooclose(j)=saccadesample;
                tooclosestart=tooclosestart+1
            elseif min(abs(sacends-samplesaroundsac(end)))<tooclosethresh
                sactooclose(j)=saccadesample;
                toocloseend=toocloseend+1
            end
            end
            sacstarts(j)=samplesaroundsac(1);
            sacends(j)=samplesaroundsac(end);
            
            %find the peak velocity between the starts and ends of the
            %saccades then find the amplitude and direction of the saccade
            %using the starts and ends
            peakvelateach(j)=max(abs(eyevelmedsmooth(samplesaroundsac(1):samplesaroundsac(end))));
            ampsateach(j)=sqrt((eyepositiondata(samplesaroundsac(end),1)-eyepositiondata(samplesaroundsac(1),1)).^2+(eyepositiondata(samplesaroundsac(end),2)-eyepositiondata(samplesaroundsac(1),2)).^2);
            dirateach(j)=atan2((eyepositiondata(samplesaroundsac(end),2)-eyepositiondata(samplesaroundsac(1),2)),(eyepositiondata(samplesaroundsac(end),1)-eyepositiondata(samplesaroundsac(1),1)));
            ampxateach(j)=ampsateach(j).*abs(cos(dirateach(j)));
            ampyateach(j)=ampsateach(j).*abs(sin(dirateach(j))); 
            accateach(j)=max(abs(eyeacc(samplesaroundsac(1):samplesaroundsac(end))));
            velxateach(j)=max(abs(eyevelxy(samplesaroundsac(1):samplesaroundsac(end),1)));
            velyateach(j)=max(abs(eyevelxy(samplesaroundsac(1):samplesaroundsac(end),2)));
            accyateach(j)=max(abs(eyeaccxy(samplesaroundsac(1):samplesaroundsac(end),2)));
            accxateach(j)=max(abs(eyeaccxy(samplesaroundsac(1):samplesaroundsac(end),1)));
            %here we throw out saccades under the minimum amplitude
            if ampsateach(j)<ampthresh
                 sacfailampthresh(j)=saccadesample;
                 toosmall=toosmall+1
                 peakvelateach(j)=[]; ampsateach(j)=[]; dirateach(j)=[];
                 sacstarts(j)=[]; sacends(j)=[];
                 ampyateach(j)=[]; ampxateach(j)=[]; 
                 velyateach(j)=[]; velxateach(j)=[]; 
                 accateach(j)=[]; accxateach(j)=[]; accyateach(j)=[];
                 
            elseif sactooclose(j)>1
                 peakvelateach(j)=[]; ampsateach(j)=[]; dirateach(j)=[];
                 sacstarts(j)=[]; sacends(j)=[];
                 ampyateach(j)=[];ampxateach(j)=[];
                 velyateach(j)=[];velxateach(j)=[];
                 accateach(j)=[]; accxateach(j)=[]; accyateach(j)=[];
            end
            
                        %plot the starts and ends on the x motion
         end
    end
       %now we are outside of the sactimes loop
    if exist('sacstarts', 'var') ~=0
        sacstarts=sacstarts(find(sacstarts));
        sacends=sacends(find(sacends));
        peakvelateach=peakvelateach(find(peakvelateach));
        ampsateach=ampsateach(find(ampsateach));
        dirateach=dirateach(find(dirateach));
        ampyateach=ampyateach(find(ampyateach));
        ampxateach=ampxateach(find(ampxateach));
        accateach=accateach(find(accateach));
        accxateach=accxateach(find(accxateach));
        accyateach=accyateach(find(accyateach));
        velyateach=velyateach(find(velyateach));
        velxateach=velxateach(find(velxateach));
%         figure(amps)
%         hold on
%         plot(timeaxis_secs(sacstarts),eyepositiondata(sacstarts,1),'go')
%         plot(timeaxis_secs(sacends),eyepositiondata(sacends,1),'ro')
    end

    %turn them into collumns for writing to excel sheet
    sacpeakvel=transpose(peakvelateach(1:size(sacstarts,2)));
    sacamps=transpose(abs(ampsateach(1:size(sacstarts,2))));
    sacampx=transpose(abs(ampxateach(1:size(sacstarts,2))));
    sacampy=transpose(abs(ampyateach(1:size(sacstarts,2))));
    sacpeakacc=transpose(abs(accateach(1:size(sacstarts,2))));
    sacdir=transpose(abs(dirateach(1:size(sacstarts,2))));
    sacvelx=transpose(abs(velxateach(1:size(sacstarts,2))));
    sacvely=transpose(abs(velyateach(1:size(sacstarts,2))));
    sacaccx=transpose(abs(accxateach(1:size(sacstarts,2))));
    sacaccy=transpose(abs(accyateach(1:size(sacstarts,2))));
    %% all the square waves!
    
    xampabovethresh=find(sacamps>0.1 & sacampy<sacampx);
    diffdirind=find(abs(diff(cos(dirateach(xampabovethresh))))>1);%find all the ones where the sign flips in the x direction
     if size(diffdirind,2)~=0
           diffdirind=[diffdirind, diffdirind(end)+1];%add one to the end so that it has the same number of elements
           timebdiffdir=abs(diff(timeaxis_secs(sacstarts(diffdirind))));%find the time between the microsaccades that swtich direction
           refractorytimes=find(timebdiffdir<0.5);
           squarewaves=diffdirind(refractorytimes);%find all the ones that have amplitudes and times that meet the previous criteria     
           numsquarewaves=size(squarewaves,2);
           averefractory=mean(refractorytimes);
           sqwjtrainnum=sum(diff(find(diff(squarewaves)==1))==1);
       else
           timebdiffdir=[];
           refractoryind=[];
           sameampdiffdir=[];
           squarewaves=[];
           numsquarewaves=0;               
           sqwjtrains=[];
           sqwjtrainnum=0;
     end
    %now lets calculate square wave stats
    sacpeakvelall{i}=sacpeakvel; sacpeakampall{i}=sacamps; sacdirall{i}=sacdir;
    numsac(i)=size(sacdir, 1);    numsacraw(i)=size(sac_click,2);
    meansacpeakvel(i)=mean(sacpeakvel); meansacamp(i)=mean(sacamps);
    meansacampx(i)=mean(sacampx); meansacampy(i)=mean(sacampy);
    meansacvelx(i)=mean(sacvelx); meansacvely(i)=mean(sacvely);
    meanacc(i)=mean(sacpeakacc); 
    meanaccx(i)=mean(sacaccx); meanaccy(i)=mean(sacaccy);
    numsqw(i)=numsquarewaves;
    numsqwtrain(i)=sqwjtrainnum;
    slope(i)=mean(sacampy./sacampx);
    
    %save the number of blinkframes       
    if isempty(blinkframes)
        numblinks(i)=0;
    elseif isempty(find(diff(blinkframes)>1))~=0
        numblinks(i)=1;
    else
        numblinks(i)=1+size(find(diff(blinkframes)>1),1);
    end
    %save the starts and ends figure to a different directory
    %saveas(amps,['C:\Users\Lenovo\Documents\matfigs\',fnames{i},'.tiff'],'tiff')

    %save frameshifts_strips_spline, eyepositiondata, timeaxis_secs, blinkframes, 
    %sacclick, sacstarts, sacends to new mat file for easy access and read
    idformat=idnum(i); eyeformat=eye(i);tracenumformat=tracenum(i); msformat=ms(i); ageformat=age(i);
    save(['C:\Users\Lenovo\Documents\matfixed\',fnames{i},'.mat'], 'idformat', 'xyfixed','eyevelxy','eyeaccxy', 'eyeformat', 'tracenumformat', 'msformat', 'ageformat', 'frameshifts_strips_spline', 'eyepositiondata', 'eyevelmedsmooth', 'eyevelmedgausssmooth', 'timeaxis_secs','blinkframes','sac_click','sacstarts','sacends','sacfailvelthresh','sacfailampthresh','sactooclose','numsqw')

end
%write this all to the second sheet of the initial excel file
idnum(size(numsac,2)+1:end)=[]; eye(size(numsac,2)+1:end)=[]; tracenum(size(numsac,2)+1:end)=[];
idnum=idnum(find(numsac)); eye=eye(find(numsac)); tracenum=tracenum(find(numsac));ms=ms(find(numsac));age=age(find(numsac));
meansacpeakvel=meansacpeakvel(find(numsac));meansacamp=meansacamp(find(numsac));numblinks=numblinks(find(numsac));
meansacampx=meansacampx(find(numsac));meansacampy=meansacampy(find(numsac));meanacc=meanacc(find(numsac));
meansacvelx=meansacvelx(find(numsac));meansacvely=meansacvely(find(numsac));
meansacaccx=meanaccx(find(numsac));meansacaccy=meanaccy(find(numsac));
numsacraw=numsacraw(find(numsac)); numsac=numsac(find(numsac)); 
slope=slope(find(numsac));
numsqw=numsqw(find(numsac)); numsqwtrain=numsqwtrain(find(numsac));
xlswrite(filenameforexcel,{'id','eye','tracenum','ms','age','numsac','meanvel','meanamp','numblinks','meanacc','meanvelx','meanvely','meanampx','meanampy','meanaccx','meanaccy','sqwnum','sqwtrainmax','slope'},2,'A1')
xlswrite(filenameforexcel,[idnum',eye',tracenum',ms',age',numsac',meansacpeakvel',meansacamp',numblinks',meanacc',meansacvelx',meansacvely',meansacampx',meansacampy',meansacaccx',meansacaccy',numsqw',numsqwtrain',slope'],2,'A2')