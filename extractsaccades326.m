function sacstats = extractsaccades218(eyepositiondata, sampletimedata, sacthresh, verbose)
% function sacstats = extractsaccades(eyepositiondata, sampletimedata);
% extractsaccades finds locations and amplitudes of saccades from an eye
% trace

%% defaults
% if verbose is not supplied, set to zero for "silent" operation
if nargin < 4;
    verbose = 0;
end
%we override this and set it to 10 in findsaccades
% set sacthresh to 5 deg/sec as default
if nargin < 3;
    sacthresh(1) = 5;
    sacthresh(2) = 1000;
end
% provide an upper limit for saccade velocity as well.
if numel(sacthresh) == 1;
    sacthresh(2) = 400;
end
% if sampletimedata is not supplied, assume 30hz continuous sampling
if nargin < 2;
    sampletimedata = [1:length(eyepositiondata(:,1))]./30;
end

%% fix the two ends so that the velocity crosses zero
eyepositiondata(2,1) = eyepositiondata(3,1) + 1;
eyepositiondata(1,1) = eyepositiondata(2,1) - 1;
eyepositiondata(end-1,1) = eyepositiondata(end-2,1) + 1;
eyepositiondata(end,1) = eyepositiondata(end-1,1) - 1;
eyepositiondata(2,2) = eyepositiondata(3,2) + 1;
eyepositiondata(1,2) = eyepositiondata(2,2) - 1;
eyepositiondata(end-1,2) = eyepositiondata(end-2,2) + 1;
eyepositiondata(end,2) = eyepositiondata(end-1,2) - 1;
%% calculate velocity and find saccades
% if diff(sampletimedata)>0.003
eyevelx= [0; diff(eyepositiondata(:,1))./diff(sampletimedata(:))]; %adding in this extra 0, give you back a place you lost taking the diff
eyevely= [0; diff(eyepositiondata(:,2))./diff(sampletimedata(:))];
%here we combine x and y
eyevel = [0; sqrt(diff(eyepositiondata(:,1)).^2+diff(eyepositiondata(:,2)).^2)./diff(sampletimedata(:))];
%here we find the unique indices
sacindices = find((abs(smoothgaussr(eyevel,3)) > sacthresh(1)) .*( abs(eyevel) < sacthresh(2)));

try
    if ~isempty(sacindices);
        %get one index for each saccade
        uniquesacindicesnew = sacindices(find(diff([-1; sacindices])>1 )); %& diff([-1; sacindices])~=15) & diff([-1; sacindices])~=16 & diff([-1; sacindices])~=15 & diff([-1; sacindices])~=32)
        if verbose > 0;
            sampleindices=1:size(eyepositiondata(:,1));
            newhandle=figure(300);hold off;
            plot(sampletimedata,eyepositiondata(:,1),'r.-');hold on;
            plot(sampletimedata,eyepositiondata(:,2),'b.-');
            plot(sampletimedata(sampleindices(uniquesacindicesnew)),eyepositiondata(uniquesacindicesnew),'bo');hold off;
            %plot where it thinks the indexes of saccades are
            k=0;
            l=0;
            button=0;
            while button ~= 13
            [x,y,button]=ginputc(1);
                    if button==1 %left click add a saccade where you clicked
                        [~,newindices]=min(abs(sampletimedata-x));
                        l=l+1;
                        uniquesacindicesnew(length(uniquesacindicesnew)+1)=newindices;
                        figure(newhandle); hold off;
                        plot(sampletimedata,eyepositiondata(:,1),'r.-');hold on;
                        plot(sampletimedata,eyepositiondata(:,2),'b.-');
                        plot(sampletimedata(sampleindices(uniquesacindicesnew)),eyepositiondata(uniquesacindicesnew),'bo');hold off;
                    end
                    if button==3 %right click remove the nearest saccade to where you clicked
                        k=k+1;
                        [~,removeindices]=min(abs(sampletimedata-x));
                        [~, index] = min(abs(uniquesacindicesnew-removeindices));
                        uniquesacindicesnew(index)=[];
                        figure(newhandle); hold off;
                        plot(sampletimedata,eyepositiondata(:,1),'r.-');hold on;
                        plot(sampletimedata,eyepositiondata(:,2),'b.-');
                        plot(sampletimedata(sampleindices(uniquesacindicesnew)),eyepositiondata(uniquesacindicesnew),'bo');hold off;
                    end
                    if button==101 %clear all the saccade indexes to start over by hitting 'e'
                        uniquesacindicesnew=[];
                        figure(newhandle); hold off;
                        plot(sampletimedata,eyepositiondata(:,1),'r.-');hold on;
                        plot(sampletimedata,eyepositiondata(:,2),'b.-');
                        plot(sampletimedata(sampleindices(uniquesacindicesnew)),eyepositiondata(uniquesacindicesnew),'bo');hold off;
                    end
            end
            positionhandle = figure(301);hold off;
            %press return when done selecting your positions
            velocityhandle = figure(302);hold off;
            directionalhandle=figure(303);hold off;
            positioncleaned=figure(304); hold off;
        end
        %% find start and ends of saccades
        %this does it very simply it's designed to go over the end of the
        %saccade to avoid lens wobble
        sacstarts = uniquesacindicesnew-20; %start 10 samples before
        sacends=uniquesacindicesnew+25; %end 50 samples after 
        for l=1:size(sacends,1)
            if sacends(l)>=size(sampletimedata,1)
                sacends(l)=size(sampletimedata,1);
            end
            if sacstarts(l)<1
                sacstarts(l)=1;
            end
        end
         %% find frequency velocity and direction of saccades
        sacfreq = length(uniquesacindicesnew) ./ (sampletimedata(end) - sampletimedata(1));
        sacamps = sqrt((eyepositiondata(sacends,2)-eyepositiondata(sacstarts,2)).^2 + (eyepositiondata(sacends,1)-eyepositiondata(sacstarts,1)).^2);
        sacvel=zeros(size(sacamps));       
        %now we need to loop through to determine the velocity of all the
        %individual microsaccades
        for i = 1:size(sacamps)
            sacvel(i)=max(abs(eyevel(sacstarts(i):sacends(i))));
            %find the maximum velocity of each saccade
            sacdirection(i)=atan2((eyepositiondata(sacends(i),2)-eyepositiondata((sacstarts(i)),2)),(eyepositiondata(sacends(i),1)-eyepositiondata((sacstarts(i)),1)));
            %find the direciton of each saccade
        end
        %% find drift amplitude and velocity
        shiftinmoviex=eyepositiondata(end-5,1)-eyepositiondata(5,1);
        shiftinmoviey=eyepositiondata(end-5,2)-eyepositiondata(5,2);
        sacx=eyepositiondata(sacends,1)- eyepositiondata(sacstarts,1);
        sacy=eyepositiondata(sacends,2)- eyepositiondata(sacstarts,2);
        driftx=shiftinmoviex-sum(sacx);
        drifty=shiftinmoviey-sum(sacy);
        driftdirection=atan2(drifty,driftx);
        driftamp=sqrt(driftx.^2+drifty.^2);
        eyevelindices=1:size(eyevel,1);
        driftindices=setdiff(eyevelindices,sacindices);
        driftvelocity=mean(abs(eyevel(driftindices)));
        %driftvelocity doesn't reveal much yet
                %% make the plots
        if verbose > 0;
            figure(positionhandle);hold on;
            plot(sampletimedata, eyepositiondata(:,1),'r.-')
            plot(sampletimedata,eyepositiondata(:,2),'b.-')
            plot(sampletimedata(sacstarts),eyepositiondata(sacstarts,1),'gs');
            plot(sampletimedata(sacends),eyepositiondata(sacends,1),'ms');%hold off
            plot(sampletimedata(sacstarts),eyepositiondata(sacstarts,2),'gs');
            plot(sampletimedata(sacends),eyepositiondata(sacends,2),'ms'); hold off;
            %plot the starts in the x and y in green and the ends in
            %magenta
            figure(positioncleaned);hold on;
            plot(sampletimedata, eyepositiondata(:,1),'r.-');
            plot(sampletimedata,eyepositiondata(:,2),'b.-');hold off
            %             plot(sampletimedata(sacstarts),eyepositiondata(sacstarts,2),'gs');
            %             plot(sampletimedata(sacends),eyepositiondata(sacends,2),'ms');
            %             figure(velocityhandle);
            %             plot(sampletimedata(sacstarts),eyevel(sacstarts),'gs');
            %             plot(sampletimedata(sacends),eyevel(sacends),'ms');
            figure(directionalhandle)
            polar(sacdirection(:),sacamps(:),'ro')
            figure(velocityhandle); hold on
            plot(sampletimedata,eyevel,'r.-');
            %polar plot of the directions of each saccade
            hold on;
            %plot(sampletimedata(sacindices),eyevel(sacindices),'bo');
            plot(sampletimedata(uniquesacindicesnew),eyevel(uniquesacindicesnew),'gx'); hold off; %latency is in uniquesacindicesnew
            figure(798);
            plot(sacx,sacy,'bs');
            figure(999);
            plot(sampletimedata,eyevelx,'r');hold on; 
            plot(sampletimedata,eyevely,'b')
            %velocity plot
        end
    else
        sacfreq = 0; sacamps = 0; uniquesacindicesnew = 1; sacstarts = 1; sacends = 1;
    end
    disp('it worked');
catch
    rethrow(lasterror);
    keyboard;
end
sacstats = {size(sacamps,1),mean(abs(sacvel)),max(abs(sacvel)),mean(abs(sacamps)),sacfreq,driftamp,driftdirection,driftvelocity,sacdirection,sacamps,uniquesacindicesnew,k,l,sacvel};
%output your sacstats
