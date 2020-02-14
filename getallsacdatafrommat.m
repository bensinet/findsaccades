clearvars; close all;
%select all the matvelocity figures
[fnames,pname] = uigetfile('*mat','Open AVI files','MultiSelect','On');
filenameforexcel='C:\Users\Lenovo\Documents\both patients_updated2.21.19 -findendsac.xlsx';
[num,text,raw]=xlsread(filenameforexcel,1);
ids=num(:,1);
agecollumn=6;
mscollumn=2;
    velx=[];vely=[]; accx=[]; accy=[];
    amps=[]; vel=[]; acc=[]; dir=[];id=[];msfull=[]; agefull=[];sacstartall=[];sacendall=[];tracenum=[];eye=[];
for i=1:size(fnames,2)
    idi=str2num(fnames{i}(1:5));
    eyel=fnames{i}(6);
    if eyel=='L' || eyel=='l'
        eyei=1;
    elseif eyel=='R' || eyel=='r'
        eyei=2;
    end
    tracenumi=str2double(fnames{i}(10:11));
    findplace=find(ids==idi);
    load([pname,fnames{i}])
    sacsize=size(sacstarts,2);
    if isempty(findplace)
        continue
    else
    end
    ms=num(findplace,mscollumn);
    age=num(findplace,agecollumn);
    for j=1:sacsize
        vel1(j)=max(abs(eyevelmedsmooth(sacstarts(j):sacends(j))));
        acc1(j)=max(abs(eyeaccxy(sacstarts(j):sacends(j))));
        velx1(j)=max(abs(eyevelxy(sacstarts(j):sacends(j),1)));
        vely1(j)=max(abs(eyevelxy(sacstarts(j):sacends(j),2)));
        accx1(j)=max(abs(eyeaccxy(sacstarts(j):sacends(j),1)));
        accy1(j)=max(abs(eyeaccxy(sacstarts(j):sacends(j),2)));
    end
    if sacsize>1
        amps(end+1:end+sacsize)=sqrt((eyepositiondata(sacends,1)-eyepositiondata(sacstarts,1)).^2+(eyepositiondata(sacends,2)-eyepositiondata(sacstarts,2)).^2);
        vel(end+1:end+sacsize)=vel1(1:sacsize);
        acc(end+1:end+sacsize)=acc1(1:sacsize);
        dir(end+1:end+sacsize)=atan2((eyepositiondata(sacends,2)-eyepositiondata(sacstarts,2)),(eyepositiondata(sacends,1)-eyepositiondata(sacstarts,1)));
        id(end+1:end+sacsize)=ones(1,sacsize).*ids(findplace);
        tracenum(end+1:end+sacsize)=ones(1,sacsize).*tracenumi;
        eye(end+1:end+sacsize)=ones(1,sacsize).*eyei;
        msfull(end+1:end+sacsize)=ones(1,sacsize).*ms;
        agefull(end+1:end+sacsize)=ones(1,sacsize).*age;
        sacstartall(end+1:end+sacsize)=timeaxis_secs(sacstarts);
        sacendall(end+1:end+sacsize)=timeaxis_secs(sacends);
        velx(end+1:end+sacsize)=velx1(1:sacsize);
        vely(end+1:end+sacsize)=vely1(1:sacsize);
        accx(end+1:end+sacsize)=accx1(1:sacsize);
        accy(end+1:end+sacsize)=accy1(1:sacsize);
    elseif sacsize==1
        amps(end+1)=sqrt((eyepositiondata(sacends,1)-eyepositiondata(sacstarts,1)).^2+(eyepositiondata(sacends,2)-eyepositiondata(sacstarts,2)).^2);
        dir(end+1)=atan2((eyepositiondata(sacends,2)-eyepositiondata(sacstarts,2)),(eyepositiondata(sacends,1)-eyepositiondata(sacstarts,1)));
        vel(end+1)=vel1(1);
        acc(end+1)=acc1(1);
        msfull(end+1)=ms;
        agefull(end+1)=age;
        id(end+1)=ids(findplace);
        tracenum(end+1)=tracenumi;
        eye(end+1)=eyei;
        sacstartall(end+1)=sacstarts;
        sacendall(end+1)=sacends;
        velx(end+1)=velx1(1);
        vely(end+1)=vely1(1);
        accx(end+1)=accx1(1);
        accy(end+1)=accy1(1);
    end
end
xlswrite(filenameforexcel,{'id','eye','tracenum','ms','age','amps','vel','acc','dir','sacstartall','sacendall','velx','vely','accx','accy'},3,'A1')
xlswrite(filenameforexcel,[id',eye',tracenum',msfull',agefull',amps',vel',acc',dir',sacstartall',sacendall',velx',vely',accx',accy'],3,'A2')