clearvars;
filenameforexcel='C:\Users\Lenovo\Documents\both patients_updated2.21.19 -findendsac.xlsx';
% filenamenosac='C:\Users\Lenovo\Documents\no microsaccades.xlsx';
[numms,textms,rawms]=xlsread(filenameforexcel,2);
% [nomic,crap,crapola]=xlsread(filenamenosac);
% nomicid=nomic(:,1);nomiceye=nomic(:,2); nomictrace=nomic(:,3);

idnmall=numms(:,1);
eyeall=numms(:,2);
idfull=idnmall*10+eyeall;
[UA,~,idx] = unique(idfull);
for i=1:size(UA,1)
    numtraces(i)=size(find(idfull==UA(i)),1);
end
numsacall = accumarray(idx,numms(:,6),[],@mean);
meansacampall = accumarray(idx,numms(:,8),[],@mean); meansacpeakvelall = accumarray(idx,numms(:,7),[],@mean);
idacc = accumarray(idx,numms(:,1),[],@mean); eyeacc= accumarray(idx,numms(:,2),[],@mean); msacc=accumarray(idx,numms(:,4),[],@mean); 
accacc=accumarray(idx,numms(:,10),[],@mean); meansacaccx=accumarray(idx,numms(:,15),[],@mean); meansacaccy=accumarray(idx,numms(:,16),[],@mean);
meansacampx=accumarray(idx,numms(:,13),[],@mean); meansacampy=accumarray(idx,numms(:,14),[],@mean);
meansacvelx=accumarray(idx,numms(:,11),[],@mean); meansacvely=accumarray(idx,numms(:,12),[],@mean);
numtraces=transpose(numtraces);
right=find(eyeacc==2);left=find(eyeacc==1);

idrall=idacc(right);idlall=idacc(left);
msrall=msacc(right);mslall=msacc(left);
trall=numtraces(right);tlall=numtraces(left);
numrall=numsacall(right);numlall=numsacall(left);
amprall=meansacampall(right);amplall=meansacampall(left);
ampxrall=meansacampx(right);ampxlall=meansacampx(left);
ampyrall=meansacampy(right);ampylall=meansacampy(left);
velxrall=meansacvelx(right);velxlall=meansacvelx(left);
velyrall=meansacvely(right);velylall=meansacvely(left);
accxrall=meansacaccx(right);accxlall=meansacaccx(left);
accyrall=meansacaccy(right);accylall=meansacaccy(left);
accrall=accacc(right);acclall=accacc(left);
velrall=meansacpeakvelall(right);vellall=meansacpeakvelall(left);
idxr=ismember(idrall,idlall); idxl=ismember(idlall,idrall);
idr=idrall(idxr);idl=idlall(idxl);
msr=msrall(idxr);
numtracesr=trall(idxr);numtracesl=tlall(idxl);
numsacallr=numrall(idxr);numsacalll=numlall(idxl);
meansacampallr=amprall(idxr);meansacampalll=amplall(idxl);
velxr=velxrall(idxr);velxl=velxlall(idxl);
velyr=velyrall(idxr);velyl=velylall(idxl);
ampxr=ampxrall(idxr);ampxl=ampxlall(idxl);
ampyr=ampyrall(idxr);ampyl=ampylall(idxl);
accr=accrall(idxr);accl=acclall(idxl);
accxr=accxrall(idxr);accxl=accxlall(idxl);
accyr=accyrall(idxr);accyl=accylall(idxl);
meansacpeakvelallr=velrall(idxr);meansacpeakvelalll=vellall(idxl);
aveind=find(numtracesr<2 | numtracesl<2);
avevel=mean([meansacpeakvelallr,meansacpeakvelalll],2);
aveamp=mean([meansacampallr,meansacampalll],2);
avenummicro=mean([numsacallr,numsacalll],2);
aveacc=mean([accr,accl],2);
ampx=mean([ampxr,ampxl],2); ampy=mean([ampyr,ampyl],2);
velx=mean([velxr,velxl],2); vely=mean([velyr,velyl],2);
accx=mean([accxr,accxl],2); accy=mean([accyr,accyl],2);
avevel(aveind)=[];aveamp(aveind)=[];avenummicro(aveind)=[];idr(aveind)=[];msr(aveind)=[];ampx(aveind)=[];
ampy(aveind)=[]; aveacc(aveind)=[]; velx(aveind)=[];vely(aveind)=[];accx(aveind)=[];accy(aveind)=[];
% xlswrite(filenameforexcel,{'id','MS','Leftnumtraces','Leftnumsac','Leftmeanamp','Leftmeanvel','Rightnumtraces','Rightnumsac','Rightmeanamp','Rightmeanvel','Meannumsac','Meanamp','Meanvel'},4,'A1')
% xlswrite(filenameforexcel,[idr,msr,numtracesl,numsacalll,meansacampalll,meansacpeakvelalll,numtracesr,numsacallr,meansacampallr,meansacpeakvelallr,avenummicro,aveamp,avevel],4,'A2')

[num,text,raw]=xlsread(filenameforexcel,1);
idnmall=num(:,1);
for i=1:size(idr,1)
idx=find(idnmall==idr(i));
velallp{idx}=avevel(i); velxp{idx}=velx(i); velyp{idx}=vely(i);
ampallp{idx}=aveamp(i); ampxp{idx}=ampx(i); ampyp{idx}=ampy(i);
idnmp{idx}=idr(i); 
accp{idx}=aveacc(i);accxp{idx}=accx(i); accyp{idx}=accy(i);



end
xlswrite(filenameforexcel,{'id','meanamp','meanvel','meanacc','ampx','ampy','velx','vely','accx','accy'},4,'A1');
xlswrite(filenameforexcel,[idnmp',ampallp',velallp',accp',ampxp',ampyp',velxp',velyp',accxp',accyp'],4,'A2');
    