%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Code written by Ovidiu Radulescu, University of Montpellier, June 2019
%%%%% reads the short movie data, filters large intensity and/or uncalibrated cells, applies a 
%%%%% low-to-zero thresholding filter to the signal, eliminates cells with too large intensities
%%%%% output :  
%%%%%  >> .mat file containing variables 'DataExp', 'Frames', 'Samples', 
%%%%% representing intensity data, time data, and indices of
%%%%% selected cells
%%%%%% >> a number of pdf files showing the data
%%%%%% >> a text file description.txt containing the selected cell file
%%%%%% names and indices
%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% USER PROVIDED INFORMATION
%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DataFilePath1 = 'images_short/'; %%% where to write images 
thr=0.5; %%%%% threshold intensity value; data smaller than thr is set to zero
MaxInt=150; %%%% maximum intensity value; cells with larger intensities are excluded
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


DataFilePath0 = 'data_short/'; %%%% where is data
nthr=strrep(num2str(thr),'.',''); %%%% string to be added to output file, representing threshold value


data_list=struct2cell(dir(fullfile(DataFilePath0)))
dir_name_list = data_list(1,3:end); %%% list of subdirectories containing data from different genotypes


for iii=1:length(dir_name_list)
    
close all    
%%%%%% where the data is %%%%%%%%%%%%%%%%%%%
dirname = dir_name_list{iii}; 
list=[];

DataFilePath = [DataFilePath0,dirname];

WriteTo =      [DataFilePath1,dirname,'_images'];



mkdir(WriteTo)
delete([WriteTo,'/*'])

data_list=struct2cell(dir(fullfile(DataFilePath)));
file_name_list = data_list(1,3:end);
nexp=length(file_name_list); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% create DataExp file   
DataExp=[]; Samples=[];



for data_i=1:nexp
    
    DataFileName = [DataFilePath,'/',file_name_list{data_i}];
    T = readtable(DataFileName,'Delimiter',';');
    M = table2array(T);
    Data = M(:,2);
    Frames = M(:,1);

   
   
%%%%% other parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frame_num = length(Frames);     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

%%%%% truncate data (all data points below threshold are set to zero)
DataTrunc=Data;
DataTrunc(DataTrunc<thr)=0;


%%%%% eliminate spikes (very large derivative) in the data %%%%%%%
Data_derivative = diff(DataTrunc);
ind_der_spike=find( abs(Data_derivative) > 20 * mean(abs(Data_derivative))); 
ind_spike=ind_der_spike(2:2:length(ind_der_spike));
for ii=1:length(ind_spike)
   DataTrunc(ind_spike(ii))= (DataTrunc(ind_spike(ii)-1)+DataTrunc(ind_spike(ii)+1))/2;
end



if isempty(list)
if max(DataTrunc) < MaxInt & (length(find(DataTrunc == 0)) < length(DataTrunc))  %%%% selects only cells with intensities smaller than Maxint and generating a signal different from zero   
    DataExp=[DataExp,DataTrunc]; %%% intensities on columns
    Samples=[Samples,data_i]; %%%% stores the selected cell number
end 
    else
if ~ismember(data_i,list)  %%%% selects cells not in list
    DataExp=[DataExp,DataTrunc]; %%% intensities on columns
    Samples=[Samples,data_i]; %%%% stores the selected cell number
end    
end
    
ifig=floor((data_i-1)/36)+1;
    
h=figure(ifig)
hold off
subplot(6,6,mod((data_i-1),36)+1)
plot(1:frame_num,Data,'k','linewidth',0.1)
hold on
plot(1:frame_num,DataTrunc,'r','linewidth',0.1)
%axis([1 500 0 100])
if ismember(data_i,list)
    title(num2str(data_i),'color','r')
else
    title(num2str(data_i),'color','k')
end
    

if mod(data_i,36)==0 || data_i == nexp
figfile = [WriteTo,'/','fig',num2str(ifig),'.pdf']
print(h,'-dpdf',figfile) 
end


end %%% nexp


%%%%%% visualise the data as it will be analysed
for i=1:length(Samples)
   ifig=floor((i-1)/36)+1;
h=figure(ifig+10)
hold off
subplot(6,6,mod((i-1),36)+1)
plot(1:frame_num,DataExp(:,i),'r','linewidth',0.1)
title(num2str(Samples(i)))
    if mod(i,36)==0 || i == length(Samples)
        figfile = [WriteTo,'/','fig',num2str(ifig+10),'.pdf']
    print(h,'-dpdf',figfile) 
    end     
end



name_list=file_name_list;
for i=1:nexp %%%% for all the cells (not only for those selected)
    name_list{i}=[num2str(i,3),': ',name_list{i}]; %%% store index and name of the file for all the cells
end
%%%%% write text file
fname = [WriteTo,'/','description.txt'] %%%% description will contain the list of cell indices and files
tfile = fopen(fname,'w')

fprintf( tfile, '%s \n',name_list{:})
fclose(tfile)




sd= size(DataExp);
Frames=0:sd(1)-1;

fname = ['data_',dirname,'_short_',nthr,'.mat'];
save(fname, 'DataExp', 'Frames', 'Samples') %%% Samples will contain the indices of selected cells




end %%% for iii
