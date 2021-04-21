
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Code written by Ovidiu Radulescu, University of Montpellier, June 2019
%%%%% reads the raw long movie data
%%%%% for inhomogeneous length of times series sort the
%%%%% data lines according to length and fill with value 2 up to the end,
%%%%% representing the end of the longest time series
%%%%% output :  
%%%%%  >> .mat file containing variables 'DataExpLong','DataExpLongRaw','Time' 
%%%%% representing raw intensity data, discretized (0,1,2) intensity data,
%%%%% and Time values
%%%%%% >> a number of pdf files showing the data
%%%%%% >> a text file description.txt containing the cell file
%%%%%% names and indices
%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

warning('off','all')

fsz=16; %%%% fontsize

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% USER PROVIDED INFORMATION
%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
threshold = 500; %%% threshold value intensity( intensity > threshold )=1; intensity( intensity<threshold ) =0
DataFilePath1 = 'images_long/'; %%% where to write images 
header = 21; %%% number of lines in the data file header
intensity_col = 5; %%% which column contains the intensity data
time_col = 1; %%% which column contains the time data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DataFilePath0 = 'data_long/'; %%%% where is data
data_list=struct2cell(dir(fullfile(DataFilePath0)));
dir_name_list = data_list(1,3:end); %%% list of subdirectories containing data from different genotypes

for iii=1:length(dir_name_list)


%%%% genotype iii    
    
dirname = dir_name_list{iii}; 





%%%% create DataExp file   
DataExp=[];

WriteTo = [DataFilePath1,dirname,'_long_raw_images',num2str(threshold)];
mkdir(WriteTo)

DataFilePath = [DataFilePath0,dirname]; %%%% where data is

data_list=struct2cell(dir(fullfile(DataFilePath))); %%%% look for available files

file_name_list = data_list(1,3:end); %%%% list of file names without paths

File_names_short= file_name_list; %%%% 

for jj=1:length(file_name_list)
    file_name_list{jj} = [DataFilePath,'/',file_name_list{jj}]; % add path to filenames
end

File_names= file_name_list; %%%% list of file names with paths




nexp=length(File_names); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%% first read to get L, Lmax 
Lmax=0;L=[];
for data_i=1:nexp
    T = readtable(File_names{data_i},'Delimiter','\t','HeaderLines',header);
    M = table2array(T);
    Time = M(:,time_col); %%%% time in minutes
    l=length(Time);
    L=[L,l];
    if l > Lmax
        Lmax=l;
        Time_max=Time;
    end
end


for data_i=1:nexp
    

    T = readtable(File_names{data_i},'Delimiter','\t','HeaderLines',header);
    M = table2array(T);
    Data = M(:,intensity_col);
    
    
    
    
    
    Time = M(:,time_col); %%%% time in minutes

    
    l=length(Time);
    
 if ~(l == Lmax) %%% time points less than Lmax
      Data=[Data;ones(Lmax-l,1)*(-1)]; %%% fill with -1 up to Lmax
      Time=Time_max;
 end
  
 DataTrunc=Data;
 DataTrunc(DataTrunc<threshold)=0;
 nn=size(Data);
 frame_num=nn(1);
 
 
% if max(Data) < MaxInt
 DataExp=[DataExp,Data]; %%% intensities on columns
%     Samples=[Samples,data_i];
% end 
    
 ifig=floor((data_i-1)/36)+1;
     
 h=figure(ifig)
 hold off
 subplot(6,6,mod((data_i-1),36)+1)
 plot(Time,Data,'k','linewidth',0.1)
 hold on
 %plot(1:frame_num,DataTrunc,'r','linewidth',0.1)
 ind=find(DataTrunc==0);
 XX=1:frame_num;
 plot(Time(ind),DataTrunc(ind),'kx','linewidth',0.2)
    axis([0, max(Time), 0, max(Data)])
 title(num2str(data_i))
% 
% 




 if mod(data_i,36)==0 || data_i == nexp
 figfile = [WriteTo,'/','fig',num2str(ifig),'.pdf']
 print(h,'-dpdf',figfile) 
 end

end %%% data_i




%end %%% jjj

DataExpLongRaw = DataExp;






%%%% apply thresholding
DataExp(DataExp <= threshold & DataExp > -1)=0;
DataExp(DataExp > threshold )=1;
DataExp(DataExp == -1)=2; %%% change -1 to 2; these are non-measured values to fill the data matrix when cell timeseries have uneven lengths 


X=DataExp'+1; %%%%%% used for vizualization, X \in {1,2,3}

DataExpLong = DataExp;
%%%%%% reorder according to number of time points
ind=find(L<Lmax+1); 
if ~isempty(ind)
   [Lsorted,isorted]=sort(L,'descend');
   X=X(isorted,:); %%%% is used for vizualization
   DataExpLong=DataExp(:,isorted); %%%% DataExpLong \in {0,1,2}
   DataExpLongRaw=DataExpLongRaw(:,isorted); %%%% 
   %File_names=File_names(isorted);
end


%%%%%% show thresholded and ordered data
h=figure(2000)
if ~isempty(ind)
map=[[0 0 1];[1 0 0];[1 1 1]];
else
map=[[0 0 1];[1 0 0]];    
end
colormap(map)
image(Time,1:length(File_names),X)
xlabel('Time [min]','fontsize',fsz)
ylabel('Transcription site','fontsize',fsz)
set(gca,'fontsize',fsz)

figfile = [WriteTo,'/','fig_binarized_',dirname,'.pdf']
print(h,'-dpdf',figfile) 

%%%% show raw data 
h=figure(3000)
colormap jet
%%%% upper truncation
imagesc(Time,1:length(File_names),DataExpLongRaw')
colorbar
xlabel('Time [min]','fontsize',fsz)
ylabel('Transcription site','fontsize',fsz)
set(gca,'fontsize',fsz)

figfile = [WriteTo,'/','fig_raw_',dirname,'.pdf']
print(h,'-dpdf',figfile)


for i=1:length(File_names_short)
    File_names_short{i}=[num2str(i,3),': ',File_names_short{i}];
end
% %%%%% write text file
 fname = [WriteTo,'/','description.txt']
 tfile = fopen(fname,'w')
 
 fprintf( tfile, '%s \n',File_names_short{:})
 fclose(tfile)
 
 fname = ['data_',dirname,'_long_raw',num2str(threshold),'.mat']
 save(fname, 'DataExpLong','DataExpLongRaw','Time')

 
 
 
 
end %%% for iii
