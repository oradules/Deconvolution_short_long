%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Code written by Ovidiu Radulescu, University of Montpellier, June 2019
%%%%% reads results of genetic algorithm for short movies, data from long
%%%%% movies
%%%%% needs 1) short movies decomvolution results  result_tat_name_short.mat
%%%%% 2) long movie data name_long.mat or long movie data name_long_raw.mat
%%%%% performs unconstrainted three exponentials fit, compute parameters of
%%%%% the three states model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all

%%%% file names %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
list_short = {{'result_tat_latency_cPosPred10.mat'}};
files_long  = {'data_tat_latency_long_raw500.mat'};
names = {'tat_latency'};

%%%%%%%%%%%%%%%%%%%%%%% Parameters
thresh = 500;
isel_vis=[1,2,3];
isel_vis_long=[1,2,3]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%% Parameters 
TaillePreMarq = 700; % 700 bases (ref publi)
TailleSeqMarq = 5800; % 2900 bases (ref publi)
TailleSeqMarq = 2900; % 2900 bases (ref publi)
EspaceInterPolyMin = 10; % space between polymerase in bases 
Polym_speed = 67; % average speed bases per second (Ref publi)
TaillePostMarq = 1600+100*Polym_speed; % 1600 bases (ref publi)
TaillePostMarq = 5300 + 67*100; % 
tstep = 3;
FreqEchImg = 1/tstep; % 1/3 image per second data time sampling
FreqEchSimu = 1/(EspaceInterPolyMin/Polym_speed); 

fsz=16; msz=5;lw=1;
outliers_long=1;
outliers_short=0;
visualize = 1;

for ibig = 2:2
    
    alpha = 0.9 - ibig * 0.3;
    


close all

 
DataFilePath0 = ['results_3exp_no_constraint_alpha',num2str(100*alpha)]; %%%% where to write results
mkdir(DataFilePath0);
%delete([DataFilePath0,'/*']); %%% delete previous result files

xlsfilename = [DataFilePath0,'/results',num2str(100*alpha),'.xlsx'];

line=1;
xlswrite(xlsfilename,{'Fname','OBJ','OBJM3','OBJONOFF','NS','NL','shift','cens','k34f','k23f','k32f','k21f','k12f','K3','K1p','K2p','K1m','K2m','lambda1','lambda2','lambda3',...
    'A1','A2','A3','mRNA'},1,'A1');

for iname=1:length(names)

    
name=names{iname};    
files_short=list_short{iname};




      
dirwrite=[DataFilePath0,'/',name,'_result']; %%%%% where to write result
mkdir(dirwrite);

%%%%%%%%%%%%%%%%% load deconvolution result; concatenate if there are several files %%%%
Gexp=[]; Gpred=[]; Gpos=[]; Gfit = [];
for jjj=1:length(files_short)
    fname = files_short{jjj};
    load(fname)
    Gexp=[Gexp,DataExp];
    Gpred=[Gpred,DataPred];
    Gpos=[Gpos,cPosPred];
    Gfit=[Gfit;Fit];
end
DataExp = Gexp;
DataPred = Gpred;
cPosPred = Gpos;
Fit = Gfit;
%%%%%%%%%%%%%%%% load long movie data

fname = files_long{iname};
load(fname); %% loads DataExpLong and Time


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sd=size(DataPred);
n_cells=sd(2); %%%% number of cells


%%%%%%%%%%%%% eliminate outliers handling short movie %%%%%%%%%%%%%
if outliers_short
    Ev=zeros(n_cells,1); %%% based on quartile of the numbers of polymerases 
for i=1:n_cells
    Ev(i) = length(cPosPred{i});
end
Q = quantile(Ev,3);
isel=find(Q(1)-1.5*(Q(3)-Q(1)) < Ev & Ev < Q(3) + 1.5*(Q(3)-Q(1)));
DataExp=DataExp(:,isel);
DataPred=DataPred(:,isel);
cPosPred=cPosPred(isel);
Fit=Fit(isel);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if outliers_short
n_cells=length(isel);
end 

%%%%% compute distribution of spacings short movie %%%%%%
dt=[]; dtc=[];

tmax=0;
for i=1:n_cells
    indices = (cPosPred{i})';
    times = indices / FreqEchSimu;
    Mtimes=max(times);
    if Mtimes > tmax
        tmax = Mtimes;
    end
end

for i=1:n_cells
    indices = (cPosPred{i})';
    times = indices / FreqEchSimu;
    lt = length(times);
    switch lt
        case 0
            dtc=[dtc;tmax]; 
        case 1
            dtc=[dtc;tmax-times(1)];
        otherwise
        dtimes = diff(times); %%%% uncensured intervals
        %dtimes=dtimes(dtimes>0);
        dt=[dt;dtimes(1:end)];
        if tmax > times(end)
            dtc=[dtc;tmax-times(end)];
        end
        if times(1) > 0
            dtc=[dtc;times(1)];
        end
    end
end
dtg=[dt;dtc];
censored_short=[zeros(size(dt));ones(size(dtc))];

%%%%%% visualize data 
if visualize
h=figure(1);
colormap(jet)
n=size(DataExp); 
stime=(0:n(1)-1)*3/60;
ncell=1:n(2);
imagesc(stime,ncell,DataExp')
colorbar
xlabel('Time [min]','fontsize',fsz)
ylabel('Transcription site','fontsize',fsz)
%title('Short movies data','fontsize',fsz)
set(gca,'fontsize',fsz)
    figfile=[dirwrite,'/figure_data_short_b.pdf'];
    print(h,'-dpdf',figfile)

h=figure(11);

for i=1:3
    subplot(3,1,i)
    hold off
 plot(stime,DataExp(:,isel_vis(i)),'r')
 hold on
 plot(stime,DataPred(:,isel_vis(i)),'k')
 
positions   =  cPosPred{isel_vis(i)};  
bppositions =  positions/FreqEchSimu*Polym_speed-(TaillePreMarq+TailleSeqMarq+TaillePostMarq); %%% positions in bp
fpositions  =  bppositions /Polym_speed/60; %%%% positions in minutes
m= max(DataPred(:,isel_vis(i)));
for j=1:length(positions)
  plot([fpositions(j),fpositions(j)],[0,m/10],'c')
  hold on
end
 axis([0, max(stime), 0, m]);
 if i==3
     xlabel('Time [min]','fontsize',fsz)
 end
 if i==2
     ylabel('Florescence intensity','fontsize',fsz)
 end
 set(gca,'fontsize',fsz)
end
    figfile=[dirwrite,'/figure_data_short_a.pdf'];
    print(h,'-dpdf',figfile)


h=figure(111);
jsel=isel_vis(2);
plot(stime,DataExp(:,jsel),'r','linewidth',2)
m= max(DataPred(:,jsel));
  xlabel('Time [min]','fontsize',fsz)
  ylabel('Fluorescence Intensity','fontsize',fsz)
  axis([0, max(stime), 0, m]);
set(gca,'fontsize',fsz)
    figfile=[dirwrite,'/figure_data_short_c.pdf'];
    print(h,'-dpdf',figfile)

h=figure(112);
positions   =  cPosPred{jsel};  
bppositions =  positions/FreqEchSimu*Polym_speed-(TaillePreMarq+TailleSeqMarq+TaillePostMarq); %%% positions in bp
fpositions  =  bppositions /Polym_speed/60; %%%% positions in minutes
hold off
for j=1:length(positions)
  plot([fpositions(j),fpositions(j)],[0,1],'k','linewidth',lw)
  hold on
end
 axis([0, max(stime), 0, 1]);
set(gca,'fontsize',fsz,'YTick',[],'YtickLabel',[])
xlabel('Time [min]','fontsize',fsz)
    figfile=[dirwrite,'/figure_data_short_d.pdf'];
    print(h,'-dpdf',figfile)


end %%% if visualize





%%%%%%%%%%%%% outliers handling long movie %%%%%%%%%%%%%
if outliers_long
n = size(DataExpLong);
FI=zeros(n(2),1);
for i=1:n(2)
    N0=length(find(DataExpLong(:,i)==0));
    N1=length(find(DataExpLong(:,i)==1));
    FI(i) = N0/(N0+N1);
end


Q = quantile(FI,3);
isel=find(Q(1)-1.5*(Q(3)-Q(1)) < FI & FI < min([Q(3) + 1.5*(Q(3)-Q(1)),1]) & FI > 0 );
%isel=find(Q(1)-5*(Q(3)-Q(1)) < FI & FI < 1 & FI < Q(3) + 5*(Q(3)-Q(1)));
DataExpLong=DataExpLong(:,isel);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if visualize
h=figure(21);
colormap(jet)
n=size(DataExpLongRaw); 
stime=(0:n(1)-1)*3;
ncell=1:n(2);
img = DataExpLongRaw';
img(img>10000)=10000;
imagesc(stime,ncell,img)
colorbar
xlabel('Time [min]','fontsize',fsz)
ylabel('Transcription site','fontsize',fsz)
set(gca,'fontsize',fsz)
    figfile=[dirwrite,'/figure_data_long_b.pdf'];
    print(h,'-dpdf',figfile)

h=figure(22);


for i=1:3
    subplot(3,1,i)
 y =   DataExpLongRaw(:,isel_vis_long(i)); 
 hold off
 plot(stime,y,'r')
 hold on
 y(y<thresh)=0;
 y(y>=thresh)=thresh;
 plot(stime,y,'k.')
 if i==3
     xlabel('Time [min]','fontsize',fsz)
 end
 if i==2
     ylabel('Florescence intensity','fontsize',fsz)
 end
 set(gca,'fontsize',fsz)
end
    figfile=[dirwrite,'/figure_data_long_a.pdf'];
    print(h,'-dpdf',figfile)    

h=figure(2);
map=[[0 0 1];[1 0 0];[1 1 1]];
colormap(map)
image(stime,ncell,DataExpLong'+1)   
xlabel('Time [min]','fontsize',fsz)
ylabel('Transcription site','fontsize',fsz)
set(gca,'fontsize',fsz)
figfile=[dirwrite,'/figure_data_long_c.pdf'];
    print(h,'-dpdf',figfile)    
    
end %%% if visualize   
   



    
store = DataExpLong;
n = size(store);
%%%%% estimate distributions for long movies %%%%%%%%%%%%%%%%%%%%%
wt=[];wtc=[];
Ninactive=0; Total=0; tmaxlomax=0;
for i=1:n(2) %%% all cells
    ilast = find(store(:,i) ==2);
    if isempty(ilast)
        tmaxlo=Time(end)*60;
    else
        tmaxlo=Time(min(ilast))*60;
    end

    
    ind=  find(store(:,i)==1);
    
    if ~isempty(ind)
    laps = diff(ind);
    Ninactive = Ninactive+length(find(laps>1))+1;
    wtimes=(laps(laps>1)-1)*tstep*60;
    wlast=tmaxlo - ind(end)*tstep*60; %%%% tmaxlo depends on cell
    wt=[wt;wtimes]; %%% waiting times in seconds
    if wlast > 0
    wtc=[wtc;wlast]; %%% add last interval
    end
    else
       %wtc=[wtc;tmaxlo]; %%% exclude 
    end 
    Total = Total + tmaxlo;
    if tmaxlo > tmaxlomax
        tmaxlomax = tmaxlo;
    end
end

ts=log(tmaxlomax)/log(10)/100;
time=10.^(-1:ts:log(tmaxlomax)/log(10)); %%%% time in seconds exponential sampling
%time=0:0.1:tmaxlomax; %%%% time in seconds

censored=[zeros(size(wt));ones(size(wtc))]; %%%% censored long movie


if visualize
%%%%% analyse noise (difference between data and prediction) in short movie
h=figure(3);
noise = DataExp - DataPred;
no=size(noise);
noise1= noise(:);
data1= DataPred(:);
%subplot(3,1,1)
%hist(noise1,50)
subplot(2,1,1)
plot(data1(1:50:no(1)*no(2)),noise1(1:50:no(1)*no(2)),'.','markersize',4)

xlabel('Data','fontsize',fsz)
ylabel('Residuals','fontsize',fsz)
set(gca,'fontsize',fsz)
    figfile=[dirwrite,'/residuals.pdf'];
    print(h,'-dpdf',figfile)

M = max(data1);
N = 10;
x=[0];y=[0];
for i=1:N
    ind=find((i-1)*M/N <= data1 & data1 < i*M/N);
    x=[x,mean(data1(ind))];
    y=[y,var(noise(ind))];
end

subplot(2,1,2)
hold off
plot(x,y,'xk','markersize',12)
hold on
p = polyfit(x,y,3);
plot(x,polyval(p,x),'r')

xlabel('mean data','fontsize',fsz)
ylabel('Variance residuals','fontsize',fsz)
set(gca,'fontsize',fsz)
title([num2str(p(1),2),'x^3 + ',num2str(p(2),2),'x^2 + ',num2str(p(3),2),'x + ',num2str(p(4),2)],'fontsize',fsz) 
    figfile=[dirwrite,'/residuals.pdf'];
    print(h,'-dpdf',figfile)
end %%%% if visualize
    
        
name= [name,num2str(100*alpha)] ;   
    
 

for cens=0:1
    
    
    
Ominmin=1e6;    
Kstore=[]; %%% will store parameter and objective function values
    
    
%%%% short movie    
if cens
    [fsg,xsg,los,ups]=ecdf(dtg,'censoring',censored_short);
else
    [fsg,xsg,los,ups]=ecdf(dtg);
end



for nnn=0:2:6
    
wwt=[wt;wtc]+3*(nnn)*60/2; %%%% shifted distribution from 0 t0 6 min 
 shift=3*nnn/2; %%% shift in minutes
 
 
%%%% estimate of integral for computing p1 
a=3*60;
Tinactive =sum(wwt);
Tactive = Total-Tinactive;
Nactive = Ninactive - n(2); %%%% number of active periods
imax=max(find(xsg<a));
imin=min(find(xsg>a));
x1=xsg(imin);f1=fsg(imin);
x2=xsg(imax);f2=fsg(imax);
f = f1 + (f2-f1)*(a-x1)/(x2-x1);
esp=(-a*(1-f)+trapz(xsg(xsg<a),1-fsg(xsg<a)))/f;
estp1= Ninactive/(Ninactive + Tactive/esp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
 

%%%%% long movie
if cens
[fl,xl,lol,upl]=ecdf(wwt,'censoring',censored);
else
[fl,xl,lol,upl]=ecdf(wwt);
end





%%%% fit p2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fact1=1;fact2=1;
[uxsg,iu]=unique(xsg);
ufsg=fsg(iu);
[uxl,iu]=unique(xl);
ufl=fl(iu);
M=max(uxsg)/fact2;
m=min(uxl(2:end))*fact1;
%%%%% interpolate on [m,M];
nsteps=100;
step = (M-m)/nsteps;
x=m:step:M;
y1=interp1(uxsg,ufsg,x,'cubic'); %%% fast 
y2=interp1(uxl,ufl,x,'cubic'); %%% slow
%%%%%%%%%%%%%%%%%%%%%%%%%%
p1=estp1;
opts = optimoptions(@lsqnonlin,'TolFun', 1e-15,'MaxIter',1e3, ...
        'MaxFunEvals',1e6,'TolX', 1e-10,'Display','off');
objective = @(k)(  log((1-y2)*p1) - log( (1-y1)*(1-k(1)) + k(1) ) );
     k0=0.001;    
     [k, obj] = lsqnonlin(objective,k0,0,1,opts);
p2=k;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if visualize
h=figure(200+nnn);
hold off
loglog(xsg,(1-fsg)*(1-p2)+p2,'og')
hold on
loglog(xl,(1-fl)*p1,'ok')
legend({'short movie','long movie'})
%loglog(time,s,'k')
axis([1e-1, 1e4, 1e-6, 1e0])
xlabel('Time [s]','fontsize',fsz)
ylabel('Frequency','fontsize',fsz)
set(gca,'fontsize',fsz)
title(['\Delta_0 =',num2str(3*nnn)],'fontsize',fsz)
    figfile=[dirwrite,'/figure2_shift',num2str(nnn),'_cens',num2str(cens),'.pdf'];
    print(h,'-dpdf',figfile)
end% if visualize

%%%% fit distribution of spacings using combination of 3 exponentials. 5 
%%%% params
xl=xl(1:end-1);
fl=fl(1:end-1);
lol=lol(1:end-1);
upl=upl(1:end-1);




NS = length(xsg); NL = length(xl);
sNS=sqrt(NS);sNL=sqrt(NL); 
fact1=sqrt(1-alpha);
fact2=sqrt(alpha);
 
eshort=(1-fsg)*(1-p2)+p2;
elong = (1-fl)*p1;


%%%%%%% compute average waiting time from AUC of the distribution
mean_w1=trapz([xsg;xl],[(1-fsg)*(1-p2)+p2;(1-fl)*p1]);


if alpha == 1 %%% linear 
    exp_fitness = @(k)( [(k(4)*exp(k(1)*xsg)+k(5)*exp(k(2)*xsg)+(1-k(4)-k(5))*exp(k(3)*xsg)-eshort)/sNS;...
                    (k(4)*exp(k(1)*xl)+k(5)*exp(k(2)*xl)+(1-k(4)-k(5))*exp(k(3)*xl)-elong)/sNL]); %%%%% mixed
else
    exp_fitness = @(k)( [log((k(4)*exp(k(1)*xsg)+k(5)*exp(k(2)*xsg)+(1-k(4)-k(5))*exp(k(3)*xsg))./eshort)/sNS*fact1;...
                     log((k(4)*exp(k(1)*xl)+k(5)*exp(k(2)*xl)+(1-k(4)-k(5))*exp(k(3)*xl))./elong)/sNL*fact1;...
                    (k(4)*exp(k(1)*xsg)+k(5)*exp(k(2)*xsg)+(1-k(4)-k(5))*exp(k(3)*xsg)-eshort)/sNS*fact2;...
                    (k(4)*exp(k(1)*xl)+k(5)*exp(k(2)*xl)+(1-k(4)-k(5))*exp(k(3)*xl)-elong)/sNL*fact2]); %%%%% mixed

end

opts = optimoptions(@lsqnonlin,'TolFun', 1e-8,'MaxIter',1e3, 'MaxFunEvals',1e6,'TolX', 1e-10,'Display','off');

%%%%%%% initial guess    
k00=[-0.1,-0.01,-0.001,0.25,0.25]; 

amp = [log(100),log(100),log(100)]; 




NbIterationinFit=100;


Omin=1e6;
%dofit=1;
%if dofit

disp('3 exp no constraint fit')
    for mc = 1:NbIterationinFit
        
       
        %%%% first try
        %%%% Change k00, preserve order lambda1 < lambda2 < lambda3
        factor=exp(amp.*(2*rand(1,3)-1)); 
        k0(1:3) = k00(1:3).*factor;  %%%% lambda_i
        if ~(k0(1) < k0(2) && k0(2) < k0(3))
           while ~(k0(1) < k0(2) && k0(2) < k0(3))
                factor=exp(amp.*(2*rand(1,3)-1)); 
                k0(1:3) = k00(1:3).*factor;  %%%% lambda_i 
           end
        end
        k0(4:5)=2*rand(1,2)-1;  %%%% A1, A2   
        if ~constraint0(k0)
           while ~constraint0(k0)
               k0(4:5)=2*rand(1,2)-1;  %%%% try until condition satisfied
           end
        end
        
        
        e0 = exp_fitness(k0);
        if all( [ isfinite(e0) ;  abs(imag(e0)) == 0] )
        %%%% Use the fcn lsqnonlin
            [k, obj] = lsqnonlin(exp_fitness,k0,[],[],opts);
        %%%% 
        if mod(mc,10)==0
            disp(mc)
        end
        if obj < Omin && sum(imag(k))==0 && constraint0(k)
           Omin=obj; %%%% optimal objective function
           kmin=k;   %%%% optimal parameters
           shift=3*nnn/2; %%%% shift in minutes
           p1min=p1; %%%% p1 value
           p2min=p2; %%%% p2 value
        end 
        if sum(imag(k))==0 && constraint0(k)
            Kstore=[Kstore;[k,obj]]; %%%% consider suboptimal values
        end
        end
     
    end %%%%% for mc

    
if Omin < Ominmin
           Ominmin=Omin; %%%% optimal optimum
           kminmin=kmin; %%%% optimal optimal parameters
           shiftmin=shift;
           p1min=p1;
           p2min=p2;
           fsgmin=fsg; losmin=los;upsmin=ups;
           flmin=fl; lolmin=lol;uplmin=upl;
           xlmin=xl;
           xsgmin=xsg;
           mw1opt=mean_w1;
           censmin=cens;
end     
    

OBJ = Omin ;    
   
if visualize
h=figure(105+nnn);
hold off
loglog(xsg,(1-fsg)*(1-p2)+p2,'og')
hold on
loglog(xl,(1-fl)*p1,'xk')
sfit = kmin(4)*exp(kmin(1)*time)+kmin(5)*exp(kmin(2)*time)+(1-kmin(4)-kmin(5))*exp(kmin(3)*time);
loglog(time,sfit,'r'); %%%% fitted 
xlabel('Time [s]','Fontsize',fsz)
ylabel('Freq','Fontsize',fsz) 
legend({'short movie','long movie','fitted'},'Location','southwest')
end %% if visualize

l1=kmin(1);
l2=kmin(2);
l3=kmin(3);
A1=kmin(4);
A2=kmin(5);
A3=1-A1-A2;
mean_w2=-A1/l1-A2/l2-A3/l3; %%%% average waiting time from parameters
L1=l1+l2+l3;
L2=l1*l2+l1*l3+l2*l3;
L3=l1*l2*l3;
S1=A1*l1+A2*l2+A3*l3;
S2=A1*l1^2+A2*l2^2+A3*l3^2;
S3=A1*l1^3+A2*l2^3+A3*l3^3;
k1=-L3*(S1^2-S2)/(S2^2-S1*S3);
k2=-(S2^2-S1*S3)/S1/(S1^2-S2);
k3=-S1;
k4=(S1^2-S2)/S1;
k5=-A1*A2*A3*(l1-l2)^2*(l1-l3)^2*(l2-l3)^2*S1/(S1^2-S2)/(S2^2-S1*S3);
K3= -S1;
K1p = 1/2 * ( -L1+S2/S1 + sqrt((S1*L1-S2)^2-4*L3*S1)/S1 );
K2p = 1/2 * ( -L1+S2/S1 - sqrt((S1*L1-S2)^2-4*L3*S1)/S1 );
K1m = 1/2 * (S1-S2/S1 - (-S1^2*L1+S1*S2+S1*L2-L3+S2^2/S1-S3)/sqrt((S1*L1-S2)^2-4*L3*S1));
K2m = 1/2 * (S1-S2/S1 + (-S1^2*L1+S1*S2+S1*L2-L3+S2^2/S1-S3)/sqrt((S1*L1-S2)^2-4*L3*S1));


% text(0.2,2e-1,['k1=',num2str(k1,2)])
% text(0.2,2e-2,['k2=',num2str(k2,2)])
% text(0.2,2e-3,['k3=',num2str(k3,2)])
% text(0.2,2e-4,['k4=',num2str(k4,2)])
% text(0.2,2e-5,['k5=',num2str(k5,2)])

if visualize
%text(0.2,2e-4,['mean_w1=',num2str(mean_w1,2)])
%text(0.2,2e-5,['mean_w2=',num2str(mean_w2,2),' mRNA=',num2str(45/mean_w2*60,2)])



title(['cens=',num2str(cens),'\Delta_0 =',num2str(shift),'Obj=',num2str(Omin)],'fontsize',fsz);
set(gca,'fontsize',fsz);
axis([1e-1,1e5, 1e-6, 1])
    figfile=[dirwrite,'/figure3_shift',num2str(nnn),'_cens',num2str(cens),'.pdf'];
    print(h,'-dpdf',figfile)
end %%%% if visualize   
%line=line+1;

% xlswrite(xlsfilename,{name_long,num2str(cens),num2str(nnn*3),num2str(OBJ),num2str(NS),num2str(NL),num2str(k3,2),num2str(k2,2),num2str(k4,2),...
%     num2str(k5,2),num2str(k1,2),num2str(K3,2),num2str(K1p,2),num2str(K2p,2),num2str(K1m,2),num2str(K2m,2),...
%     num2str(l1,3),num2str(l2,3),num2str(l3,3),num2str(A1,3),num2str(A2,3),num2str(A3,3)},1,['A',num2str(line)]);
    
    
%end %%%% dofit 
    
end %%% for nnn




disp('3 exp constrained fit')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% constrained 3-exp fit %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xsg=xsgmin;
xl=xlmin;
fsg=fsgmin;
fl=flmin;
NS = length(xsg); NL = length(xl);
sNS=sqrt(NS);sNL=sqrt(NL); 
eshort=(1-fsg)*(1-p2)+p2;
elong = (1-fl)*p1;
fact1=sqrt(1-alpha);
fact2=sqrt(alpha);
if alpha==1
    exp_fitness = @(k)( [(k(4)*exp(k(1)*xsg)-(k(3)+k(4)*(k(1)-k(3)))/(k(2)-k(3))*exp(k(2)*xsg)+(k(2)+k(4)*(k(1)-k(2)))/(k(2)-k(3)) *exp(k(3)*xsg)-eshort)/sNS;...
                    (k(4)*exp(k(1)*xl)-(k(3)+k(4)*(k(1)-k(3)))/(k(2)-k(3))*exp(k(2)*xl)+(k(2)+k(4)*(k(1)-k(2)))/(k(2)-k(3)) *exp(k(3)*xl)-elong)/sNL]); %%%%% mixed
else
exp_fitness = @(k)( [log((k(4)*exp(k(1)*xsg)-(k(3)+k(4)*(k(1)-k(3)))/(k(2)-k(3))*exp(k(2)*xsg)+(k(2)+k(4)*(k(1)-k(2)))/(k(2)-k(3)) *exp(k(3)*xsg))./eshort)/sNS*fact1;...
                     log((k(4)*exp(k(1)*xl)-(k(3)+k(4)*(k(1)-k(3)))/(k(2)-k(3))*exp(k(2)*xl)+(k(2)+k(4)*(k(1)-k(2)))/(k(2)-k(3)) *exp(k(3)*xl))./elong)/sNL*fact1;...
                    (k(4)*exp(k(1)*xsg)-(k(3)+k(4)*(k(1)-k(3)))/(k(2)-k(3))*exp(k(2)*xsg)+(k(2)+k(4)*(k(1)-k(2)))/(k(2)-k(3)) *exp(k(3)*xsg)-eshort)/sNS*fact2;...
                    (k(4)*exp(k(1)*xl)-(k(3)+k(4)*(k(1)-k(3)))/(k(2)-k(3))*exp(k(2)*xl)+(k(2)+k(4)*(k(1)-k(2)))/(k(2)-k(3)) *exp(k(3)*xl)-elong)/sNL*fact2]); %%%%% mixed
end
opts = optimoptions(@lsqnonlin,'TolFun', 1e-8,'MaxIter',1e3, 'MaxFunEvals',1e6,'TolX', 1e-10,'Display','off');

%%%%%%% initial guess    
k00=[-0.1,-0.01,-0.001,0.25]; 
k0=k00;
amp = [log(10),log(10),log(10)]; 
NbIterationinFit=100;    
O3min=1e6;
    for mc = 1:NbIterationinFit
%%%% first try
        %%%% Change k00, preserve order lambda1 < lambda2 < lambda3
        factor=exp(amp.*(2*rand(1,3)-1)); 
        k0(1:3) = k00(1:3).*factor;  %%%% lambda_i
        if ~(k0(1) < k0(2) && k0(2) < k0(3))
           while ~(k0(1) < k0(2) && k0(2) < k0(3))
                factor=exp(amp.*(2*rand(1,3)-1)); 
                k0(1:3) = k00(1:3).*factor;  %%%% lambda_i 
           end
        end
        k0(4)=2*rand-1;  %%%% A1   
        if ~( k0(2) + k0(4)*(k0(1)-k0(2)) < 0 )
           while ~( k0(2) + k0(4)*(k0(1)-k0(2)) < 0 )
               k0(4)=2*rand-1;  %%%% try until condition satisfied
           end
        end
        
        e0 = exp_fitness(k0);
        if all( [ isfinite(e0) ;  abs(imag(e0)) == 0] )
        %%%% Use the fcn lsqnonlin
         
            [k, obj] = lsqnonlin(exp_fitness,k0,[],[],opts);

        %%%% 
        if mod(mc,10)==0
        disp(mc)
        end
        if obj < O3min && sum(imag(k))==0
           O3min=obj; %%%% optimal objective function
           kmin=k;   %%%% optimal parameters
        end 
        end   
    end %%%%% for mc
 
sfit3expconstrained = kmin(4)*exp(kmin(1)*time)-(kmin(3)+kmin(4)*(kmin(1)-kmin(3)))/(kmin(2)-kmin(3))*exp(kmin(2)*time)+(kmin(2)+kmin(4)*(kmin(1)-kmin(2)))/(kmin(2)-kmin(3))*exp(kmin(3)*time);

disp('2 exp  fit')
%%%%%%%%%%%%%%%%% 2-exp fit%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xsg=xsgmin;
xl=xlmin;
fsg=fsgmin;
fl=flmin;

NS = length(xsg); NL = length(xl);
sNS=sqrt(NS);sNL=sqrt(NL); 
fact1=sqrt(1-alpha);
fact2=sqrt(alpha); 
eshort=(1-fsg)*(1-p2)+p2;
elong = (1-fl)*p1;

if alpha == 1
    exp_fitness = @(k)( [  (k(3)*exp(k(1)*xsg)+(1-k(3))*exp(k(2)*xsg)-eshort)/sNS;...
                    (k(3)*exp(k(1)*xl)+(1-k(3))*exp(k(2)*xl)-elong)/sNL]); %%%%% linear
else
exp_fitness = @(k)( [log((k(3)*exp(k(1)*xsg)+(1-k(3))*exp(k(2)*xsg))./eshort)/sNS*fact1;...
                     log((k(3)*exp(k(1)*xl)+(1-k(3))*exp(k(2)*xl))./elong)/sNL*fact1;...
                    (k(3)*exp(k(1)*xsg)+(1-k(3))*exp(k(2)*xsg)-eshort)/sNS*fact2;...
                    (k(3)*exp(k(1)*xl)+(1-k(3))*exp(k(2)*xl)-elong)/sNL*fact2]); %%%%% mixed
end
opts = optimoptions(@lsqnonlin,'TolFun', 1e-8,'MaxIter',1e3, 'MaxFunEvals',1e6,'TolX', 1e-10);

%%%%%%% initial guess    
k00=[-0.1,-0.01,0.25]; 
k01=[-0.1,-0.001,0.25];
k0=k00;
amp = [log(100),log(100)]; 

NbIterationinFit=100;    

O2min=1e6;
    for mc = 1:NbIterationinFit
        
        %%%% first try
        %%%% Change k00
        factor=exp(amp.*(2*rand(1,2)-1)); 
        k0(1:2) = k00(1:2).*factor;
        %%%% sort %%%%        
        k0(1:2)=sort(k0(1:2),'ascend');
        
        k0min = k0(2)/(k0(2)-k0(1));
        k0max = 1;
        k0(3)=k0min + rand*(k0max-k0min);
        
        e0 = exp_fitness(k0);
        if all( [ isfinite(e0) ;  abs(imag(e0)) == 0] )
        %%%% Use the fcn lsqnonlin
            [k, obj] = lsqnonlin(exp_fitness,k0,[],[],opts);
        %%%% 
        if mod(mc,10)==0
        disp(mc)
        end
        if obj < O2min && sum(imag(k))==0
           O2min=obj; %%%% optimal objective function
           kmin=k;   %%%% optimal parameters
        end 
        end
      
    end %%%%% for mc

     
sfit2exp = kmin(3)*exp(kmin(1)*time)+(1-kmin(3))*exp(kmin(2)*time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if visualize
    
    
%%%%%%%%%%%%%%%% visualize optimal optimal fit 
h=figure(1001);
hold off
plot(xsgmin/60,(1-fsgmin),'or','markersize',msz,'linewidth',lw)  %%%% short movie
hold on
plot(xsgmin/60,(1-fsgmin),'r','linewidth',lw) %%%% long movie
xlabel('\Delta t [min]','Fontsize',fsz)
ylabel('Survival function short','Fontsize',fsz) 
set(gca,'Fontsize',fsz)
axis([0,600/60, 1e-4, 1])
    figfile=[dirwrite,'/figure1_optimal_',name,'_cens',num2str(cens),'.pdf'];
    print(h,'-dpdf',figfile)


h=figure(1002);
hold off
plot(xlmin/3600,(1-flmin),'og','markersize',msz,'linewidth',lw) %%%% long movie
hold on
plot(xlmin/3600,(1-flmin),'g','linewidth',lw) %%%% long movie
xlabel('\Delta t [h]','Fontsize',fsz)
ylabel('Survival function long','Fontsize',fsz) 
set(gca,'Fontsize',fsz)
axis([100/3600,1.2, 1e-4, 1])
    figfile=[dirwrite,'/figure2_optimal_',name,'_cens',num2str(cens),'.pdf'];
    print(h,'-dpdf',figfile)

end %%% if visualize







if visualize
%%%%%%%%%%%%%%%% visualize optimal optimal fit 
h=figure(1000);
hold off
loglog(xsgmin,(1-fsgmin)*(1-p2min)+p2min,'or')  %%%% short movie
hold on
loglog(xlmin,(1-flmin)*p1min,'og') %%%% long movie

kmin=kminmin;
sfit = kmin(4)*exp(kmin(1)*time)+kmin(5)*exp(kmin(2)*time)+(1-kmin(4)-kmin(5))*exp(kmin(3)*time);
loglog(time,sfit,'k','linewidth',2); %%%% fitted 




ind1=1:2:length(time);

ind2=1+1:2:length(time-1);

loglog(time(ind1),sfit3expconstrained(ind1),'k+','markersize',msz); %%%% fitted 
loglog(time(ind2),sfit2exp(ind2),'kx','markersize',msz); %%%% fitted 

%loglog(xsgmin,(1-losmin)*(1-p2min)+p2min,'--r')
%loglog(xsgmin,(1-upsmin)*(1-p2min)+p2min,'--r')
%loglog(xlmin,(1-lolmin)*p1min,'--g')
%loglog(xlmin,(1-uplmin)*p1min,'--g')

xlabel('\Delta t [s]','Fontsize',fsz)
ylabel('Survival function','Fontsize',fsz) 
%legend({'short movie','long movie','fitted 3exp','fitted 3exp constrained','fitted 2exp'},'Fontsize',fsz,'Location','southwest')

set(gca,'Fontsize',fsz)
%title(['Obj=',num2str(Ominmin)],'fontsize',fsz);
axis([1e-1,1e5, 1e-6, 1])



set(gca,'fontsize',fsz);
    figfile=[dirwrite,'/figure3_optimal_',name,'_cens',num2str(cens),'.pdf'];
    print(h,'-dpdf',figfile)


%%%%% comparison linear scale 
h=figure(1003);
hold off
plot(xsgmin/60,(1-fsgmin)*(1-p2min)+p2min,'or')  %%%% short movie
hold on
plot(time/60,sfit,'k','linewidth',2); %%%% fitted 
plot(time(ind1)/60,sfit3expconstrained(ind1),'k+','markersize',msz); %%%% fitted 
plot(time(ind2)/60,sfit2exp(ind2),'kx','markersize',msz); %%%% fitted 
xlabel('\Delta t [min]','Fontsize',fsz)
ylabel('Survival function','Fontsize',fsz) 
%legend({'short movie','fitted 3exp','fitted 3exp constrained','fitted 2exp'},'Fontsize',fsz,'Location','northeast')
set(gca,'Fontsize',fsz)
axis([0, 3,0, 1])


set(gca,'fontsize',fsz);
    figfile=[dirwrite,'/figure4_optimal_',name,'_cens',num2str(cens),'.pdf'];
    print(h,'-dpdf',figfile)


fname=[dirwrite,'/figure3_optimal_',name,'_cens',num2str(cens),'.mat'];
save(fname,'xsgmin','fsgmin','xlmin','flmin','time','sfit','sfit3expconstrained','sfit2exp','ind1','ind2');

end %%% if visualize




%%%%%%%%% optimal parameters
l1=kmin(1);
l2=kmin(2);
l3=kmin(3);
A1=kmin(4);
A2=kmin(5);
A3=1-A1-A2;
mean_w2=-A1/l1-A2/l2-A3/l3; %%%% average waiting time from parameters
mean_w1=mw1opt;
L1=l1+l2+l3;
L2=l1.*l2+l1.*l3+l2.*l3;
L3=l1.*l2.*l3;
S1=A1.*l1+A2.*l2+A3.*l3;
S2=A1.*l1.^2+A2.*l2.^2+A3.*l3.^2;
S3=A1.*l1.^3+A2.*l2.^3+A3.*l3.^3;
k1min=-L3.*(S1.^2-S2)./(S2.^2-S1.*S3);
k2min=-(S2.^2-S1.*S3)./S1./(S1.^2-S2);
k3min=-S1;
k4min=(S1.^2-S2)./S1;
k5min=-A1.*A2.*A3.*(l1-l2).^2.*(l1-l3).^2.*(l2-l3).^2.*S1./(S1.^2-S2)./(S2.^2-S1.*S3);
K3min= -S1;
K1pmin = 1/2 * ( -L1+S2/S1 + sqrt((S1*L1-S2)^2-4*L3*S1)/S1 );
K2pmin = 1/2 * ( -L1+S2/S1 - sqrt((S1*L1-S2)^2-4*L3*S1)/S1 );
K1mmin = 1/2 * (S1-S2/S1 - (-S1^2*L1+S1*S2+S1*L2-L3+S2^2/S1-S3)/sqrt((S1*L1-S2)^2-4*L3*S1));
K2mmin = 1/2 * (S1-S2/S1 + (-S1^2*L1+S1*S2+S1*L2-L3+S2^2/S1-S3)/sqrt((S1*L1-S2)^2-4*L3*S1));
l1min=l1;l2min=l2;l3min=l3;A1min=A1;A2min=A2;A3min=A3;


%%%%%%%%%%%%%%%% compute parameters with error bars
overflow=1; %%%%% 100% overflow
%%%% select near-optimal parameters
%%%% O between Ominmin and Ominmin*(1+overflow)  %%%  
Ksel = Kstore( Kstore(:,6) < Ominmin*(1+overflow) , :  );
%%%% compute intervals for parameters
l1=Ksel(:,1);
l2=Ksel(:,2);
l3=Ksel(:,3);
A1=Ksel(:,4);
A2=Ksel(:,5);
A3=1-A1-A2;


MRNA = 45./(-A1./l1-A2./l2-A3./l3)*60;



L1=l1+l2+l3;
L2=l1.*l2+l1.*l3+l2.*l3;
L3=l1.*l2.*l3;
S1=A1.*l1+A2.*l2+A3.*l3;
S2=A1.*l1.^2+A2.*l2.^2+A3.*l3.^2;
S3=A1.*l1.^3+A2.*l2.^3+A3.*l3.^3;
%%%% model M1
K1=-L3.*(S1.^2-S2)./(S2.^2-S1.*S3); %%% k1p
K2=-(S2.^2-S1.*S3)./S1./(S1.^2-S2); %%% k2p
K3=-S1; %%%% k3
K4=(S1.^2-S2)./S1; %%% k2m
K5=-A1.*A2.*A3.*(l1-l2).^2.*(l1-l3).^2.*(l2-l3).^2.*S1./(S1.^2-S2)./(S2.^2-S1.*S3); %%% k1m

K3p = -S1;
K1p = 1/2 * ( -L1+S2./S1 + sqrt((S1.*L1-S2).^2-4*L3.*S1)./S1 );
K2p = 1/2 * ( -L1+S2./S1 - sqrt((S1.*L1-S2).^2-4*L3.*S1)./S1 );
K1m = 1/2 * (S1-S2./S1 - (-S1.^2.*L1+S1.*S2+S1.*L2-L3+S2.^2./S1-S3)./sqrt((S1.*L1-S2).^2-4*L3.*S1));
K2m = 1/2 * (S1-S2./S1 + (-S1.^2.*L1+S1.*S2+S1.*L2-L3+S2.^2./S1-S3)./sqrt((S1.*L1-S2).^2-4*L3.*S1));

line=line+1;

 xlswrite(xlsfilename,{name,num2str(Ominmin),num2str(O3min),num2str(O2min),num2str(NS),num2str(NL),num2str(shiftmin),num2str(cens),num2str(k3min,2),num2str(k2min,2),num2str(k4min,2),...
     num2str(k5min,2),num2str(k1min,2),num2str(K3min,2),num2str(K1pmin,2),num2str(K2pmin,2),num2str(K1mmin,2),num2str(K2mmin,2),...
    num2str(l1min,3),num2str(l2min,3),num2str(l3min,3),num2str(A1min,3),num2str(A2min,3),num2str(A3min,3),45/mean_w2*60},1,['A',num2str(line)]);
 
line = line+1;

 xlswrite(xlsfilename,{name,'','','','','','','',num2str(min(K3),2),num2str(min(K2),2),num2str(min(K4),2),...
     num2str(min(K5),2),num2str(min(K1),2),num2str(min(K3p),2),num2str(min(K1p),2),num2str(min(K2p),2),num2str(min(K1m),2),num2str(min(K2m),2),...
    num2str(min(l1),3),num2str(min(l2),3),num2str(min(l3),3),num2str(min(A1),3),num2str(min(A2),3),num2str(min(A3),3),num2str(min(MRNA),3)},1,['A',num2str(line)]);

line = line+1;

 xlswrite(xlsfilename,{name,'','','','','','','',num2str(max(K3),2),num2str(max(K2),2),num2str(max(K4),2),...
     num2str(max(K5),2),num2str(max(K1),2),num2str(max(K3p),2),num2str(max(K1p),2),num2str(max(K2p),2),num2str(max(K1m),2),num2str(max(K2m),2),...
    num2str(max(l1),3),num2str(max(l2),3),num2str(max(l3),3),num2str(max(A1),3),num2str(max(A2),3),num2str(max(A3),3),num2str(max(MRNA),3)},1,['A',num2str(line)]);

end %%% cens


end %%%% iname


%fclose(fid);



end %%% ibig

