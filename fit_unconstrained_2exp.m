%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Code written by Ovidiu Radulescu, University of Montpellier, June 2019
%%%%% reads results of genetic algorithm for short movies, data from long
%%%%% movies
%%%%% needs 1) short movies decomvolution results  result_namecPosPred.mat
%%%%% 2) long movie data name_long_raw.mat
%%%%% performs unconstrained two exponentials fit, compute parameters of
%%%%% the two states model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all


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



%%%% file names %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
list_short = {{'result_tat_latency_cPosPred10.mat'}};
files_long  = {'data_tat_latency_long_raw500.mat'};
names = {'tat_latency'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ibig = 0:2
    
    alpha = 0.9 - ibig * 0.3;
    


DataFilePath0 = ['results_2exp_no_constraint_alpha',num2str(100*alpha)]; %%%% where to write results
mkdir(DataFilePath0);
delete([DataFilePath0,'/*']); %%% delete previous result files

xlsfilename = [DataFilePath0,'/results_2exp_',num2str(100*alpha),'.xlsx'];


line=1;
xlswrite(xlsfilename,{'Fname','OBJ','shift','cens','lambda1','lambda2','A1','A2','k2','k1m','k1p','mRNA'},1,'A1');

for iname=1:length(names)
line=line+1;
    

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
cPosPred=cPosPred(:,isel);
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
        dt=[dt;dtimes(1:end)];
        %%%% censured intervals
        if tmax > times(end)
            dtc=[dtc;tmax-times(end)];
        end
        if times(1) > 0
            dtc=[dtc;times(1)];
        end
    end
end
dt(dt==0)=[];%%%% eliminate zero intervals
dtg=[dt;dtc];
censored_short=[zeros(size(dt));ones(size(dtc))];



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
isel=find(Q(1)-1.5*(Q(3)-Q(1)) < FI & FI < min([Q(3) + 1.5*(Q(3)-Q(1)),1]) & FI>0);
%isel=find(Q(1)-5*(Q(3)-Q(1)) < FI & FI < 1 & FI < Q(3) + 5*(Q(3)-Q(1)));
DataExpLong=DataExpLong(:,isel);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  
   
    
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

time=0:0.1:tmaxlomax; %%%% time in seconds

censored=[zeros(size(wt));ones(size(wtc))]; %%%% censored long movie



    
    
    
Ominmin=1e6;    
Kstore=[]; %%% will store parameter and objective function values



    
%%%% short movie    
if cens
    [fsg,xsg,logg,upg]=ecdf(dtg,'censoring',censored_short);
else
    [fsg,xsg,logg,upg]=ecdf(dtg);
end





name= [name,num2str(100*alpha)] ;   

    
 


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
[fl,xl]=ecdf(wwt,'censoring',censored);
else
[fl,xl]=ecdf(wwt);
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
        'MaxFunEvals',1e6,'TolX', 1e-10);
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
    figfile=[dirwrite,'/figure2_2exp_shift',num2str(nnn),'_cens',num2str(cens),'.pdf'];
    print(h,'-dpdf',figfile)
end% if visualize 

%%%% fit distribution of spacings using combination of 3 exponentials. 5 
%%%% params
xl=xl(1:end-1);
fl=fl(1:end-1);
NS = length(xsg); NL = length(xl);
sNS=sqrt(NS);sNL=sqrt(NL); 
fact1=sqrt(1-alpha);
fact2=sqrt(alpha);
 
eshort=(1-fsg)*(1-p2)+p2;
elong = (1-fl)*p1;


%%%%%%% compute average waiting time from AUC of the distribution
mean_w1=trapz([xsg;xl],[eshort;elong]);


exp_fitness = @(k)( [log((k(3)*exp(k(1)*xsg)+(1-k(3))*exp(k(2)*xsg))./eshort)/sNS*fact1;...
                     log((k(3)*exp(k(1)*xl)+(1-k(3))*exp(k(2)*xl))./elong)/sNL*fact1;...
                    (k(3)*exp(k(1)*xsg)+(1-k(3))*exp(k(2)*xsg)-eshort)/sNS*fact2;...
                    (k(3)*exp(k(1)*xl)+(1-k(3))*exp(k(2)*xl)-elong)/sNL*fact2]); %%%%% mixed


opts = optimoptions(@lsqnonlin,'TolFun', 1e-8,'MaxIter',1e3, 'MaxFunEvals',1e6,'TolX', 1e-10);

%%%%%%% initial guess    
k00=[-0.1,-0.01,0.25]; 
k01=[-0.1,-0.001,0.25];

amp = [log(100),log(100)]; 

NbIterationinFit=100;    

Omin=1e6;
dofit=1;
if dofit
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
        disp(mc)
        if obj < Omin && sum(imag(k))==0
           Omin=obj; %%%% optimal objective function
           kmin=k;   %%%% optimal parameters
           shift=3*nnn/2; %%%% shift in minutes
           p1min=p1; %%%% p1 value
           p2min=p2; %%%% p2 value
        end 
        if sum(imag(k))==0
            Kstore=[Kstore;[k,obj]]; %%%% consider suboptimal values
        end
        end
      
    end %%%%% for mc

    
if Omin < Ominmin
           Ominmin=Omin; %%%% optimal optimum
           kminmin=kmin; %%%% optimal optimal parameters
           p1min=p1;
           p2min=p2;
           fsgmin=fsg;
           flmin=fl;
           xlmin=xl;
           xsgmin=xsg;
           mw1opt=mean_w1;
           shiftmin = shift;
end     
    

OBJ = Omin ;    
   
if visualize
h=figure(105+nnn);
hold off
loglog(xsg,(1-fsg)*(1-p2)+p2,'og')
hold on
loglog(xl,(1-fl)*p1,'xk')
sfit = kmin(3)*exp(kmin(1)*time)+(1-kmin(3))*exp(kmin(2)*time);
loglog(time,sfit,'r'); %%%% fitted 
xlabel('Time [s]','Fontsize',fsz)
ylabel('Freq','Fontsize',fsz) 
legend({'short movie','long movie','fitted'})
end %% if visualize

l1=kmin(1);
l2=kmin(2);
A1=kmin(3);
A2=1-A1;
mean_w2=-A1/l1-A2/l2; %%%% average waiting time from parameters


if visualize
text(0.2,2e-4,['mean_w1=',num2str(mean_w1,2)])
text(0.2,2e-5,['mean_w2=',num2str(mean_w2,2),' mRNA=',num2str(45/mean_w2*60,2)])



title(['cens=',num2str(cens),'\Delta_0 =',num2str(shift),'Obj=',num2str(Omin)],'fontsize',fsz);
set(gca,'fontsize',fsz);
axis([1e-1,1e5, 1e-6, 1])
    figfile=[dirwrite,'/figure3_2exp_shift',num2str(nnn),'_cens',num2str(cens),'.pdf'];
    print(h,'-dpdf',figfile)
end %%%% if visualize   

    
end %%%% dofit 
    
end %%% for nnn



if visualize
%%%%%%%%%%%%%%%% visualize optimal optimal fit 
h=figure(1000);
hold off
loglog(xsgmin,(1-fsgmin)*(1-p2min)+p2min,'or')
hold on
loglog(xlmin,(1-flmin)*p1min,'og')
kmin=kminmin;
sfit = kmin(3)*exp(kmin(1)*time)+(1-kmin(3))*exp(kmin(2)*time);
loglog(time,sfit,'k'); %%%% fitted 
xlabel('\Delta t [s]','Fontsize',fsz)
ylabel('Survival function','Fontsize',fsz) 
legend({'short movie','long movie','fit 2exp'})
set(gca,'Fontsize',fsz)
%title(['Obj=',num2str(Ominmin)],'fontsize',fsz);
axis([1e-1,1e5, 1e-6, 1])
    figfile=[dirwrite,'/figure3_optimal_2exp_',name,'_cens',num2str(cens),'.pdf'];
    print(h,'-dpdf',figfile)
    
    
    
h=figure(1003);
hold off
plot(xsgmin/60,(1-fsgmin)*(1-p2min)+p2min,'or')  %%%% short movie
hold on
plot(time/60,sfit,'k','linewidth',2); %%%% fitted 

xlabel('\Delta t [min]','Fontsize',fsz)
ylabel('Survival function','Fontsize',fsz) 
%legend({'short movie','fitted 3exp','fitted 3exp constrained','fitted 2exp'},'Fontsize',fsz,'Location','northeast')
set(gca,'Fontsize',fsz)
axis([0, 3,0, 1])

set(gca,'fontsize',fsz);
    figfile=[dirwrite,'/figure4_optimal_2exp_',name,'_cens',num2str(cens),'.pdf'];
    print(h,'-dpdf',figfile)
  
    
 fname=[dirwrite,'/figure3_optimal_2exp_',name,'_cens',num2str(cens),'.mat'];
save(fname,'xsgmin','fsgmin','xlmin','flmin','time','sfit');
   
    

end %%% if visualize

%%%%%%%%% optimal parameters
l1=kmin(1);
l2=kmin(2);
A1=kmin(3);
A2=1-A1;
mean_w2=-A1/l1-A2/l2; %%%% average waiting time from parameters
mean_w1=mw1opt;
% L1=l1+l2+l3;
% L2=l1.*l2+l1.*l3+l2.*l3;
% L3=l1.*l2.*l3;
S1=A1.*l1+A2.*l2;
S2=A1.*l1.^2+A2.*l2.^2;
S3=A1.*l1.^3+A2.*l2.^3;
k2min=-S1;
k1mmin=S1-S2/S1;
k1pmin=(S3*S1-S2^2)/S1/(S1^2-S2);

%%%%%%% sort %%%%%%%%%%%%%%%%%%%%%%
KK  = [A1,A2];
K = [l1,l2];
[K_sorted,isort] = sort(K,2,'ascend');
KK=KK(isort);
l1=K_sorted(1);
l2=K_sorted(2);
A1=KK(1);
A2=KK(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


l1min=l1;l2min=l2;A1min=A1;A2min=A2;



%%%%%%%%%%%%%%%% compute parameters with error bars
overflow=1; %%%%% 100% overflow
%%%% select near-optimal parameters
%%%% O between Ominmin and Ominmin*(1+overflow)  %%%  
Ksel = Kstore( Kstore(:,4) < Ominmin*(1+overflow) , :  );
%%%% compute intervals for parameters
l1=Ksel(:,1);
l2=Ksel(:,2);
A1=Ksel(:,3);
A2=1-A1;
S1=A1.*l1+A2.*l2;
S2=A1.*l1.^2+A2.*l2.^2;
S3=A1.*l1.^3+A2.*l2.^3;
K2=-S1;
K1m=S1-S2./S1;
K1p=(S3.*S1-S2.^2)./S1./(S1.^2-S2);



%% sort Ksel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ksel_ext=[A1,A2];
[Ksel_sorted,isort] = sort(Ksel(:,1:2),2,'ascend');
nn=size(isort);
for i1=1:nn(1)
    Ksel_ext(i1,:) = Ksel_ext(i1,isort(i1,:));
end
Ksel_ext=[Ksel_sorted,Ksel_ext];
l1=Ksel_ext(:,1);
l2=Ksel_ext(:,2);
A1=Ksel_ext(:,3);
A2=Ksel_ext(:,4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








xlswrite(xlsfilename,{name,num2str(Ominmin),num2str(shiftmin),num2str(cens),num2str(l1min,3),num2str(l2min,3),num2str(A1min,3),num2str(A2min,3),num2str(k2min,3),num2str(k1mmin,3),num2str(k1pmin,3),45/mean_w2*60},1,['A',num2str(line)]);
 
line = line+1;

xlswrite(xlsfilename,{name,'','','',num2str(min(l1),3),num2str(min(l2),3),num2str(min(A1),3),num2str(min(A2),3),num2str(min(K2),3),num2str(min(K1m),3),num2str(min(K1p),3),},1,['A',num2str(line)]);

line = line+1;

xlswrite(xlsfilename,{name,'','','',num2str(max(l1),3),num2str(max(l2),3),num2str(max(A1),3),num2str(max(A2),3),num2str(max(K2),3),num2str(max(K1m),3),num2str(max(K1p),3),},1,['A',num2str(line)]);


end %%%% iname


end %%% alpha




