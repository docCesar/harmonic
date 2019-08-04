clear;close all;clc

% Next update
% 0. (IMPORTANT!!!) TO solve the problem that the number of lines is
% different.
% 1. Add 3D samples cases.
% 2. Seperate Hall data and WAL data automatically.

% Constants
h=6.6260699*10^-34;
hbar=1.0545718e-34;
m=9.10938356*10^(-31);
e=1.60217662*10^-19;

% Parameters of samples
% asp=4;          % Old mask
% asp=2.5;        % Small devices of new mask
asp=24;         % Large devices of new mask
thick=5*10^-9; % Thickness of sample layer

%% Import data
fileName=dir(fullfile('*.txt'));
%fileName= dir(fullfile('d:/datafile','*.txt'));             
%The path should be changed to what you want.
number=length(fileName);
% Claim variables.
temForCheck=[];
thicknessForCheck=[];
hall=[];
bH=[];
rH=[];
gHall=[];
bWAL=[];
rWAL=[];
gWAL=[];
mu=[];
taup=[];
vf=[];
dif=[];
be=[];
ne=[];
kf=[];

bi=[];
bso=[];

for i=1:number
    name=convertCharsToStrings(fileName(i).name);
    % Get the temperature for check.
    patternTem='(?<=_)\w*(?=K_)';
    temRegTem=regexp(name,patternTem,'match');
    temForCheck(end+1)=str2num(temRegTem{1});
    
    % Get the thickness for check.
    patternThickness='(?<=K_)\w*(?=nm_Vg)';
    temRegThickness=regexp(name,patternThickness,'match');
    thicknessForCheck(end+1)=str2num(temRegThickness{1});
    %% Initialization

    delimiter = '\t';
    startRow = 2;
    
    %% Format
    formatSpec = '%f%f%*s%[^\n\r]';
    
    %% Open the data.
    fileID = fopen(name,'r');

    %% Read the data column.

    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    [l,c]=size(dataArray{1});
    
	%% Close the data.
    fclose(fileID);
    
    %% Name the data column.

%     
%     if contains(name,'Rxy_')==1    % 稍后应添加一个匹配判断
%         bH(:,end+1) = dataArray{:, 1}./1000;
%         rH(:,end+1) = dataArray{:, 2};
%         gHall(:,end+1)=asp*h/(e^2)./dataArray{:, 2};
%         bWAL(1:l,end+1)=0;
%         rWAL(1:l,end+1)=0;
%         gWAL(1:l,end+1)=0;
%     elseif contains(name,'Rxx_')==1
%         bWAL(:,end+1) = dataArray{:, 1}./1000;
%         rWAL(:,end+1) = dataArray{:, 2};
%         gWAL(:,end+1)=asp*h/(e^2)./dataArray{:, 2};
%         bH(1:l,end+1)=0;
%         rH(1:l,end+1)=0;
%         gHall(1:l,end+1)=0;
%     else
%         "Error code 01 (Wrong data name)"
%         return
%     end

    
    if contains(name,'Rxy_')==1    % 稍后应添加一个匹配判断
        'a'
        eval(['bH_' , temRegThickness{1} , 'nm = dataArray{:, 1}./1000;'])
        eval(['rH_' , temRegThickness{1} , 'nm = dataArray{:, 2};'])
        eval(['gHall_' , temRegThickness{1} , 'nm = asp*h/(e^2)./dataArray{:, 2};'])
    elseif contains(name,'Rxx_')==1
        'b'
        eval(['bWAL_' , temRegThickness{1} , 'nm = dataArray{:, 1}./1000;'])
        eval(['rWAL_' , temRegThickness{1} , 'nm = dataArray{:, 2};'])
        eval(['gWAL_' , temRegThickness{1} , 'nm =asp*h/(e^2)./dataArray{:, 2};'])
    else
        "Error code 01 (Wrong data name)"
        return
    end

    %% Clean the temporary variation.
    clearvars temRegTem patternTem temRegThickness patternThickness filename delimiter startRow formatSpec fileID dataArray ans;
end


%% rHall vs B
figure
for i=number/2+1:number
   eval(['bH=bH_', num2str(thicknessForCheck(i),'%.2d') , 'nm;'])
   eval(['rH=rH_', num2str(thicknessForCheck(i),'%.2d') , 'nm;'])
   plot(bH,rH,'Linewidth',2)
   xlabel {B (T)}
   ylabel {R_{xy} (\Omega)}
   grid on
   hold on    
end
% temForLegend=transpose(temForCheck(1,1:number/2));
% temForLegend=num2str(temForLegend);
thicknessForLegend=transpose(thicknessForCheck(1,1:number/2));
thicknessForLegend=num2str(thicknessForLegend);
title('Hall resistance')
legend(thicknessForLegend)

%% Hall fitting



for i=1:number/2
    eval(['bH=bH_', num2str(thicknessForCheck(i),'%.2d') , 'nm;'])
    eval(['rH=rH_', num2str(thicknessForCheck(i),'%.2d') , 'nm;'])
    eval(['gWAL=abs(gWAL_', num2str(thicknessForCheck(i),'%.2d') , 'nm);'])
    [xData, yData] = prepareCurveData( bH, rH );

    % Set up fittype and options.
    ft = fittype( 'p1*x+p2', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Robust = 'Bisquare';
    opts.StartPoint = [0.9134 0.0975404049994095];
    opts.TolFun = 1e-08;
    opts.TolX = 1e-08;

    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );

    % Extract fitting parameters.
    hall(end+1)=abs(fitresult.p1);
%     'hall'
%     fitresult.p1
    
    % Plot fit with data.
    figure( 'Name', 'untitled fit 1' );
    plotH = plot( fitresult, xData, yData );
    legend( plotH, 'rH vs. bH', 'untitled fit 1', 'Location', 'NorthEast' );
    % Label axes
    xlabel bH
    ylabel rH
    grid on

% Properties from Hall fitting.
ne(i)=abs(1/e/hall(i));
mu(i)=e^2/h*max(gWAL)/e/ne(i);
taup(i)=m/e*mu(i);
kf(i)=(2*pi*ne(i))^(1/2);  % 2D case.
%kf(i)=(3*pi^2*ne)^(1/3);  % 3D case.
vf(i)=hbar/m*kf(i);
dif(i)=vf(i)^2*taup(i)/2;
be(i)=hbar/(4*e*dif(i)*taup(i));

end
'Hall finish'
% return

%% Rxx vs B
figure
for i=1:number/2
%     figure
   eval(['bWAL=bWAL_', num2str(thicknessForCheck(i),'%.2d') , 'nm;'])
   eval(['rWAL=abs(rWAL_', num2str(thicknessForCheck(i),'%.2d') , 'nm);'])
   plot(bWAL,rWAL,'Linewidth',2)
   xlabel {B (T)}
   ylabel {R_{xx} (\Omega)}
   grid on
   hold on
    
end
title('Longitude resistance')
legend(thicknessForLegend)

%% Rxx(Normalized) vs B
figure
for i=1:number/2
%     figure
   eval(['bWAL=bWAL_', num2str(thicknessForCheck(i),'%.2d') , 'nm;'])
   eval(['rWAL=abs(rWAL_', num2str(thicknessForCheck(i),'%.2d') , 'nm);'])
   rWALNor=(rWAL-min(rWAL))./min(rWAL);
   plot(bWAL,rWALNor,'Linewidth',2)
   xlabel {B (T)}
   ylabel {R_{xx}(Normalized)}
   grid on
   hold on
    
end
title('Longitude resistance(Normalized)')
legend(thicknessForLegend)

%% Gxx(Normalized) vs B
figure
for i=1:number/2
%     figure
   eval(['bWAL=bWAL_', num2str(thicknessForCheck(i),'%.2d') , 'nm;'])
   eval(['gWAL=gWAL_', num2str(thicknessForCheck(i),'%.2d') , 'nm;'])
   gWALNor=abs((gWAL-min(gWAL))./min(gWAL));
   plot(bWAL,gWALNor,'Linewidth',2)
   xlabel {B (T)}
   ylabel {G_{xx}(Normalized)}
   grid on
   hold on
    
end
title('Longitude conductance(Normalized)')
legend(thicknessForLegend)
%% WAL fitting
% Pretreatment of data.
% bWAL=abs(bWAL);     % Get the abstract of bWAL.

% [bWALmin,bWALpos]=min(bWAL)
% bWAL=bWAL-bWALmin;
% bWAL(1800,:)=[];
% gWAL(1800,:)=[];



for i=1:number/2
    eval(['bWAL=abs(bWAL_', num2str(thicknessForCheck(i),'%.2d') , 'nm);'])
%     gWAL=asp.*h./(e^2)./rWAL(:,i);
    eval(['gWAL=abs(gWAL_', num2str(thicknessForCheck(i),'%.2d') , 'nm);'])
    dgWAL=gWAL-max(gWAL);
    [xData, yData] = prepareCurveData( bWAL, dgWAL );
    

    % Set up fittype and options.
    ft = fittype( '1/3.14159*(-psi((bso+be)/x+0.5)+log((bso+be)/x)+1.5*psi((bi+4*bso/3)/x+0.5)-1.5*log((bi+4*bso/3)/x)-0.5*psi(bi/x+0.5)+0.5*log(bi/x))+kFactor*x^2', 'independent', 'x', 'dependent', 'y' , 'problem' , 'be' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.DiffMinChange = 1e-12;
    opts.Display = 'Off'; 
    opts.Lower = [0 0];
    opts.MaxFunEvals = 3000;
    opts.MaxIter = 2000;
    opts.Robust = 'LAR';
    opts.StartPoint = [0.02 -0.05 0];
    opts.TolFun = 1e-08;
    opts.TolX = 1e-08;
    
    % Fit model to data.
    i
    [fitresult, gof] = fit( xData, yData, ft, opts , 'problem' , be(i))
    
    bso(end+1)=fitresult.bso;
    bi(end+1)=fitresult.bi;
    
    % Plot fit with data.
    figure( 'Name', 'untitled fit 1' );
    scatter( xData, yData,'o' );
    box on
    hold on
    plot( fitresult, xData, yData );
    legend( 'gWALhalf vs. bWALhalf', 'untitled fit 1', 'Location', 'NorthEast' );
    % Label axes
    xlabel bWALhalf
    ylabel gWALhalf
    grid on
    
    
    
end

lso=sqrt(hbar./(4.*e.*bso));
lphi=sqrt(hbar./(4.*e.*bi));
ltr=(hbar/e*sqrt(2*pi)).*mu;
tauso=hbar./(4*e.*bso.*dif);

%% Change the units
dif=dif.*10000;
mu=mu.*10000;


%% Ne vs thickness
figure
scatter(thicknessForCheck(1:number/2),ne,'o');
xlim([0 1.1*max(thicknessForCheck)])
% upLim=max(ne)+0.1*(max(ne)-min(ne));
% downLim=min(ne)-0.1*(max(ne)-min(ne));
% ylim([downLim upLim])
% clearvars upLim downLim

xlabel {Thickness (nm)}
ylabel {n_e (cm^{-3})}
title('N_e vs thickness')
grid on
box on

%% D vs thickness
figure
scatter(thicknessForCheck(1:number/2),dif,'o');
xlim([0 1.1*max(thicknessForCheck)])
% ylim([0.999*min(dif) 1.001*max(dif)])
xlabel {Thickness (nm)}
ylabel {D (cm^2/s)}
title('D vs thickness')
grid on
box on

%% L_SO vs thickness
figure
scatter(thicknessForCheck(1:number/2),lso,'o');
xlim([0 1.1*max(thicknessForCheck)])
% ylim([0.999*min(dif) 1.001*max(dif)])
xlabel {Thickness (nm)}
ylabel {L_{SO} (m)}
title('L_{SO} vs thickness')
grid on
box on

%% L_Phi vs thickness
figure
scatter(thicknessForCheck(1:number/2),lphi,'o');
xlim([0 1.1*max(thicknessForCheck)])
% ylim([0.999*min(dif) 1.001*max(dif)])
xlabel {Thickness (nm)}
ylabel {L_{\phi} (m)}
title('L_{\phi} vs thickness')
grid on
box on

%% L_SO vs D
figure
scatter(dif,lso,'o');
% xlim([0 1.1*max(temForCheck)])
% ylim([0.999*min(dif) 1.001*max(dif)])
xlabel {D (m^2/s)}
ylabel {L_{SO} (m)}
title('L_{SO} vs D')
grid on
box on

%% L_phi vs D
figure
scatter(dif,lphi,'o');
% xlim([0 1.1*max(temForCheck)])
% ylim([0.999*min(dif) 1.001*max(dif)])
xlabel {D (m^2/s)}
ylabel {L_{\phi} (m)}
title('L_{\phi} vs D')
grid on
box on

%% Tau_SO vs Tau_p
figure
scatter(taup,tauso,'o');
xlabel {\tau_p (s)}
ylabel {\tau_{SO} (s)}
title('\tau_{SO} vs tau_p')
grid on
box on

%% Tau_SO vs D
figure
scatter(dif,tauso,'o');
xlabel {D (cm^2/s)}
ylabel {\tau_{SO} (s)}
title('\tau_{SO} vs D')
grid on
box on

%% mu vs thickness
figure
scatter(thicknessForCheck(1:number/2),mu,'o');
xlim([0 1.1*max(thicknessForCheck)])
xlabel {Thickness (nm)}
ylabel {\mu (cm^2/V/s)}
title('\mu vs thickness')
grid on
box on

%% Tau_SO vs thickness
figure
scatter(thicknessForCheck(1:number/2),tauso,'o');
xlim([0 1.1*max(thicknessForCheck)])
% ylim([0.999*min(dif) 1.001*max(dif)])
xlabel {Thickness (nm)}
ylabel {\tau_{SO} (s)}
title('\tau_{SO} vs thickness')
grid on
box on

%% Tau_p vs thickness
figure
scatter(thicknessForCheck(1:number/2),taup,'o');
xlim([0 1.1*max(thicknessForCheck)])
% ylim([0.999*min(dif) 1.001*max(dif)])
xlabel {Thickness (nm)}
ylabel {\tau_p (s)}
title('\tau_p vs thickness')
grid on
box on

%% vf vs thickness
figure
scatter(thicknessForCheck(1:number/2),vf,'o');
xlim([0 1.1*max(thicknessForCheck)])
xlabel {Thickness (nm)}
ylabel {v_f }
title('v_f vs thickness')
grid on
box on

%% Bso vs thickness
figure
scatter(thicknessForCheck(1:number/2),bso,'o');
xlim([0 1.1*max(thicknessForCheck)])
xlabel {Thickness (nm)}
ylabel {B_{SO} (T)}
title('B_{SO} vs thickness')
grid on
box on


'Finished'