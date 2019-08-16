clear;close all;clc;

fileName=dir(fullfile('*.txt'));
%fileName= dir(fullfile('d:/datafile','*.txt'));             
%The path should be changed to what you want.

number=length(fileName);
thicknessForTag=string;
magneticForTag=string;

% Parameters
rRef = 20;
bk = 1;

magneticForNum = [];
maxR2ndHarmonic = [];
rad = [];
rfl = [];
off = [];
c = [];
rphe = [];
off1 = [];
c1 = [];
kUnknow = [];

tagN = 0;
tag = [];

for i=1:number
    isDeg=~isempty(strfind(fileName(i).name,'deg'));
    if isDeg==1
        y = importfile(fileName(i).name);
        
        v0 = y(:,1);
        i0=median(v0)/rRef;
        y(:,1) = y(:,1)./i0;
        y(:,2) = y(:,2)./i0;
        har45Deg = y;
        tag=i;
        tagN=tagN+1;
    end
    clearvars y isDeg
end

fileName(tag)=[];
if tagN~=0
    number=number-tagN;
    i0=median(v0)/rRef;
else
    "Lack of the value of current. Please enter it."
    i0=input('In (A)\n')
end
clearvars tag tagN    


for i=1:number    
    % Extract the thickness and magnetic field.
    patternTh = '(?<=Py_)\w*(?=nm_)';
    patternMag = '(?<=harmonic)\-?\w*(?=\(mT\))';
    thickness=regexp(fileName(i).name,patternTh, 'match');
    thicknessForTag(end+1)=thickness{1};
    magneticField=regexp(fileName(i).name,patternMag, 'match');
    
    magneticForNum(end+1)=str2num(magneticField{1});
    magneticForTag(end+1)=magneticField{1};
    
    y = importfile(fileName(i).name);
    y(:,1)=y(:,1).*pi./180;
    y(:,2)=y(:,2)./i0;
    y(:,3)=y(:,3)./i0;
    
    isNegative=~isempty(strfind(magneticField{1},'-'));
    
    if isNegative==1
        magneticField{1}=strrep(magneticField{1},'-','N');
        maxR2ndHarmonic(end+1)=min(y(:,3));
    else
        maxR2ndHarmonic(end+1)=max(y(:,3));
    end
    
    eval(['b4t3_10nm_Py_' , thickness{1} , 'nm_' , magneticField{1} ,...
        'mT=y;'])
    
    
    
    % Plot
    figure(1)
    plot(y(:,1),y(:,2))
    hold on
    
    figure(2)
    plot(y(:,1),y(:,3))
    hold on
    
    % Fitting
    [rphe(end+1), off1(end+1), c1(end+1), kUnknow(end+1)] =  create1stFit(y(:,1),y(:,2),char(magneticForTag(i+1)));
    [rad(end+1), rfl(end+1), off(end+1), c(end+1)] = create2ndFit(y(:,1),y(:,3),char(magneticForTag(i+1)));
%     create2ndFit(y(:,1),y(:,3),char(magneticForTag(i+1)));
    
    clearvars y isNegative magneticField thickness
end

thicknessForTag(1)=[];
magneticForTag(1)=[];

% Show raw data plots
figure(1)
title("The 1st harmonic result");
legend(magneticForTag);
xlabel('Angle')
ylabel('Resistance (\Omega)')
grid on

figure(2)
title("The 2nd harmonic result");
legend(magneticForTag);
xlabel('Angle')
ylabel('Resistance (\Omega)')
grid on

figure
scatter(magneticForNum,maxR2ndHarmonic)
xlabel('Magnetic field (mT)')
ylabel('R^{2\omega}_{xy,max} (\Omega)')
box on
grid on

bInvK = 1./(magneticForNum./1000+bk);
bInv = 1./(magneticForNum./1000);
rflSub = rfl./rphe;

%% Fit: For B_AD
[xData, yData] = prepareCurveData( bInvK, rad );

% Set up fittype and options.
ft = fittype( 'poly1' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );

% Plot fit with data.
figure( 'Name','B_AD' );
plot(xData, yData ,'o');
hold on 
h = plot(fitresult);
set(h, 'LineStyle',':', 'LineWidth',2)
xlabel('1/(B_{ext}+B_k) (T^{-1})')
ylabel('R_{AD+\nablaT} (\Omega)')
grid on

bAD=fitresult.p1;

clearvars xData yData

%% Fit: For B_FL
[xData, yData] = prepareCurveData( bInv, rflSub );

% Set up fittype and options.
ft = fittype( 'poly1' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );

% Plot fit with data.
figure( 'Name','B_FL' );
plot(xData, yData ,'o');
hold on 
h = plot(fitresult);
set(h, 'LineStyle',':', 'LineWidth',2)
xlabel('1/B_{ext} (T^{-1})')
ylabel('R_{FL+Oe}/R_{PHE} (\Omega)')
grid on

bFLOE=fitresult.p1;

clearvars xData yData

%% Plot of 45deg
figure('Name','45deg')
plot(har45Deg(:,3),har45Deg(:,2))
xlabel 'H (Oe)'
ylabel 'R^{2\omega}_{xy} {\Omega}'
grid on

function y = importfile(filename, startRow, endRow)
%% Initialization
delimiter = '\t';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%% Format

formatSpec = '%f%f%f%[^\n\r]';

%% Open the data
fileID = fopen(filename,'r');

%% Read the data column.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the data.
fclose(fileID);

%% Name the data column.

y=[dataArray{1:end-1}];


end

function [rphe, off1, c1, kUnknow] = create1stFit(x, y, field)
%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( x, y );

% Set up fittype and options.
ft = fittype( 'rphe*sin(2*x+2*off1)+c1+kUnknow*cos(x+off1)', 'independent', 'x', 'dependent', 'y' );
% ft = fittype( 'rphe*sin(2*x+2*off1)+c1', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.DiffMinChange = 1e-10;
opts.Display = 'Off';
opts.MaxFunEvals = 600;
opts.MaxIter = 400;
opts.Robust = 'LAR';
opts.StartPoint = [0.1 0.278498218867048 0.2 0];
opts.TolFun = 1e-07;
opts.TolX = 1e-07;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name',field );
rphe = fitresult.rphe;
off1 = fitresult.off1;
c1 = fitresult.c1;
kUnknow = fitresult.kUnknow;
plot(xData, yData ,'o');
hold on 
h = plot(fitresult);
set(h, 'LineStyle',':', 'LineWidth',2)
xlabel('Angle')
ylabel('Resistance (\Omega)')
title(['Fitting result(1st) of ', field, 'mT case'])
grid on

end

function [rad, rfl, off, c] = create2ndFit(x, y, field)
[xData, yData] = prepareCurveData( x, y );

% Set up fittype and options.
ft = fittype( 'rad*cos(x+off)+rfl*(2*cos(x+off)^3-cos(x+off))+c', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.DiffMinChange = 1e-10;
opts.Display = 'Off';
opts.MaxFunEvals = 6000;
opts.MaxIter = 4000;
opts.Robust = 'Bisquare';
opts.StartPoint = [0 0 0 0];
opts.TolFun = 1e-08;
opts.TolX = 1e-08;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name',field );
rad = fitresult.rad;
rfl = fitresult.rfl;
off = fitresult.off;
c = fitresult.c;
radCurve = rad.*cos(xData+off)+c;
rflCurve = rfl.*(2.*(cos(xData+off).^3)-cos(xData+off))+c;
plot(xData, yData ,'o');
hold on 
h = plot(fitresult);
set(h, 'LineStyle',':', 'LineWidth',2)
hold on 
plot(  xData , radCurve , xData , rflCurve,'LineWidth',1.5)
legend('Data points', 'Fitting curve', 'AD contribution', 'FL contribution', 'Location', 'NorthEast' );
xlabel('Angle')
ylabel('Resistance (\Omega)')
title(['Fitting result of ', field, 'mT case'])
grid on
end

