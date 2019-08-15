clear;close all;clc;

fileName=dir(fullfile('*.txt'));
%fileName= dir(fullfile('d:/datafile','*.txt'));             
%The path should be changed to what you want.

number=length(fileName);
thicknessForTag=string;
magneticForTag=string;

magneticForNum=[];
maxR2ndHarmonic=[];
rad = [];
rfl = [];
off = [];
c = [];
rphe = [];
off1 = [];
c1 = [];

for i=1:number
    % Here we can add an if to separate 2 types of measurements.
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
    figure
    plot(y(:,1),y(:,2))
    title(['The 1st harmonic result of ',char(magneticForTag(i+1)),'mT'])
    xlabel('Angle')
    ylabel('Votage (V)')
    hold on
    grid on
    
    figure
    plot(y(:,1),y(:,3))
    title(['The 2nd harmonic result of ',char(magneticForTag(i+1)),'mT'])
    xlabel('Angle')
    ylabel('Votage (V)')
    grid on
    
    % Fitting
    [rphe(end+1), off1(end+1), c1(end+1)] =  create1stFit(y(:,1),y(:,2),char(magneticForTag(i+1)));
    [rad(end+1), rfl(end+1), off(end+1), c(end+1)] = create2ndFit(y(:,1),y(:,3),char(magneticForTag(i+1)));
%     create2ndFit(y(:,1),y(:,3),char(magneticForTag(i+1)));
    
    clearvars y isNegative magneticField thickness
end

thicknessForTag(1)=[];
magneticForTag(1)=[];

% figure(1)
% title("2nd harmonic");
% legend(magneticForTag);
% xlabel('Degree')
% ylabel('Votage (V)')
% 
% figure(2)
% title("1st harmonic");
% legend(magneticForTag);
% xlabel('Degree')
% ylabel('Votage (V)')

figure
scatter(magneticForNum,maxR2ndHarmonic)
xlabel('Magnetic field (mT)')
ylabel('Votage (V)')
box on
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

function [rphe, off1, c1] = create1stFit(x, y, field)
%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( x, y );

% Set up fittype and options.
ft = fittype( 'rphe*sin(2*x+2*off1)+c1', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.DiffMinChange = 1e-10;
opts.Display = 'Off';
opts.MaxFunEvals = 600;
opts.MaxIter = 400;
opts.Robust = 'LAR';
opts.StartPoint = [0.0975404049994095 0.278498218867048 0.546881519204984];
opts.TolFun = 1e-06;
opts.TolX = 1e-06;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name',field );
rphe = fitresult.rphe
off1 = fitresult.off1
c1 = fitresult.c1
plot(xData, yData ,'o');
hold on 
h = plot(fitresult);
set(h, 'LineStyle',':', 'LineWidth',2)
xlabel x
ylabel y
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
xlabel x
ylabel y
title(['Fitting result of ', field, 'mT case'])
grid on
end

