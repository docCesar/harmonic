clear;close all;clc;

fileName=dir(fullfile('*.txt'));
%fileName= dir(fullfile('d:/datafile','*.txt'));             
%The path should be changed to what you want.

number=length(fileName);
thicknessForTag=string;
magneticForTag=string;

magneticForNum=[];
maxR2ndHarmonic=[];

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
    plot(y(:,1),y(:,3))
    hold on
    grid on
    
    figure(2)
    plot(y(:,1),y(:,2))
    hold on
    grid on
    
    figure
    plot(y(:,1),y(:,3))
    title(['The 2nd harmonic result of ',char(magneticForTag(i+1)),'mT'])
    xlabel('H (mT)')
    ylabel('Votage (V)')
    grid on
    
    clearvars y isNegative magneticField thickness
end

thicknessForTag(1)=[];
magneticForTag(1)=[];

figure(1)
title("2nd harmonic");
legend(magneticForTag);
xlabel('H (mT)')
ylabel('Votage (V)')

figure(2)
title("1st harmonic");
legend(magneticForTag);
xlabel('H (mT)')
ylabel('Votage (V)')

figure
scatter(magneticForNum,maxR2ndHarmonic)
xlabel('H (mT)')
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
