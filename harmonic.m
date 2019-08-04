clear;close all;clc;

b=importfile("B4T3_10nm_Py_5nm_harmonic900(mT)_300(K)_7(Vrms)_.txt",5,900);



function y = importfile(filename, thickness, magneticField, startRow, endRow)
%% Initialization
delimiter = '\t';
if nargin<=4
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
