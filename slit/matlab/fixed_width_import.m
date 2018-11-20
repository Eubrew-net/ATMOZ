function [data] = fixed_width_import(filename,startline,number_of_lines,columns_width)
% FIXED_WIDTH_IMPORT - Imports data from fixed width text-files
%
% Reads data from a fixed width textfile (i.e. the numbers are arranged in
% columns that are a given number of characters wide). Non-numerical data 
% is converted in NaN. You can use this function for reading very large
% files in chunks, because you have to specify the line where to start 
% reading (first line = 1), and the number of lines you wish to read.
%
% Input:
%   - filename (string): filename to open (you have to include the extension)
%   - startline (int): line where to start reading
%   - number_of_lines (int): number of lines you want to read
%   - columns_width (vector): A vector containing for each column it's
%   width. (The width of a column is the number of characters it is wide.)
%
% Output: a matrix containing the read data.
%
% Example:
%   
% The file example.dat contains:
%--------------------------------------------------------------------------------- 
% This header will be changed into NaN!
% 1
% Here is some text
%    ROW   1   NODE    18     DEG. OF. FR. =  UX  
%     1 0.42781213E+01    2-0.31885678E+00    3 0.16586227E+00    4-0.15326625E+01
%     5-0.52177833E+00    6 0.10979865E+00    7 0.22749284E+00    8-0.51473827E+00
%     9-0.15022180E-01   10 0.81915685E+00   11 0.91007242E+00   12 0.78242910E+00
% The next line is a white line.
% 
%     1-0.31885678E+00    2 0.70635274E+01    3-0.19959239E-03    4 0.41480462E+00
%     5-0.26786897E-02    6-0.42934296E+01    7 0.84337345E+00    8-0.52414042E+00
%     9-0.13729736E-03   10 0.94989062E+01   11 0.33425710E+00   12 0.78303337E+00
% The next line is the last line.
%     1-0.31885678E+00    2 0.70635274E+01    3-0.19959239E-03    4 0.41480462E+00
%---------------------------------------------------------------------------------
%
% [data] = fixed_width_import('example.dat',2,14,[5 15 5 15 5 15 5 15])
%
%   will then return one warning for reading too many lines, and
%
% data =
% 
%     1.0000       NaN       NaN       NaN       NaN       NaN       NaN       NaN
%        NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN
%        NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN
%        NaN       NaN   18.0000       NaN       NaN       NaN       NaN       NaN
%     1.0000    4.2781    2.0000   -0.3189    3.0000    0.1659    4.0000   -1.5327
%     5.0000   -0.5218    6.0000    0.1098    7.0000    0.2275    8.0000   -0.5147
%     9.0000   -0.0150   10.0000    0.8192   11.0000    0.9101   12.0000    0.7824
%        NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN
%        NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN
%     1.0000   -0.3189    2.0000    7.0635    3.0000   -0.0002    4.0000    0.4148
%     5.0000   -0.0027    6.0000   -4.2934    7.0000    0.8434    8.0000   -0.5241
%     9.0000   -0.0001   10.0000    9.4989   11.0000    0.3343   12.0000    0.7830
%        NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN
%     1.0000   -0.3189    2.0000    7.0635    3.0000   -0.0002    4.0000    0.4148
%        NaN       NaN       NaN       NaN       NaN       NaN       NaN       NaN
%
% 01/05/2006, Adriaan Van Nuffel

% Open file for reading only, in text mode.
fid = fopen(filename,'rt');

if fid == -1
    % If file doesn't exist: exit from the program with an error message.
    error(['File ' filename ' not found. Check if the file exists in the current work directory (' cd ').'])
end

nr_cols = length(columns_width);
% Number of columns to be read
total_width = sum(columns_width);
% Total width of columns to be read. If the file has lines that aren't as
% wide as total_width, the rest of the data will be read as white space (NaN)

if startline > 1
    % Go to the wanted line (lines are read but nothing is done)
    for i = 1:(startline - 1)
        fgetl(fid); % Read (skip) a line
    end
end

% Now that the starting line has been reached: read the wanted data
data = zeros(number_of_lines,nr_cols); % Preallocating memory to speed things up.
for i = 1:number_of_lines

    % Read a line
    textline = fgetl(fid);
    
    % Check if the line is long enough, if not, add whitespace.
    if length(textline) < total_width
        textline(end+1:end+total_width)=' ';
    end
    
    % Converting the read line in numerical data
    startcol = 1;              % Start reading from character nº startcol
    endcol = columns_width(1); % ... up to character nº endcol
    for j = 1:nr_cols
        endcol = startcol + columns_width(j) - 1;
        data_txt = textline(startcol:endcol);
        startcol = endcol + 1;
        % Convert the text data in numbers or NaN
        data(i,j) = str2double(data_txt);
    end
    
end

% Close the file.
fclose(fid);