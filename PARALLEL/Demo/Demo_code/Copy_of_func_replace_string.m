
    
function f = Copy_of_func_replace_string(InputFile, OutputFile, SearchString, ReplaceString)
    %% change data [e.g. initial conditions] in model file 
% InputFile
fid = fopen(InputFile); data = textscan(fid, '%s', 'Delimiter', '\n', 'CollectOutput', true); fclose(fid);
    
    
for I = 1:length(data{1}) tf = strcmp(data{1}{I}, SearchString); % search for this string in the array 
    if tf == 1 data{1}{I} = ReplaceString; % replace with this string 
    end
end
    
    
% read whole model file data into cell array 

% modify the cell array 
% find the position where changes need to be applied and insert new data 
    

% write the modified cell array into the text file 
    

fid = fopen(OutputFile, 'w'); 
for I = 1:length(data{1}) 
    fprintf(fid, '%s\n', char(data{1}{I})); 
end
fclose(fid);

    
