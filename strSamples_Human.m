%%% 28 May 2014
%%% return the corresponding samples of a structure (defined by the acronym)
%%% the first ID in the output list is the id of the input structure
%%% NOTE:  THE FUNCTION RETURNS STRUCTURE ID (THE REFERENCE VOLUME INCLUDES THESE IDs)
%%%
%%% EXAMPLE: 
%%% inputStr = 'Br';
%%% ontologyFile = 'C:\Users\amahfouz\Documents\MATLAB\Data\ABA_Human_Data_Analysis\rawData_25Feb2014\Ontology.xlsx';
%%% strSamples_Human(inputStr, ontologyFile)

function subStructures = strSamples_Human(inputStr, ontologyFile, level)

% fID = fopen(ontologyFile);
% data = textscan(fID,'%d %q %q %d %c %d %q %q', ...
%     'delimiter', ',', 'HeaderLines', 1);
% fclose(fID);
% strIDs = data{1};
% strAcronyms = data{2};
% parentIDs = data{4};

% load the ontology data
[num txt] = xlsread(ontologyFile);
% extract structures acronyms
strAcronyms = txt(2:end,2);
strIDs = num(:,1);
parentIDs = num(:,4);
clear num; clear txt;
% find the index of the input structure
inputIndex = find(strcmpi(strAcronyms, inputStr) == 1, 1, 'first');
% check if the structure acronym is found
if isempty(inputIndex)
    disp([inputStr ' not found'])
    samplesInd = [];
    subStructures = [];
    return
else
    % first, add the main structure to the list
    samplesInd(1) = inputIndex;
    samplesInd = getChildren(samplesInd, inputIndex, strIDs, parentIDs, level);
    subStructures = strAcronyms(samplesInd);
end
end

% a function to recursively retrieve the children info
function sInd = getChildren(sInd, ind, strIDs, parentIDs, level)
    % get the current structure ID
    currStrID = strIDs(ind);
    % find children of the input structure
    chil = find(parentIDs == currStrID);
    if ~isempty(chil)
        sInd = [sInd; chil];
        if level == 0
            for i = 1 : length(chil)
                sInd = getChildren(sInd, chil(i), strIDs, parentIDs, level);
            end
        end
    end
end






