%%% 5 March 2015
%%% Read the gene expression data of DMD in LMD prenatal samples

%% select gene of interest
geneName = 'DMD';
probeNames = {'A\_23\_P321860','A\_32\_P199796','A\_24\_P342388','A\_23\_P113453','A\_24\_P185854','A\_24\_P34186'};
%%% Probe#5 has the highest connectivity
probe_id = 5;
donor_age = {'15 pcw', '16 pcw', '21 pcw', '21 pcw'};
donor_name = {'H376.IIIA.02', 'H376.IIIB.02', 'H376.IV.02', 'H376.IV.03'};
resDir = 'C:\Users\amahfouz\SURFdrive\Projects\DMD\Results\Prenatal_LMD_Microarray\';

%% get the expression of a gene from the AHBA
dataDir = 'C:\Users\amahfouz\SURFdrive\Projects\DMD\Data\Prenatal_LMD_Microarray_DMD\';
% read samples meta data
[num txt] = xlsread([dataDir 'Columns.xls']);
sample.donor_name = txt(2:end,2);
sample.donor_age = txt(2:end,3);
sample.structure_id = num(:,5);
sample.structure_name = txt(2:end,6);
sample.structure_acronym = txt(2:end,7);
sample.structure_color = txt(2:end,8);
sample.top_level_structure_name = txt(:,10);
sample.top_level_structure_acronym = txt(2:end,11);
sample.top_level_structure_color = txt(2:end,12);
clear num; clear txt;
% read probe meta data
[num txt] = xlsread([dataDir 'Rows.xls']);
probe.name = txt(2:end,2);
clear num; clear txt;
% read expression data
expData = csvread([dataDir 'Expression.csv']);
expData(:,1) = [];

%% visualize common regions across the 4 specimens 
% find unique regions per specimn and get the average expression
for i = 1 : 4
    specimen_samples_idx = find(ismember(sample.donor_name, donor_name{i}));
    uniqueSamples = unique(sample.structure_id(specimen_samples_idx));
    if i == 1
        allSamples = uniqueSamples;
    else
        allSamples = [allSamples; uniqueSamples];
    end
    for s = 1 : length(uniqueSamples)
        specimen_str_idx = find(sample.structure_id((specimen_samples_idx)) == uniqueSamples(s));
        specimen_samples{i}.sampleNum(s,1) = numel(specimen_str_idx);
        specimen_samples{i}.avgExp(s,1) = mean(expData(probe_id,specimen_samples_idx(specimen_str_idx)),2);
        specimen_samples{i}.stucture_id(s,1) = uniqueSamples(s);
    end
end
% find th set of all regions analyzed across the 4 donors
uniqueAllSamples = unique(allSamples);
% for each of the "all regions", get the expression from all donors 
for s = 1 : length(uniqueAllSamples)
    count = 0;
    for D = 1 : 4
        IX = find(specimen_samples{D}.stucture_id == uniqueAllSamples(s));
        if ~isempty(IX)
            count = count+1;
            if count == 1
%                 IX2 = find(donor_sample{D}.structure_id == uniqueAllSamples(s),1);
                Region.structure_name{s,1} = sample.structure_name{find(sample.structure_id==uniqueAllSamples(s),1)};
                Region.structure_acronym{s,1} = sample.structure_acronym{find(sample.structure_id==uniqueAllSamples(s),1)};
                Region.structure_color{s,1} = sample.structure_color{find(sample.structure_id==uniqueAllSamples(s),1)};
                Region.structure_id(s,1) = uniqueAllSamples(s);
                Region.top_level_structure_acronym{s,1} = sample.top_level_structure_acronym{find(sample.structure_id==uniqueAllSamples(s),1)};
                Region.top_level_structure_color{s,1} = sample.top_level_structure_color{find(sample.structure_id==uniqueAllSamples(s),1)};
            end
            Region.expression{s,1}(D) = specimen_samples{D}.avgExp(IX);
            Region.count{s,1} = count;
        else
            Region.expression{s,1}(D) = NaN;
        end
        clear IX;
    end
end
% find structures analyzed in all 4 specimens
IX_4 = find(cell2mat(Region.count) == 4);
% plot
expMat4 = cell2mat(Region.expression(IX_4));
barColors = Region.structure_color(IX_4);
for i = 1 : 4
    [~,sortingIX_4] = sort(expMat4(:,i), 'descend');
    % boxplot of the average expression of each region across the 6 donors
    f = figure('units','normalized','outerposition',[0 0 1 1]) 
    hold on
    for S = 1 : numel(IX_4)
        b = bar(S, expMat4(sortingIX_4(S),i));
        set(b,'EdgeColor','w','FaceColor',hex2rgb(Region.top_level_structure_color{IX_4(sortingIX_4(S))})/255);
        combinedLabel{S} = [Region.structure_acronym{IX_4(sortingIX_4(S))} ' ('  Region.top_level_structure_acronym{IX_4(sortingIX_4(S))} ')'];
    end
    hold off
    set(gca, 'ylim', [-3 3], 'xlim', [0 numel(IX_4)+1])%, 'XTick', 1:numel(IX_4))
    xticklabel_rotate(1:numel(IX_4),90,combinedLabel, 'FontSize', 8)
%     boxplot(expMat4(sortingIX_4,i)', Region.structure_acronym(IX_4(sortingIX_4)), ...
%         'colorgroup',Region.top_level_structure_acronym(IX_4(sortingIX_4)),...
%         'factorgap',5, 'factorseparator',1, 'labelorientation', 'inline', 'plotstyle', 'compact')
    grid on
    ylabel('Expression', 'FontWeight', 'bold', 'FontSize', 15)
    title([{geneName}; {['Donor ' donor_name{i} ' (' donor_age{i} ')']}], 'FontWeight', 'bold', 'FontSize', 15)
%     saveas(f, [resDir donor_name{i} '_' donor_age{i} '.png']);
%     saveas(f, [resDir donor_name{i} '_' donor_age{i} '.fig']);
end

%% boxplot of the expression
figure,
boxplot(expData(probe_id,:), ...
    {sample.donor_name,sample.top_level_structure_acronym}, ...
    'colorgroup',sample.top_level_structure_acronym, 'factorgap',5, ...
    'factorseparator',1, 'labelorientation', 'inline', 'plotstyle', 'compact')
grid on
ylabel('Expression (log_2)', 'FontWeight', 'bold', 'FontSize', 15)
line([])
title(['Probes: ' num2str(transcript_ranges(i,1)) ' - ' num2str(transcript_ranges(i,2))],...
    'FontWeight', 'bold', 'FontSize', 15);

