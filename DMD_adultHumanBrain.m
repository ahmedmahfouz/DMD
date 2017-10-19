%%% 26 Jan 2015
%%% Read the gene expression data of DMD in the adult human brain

%% select gene of interest
geneName = 'DMD';
donors = {'10021','12876','14380','15496','15697','9861'};
probeNames = {'A\_23\_P321860','A\_32\_P199796','A\_24\_P342388','A\_23\_P113453','A\_24\_P185854','A\_24\_P34186'};
% mainStructures = {'N/A', 'FL', 'PL', 'TL', 'OL', 'HiF', 'Str', 'GP', 'Amg', 'TH', ...
%     'Hy', 'MES', 'Pons', 'MY', 'Cb', 'WM'};
mainStructures = {'FL','OL','PL','TL','Ins','CgG','HiF','PHG','Amg',...
    'BF','GP','Str','Cl','ET','Hy','SbT','TH','MES','Cb','Pons',...
    'MY','WM','SS'};
%%% Probe#5 has the highest connectivity
probe_id = 5;

%% get the expression of a gene from the AHBA
dataDir = '/tudelft.net/staff-bulk/ewi/insy/VisionLab/mvandegiessen/tSNE_ABA/rawData_25Feb2014/';
load([dataDir 'normalized_microarray_donor' donors{1} '/probe.mat']);
gene_idx = find(strcmpi(probe.gene_symbol,geneName));
for i = 1 : length(donors)
    load([dataDir 'normalized_microarray_donor' donors{i} '/sample.mat']);
    load([dataDir 'normalized_microarray_donor' donors{i} '/sampleLables.mat']);
    donor_sample{i} = sample;
    donor_sampleLabels{i} = sampleLables;
    load([dataDir 'normalized_microarray_donor' donors{i} '/MicroarrayExpression.mat']);
    gene_expression_donors{i} = MicroarrayExpression(gene_idx,:);
end
save('donor_sample.mat','donor_sample');
save('donor_sampleLabels.mat','donor_sampleLabels');
save('gene_expression_donors.mat','gene_expression_donors');

%% get the expression of a gene from the AHBA (Sjoerd's data version)
% dataDir = '/tudelft.net/staff-bulk/ewi/insy/DBL/sjoerdhuisman/ABA_human_brain_probegene/';
% T = readtable([dataDir 'probe_info_2014-11-11.csv']);
% probe.gene_symbol = T.gene_symbol;
% probe.entrez_id = T.entrez_id;
% gene_idx = find(strcmpi(probe.gene_symbol,geneName));
% clear T;
% for i = 1 : length(donors)
%     T = readtable([dataDir 'sample_info_normalized_microarray_donor' donors{i} '_2014-11-11.txt'],'Delimiter','\t');
%     sample{i}.structure_id = T.structure_id;
%     sample{i}.structure_acronym = T.structure_acronym;
%     sample{i}.structure_name = T.structure_name;
%     clear T;
%     MicroarrayExpression = csvread([dataDir 'gene_expr_normalized_microarray_donor' donors{i} '_2014-11-11.csv']);
%     gene_expression_donors{i} = MicroarrayExpression(gene_idx,:);
% end
% save('sample.mat','sample');
% save('probe.mat','probe');
% save('DMD_expression_donors.mat','gene_expression_donors');

%% Define directories and load data
if ispc
    dataDir = 'C:/Users/amahfouz/SURFdrive/Projects/DMD/Data/';
end
load([dataDir 'donor_sample.mat']);
load([dataDir 'donor_sampleLabels.mat']);
load([dataDir 'gene_expression_donors.mat']);

%% Define directories and load data (Sjoerd's data version)
if ispc
    dataDir = 'C:/Users/amahfouz/SURFdrive/Projects/DMD/Data/';
end
load([dataDir 'sample.mat']);
load([dataDir 'probe.mat']);
load([dataDir 'DMD_expression_donors.mat']);

%% Plot the expression of the gene
% concatinate expression data & generate labels
for i = 1 : length(gene_expression_donors)
    if i == 1
        expData = gene_expression_donors{i};
        donorLabel = ones(1,length(gene_expression_donors{i}))*i;
        strLabel = donor_sampleLabels{i};
    else
        expData = [expData, gene_expression_donors{i}];
        donorLabel = [donorLabel, ones(1,length(gene_expression_donors{i}))*i];
        strLabel = [strLabel,donor_sampleLabels{i}];
    end
end
% generate string labels
strLabel_name = cell(1,numel(strLabel));
donorLabel_name = cell(1,numel(strLabel));
for i = 1 : numel(unique(strLabel))
    strLabel_name(strLabel==i-1) = mainStructures(i);
end
for i = 1 : numel(unique(donorLabel))
    donorLabel_name(donorLabel==i) = donors(i);
end
% select probe to plot
probe = 5;
figure,
boxplot(expData(probe,:), {strLabel_name,donorLabel_name}, ...
    'colorgroup',strLabel_name, 'factorgap',5, ...
    'factorseparator',1, 'labelorientation', 'inline')
grid on
ylabel('Expression', 'FontWeight', 'bold', 'FontSize', 15)
title([geneName ' - Probe#' num2str(probe)], 'FontWeight', 'bold', 'FontSize', 15)
figure,
boxplot(expData(probe,:), {donorLabel_name,strLabel_name}, ...
    'colorgroup',strLabel_name, 'factorgap',5, ...
    'factorseparator',1, 'labelorientation', 'inline')
grid on
ylabel('Expression', 'FontWeight', 'bold', 'FontSize', 15)
title([geneName ' - Probe#' num2str(probe)], 'FontWeight', 'bold', 'FontSize', 15)

%% Plot the expression of the gene averaged per structre
% calculate the average per structure for each probe and donor
ontologyFile = 'C:\Ahmed\Work\Data\ABA_Human_Data_Analysis\Ontology.xlsx';
for i = 1 : length(gene_expression_donors)
    tempExp = zscore(gene_expression_donors{i},[],2);
    for j = 1 : length(mainStructures)
        str_children = strSamples_Human(mainStructures{j}, ontologyFile);
        region_idx = ismember(sample{i}.structure_acronym, str_children);
%         avgExpMat(:,i,j) = mean(gene_expression_donors{i}(:,region_idx),2);
        avgExpMat(:,i,j) = mean(tempExp(:,region_idx),2);
        regionSize(i,j) = sum(region_idx);
        clear region_idx;
%         avgExpMat(:,i,j) = mean(gene_expression_donors{i}(:,donor_sampleLabels{i} == j-1),2);
    end
end
% probe = 5;
figure,
boxplot(squeeze(avgExpMat(:,:,:)), mainStructures, ...
    'colorgroup',mainStructures, 'factorgap',5, ...
    'factorseparator',1, 'labelorientation', 'inline', 'plotstyle', 'compact')
grid on
ylabel('Expression', 'FontWeight', 'bold', 'FontSize', 15)
title([{geneName}; {['Average Expression - Probe: ' probeNames{probe}]}], 'FontWeight', 'bold', 'FontSize', 15)

%% Calculate the correlation of all genes with DMD for each donor separately
strOfInterest = {'Br'};
% get the expression of a gene from the AHBA (probes mapped to genes)
dataDir = '/tudelft.net/staff-bulk/ewi/insy/DBL/sjoerdhuisman/ABA_human_brain_probegene/';
resDir = '/tudelft.net/staff-bulk/ewi/insy/DBL/amahfouz/MATLAB/Results/DMD/';
% read probe info
T = readtable([dataDir 'probe_info_2014-11-11.csv']);
probe.gene_symbol = T.gene_symbol;
probe.entrez_id = T.entrez_id;
gene_idx = find(strcmpi(probe.gene_symbol,geneName));
clear T;
% calculate correlations wih the 5th prob of DMD in each donor separately
for i = 1 : length(donors)
    % gene expression data
    MicroarrayExpression = csvread([dataDir 'gene_expr_normalized_microarray_donor' donors{i} '_2014-11-11.csv']);    
%     MicroarrayExpression = zscore(MicroarrayExpression,[],2);
    % select samples
    for S = 1 : length(strOfInterest)
        sample_children = strSamples_Human(strOfInterest{S}, [dataDir 'Ontology.csv']);
        % read sample info
        dataDir2 = '/tudelft.net/staff-bulk/ewi/insy/VisionLab/mvandegiessen/tSNE_ABA/rawData_25Feb2014/';
        load([dataDir2 'normalized_microarray_donor' donors{i} '/sample.mat']);
        if S == 1
            sample_idx = find(ismember(sample.structure_acronym, sample_children));
        else
            sample_idx = [sample_idx; find(ismember(sample.structure_acronym, sample_children))];
        end
    end
    % calculate the correlation
    corrMat(i,:) = corr(MicroarrayExpression(gene_idx,sample_idx)', MicroarrayExpression(:,sample_idx)');
    % keep track of structure size
    strSize = numel(sample_idx);
    clear MicroarrayExpression; clear sample_idx;
end
% create file name
fileNameExt = strOfInterest{1};
for i = 2 : length(strOfInterest)
    fileNameExt = [fileNameExt '_' strOfInterest{i}];
end
save([resDir fileNameExt '_corrMat.mat'],'corrMat','probe','strSize');

%% Analyze the correlation lists
dataDir = 'C:/Ahmed/Work/Data/DMD/';
resultsDir = 'C:/Ahmed/Work/Results/DMD/';
strOfInterest = {'Cb'};
fileNameExt = strOfInterest{1};
for i = 2 : length(strOfInterest)
    fileNameExt = [fileNameExt '_' strOfInterest{i}];
end
load([dataDir fileNameExt '_corrMat.mat']);
for i = 1 : length(donors)
    [sortedCorrMat(i,:) IX(i,:)] = sort(corrMat(i,:),2,'descend');
    xlswrite([resultsDir geneName '_' fileNameExt '_possitively_coexpressed_genes.xlsx'], probe.gene_symbol(IX(i,2:end)), i, 'A1');
    [sortedCorrMat(i,:) IX(i,:)] = sort(corrMat(i,:),2);
    xlswrite([resultsDir geneName '_' fileNameExt '_negatively_coexpressed_genes.xlsx'], probe.gene_symbol(IX(i,1:end-1)), i, 'A1');
end
xlsheets(donors, [resultsDir geneName '_' fileNameExt '_possitively_coexpressed_genes.xlsx']);
xlsheets(donors, [resultsDir geneName '_' fileNameExt '_negatively_coexpressed_genes.xlsx']);

%% Calculate the correlation of all genes with DMD across ALL donor 
% get the expression of a gene from the AHBA (probes mapped to genes)
dataDir = '/tudelft.net/staff-bulk/ewi/insy/DBL/sjoerdhuisman/ABA_human_brain_probegene/';
resDir = '/tudelft.net/staff-bulk/ewi/insy/DBL/amahfouz/MATLAB/Results/DMD/';
% read probe info
T = readtable([dataDir 'probe_info_2014-11-11.csv']);
probe.gene_symbol = T.gene_symbol;
probe.entrez_id = T.entrez_id;
gene_idx = find(strcmpi(probe.gene_symbol,geneName));
clear T;
for i = 1 : length(donors)
    % gene expression data
    if i == 1
        MicroarrayExpression = csvread([dataDir 'gene_expr_normalized_microarray_donor' donors{i} '_2014-11-11.csv']);
        MicroarrayExpression = zscore(MicroarrayExpression,[],2);
    else
        tempExpMat = csvread([dataDir 'gene_expr_normalized_microarray_donor' donors{i} '_2014-11-11.csv']);
        tempExpMat = zscore(tempExpMat,[],2);
        MicroarrayExpression = [MicroarrayExpression, tempExpMat];
        clear tempExpMat;
    end
end
% calculate the correlation
[corrMat_combined pVal] = corr(MicroarrayExpression(gene_idx,:)', MicroarrayExpression');
save([resDir 'corrMat_combined.mat'],'corrMat_combined','pVal','probe');

%% Analyze the combined donor correlation list
dataDir = 'C:/Users/amahfouz/SURFdrive/Projects/DMD/Data/';
resultsDir = 'C:/Users/amahfouz/SURFdrive/Projects/DMD/Results/';
load([dataDir 'corrMat_combined.mat']);
% sort corelations
[sortedCorrMat IX] = sort(corrMat_combined,2,'descend');
% save data to excel
xlswrite([resultsDir geneName '_combined_possitively_coexpressed_genes.xlsx'], probe.gene_symbol(IX(2:end)), 1, 'A1');
xlswrite([resultsDir geneName '_combined_possitively_coexpressed_genes.xlsx'], probe.entrez_id(IX(2:end)), 1, 'B1');
xlswrite([resultsDir geneName '_combined_possitively_coexpressed_genes.xlsx'], sortedCorrMat(2:end)', 1, 'C1');
% [sortedCorrMat IX] = sort(corrMat_combined,2);
xlswrite([resultsDir geneName '_combined_negatively_coexpressed_genes.xlsx'], probe.gene_symbol(IX(1:end-1)), 1, 'A1');
xlsheets(donors, [resultsDir geneName '_combined_possitively_coexpressed_genes.xlsx']);
xlsheets(donors, [resultsDir geneName '_combined_negatively_coexpressed_genes.xlsx']);
% plot correlation and p-value
figure, hold on
subplot(2,1,1), hold on
line([0 length(sortedCorrMat)-1], [0 0],'LineStyle','--','Color',[0.5,0.5,0.5],'LineWidth',2)
plot(2:201,sortedCorrMat(2:201),'LineWidth',3,'Color','r'),
plot(202:numel(sortedCorrMat)-200,sortedCorrMat(202:end-200),'LineWidth',3),
plot(numel(sortedCorrMat)-199:numel(sortedCorrMat),sortedCorrMat(end-199:end),'LineWidth',3,'Color','r'),
grid on, hold off
ylabel('Correlation', 'FontWeight', 'bold', 'FontSize', 15)
xlabel('Genes sorted on corrleation')
set(gca,'XTick',[],'XTickLabel',[])
title([{['Correlation to ' geneName]}; {['Probe: ' probeNames{probe_id}]}], 'FontWeight', 'bold', 'FontSize', 15)
subplot(2,1,2), bar(-log10(pVal(IX(2:end))),'b','EdgeColor','w'), grid on
ylabel('-log_1_0 (p-value)', 'FontWeight', 'bold', 'FontSize', 15)
xlabel('Genes sorted on corrleation')
set(gca,'XTick',[],'XTickLabel',[])
hold off

figure, hold on
line([0 length(sortedCorrMat)-1], [0 0],'LineStyle','--','Color',[0.5,0.5,0.5],'LineWidth',2)
plot(2:201,sortedCorrMat(2:201),'LineWidth',3,'Color',[0.8 0.2 0]),
plot(202:numel(sortedCorrMat)-200,sortedCorrMat(202:end-200),'LineWidth',3,'Color',[0.2 0.2 0.2]),
plot(numel(sortedCorrMat)-199:numel(sortedCorrMat),sortedCorrMat(end-199:end),'LineWidth',3,'Color',[0.8 0.2 0]),
grid on
grid minor
ylabel('Correlation', 'FontWeight', 'bold', 'FontSize', 15)
xlabel('All 19,992 genes sorted on corrleation', 'FontWeight', 'bold', 'FontSize', 15)
% xlabel('Genes sorted on corrleation')
set(gca,'XTick',[],'XTickLabel',[])
hold off

%% plot a heatmap of the most/least correlated genes
N = 25;
donorNo = 1;
resultsDir = 'C:/Users/amahfouz/SURFdrive/Projects/DMD/Results/';
load([dataDir 'corrMat_combined.mat']);
load([dataDir 'corrGeneExp.mat']);
[XXX,ix_probes] = sort(corrMat_combined,'descend');
% % plot the correlation of donor#1 ()
% figure, imagesc([topCorr(:,1:893);bottomCorr(:,1:893)]), colormap('redbluecmap')
% caxis([-3 3])
% set(gca, 'XTick', 1:893, 'XTickLabel', donor_sample{1,1}.structure_acronym)

% find higher level annotation of a set of samples
ontologyFile = 'C:\Users\amahfouz\SURFdrive\Data\ABA_adult_human_brain\Ontology.xlsx';
[num txt] = xlsread(ontologyFile);
structure.id = num(:,1);
structure.acronym = txt(2:end,2);
structure.parent_structue = num(:,4);
structure.hemisphere = txt(2:end,5);
structure.structure_id_path = txt(2:end,7);
structure.hexCOLOR = txt(2:end,8);
xx = find(ismember(structure.hexCOLOR,'')==1);
yy = num(~isnan(num(:,8)),8);
for i = 1 : numel(xx)
    if yy(i) > 99999
        structure.hexCOLOR(xx(i)) = str2cell(num2str(yy(i)));
    else
        structure.hexCOLOR(xx(i)) = str2cell(['0' num2str(yy(i))]);
    end
end
structure.Order = num(:,6);
mainStructures = {'FL','OL','PL','TL','Ins','CgG','HiF','PHG','Pir','Amg',...
    'AO','BF','GP','Str','Cl','ET','Hy','SbT','TH','MES','Cb','Pons',...
    'MY','WM','SS'};
mainStrID = structure.id(ismember(structure.acronym, mainStructures));
% find the parent of all samples of selected donor
for i = 1 : length(donor_sample{donorNo}.structure_id)
    sampleOrder(i) = structure.Order(structure.id==donor_sample{donorNo}.structure_id(i));
    sampleColor(i) = structure.hexCOLOR(structure.id==donor_sample{donorNo}.structure_id(i));
    sampleColor_RGB(i,:) =  hex2rgb(sampleColor{i});
    for j = 1 : length(mainStrID)
        if findstr(structure.structure_id_path{structure.id==donor_sample{donorNo}.structure_id(i)}, num2str(mainStrID(j)));
            leafStr_parent(i,1) = mainStrID(j);
            leafStr_parent_acronym{i,1} = structure.acronym{structure.id==mainStrID(j)};
            leafStr_parent_hexCOLOR{i,1} = structure.hexCOLOR{structure.id==mainStrID(j)};
        end
    end    
end
% sort the samples based on parent structures
rowNames = [probe.gene_symbol(ix_probes(1:N+1)); probe.gene_symbol(ix_probes(end-N+1:end))];
% [colNames IX] = sort(leafStr_parent_acronym);
% [colLabels.Labels, m] = unique(leafStr_parent_acronym);
% % C = jet(length(unique(colNames)));
% C = hex2rgb(leafStr_parent_hexCOLOR(m))/255;
% % colLabels.Colors = C;
% for i = 1 : length(unique(colLabels.Labels))
%     colLabels.Colors(i) = mat2cell(C(i,:));
% end
% H = HeatMap([topCorr(:,IX);bottomCorr(:,IX)], ...
%     'RowLabels', rowNames, 'ColumnLabels', colNames,...
%     'Colormap', 'redbluecmap', 'DisplayRange', 3,...
%     'LabelsWithMarkers', true, 'ColumnLabelsColor', colLabels);

[~,IX2] = sort(sampleOrder);
colNames = donor_sample{donorNo}.structure_acronym(IX2);
colNames2 = leafStr_parent_acronym(IX2);

[colLabels.Labels, m] = unique(donor_sample{donorNo}.structure_acronym);
C = hex2rgb(sampleColor(m))/255;
for i = 1 : length(colLabels.Labels)
%     colLabels.Colors(i) = mat2cell(C(i,:));
    temp = mat2cell(C(i,:),1,3);
    colLabels.Colors{i,1} = temp{:};
end
H = HeatMap([topCorr(1:N+1,IX2);bottomCorr(end-N+1:end,IX2)]', ...
    'RowLabels', colNames, 'ColumnLabels', rowNames,...
    'Colormap', 'redbluecmap', 'DisplayRange', 3,...
    'LabelsWithMarkers', true, 'RowLabelsColor', colLabels);

hF = plot(H);
cF = get(0,'CurrentFigure');
set(cF,'CurrentAxes',findobj(hF,'Type','Axes'));
X = flipud([topCorr(1:N+1,IX2);bottomCorr(end-N+1:end,IX2)]');
hold on
for i = 1:size(X,2)
    plot([i-.5,i-.5],[.5,size(X,1)+.5],'k-');
end
hold off

% figure
% imagesc(X,'CDataMapping','scaled'), colormap('redbluecmap')
% set(gca, 'CLim', [-3 3])
% set(gca,'YTick',[])
% set(gca,'XTick',1:length(rowNames),'XTickLabel',rowNames)
% rotateXLabels(gca,90)

%% rank regions based on their DMD expression and plot the average expression
% select Probe#5 (A_24_P185854) has the highest connectivity
probe_id = 1; % probe =1 if using Sjoerd's data
% for each donor, calculate the average expression per unique samples
for D = 1 : length(donors)
    MicroarrayExpression = zscore(gene_expression_donors{D}(probe_id,:),[],2);
    uniqueSamples = unique(donor_sample{D}.structure_id);
    if D == 1
        allSamples = uniqueSamples;
    else
        allSamples = [allSamples; uniqueSamples];
    end
    for s = 1 : length(uniqueSamples)
        donorSamples{D}.sampleNum(s,1) = numel(find(donor_sample{D}.structure_id == uniqueSamples(s)));
        donorSamples{D}.avgExp(s,1) = mean(MicroarrayExpression(1,donor_sample{D}.structure_id == uniqueSamples(s)));
        donorSamples{D}.sample_id(s,1) = uniqueSamples(s);
    end
end
% list all samples analyzed
resDir = 'C:/Users/amahfouz/SURFdrive/Projects/DMD/Results/AdultHumanBrain/';
uniqueAllSamples = unique(allSamples);
% save([dataDir 'uniqueAllSamples.mat'], 'uniqueAllSamples');
for s = 1 : length(uniqueAllSamples)
    count = 0;
    for D = 1 : length(donors)
        IX = find(donorSamples{D}.sample_id == uniqueAllSamples(s));
        if ~isempty(IX)
            count = count+1;
            if count == 1
                IX2 = find(donor_sample{D}.structure_id == uniqueAllSamples(s),1);
                Region.name{s,1} = donor_sample{D}.structure_name{IX2};
                Region.structure_acronym{s,1} = donor_sample{D}.structure_acronym{IX2};
                Region.structure_id{s,1} = donor_sample{D}.structure_id(IX2);
            end
            Region.donorExp{s,1}(D) = donorSamples{D}.avgExp(IX);
            Region.sampelCountPeDonor{s,1}(D) = donorSamples{D}.sampleNum(IX);
            Region.count{s,1} = count;
        else
            Region.donorExp{s,1}(D) = NaN;
            Region.sampelCountPeDonor{s,1}(D) = NaN;
        end
        clear IX;
    end
    Region.avgExp(s,1) = nanmean(Region.donorExp{s,1});
    Region.totalNoSamples(s,1) = nansum(Region.sampelCountPeDonor{s,1});
end
% save([dataDir 'Region.mat'], 'Region');
% load([dataDir 'Region.mat']);
% add higher level annotation
ontologyFile = 'C:/Users/amahfouz/SURFdrive/Data/ABA_adult_human_brain/Ontology.xlsx';
[num txt] = xlsread(ontologyFile);
structure.id = num(:,1);
structure.acronym = txt(2:end,2);
structure.parent_structue = num(:,4);
structure.hemisphere = txt(2:end,5);
structure.structure_id_path = txt(2:end,7);
mainStructures = {'FL','OL','PL','TL','Ins','CgG','HiF','PHG','Pir','Amg',...
    'AO','BF','GP','Str','Cl','ET','Hy','SbT','TH','MES','Cb','Pons',...
    'MY','WM','SS'};
mainStrID = structure.id(ismember(structure.acronym, mainStructures));
for i = 1 : length(uniqueAllSamples)
    for j = 1 : length(mainStrID)
        if findstr(structure.structure_id_path{structure.id==uniqueAllSamples(i)}, num2str(mainStrID(j)));
            leafStr_parent(i,1) = mainStrID(j);
            leafStr_parent_acronym{i,1} = structure.acronym{structure.id==mainStrID(j)};
        end
    end    
end
% sort the regions
[~,sortingIX] = sort(Region.avgExp, 'descend');
xlswrite([resDir 'DMD_AdultHumanBrain.xlsx'], [{'structure_name'},...
    {'structure_acronym'},{'higher_order_annotation'},{'bains_used'},{'mean_DMD'},...
    strcat(donors,'_mean_DMD'),strcat(donors,'_numberOfSamples'),{'total_number_of_samples'},...
    {'structure_id'}], 1, 'A1');
xlswrite([resDir 'DMD_AdultHumanBrain.xlsx'], [Region.name(sortingIX), ...
    Region.structure_acronym(sortingIX)], 1, 'A2');
xlswrite([resDir 'DMD_AdultHumanBrain.xlsx'], leafStr_parent_acronym(sortingIX), 1, 'C2');
xlswrite([resDir 'DMD_AdultHumanBrain.xlsx'], Region.count(sortingIX), 1, 'D2');
xlswrite([resDir 'DMD_AdultHumanBrain.xlsx'], Region.avgExp(sortingIX,:), 1, 'E2');
xlswrite([resDir 'DMD_AdultHumanBrain.xlsx'], cell2mat(Region.donorExp(sortingIX)), 1, 'F2');
xlswrite([resDir 'DMD_AdultHumanBrain.xlsx'], cell2mat(Region.sampelCountPeDonor(sortingIX)), 1, 'L2');
xlswrite([resDir 'DMD_AdultHumanBrain.xlsx'], Region.totalNoSamples(sortingIX), 1, 'R2');
xlswrite([resDir 'DMD_AdultHumanBrain.xlsx'], Region.structure_id(sortingIX), 1, 'S2');
% select and sort regions sampled in all 6 brains
IX_6 = find(cell2mat(Region.count) == 6);
[~,sortingIX_6] = sort(Region.avgExp(IX_6), 'descend');
% create labels
for i = 1 : length(sortingIX_6)
    combinedLabel{i} = [Region.structure_acronym{IX_6(sortingIX_6(i))} ' ('  leafStr_parent_acronym{IX_6(sortingIX_6(i))} ')'];
end
% boxplot of the average expression of each region across the 6 donors
figure,
boxplot(cell2mat(Region.donorExp(IX_6(sortingIX_6),:))', combinedLabel, ...
    'colorgroup',leafStr_parent_acronym(IX_6(sortingIX_6)),...
    'factorgap',5, 'factorseparator',1, 'labelorientation', 'inline', 'plotstyle', 'compact')
grid on
ylabel('Expression (z-score)', 'FontWeight', 'bold', 'FontSize', 15)
title([{geneName}; {['Average Expression - Probe: ' probeNames{probe_id}]}], 'FontWeight', 'bold', 'FontSize', 15)
% plot the number of samples per reion
figure,
bar(Region.totalNoSamples(IX_6(sortingIX_6)),'b','EdgeColor','w')
xlim([0 numel(IX_6)+1])
ylabel('Number of samples', 'FontWeight', 'bold', 'FontSize', 15)
title('Number of samples per region in all 6 donors', 'FontWeight', 'bold', 'FontSize', 15)
grid on
xticklabel_rotate([1:numel(IX_6)],90,combinedLabel) 

% analyze significance of brain regions
sortedRegions = leafStr_parent_acronym(sortingIX);
uniqueRegions = unique(leafStr_parent_acronym);
for i = 1 : length(uniqueRegions)
    ranks = find(ismember(sortedRegions,uniqueRegions{i})==1);
    allRanks = setdiff(1:length(sortedRegions),ranks);
    [p(i,1),~,stats] = ranksum(ranks, allRanks);
    rSum(i,1) = stats.ranksum;
    nSamples(i,1) = numel(ranks);
end
p_corrected = multtest(p, 'method', 'holm');
xlswrite([resDir 'regionStats.xlsx'], uniqueRegions, 1, 'A2');
xlswrite([resDir 'regionStats.xlsx'], p_corrected, 1, 'B2');
xlswrite([resDir 'regionStats.xlsx'], rSum, 1, 'C2');
xlswrite([resDir 'regionStats.xlsx'], nSamples, 1, 'D2');
xlswrite([resDir 'regionStats.xlsx'], [{'Region'}, {'p-value'}, {'ranksum'}, ...
    {'number of samples'}], 1, 'A1');

%% calculate p-value of expression 
ontologyFile = 'C:/Users/amahfouz/SURFdrive/Data/ABA_adult_human_brain/Ontology.xlsx';
for i = 1 : length(donors)
%     [~,IXexp] = sort(zscore(gene_expression_donors{i}(probe_id,:),[],2),'descend');
    pvalue = 2*(1-normcdf(abs(zscore(gene_expression_donors{i}(probe_id,:),[],2)),0,1));
    for j = 1 : length(mainStructures)-1
            str_children = strSamples_Human(mainStructures{j}, ontologyFile);
            region_idx = ismember(donor_sample{i}.structure_acronym, str_children); 
            if sum(region_idx) ~= 0
                minP(j,i) = min(pvalue(region_idx));
                maxP(j,i) = max(pvalue(region_idx));
            else
                minP(j,i) = 1;
                maxP(j,i) = 1;
            end
%             ranks = IXexp(region_idx);
%             allRanks = setdiff(1:length(IXexp),ranks);
%             [p(j,i),~,stats] = ranksum(ranks, allRanks);
%             rSum(j,i) = stats.ranksum;
            regionSize(j,i) = sum(region_idx);
            clear region_idx;
    end
%     p_corrected(:,i) = multtest(p(:,i), 'method', 'holm');
end
resDir = 'C:\Ahmed\Work\Results\DMD\AdultHumanBrain\';
xlswrite([resDir 'regions_stats_pVal.xlsx'], mainStructures(1:end-1)', 1, 'A2');
xlswrite([resDir 'regions_stats_pVal.xlsx'], minP, 1, 'B2');
xlswrite([resDir 'regions_stats_pVal.xlsx'], [{'region'} donors], 1, 'A1');
% save region sizes
xlswrite([resDir 'regions_no_of_samples.xlsx'], mainStructures(1:end-1)', 1, 'A2');
xlswrite([resDir 'regions_no_of_samples.xlsx'], regionSize, 1, 'B2');
xlswrite([resDir 'regions_no_of_samples.xlsx'], sum(regionSize')', 1, 'H2');
xlswrite([resDir 'regions_no_of_samples.xlsx'], [{'region'} donors {'total number of samples'}], 1, 'A1');

figure, hold on
for i = 1 : length(donors)
    subplot(2,3,i), hist(zscore(gene_expression_donors{i}(probe_id,:),[],2),30)
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','b','EdgeColor','w')
    grid on
    xlabel('expression (z-score)')
    ylabel('frequency')
    title(['Donor ' donors{i}], 'FontWeight', 'bold')
end
hold off

%% differential expression between each pair of regions (only the 105 analyzed in all 6 brains)
addpath('/tudelft.net/staff-bulk/ewi/insy/DBL/amahfouz/MATLAB/Libraries')
dataDir = '/tudelft.net/staff-bulk/ewi/insy/DBL/sjoerdhuisman/ABA_human_brain_probegene/';
resDir = '/tudelft.net/staff-bulk/ewi/insy/DBL/amahfouz/MATLAB/Results/DMD/';
donors = {'10021','12876','14380','15496','15697','9861'};
% read probe info
T = readtable([dataDir 'probe_info_2014-11-11.csv']);
probe.gene_symbol = T.gene_symbol;
probe.entrez_id = T.entrez_id;
clear T;
% find unique samples across the 6 brains
load('/tudelft.net/staff-bulk/ewi/insy/DBL/amahfouz/MATLAB/Results/DMD/uniqueAllSamples.mat');
load('/tudelft.net/staff-bulk/ewi/insy/DBL/amahfouz/MATLAB/Results/DMD/donor_sample.mat');
load('/tudelft.net/staff-bulk/ewi/insy/DBL/amahfouz/MATLAB/Results/DMD/Region.mat');
% select regions analyzed in all 6 brains (105 regions)
IX_6 = find(cell2mat(Region.count) == 6);
% initialize vectors
DE.pVal = zeros(length(donors),numel(IX_6),numel(IX_6),numel(probe.entrez_id));
DE.tStat = zeros(length(donors),numel(IX_6),numel(IX_6),numel(probe.entrez_id));
DE.pVal_corrected = zeros(length(donors),numel(IX_6),numel(IX_6),numel(probe.entrez_id));
for D = 1 : length(donors)
    D
    % gene expression data
    MicroarrayExpression = csvread([dataDir 'gene_expr_normalized_microarray_donor' donors{D} '_2014-11-11.csv']);
    % sample data
    for s1 = 1 : numel(IX_6)
        idx1 = find(donor_sample{D}.structure_id == uniqueAllSamples(IX_6(s1)));
        for s2 = 1 : numel(IX_6)
            if s2 ~= s1
                idx2 = find(donor_sample{D}.structure_id == uniqueAllSamples(IX_6(s2)));
%                 % find differentially expressed genes (t-test)
%                 [DE.pVal(D,:), DE.tStat(D,:)] = mattest(MicroarrayExpression(:,idx1), MicroarrayExpression(:,idx2));
%                 % correct for multiple testing
%                 DE.pVal_corrected(D,:) = multtest(squeeze(DE.pVal(D,:)),'method','BH');
                % find differentially expressed genes (t-test)
                [DE.pVal(D,s1,s2,:), DE.tStat(D,s1,s2,:)] = mattest(MicroarrayExpression(:,idx1), MicroarrayExpression(:,idx2));
                % correct for multiple testing
                DE.pVal_corrected(D,s1,s2,:) = multtest(squeeze(DE.pVal(D,s1,s2,:)),'method','BH');
            end
        end
    end
end
save('/tudelft.net/staff-bulk/ewi/insy/DBL/amahfouz/MATLAB/Results/DMD/DE.mat','DE','-v7.3')

%% %% Plot the expression of the gene in the cortex
% find first level substructures of cortical substructures
ontologyFile = 'C:\Users\amahfouz\SURFdrive\Data\ABA_adult_human_brain\Ontology.xlsx';
structures = {'FL', 'Ins', 'LL', 'OL', 'PL', 'TL'};
for i = 1 : length(gene_expression_donors)
    tempExp = zscore(gene_expression_donors{i},[],2);
    X = 0;
    for j = 1 : length(structures)
        % extract the first level children
        level = 1;
        str_children = strSamples_Human(structures{j}, ontologyFile, level);
        for k = 1 : length(str_children)
            level = 0;
            str_children_2 = strSamples_Human(str_children{k}, ontologyFile, level);
            region_idx = ismember(sample{i}.structure_acronym, str_children_2);
            avgExpMat(i,X+k) = mean(tempExp(:,region_idx),2);
%             avgExpMat(:,i,j,k) = mean(tempExp(:,region_idx),2);
            regionSize(i,X+k) = sum(region_idx);
            label_vector{X+k} = [structures{j} '_' str_children{k}];
        end
        X = X + k;
    end
end

figure,
boxplot(avgExpMat, label_vector, ...
    'factorgap',5, ...
    'factorseparator',1, 'labelorientation', 'inline', 'plotstyle', 'compact')
grid on
ylabel('Expression', 'FontWeight', 'bold', 'FontSize', 15)
title([{geneName}; {['Average Expression - Probe: ' probeNames{5}]}], 'FontWeight', 'bold', 'FontSize', 15)

%% save data for plotting in R
% for each donor, read the data and zscore and concatinate
for D = 1 : length(donors)
    if D == 1
        concatExpression = zscore(gene_expression_donors{D},[],2);
        concatIDs = sample{D}.structure_id';
    else
        concatExpression = [concatExpression zscore(gene_expression_donors{D},[],2)];
        concatIDs = [concatIDs sample{D}.structure_id'];
    end    
end
T = table(concatExpression', concatIDs', 'VariableNames',{'exp','structure_id'});
writetable(T, 'C:\Users\amahfouz\SURFdrive\MATLAB_scripts\Human_brain_expression_visualization\DMD_6donors_concat_zscore.csv')

%% check PPI connections among a set of genes
% read PPI network data
PPI = readtable('C:\Users\amahfouz\SURFdrive\Data\StringNew_HPRD.txt',...
    'Delimiter','\t','ReadVariableNames',false);
% read the top 200 genes correlated with DMD in the adult brain
load([dataDir 'corrMat_combined.mat']);
[sortedCorrMat IX] = sort(corrMat_combined,2,'descend');
sortedGenes = probe.gene_symbol(IX);
top200_genes = sortedGenes(1:201); % 200 + DMD
% retrive all interactions between the 200 genes
ind1 = ismember(PPI.Var1, top200_genes);
ind2 = ismember(PPI.Var2, top200_genes);
top_interact = find((ind1+ind2) == 2); % only 6 interactions returned

%% build a coexpression netwok between the top 200 genes correalted with DMD in the adult brain
N = 25;
dataDir = '/tudelft.net/staff-bulk/ewi/insy/DBL/sjoerdhuisman/ABA_human_brain_probegene/';
resDir = '/tudelft.net/staff-bulk/ewi/insy/DBL/amahfouz/MATLAB/Results/DMD/';
% read the top 200 genes correlated with DMD in the adult brain
load([resDir 'corrMat_combined.mat']);
[sortedCorrMat IX] = sort(corrMat_combined,2,'descend');
sortedGenes = probe.gene_symbol(IX);
top200_genes = sortedGenes(1:N+1); % 200 + DMD
% generate the correlation matrix
% get the expression of genes from the AHBA (probes mapped to genes)
% read probe info
T = readtable([dataDir 'probe_info_2014-11-11.csv']);
probe.gene_symbol = T.gene_symbol;
probe.entrez_id = T.entrez_id;
gene_idx = find(ismember(probe.gene_symbol,top200_genes));
clear T;
% calfor i = 1 : length(donors)
    % gene expression data
    if i == 1
        MicroarrayExpression = csvread([dataDir 'gene_expr_normalized_microarray_donor' donors{i} '_2014-11-11.csv']);
        MicroarrayExpression = zscore(MicroarrayExpression,[],2);
    else
        tempExpMat = csvread([dataDir 'gene_expr_normalized_microarray_donor' donors{i} '_2014-11-11.csv']);
        tempExpMat = zscore(tempExpMat,[],2);
        MicroarrayExpression = [MicroarrayExpression, tempExpMat];
        clear tempExpMat;
    end
end
culate the correlation
[corrMat_combined pVal] = corr(MicroarrayExpression(gene_idx,:)', MicroarrayExpression(gene_idx,:)');
save([resDir 'corrMat_combined_top' num2str(N) '.mat'],'corrMat_combined','pVal','probe');

%% load the top N corelation matrix and save to excel
load([dataDir 'corrMat_combined_top' num2str(N) '.mat'])
corrMat_combined = corrMat_combined .* abs(eye(size(corrMat_combined,1))-1);
corr_top200 = squareform(corrMat_combined,'tovector');
% get the gene names
load([dataDir 'top200_gene_symbol.mat'])
count = 0;
for i = 1 : N+1%length(top200_gene_symbol)
    for j = i+1 : N+1%length(top200_gene_symbol)
        count = count + 1;
        top200_pair(count,1) = top200_gene_symbol(j);
        top200_pair(count,2) = top200_gene_symbol(i);
    end
end
T = table(top200_pair(:,1), top200_pair(:,2), corr_top200', 'VariableNames',{'Gene1','Gene2','correlation'});
writetable(T, [dataDir 'DMD_adult_top' num2str(N) '.xls'])
% figure, hist(corr_top200,100)

X = sort(corr_top200, 'descend');
Y = X(1:1000);
Y(end)




