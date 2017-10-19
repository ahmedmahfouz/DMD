
%%% 26 Jan 2015
%%% Read the gene expression data of DMD

%% Define directories
if ispc
%     atlasDir = 'C:/Ahmed/Work/Data/ALLEN BRAIN ATLASES/genes_matrix_csv_genCodeV10/';
    dataDir = 'C:\Users\amahfouz\SURFdrive\Projects\DMD\Data\'; 
end
geneName = 'DMD';

%% read data
% read rows meta data
[num,txt] = xlsread([atlasDir 'rows_metadata.xlsx']);
gene.gene_symbol = txt(2:end,4);
gene.ensembl_gene_id = txt(2:end,3);
gene.entrez_id = num(:,5);
clear num; clear txt;
save([dataDir 'gene.mat'],'gene');
% read columns meta data
[num,txt] = xlsread([atlasDir 'columns_metadata.xlsx']);
sample.donor_id = num(:,2);
sample.donor_name = txt(2:end,3);
sample.age = txt(2:end,4);
sample.gender = txt(2:end,5);
sample.structure_acronym = txt(2:end,7);
clear num; clear txt;
save([dataDir 'gene_sample.mat'],'sample');
% select only DMD data
gene_idx = find(strcmpi(gene.gene_symbol,geneName));
data = csvread([atlasDir 'expression_matrix.csv']);
dmd_data = data(gene_idx,2:end);
clear data;
save([dataDir 'dmd_gene_data.mat'],'dmd_data','gene_idx');

%% load data
load([dataDir 'dmd_gene_data.mat']);
load([dataDir 'gene.mat']);
load([dataDir 'gene_sample.mat']);

%% group ages into developmental stages and define structures
devStages = {'8 pcw','9 pcw','12 pcw',... % early fetal
    '13 pcw','16 pcw','17 pcw',... % early-mid fetal
    '19 pcw','21 pcw',... % late-mid fetal
    '24 pcw','25 pcw','26 pcw','35 pcw','37 pcw',... % late fetal
    '4 mos',... % Neonatal & early enfancy
    '10 mos',... % late enfancy
    '1 yrs','2 yrs','3 yrs','4 yrs',... % early childhood
    '8 yrs','11 yrs',... % middle & late childhood
    '13 yrs','15 yrs','18 yrs','19 yrs',... % adolescence
    '21 yrs','23 yrs','30 yrs','36 yrs','37 yrs',... % young adulthood
    '40 yrs'}; % middle adulthood
devStagesNum = [1,1,1,2,2,2,3,3,4,4,4,4,4,5,6,7,7,7,7,8,8,9,9,9,9,10,10,10,10,10,...
    11];
devStagesName = {'early_fetal','early_fetal','early_fetal',...
    'early_mid_fetal','early_mid_fetal','early_mid_fetal',...
    'late_mid_fetal','late_mid_fetal',...
    'late_fetal','late_fetal','late_fetal','late_fetal','late_fetal',...
    'Neonatal_and_early_enfancy',...
    'late_enfancy',...
    'early_childhood','early_childhood','early_childhood','early_childhood',...
    'middle_and_late_childhood','middle_and_late_childhood',...
    'adolescence','adolescence','adolescence','adolescence',...
    'young_adulthood','young_adulthood','young_adulthood','young_adulthood','young_adulthood',...
    'middle_adulthood'};
% assign developmental stage label to each sample
for i = 1 : length(sample.age)
    age_idx = find(strcmpi(devStages,sample.age{i})==1);
    devStageLabel(i) = devStagesNum(age_idx);
    devStageLabel_name{i} = devStagesName{age_idx};
end
structures = {'AMY','HIP','STR','MD','CBC',...
    'VFC','DFC','ITC','STC','MFC','OFC','M1C','IPC','A1C','V1C','S1C'};

%% plot the expression of the gene across all samples
% create sample labels
for i = 1 : length(sample.age)
    sampleLabel{i} = [sample.structure_acronym{i} ' _ ' sample.age{i}];
end
figure,
bar(dmd_data,'b','EdgeColor','w');
set(gca, 'XTickLabel', sampleLabel, 'XTick', 1:numel(sampleLabel), ...
    'xlim', [0 numel(sampleLabel)+1]);
ylabel('Expression (RPKM)', 'FontWeight', 'bold', 'FontSize', 15)
title([geneName ' Expression'], 'FontWeight', 'bold', 'FontSize', 15)
rotateXLabels(gca, 45)

%% plot different gene expression across age for different stuctures
structures = {'AMY','HIP','STR','MD','CBC',...
    'VFC','DFC','ITC','STC','MFC','OFC','M1C','IPC','A1C','V1C','S1C'};
ages = {'12 pcw','13 pcw','16 pcw','17 pcw','19 pcw','21 pcw','24 pcw','37 pcw',...
    '4 mos','1 yrs','2 yrs','3 yrs','4 yrs','8 yrs','11 yrs','13 yrs','18 yrs',...
    '19 yrs','21 yrs','23 yrs','30 yrs','36 yrs','37 yrs','40 yrs'};
% ages = unique(sample.age);
for i = 1 : length(structures)
    str_samples = find(strcmpi(sample.structure_acronym,structures{i})==1);
    for j = 1 : length(ages)
        str_age_samples = find(strcmpi(sample.age(str_samples),ages{j})==1);
        if isempty(str_age_samples)
            dmd_data_reshaped(i,j) = NaN;
        else
            dmd_data_reshaped(i,j) = nanmean(dmd_data(1,str_samples(str_age_samples)),2);
        end
    end
end
% average the cortical samples
dmd_ctx = squeeze(nanmean(dmd_data_reshaped(6:end,:),1));
% colors = jet(length(structures));
colors = lines(6);
figure, hold on
for j = 1 : 5
    plot(dmd_data_reshaped(j,:), ...
        'LineWidth', 2, 'Color', colors(j,:))
end
plot(dmd_ctx, ...
        'LineWidth', 2, 'Color', colors(6,:))
hold off
grid on
legend([structures(1:5) 'CTX'])
set(gca, 'XTickLabel', ages, 'XTick', 1:numel(ages), 'FontWeight', 'bold');
ylabel('Expression (RPKM)', 'FontWeight', 'bold', 'FontSize', 15)
title([geneName ' Expression Per Structure'], 'FontWeight', 'bold', 'FontSize', 15)

%% plot the expression of each structure across developmental stages
strLabelsModified = sample.structure_acronym;
strLabelsModified(find(ismember(sample.structure_acronym,structures(6:end))==1)) = {'CTX'};
strOfInterest = find(ismember(sample.structure_acronym,structures)==1);
figure,
boxplot(dmd_data(1,strOfInterest), {strLabelsModified(strOfInterest),devStageLabel_name(strOfInterest)}, ...
    'colorgroup',strLabelsModified(strOfInterest), 'factorgap',5, ...
    'factorseparator',1, 'labelorientation', 'inline', 'plotstyle', 'compact')
grid on
ylabel('Expression (RPKM)', 'FontWeight', 'bold', 'FontSize', 15)
title([geneName ' Expression'], 'FontWeight', 'bold', 'FontSize', 15)
figure,
boxplot(dmd_data(1,strOfInterest), {devStageLabel_name(strOfInterest),strLabelsModified(strOfInterest)}, ...
    'colorgroup',strLabelsModified(strOfInterest), 'factorgap',5, ...
    'factorseparator',1, 'labelorientation', 'inline', 'plotstyle', 'compact')
grid on
ylabel('Expression (RPKM)', 'FontWeight', 'bold', 'FontSize', 15)
title([geneName ' Expression'], 'FontWeight', 'bold', 'FontSize', 15)
% plot the expression across cortical structures
strOfInterest = find(ismember(sample.structure_acronym,structures(6:end))==1);
figure,
boxplot(dmd_data(1,strOfInterest), {sample.structure_acronym(strOfInterest),devStageLabel_name(strOfInterest)}, ...
    'colorgroup',sample.structure_acronym(strOfInterest), 'factorgap',5, ...
    'factorseparator',1, 'labelorientation', 'inline', 'plotstyle', 'compact')
grid on
ylabel('Expression (RPKM)', 'FontWeight', 'bold', 'FontSize', 15)
title([geneName ' Expression'], 'FontWeight', 'bold', 'FontSize', 15)
figure,
boxplot(dmd_data(1,strOfInterest), {devStageLabel_name(strOfInterest),sample.structure_acronym(strOfInterest)}, ...
    'colorgroup',sample.structure_acronym(strOfInterest), 'factorgap',5, ...
    'factorseparator',1, 'labelorientation', 'inline', 'plotstyle', 'compact')
grid on
ylabel('Expression (RPKM)', 'FontWeight', 'bold', 'FontSize', 15)
title([geneName ' Expression'], 'FontWeight', 'bold', 'FontSize', 15)

%% plot the expression of each structure across developmental stages (not BOXPLOT)
strLabelsModified = sample.structure_acronym;
strLabelsModified(find(ismember(sample.structure_acronym,structures(6:end))==1)) = {'CTX'};
structures_6 = {'AMY','HIP','STR','MD','CBC','CTX'};
for S = 1 : length(structures_6)
    str_idx = find(ismember(strLabelsModified,structures_6{S})==1);
    if ~isempty(str_idx)
        for T = 1 : max(devStageLabel)
            str_devStage_idx = str_idx(devStageLabel(str_idx) == T);
            if ~isempty(str_devStage_idx)
                expMean(S,T) = mean(dmd_data(1,str_devStage_idx));
                expMax(S,T) = max(dmd_data(1,str_devStage_idx)) - expMean(S,T);
                expMin(S,T) = expMean(S,T) - min(dmd_data(1,str_devStage_idx));
            else
                expMean(S,T) = NaN;
                expMax(S,T) = NaN;
                expMin(S,T) = NaN;
            end
        end
    else
        expMean(S,T) = 0;
        expMax(S,T) = 0;
        expMin(S,T) = 0;
    end
    clear str_idx;
end

for i = 1 : length(unique(devStagesNum))
    devS{i} = devStageLabel_name{find(devStageLabel == i, 1, 'first')};
end

C = jet(6);
figure, hold on
for i = 1 : size(expMean,1)
    h = errorbar([1:max(devStageLabel)], expMean(i,:), expMin(i,:), expMax(i,:));
    hc = get(h, 'Children');
    set(hc(1), 'lineWidth', 2)
    set(hc(1), 'color', C(i,:))
    set(hc(2), 'lineWidth', 2)
    set(hc(2), 'color', C(i,:))
end
hold off
grid on
ylabel('Expression (RPKM)', 'FontWeight', 'bold', 'FontSize', 15)
set(gca, 'XTick', 1:max(devStageLabel), 'XTickLabel', devS)
XYrotalabel

%% old plotting stuff
colors = lines(16);
figure, hold on
for i = 1 : length(structures)
    str_samples = find(strcmpi(sample.structure_acronym,structures{i})==1);
    if i == 1
        currData = dmd_data(1,str_samples);
        strLabels = devStageLabel(str_samples);
    else
        currData = [currData, dmd_data(1,str_samples)];
    end
%     currStrData = dmd_data(1,str_samples);
%     currStrLabels = devStageLabel(str_samples);
    subplot(5,1,i), 
    boxplot(currStrData,currStrLabels,'plotstyle','compact','colors',colors(i,:));
end
hold off

for i = 1 : length(structures)
    str_samples = find(strcmpi(sample.structure_acronym,structures{i})==1);
    for j = 1 : length(devStage)
        str_age_samples = find(ismember(sample.age(str_samples),devStage{j})==1);
        if isempty(str_age_samples)
            dmd_data_reshaped{i,j} = NaN;
        else
            dmd_data_reshaped{i,j} = dmd_data(1,str_samples(str_age_samples));
        end
    end
end
% average the cortical samples
dmd_ctx = squeeze(nanmean(dmd_data_reshaped(6:end,:),1));
% colors = jet(length(structures));
colors = lines(6);
figure, hold on
for j = 1 : 5
    plot(dmd_data_reshaped(j,:), ...
        'LineWidth', 2, 'Color', colors(j,:))
end
plot(dmd_ctx, ...
        'LineWidth', 2, 'Color', colors(6,:))
hold off
grid on
legend([structures(1:5) 'CTX'])
set(gca, 'XTickLabel', ages, 'XTick', 1:numel(ages), 'FontWeight', 'bold');
ylabel('Expression (RPKM)', 'FontWeight', 'bold', 'FontSize', 15)
title([geneName ' Expression Per Structure'], 'FontWeight', 'bold', 'FontSize', 15)

%% co-expression analysis
dataDir = '/tudelft.net/staff-bulk/ewi/insy/DBL/amahfouz/MATLAB/Data/genes_matrix_csv_genCodeV10/';
resDir = '/tudelft.net/staff-bulk/ewi/insy/DBL/amahfouz/MATLAB/Results/DMD/';geneName = 'DMD';
% read probe info
T = readtable([dataDir 'rows_metadata.csv']);
probe.gene_symbol = T.gene_symbol;
probe.entrez_id = T.entrez_id;
% remove pobes with no entrez_id
removedProbes = find(isnan(probe.entrez_id));
probe.gene_symbol(removedProbes) = [];
probe.entrez_id(removedProbes) = [];
gene_idx = find(strcmpi(probe.gene_symbol,geneName));
clear T;
% read expression data
data = csvread([dataDir 'expression_matrix.csv']);
data(removedProbes,:) = [];
data(:,1) = [];
% calculate the correlation
[corrMat pVal] = corr(data(gene_idx,:)', data');
save([resDir 'corrMat_brainSpan.mat'],'corrMat','pVal','probe');

%% Analyze the combined donor correlation list
dataDir = 'C:/Ahmed/Work/Data/DMD/';
resultsDir = 'C:/Ahmed/Work/Results/DMD/';
load([dataDir 'corrMat_brainSpan.mat']);
% remove NaN corelations
removedRows = find(isnan(corrMat));
corrMat(removedRows) = [];
probe.gene_symbol(removedRows) = [];
probe.entrez_id(removedRows) = [];
pVal(removedRows) = [];
% sort correlations
[sortedCorrMat IX] = sort(corrMat,2,'descend');
% save to excel
xlswrite([resultsDir geneName '_BrainSpan_coexpressed_genes.xlsx'], probe.gene_symbol(IX(2:end)), 1, 'A1');
xlswrite([resultsDir geneName '_BrainSpan_coexpressed_genes.xlsx'], probe.entrez_id(IX(2:end)), 1, 'B1');
xlswrite([resultsDir geneName '_BrainSpan_coexpressed_genes.xlsx'], sortedCorrMat(2:end)', 1, 'C1');
[sortedCorrMat IX] = sort(corrMat,2);
xlswrite([resultsDir geneName '_BrainSpan_coexpressed_genes_neg.xlsx'], probe.gene_symbol(IX(2:end)), 1, 'A1');
xlswrite([resultsDir geneName '_BrainSpan_coexpressed_genes_neg.xlsx'], probe.entrez_id(IX(2:end)), 1, 'B1');
xlswrite([resultsDir geneName '_BrainSpan_coexpressed_genes_neg.xlsx'], sortedCorrMat(2:end)', 1, 'C1');
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
set(gca,'XTick',[],'XTickLabel',[],'xlim',[0 numel(corrMat)])
title(['Correlation to ' geneName], 'FontWeight', 'bold', 'FontSize', 15)
subplot(2,1,2), bar(-log10(pVal(IX(2:end))),'b','EdgeColor','w'), grid on
ylabel('-log_1_0 (p-value)', 'FontWeight', 'bold', 'FontSize', 15)
xlabel('Genes sorted on corrleation')
set(gca,'XTick',[],'XTickLabel',[],'xlim',[0 numel(corrMat)])
hold off

%% save data to xls
xlswrite('dmd_gene_BrainSpan.xlsx', dmd_data, 1, 'B4')
xlswrite('dmd_gene_BrainSpan.xlsx', sample.donor_name', 1, 'B1');
xlswrite('dmd_gene_BrainSpan.xlsx', sample.age', 1, 'B2');
xlswrite('dmd_gene_BrainSpan.xlsx', sample.structure_acronym', 1, 'B3');
xlswrite('dmd_gene_BrainSpan.xlsx', [{'Name'},{'Age'},{'Structure'},...
    {'DMD'}]', 1, 'A1');



