%%% 19 Aug 2015
%%% Read the exon expression (microarray) data of DMD

%% Define directories
if ispc
    atlasDir = 'C:/Users/amahfouz/SURFdrive/Projects/DMD/Data/exon_array_matrix_csv/'; 
    dataDir = 'C:/Users/amahfouz/SURFdrive/Projects/DMD/Data/';
end

%% read data
dmdGene = 'DMD';
% read rows meta data
T = readtable([atlasDir 'rows_metadata.csv']);
transcript.gene_symbol = T.gene_symbol;
transcript.ensembl_gene_id = T.ensembl_gene_id;
transcript.entrez_id = T.entrez_id;
transcript.start = T.start;
transcript.end = T.xEnd;
clear num; clear txt;
save([dataDir 'transcript_mircoarray.mat'],'transcript');
% read columns meta data
T = readtable([atlasDir 'columns_metadata.csv']);
sample.donor_id = T.donor_id;
sample.donor_name = T.donor_name;
sample.age = T.age;
sample.gender = T.gender;
sample.structure_acronym = T.structure_acronym;
clear num; clear txt;
save([dataDir 'exon_sample_microarray.mat'],'sample');
% select only DMD data
gene_idx = find(strcmpi(transcript.gene_symbol,dmdGene));
data = csvread([atlasDir 'expression_matrix.csv']);
dmd_data = data(gene_idx,2:end);
clear data;
save([dataDir 'dmd_exon_data_miicroarray.mat'],'dmd_data','gene_idx');

%% load data
load([dataDir 'dmd_exon_data_microarray.mat']);
load([dataDir 'transcript_microarray.mat']);
load([dataDir 'exon_sample_microarray.mat']);

%% save data to xlsx
xlswrite([dataDir 'dmd_exon_data.xlsx'], sample.donor_name', 1, 'C1');
xlswrite([dataDir 'dmd_exon_data.xlsx'], sample.gender', 1, 'C2');
xlswrite([dataDir 'dmd_exon_data.xlsx'], sample.age', 1, 'C3');
xlswrite([dataDir 'dmd_exon_data.xlsx'], sample.structure_acronym', 1, 'C4');
xlswrite([dataDir 'dmd_exon_data.xlsx'], dmd_data, 1, 'C5');
transcript_info(1,:) = transcript.start(gene_idx);
transcript_info(2,:) = transcript.end(gene_idx);
xlswrite([dataDir 'dmd_exon_data.xlsx'], transcript_info', 1, 'A5');
xlswrite([dataDir 'dmd_exon_data.xlsx'], [{'transcript_start'},{'transcript_end'}], 1, 'A4');

%% developmental stages
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
    'Neonatal/early_enfancy',...
    'late_enfancy',...
    'early_childhood','early_childhood','early_childhood','early_childhood',...
    'middle/late_childhood','middle/late_childhood',...
    'adolescence','adolescence','adolescence','adolescence',...
    'young_adulthood','young_adulthood','young_adulthood','young_adulthood','young_adulthood',...
    'middle_adulthood'};

%% select transcripts of each isoform
transcript_ranges = [31150000,31300000; 31600000,32000000; 32200000,33000000];
for i = 1 : size(transcript_ranges,1)
    % find transcripts within the curent range
    R1 = find(transcript.start(gene_idx) > transcript_ranges(i,1));
    R2 = find(transcript.end(gene_idx) < transcript_ranges(i,2));
    isoform_transcripts{i} = intersect(R1,R2);
    trans_exp(i,:) = mean(dmd_data(isoform_transcripts{i},:));
    L{i} = ['exons: ' num2str(transcript_ranges(i,1)) '-' num2str(transcript_ranges(i,2))];
    clear R1; clear R2;
end

%% plot the average expression of exon groups
C = [1,0,0;0,1,0;0,0,1];
for i = 1 : length(devStages)
    xTicks(i) = find(ismember(sample.age, devStages{i})==1, 1, 'first');
end
figure, hold on
for i = 1 : length(isoform_transcripts)
    plot(trans_exp(i,:),'color',C(i,:),'linewidth',2);
end
hold off
grid on
ylabel('Expression (RPKM)', 'FontWeight', 'bold', 'FontSize', 15);
set(gca, 'XTickLabel', devStages, 'XTick', xTicks, 'xlim', [0 size(trans_exp,2)]);
rotateXLabels(gca, 45)
legend(L);

%% plot an image of the probe 
% create sample labels
for i = 1 : length(sample.age)
    sampleLabel{i} = [sample.structure_acronym{i} ' _ ' sample.age{i}];
end
for i = 1 : size(dmd_data,1)
    transcriptLabel{i} = [num2str(i) ':' num2str(transcript.start(gene_idx(i))) ' - ' num2str(transcript.end(gene_idx(i)))];
end

for i = 1 : length(devStages)
    xTicks(i) = find(ismember(sample.age, devStages{i})==1, 1, 'first');
end

% figure,
% imagesc(log10(dmd_data+1)), colormap('hot'), colorbar
% set(gca, 'XTickLabel', sampleLabel, 'XTick', 1:numel(sampleLabel));
% set(gca, 'YTickLabel', transcriptLabel, 'YTick', 1:numel(transcriptLabel));
% rotateXLabels(gca, 45)

figure,
imagesc(fliplr(log10(dmd_data+1)')), colormap('hot'), colorbar
set(gca, 'XTickLabel', fliplr(transcriptLabel), 'XTick', 1:numel(transcriptLabel));
set(gca, 'YTickLabel', devStages, 'YTick', xTicks);
% for i = 1 : length(isoform_transcripts)
%     line([size(dmd_data,1)-isoform_transcripts{i}(end) size(dmd_data,1)-isoform_transcripts{i}(end)],...
%         [0 524], 'Color', 'w', 'LineWidth', 5);
% end
% XYrotalabel
% set(gca, 'YTickLabel', sample.age, 'YTick', 1:numel(sampleLabel), 'FontSize', 5);
rotateXLabels(gca, 45)

figure, 
boxplot(fliplr(log2(dmd_data+1)'))
ylabel('Expression (RPKM)', 'FontWeight', 'bold', 'FontSize', 15)
set(gca, 'XTickLabel', fliplr(1:length(transcriptLabel)), 'XTick', 1:numel(transcriptLabel));
for i = 1 : length(isoform_transcripts)
    line([size(dmd_data,1)-isoform_transcripts{i}(end) size(dmd_data,1)-isoform_transcripts{i}(end)],...
        [0 524], 'Color', 'black', 'LineWidth', 5);
end
% ylim([-0.1 6])
XYrotalabel
rotateXLabels(gca, 45)

%% exclude some exons
exon_exclude = [92,91,89,88,87,46,40,22];
dmd_data_updated = dmd_data;
dmd_data_updated(exon_exclude,:) = [];
transcriptLabel_updated = transcriptLabel;
transcriptLabel_updated(exon_exclude) = [];

custom_map = csvread('redblue_256_rgb.txt');
custom_map = custom_map / 255;
custom_map = flipud(custom_map);
figure,
imagesc(fliplr(log2(dmd_data_updated+1)')), colormap(custom_map), colorbar
imagesc(fliplr(dmd_data_updated')), colormap(custom_map), colorbar
X = fliplr(log10(dmd_data_updated+1)');
hold on;
% for i = 1:size(X,1)
%    plot([.5,size(X,2)+.5],[i-.5,i-.5],'k-');
% end
for i = 1:size(X,2)
    plot([i-.5,i-.5],[.5,size(X,1)+.5],'k-');
end
hold off
set(gca, 'XTickLabel', fliplr(transcriptLabel_updated), 'XTick', 1:numel(transcriptLabel_updated));
set(gca, 'YTickLabel', devStages, 'YTick', xTicks);
rotateXLabels(gca, 45)

%% probe selection
exon_select = [19,41,86];
dmd_sub = dmd_data(exon_select,:);
xlswrite('individual_exons.xlsx', dmd_sub, 1, 'B4')
xlswrite('individual_exons.xlsx', sample.donor_name', 1, 'B1');
xlswrite('individual_exons.xlsx', sample.age', 1, 'B2');
xlswrite('individual_exons.xlsx', sample.structure_acronym', 1, 'B3');
xlswrite('individual_exons.xlsx', [{'Name'},{'Age'},{'Structure'},...
    {'probe#19'},{'probe#41'},{'probe#86'}]', 1, 'A1');

%% Plot the expression of selected exons across development for different structures
exon_select = [19,41,86];
dmd_sub = dmd_data(exon_select,:);
% dmd_sub = log2(sum(dmd_data)+1);
isoforms = {'Dp71-Dp40', 'Dp140', 'Dp427'};

figure, hold on
% cortex
subplot(2,2,1), plot(dmd_sub(:,ismember(sample.structure_acronym,{'VFC','DFC','ITC','STC','MFC','OFC','M1C','IPC','A1C','V1C','S1C'}))')
legend(isoforms)
% cerebellum
subplot(2,2,2), plot(dmd_sub(:,ismember(sample.structure_acronym,{'CB','CBC'}))')
legend(isoforms)
% hippocampus
subplot(2,2,3), plot(dmd_sub(:,ismember(sample.structure_acronym,{'HIP'}))')
legend(isoforms)
% Amygdala
subplot(2,2,4), plot(dmd_sub(:,ismember(sample.structure_acronym,{'AMY'}))')
legend(isoforms)
hold off

figure, plot(sum(dmd_data))

%% cluster the DMD transcripts
% clustergram((log10(dmd_data+1)'))

%% group ages into developmental stages and define structures
% assign developmental stage label to each sample
for i = 1 : length(sample.age)
    age_idx = find(strcmpi(devStages,sample.age{i})==1);
    devStageLabel(i) = devStagesNum(age_idx);
    devStageLabel_name{i} = devStagesName{age_idx};
end
structures = {'AMY','HIP','STR','MD','CBC',...
    'VFC','DFC','ITC','STC','MFC','OFC','M1C','IPC','A1C','V1C','S1C'};

%% plot different transcripts across age for different stuctures 
strLabelsModified = sample.structure_acronym;
strLabelsModified(find(ismember(sample.structure_acronym,structures(6:end))==1)) = {'CTX'};
strOfInterest = find(ismember(sample.structure_acronym,structures)==1);
for i = 1 : length(isoform_transcripts)
    figure,
    boxplot(nanmean(dmd_data(isoform_transcripts{i},strOfInterest)), ...
        {strLabelsModified(strOfInterest),devStageLabel_name(strOfInterest)}, ...
        'colorgroup',strLabelsModified(strOfInterest), 'factorgap',5, ...
        'factorseparator',1, 'labelorientation', 'inline', 'plotstyle', 'compact')
    grid on
    ylabel('Expression (RPKM)', 'FontWeight', 'bold', 'FontSize', 15)
    ylim([0 20])
    title(['Probes: ' num2str(transcript_ranges(i,1)) ' - ' num2str(transcript_ranges(i,2))],...
        'FontWeight', 'bold', 'FontSize', 15);
    figure,
    boxplot(nanmean(dmd_data(isoform_transcripts{i},strOfInterest)), ...
        {devStageLabel_name(strOfInterest),strLabelsModified(strOfInterest)}, ...
        'colorgroup',strLabelsModified(strOfInterest), 'factorgap',5, ...
        'factorseparator',1, 'labelorientation', 'inline', 'plotstyle', 'compact')
    grid on
    ylabel('Expression (RPKM)', 'FontWeight', 'bold', 'FontSize', 15)
    ylim([0 20])
    title(['Probes: ' num2str(transcript_ranges(i,1)) ' - ' num2str(transcript_ranges(i,2))],...
        'FontWeight', 'bold', 'FontSize', 15);
end
% plot the expression across cortical structures
strOfInterest = find(ismember(sample.structure_acronym,structures(6:end))==1);
for i = 1 : length(isoform_transcripts)
    figure,
%     boxplot(nanmean(dmd_data(isoform_transcripts{i},strOfInterest)), ...
%         {sample.structure_acronym(strOfInterest),devStageLabel_name(strOfInterest)}, ...
%         'colorgroup',sample.structure_acronym(strOfInterest), 'factorgap',5, ...
%         'factorseparator',1, 'labelorientation', 'inline')
    boxplot(nanmean(dmd_data(isoform_transcripts{i},strOfInterest)), ...
        {devStageLabel_name(strOfInterest),sample.structure_acronym(strOfInterest)}, ...
        'colorgroup',sample.structure_acronym(strOfInterest), 'factorgap',5, ...
        'factorseparator',1, 'labelorientation', 'inline')
    grid on
    ylabel('Expression (RPKM)', 'FontWeight', 'bold', 'FontSize', 15)
    ylim([0 20])
    title(['Probes: ' num2str(transcript_ranges(i,1)) ' - ' num2str(transcript_ranges(i,2))],...
        'FontWeight', 'bold', 'FontSize', 15);
end

%% OLD plot different transcripts across age for different stuctures 
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
            dmd_data_reshaped(:,i,j) = NaN;
        else
            dmd_data_reshaped(:,i,j) = nanmean(dmd_data(:,str_samples(str_age_samples)),2);
        end
    end
end
% average the cortical samples
dmd_ctx = squeeze(nanmean(dmd_data_reshaped(:,6:end,:),2));
% colors = jet(length(structures));
colors = lines(6);
for i = 1 : length(isoform_transcripts)
    figure, hold on
    for j = 1 : 5
        plot(nanmean(squeeze(dmd_data_reshaped(isoform_transcripts{i},j,:)))', ...
            'LineWidth', 2, 'Color', colors(j,:))
    end
    plot(nanmean(dmd_ctx(isoform_transcripts{i},:)), ...
            'LineWidth', 2, 'Color', colors(6,:))
    hold off
    grid on
    legend([structures(1:5) 'CTX'])
    set(gca, 'XTickLabel', ages, 'XTick', 1:numel(ages), 'FontWeight', 'bold');
    ylim([0 18])
    ylabel('Expression (RPKM)', 'FontWeight', 'bold', 'FontSize', 15)
    title(['Probes: ' num2str(transcript_ranges(i,1)) ' - ' num2str(transcript_ranges(i,2))],...
        'FontWeight', 'bold', 'FontSize', 15);
end

%% plot the expression of each isofom's transcripts
for i = 1 : 3
    figure, plot(dmd_data(isoform_transcripts{i},:)');
end

%% co-expression analysis
dataDir = '/tudelft.net/staff-bulk/ewi/insy/DBL/amahfouz/MATLAB/Data/exons_matrix_csv/';
resDir = '/tudelft.net/staff-bulk/ewi/insy/DBL/amahfouz/MATLAB/Results/DMD/';
geneName = 'DMD';
% read probe info
T = readtable([dataDir 'rows_metadata.xlsx']);
transcript.gene_symbol = T.gene_symbol;
transcript.entrez_id = T.entrez_id;
transcript.start = T.start;
transcript.end = T.xEnd;
clear T;
% remove pobes with no entrez_id
removedProbes = find(isnan(transcript.entrez_id));
transcript.gene_symbol(removedProbes) = [];
transcript.entrez_id(removedProbes) = [];
transcript.start(removedProbes) = [];
transcript.end(removedProbes) = [];
gene_idx = find(strcmpi(transcript.gene_symbol,geneName));
% read expression data
data = csvread([dataDir 'expression_matrix.csv']);
data(removedProbes,:) = [];
data(:,1) = [];
% create expression profile of isoforms
transcript_ranges = [31150000,31300000; 31600000,32000000; 32200000,33000000];
for i = 1 : size(transcript_ranges,1)
    % find transcripts within the curent range
    R1 = find(transcript.start(gene_idx) > transcript_ranges(i,1));
    R2 = find(transcript.end(gene_idx) < transcript_ranges(i,2));
    isoform_transcripts{i} = intersect(R1,R2);
    clear R1; clear R2;
end
trans1_exp = mean(data(gene_idx(isoform_transcripts{1}),:));
trans2_exp = mean(data(gene_idx(isoform_transcripts{2}),:));
trans3_exp = mean(data(gene_idx(isoform_transcripts{3}),:));
% calculate the correlation
[corrMat pVal] = corr([trans1_exp;trans2_exp;trans3_exp]', data');
save([resDir 'corrMat_exonArray.mat'],'corrMat','pVal','transcript');

%% coexpression of individual exons
exon_select = [19,41,86];
dataDir = '/tudelft.net/staff-bulk/ewi/insy/DBL/amahfouz/MATLAB/Data/exons_matrix_csv/';
resDir = '/tudelft.net/staff-bulk/ewi/insy/DBL/amahfouz/MATLAB/Results/DMD/';
geneName = 'DMD';
% read probe info
T = readtable([dataDir 'rows_metadata.xlsx']);
transcript.gene_symbol = T.gene_symbol;
transcript.entrez_id = T.entrez_id;
transcript.start = T.start;
transcript.end = T.xEnd;
clear T;
% remove pobes with no entrez_id
removedProbes = find(isnan(transcript.entrez_id));
transcript.gene_symbol(removedProbes) = [];
transcript.entrez_id(removedProbes) = [];
transcript.start(removedProbes) = [];
transcript.end(removedProbes) = [];
gene_idx = find(strcmpi(transcript.gene_symbol,geneName));
% read expression data
data = csvread([dataDir 'expression_matrix.csv']);
data(removedProbes,:) = [];
data(:,1) = [];
% create expression profile of isoforms
trans1_exp = (data(gene_idx(exon_select(1)),:));
trans2_exp = (data(gene_idx(exon_select(2)),:));
trans3_exp = (data(gene_idx(exon_select(3)),:));
% calculate the correlation
[corrMat pVal] = corr([trans1_exp;trans2_exp;trans3_exp]', data');
save([resDir 'corrMat_exonArray_individual_exons.mat'],'corrMat','pVal','transcript');

%% Analyze the combined donor correlation list
geneName = 'DMD';
transcript_ranges = [31150000,31300000; 31600000,32000000; 32200000,33000000];
dataDir = 'C:/Users/amahfouz/SURFdrive/Projects/DMD/Data/';
resultsDir = 'C:/Users/amahfouz/SURFdrive/Projects/DMD/Results/';
% load([dataDir 'corrMat_exonArray.mat']);
load([dataDir 'corrMat_exonArray_individual_exons.mat']);
removedRows = logical(prod(double(isnan(corrMat))));
corrMat(:,removedRows) = [];
transcript.gene_symbol(removedRows) = [];
transcript.entrez_id(removedRows) = [];
transcript.start(removedRows) = [];
transcript.end(removedRows) = [];
pVal(:,removedRows) = [];
% f = figure; hold on
for i = 1 : size(corrMat,1)
    % sort the data
    [sortedCorrMat IX] = sort(corrMat(i,:),2,'descend');
    rankedGenes = transcript.gene_symbol(IX);
    [b,m,~] = unique(rankedGenes, 'stable');
    [~,IX2] = sort(m);
    sortedG = b(IX2);
    sortedCorr = sortedCorrMat(m(IX2));
    sortedPval = pVal(i,IX);
    sortedPval = sortedPval(m(IX2));
    sotredEntrez = transcript.entrez_id(IX);
    sotredEntrez = sotredEntrez(m(IX2));
    % save to excel
    xlswrite([resultsDir geneName '_exonArray_coexpressed_genes.xlsx'], sortedG(2:end), i, 'A1');
    xlswrite([resultsDir geneName '_exonArray_coexpressed_genes.xlsx'], sotredEntrez(2:end), i, 'B1');
    xlswrite([resultsDir geneName '_exonArray_coexpressed_genes.xlsx'], sortedCorr(2:end)', i, 'C1');
    clear IXl clear IX2; clear b; clear m;
%     % flip the list for negaive correlations
%     [sortedCorrMat IX] = sort(corrMat(i,:),2);
%     rankedGenes = transcript.gene_symbol(IX);
%     [b,m,~] = unique(rankedGenes, 'first');
%     [~,IX2] = sort(m);
%     sortedG = b(IX2);
%     sortedCorr = sortedCorrMat(m(IX2));
%     sotredEntrez = transcript.entrez_id(IX);
%     sotredEntrez = sotredEntrez(m(IX2));
%     xlswrite([resultsDir geneName '_exonArray_coexpressed_genes.xlsx'], sortedG(1:end-1), i, 'A1');
%     xlswrite([resultsDir geneName '_exonArray_coexpressed_genes.xlsx'], sotredEntrez(1:end-1), i, 'B1');
%     xlswrite([resultsDir geneName '_exonArray_coexpressed_genes.xlsx'], sortedCorr(1:end-1)', i, 'C1');
%     clear IXl; clear IX2; clear b; clear m;
%     % plot correlation and p-value    
%     subplot(2,3,i), hold on
%     line([0 length(sortedCorr)-1], [0 0],'LineStyle','--','Color',[0.5,0.5,0.5],'LineWidth',2)
%     plot(2:201,sortedCorr(2:201),'LineWidth',3,'Color','r'),
%     plot(202:numel(sortedCorr)-200,sortedCorr(202:end-200),'LineWidth',3),
%     plot(numel(sortedCorr)-199:numel(sortedCorr),sortedCorr(end-199:end),'LineWidth',3,'Color','r'),
%     grid on, hold off
%     ylabel('Correlation', 'FontWeight', 'bold', 'FontSize', 15)
%     xlabel('Genes sorted on corrleation')
%     set(gca,'XTick',[],'XTickLabel',[],'xlim',[0 numel(sortedCorr)],'ylim',[-1 1])
%     title([{['Correlation to ' geneName]},...
%         {['Probes: ' num2str(transcript_ranges(i,1)) ' - ' num2str(transcript_ranges(i,2))]}],...
%         'FontWeight', 'bold', 'FontSize', 15)
%     subplot(2,3,i+3), bar(-log10(sortedPval(2:end)),'b','EdgeColor','w'), grid on
%     ylabel('-log_1_0 (p-value)', 'FontWeight', 'bold', 'FontSize', 15)
%     xlabel('Genes sorted on corrleation')
%     set(gca,'XTick',[],'XTickLabel',[],'xlim',[0 numel(sortedCorr)],'ylim',[0 200])
end
% hold off
% saveas(f, [resultsDir 'corelationPlot_BrainSpan_exons_individual_exons.fig']);
% saveas(f, [resultsDir 'corelationPlot_BrainSpan_exons_individual_exons.png']);

%% build a coexpression netwok between the top 25 genes correalted with DMD in the adult brain
N = 25;
exon_select = [19,41,86];
dataDir = '/tudelft.net/staff-bulk/ewi/insy/DBL/amahfouz/MATLAB/Data/exons_matrix_csv/';
resDir = '/tudelft.net/staff-bulk/ewi/insy/DBL/amahfouz/MATLAB/Results/DMD/';
geneName = 'DMD';
% read probe info
T = readtable([dataDir 'rows_metadata.xlsx']);
origTranscript.gene_symbol = T.gene_symbol;
origTranscript.entrez_id = T.entrez_id;
origTranscript.start = T.start;
origTranscript.end = T.xEnd;
clear T;
% remove pobes with no entrez_id
removedProbes = find(isnan(origTranscript.entrez_id));
origTranscript.gene_symbol(removedProbes) = [];
origTranscript.entrez_id(removedProbes) = [];
origTranscript.start(removedProbes) = [];
origTranscript.end(removedProbes) = [];
% read expression data
data = csvread([dataDir 'expression_matrix.csv']);
data(removedProbes,:) = [];
data(:,1) = [];
% read the coexpression results and clean the data
load([resDir 'corrMat_exonArray_individual_exons.mat']);
removedRows = logical(prod(double(isnan(corrMat))));
corrMat(:,removedRows) = [];
transcript.gene_symbol(removedRows) = [];
transcript.entrez_id(removedRows) = [];
transcript.start(removedRows) = [];
transcript.end(removedRows) = [];
pVal(:,removedRows) = [];
data(removedRows,:) = [];
% loop on the different isoform-specific exons
for i = 1 : size(corrMat,1)
    % sort the data
    [sortedCorrMat IX] = sort(corrMat(i,:),2,'descend');
    rankedGenes = transcript.gene_symbol(IX);
    [b,m,~] = unique(rankedGenes, 'stable');
    [~,IX2] = sort(m);
    sortedG = b(IX2);
    sortedData = data(IX,:);
    sortedData = sortedData(m(IX2),:);
    % select the top N genes and calculate the correlation between them
    topN_genes{i} = sortedG(1:N+1); % N + DMD
    [corrMat_topN(:,:,i),~] = corr(sortedData(1:N+1,:)');
end
save([resDir 'corrMat_exonArray_individual_exons_topN.mat'],'corrMat_topN','topN_genes');

%% load the topN corelations and save to excel
load([dataDir 'corrMat_exonArray_individual_exons_topN.mat'])
for k = 1 : size(corrMat_topN,3)
    corrMat_curr = corrMat_topN(:,:,k) .* abs(eye(size(corrMat_topN,1))-1);
    corr_topN = squareform(corrMat_curr,'tovector');
    count = 0;
    for i = 1 : length(topN_genes{k})
        for j = i+1 : length(topN_genes{k})
            count = count + 1;
            top200_pair(count,1) = topN_genes{k}(j);
            top200_pair(count,2) = topN_genes{k}(i);
        end
    end
    T = table(top200_pair(:,1), top200_pair(:,2), corr_topN', 'VariableNames',{'Gene1','Gene2','correlation'});
    writetable(T, [dataDir 'DMD_exons_individual_topN_' num2str(k) '.xls']);
end

%% read the GOElite results and plot a graph
resDir = 'C:\Users\amahfouz\SURFdrive\Projects\DMD\Results\BrainSpan_exons\GO_elite\';
T_Dp71_Dp40 = readtable([resDir 'pruned-results_combination_elite_Dp71_Dp40.txt'],'Delimiter','\t');
T_Dp140 = readtable([resDir 'pruned-results_combination_elite_Dp140.txt'],'Delimiter','\t');
T_Dp427 = readtable([resDir 'pruned-results_combination_elite_Dp427.txt'],'Delimiter','\t');
tableList = {T_Dp71_Dp40, T_Dp140, T_Dp427};
for ii = 1 : length(tableList)
    currT = tableList{ii};
    % select GO terms related to "biological_processes" and "molecular_function"
    GOtype_idx = find(ismember(currT.GOType,{'biological_process','molecular_function'})==1);
    % select GO terms with a Z-score > 4
    zScore_idx = find(currT.ZScore >= 4);
    % select the GO terms fulfilling bith criteria
    selection_idx = intersect(GOtype_idx, zScore_idx);
    go_selected = currT.GOName_GOID_(selection_idx);
    
    % plot the GO terms enrichment
    f = figure;
    hold on
    set(f,'Position',[200, 200, 1024, 300])
    barh(flipud(currT.ZScore(selection_idx)),'FaceColor',[0 0 153/255],'EdgeColor',[0 0 153/255])
    set(gca,'YTick',1:length(selection_idx),'YTickLabel',flipud(go_selected))
    ylim([0 9])
    xlim([0 12])
    line([4 4],[0.5 length(selection_idx)+0.5],'Color','r','LineWidth',2,'LineStyle','--');
    grid on
    grid minor
    hold off
end







