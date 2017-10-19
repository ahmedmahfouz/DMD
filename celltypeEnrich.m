%%% 4 March 2015
%%% test gene sets for enrichment of cell-type specific genes

%% read cell-type marker data
dataDir = 'C:\Ahmed\Work\Data\CellTypeSpecificGenes\';
% Cahoy_2008: oligodendrocytes, astrocytes, neurons 
fe = 10; % fold enrichment threshold
% astrocytes
[num txt] = xlsread([dataDir 'Cahoy2008\340_Cahoy_S_Table_S4_AvX_2007-09-11.xls']);
cellTypeMarker.data{1} = txt(3:end,3);
cellTypeMarker.name{1} = 'Cahoy2008_astro';
foldEnrich = num(:,4);
belowThresh = find(foldEnrich < fe);
cellTypeMarker.data{1}(belowThresh) = [];
% oligodendrocytes
[num txt] = xlsread([dataDir 'Cahoy2008\350_Cahoy_S_Table_S5_OvX_2007-09-11.xls']);
cellTypeMarker.data{2} = txt(3:end,3);
cellTypeMarker.name{2} = 'Cahoy2008_oligo';
foldEnrich = num(:,4);
belowThresh = find(foldEnrich < fe);
cellTypeMarker.data{2}(belowThresh) = [];
% neurons
[num txt] = xlsread([dataDir 'Cahoy2008\360_Cahoy_S_Table_S6_NvX_2007-09-11.xls']);
cellTypeMarker.data{3} = txt(3:end,3);
cellTypeMarker.name{3} = 'Cahoy2008_neuro';
foldEnrich = num(:,4);
belowThresh = find(foldEnrich < fe);
cellTypeMarker.data{3}(belowThresh) = [];

% Lein_2007: oligodendrocytes, astrocytes, neurons
[num txt] = xlsread([dataDir 'Lein_Nature2007_TableS1.xls']);
cellTypeMarker.data{4} = txt(2:74,3);
cellTypeMarker.name{4} = 'Lein2007_neuro';
cellTypeMarker.data{5} = txt(2:49,5);
cellTypeMarker.name{5} = 'Lein2007_astro';
cellTypeMarker.data{6} = txt(2:80,7);
cellTypeMarker.name{6} = 'Lein2007_oligo';

% Bachoo_2004: astrocytes
[num txt] = xlsread([dataDir 'Bachoo_PNAS2004_TableS7.xls']);
cellTypeMarker.data{7} = unique(txt(2:end,4));
cellTypeMarker.name{7} = 'Bachoo2004_astro';

% Rong_2004: Purkinje neurons
[num txt] = xlsread([dataDir 'Rong_MolBrainRes2004_Table3.xls']);
cellTypeMarker.data{8} = unique(txt(4:end,1));
cellTypeMarker.name{8} = 'Rong_2004_purkinje';

% Hawrylycz_2012: post-synaptic density (psd)
[num txt] = xlsread([dataDir 'Hawrylycz_Nature2012_TableS3.xls']);
cellTypeMarker.data{9} = txt(2:end,2);
cellTypeMarker.name{9} = 'Hawrylycz_2012_psd';
cellTypeMarker.data{9}(74) = []; % remove a gap

% Nagata_2004: Ischemia and reperfusion in hippocampus
[num txt] = xlsread([dataDir 'Nagata_MolBrainRes2004_Table1.xls']);
cellTypeMarker.data{10} = txt(4:end,2);
cellTypeMarker.name{10} = 'Nagata_2004_isc_reperf_hip';

% Morciano_2005: synaptic vesicles & presynaptic membrane compartment

% save cell-type marker data
save('C:\Ahmed\Work\Data\DMD\cellTypeMarker.mat','cellTypeMarker');

%% test sets of genes for enrichment in the top N genes
N = 200;
load('C:\Ahmed\Work\Data\DMD\cellTypeMarker.mat','cellTypeMarker');
% Adult Data
[num txt] = xlsread('C:\Ahmed\Work\Results\DMD\Table~S2 - DMD_combined_coexpressed_genes.xlsx');
topGenes = txt(2:N+1,1);
bottomGenes = txt(end-N+1:end,1);
allGenes = txt(2:end,1);
for i = 1 : length(cellTypeMarker.data)
    overlap_top(i,1) = sum(ismember(lower(topGenes),lower(cellTypeMarker.data{i})));
    p_top(i,1) = hygepdf(overlap_top(i,1), length(allGenes), length(cellTypeMarker.data{i}), N);
    overlap_bottom(i,1) = sum(ismember(lower(bottomGenes),lower(cellTypeMarker.data{i})));
    p_bottom(i,1) = hygepdf(overlap_bottom(i,1), length(allGenes), length(cellTypeMarker.data{i}), N);
end
% BrainSpan
[num txt] = xlsread('C:\Ahmed\Work\Results\DMD\Table~S6 - DMD_BrainSpan_coexpressed_genes.xlsx');
topGenes = txt(2:N+1,1);
bottomGenes = txt(end-N+1:end,1);
allGenes = txt(2:end,1);
for i = 1 : length(cellTypeMarker.data)
    overlap_top(i,1) = sum(ismember(lower(topGenes),lower(cellTypeMarker.data{i})));
    p_top(i,1) = hygepdf(overlap_top(i,1), length(allGenes), length(cellTypeMarker.data{i}), N);
    overlap_bottom(i,1) = sum(ismember(lower(bottomGenes),lower(cellTypeMarker.data{i})));
    p_bottom(i,1) = hygepdf(overlap_bottom(i,1), length(allGenes), length(cellTypeMarker.data{i}), N);
end
% ExonArray
[num txt] = xlsread('C:\Ahmed\Work\Results\DMD\Table~S10 - DMD_exonArray_coexpressed_genes.xlsx');
topGenes = txt(2:N+1,1);
bottomGenes = txt(end-N+1:end,1);
allGenes = txt(2:end,1);
for i = 1 : length(cellTypeMarker.data)
    overlap_top(i,1) = sum(ismember(lower(topGenes),lower(cellTypeMarker.data{i})));
    p_top(i,1) = hygepdf(overlap_top(i,1), length(allGenes), length(cellTypeMarker.data{i}), N);
    overlap_bottom(i,1) = sum(ismember(lower(bottomGenes),lower(cellTypeMarker.data{i})));
    p_bottom(i,1) = hygepdf(overlap_bottom(i,1), length(allGenes), length(cellTypeMarker.data{i}), N);
end

%% save expression data for next step
atlasDir = '/tudelft.net/staff-bulk/ewi/insy/DBL/sjoerdhuisman/ABA_human_brain_probegene/';
dmdDir = '/tudelft.net/staff-bulk/ewi/insy/DBL/amahfouz/MATLAB/Results/DMD/';
donors = {'10021','12876','14380','15496','15697','9861'};
for i = 1 : length(donors)
    % gene expression data
    if i == 1
        MicroarrayExpression = csvread([atlasDir 'gene_expr_normalized_microarray_donor' donors{i} '_2014-11-11.csv']);
        MicroarrayExpression = zscore(MicroarrayExpression,[],2);
    else
        tempExpMat = csvread([atlasDir 'gene_expr_normalized_microarray_donor' donors{i} '_2014-11-11.csv']);
        tempExpMat = zscore(tempExpMat,[],2);
        MicroarrayExpression = [MicroarrayExpression, tempExpMat];
        clear tempExpMat;
    end
end
save([dmdDir 'combinedExpression.mat'],'MicroarrayExpression')

%% test sets of genes for for correlation with the top N genes
atlasDir = '/tudelft.net/staff-bulk/ewi/insy/DBL/sjoerdhuisman/ABA_human_brain_probegene/';
dmdDir = '/tudelft.net/staff-bulk/ewi/insy/DBL/amahfouz/MATLAB/Results/DMD/';
N = 200;
nPerm = 1000;
load([dmdDir 'cellTypeMarker.mat']);
% Adult Data
T = readtable([dmdDir 'Table~S2 - DMD_combined_coexpressed_genes.xlsx']);
rankedGenes = T.Gene_symbol(1:N);
clear T;
% read probe info
T = readtable([atlasDir 'probe_info_2014-11-11.csv']);
probe.gene_symbol = T.gene_symbol;
probe.entrez_id = T.entrez_id;
geneSet_idx = find(ismember(lower(probe.gene_symbol),lower(rankedGenes)));
clear T;
% load expression data
load([dmdDir 'combinedExpression.mat']);
for i = 1 : length(cellTypeMarker.data)
    currMarkerSet_idx = find(ismember(lower(probe.gene_symbol),lower(cellTypeMarker.data{i})));
    corrVals = corr(mean(MicroarrayExpression(currMarkerSet_idx,:))', MicroarrayExpression(geneSet_idx,:)');
    % background corr
    bckCorr = corr(mean(MicroarrayExpression(currMarkerSet_idx,:))', MicroarrayExpression');
%     % determine significance using permutations
%     for R = 1 : nPerm
%         randGenes = randi(size(MicroarrayExpression,1),numel(geneSet_idx),1);
%         randCorr(R,:) = corr(mean(MicroarrayExpression(currMarkerSet_idx,:))', MicroarrayExpression(randGenes,:)');
%     end
    save([dmdDir 'corrData_' cellTypeMarker.name{i} '.mat'], 'corrVals', 'bckCorr')
end




