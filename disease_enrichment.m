%%% 25 June 201
%%% check for asd/id enrichment
%%% 5 Aug: added enrichment analysis for ADD, OCD, and Dyslexia

%% 
%select gene of interest
geneName = 'DMD';
if ispc
    dataDir = 'C:/Users/amahfouz/SURFdrive/Projects/DMD/Data/';
    resDir = 'C:/Users/amahfouz/SURFdrive/Projects/DMD/Results/';
end
% read the ASD/ID genes (Fereydoun Genome Res 2015)
T = readtable([dataDir 'ID_2_Autism_4_Severe_Missense.Clean_WithNew.txt'],...
    'Delimiter','\t','ReadVariableNames',false);
geneSet{1} = unique(T.Var1);
geneSetName{1} = 'ASD_ID';
clear T;
% read the SFARI ASD genes (downloaded on 25 June 2015)
T = readtable([dataDir 'SFARI_ASD_25JUNE2015.txt']);
geneSet{2} = unique(T.GeneSymbol);
geneSetName{2} = 'SFARI_ASD';
clear T;
% read the ADD, OCD, and Dyslexia genes (downloaded from DisGeNet on 30 July 2015)
% T = readtable([dataDir 'DisGeNET_ADD.xls']);
T = readtable([dataDir 'DisGeNET_ADD.csv']);
geneSet{3} = unique(T.Symbol);
geneSetName{3} = 'ADD';
clear T;
% T = readtable([dataDir 'DisGeNET_OCD.xls']);
T = readtable([dataDir 'DisGeNET_OCD.csv']);
geneSet{4} = unique(T.Symbol);
geneSetName{4} = 'OCD';
clear T;
% T = readtable([dataDir 'DisGeNET_Dyslexia.xls']);
T = readtable([dataDir 'DisGeNET_Dyslexia.csv']);
geneSet{5} = unique(T.Symbol);
geneSetName{5} = 'Dyslexia';
clear T;
% isoform names
isoForms = {'Dp71_Dp40', 'Dp140', 'Dp427'};

%% Adult Human Brain
% read the genes correlated with DMD in the adult brain
load([dataDir 'corrMat_combined.mat']);
[sortedCorrMat IX] = sort(corrMat_combined,2,'descend');
sortedGenes = probe.gene_symbol(IX);
top200_genes = sortedGenes(1:200);
% test the enrichment of all gene sets
for i = 1 : length(geneSet)
    I1 = find(ismember(sortedGenes,geneSet{i})==1);
    p(i) = ranksum(1:length(sortedGenes), I1, 'tail', 'right');
    % save gene ranks to file
    T = table(sortedGenes(I1), I1, 'VariableNames', {'Gene_symbol','rank'});
%     writetable(T, [resDir 'disease_enrichment_adult.xlsx'], 'Sheet', i)
    writetable(T, [resDir 'disease_enrichment_adult_' geneSetName{i} '.csv'])
    clear T;
end
% correct for multiple testing
p_corrected_adult = multtest(p,'method','BH');
% % change sheet names
% xlsheets(geneSetName, [resDir 'disease_enrichment_adult.xlsx']);

%% BrainSpan exons
% load exon data
load([dataDir 'dmd_exon_data.mat']);
load([dataDir 'transcript.mat']);
load([dataDir 'exon_sample.mat']);
% read the top 200 genes correlated with DMD in the exon brainSpan data
load([dataDir 'corrMat_exonArray_individual_exons_log2.mat']);
removedRows = logical(prod(double(isnan(corrMat))));
corrMat(:,removedRows) = [];
transcript.gene_symbol(removedRows) = [];
transcript.entrez_id(removedRows) = [];
transcript.start(removedRows) = [];
transcript.end(removedRows) = [];
pVal(:,removedRows) = [];
for i = 1 : size(corrMat,1)
    % sort the data
    [sortedCorrMat IX] = sort(corrMat(i,:),2,'descend');
    rankedGenes = transcript.gene_symbol(IX);
    [b,m,~] = unique(rankedGenes, 'stable');
    [~,IX2] = sort(m);
    sortedGenes = b(IX2);
    % test the enrichment of all gene sets
    for j = 1 : length(geneSet)
        I1 = find(ismember(sortedGenes,geneSet{j})==1);
        p(j) = ranksum(1:length(sortedGenes), I1, 'tail', 'right');
        % save gene ranks to file
        T = table(sortedGenes(I1), I1, 'VariableNames', {'Gene_symbol','rank'});
%         writetable(T, [resDir 'disease_enrichment_developing_' isoForms{i} '.xlsx'], 'Sheet', j)
        writetable(T, [resDir 'disease_enrichment_developing_' isoForms{i} '_' geneSetName{j} '.csv'])
        clear T;
    end
    % correct for multiple testing
    p_corrected(i,:) = multtest(p,'method','BH');
    % change sheet names
%     xlsheets(geneSetName, [resDir 'disease_enrichment_developing_' isoForms{i} '.xlsx']);
end

%% save p-vales
T = table([p_corrected_adult(1);p_corrected(:,1)], [p_corrected_adult(2);p_corrected(:,2)],...
    [p_corrected_adult(3);p_corrected(:,3)], [p_corrected_adult(4);p_corrected(:,4)],...
    [p_corrected_adult(5);p_corrected(:,5)], 'VariableNames', geneSetName, 'RowNames', [{'Adult'},isoForms]);
writetable(T, [resDir 'disease_enrichment_pValues_corrected.csv'])
clear T;
        
%% make figure
custom_map = csvread('redblue_256_rgb.txt');
custom_map = custom_map / 255;
custom_map = flipud(custom_map);

isoForms_names = {'Dp71\_Dp40', 'Dp140'    'Dp427'};
figure,
imagesc(-log10([p_corrected_adult;p_corrected])), colormap(custom_map), colorbar
% imagesc(fliplr(dmd_data_updated')), colormap(custom_map), colorbar
X = -log10([p_corrected_adult;p_corrected]);
hold on;
for i = 1:size(X,1)
   plot([.5,size(X,2)+.5],[i-.5,i-.5],'k-');
end
for i = 1:size(X,2)
    plot([i-.5,i-.5],[.5,size(X,1)+.5],'k-');
end
hold off
axis image
set(gca, 'XTickLabel', geneSetName, 'XTick', 1:length(geneSetName));
set(gca, 'YTickLabel', [{'Adult'},isoForms], 'YTick', 1:4);



