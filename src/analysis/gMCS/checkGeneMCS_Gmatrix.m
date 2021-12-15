function [IsCutSet, IsMinimal, geneNotMinimal] = checkGeneMCS_Gmatrix(model, G_matrix, gMCSs)
% Check if the gMCS is Minimal and Cut Set
% gMCS = gene Minimal Cut Set
% Minimal = the knock-off of all the genes of the gMCS is neccesary to
% eliminate the objective function
% Cut Set = The elimination of all the genes ensures the elimination of the
% objective reaction
%
% USAGE:
%
%    [IsCutSet, IsMinimal, geneNotMinimal] = checkGeneMCS(model, gMCSs, isoform_separator)
%
% INPUTS:
%    model:                 Metabolic model in COBRA format.
%    G_matrix:              Character indicating the MAT object which
%                           contains the G matrix information or an struct
%                           object with all the information of the
%                           G matrix.
%    gMCSs:                 cell array which stores all the gMCSs previously
%                           calculated. Each position contains a gMCS.
%
% OUTPUT:
%    IsCutSet:              array which indicates the biomass producion if
%                           all genes of the gMCSs are deleted
%    IsMinimal:             array which indicates the maximum of biomass if
%                           any gene is deleted
%    geneNotMinimal:        cell array which stores all the genes which
%                           make the gMCS not minimal (this errors are
%                           produced by the restriction of coputation
%                           time).
%
% NOTE:
%    The isoform separator is neccesary to consider some genes as the same
%    gene. Example: genes 555.1, 555.2 and 555.3 are considered as gene 555
%
% .. Author: - Luis V. Valcarcel 2017-11-13

if (nargin < 2)
    error('No G matrix introduced');
end

if nargin < 3 || isempty(gMCSs)
    % Nothing to check
    BiomassCutSet = nan;
    BiomassMinimal = nan;
    geneNotMinimal = nan;
    warning('There is no gMCS to check');
    return
end

% Load the information from the G matrix
if ischar(G_matrix)
    G_file = [G_matrix '.mat'];
    if exist(G_file, 'file')
        load(G_file, 'G', 'G_ind', 'G_rxns', 'separate_isoform');
    else
        error(['No file with this name: ' G_file])
    end
    G_matrix = struct();
    G_matrix.G = G;
    G_matrix.G_ind = G_ind;
    if exist('G_rxns', 'var')
        G_matrix.G_rxns = G_rxns;
    end
    if exist('separate_isoform', 'var')
        G_matrix.separate_isoform = separate_isoform;
    end
elseif isstruct(G_matrix)
    G = G_matrix.G;
    G_ind = G_matrix.G_ind;
    if isfield(G_matrix, 'G_rxns')
        G_rxns = G_matrix.G_rxns;
    end
    if isfield(G_matrix, 'separate_isoform')
        separate_isoform = G_matrix.separate_isoform;
    end
else
    error('No format recognised for G matrix (char or struct)');
end

% check the fields to ensure is you can use it in the function
assert(size(model.S,2) == size(G,2), 'G matrix and model differ in number of reactions.');
if exist('G_rxns', 'var')
    assert(isequal(model.rxns, G_rxns), 'G matrix is not the same, different reaction order in the G matrix and the model.');
end
assert(all(ismember([G_ind{:}], strtok(model.genes,separate_isoform))), 'G matrix is not the same, different reaction order in the G matrix and the model.');


% Expand information of the G_matrix for easy search of G_ind
Gind2genes_genes = unique([G_ind{:}]);
Gind2genes_mat = spalloc(length(G_ind),length(Gind2genes_genes), sum(cellfun(@length,G_ind)));
for i = 1:length(G_ind)
    Gind2genes_mat(i,ismember(Gind2genes_genes, G_ind{i})) = 1;
end
% make the matrix to sum 1 by row
Gind2genes_mat = Gind2genes_mat./repmat(sum(Gind2genes_mat,2), 1,size(Gind2genes_mat,2));
assert(all(abs(sum(Gind2genes_mat,2)-1)<1e-4));
% store the results
G_matrix.Gind2genes_genes = Gind2genes_genes;
G_matrix.Gind2genes_mat = Gind2genes_mat;


% Define the variables
BiomassCutSet = zeros(size(gMCSs));
BiomassMinimal = zeros(size(gMCSs));
geneNotMinimal = cell(size(gMCSs));

% Loop for all the gMCSs
showprogress(0,'Check all gMCS');
for i = 1:numel(gMCSs)
    showprogress(i/numel(gMCSs));
    gmcs = gMCSs{i};
    % Check if Cut Set
    BiomassCutSet(i) = checkGeneCutSet_Gmat(model, G_matrix, gmcs);
    % Check if Minimal
    if nargout > 1
        %     if true
        if length(gmcs) == 1
            sol = optimizeCbModel(model);
            BiomassMinimal(i) = sol.f;
        else
            check_Minimal = zeros(size(gmcs));
            % Loop to check if a subset of the gMCS is also gMCS
            for g = 1:numel(gmcs)
                idx = true(size(gmcs));
                idx(g) = 0; % not include the gene in pos "g"
                gmcs_aux = gmcs(idx);
                check_Minimal(g) = checkGeneCutSet_Gmat(model, G_matrix, gmcs_aux);
            end
            % Select a threshold for the objective function similar to the MCS target
            threshold_objective = 1e-3;
            geneNotMinimal{i} = gmcs(check_Minimal>threshold_objective);
            BiomassMinimal(i) = min(check_Minimal); % if any of them is a Cut Set,
            % the sum is > 0, not minimal
        end
    end
    
end


th = getCobraSolverParams('LP', 'feasTol')*100; % same as default target b for gMCS

IsCutSet = BiomassCutSet < th;
IsMinimal = BiomassMinimal > th;

end


function biomassCS = checkGeneCutSet_Gmat(model_raw, G_matrix, gMCS)
% Check if the gMCS is a Cut Set
% Cut Set = The elimination of all the genes ensures the elimination of the
% objective reaction
%
% USAGE:
%
%    checkCS = checkGeneCutSet(model_raw,gMCS,isoform_separator)
%
% INPUTS:
%    model:                 Metabolic model in COBRA format.
%    G_matrix:              Character indicating the MAT object which
%                           contains the G matrix information or an struct
%                           object with all the information of the
%                           G matrix.
%
% OUTPUT:
%    biomassCS:             biomass value which indicates if the set of
%                           genes are a Cut Sets
%
% .. Author: - Luis V. Valcarcel 2021-11-29


% Generate a model in which we delete the genes of the Cut Set
[model, ~, ~, deletedGenes] = fromGeneMCStoRxnKO(model_raw, G_matrix,  gMCS);

% Check that eliminated genes and geneList are the same
if ~isequal(gMCS,deletedGenes)
    error('Elimated genes and Cut Set are not the same');
end

% check if the elimination of all genes is effective
sol = optimizeCbModel(model);
biomassCS = sol.f ;

end


function [model, hasEffect, constrRxnNames, deletedGenes] = fromGeneMCStoRxnKO(model, G_matrix, gmcs)
% Deletes one or more genes and constrain the reactions
% affected to zero and appends '_deleted' to the gene(s)
%
% USAGE:
%
%    [model, hasEffect, constrRxnNames, deletedGenes] = fromGeneMCStoRxnKO(model, G_matrix, gmcs, downRegFraction)
%
% INPUT:
%    model:             COBRA model with the appropriate constrains for a
%                       particular condition
%
% OPTIONAL INPUTS:
%    G_matrix:          Character indicating the MAT object which contains
%                       the G matrix information or an struct object with
%                       all the information of the G matrix.
%    gmcs:              List of genes to be deleted, one simple gMCS (Default =  NULL)
%
% OUTPUTS:
%    model:             COBRA model with the selected genes deleted
%    hasEffect:         True if the gene deletion has an effect on the model
%    constrRxnNames:    Reactions that are associated to the genes in `deletedGenes`
%    deletedGenes:      The list of genes removed from the model.
%
% .. Authors:
%       - Luis V. Valcárcel 2021/11/24 (function based on deleteModelGenes)
%

if (nargin < 2)
    error('No G matrix introduced');
else
    G = G_matrix.G;
    G_ind = G_matrix.G_ind;
%     if isfield(G_matrix, 'G_rxns')
%         G_rxns = G_matrix.G_rxns;
%     end
    if isfield(G_matrix, 'separate_isoform')
        separate_isoform = G_matrix.separate_isoform;
    end
    Gind2genes_mat = G_matrix.Gind2genes_mat;
    Gind2genes_genes = G_matrix.Gind2genes_genes;
    
end
if (nargin < 3)
    gmcs = {''};
end

if (~iscell(gmcs))
    geneName = gmcs;
    clear geneList;
    gmcs{1} = geneName;
end

if (~isfield(model,'genes'))
    error('Gene-reaction associations not included with the model');
end


%deletedGenes is a cell array for returning the genes that are
%eliminated from the model.
deletedGenes = {};

% Find gene indices in G_mat
[isInModel,geneInd] = ismember(gmcs,strtok(model.genes, separate_isoform));
deletedGenes = gmcs( find(geneInd) );

% mark rows of G matrix that are fully included in the gMCS
rowsInd = sum(Gind2genes_mat(:,ismember(Gind2genes_genes, deletedGenes)),2)==1;
num_rowsInd = sum(rowsInd);
% Find rxns associated with this gene

if num_rowsInd == 0
    constrRxnNames = {''};
    hasEffect = false;
elseif num_rowsInd == 1
    rxnInd = G(rowsInd,:)>0;
    hasEffect = true;
else
    rxnInd = sum(G(rowsInd,:)>0)>0;
    hasEffect = true;
end
if hasEffect
    % Full deletion
    constrRxnNames = model.rxns(rxnInd);
    %     model = changeRxnBounds(model,constrRxnNames,0,'b');
    model.lb(rxnInd) = 0;
    model.ub(rxnInd) = 0;
end

if any(~isInModel)
    warning(['Gene',' ',gmcs{~isInModel}, ' not in model!']);
end

end

