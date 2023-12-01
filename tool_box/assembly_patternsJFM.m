function [Members , AssemblyTemplates] = assembly_patterns(SpikeCount,opts)

% Patterns = assembly_patterns(Activitymatrix,opts): extracts assembly patterns from the spike matrix.
% 
% Description of inputs: 	
%   Activitymatrix: spike matrix. Rows represent neurons, columns represent time bins. Thus, each element of the matrix carries the spike count of a given neuron at a given bin.
%       opts: set parameters. All fields described below.
%       opts.threshold.method: defines the method to compute the threshold for assembly detection. Options are:
%           'MarcenkoPastur': uses the analytical bound. 
%           'binshuffling': estimate eigenvalue distribution for independent activity from surrogate matrices generated by shuffling time bins.
%           'circularshift': estimate eigenvalue distribution for independent activity from surrogate matrices generated by random circular shifts of original spike matrix.
%       opts.threshold.permutations_percentile: defines which percentile of the surrogate distribution of maximal eigenvalues is used as statistical threshold. It must be a number between 0 and 100 (95 or larger recommended). Not used when 'MarcenkoPastur' is chosen.
%       opts.threshold.number_of_permutations: defines how many surrogate matrices are generated (100 or more recommended). Not used when 'MarcenkoPastur' is chosen.
%       opts.Patterns.method: defines which method is used to extract assembly patterns. Options are: 'PCA' or 'ICA' (recommended). 
%       opts.Patterns.number_of_iterations: number of iterations for fastICA algorithm (100 or more recommended). Not used when 'PCA' is chosen.
%       opts.Members.method: To define the way to detect the members of the assemblie.
%                           If 'Otsu', it will aply a linear Otsu threshold.
%                           If 'Mean', it will take mean+std as a threshold. (less restrictive)
%                           If 'Sqrt', it will take 1/(sqrt(N)
%       
%
% Description of outputs:
%   Patterns: assembly patterns. Columns denote assembly # and rows neuron #.
% 
% By simply running Patterns = assembly_patterns(Activitymatrix); The following default options will be used: 
%   opts.Patterns.method: 'ICA'
%   opts.threshold.method: 'MarcenkoPastur'
%   opts.Patterns.number_of_iterations: 500.
% 
% This framework is described in: 
% Lopes-dos-Santos V, Ribeiro S, Tort ABL 
% (2013) Detecting cell assemblies in large neuronal populations, 
% Journal of Neuroscience Methods.
%
% Please send bug reports to vitor@neuro.ufrn.br (V�tor)
%
% Adapted by Morici Juan Facundo 07/2023
% - Otsu, Mean and 1/Sqrt(N) thresholding for assemblie's members detection
% - The sign of the AssemblyTemplates was defined using the max(abs(template))
%   Example: If the max(abs(template)) was a negative value, the templated was
%   multiplied by -1

if nargin<2
    opts.threshold.method = 'MarcenkoPastur';
    opts.Patterns.method = 'ICA';
    opts.Patterns.number_of_iterations = 500;
    opts.Members.method = 'Mean'
end
zSpikeCount = zscore(SpikeCount');

CorrMatrix = corr(zSpikeCount);
CorrMatrix(isnan(CorrMatrix))=0;
[eigenvectors,d] = eig(CorrMatrix);
eigenvalues=diag(d);

q = size(zSpikeCount,1)/size(zSpikeCount,2);

if q<1
    fprintf('Error: Number of time bins must be larger than number of neurons \n')
    return
end

%% Finding number of assemblies

switch opts.threshold.method
    case 'MarcenkoPastur'
        fprintf('Using Marcenko-Pastur distribution for estimating number of assemblies \n')
        lambda_max = ((1+sqrt(1/q))^2);
    case 'binshuffling'
        fprintf('Generating control spike count matrix for estimating number of assemblies \n')
        if ~isfield(opts.threshold,'number_of_permutations')
            auxmsdg = 'Please enter a number of surrogates larger than zero: ';
            opts.threshold.number_of_permutations = input(auxmsdg);
        end
        if ~isfield(opts.threshold,'permutations_percentile')
            auxmsdg = 'Please enter percentile for statistical threshold: ';
            opts.threshold.permutations_percentile = input(auxmsdg);
        end
        while opts.threshold.number_of_permutations<=0
            auxmsdg = 'Please enter a number of surrogates larger than zero: ';
            opts.threshold.number_of_permutations = input(auxmsdg);
        end
        fprintf(['Number of permutations:  ' num2str(opts.threshold.number_of_permutations) '\n'])
        control_max_eig = bin_shuffling(SpikeCount,opts.threshold.number_of_permutations);
        lambda_max = prctile(control_max_eig,opts.threshold.permutations_percentile);
    case 'circularshift'
        fprintf('Generating control spike count matrix for estimating number of assemblies \n')
        if ~isfield(opts.threshold,'number_of_permutations')
            auxmsdg = 'Please enter a number of surrogates larger than zero: ';
            opts.threshold.number_of_permutations = input(auxmsdg);
        end
        while opts.threshold.number_of_permutations<=0
            auxmsdg = 'Please enter a number of surrogates larger than zero: ';
            opts.threshold.number_of_permutations = input(auxmsdg);
        end
        if ~isfield(opts.threshold,'permutations_percentile')
            auxmsdg = 'Please enter percentile for statistical threshold: ';
            opts.threshold.permutations_percentile = input(auxmsdg);
        end
        fprintf(['Number of permutations:  ' num2str(opts.threshold.number_of_permutations) '\n'])
        control_max_eig = circular_shift(SpikeCount,opts.threshold.number_of_permutations);
        lambda_max = prctile(control_max_eig,opts.threshold.permutations_percentile);
    case 'thetacycleshift'
        control_max_eig = shift_thetacycle(SpikeCount,opts.threshold.number_of_permutations,opts.binspercycle);
        lambda_max = prctile(control_max_eig,opts.threshold.permutations_percentile);
end
NumberOfAssemblies = sum(eigenvalues>lambda_max);
fprintf(['Number of assemblies detected: ' num2str(NumberOfAssemblies) '\n'])
if NumberOfAssemblies<1
    AssemblyTemplates=[];
    return
end
%% Finding co-activation patterns

switch opts.Patterns.method
    case 'PCA'
        [garbage,PC_position] = sort(-eigenvalues);
        AssemblyTemplates = eigenvectors(:,PC_position(1:NumberOfAssemblies));

    case 'ICA'
        AssemblyTemplates=...
            fast_ica(zSpikeCount,NumberOfAssemblies,opts.Patterns.number_of_iterations);
        
    case 'ICA2'
        AssemblyTemplates=...
            fast_ica(zSpikeCount,length(eigenvalues),opts.Patterns.number_of_iterations);
        prjs = AssemblyTemplates'*zSpikeCount';
        sigassemblies = var(prjs,[],2)>lambda_max;
        AssemblyTemplates = AssemblyTemplates(:,sigassemblies);
end


for i = 1 : NumberOfAssemblies %correction of the sign
   [iii , ii] = max(abs(AssemblyTemplates(:,i)));
   AssemblyTemplates(:,i) = AssemblyTemplates(:,i) * sign(AssemblyTemplates(ii,i));
   clear ii iii
end


switch opts.Members.method
    case 'Otsu'
        Members = [];
        for i = 1 : size(AssemblyTemplates,2)
            Members = [Members , abs(AssemblyTemplates(:,i)) > otsuthresh(abs(AssemblyTemplates(:,i)))];
        end
    case 'Mean'
        Members = [];
        for i = 1 : size(AssemblyTemplates,2)
            Members = [Members , abs(AssemblyTemplates(:,i)) > mean(abs(AssemblyTemplates(:,i))) + std(abs(AssemblyTemplates(:,i)))];
        end
        if isempty(Members)
            Members = [];
        end
    case 'Sqrt'
        Members =  abs(AssemblyTemplates) > 1/(sqrt(size(AssemblyTemplates,1)));
end
end

