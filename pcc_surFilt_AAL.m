% Statistical filtering of functional connectivity graph using Iterated Amplitude Adjusted Fourier Transform surrogate time series
% Only connections that reach significance are retained in the final network i.e. FDR-corrected one sided p-value for that connection is < q

% Reference: 
% Simos, N.J., Dimitriadis, S.I., Kavroulakis, E., Manikis, G.C., Bertsias, G., Simos, P., Maris, T.G., & Papadaki, E. (2020). 
% Quantitative identification of functional connectivity disturbances in neuropsychiatric lupus based on resting-state fMRI: A robust machine learning approach. 
% Brain Sciences, 10(11), 1–18. https://doi.org/10.3390/brainsci10110777

% N.J. Simos 2020

clear all
close all
clc

numOfRois = 90;  % Number of ROIs
numOfTp = 147;   % Number of timepoints
numOfsur = 1000; % Number of surrogate timeseries to compute
q = 0.02;        % False discovery rate used as significance threshold in FDR correction for multiple comparisons
% numOfSubs = number of subjects

for subNum = 1:numOfSubs
    
    rng('default')
    
    % load mean ROI timeseries for each subject into 'roi_tss', of dimensions numOfRois x numOfTp
    % AAL atlas was used in our case but this can be of course adjusted for any number of regions (numOfRois)
    
    fcg = squareform(pdist(roi_tss,@pearson_cc)); % calculate functional connectivity graph as pearson correlations between all pairs of ROI timeseries
    pos_fcg = pos_fcg_func(fcg); % Positive connections only, set negative to 0
    abs_fcg = abs(fcg); % Alternatively get absolute of fcg
    
    % numOfsur timeseries for each original ROI timeseries
    surr_tss = zeros(numOfRois, numOfsur, numOfTp);
    for roiNum = 1:numOfRois
        surr_tss(roiNum, :, :) = IAAFT_sur(roi_tss(roiNum, :), numOfsur)';
    end
    
    % numOfsur fcgs, one graph (dimensions: numOfRois x numOfRois) for each set of surrogate timeseries
    surr_fcgs = zeros(numOfsur, numOfRois, numOfRois);
    for surCount = 1:numOfsur
        surr_fcgs(surCount,:,:) = pos_fcg_func(squareform(pdist(squeeze(surr_tss(:,surCount,:)), @pearson_cc)));
    end
    
    % one-sided p-values for each connection (each graph edge)
    % calculated as the number of connection among surrogate roi timeseries greater than the observed connection / number of surrogates
    pval_fcg = zeros(numOfRois, numOfRois);
    for k = 1:numOfRois
        for l = (k+1):numOfRois
            pval_fcg(k,l) = length( find(pos_fcg(k,l) < surr_fcgs(:,k,l)) ) / numOfsur;
        end
    end
    
    [h, ~] = fdr(pval_fcg, q);
    pos_fcg_surr_filt = h.*pos_fcg; % final matrix containing only significant connections
    pos_fcg_surr_filt = triu(pos_fcg_surr_filt) + triu(pos_fcg_surr_filt,1)'; % Conversion to full symmetric matrix, previous version onyly had values in upper triangular portion
    
    globEff(subNum) = efficiency_wei(pos_fcg_surr_filt); % global efficiency
    locEff(subNum, :) = efficiency_wei(pos_fcg_surr_filt, 2)'; % local efficiency
    degree(subNum, :) = degrees_und(pos_fcg_surr_filt); % nodal degree
    BetCentr(subNum, :) = ( betweenness_wei(1./pos_fcg_surr_filt) / ((numOfRois-1)*(numOfRois-2)) )'; % scaling of betweenness centrality as proposed by BCT
    EigCentr(subNum, :) = eigenvector_centrality_und(pos_fcg_surr_filt)'; % eigenvector centrality
    
    % All graph metrics combined (concatenated) in a single feature vector for each subject
    % Exactly as described in the paper, used in our case as input to a Machine Learning classifier
    all_metrics_surr_filt_fcg(subNum, :) = cat(2, globEff(subNum, :), locEff(subNum, :), degree(subNum, :), BetCentr(subNum, :), EigCentr(subNum, :)); 
    
end



