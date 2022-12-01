% fdr_bh() - Executes the Benjamini & Hochberg (1995) procedure for
%            controlling the false discovery rate (FDR) of a family of 
%            hypothesis tests. FDR is the expected proportion of rejected
%            hypotheses that are mistakenly rejected (i.e., the null
%            hypothesis is actually true for those tests). FDR is a
%            somewhat less conservative/more powerful method for correcting
%            for multiple comparisons than methods like Bonferroni
%            correction that provide strong control of the family-wise
%            error rate (i.e., the probability that one or more null
%            hypotheses are mistakenly rejected).
%
% Usage:
%  >> [h, crit_p]=fdr(pvals,q,method,report);
%
% Required Input:
%   pvals - A vector or matrix (two dimensions or more) containing the
%           p-value of each individual test in a family of tests.
%
% Optional Inputs:
%   q       - The desired false discovery rate. {default: 0.05}
%   method  - ['pdep' or 'dep'] If 'pdep,' the original Bejnamini & Hochberg
%             FDR procedure is used, which is guaranteed to be accurate if
%             the individual tests are independent or positively dependent
%             (e.g., positively correlated).  If 'dep,' the FDR procedure
%             described in Benjamini & Yekutieli (2001) that is guaranteed
%             to be accurate for any test dependency structure 

%
%
% Outputs:
%   h       - A binary vector or matrix of the same size as the input "pvals."
%             If the ith element of h is 1, then the test that produced the 
%             ith p-value in pvals is significant (i.e., the null hypothesis
%             of the test is rejected).
%   crit_p  - All p-values less than or equal to crit_p are significant
%             (i.e., their null hypotheses are rejected).  If no p-values are
%             significant, crit_p=0.
%
%
% References:
%   Benjamini, Y. & Hochberg, Y. (1995) Controlling the false discovery
%     rate: A practical and powerful approach to multiple testing. Journal
%     of the Royal Statistical Society, Series B (Methodological). 57(1),
%     289-300.
%
%   Benjamini, Y. & Yekutieli, D. (2001) The control of the false discovery
%     rate in multiple testing under dependency. The Annals of Statistics.
%     29(4), 1165-1188.

% For complex brain networks see also :
% J. Toppi, F. De Vico Fallani, G. Vecchiato, A. G. Maglione, F. Cincotti, D. Mattia, S. Salinari, F. Babiloni, and L. Astolfi 
% How the Statistical Validation of Functional Connectivity Patterns Can Prevent Erroneous Definition of Small-World Properties 
%of a Brain Connectivity Network,  Volume 2012 (2012), Article ID 130985, 13 pages


%DIMITRIADIS STAVROS 7/2011

function [h crit_p]=fdr(pvals,q,method)

if nargin<1,
    error('You need to provide a vector or matrix of p-values.');
else
    if ~isempty(find(pvals<0,1) & find(pvals>1,1) ),
        error(' p-values are out of range');
    end
end

if nargin<2,
    q=.01;
end

if nargin<3,
    method='pdep';
end


s=size(pvals);

if (length(s)>2) || s(1)>1,
    p_sorted=sort(reshape(pvals,1,prod(s)));
else
    %p-values are already a row vector
    p_sorted=sort(pvals);
end
m=length(p_sorted); %number of tests

if strcmpi(method,'pdep'),
    %BH procedure for independence or positive dependence
    thresh=[1:m]*q/m;
elseif strcmpi(method,'dep')
    %BH procedure for any dependency structure
    denom=m*sum(1./[1:m]);
    thresh=[1:m]*q/denom;
else
    error('Argument ''method'' needs to be ''pdep'' or ''dep''.');
end

rej=p_sorted<=thresh;

max_id=find(rej,1,'last'); 

if isempty(max_id),
    crit_p=0;
    h=pvals*0;
else
    crit_p=p_sorted(max_id);
    h=pvals<=crit_p;
end





