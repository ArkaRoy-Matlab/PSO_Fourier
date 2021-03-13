% ***************************************************************
% *** Matlab function for pca for reducing data dimension
% *** Source Code is mainly written for research purposes. The codes are
% *** having copyrights and required proper citations whenever it is used.
% *** Originated by:
% ***       Mr. Arka Roy (email: arka.phy@gmail.com)
% ***       Dr. Chandra Prakash Dubey (email:p.dubey48@gmail.com)
% ***       Mr. M. Prasad (email:prasadgoud333@gmail.com)
% ***       Crustal Processes Group, National Centre for Earth Science Studies,
% ***       Ministry of Earth Sciences, Government of India
% ***       Thiruvanthapuram, Kerala, India
% ***************************************************************
%% Matlab function for pca for reducing data dimention

% consider an data set of m parameters and n models
%   For example : data=rand(m,n);
%   data dimension should be (mxn) order
 %% Matlab function for Principal Component Analysis
function  [pc,Evalues,W] =pca_reduction(data)

    % remove the mean variable-wise (row-wise)
    data=data-repmat(mean(data,2),1,size(data,2));

% calculate eigenvectors (loadings) W, and eigenvalues of the covariance matrix
    [W, EvalueMatrix] = eig(cov(data'));
    Evalues = diag(EvalueMatrix);

% order by largest eigenvalue
    Evalues = Evalues(end:-1:1);
    W = W(:,end:-1:1); W=W';  

% generate PCA component space (PCA scores)
    pc = W * data;
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%% End of Code %%%%%%%%%%%%%%%%%%%%%%%%%%%
