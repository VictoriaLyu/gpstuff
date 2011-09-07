function [ps, fm, fprcs, m, prcs]=gp_avpredcomp(gp, x, y, varargin)
%GP_AVPREDCOMP  Average predictive comparison for Gaussian process model
%
%  Description
%    [PS, FM, FPRCS, M, PRCS]=GP_AVPREDCOMP(GP, X, Y, OPTIONS)
%    Takes a Gaussian process structure GP together with a matrix X
%    of training inputs and vector Y of training targets, and
%    returns average predictive comparison estimates for each
%    input. PS is the probability of knowing the sign of the
%    average change in the latent outcome for each input variable. 
%    FM and FPRCS are estimated predictive relevances for means and
%    percentiles for each input variable when latent outcome
%    variable is computed through the inverse link function. M and
%    PRCS are estimated predictive relevances for mean and
%    percentiles for each input variable when outcome variable is
%    the latent variable.
%
%    OPTIONS is optional parameter-value pair
%      z         - optional observed quantity in triplet (x_i,y_i,z_i)
%                  Some likelihoods may use this. For example, in
%                  case of Poisson likelihood we have z_i=E_i, that
%                  is, expected value for ith case.
%      nsamp     - determines the number of samples used (default=500).
%      prctiles  - determines percentiles that are computed from 
%                  0 to 100 (default=[5 95]).
%
%  See also
%    GP_PRED

% Copyright (c) 2011      Jaakko Riihimäki
% Copyright (c) 2011      Aki Vehtari

% This software is distributed under the GNU General Public
% License (version 3 or later); please refer to the file
% License.txt, included with the software, for details.

ip=inputParser;
ip.FunctionName = 'GP_AVPREDCOMP';
ip.addRequired('gp',@isstruct);
ip.addRequired('x', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
ip.addRequired('y', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
ip.addParamValue('z', [], @(x) isreal(x) && all(isfinite(x(:))))
ip.addParamValue('nsamp', 500, @(x) isreal(x) && isscalar(x))
ip.addParamValue('prctiles', [5 95], @(x) isreal(x) && isvector(x) && all(x>=0) && all(x<=100))

ip.parse(gp, x, y, varargin{:});
z=ip.Results.z;
nsamp=ip.Results.nsamp;
prctiles=ip.Results.prctiles;

[n, nin]=size(x);

ps=zeros(nin,1);
fm=zeros(nin,1);
fprcs=zeros(nin,length(prctiles));
m=zeros(nin,1);
prcs=zeros(nin,length(prctiles));

covx=cov(x);

stream = RandStream('mrg32k3a');
prevstream=RandStream.setDefaultStream(stream);

% loop through the input variables
for k1=1:nin
  fprintf('k1=%d\n',k1)
    %- Compute the weight matrix based on Mahalanobis distances:
    x_=x; x_(:,k1)=[];
    covx_=covx; covx_(:,k1)=[]; covx_(k1,:)=[];
    % weight matrix:
    W=zeros(n);
    for i1=1:n
        x_diff=bsxfun(@minus,x_(i1,:),x_((i1+1):n,:))';
        W(i1,(i1+1):n)=1./(1+sum(x_diff.*(covx_\x_diff)));
    end
    W=W+W'+eye(n);
    %-
    
    rsubstream=round(rand*10e9);
    
    num=zeros(1,nsamp); den=0;
    numf=zeros(1,nsamp); pp=0;
    for i1=1:n
        % inputs of interest
        ui=x(i1, k1);
        ujs=x(:, k1);
        
        % replicate same values for other inputs
        xrep=repmat(x(i1,:),n,1); xrep(:,k1)=ujs;
        
        Udiff=ujs-ui;
        Usign=sign(Udiff);
        stream.Substream = rsubstream;
        % draw random samples from the posterior
        fs = gp_rnd(gp, x, y, xrep, 'z', z, 'nsamp', nsamp);
        
        % compute latent values through the inverse link function
        if isfield(gp.lik.fh, 'invlink')
            ilfs = gp.lik.fh.invlink(gp.lik, fs, repmat(z,1,nsamp));
            % average change in outcome
            num=num+sum(bsxfun(@times,W(:,i1).*Usign,bsxfun(@minus,ilfs,ilfs(i1,:))));
        end
        
        % average change in latent outcome
        numf=numf+sum(bsxfun(@times,W(:,i1).*Usign,bsxfun(@minus,fs,fs(i1,:))));

        % average change in input
        den=den+sum(W(:,i1).*Udiff.*Usign);
        
    end
    
    % means and percentiles when outcome is the latent function
    numfden=numf./den;
    fm(k1,1)=mean(numfden);
    fprcs(k1,:)=prctile(numfden,prctiles);
    
    % means and percentiles when outcome is computed through the inverse
    % link function
    numden=num./den;
    m(k1,1)=mean(numden);
    prcs(k1,:)=prctile(numden,prctiles);

    % probability of knowing the sign of the change in
    % latent function 
    ps(k1,1)=mean(numfden>0);
    if ps(k1,1)<0.5
      ps(k1,1)=1-ps(k1,1);
    end
    
end

if ~isfield(gp.lik.fh, 'invlink')
    m=fm;
    prcs=fprcs;
end

RandStream.setDefaultStream(prevstream);
