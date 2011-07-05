function [Eft, Varft, lpyt, Eyt, Varyt] = gpla_loopred(gp, x, y, varargin)
%GPLA_LOOPRED  Leave-one-out predictions with Laplace approximation
%
%  Description
%    [EFT, VARFT, LPYT, EYT, VARYT] = GPLA_LOOPRED(GP, X, Y,
%    OPTIONS) takes a Gaussian process structure GP together with a
%    matrix XT of input vectors, matrix X of training inputs and
%    vector Y of training targets, and evaluates the leave-one-out
%    predictive distribution at inputs X. Returns a posterior mean
%    EFT and variance VARFT of latent variables, the posterior
%    predictive mean EYT and variance VARYT of observations, and logarithm
%    of posterior predictive density PYT at input locations X.
%
%    Laplace leave-one-out is approximated by computing a Laplace
%    site-term and computing a Laplace cavity distribution as
%    leave-one-out posterior for the ith latent value. Since the
%    ith likelihood has influenced other site terms through the
%    prior, this estimate can be over-optimistic.
%
%    OPTIONS is optional parameter-value pair
%      z  - optional observed quantity in triplet (x_i,y_i,z_i)
%           Some likelihoods may use this. For example, in case of
%           Poisson likelihood we have z_i=E_i, that is, expected
%           value for ith case.
%
%  See also
%    GPLA_E, GPLA_G, GP_PRED, DEMO_SPATIAL, DEMO_CLASSIFIC
  
% Copyright (c) 2011  Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 2 or later); please refer to the file 
% License.txt, included with the software, for details.

  ip=inputParser;
  ip.FunctionName = 'GPLA_LOOPRED';
  ip.addRequired('gp', @(x) isstruct(x) || iscell(x));
  ip.addRequired('x', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
  ip.addRequired('y', @(x) ~isempty(x) && isreal(x) && all(isfinite(x(:))))
  ip.addParamValue('z', [], @(x) isreal(x) && all(isfinite(x(:))))
  ip.parse(gp, x, y, varargin{:});
  z=ip.Results.z;

  if ~iscell(gp)
    % Single GP
    
    % latent posterior
    [f, sigm2ii] = gpla_pred(gp, x, y, x);
    
    % "site parameters"
    W        = -gp.lik.fh.llg2(gp.lik, y, f, 'latent', z);
    deriv    = gp.lik.fh.llg(gp.lik, y, f, 'latent', z);
    sigm2_t  = 1./W;
    mu_t     = f + sigm2_t.*deriv;
    
    % "cavity parameters"
    sigma2_i = 1./(1./sigm2ii-1./sigm2_t);
    myy_i    = sigma2_i.*(f./sigm2ii-mu_t./sigm2_t);
    
    % leave-one-out predictions
    Eft=myy_i;
    Varft=sigma2_i;
    if nargin==3
      lpyt = gp.lik.fh.predy(gp.lik, myy_i, sigma2_i, y, z);
    elseif nargin>3
      [lpyt,Eyt,Varyt] = gp.lik.fh.predy(gp.lik, myy_i, sigma2_i, y, z);
    end

  else
    % Cell array of GPs
%    nGP = numel(gp);
%    for j = 1:nGP
%      [~,~,~,~,~,~,~,~,muvec_i,sigm2vec_i,Z_i] = gpep_e(gp_pak(gp{j}), ...
%                                                        gp{j}, x, y, 'z', z);
      
%      P_TH(j,:) = gp{j}.ia_weight;
%      Eft_grid(j,:)=muvec_i;
%      Varft_grid(j,:)=sigm2vec_i;
%      pyt_grid(j,:)=Z_i;
%      n=length(y);
%      if nargout > 3
%        for cvi=1:n
%          if isempty(z)
%            [~, Eyt_grid(j,cvi), Varyt_grid(j,cvi)] = ...
%                gp{j}.lik.fh.predy(gp{j}.lik, muvec_i(cvi), ...
%                      sigm2vec_i(cvi), [], []);
%          else
%            [~, Eyt_grid(j,cvi), Varyt_grid(j,cvi)] = ...
%                gp{j}.lik.fh.predy(gp{j}.lik, muvec_i(cvi), ...
%                      sigm2vec_i(cvi), [], z(cvi));
%          end
%        end
%      end
%    end
    
%    ft = zeros(size(Eft_grid,2),501);
%    for j = 1 : size(Eft_grid,2);
%        ft(j,:) = Eft_grid(1,j)-10*sqrt(Varft_grid(1,j)) : 20*sqrt(Varft_grid(1,j))/500 : Eft_grid(1,j)+10*sqrt(Varft_grid(1,j));  
%    end
    
%    % Calculate the density in each grid point by integrating over
%    % different models
%    pft = zeros(size(Eft_grid,2),501);
%    for j = 1 : size(Eft_grid,2)
%        pft(j,:) = sum(norm_pdf(repmat(ft(j,:),size(Eft_grid,1),1), repmat(Eft_grid(:,j),1,size(ft,2)), repmat(sqrt(Varft_grid(:,j)),1,size(ft,2))).*repmat(P_TH,1,size(ft,2)),1); 
%    end

%    % Normalize distributions
%    pft = bsxfun(@rdivide,pft,sum(pft,2));

%    % Widths of each grid point
%    dft = diff(ft,1,2);
%    dft(:,end+1)=dft(:,end);

%    % Calculate mean and variance of the distributions
%    Eft = sum(ft.*pft,2)./sum(pft,2);
%    Varft = sum(pft.*(repmat(Eft,1,size(ft,2))-ft).^2,2)./sum(pft,2);
    
%    if nargout > 3
%      Eyt = sum(Eyt_grid.*repmat(P_TH,1,size(Eyt_grid,2)),1);
%      Varyt = sum(Varyt_grid.*repmat(P_TH,1,size(Eyt_grid,2)),1) + sum((Eyt_grid - repmat(Eyt,nGP,1)).^2, 1);
%      Eyt=Eyt';
%      Varyt=Varyt';
%    end
%    lpyt = log(sum(bsxfun(@times,pyt_grid,P_TH),1)');

%  end
end