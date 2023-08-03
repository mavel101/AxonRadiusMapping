function sphericalMean = fitSphericalMean_mle(dirs, y, sigma, order)

% Estimation of spherical harmonic coefficients of a given order from
% single-shell Rician distributed data using one of the following two estimators:
%  1. 'mle': maximum likelihood estimator
%  2. 'cls': conditional least squares estimator (see https://onlinelibrary.wiley.com/doi/10.1002/mrm.24529)
%
% The Gaussian noise level must estimated a priori. 
%
% input arguments:
% dirs:  Gradient directions for a single b-shell [Nqx3]
% y:     Vectorized diffusion-weighted signals    [Nq x Nvoxels]
% sigma: Vectorized Gaussian noise level                     [1 x Nvoxels]
% order: Maximal spherical harmonics order
% estimator: 'mle' or 'cls'
%
% Author: Jelle Veraart (jelle.veraart@nyulangone.org)
%
% copyright NYU School of Medicine, 2023

X = getSH(order, dirs);
start = X\y;
coef = start;
y(y<eps) = eps;


options = optimset('fminunc'); options = optimset(options,'GradObj','off','Hessian','off','Display','off', 'MaxFunEvals', 20000, 'MaxIter', 2000);
parfor i = 1:size(y, 2)  
    coef(:,i) = fminunc(@(x)LogLikelihood(x,double(X),double(y(:, i)),double(sigma(:,i))),double(start(:,i)),options);
end
   
sphericalMean = coef(1, :) ./ sqrt(pi*4);

end


function f = LogLikelihood(coef, X, s, sigma)

    s_hat = X*coef(:);

    f = logRicianpdf(s, s_hat, sigma);
    f = -sum(f,1)'; % 1 x 1

   
end
        
function [P] = logRicianpdf(x, A, sigma)

    L = 1;
    Arg = (A .* x) ./ (sigma.^2);
    besL1 = besseli(L-1, Arg, 1);
    P = log(x.^L) - log(sigma.^2) + log(A.^(L-1)) - (A.^2 + x.^2)./(2*sigma.^2) + log(besL1) + abs(real(Arg));
   
end

function Y = getSH (Lmax, dirs, CS_phase)

            % Ylm_n = get_even_SH(dirs,Lmax,CS_phase)
            %
            % if CS_phase=1, then the definition uses the Condon-Shortley phase factor
            % of (-1)^m. Default is CS_phase=0 (so this factor is ommited)
            %
            % By: Santiago Coelho (https://github.com/NYU-DiffusionMRI/SMI/blob/master/SMI.m)
            
            if size(dirs,2)~=3
                dirs=dirs';
            end
            Nmeas=size(dirs,1);
            [PHI,THETA,~]=cart2sph(dirs(:,1),dirs(:,2),dirs(:,3)); THETA=pi/2-THETA;
            l=0:2:Lmax;
            l_all=[];
            m_all=[];
            for ii=1:length(l)
                l_all=[l_all, l(ii)*ones(1,2*l(ii)+1)];
                m_all=[m_all -l(ii):l(ii)];
            end
            K_lm=sqrt((2*l_all+1)./(4*pi) .* factorial(l_all-abs(m_all))./factorial(l_all+abs(m_all)));
            if nargin==2 || isempty(CS_phase) || ~exist('CS_phase','var') || ~CS_phase
                extra_factor=ones(size(K_lm));
                extra_factor(m_all~=0)=sqrt(2);
            else
                extra_factor=ones(size(K_lm));
                extra_factor(m_all~=0)=sqrt(2);
                extra_factor=extra_factor.*(-1).^(m_all);
            end
            P_l_in_cos_theta=zeros(length(l_all),Nmeas);
            phi_term=zeros(length(l_all),Nmeas);
            id_which_pl=zeros(1,length(l_all));
            for ii=1:length(l_all)
                all_Pls=legendre(l_all(ii),cos(THETA));
                P_l_in_cos_theta(ii,:)=all_Pls(abs(m_all(ii))+1,:);
                id_which_pl(ii)=abs(m_all(ii))+1;
                if m_all(ii)>0
                    phi_term(ii,:)=cos(m_all(ii)*PHI);
                elseif m_all(ii)==0
                    phi_term(ii,:)=1;
                elseif m_all(ii)<0
                    phi_term(ii,:)=sin(-m_all(ii)*PHI);
                end
            end
            Y_lm=repmat(extra_factor',1,Nmeas).*repmat(K_lm',1,Nmeas).*phi_term.*P_l_in_cos_theta;
            Y=Y_lm';
end