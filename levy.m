% Input parameters
% n     -> Number of steps 
% m     -> Number of Dimensions 
% beta  -> Power law index  % Note: 1 < beta < 2
% Output 
% z     -> 'n' levy steps in 'm' dimension
%_________________________________________________________________________
%  Lioness Optimization Algorithm source code (Developed in MATLAB R2020a)
%
%
% paper:
%  Dong Li, Xiaofei Ren, Xiujuan Lei.  
%  Lioness Optimization Algorithm：A meta-heuristic algorithm inspired by lioness hunting
%  
%  E-mails: ddli1009@xupt.edu.cn            (Dong Li)
%           171506694@qq.com                (Xiaofei Ren)
%           xjlei@snnu.edu.cn               (Xiujuan Lei) 
%_________________________________________________________________________

function [z] = levy(n,m,beta)

    num = gamma(1+beta)*sin(pi*beta/2); % used for Numerator/分子
    
    den = gamma((1+beta)/2)*beta*2^((beta-1)/2); % used for Denominator/分母

    sigma_u = (num/den)^(1/beta);% Standard deviation/标准偏差

    u = random('Normal',0,sigma_u,n,m); 
    
    v = random('Normal',0,1,n,m);

    z =u./(abs(v).^(1/beta));

  
  end

