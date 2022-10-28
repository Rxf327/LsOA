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

% This function initialize the first population of search agents/初始化
function Positions=initialization(SearchAgents_no,dim,ub,lb)

Boundary_no= size(ub,2); % numnber of boundaries

% If the boundaries of all variables are equal and user enter a signle
% number for both ub and lb 如果所有变量上下界都相同
if Boundary_no==1
    Positions=rand(SearchAgents_no,dim).*(ub-lb)+lb;
end

% If each variable has a different lb and ub/如果每个变量都有不同的lb和ub 
if Boundary_no>1
    for i=1:dim
        ub_i=ub(i);
        lb_i=lb(i);
        Positions(:,i)=rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;
    end
end