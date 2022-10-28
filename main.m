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

% --------------------------------------------
% fobj = @YourCostFunction
% dim = number of your variables
% Max_iteration = maximum number of iterations
% SearchAgents_no = number of search agents
% lb=[lb1,lb2,...,lbn] where lbn is the lower bound of variable n
% ub=[ub1,ub2,...,ubn] where ubn is the upper bound of variable n
% ---------------------------------------------------------

clear all
clc
format long

SearchAgents_no=30; % Number of search agents/搜索代理的数量

runs=30;            % Number of independent runs/独立运行次数

Function_name='F1'; % Name of the test function that can be from F1 to F23(Table 1 in the paper)
                    % 测试函数的名称，可以从F1到F23（本文的表1） 
   
Max_iteration=250;  % Maximum number of iterations/最大迭代次数

Best=zeros(1,runs);

for i=1:23
    
   Function_name=['F',num2str(i)] 
   
   % Load details of the selected benchmark function/调入所选基准函数的详细信息 
   [lb,ub,dim,fobj]=Get_Functions_details(Function_name);
   
    for j=1:runs

       [Best_score,Best_pos,Convergence_curve]=LsOA(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
        Best(i,j)=Best_score;
       
    end
    
end

% Draw search space/绘制搜索空间 
figure('Position',[500 400 700 290])
subplot(1,2,1);
func_plot(Function_name);
title('Function Topology')
xlabel('x_1');
ylabel('x_2');
zlabel([Function_name,'( x_1 , x_2 )'])

% Convergence curve/收敛曲线
subplot(1,2,2);
semilogy(Convergence_curve,'Color','r')
title('Objective space')
xlabel('Iteration');
ylabel('Best score obtained so far');


display(['The best solution obtained by LsOA is : ', num2str(Best_pos)]);
display(['The best optimal value of the objective function found by LsOA is : ', num2str(Best_score)]);

