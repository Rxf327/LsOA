%_________________________________________________________________________
%  Lioness Optimization Algorithm source code (Developed in MATLAB R2020a)
%
%
% paper:
%  Dong Li, Xiaofei Ren, Xiujuan Lei.  
%  Lioness Optimization Algorithm��A meta-heuristic algorithm inspired by lioness hunting
%  
%  E-mails: ddli1009@xupt.edu.cn            (Dong Li)
%           171506694@qq.com                (Xiaofei Ren)
%           xjlei@snnu.edu.cn               (Xiujuan Lei) 
%_________________________________________________________________________


function [Top_lioness_fit,Top_lioness_pos,Convergence_curve]=LsOA(SearchAgents_no,Max_iter,lb,ub,dim,fobj)

% Initialize position and fitness of the top_lioness/��ʼ��top_lioness��λ�ú���Ӧ��ֵ 
Top_lioness_pos=zeros(1,dim);
Top_lioness_fit=inf;  % change this to -inf for maximization problems/Ҫ���������⣬Ҫ�������Ϊ-inf

Convergence_curve=zeros(1,Max_iter);
step=zeros(SearchAgents_no,dim);
fitness=inf(SearchAgents_no,1);  % change this to -inf for maximization problems/Ҫ���������⣬Ҫ�������Ϊ-inf

% Initialize the positions of search agents/��ʼ�����������λ��
position=initialization(SearchAgents_no,dim,ub,lb);
  
Xmin=repmat(ones(1,dim).*lb,SearchAgents_no,1);
Xmax=repmat(ones(1,dim).*ub,SearchAgents_no,1);
         

Iter=0;  % Loop counter/ѭ��������
FADs=0.2;   
P=0.5;     

%%% Initialization of the "Central Circle" /������Ȧ�����Ƴ�ʼ��
A_pos=zeros(1,dim); % Lioness A/��ʨA
A_score=inf; % change this to -inf for maximization problems/Ҫ���������⣬Ҫ�������Ϊ-inf

B_pos=zeros(1,dim); % Lioness B/��ʨB
B_score=inf; % change this to -inf for maximization problems/Ҫ���������⣬Ҫ�������Ϊ-inf

C_pos=zeros(1,dim); % Lioness C/��ʨC
C_score=inf; % change this to -inf for maximization problems/Ҫ���������⣬Ҫ�������Ϊ-inf

D_pos=zeros(1,dim); % Lioness D/��ʨD
D_score=inf; % change this to -inf for maximization problems/Ҫ���������⣬Ҫ�������Ϊ-inf
%%
% Main loop/��ѭ��
while Iter<Max_iter    
 
 for i=1:size(position,1)  
     
    % Check boundries/���߽�     
    Flag4ub=position(i,:)>ub;
    Flag4lb=position(i,:)<lb;    
    position(i,:)=(position(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;   
    
    % Calculate fitness/������Ӧ��ֵ    
    fitness(i,1)=fobj(position(i,:));
    
    % Update the position of Top_lioness/����Top_lioness��λ�� 
    if fitness(i,1)<Top_lioness_fit 
        Top_lioness_fit=fitness(i,1); 
        Top_lioness_pos=position(i,:);
    end 
    
 end
 
     %--------------------- Diversity  --------------------- 
     
%   mn=0;  % mn is used to record all distances in the population/mn���ڼ�¼��Ⱥ�����о���
%   for i=1:size(Prey,1)
%       for k=i+1:size(Prey,1)
%           mn=mn+1;
%           dis(mn)=distance(Prey(i,:),Prey(k,:));
%       end
%   end
%   dis_group(Iter+1)=mean(dis);
     
     %------------------- Memory saving ------------------- 
    
 if Iter==0
   fit_old=fitness;    Prey_old=position;
 end
     
  Inx=(fit_old<fitness);   
  Indx=repmat(Inx,1,dim);
  position=Indx.*Prey_old+~Indx.*position;     % Update the location of some search agents/���²������������λ��
  fitness=Inx.*fit_old+~Inx.*fitness;  % Update the fitness corresponding to the above positions/��������λ������Ӧ����Ӧ��ֵ
        
  fit_old=fitness;    Prey_old=position;

  %%% The structure of the "Central Circle"/������Ȧ���Ĺ���
  for i=1:size(position,1)  
        if fitness(i)<A_score 
            A_score=fitness(i); % Update Lioness A/���´�ʨA
            A_pos=position(i,:);
        end
        
        if fitness(i)>A_score && fitness(i)<B_score 
            B_score=fitness(i); % Update Lioness B/���´�ʨB
            B_pos=position(i,:);
        end
        
        if fitness(i)>A_score && fitness(i)>B_score && fitness(i)<C_score 
            C_score=fitness(i); % Update Lioness C/���´�ʨC
            C_pos=position(i,:);
        end
        
        if fitness(i)>A_score && fitness(i)>B_score && fitness(i)>C_score && fitness(i)<D_score
            D_score=fitness(i); % Update Lioness D/���´�ʨD
            D_pos=position(i,:);
        end
      
  end
  
 % The construction of the elite matrix/��Ӣ����Ĺ���
 r1_pos=[A_pos;B_pos;C_pos;D_pos];
 Top_lioness_pos1=mean(r1_pos(1:3,:),1);        % Eq.��7��     
 Elite_lioness=repmat(Top_lioness_pos1,SearchAgents_no,1);  
 % The construction of the elite matrix is completed/��Ӣ����Ĺ�������
 
 CF=(1-Iter/Max_iter)^(2*Iter/Max_iter);  
 a=2-Iter*((2)/Max_iter);                             
 RL=0.05*levy(SearchAgents_no,dim,1.5);   % Levy random number vector/Levy�����
 RB=randn(SearchAgents_no,dim);           % Brownian random number vector/Brownian�����
 q=rand();
 Q=0.5;

 if q <= Q   
 %% Teaming hunting/�ŶӺ������� 
      for i=1:size(position,1)
          
        for j=1:size(position,2)     
                 
            r1=rand(); % r1 is a random number in [0,1]
            r2=rand(); % r2 is a random number in [0,1]
            
            A1=2*a*r1-a; % Eq.��4��
            C1=2*r2;     % Eq.��5��

            D_A=abs(C1*A_pos(j)-position(i,j)); % Eq.��1��
            Xa=A_pos(j)-A1*D_A;             % Eq.��2��
                       
            r1=rand();
            r2=rand();
            
            A2=2*a*r1-a; % Eq.��4��                 
            C2=2*r2;     % Eq.��5��                   

            D_B=abs(C2*B_pos(j)-position(i,j)); % Eq.��1��
            Xb=B_pos(j)-A2*D_B;             % Eq.��2��
           
            r1=rand();
            r2=rand(); 
            
            A3=2*a*r1-a; % Eq.��4��
            C3=2*r2;     % Eq.��5�� 
            
            D_C=abs(C3*C_pos(j)-position(i,j));  % Eq.��1��
            Xc=C_pos(j)-A3*D_C;              % Eq.��2��
           
            r1=rand();
            r2=rand(); 
            
            A4=2*a*r1-a; % Eq.��4��
            C4=2*r2;     % Eq.��5�� 
           
            D_D=abs(C4*D_pos(j)-position(i,j));  % Eq.��1��
            Xd=D_pos(j)-A4*D_D;              % Eq.��2��
           
            X_ave=(Xa+Xb+Xc+Xd)/4;  % Eq.��3��
            E_Team=[Xa; Xb; Xc; Xd; X_ave];  % Eq.��6��
  %%% The construction of the "Center Circle" is completed/������Ȧ����������   
  
            M=0.9;
            if r2>M
                position(i,j)= E_Team(randi(size(E_Team,1)),:); % Eq.��8��
            else
                position(i,:)= E_Team(randi(size(E_Team,1)),:); % Eq.��9��
            end  
            
        end
      end
   
 else
 %% Elite hunting/��Ӣ����   
  for i=1:size(position,1)
     for j=1:size(position,2)        
       R=rand();
          %------------------ the early iteration (Eq.10) -----------------
       if Iter<Max_iter/3 
          step(i,j)=RB(i,j)*(Elite_lioness(i,j)-RB(i,j)*position(i,j));         
          position(i,j)=position(i,j)+P*R*step(i,j);                       
             
          %--------------- the middle iteration (Eqs.11 & 12)--------------
       elseif Iter>Max_iter/3 && Iter<2*Max_iter/3 
          
         if i>size(position,1)/2
            step(i,j)=RB(i,j)*(RB(i,j)*Elite_lioness(i,j)-position(i,j));
            position(i,j)=Elite_lioness(i,j)+P*CF*step(i,j); 
         else
            step(i,j)=RL(i,j)*(Elite_lioness(i,j)-RL(i,j)*position(i,j));                     
            position(i,j)=position(i,j)+P*R*step(i,j);  
         end  
         
         %----------------- the late iteration (Eq.13)---------------------
       else 
           
           step(i,j)=RL(i,j)*(RL(i,j)*Elite_lioness(i,j)-position(i,j)); 
           position(i,j)=Elite_lioness(i,j)+P*CF*step(i,j);  
    
       end 
       
     end 
  end 
  
 end
 %%             
  for i=1:size(position,1)  
        
    Flag4ub=position(i,:)>ub;  
    Flag4lb=position(i,:)<lb;  
    position(i,:)=(position(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
  
    fitness(i,1)=fobj(position(i,:));  
    
    if fitness(i,1)<Top_lioness_fit 
       Top_lioness_fit=fitness(i,1);
       Top_lioness_pos=position(i,:);
    end     
  end
        
     %----------------------  Memory saving ----------------------
    
  if Iter==0
    fit_old=fitness;    Prey_old=position;
  end
     
    Inx=(fit_old<fitness);
    Indx=repmat(Inx,1,dim);
    position=Indx.*Prey_old+~Indx.*position;
    fitness=Inx.*fit_old+~Inx.*fitness;
        
    fit_old=fitness;    Prey_old=position;

     %---------------- the impact of FADs (Eq.14) ----------------- 
                             
  if rand()<FADs  
     Ma=rand(SearchAgents_no,dim)<FADs;                                                                                              
     position=position+CF*((Xmin+rand(SearchAgents_no,dim).*(Xmax-Xmin)).*Ma);
  else
     r=rand();  Rn=size(position,1);
     step=(FADs*(1-r)+r)*(position(randperm(Rn),:)-position(randperm(Rn),:));
     position=position+step;
  end
                                                        
  Iter=Iter+1;  
  Convergence_curve(Iter)=Top_lioness_fit; 
%   mean(dis)    
end
  
   display(['The best optimal value of the objective funciton found by LsOA is : ', num2str(Top_lioness_fit)]);


