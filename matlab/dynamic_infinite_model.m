clear; close all;

%%%% 参数设定
% utility function
u=@(c) log(c);      
up=@(c) 1./c;  
inv_up=@(mc) 1./mc;

% production function
alpha=0.5;  
f=@(k) k.^alpha; 
fp=@(k) alpha*k.^(alpha-1); 
inv_fp=@(mk) (mk/alpha).^(1/(alpha-1));
inv_fk=@(fk) fk.^(1/(alpha-1));

beta=0.9;  % discount factor
delta=0.5; % depreciation rate


%%  无穷期模型

%%%%%% 求解稳态
ks=inv_fp(1/beta-1+delta);
cs=f(ks)-delta*ks;
% kbar=fsolve(@(k) f(k)-delta*k,1000);
kbar=inv_fk(delta);



%% 根据任意初始值，画出路径图
close all;
T=100;
k0=2;
c1=0.9482;
c=NaN(T+1,1);
k=NaN(T+1,1);
k(1)=k0;   c(1)=c1; c1_high=0; 
for t=1:T
    if c1_high==0
        k(t+1)=f(k(t)) +(1-delta)*k(t) - c(t);
        if k(t+1)>0
            c(t+1)=inv_up( up(c(t)) / ( beta*( fp(k(t+1)) + 1-delta )) );
        else
            c1_high=1;
        end
    end
end

k_grid=linspace(0,kbar,1000);
cs_k=f(k_grid)-delta*k_grid;
plot(k,c,'b','linewidth',2);  hold on;
plot(k_grid,cs_k,'r','linewidth',2); hold on;
line([ks ks],[0,10],'linewidth',2); hold off;
axis([0 kbar 0 1]);
xlabel('k_t');  ylabel('c_t'); 


%% 寻找最优的路径
T=100;
k0=2;

eps=1; 
step=0;
c1_min=0.0001;
c1_max=f(k0)+(1-delta)*k0-0.0001;
while eps>0.0001 && step<1000
    c=NaN(T+1,1);
    k=NaN(T+1,1);
    k(1)=k0;
    % 选择c1，利用动态系统方程递推，使得k(T+1)=ks; 以及c(T+1)=cs;
    c1=(c1_min+c1_max)/2;
    c1_high=0;
    
    k(1)=k0;   c(1)=c1;
    for t=1:T
        if c1_high==0
            k(t+1)=f(k(t)) +(1-delta)*k(t) - c(t);
            if k(t+1)>0
                c(t+1)=inv_up( up(c(t)) / ( beta*( fp(k(t+1)) + 1-delta )) );
            else
                c1_high=1;
            end
        end
    end
    
    if c1_high==1
        c1_max=c1;  %%%%%% 初始c1太高了，调低c1
    else
        c1_min=c1;  %%%%%%   初始的c1太低了，调高c1
    end
    if isnan(k(T+1))
        eps=1;
    else
        eps=abs(k(T+1)/ks-1);
    end
    step=step+1;
end

k_grid=linspace(0,kbar,1000);
cs_k=f(k_grid)-delta*k_grid;
plot(k,c,'b','linewidth',2);  hold on;
plot(k_grid,cs_k,'r','linewidth',2); hold on;
line([ks ks],[0,10],'linewidth',2); hold off;
axis([0 kbar 0 1]);
xlabel('k_t');  ylabel('c_t'); 




