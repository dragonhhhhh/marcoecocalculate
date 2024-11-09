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


%%  值函数(value function)迭代

%%%%%% 求解稳态
ks=inv_fp(1/beta-1+delta);
cs=f(ks)-delta*ks;
% kbar=fsolve(@(k) f(k)-delta*k,1000);
kbar=inv_fk(delta);


kmin=0.1;
kmax=kbar;
N=200;
k=linspace(kmin,kmax,N)';

%%%%%%% 数值解
%%%%%%% 初始化值函数
v0=zeros(N,1);
v1=zeros(N,1);

figure(1);
%%%%%%% 值函数迭代
step=0;
eps=10;
crit=10^(-6);
while eps>crit && step<1000
    
    if step<=0
        plot(k,v0,'b--','linewidth',2); hold on;
    end
    
    if step<=5 && step>=1
        plot(k,v0,'--','linewidth',1); hold on;
    end  

    if step<=25 && step>=21
        plot(k,v0,'--','linewidth',1); hold on;
    end     
    
    for i=1:N
        c=f(k(i))+(1-delta)*k(i) - k;
        index=c>0;
        max_v=max(u(c(index))+beta*v0(index));
        v1(i)=max_v;
    end
    eps=sqrt((sum((v0-v1).^2)));
    step=step+1;
    v0=v1;     
end
plot(k,v0,'r','linewidth',2); hold on;
xlabel('k'); ylabel('v(k)');


%%%%%%% 根据value function 算出policy function： c()和kp()
opt_c=zeros(N,1);
opt_kp=zeros(N,1);
for i=1:N
    c=f(k(i))+(1-delta)*k(i) - k;
    c(c<=0)=NaN;
    max_v=u(c)+beta*v0;
    index=find(max_v>=max(max_v),1);
    opt_c(i)=c(index);
    opt_kp(i)=index;
end



%%  根据policy function，画出最优路径
figure(2);
T=100;
ct=NaN(T+1,1);
kt=NaN(T+1,1);
kt_value=NaN(T+1,1);

%%%%%%%%%%%%  较低的初始k0
kt(1)=find(k>=2,1);
for t=1:T
    kt_value(t)=k(kt(t));
    ct(t)=opt_c(kt(t));
    kt(t+1)=opt_kp(kt(t));
end
k_grid=linspace(0,kbar,1000);
cs_k=f(k_grid)-delta*k_grid;
plot(kt_value,ct,'k--','linewidth',2);  hold on;
plot(k_grid,cs_k,'r','linewidth',2); hold on;
line([ks ks],[0,kbar],'linewidth',2); hold off;
axis([0 kbar 0 1]); grid on;
xlabel('k_t');  ylabel('c_t'); 

