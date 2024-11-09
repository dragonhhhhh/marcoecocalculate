**Unlike the previous value function solving problem, the technical change in the problem solved in this document is a three-state Markov process.This document is a synthesis of the previous three documents.

mata
mu = 0
rho = 0.9
sigma_epslilon = 0.1
sigma_z = sqrt(sigma_epslilon^2/(1-rho^2))
N=3
sigma_z
z = J(1,N,.)
z[1,1] = mu - sigma_z
z[1,N] = mu + sigma_z

for (i=1; i<=N; i++) {
z[1,i] = z[1,1] + (z[1,N] - z[1,1]) * (i - 1) / (N - 1)
}

P = J(N,N,.)
for (i=1; i<=N; i++) {
P[i,1] =  normal( ((z[1,1]+z[1,2])/2 - (1-rho)*mu - rho*z[1,i])/sigma_epslilon )
}

for (i=1; i<=N; i++) {
P[i,2] =  normal( ((z[1,2]+z[1,3])/2 - (1-rho)*mu - rho*z[1,i])/sigma_epslilon ) - normal( ((z[1,1]+z[1,2])/2 - (1-rho)*mu - rho*z[1,i])/sigma_epslilon )
}

for (i=1; i<=N; i++) {
P[i,3] =  1 - normal( ((z[1,2]+z[1,3])/2 - (1-rho)*mu - rho*z[1,i])/sigma_epslilon )
}
P

end


mata
Pi0 = (0\0.5\0.5)
Pi0_T = Pi0'
Pi1_T = Pi0_T*P
while(max(abs(Pi1_T-Pi0_T))>0.00001){	
	Pi0_T = Pi1_T
	Pi1_T = Pi0_T*P
}
Pi1_T
end


mata
    // 模型参数
    alpha = 0.3    // 资本系数
    beta = 0.96     // 折现因子
	delta = 0.08
	gamma = 2
    tol = 1e-6     // 收敛容忍度
    max_iter = 500 // 最大迭代次数
    // 定义资本存量网格
    k_min = 0.02   // 最小资本存量
    k_max = 10   // 最大资本存量
    k_grid = range(k_min, k_max, 0.02)
    n_k = rows(k_grid)     // 网格点数量
	V1 = J(n_k, 1, 0) // 初始价值函数为0
	V_new1 = V1       // 用于更新的价值函数
	V2 = J(n_k, 1, 0) // 初始价值函数为0
	V_new2 = V2       // 用于更新的价值函数
	V3 = J(n_k, 1, 0) // 初始价值函数为0
	V_new3 = V3       // 用于更新的价值函数
	policy1 = J(n_k, 1, 0) // 策略函数，用于记录最优 k'
	policy2 = J(n_k, 1, 0) // 策略函数，用于记录最优 k'
	policy3 = J(n_k, 1, 0) // 策略函数，用于记录最优 k'
end




mata
for (iter = 1; iter <= max_iter; iter++){
		
	for (i = 1; i <= n_k; i++) {
		k = k_grid[i] // 当前资本存量
        utility1 = J(n_k, 1, -9999) // 初始化效用函数
		utility2 = J(n_k, 1, -9999) // 初始化效用函数
		utility3 = J(n_k, 1, -9999) // 初始化效用函数

        for (j = 1; j <= n_k; j++) {
			k_prime = k_grid[j] // 下期资本存量
            
			consumption1 = exp(z[1,1])*k^alpha +(1-delta)*k- k_prime // 计算消费
            if (consumption1 > 0) {
				utility1[j] = (consumption1^(1-gamma)-1)/(1-gamma) + beta*(P[1,1]*V1[j] + P[1,2]*V2[j] + P[1,3]*V3[j]) // 计算效用
            }
			
            consumption2 = exp(z[1,2])*k^alpha +(1-delta)*k- k_prime // 计算消费
            if (consumption2 > 0) {
				utility2[j] = (consumption2^(1-gamma)-1)/(1-gamma) + beta*(P[2,1]*V1[j] + P[2,2]*V2[j] + P[2,3]*V3[j]) // 计算效用
            }
			
            consumption3 = exp(z[1,3])*k^alpha +(1-delta)*k- k_prime // 计算消费            
            if (consumption3 > 0) {
				utility3[j] = (consumption3^(1-gamma)-1)/(1-gamma) + beta*(P[3,1]*V1[j] + P[3,2]*V2[j] + P[3,3]*V3[j]) // 计算效用
            }			
			
			
        }	

        // 寻找最大效用对应的价值
        V_new1[i] = max(utility1)
		maxindex1 = .
		w1 = .
		maxindex(utility1,1,maxindex1,w1)
		policy1[i] = maxindex1
		
        V_new2[i] = max(utility2)
		maxindex2 = .
		w2 = .
		maxindex(utility2,1,maxindex2,w2)
		policy2[i] = maxindex2		

        V_new3[i] = max(utility3)
		maxindex3 = .
		w3 = .
		maxindex(utility3,1,maxindex3,w3)
		policy3[i] = maxindex3		
		
     }

	if (max(abs(V_new1 - V1)) < tol & max(abs(V_new2 - V3)) < tol & max(abs(V_new3 - V3)) < tol ) {
    printf("在第 %f 次迭代后收敛\n", iter)
	break
    }
		
	V1 = V_new1
	V2 = V_new2
	V3 = V_new3
}
end



//mata语言转为stata语言
mata
k_policy1 = k_grid[policy1]
k_policy2 = k_grid[policy2]
k_policy3 = k_grid[policy3]
st_matrix("kpolicy1",k_policy1)
st_matrix("kpolicy2",k_policy2)
st_matrix("kpolicy3",k_policy3)
st_matrix("kgrid",k_grid)
st_matrix("V1", V1)
st_matrix("V2", V2)
st_matrix("V3", V3)

end

clear
svmat kgrid,names(k)
svmat kpolicy1,names(kpolicy1)
svmat kpolicy2,names(kpolicy2)
svmat kpolicy3,names(kpolicy3)
svmat V1,names(v1)
svmat V2,names(v2)
svmat V3,names(v3)

twoway (scatter kpolicy11 k1) (scatter kpolicy21 k1) (scatter kpolicy31 k1)
twoway (scatter v11 k1) (scatter v21 k1) (scatter v31 k1)