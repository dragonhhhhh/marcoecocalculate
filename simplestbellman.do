**bellman方程数值求解
**bellman方程形如：V(k) = ln(k^alpha) + V(k')

clear
mata: mata clear

mata
    // 模型参数
    alpha = 0.3    		// 资本系数
    beta = 0.9     		// 折现因子
    tol = 1e-6     		// 收敛容忍度
    max_iter = 1000 	// 最大迭代次数

    k_min = 0.01   		// 最小资本存量
    k_max = 10  		// 最大资本存量
    k_grid = range(k_min, k_max, 0.01)
    n_k = rows(k_grid)  // 网格点数量
	V = J(n_k, 1, 0) 	// 初始价值函数为0
	V_new = V       	// 用于更新的价值函数
	policy = J(n_k, 1, 0) // 策略函数，用于记录最优 k'
end




mata
for (iter = 1; iter <= max_iter; iter++){
		
	for (i = 1; i <= n_k; i++) {
		k = k_grid[i] // 当前资本存量
        utility = J(n_k, 1, -9999) // 初始化效用函数

        for (j = 1; j <= n_k; j++) {
			k_prime = k_grid[j] // 下期资本存量
            consumption = k^alpha - k_prime // 计算消费

            if (consumption > 0) {
				utility[j] = log(consumption) + beta * V[j] // 计算效用
            }
        }

        // 寻找最大效用对应的价值
        V_new[i] = max(utility)
		maxindex = .
		w = .
		maxindex(utility,1,maxindex,w)
		policy[i] = maxindex
        }

		if (max(abs(V_new - V)) < tol) {
        printf("在第 %f 次迭代后收敛\n", iter)
		break
        }
		
		V = V_new
}
end



//mata语言转为stata语言
mata
k_policy = k_grid[policy]

st_matrix("kpolicy",k_policy)
st_matrix("kgrid",k_grid)
st_matrix("V", V)

end


*matrix list kpolicy
*matrix list kgrid
*clear
svmat kgrid,names(k)
svmat kpolicy,names(kpolicy)
svmat V,names(v)


twoway (scatter kpolicy1 k1)
twoway (scatter v1 k1)
