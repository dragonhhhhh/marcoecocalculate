********* Finite time *********

*The code in this document can implement Professor Jie Luo's MATLAB script `dynamic_finite_model.m`.





**# 函数、参数设定
program drop _all

**# Bookmark #4
mata:mata clear
mata
	
	alpha = 0.5
	beta  = 0.9
	delta = 0.1
	T 	  = 50
	k0 	  = 5
	eps   = 1
	step  = 0
	c1_min= 0.0001
	max_step=1000
	toleration = 0.0001
end

mata
	real scalar utility(real scalar c){
		return(log(c))
	}

	real scalar focu(real scalar c){
	return(1/c)
	}
	
    real scalar production(real scalar k, real scalar alpha){
		return(k^alpha)
	}
	
    real scalar focp(real scalar k, real scalar alpha){
		return(alpha*k^(alpha-1))
	}	
	
	
    real matrix inv_focp(real scalar y, real scalar alpha){
		return((y/alpha)^(1/(alpha-1)))
	}
	
	
	c1_max= production(k0,alpha) + (1-delta)*k0-0.0001
end


**# 稳态求解
mata
	ks = inv_focp(1/beta-1+delta,alpha) 
	ks
	cs = production(ks,alpha)-delta*ks
	cs
end


**# 有限期模型
mata
C = J(T+1, 1, .)
K = J(T+1, 1, .)
for (step = 1; step <= max_step; step++){
	K[1] = k0
	C[1] = (c1_min+c1_max)/2
    c1_high = 0
	
	for (t=1; t<=T; t++){
		if (c1_high==0){
			K[t+1] = production(K[t],alpha)+(1-delta)*K[t]- C[t]
			if (K[t+1]>0){
				C[t+1] = beta*C[t]*(focp(K[t+1],alpha) + 1-delta)
			}
			else{
				c1_high=1
			}
		}
	}
        
    if (c1_high==1){
		c1_max = C[1]
	}				 //初始c1太高了，调低c1
    else{
		c1_min = C[1]
	}				 //初始的c1太低了，调高c1
	step = step + 1
}

if (abs(K[51])<toleration = 0.0001){
	C[51]=.
	K[51]=.
}

C
K
end

matrix drop _all
mata
st_matrix("Consumption",C)
st_matrix("Kapital",K)
end

clear
svmat Consumption,names(C)
svmat Kapital,names(K)
gen t = _n

twoway (scatter C1 t) (line K1 t), ///
yline(1.807479224, lcolor(blue) lpattern(dash)) ///
yline(5.609418283, lcolor(red) lpattern(dash))
