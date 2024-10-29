**Find the steady-state probability distribution of the Markov process
**找马尔科夫过程的稳态概率分布

mata
P = (0.9, 0.05, 0.05 \ 0.25, 0.65, 0.1\ 0, 0.3, 0.7)
Pi0 = (0\0.5\0.5)
Pi0_T = Pi0'
Pi1_T = Pi0_T*P
while(max(abs(Pi1_T-Pi0_T))>0.00001){	
	Pi0_T = Pi1_T
	Pi1_T = Pi0_T*P
}
Pi1_T
end