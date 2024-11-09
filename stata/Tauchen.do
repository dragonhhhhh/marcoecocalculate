**Discretizing a first-order autoregressive process using the Tauchen method, for simplicity, only three states are distinguished.
**利用Tauchen方法对一阶自回归过程离散化, 为了简便起见，只区分3种状态

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

