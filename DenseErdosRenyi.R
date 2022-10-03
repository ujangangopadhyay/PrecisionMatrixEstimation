set.seed(5)
n = 200
probability = 0.5
d = 20
Adjacency = matrix(0,nr=d,nc=d)
for(i in 1:d){
	for(j in i:d){
		if(i != j){
		Adjacency[i,j] = rbinom(1,1,probability)
		Adjacency[j,i] = Adjacency[i,j]
		}
	}
}
Degree = matrix(0,nr=d,nc=d)
for(i in 1:d){
	Degree[i,i] = sum(Adjacency[i,])
}
Laplacian = Degree - Adjacency
Sigma = ginv(Laplacian)
eta = 0.01
Identity = diag(1,d)
Leta = solve(Sigma + (1/eta)*Identity)

X = mvrnorm(n,numeric(d),Sigma)
Y = mvrnorm(n,numeric(d),eta*Identity)

f = function(YY,XX,tt,nn){
	ss = 0
	for(kk in 1:nn){
	ss = ss + exp(1i*(( t(YY[kk,])%*%(XX[kk,]+tt) )))
	}
	if(Re(ss)<0){
		print(c(ss/nn,tt))
	}	
	ss = Re(ss)
	return(log(ss/nn))
}

e = function(ii,dd){
ee = numeric(dd)
ee[ii] = 1
return(ee)
}

N = 50

est = matrix(0,nr=d,nc=d)

for(i in 1:d){
	eid = e(i,d)
	if(i<d){
		for(j in (i+1):d){
		ejd = e(j,d)
		temp = 0
		for(NN in 1:N){
			tempY = Y = mvrnorm(n,numeric(d),eta*Identity)
			temp = temp - 0.5*(f(tempY,X,eid+ejd,n)-f(tempY,X,eid-ejd,n))
			}
			temp = temp/N
			est[i,j] = temp/eta
			est[j,i] = est[i,j]
		}
	}
}

## Estimating diagonal entries
for(i in 1:d){
	eid = e(i,d) 
	temp = 0
	for(NN in 1:N){
		tempY = mvrnorm(n,numeric(d),eta*Identity)
		temp = temp + 2*(f(tempY,X,numeric(d),n)-f(tempY,X,eid,n))
		}
	temp = temp/N
	est[i,i] = temp/eta
}


