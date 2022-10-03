set.seed(5)
n = 200
probability = 0.5
d = 20
N = 100

## Simulating a Erdos-Renyi Graph Adjacency matrix

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

eta = 0.05

Identity = diag(1,d)

Leta = solve(Sigma + (1/eta)*Identity)

Leta_over_eta = solve(eta*Sigma + Identity)


## The Data Mx

X = mvrnorm(n,numeric(d),Sigma)

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

est = matrix(0,nr=d,nc=d) ## Estimate of Leta_over_eta

## Estimating off-diagonal entries

for(i in 1:d){
	eid = e(i,d)
	if(i<d){
		for(j in (i+1):d){
			ejd = e(j,d)
			temp = 0
			for(NN in 1:N){
				tempY = mvrnorm(n,numeric(d),eta*Identity)
				temp = temp - 0.5*(f(tempY,X,eid+ejd,n)-f(tempY,X,eid-ejd,n))
			}
			temp = temp/N
			temp = temp/eta
			est[i,j] = temp
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
	temp = temp/eta
	est[i,i] = temp
}

## Distance between Leta_over_eta and est

DistanceF = norm(Leta_over_eta - est,"F")/d

## Finding distance between eigenvectors

dist1 = function(x,y){
		return(min(mean((x-y)^2),mean((x+y)^2))) # Because if x is an eigenvector so is -x
}

dist2 = function(X,Y,d_var=d){
	temp = 0
	for(i in 1:d_var){
		temp = temp + dist1(X[,i],Y[,i])
	}
	temp = temp / d_var
	temp = sqrt(temp)
	return(temp)
}

distance = dist2(eigen(est)$vectors,eigen(Leta_over_eta)$vectors,d)
