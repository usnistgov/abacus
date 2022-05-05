#### R program for design of calibration and quantitation experiment given in the write-up from Aug 28.
### It selects the best designs among all possible calibration designs up to size caln. Balanced designs of size 4, 
### 6,8,10,12,14,15,16,18,20,22,24, 28,30, 36, 40, 44, and 48 are considered. 


expdesign2<-function(caln,totn,a,b,vy,vx,beta,eta,vs,xnew,vi){
library(Matrix)
library(magic)
 
divisors <- function(x){
 #  Vector of numberes to test against
 y <- seq_len(x)
 #  Modulo division. If remainder is 0 that number is a divisor of x so return it 
 y[ x%%y == 0 ]}

caldesign104<-function(n,a,b,vy,vx,beta,eta,ny,vs,xnew,vi){
nq<-1
fy<-divisors(ny)
f<-divisors(n)
r<-1
if(vs>0 |vi>0){r=ny}
c<-diag(1,n)

#########
# design 1 (all I = 4 values are different, replicated J = 1 time)
I<-f[length(f)]
J<-1

th<-c(a,a+(b-a)/(I-1),a+2*(b-a)/(I-1),a+3*(b-a)/(I-1))
 

a1<-matrix(c(1,1,1,1,th[1],th[2],th[3],th[4]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),
           matrix((a+(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
           matrix((a+2*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
           matrix((a+3*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J))

d1temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),
             matrix(d11a[[2]],ncol=J,nrow=J),
             matrix(d11a[[3]],ncol=J,nrow=J),
             matrix(d11a[[4]],ncol=J,nrow=J))


d11a2 <- list(matrix(nq, nrow=J,ncol=J),
              matrix(nq,nrow=J,ncol=J),
              matrix(nq,nrow=J,ncol=J),
              matrix(nq,nrow=J,ncol=J))

d1temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J))

d1<-do.call(adiag,d1temp)
d12<-do.call(adiag,d1temp2)
sig1<-vy*c+beta^2*eta*d1+beta^2*vx*d12
sig1i<-solve(sig1)
b1<-t(a1)%*%sig1i%*%a1
b1i<-solve(b1)
i1<-I
j1<-J

#############
# design 2 (2 different values, replicated J = 2 times)
I <-f[length(f)-1]
J<-n/I
th<-c(rep(a,J),rep(a+(b-a)/(I-1),J))
a2<-matrix(c(1,1,1,1,th[1],th[2],th[3],th[4]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J))
d2temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J))
d2<-do.call(adiag,d2temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))
d2temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J))
d22<-do.call(adiag,d2temp2)


sig2<-vy*c+beta^2*eta*d2+beta^2*vx*d22
sig2i<-solve(sig2)
b2<-t(a2)%*%sig2i%*%a2
b2i<-solve(b2)
i2<-I
j2<-J
#############

eig1<-eigen(b1i)
tr1<-sum(eig1$values)
eig2<-eigen(b2i)
tr2<-sum(eig2$values)
trace<-c(tr1,tr2)
best<-which(trace==min(trace))

vb<-c((b1i[[2,2]])^0.5,(b2i[[2,2]])^0.5)
bestb<-which(vb==min(vb))

cq<-diag(1,nq)
thq<-c(rep(beta,nq))
th2q<-matrix(c(1,xnew),nrow=1)
a5q<-matrix(c(rep(beta,nq)), ncol=1)
varab1<-th2q%*%b1i%*%t(th2q)
d11aq1<-list(matrix(varab1, nrow=nq,ncol=nq))
dq1<-do.call(adiag,d11aq1)
sigq1<-dq1
sigq1i<-solve(sigq1)
bq1<-t(a5q)%*%sigq1i%*%a5q 
stuq1<-(1/bq1+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b1i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq1<-stuq1/xnew*100

varab2<-th2q%*%b2i%*%t(th2q)
d11aq2<-list(matrix(varab2, nrow=nq,ncol=nq))
dq2<-do.call(adiag,d11aq2)
sigq2<-dq2
sigq2i<-solve(sigq2)
bq2<-t(a5q)%*%sigq2i%*%a5q 
stuq2<-(1/bq2+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b2i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq2<-stuq2/xnew*100


finalrelst<-c(relstuq1,relstuq2)
finali<-c(i1,i2)
finalj<-c(j1,j2)
finbest<-which(finalrelst==min(finalrelst))
jbest<-finalj[finbest]
ibest<-finali[finbest]
finstd<-finalrelst[finbest]
totrep<-n+ny
samplrep<-ny/r
result<-list(samplrep=samplrep,nosamples=r,finstd=finstd,ibest=ibest,jbest=jbest,totrep=totrep)
return(result)
}

##################################################################################################
##################################################################################################

caldesign106<-function(n,a,b,vy,vx,beta,eta,ny,vs,xnew,vi){
f<-divisors(n)
c<-diag(1,n)
fy<-divisors(ny)
r<-1
if(vs>0 |vi>0){r=ny}
nq<-1
#########
# design 1 (all I = 6 values are different, replicated J = 1 time)
I<-f[length(f)]
J<-1

th<-c(a,a+(b-a)/(I-1),a+2*(b-a)/(I-1),a+3*(b-a)/(I-1),a+4*(b-a)/(I-1),
                  a+5*(b-a)/(I-1))
 

a1<-matrix(c(1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+2*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+4*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+5*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J))

d1temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J),matrix(d11a[[6]],ncol=J,nrow=J))
d1<-do.call(adiag,d1temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d1temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J),matrix(d11a2[[6]],ncol=J,nrow=J))
d12<-do.call(adiag,d1temp2)

sig1<-vy*c+beta^2*eta*d1+beta^2*vx*d12
sig1i<-solve(sig1)
b1<-t(a1)%*%sig1i%*%a1
b1i<-solve(b1)
i1<-I
j1<-J
#############
# design 2 (3 different values, replicated J = 2 times)
I <-f[length(f)-1]
J<-n/I
th<-c(rep(a,J),rep(a+(b-a)/(I-1),J),rep(a+2*(b-a)/(I-1),J))
a2<-matrix(c(1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J))
d2temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J))
d2<-do.call(adiag,d2temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))
d2temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J))
d22<-do.call(adiag,d2temp2)
sig2<-vy*c+beta^2*eta*d2+beta^2*vx*d22
sig2i<-solve(sig2)
b2<-t(a2)%*%sig2i%*%a2
b2i<-solve(b2)
i2<-I
j2<-J
#############

#############
# design 3 (I = 2 different values, replicated J=3 times)
I<-f[length(f)-2]
J<-n/I
th<-c(rep(a,J), rep(a+(b-a)/(I-1),J))
a3<-matrix(c(1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6]), ncol=2)
d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J))
d3temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J))
d3<-do.call(adiag,d3temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))
d3temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J))
d32<-do.call(adiag,d3temp2)

sig3<-vy*c+beta^2*eta*d3+beta^2*vx*d32
sig3i<-solve(sig3)
b3<-t(a3)%*%sig3i%*%a3
b3i<-solve(b3)
i3<-I
j3<-J
#################


deter<-c(det(b1i),det(b2i),det(b3i))
which(deter==min(deter))

eig1<-eigen(b1i)
tr1<-sum(eig1$values)
eig2<-eigen(b2i)
tr2<-sum(eig2$values)
eig3<-eigen(b3i)
tr3<-sum(eig3$values)
trace<-c(tr1,tr2,tr3)
best<-which(trace==min(trace))


vb<-c((b1i[[2,2]])^0.5,(b2i[[2,2]])^0.5,(b3i[[2,2]])^0.5)
bestb<-which(vb==min(vb))

cq<-diag(1,nq)
thq<-c(rep(beta,nq))
th2q<-matrix(c(1,xnew),nrow=1)
a5q<-matrix(c(rep(beta,nq)), ncol=1)
varab1<-th2q%*%b1i%*%t(th2q)
d11aq1<-list(matrix(varab1, nrow=nq,ncol=nq))
dq1<-do.call(adiag,d11aq1)
sigq1<-dq1
sigq1i<-solve(sigq1)
bq1<-t(a5q)%*%sigq1i%*%a5q 
stuq1<-(1/bq1+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b1i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq1<-stuq1/xnew*100

varab2<-th2q%*%b2i%*%t(th2q)
d11aq2<-list(matrix(varab2, nrow=nq,ncol=nq))
dq2<-do.call(adiag,d11aq2)
sigq2<-dq2
sigq2i<-solve(sigq2)
bq2<-t(a5q)%*%sigq2i%*%a5q 
stuq2<-(1/bq2+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b2i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq2<-stuq2/xnew*100

varab3<-th2q%*%b3i%*%t(th2q)
d11aq3<-list(matrix(varab3, nrow=nq,ncol=nq))
dq3<-do.call(adiag,d11aq3)
sigq3<-dq3
sigq3i<-solve(sigq3)
bq3<-t(a5q)%*%sigq3i%*%a5q 
stuq3<-(1/bq3/beta^2+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b3i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq3<-stuq3/xnew*100

finalrelst<-c(relstuq1,relstuq2,relstuq3)
finbest<-which(finalrelst==min(finalrelst))
finali<-c(i1,i2,i3)
finalj<-c(j1,j2,j3)
jbest<-finalj[finbest]
ibest<-finali[finbest]
totrep<-n+ny
finstd<-finalrelst[finbest]
samplrep<-ny/r
result<-list(samplrep=samplrep,nosamples=r,finstd=finstd,ibest=ibest,jbest=jbest,totrep=totrep)

return(result)
}

###############################################################################
###############################################################################
caldesign108<-function(n,a,b,vy,vx,beta,eta,ny,vs,xnew,vi){
f<-divisors(n)
fy<-divisors(ny)
r<-1
if(vs>0 |vi>0){r=ny}
c<-diag(1,n)
nq<-1
#########
# design 1 (all I = 8 values are different, replicated J = 1 time)
I<-f[length(f)]
J<-1

th<-c(a,a+(b-a)/(I-1),a+2*(b-a)/(I-1),a+3*(b-a)/(I-1),a+4*(b-a)/(I-1),
                  a+5*(b-a)/(I-1),a+6*(b-a)/(I-1),a+7*(b-a)/(I-1))
 

a1<-matrix(c(1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+2*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+4*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+5*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+6*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+7*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J))

d1temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J),matrix(d11a[[6]],ncol=J,nrow=J),matrix(d11a[[7]],ncol=J,nrow=J),matrix(d11a[[8]],ncol=J,nrow=J))

d1<-do.call(adiag,d1temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d1temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J),matrix(d11a2[[6]],ncol=J,nrow=J),matrix(d11a2[[7]],ncol=J,nrow=J),matrix(d11a2[[8]],ncol=J,nrow=J))

d12<-do.call(adiag,d1temp2)

sig1<-vy*c+beta^2*eta*d1+beta^2*vx*d12
sig1i<-solve(sig1)
b1<-t(a1)%*%sig1i%*%a1
b1i<-solve(b1)
i1<-I
j1<-J
#############
# design 2 (4 different values, replicated J = 2 times)
I <-f[length(f)-1]
J<-n/I
th<-c(rep(a,J),rep(a+(b-a)/(I-1),J),rep(a+2*(b-a)/(I-1),J),rep(a+3*(b-a)/(I-1),J))
a2<-matrix(c(1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J))

d2temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J))

d2<-do.call(adiag,d2temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d2temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J))

d22<-do.call(adiag,d2temp2)

sig2<-vy*c+beta^2*eta*d2+beta^2*vx*d22
sig2i<-solve(sig2)
b2<-t(a2)%*%sig2i%*%a2
b2i<-solve(b2)
i2<-I
j2<-J
#############
# design 3 (I = 2 different values, replicated J=4 times)
I<-f[length(f)-2]
J<-n/I
th<-c(rep(a,J), rep(a+(b-a)/(I-1),J))
a3<-matrix(c(1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8]), ncol=2)
d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J))
d3temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J))
d3<-do.call(adiag,d3temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))
d3temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J))
d32<-do.call(adiag,d3temp2)



sig3<-vy*c+beta^2*eta*d3+beta^2*vx*d32
sig3i<-solve(sig3)
b3<-t(a3)%*%sig3i%*%a3
b3i<-solve(b3)
i3<-I
j3<-J
#################


deter<-c(det(b1i),det(b2i),det(b3i))
which(deter==min(deter))

eig1<-eigen(b1i)
tr1<-sum(eig1$values)
eig2<-eigen(b2i)
tr2<-sum(eig2$values)
eig3<-eigen(b3i)
tr3<-sum(eig3$values)
trace<-c(tr1,tr2,tr3)
best<-which(trace==min(trace))


vb<-c((b1i[[2,2]])^0.5,(b2i[[2,2]])^0.5,(b3i[[2,2]])^0.5)
bestb<-which(vb==min(vb))

cq<-diag(1,nq)
thq<-c(rep(beta,nq))
th2q<-matrix(c(1,xnew),nrow=1)
a5q<-matrix(c(rep(beta,nq)), ncol=1)
varab1<-th2q%*%b1i%*%t(th2q)
d11aq1<-list(matrix(varab1, nrow=nq,ncol=nq))
dq1<-do.call(adiag,d11aq1)
sigq1<-dq1
sigq1i<-solve(sigq1)
bq1<-t(a5q)%*%sigq1i%*%a5q 
stuq1<-(1/bq1+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b1i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq1<-stuq1/xnew*100

varab2<-th2q%*%b2i%*%t(th2q)
d11aq2<-list(matrix(varab2, nrow=nq,ncol=nq))
dq2<-do.call(adiag,d11aq2)
sigq2<-dq2
sigq2i<-solve(sigq2)
bq2<-t(a5q)%*%sigq2i%*%a5q 
stuq2<-(1/bq2+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b2i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq2<-stuq2/xnew*100

varab3<-th2q%*%b3i%*%t(th2q)
d11aq3<-list(matrix(varab3, nrow=nq,ncol=nq))
dq3<-do.call(adiag,d11aq3)
sigq3<-dq3
sigq3i<-solve(sigq3)
bq3<-t(a5q)%*%sigq3i%*%a5q 
stuq3<-(1/bq3+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b3i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq3<-stuq3/xnew*100

finalrelst<-c(relstuq1,relstuq2,relstuq3)
finbest<-which(finalrelst==min(finalrelst))
finali<-c(i1,i2,i3)
finalj<-c(j1,j2,j3)
jbest<-finalj[finbest]
ibest<-finali[finbest]
totrep<-n+ny
finstd<-finalrelst[finbest]
samplrep<-ny/r
result<-list(samplrep=samplrep,nosamples=r,finstd=finstd,ibest=ibest,jbest=jbest,totrep=totrep)

return(result)
}
#######################################################################################
#######################################################################################


caldesign110<-function(n,a,b,vy,vx,beta,eta,ny,vs,xnew,vi){
f<-divisors(n)
fy<-divisors(ny)
r<-1
if(vs>0 |vi>0){r=ny}
c<-diag(1,n)
nq<-1
#########
# design 1 (all I = 10 values are different, replicated J = 1 time)
I<-f[length(f)]
J<-1

th<-c(a,a+(b-a)/(I-1),a+2*(b-a)/(I-1),a+3*(b-a)/(I-1),a+4*(b-a)/(I-1),
                  a+5*(b-a)/(I-1),a+6*(b-a)/(I-1),a+7*(b-a)/(I-1),a+8*(b-a)/(I-1),a+9*(b-a)/(I-1))
 

a1<-matrix(c(1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+2*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+4*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+5*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+6*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+7*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+8*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+9*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J))

d1temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J),matrix(d11a[[6]],ncol=J,nrow=J),matrix(d11a[[7]],ncol=J,nrow=J),matrix(d11a[[8]],ncol=J,nrow=J),matrix(d11a[[9]],ncol=J,nrow=J),
matrix(d11a[[10]],ncol=J,nrow=J))

d1<-do.call(adiag,d1temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d1temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J),
matrix(d11a2[[6]],ncol=J,nrow=J),matrix(d11a2[[7]],ncol=J,nrow=J),matrix(d11a2[[8]],ncol=J,nrow=J),matrix(d11a2[[9]],ncol=J,nrow=J),matrix(d11a2[[10]],ncol=J,nrow=J))

d12<-do.call(adiag,d1temp2)

sig1<-vy*c+beta^2*eta*d1+beta^2*vx*d12
sig1i<-solve(sig1)
b1<-t(a1)%*%sig1i%*%a1
b1i<-solve(b1)
i1<-I
j1<-J
#############
# design 2 (5 different values, replicated J = 2 times)
I <-f[length(f)-1]
J<-n/I
th<-c(rep(a,J),rep(a+(b-a)/(I-1),J),rep(a+2*(b-a)/(I-1),J),rep(a+3*(b-a)/(I-1),J),rep(a+4*(b-a)/(I-1),J))
 
a2<-matrix(c(1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),d55<-matrix((a+4*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J))

d2temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J))

d2<-do.call(adiag,d2temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),d55<-matrix(nq,nrow=J,ncol=J))

d2temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J))

d22<-do.call(adiag,d2temp2)

sig2<-vy*c+beta^2*eta*d2+beta^2*vx*d22
sig2i<-solve(sig2)
b2<-t(a2)%*%sig2i%*%a2
b2i<-solve(b2)
i2<-I
j2<-J
#############
# design 3 (I = 2 different values, replicated J = 5 times)
I<-f[length(f)-2]
J<-n/I
th<-c(rep(a,J),rep(a+(b-a)/(I-1),J))
a3<-matrix(c(1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10]), ncol=2)
d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J))
d3temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J))
d3<-do.call(adiag,d3temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))
d3temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J))
d32<-do.call(adiag,d3temp2)

sig3<-vy*c+beta^2*eta*d3+beta^2*vx*d32
sig3i<-solve(sig3)
b3<-t(a3)%*%sig3i%*%a3
b3i<-solve(b3)
i3<-I
j3<-J
##################

deter<-c(det(b1i),det(b2i),det(b3i))
which(deter==min(deter))

eig1<-eigen(b1i)
tr1<-sum(eig1$values)
eig2<-eigen(b2i)
tr2<-sum(eig2$values)
eig3<-eigen(b3i)
tr3<-sum(eig3$values)

trace<-c(tr1,tr2,tr3)
best<-which(trace==min(trace))

vb<-c((b1i[[2,2]])^0.5,(b2i[[2,2]])^0.5,(b3i[[2,2]])^0.5)
bestb<-which(vb==min(vb))

cq<-diag(1,nq)
thq<-c(rep(beta,nq))
th2q<-matrix(c(1,xnew),nrow=1)
a5q<-matrix(c(rep(beta,nq)), ncol=1)
varab1<-th2q%*%b1i%*%t(th2q)
d11aq1<-list(matrix(varab1, nrow=nq,ncol=nq))
dq1<-do.call(adiag,d11aq1)
sigq1<-dq1
sigq1i<-solve(sigq1)
bq1<-t(a5q)%*%sigq1i%*%a5q 
stuq1<-(1/bq1+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b1i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq1<-stuq1/xnew*100

varab2<-th2q%*%b2i%*%t(th2q)
d11aq2<-list(matrix(varab2, nrow=nq,ncol=nq))
dq2<-do.call(adiag,d11aq2)
sigq2<-dq2
sigq2i<-solve(sigq2)
bq2<-t(a5q)%*%sigq2i%*%a5q 
stuq2<-(1/bq2+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b2i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq2<-stuq2/xnew*100

varab3<-th2q%*%b3i%*%t(th2q)
d11aq3<-list(matrix(varab3, nrow=nq,ncol=nq))
dq3<-do.call(adiag,d11aq3)
sigq3<-dq3
sigq3i<-solve(sigq3)
bq3<-t(a5q)%*%sigq3i%*%a5q 
stuq3<-(1/bq3+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b3i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq3<-stuq3/xnew*100

finalrelst<-c(relstuq1,relstuq2,relstuq3)
finbest<-which(finalrelst==min(finalrelst))
finali<-c(i1,i2,i3)
finalj<-c(j1,j2,j3)
jbest<-finalj[finbest]
ibest<-finali[finbest]
totrep<-n+ny
finstd<-finalrelst[finbest]
samplrep<-ny/r
result<-list(samplrep=samplrep,nosamples=r,finstd=finstd,ibest=ibest,jbest=jbest,totrep=totrep)

return(result)
}
#######################################################
caldesign112<-function(n,a,b,vy,vx,beta,eta,ny,vs,xnew,vi){
f<-divisors(n)
fy<-divisors(ny)
r<-1
if(vs>0 |vi>0){r=ny}
c<-diag(1,n)
nq<-1
#########
# design 1 (all I = 12 values are different, replicated J = 1 time)
I<-f[length(f)]
J<-1

th<-c(a,a+(b-a)/(I-1),a+2*(b-a)/(I-1),a+3*(b-a)/(I-1),a+4*(b-a)/(I-1),
                  a+5*(b-a)/(I-1),a+6*(b-a)/(I-1),a+7*(b-a)/(I-1),a+8*(b-a)/(I-1),a+9*(b-a)/(I-1),a+10*(b-a)/(I-1),a+11*(b-a)/(I-1))
 

a1<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+2*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+4*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+5*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+6*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+7*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+8*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+9*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+10*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+11*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J))

d1temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J),matrix(d11a[[6]],ncol=J,nrow=J),matrix(d11a[[7]],ncol=J,nrow=J),matrix(d11a[[8]],ncol=J,nrow=J),matrix(d11a[[9]],ncol=J,nrow=J),
matrix(d11a[[10]],ncol=J,nrow=J),matrix(d11a[[11]],ncol=J,nrow=J),matrix(d11a[[12]],ncol=J,nrow=J))

d1<-do.call(adiag,d1temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d1temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J),matrix(d11a2[[6]],ncol=J,nrow=J),matrix(d11a2[[7]],ncol=J,nrow=J),matrix(d11a2[[8]],ncol=J,nrow=J),matrix(d11a2[[9]],ncol=J,nrow=J),
matrix(d11a2[[10]],ncol=J,nrow=J),matrix(d11a2[[11]],ncol=J,nrow=J),matrix(d11a2[[12]],ncol=J,nrow=J))

d12<-do.call(adiag,d1temp2)

sig1<-vy*c+beta^2*eta*d1+beta^2*vx*d12
sig1i<-solve(sig1)
b1<-t(a1)%*%sig1i%*%a1
b1i<-solve(b1)
i1<-I
j1<-J
#############
# design 2 (6 different values, replicated J = 2 times)
I <-f[length(f)-1]
J<-n/I
th<-c(rep(a,J),rep(a+(b-a)/(I-1),J),rep(a+2*(b-a)/(I-1),J),rep(a+3*(b-a)/(I-1),J),rep(a+4*(b-a)/(I-1),J),rep(a+5*(b-a)/(I-1),J))
 
a2<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),d55<-matrix((a+4*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),d66<-matrix((a+5*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J))

d2temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J),matrix(d11a[[6]],ncol=J,nrow=J))

d2<-do.call(adiag,d2temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),d55<-matrix(nq,nrow=J,ncol=J),d66<-matrix(nq,nrow=J,ncol=J))

d2temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J),matrix(d11a2[[6]],ncol=J,nrow=J))

d22<-do.call(adiag,d2temp2)


sig2<-vy*c+beta^2*eta*d2+beta^2*vx*d22
sig2i<-solve(sig2)
b2<-t(a2)%*%sig2i%*%a2
b2i<-solve(b2)
i2<-I
j2<-J
#############
# design 3 (I = 4 different values, replicated J=3 times)
I<-f[length(f)-2]
J<-n/I
th<-c(rep(a,J), rep(a+(b-a)/(I-1),J), rep(a+2*(b-a)/(I-1),J),rep(a+3*(b-a)/(I-1),J))
a3<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12]), ncol=2)
d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(I-1))^2,nrow=J,ncol=J))
d3temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J))
d3<-do.call(adiag,d3temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))
d3temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J))
d32<-do.call(adiag,d3temp2)


sig3<-vy*c+beta^2*eta*d3+beta^2*vx*d32
sig3i<-solve(sig3)
b3<-t(a3)%*%sig3i%*%a3
b3i<-solve(b3)
i3<-I
j3<-J
#################
# design 4 (I = 3 different values, replicated J = 4 times)
I<-f[length(f)-3]
J<-n/I
th<-c(rep(a,J), rep(a+(b-a)/(I-1),J), rep(a+2*(b-a)/(I-1),J))
 a4<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12]), ncol=2)
d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(I-1))^2,nrow=J,ncol=J))
d4temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J))
d4<-do.call(adiag,d4temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))
d4temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J))
d42<-do.call(adiag,d4temp2)

sig4<-vy*c+beta^2*eta*d4+beta^2*vx*d42
sig4i<-solve(sig4)
b4<-t(a4)%*%sig4i%*%a4
b4i<-solve(b4)
i4<-I
j4<-J
#################
# design 5 (I = 2 different values, replicated J = 6 times)
I<-f[length(f)-4]
J<-n/I
th<-c(rep(a,J), rep(a+(b-a)/(I-1),J))
a5<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12]), ncol=2)
d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J))
d5temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J))
d5<-do.call(adiag,d5temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))
d5temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J))
d52<-do.call(adiag,d5temp2)

sig5<-vy*c+beta^2*eta*d5+beta^2*vx*d52
sig5i<-solve(sig5)
b5<-t(a5)%*%sig5i%*%a5
b5i<-solve(b5)
i5<-I
j5<-J
##################

deter<-c(det(b1i),det(b2i),det(b3i),det(b4i),det(b5i))
which(deter==min(deter))

eig1<-eigen(b1i)
tr1<-sum(eig1$values)
eig2<-eigen(b2i)
tr2<-sum(eig2$values)
eig3<-eigen(b3i)
tr3<-sum(eig3$values)
eig4<-eigen(b4i)
tr4<-sum(eig4$values)
eig5<-eigen(b5i)
tr5<-sum(eig5$values)
trace<-c(tr1,tr2,tr3,tr4,tr5)
best<-which(trace==min(trace))

vb<-c((b1i[[2,2]])^0.5,(b2i[[2,2]])^0.5,(b3i[[2,2]])^0.5,(b4i[[2,2]])^0.5,(b5i[[2,2]])^0.5)
bestb<-which(vb==min(vb))

cq<-diag(1,nq)
thq<-c(rep(beta,nq))
th2q<-matrix(c(1,xnew),nrow=1)
a5q<-matrix(c(rep(beta,nq)), ncol=1)
varab1<-th2q%*%b1i%*%t(th2q)
d11aq1<-list(matrix(varab1, nrow=nq,ncol=nq))
dq1<-do.call(adiag,d11aq1)
sigq1<-dq1
sigq1i<-solve(sigq1)
bq1<-t(a5q)%*%sigq1i%*%a5q 
stuq1<-(1/bq1+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b1i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq1<-stuq1/xnew*100

varab2<-th2q%*%b2i%*%t(th2q)
d11aq2<-list(matrix(varab2, nrow=nq,ncol=nq))
dq2<-do.call(adiag,d11aq2)
sigq2<-dq2
sigq2i<-solve(sigq2)
bq2<-t(a5q)%*%sigq2i%*%a5q 
stuq2<-(1/bq2+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b2i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq2<-stuq2/xnew*100

varab3<-th2q%*%b3i%*%t(th2q)
d11aq3<-list(matrix(varab3, nrow=nq,ncol=nq))
dq3<-do.call(adiag,d11aq3)
sigq3<-dq3
sigq3i<-solve(sigq3)
bq3<-t(a5q)%*%sigq3i%*%a5q 
stuq3<-(1/bq3+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b3i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq3<-stuq3/xnew*100

varab4<-th2q%*%b4i%*%t(th2q)
d11aq4<-list(matrix(varab4, nrow=nq,ncol=nq))
dq4<-do.call(adiag,d11aq4)
sigq4<-dq4
sigq4i<-solve(sigq4)
bq4<-t(a5q)%*%sigq4i%*%a5q 
stuq4<-(1/bq4+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b4i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq4<-stuq4/xnew*100

varab5<-th2q%*%b5i%*%t(th2q)
d11aq5<-list(matrix(varab5, nrow=nq,ncol=nq))
dq5<-do.call(adiag,d11aq5)
sigq5<-dq5
sigq5i<-solve(sigq5)
bq5<-t(a5q)%*%sigq5i%*%a5q 
stuq5<-(1/bq5+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b5i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq5<-stuq5/xnew*100
finalrelst<-c(relstuq1,relstuq2,relstuq3,relstuq4,relstuq5)

finbest<-which(finalrelst==min(finalrelst))
finali<-c(i1,i2,i3,i4,i5)
finalj<-c(j1,j2,j3,j4,j5)
jbest<-finalj[finbest]
ibest<-finali[finbest]
totrep<-n+ny
finstd<-finalrelst[finbest]
samplrep<-ny/r
result<-list(samplrep=samplrep,nosamples=r,finstd=finstd,ibest=ibest,jbest=jbest,totrep=totrep)
return(result)
}
#######################################################
#######################################################
caldesign114<-function(n,a,b,vy,vx,beta,eta,ny,vs,xnew,vi){
f<-divisors(n)
fy<-divisors(ny)
r<-1
if(vs>0 |vi>0){r=ny}
c<-diag(1,n)
v1<-(beta^2*eta*a^2)^0.5
v2<-(beta^2*eta*b^2)^0.5
nq<-1

#########
# design 1 (all I = 14 values are different, replicated J = 1 time)


I<-f[length(f)]
J<-1
th<-c(a,a+(b-a)/(I-1),a+2*(b-a)/(I-1),a+3*(b-a)/(I-1),a+4*(b-a)/(I-1),
                  a+5*(b-a)/(I-1),a+6*(b-a)/(I-1),a+7*(b-a)/(I-1),a+8*(b-a)/(I-1),a+9*(b-a)/(I-1),a+10*(b-a)/(I-1),a+11*(b-a)/(I-1),
a+12*(b-a)/(I-1),b)
 
a1<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+2*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+4*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+5*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+6*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+7*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+8*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+9*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+10*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+11*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+12*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+13*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J))

d1temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J),matrix(d11a[[6]],ncol=J,nrow=J),matrix(d11a[[7]],ncol=J,nrow=J),matrix(d11a[[8]],ncol=J,nrow=J),matrix(d11a[[9]],ncol=J,nrow=J),
matrix(d11a[[10]],ncol=J,nrow=J),matrix(d11a[[11]],ncol=J,nrow=J),matrix(d11a[[12]],ncol=J,nrow=J),matrix(d11a[[13]],ncol=J,nrow=J),matrix(d11a[[14]],ncol=J,nrow=J))

d1<-do.call(adiag,d1temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d1temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J),matrix(d11a2[[6]],ncol=J,nrow=J),matrix(d11a2[[7]],ncol=J,nrow=J),matrix(d11a2[[8]],ncol=J,nrow=J),matrix(d11a2[[9]],ncol=J,nrow=J),
matrix(d11a2[[10]],ncol=J,nrow=J),matrix(d11a2[[11]],ncol=J,nrow=J),matrix(d11a2[[12]],ncol=J,nrow=J),matrix(d11a2[[13]],ncol=J,nrow=J),matrix(d11a2[[14]],ncol=J,nrow=J))

d12<-do.call(adiag,d1temp2)

sig1<-vy*c+beta^2*eta*d1+beta^2*vx*d12
sig1i<-solve(sig1)
b1<-t(a1)%*%sig1i%*%a1
b1i<-solve(b1)
i1<-I
j1<-J
#############
# design 2 (I=7 different values, replicated J = 2 times)
I <-f[length(f)-1]
J<-n/I
th<-c(rep(a,J),rep(a+(b-a)/(I-1),J),rep(a+2*(b-a)/(I-1),J),rep(a+3*(b-a)/(I-1),J),rep(a+4*(b-a)/(I-1),J),rep(a+5*(b-a)/(I-1),J),rep(a+6*(b-a)/(I-1),J))
 
a2<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+4*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+5*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+6*(b-a)/(I-1))^2,nrow=J,ncol=J))

d2temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J),matrix(d11a[[6]],ncol=J,nrow=J),matrix(d11a[[7]],ncol=J,nrow=J))

d2<-do.call(adiag,d2temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d2temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J),matrix(d11a2[[6]],ncol=J,nrow=J),matrix(d11a2[[7]],ncol=J,nrow=J))

d22<-do.call(adiag,d2temp2)


sig2<-vy*c+beta^2*eta*d2+beta^2*vx*d22
sig2i<-solve(sig2)
b2<-t(a2)%*%sig2i%*%a2
b2i<-solve(b2)
i2<-I
j2<-J
#############
# design 3 (I = 2 different values, replicated J=7 times)
I<-f[length(f)-2]
J<-n/I
th<-c(rep(a,J), rep(a+(b-a)/(I-1),J))
a3<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J))

d3temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J))
d3<-do.call(adiag,d3temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d3temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J))
d32<-do.call(adiag,d3temp2)


sig3<-vy*c+beta^2*eta*d3+beta^2*vx*d32
sig3i<-solve(sig3)
b3<-t(a3)%*%sig3i%*%a3
b3i<-solve(b3)
i3<-I
j3<-J
#################

deter<-c(det(b1i),det(b2i),det(b3i))
which(deter==min(deter))

eig1<-eigen(b1i)
tr1<-sum(eig1$values)
eig2<-eigen(b2i)
tr2<-sum(eig2$values)
eig3<-eigen(b3i)
tr3<-sum(eig3$values)
 
trace<-c(tr1,tr2,tr3)
best<-which(trace==min(trace))

vb<-c((b1i[[2,2]])^0.5,(b2i[[2,2]])^0.5,(b3i[[2,2]])^0.5)
bestb<-which(vb==min(vb))

cq<-diag(1,nq)
thq<-c(rep(beta,nq))
th2q<-matrix(c(1,xnew),nrow=1)
a5q<-matrix(c(rep(beta,nq)), ncol=1)
varab1<-th2q%*%b1i%*%t(th2q)
d11aq1<-list(matrix(varab1, nrow=nq,ncol=nq))
dq1<-do.call(adiag,d11aq1)
sigq1<-dq1
sigq1i<-solve(sigq1)
bq1<-t(a5q)%*%sigq1i%*%a5q 
stuq1<-(1/bq1+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b1i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq1<-stuq1/xnew*100

varab2<-th2q%*%b2i%*%t(th2q)
d11aq2<-list(matrix(varab2, nrow=nq,ncol=nq))
dq2<-do.call(adiag,d11aq2)
sigq2<-dq2
sigq2i<-solve(sigq2)
bq2<-t(a5q)%*%sigq2i%*%a5q 
stuq2<-(1/bq2+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b2i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq2<-stuq2/xnew*100

varab3<-th2q%*%b3i%*%t(th2q)
d11aq3<-list(matrix(varab3, nrow=nq,ncol=nq))
dq3<-do.call(adiag,d11aq3)
sigq3<-dq3
sigq3i<-solve(sigq3)
bq3<-t(a5q)%*%sigq3i%*%a5q 
stuq3<-(1/bq3+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b3i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq3<-stuq3/xnew*100
finalrelst<-c(relstuq1,relstuq2,relstuq3)
finbest<-which(finalrelst==min(finalrelst))
finali<-c(i1,i2,i3)
finalj<-c(j1,j2,j3)
jbest<-finalj[finbest]
ibest<-finali[finbest]
totrep<-n+ny
finstd<-finalrelst[finbest]
samplrep<-ny/r
result<-list(samplrep=samplrep,nosamples=r,finstd=finstd,ibest=ibest,jbest=jbest,totrep=totrep)


return(result)
}


#############################################################

#############################################################

#######################################################
caldesign115<-function(n,a,b,vy,vx,beta,eta,ny,vs,xnew,vi){
f<-divisors(n)
fy<-divisors(ny)
r<-1
if(vs>0 |vi>0){r=ny}
c<-diag(1,n)
v1<-(beta^2*eta*a^2)^0.5
v2<-(beta^2*eta*b^2)^0.5
nq<-1

#########
# design 1 (all I = 15 values are different, replicated J = 1 time)


I<-f[length(f)]
J<-1
th<-c(a,a+(b-a)/(I-1),a+2*(b-a)/(I-1),a+3*(b-a)/(I-1),a+4*(b-a)/(I-1),
                  a+5*(b-a)/(I-1),a+6*(b-a)/(I-1),a+7*(b-a)/(I-1),a+8*(b-a)/(I-1),a+9*(b-a)/(I-1),a+10*(b-a)/(I-1),a+11*(b-a)/(I-1),
a+12*(b-a)/(I-1),a+13*(b-a)/(I-1),b)
 
a1<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+2*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+4*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+5*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+6*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+7*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+8*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+9*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+10*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+11*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+12*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+13*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+14*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J))

d1temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J),matrix(d11a[[6]],ncol=J,nrow=J),matrix(d11a[[7]],ncol=J,nrow=J),matrix(d11a[[8]],ncol=J,nrow=J),matrix(d11a[[9]],ncol=J,nrow=J),
matrix(d11a[[10]],ncol=J,nrow=J),matrix(d11a[[11]],ncol=J,nrow=J),matrix(d11a[[12]],ncol=J,nrow=J),matrix(d11a[[13]],ncol=J,nrow=J),matrix(d11a[[14]],ncol=J,nrow=J),matrix(d11a[[15]],ncol=J,nrow=J))

d1<-do.call(adiag,d1temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J))

d1temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J),matrix(d11a2[[6]],ncol=J,nrow=J),matrix(d11a2[[7]],ncol=J,nrow=J),matrix(d11a2[[8]],ncol=J,nrow=J),matrix(d11a2[[9]],ncol=J,nrow=J),
matrix(d11a2[[10]],ncol=J,nrow=J),matrix(d11a2[[11]],ncol=J,nrow=J),matrix(d11a2[[12]],ncol=J,nrow=J),matrix(d11a2[[13]],ncol=J,nrow=J),matrix(d11a2[[14]],ncol=J,nrow=J),matrix(d11a2[[15]],ncol=J,nrow=J))

d12<-do.call(adiag,d1temp2)

sig1<-vy*c+beta^2*eta*d1+beta^2*vx*d12
sig1i<-solve(sig1)
b1<-t(a1)%*%sig1i%*%a1
b1i<-solve(b1)
i1<-I
j1<-J
#############
# design 2 (I=5 different values, replicated J = 3 times)
I <-f[length(f)-1]
J<-n/I
th<-c(rep(a,J),rep(a+(b-a)/(I-1),J),rep(a+2*(b-a)/(I-1),J),rep(a+3*(b-a)/(I-1),J),rep(a+4*(b-a)/(I-1),J))
 
a2<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+4*(b-a)/(I-1))^2,nrow=J,ncol=J))

d2temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J))

d2<-do.call(adiag,d2temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d2temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J))

d22<-do.call(adiag,d2temp2)

sig2<-vy*c+beta^2*eta*d2+beta^2*vx*d22
sig2i<-solve(sig2)
b2<-t(a2)%*%sig2i%*%a2
b2i<-solve(b2)
i2<-I
j2<-J
#############
# design 3 (I = 3 different values, replicated J=5 times)
I<-f[length(f)-2]
J<-n/I
th<-c(rep(a,J), rep(a+(b-a)/(I-1),J), rep(a+2*(b-a)/(I-1),J))
a3<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(I-1))^2,nrow=J,ncol=J))

d3temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J))

d3<-do.call(adiag,d3temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d3temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J))

d32<-do.call(adiag,d3temp2)

sig3<-vy*c+beta^2*eta*d3+beta^2*vx*d32
sig3i<-solve(sig3)
b3<-t(a3)%*%sig3i%*%a3
b3i<-solve(b3)
i3<-I
j3<-J
#################

deter<-c(det(b1i),det(b2i),det(b3i))
which(deter==min(deter))

eig1<-eigen(b1i)
tr1<-sum(eig1$values)
eig2<-eigen(b2i)
tr2<-sum(eig2$values)
eig3<-eigen(b3i)
tr3<-sum(eig3$values)
 

trace<-c(tr1,tr2,tr3)
best<-which(trace==min(trace))

vb<-c((b1i[[2,2]])^0.5,(b2i[[2,2]])^0.5,(b3i[[2,2]])^0.5)
bestb<-which(vb==min(vb))

cq<-diag(1,nq)
thq<-c(rep(beta,nq))
th2q<-matrix(c(1,xnew),nrow=1)
a5q<-matrix(c(rep(beta,nq)), ncol=1)
varab1<-th2q%*%b1i%*%t(th2q)
d11aq1<-list(matrix(varab1, nrow=nq,ncol=nq))
dq1<-do.call(adiag,d11aq1)
sigq1<-dq1
sigq1i<-solve(sigq1)
bq1<-t(a5q)%*%sigq1i%*%a5q 
stuq1<-(1/bq1+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b1i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq1<-stuq1/xnew*100

varab2<-th2q%*%b2i%*%t(th2q)
d11aq2<-list(matrix(varab2, nrow=nq,ncol=nq))
dq2<-do.call(adiag,d11aq2)
sigq2<-dq2
sigq2i<-solve(sigq2)
bq2<-t(a5q)%*%sigq2i%*%a5q 
stuq2<-(1/bq2+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b2i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq2<-stuq2/xnew*100

varab3<-th2q%*%b3i%*%t(th2q)
d11aq3<-list(matrix(varab3, nrow=nq,ncol=nq))
dq3<-do.call(adiag,d11aq3)
sigq3<-dq3
sigq3i<-solve(sigq3)
bq3<-t(a5q)%*%sigq3i%*%a5q 
stuq3<-(1/bq3+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b3i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq3<-stuq3/xnew*100
finalrelst<-c(relstuq1,relstuq2,relstuq3)
finbest<-which(finalrelst==min(finalrelst))

finali<-c(i1,i2,i3)
finalj<-c(j1,j2,j3)
jbest<-finalj[finbest]
ibest<-finali[finbest]
totrep<-n+ny
finstd<-finalrelst[finbest]
samplrep<-ny/r
result<-list(samplrep=samplrep,nosamples=r,finstd=finstd,ibest=ibest,jbest=jbest,totrep=totrep)


return(result)
}


#############################################################

#############################################################


caldesign116<-function(n,a,b,vy,vx,beta,eta,ny,vs,xnew,vi){
f<-divisors(n)
fy<-divisors(ny)
r<-1
if(vs>0 |vi>0){r=ny}
c<-diag(1,n)
v1<-(beta^2*eta*a^2)^0.5
v2<-(beta^2*eta*b^2)^0.5
nq<-1

#########
# design 1 (all I = 16 values are different, replicated J = 1 time)

I<-f[length(f)]
J<-1
th<-c(a,a+(b-a)/(I-1),a+2*(b-a)/(I-1),a+3*(b-a)/(I-1),a+4*(b-a)/(I-1),
                  a+5*(b-a)/(I-1),a+6*(b-a)/(I-1),a+7*(b-a)/(I-1),a+8*(b-a)/(I-1),a+9*(b-a)/(I-1),a+10*(b-a)/(I-1),a+11*(b-a)/(I-1),
a+12*(b-a)/(I-1),a+13*(b-a)/(I-1),a+14*(b-a)/(I-1),b)
 
a1<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+2*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+4*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+5*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+6*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+7*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+8*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+9*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+10*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+11*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+12*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+13*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+14*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+15*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J))

d1temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J),matrix(d11a[[6]],ncol=J,nrow=J),matrix(d11a[[7]],ncol=J,nrow=J),matrix(d11a[[8]],ncol=J,nrow=J),matrix(d11a[[9]],ncol=J,nrow=J),
matrix(d11a[[10]],ncol=J,nrow=J),matrix(d11a[[11]],ncol=J,nrow=J),matrix(d11a[[12]],ncol=J,nrow=J),matrix(d11a[[13]],ncol=J,nrow=J),matrix(d11a[[14]],ncol=J,nrow=J),matrix(d11a[[15]],ncol=J,nrow=J),matrix(d11a[[16]],ncol=J,nrow=J))
d1<-do.call(adiag,d1temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d1temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J),matrix(d11a2[[6]],ncol=J,nrow=J),matrix(d11a2[[7]],ncol=J,nrow=J),matrix(d11a2[[8]],ncol=J,nrow=J),matrix(d11a2[[9]],ncol=J,nrow=J),
matrix(d11a2[[10]],ncol=J,nrow=J),matrix(d11a2[[11]],ncol=J,nrow=J),matrix(d11a2[[12]],ncol=J,nrow=J),matrix(d11a2[[13]],ncol=J,nrow=J),matrix(d11a2[[14]],ncol=J,nrow=J),matrix(d11a2[[15]],ncol=J,nrow=J),matrix(d11a2[[16]],ncol=J,nrow=J))
d12<-do.call(adiag,d1temp2)


sig1<-vy*c+beta^2*eta*d1+beta^2*vx*d12
sig1i<-solve(sig1)
b1<-t(a1)%*%sig1i%*%a1
b1i<-solve(b1)
i1<-I
j1<-J
#############
# design 2 (I=8 different values, replicated J = 2 times)
I <-f[length(f)-1]
J<-n/I
th<-c(rep(a,J),rep(a+(b-a)/(I-1),J),rep(a+2*(b-a)/(I-1),J),rep(a+3*(b-a)/(I-1),J),rep(a+4*(b-a)/(I-1),J),rep(a+5*(b-a)/(I-1),J),rep(a+6*(b-a)/(I-1),J),rep(a+7*(b-a)/(I-1),J))
 
a2<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+4*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+5*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+6*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+7*(b-a)/(I-1))^2,nrow=J,ncol=J))

d2temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J),matrix(d11a[[6]],ncol=J,nrow=J),matrix(d11a[[7]],ncol=J,nrow=J),matrix(d11a[[8]],ncol=J,nrow=J))

d2<-do.call(adiag,d2temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d2temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J),matrix(d11a2[[6]],ncol=J,nrow=J),matrix(d11a2[[7]],ncol=J,nrow=J),matrix(d11a2[[8]],ncol=J,nrow=J))

d22<-do.call(adiag,d2temp2)

sig2<-vy*c+beta^2*eta*d2+beta^2*vx*d22
sig2i<-solve(sig2)
b2<-t(a2)%*%sig2i%*%a2
b2i<-solve(b2)
i2<-I
j2<-J
#############
# design 3 (I = 4 different values, replicated J=4 times)
I<-f[length(f)-2]
J<-n/I
th<-c(rep(a,J), rep(a+(b-a)/(I-1),J), rep(a+2*(b-a)/(I-1),J),rep(a+3*(b-a)/(I-1),J))
a3<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(I-1))^2,nrow=J,ncol=J))

d3temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J))

d3<-do.call(adiag,d3temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d3temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J))

d32<-do.call(adiag,d3temp2)

sig3<-vy*c+beta^2*eta*d3+beta^2*vx*d32
sig3i<-solve(sig3)
b3<-t(a3)%*%sig3i%*%a3
b3i<-solve(b3)
i3<-I
j3<-J
#################
# design 4 (I = 2 different values, replicated J = 8 times)
I<-f[length(f)-3]
J<-n/I
th<-c(rep(a,J), rep(a+(b-a)/(I-1),J))

a4<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J))

d4temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J))

d4<-do.call(adiag,d4temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d4temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J))

d42<-do.call(adiag,d4temp2)


sig4<-vy*c+beta^2*eta*d4+beta^2*vx*d42
sig4i<-solve(sig4)
b4<-t(a4)%*%sig4i%*%a4
b4i<-solve(b4)
i4<-I
j4<-J
#################
deter<-c(det(b1i),det(b2i),det(b3i),det(b4i))
which(deter==min(deter))

eig1<-eigen(b1i)
tr1<-sum(eig1$values)
eig2<-eigen(b2i)
tr2<-sum(eig2$values)
eig3<-eigen(b3i)
tr3<-sum(eig3$values)
eig4<-eigen(b4i)
tr4<-sum(eig4$values)

trace<-c(tr1,tr2,tr3,tr4)
best<-which(trace==min(trace))


vb<-c((b1i[[2,2]])^0.5,(b2i[[2,2]])^0.5,(b3i[[2,2]])^0.5,(b4i[[2,2]])^0.5)
bestb<-which(vb==min(vb))

cq<-diag(1,nq)
thq<-c(rep(beta,nq))
th2q<-matrix(c(1,xnew),nrow=1)
a5q<-matrix(c(rep(beta,nq)), ncol=1)
varab1<-th2q%*%b1i%*%t(th2q)
d11aq1<-list(matrix(varab1, nrow=nq,ncol=nq))
dq1<-do.call(adiag,d11aq1)
sigq1<-dq1
sigq1i<-solve(sigq1)
bq1<-t(a5q)%*%sigq1i%*%a5q 
stuq1<-(1/bq1+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b1i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq1<-stuq1/xnew*100

varab2<-th2q%*%b2i%*%t(th2q)
d11aq2<-list(matrix(varab2, nrow=nq,ncol=nq))
dq2<-do.call(adiag,d11aq2)
sigq2<-dq2
sigq2i<-solve(sigq2)
bq2<-t(a5q)%*%sigq2i%*%a5q 
stuq2<-(1/bq2+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b2i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq2<-stuq2/xnew*100

varab3<-th2q%*%b3i%*%t(th2q)
d11aq3<-list(matrix(varab3, nrow=nq,ncol=nq))
dq3<-do.call(adiag,d11aq3)
sigq3<-dq3
sigq3i<-solve(sigq3)
bq3<-t(a5q)%*%sigq3i%*%a5q 
stuq3<-(1/bq3+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b3i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq3<-stuq3/xnew*100

varab4<-th2q%*%b4i%*%t(th2q)
d11aq4<-list(matrix(varab4, nrow=nq,ncol=nq))
dq4<-do.call(adiag,d11aq4)
sigq4<-dq4
sigq4i<-solve(sigq4)
bq4<-t(a5q)%*%sigq4i%*%a5q 
stuq4<-(1/bq4+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b4i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq4<-stuq4/xnew*100



finalrelst<-c(relstuq1,relstuq2,relstuq3,relstuq4)
finbest<-which(finalrelst==min(finalrelst))
finali<-c(i1,i2,i3,i4)
finalj<-c(j1,j2,j3,j4)
jbest<-finalj[finbest]
ibest<-finali[finbest]
totrep<-n+ny
finstd<-finalrelst[finbest]
samplrep<-ny/r
result<-list(samplrep=samplrep,nosamples=r,finstd=finstd,ibest=ibest,jbest=jbest,totrep=totrep)


return(result)
}
#################################################################

#################################################################
#################################################################

caldesign118<-function(n,a,b,vy,vx,beta,eta,ny,vs,xnew,vi){
f<-divisors(n)
fy<-divisors(ny)
r<-1
if(vs>0 |vi>0){r=ny}
c<-diag(1,n)
v1<-(beta^2*eta*a^2)^0.5
v2<-(beta^2*eta*b^2)^0.5
nq<-1
#########
# design 1 (all I = 18 values are different, replicated J = 1 time)

J<-1
I<-f[length(f)]
th<-c(a,a+(b-a)/(I-1),a+2*(b-a)/(I-1),a+3*(b-a)/(I-1),a+4*(b-a)/(I-1),
                  a+5*(b-a)/(I-1),a+6*(b-a)/(I-1),a+7*(b-a)/(I-1),a+8*(b-a)/(I-1),a+9*(b-a)/(I-1),a+10*(b-a)/(I-1),a+11*(b-a)/(I-1),
a+12*(b-a)/(I-1),a+13*(b-a)/(I-1),a+14*(b-a)/(I-1),a+15*(b-a)/(I-1),a+16*(b-a)/(I-1),b)
 
a1<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+2*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+4*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+5*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+6*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+7*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+8*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+9*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+10*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+11*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+12*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+13*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+14*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+15*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+16*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+17*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J))

d1temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J),matrix(d11a[[6]],ncol=J,nrow=J),matrix(d11a[[7]],ncol=J,nrow=J),matrix(d11a[[8]],ncol=J,nrow=J),matrix(d11a[[9]],ncol=J,nrow=J),
matrix(d11a[[10]],ncol=J,nrow=J),matrix(d11a[[11]],ncol=J,nrow=J),matrix(d11a[[12]],ncol=J,nrow=J),matrix(d11a[[13]],ncol=J,nrow=J),matrix(d11a[[14]],ncol=J,nrow=J),matrix(d11a[[15]],ncol=J,nrow=J),matrix(d11a[[16]],ncol=J,nrow=J),matrix(d11a[[17]],ncol=J,nrow=J),
matrix(d11a[[18]],ncol=J,nrow=J))
d1<-do.call(adiag,d1temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d1temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J),matrix(d11a2[[6]],ncol=J,nrow=J),matrix(d11a2[[7]],ncol=J,nrow=J),matrix(d11a2[[8]],ncol=J,nrow=J),matrix(d11a2[[9]],ncol=J,nrow=J),
matrix(d11a2[[10]],ncol=J,nrow=J),matrix(d11a2[[11]],ncol=J,nrow=J),matrix(d11a2[[12]],ncol=J,nrow=J),matrix(d11a2[[13]],ncol=J,nrow=J),matrix(d11a2[[14]],ncol=J,nrow=J),matrix(d11a2[[15]],ncol=J,nrow=J),matrix(d11a2[[16]],ncol=J,nrow=J),matrix(d11a2[[17]],ncol=J,nrow=J),
matrix(d11a2[[18]],ncol=J,nrow=J))
d12<-do.call(adiag,d1temp2)



sig1<-vy*c+beta^2*eta*d1+beta^2*vx*d12

sig1i<-solve(sig1)
b1<-t(a1)%*%sig1i%*%a1
b1i<-solve(b1)
i1<-I
j1<-J
#############
# design 2 (I=9 different values, replicated J = 2 times)
I <-f[length(f)-1]
J<-n/I
th<-c(rep(a,J),rep(a+(b-a)/(I-1),J),rep(a+2*(b-a)/(I-1),J),rep(a+3*(b-a)/(I-1),J),rep(a+4*(b-a)/(I-1),J),rep(a+5*(b-a)/(I-1),J),rep(a+6*(b-a)/(I-1),J),rep(a+7*(b-a)/(I-1),J),rep(a+8*(b-a)/(I-1),J))
 
a2<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+4*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+5*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+6*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+7*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+8*(b-a)/(I-1))^2,nrow=J,ncol=J))

d2temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J),matrix(d11a[[6]],ncol=J,nrow=J),matrix(d11a[[7]],ncol=J,nrow=J),matrix(d11a[[8]],ncol=J,nrow=J),matrix(d11a[[9]],ncol=J,nrow=J))

d2<-do.call(adiag,d2temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d2temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J),matrix(d11a2[[6]],ncol=J,nrow=J),matrix(d11a2[[7]],ncol=J,nrow=J),matrix(d11a2[[8]],ncol=J,nrow=J),matrix(d11a2[[9]],ncol=J,nrow=J))

d22<-do.call(adiag,d2temp2)


sig2<-vy*c+beta^2*eta*d2+beta^2*vx*d22
sig2i<-solve(sig2)
b2<-t(a2)%*%sig2i%*%a2
b2i<-solve(b2)
i2<-I
j2<-J
#############
# design 3 (I = 6 different values, replicated J=3 times)
I<-f[length(f)-2]
J<-n/I
th<-c(rep(a,J), rep(a+(b-a)/(I-1),J), rep(a+2*(b-a)/(I-1),J),rep(a+3*(b-a)/(I-1),J),rep(a+4*(b-a)/(I-1),J), rep(a+5*(b-a)/(I-1),J))
a3<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+4*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+5*(b-a)/(I-1))^2,nrow=J,ncol=J))

d3temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J),matrix(d11a[[6]],ncol=J,nrow=J))

d3<-do.call(adiag,d3temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d3temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J),matrix(d11a2[[6]],ncol=J,nrow=J))

d32<-do.call(adiag,d3temp2)

sig3<-vy*c+beta^2*eta*d3+beta^2*vx*d32
sig3i<-solve(sig3)
b3<-t(a3)%*%sig3i%*%a3
b3i<-solve(b3)
i3<-I
j3<-J
#################
# design 4 (I = 3 different values, replicated J = 6 times)
I<-f[length(f)-3]
J<-n/I
th<-c(rep(a,J), rep(a+(b-a)/(I-1),J), rep(a+2*(b-a)/(I-1),J),rep(a+3*(b-a)/(I-1),J),rep(a+4*(b-a)/(I-1),J), rep(a+5*(b-a)/(I-1),J))
a4<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(I-1))^2,nrow=J,ncol=J))

d4temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J))

d4<-do.call(adiag,d4temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d4temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J))

d42<-do.call(adiag,d4temp2)

sig4<-vy*c+beta^2*eta*d4+beta^2*vx*d42
sig4i<-solve(sig4)
b4<-t(a4)%*%sig4i%*%a4
b4i<-solve(b4)
i4<-I
j4<-J
#################
# design 5 (I = 2 different values, replicated J = 8 times)
I<-f[length(f)-4]
J<-n/I
th<-c(rep(a,J), rep(a+(b-a)/(I-1),J))

a5<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J))

d5temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J))

d5<-do.call(adiag,d5temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d5temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J))

d52<-do.call(adiag,d5temp2)


sig5<-vy*c+beta^2*eta*d5+beta^2*vx*d52
sig5i<-solve(sig5)
b5<-t(a5)%*%sig5i%*%a5
b5i<-solve(b5)
i5<-I
j5<-J
##################



deter<-c(det(b1i),det(b2i),det(b3i),det(b4i),det(b5i))
which(deter==min(deter))

eig1<-eigen(b1i)
tr1<-sum(eig1$values)
eig2<-eigen(b2i)
tr2<-sum(eig2$values)
eig3<-eigen(b3i)
tr3<-sum(eig3$values)
eig4<-eigen(b4i)
tr4<-sum(eig4$values)
eig5<-eigen(b5i)
tr5<-sum(eig5$values)

 
trace<-c(tr1,tr2,tr3,tr4,tr5)
best<-which(trace==min(trace))


vb<-c((b1i[[2,2]])^0.5,(b2i[[2,2]])^0.5,(b3i[[2,2]])^0.5,(b4i[[2,2]])^0.5,(b5i[[2,2]])^0.5)
bestb<-which(vb==min(vb))

cq<-diag(1,nq)
thq<-c(rep(beta,nq))
th2q<-matrix(c(1,xnew),nrow=1)
a5q<-matrix(c(rep(beta,nq)), ncol=1)
varab1<-th2q%*%b1i%*%t(th2q)
d11aq1<-list(matrix(varab1, nrow=nq,ncol=nq))
dq1<-do.call(adiag,d11aq1)
sigq1<-dq1
sigq1i<-solve(sigq1)
bq1<-t(a5q)%*%sigq1i%*%a5q 
stuq1<-(1/bq1+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b1i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq1<-stuq1/xnew*100

varab2<-th2q%*%b2i%*%t(th2q)
d11aq2<-list(matrix(varab2, nrow=nq,ncol=nq))
dq2<-do.call(adiag,d11aq2)
sigq2<-dq2
sigq2i<-solve(sigq2)
bq2<-t(a5q)%*%sigq2i%*%a5q 
stuq2<-(1/bq2+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b2i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq2<-stuq2/xnew*100

varab3<-th2q%*%b3i%*%t(th2q)
d11aq3<-list(matrix(varab3, nrow=nq,ncol=nq))
dq3<-do.call(adiag,d11aq3)
sigq3<-dq3
sigq3i<-solve(sigq3)
bq3<-t(a5q)%*%sigq3i%*%a5q 
stuq3<-(1/bq3+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b3i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq3<-stuq3/xnew*100

varab4<-th2q%*%b4i%*%t(th2q)
d11aq4<-list(matrix(varab4, nrow=nq,ncol=nq))
dq4<-do.call(adiag,d11aq4)
sigq4<-dq4
sigq4i<-solve(sigq4)
bq4<-t(a5q)%*%sigq4i%*%a5q 
stuq4<-(1/bq4+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b4i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq4<-stuq4/xnew*100

varab5<-th2q%*%b5i%*%t(th2q)
d11aq5<-list(matrix(varab5, nrow=nq,ncol=nq))
dq5<-do.call(adiag,d11aq5)
sigq5<-dq5
sigq5i<-solve(sigq5)
bq5<-t(a5q)%*%sigq5i%*%a5q 
stuq5<-(1/bq5+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b5i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq5<-stuq5/xnew*100


finalrelst<-c(relstuq1,relstuq2,relstuq3,relstuq4,relstuq5)
finbest<-which(finalrelst==min(finalrelst))

finali<-c(i1,i2,i3,i4,i5)
finalj<-c(j1,j2,j3,j4,j5)
jbest<-finalj[finbest]
ibest<-finali[finbest]
totrep<-n+ny
finstd<-finalrelst[finbest]
samplrep<-ny/r
result<-list(samplrep=samplrep,nosamples=r,finstd=finstd,ibest=ibest,jbest=jbest,totrep=totrep)

return(result)
}

#########################################################################################
#########################################################################################
caldesign120<-function(n,a,b,vy,vx,beta,eta,ny,vs,xnew,vi){
f<-divisors(n)
fy<-divisors(ny)
r<-1
if(vs>0 |vi>0){r=ny}
c<-diag(1,n)
v1<-(beta^2*eta*a^2)^0.5
v2<-(beta^2*eta*b^2)^0.5
nq<-1
#########
# design 1 (all I = 20 values are different, replicated J = 1 time)

J<-1
I<-f[length(f)]
th<-c(a,a+(b-a)/(I-1),a+2*(b-a)/(I-1),a+3*(b-a)/(I-1),a+4*(b-a)/(I-1),
                  a+5*(b-a)/(I-1),a+6*(b-a)/(I-1),a+7*(b-a)/(I-1),a+8*(b-a)/(I-1),a+9*(b-a)/(I-1),a+10*(b-a)/(I-1),a+11*(b-a)/(I-1),
a+12*(b-a)/(I-1),a+13*(b-a)/(I-1),a+14*(b-a)/(I-1),a+15*(b-a)/(I-1),
                  a+16*(b-a)/(I-1),a+17*(b-a)/(I-1),a+18*(b-a)/(I-1),b)
 
a1<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],th[20]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+2*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+4*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+5*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+6*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+7*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+8*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+9*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+10*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+11*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+12*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+13*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+14*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+15*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+16*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+17*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+18*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+19*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J))

d1temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J),matrix(d11a[[6]],ncol=J,nrow=J),matrix(d11a[[7]],ncol=J,nrow=J),matrix(d11a[[8]],ncol=J,nrow=J),matrix(d11a[[9]],ncol=J,nrow=J),
matrix(d11a[[10]],ncol=J,nrow=J),matrix(d11a[[11]],ncol=J,nrow=J),matrix(d11a[[12]],ncol=J,nrow=J),matrix(d11a[[13]],ncol=J,nrow=J),matrix(d11a[[14]],ncol=J,nrow=J),matrix(d11a[[15]],ncol=J,nrow=J),matrix(d11a[[16]],ncol=J,nrow=J),matrix(d11a[[17]],ncol=J,nrow=J),
matrix(d11a[[18]],ncol=J,nrow=J),matrix(d11a[[19]],ncol=J,nrow=J),matrix(d11a[[20]]))
d1<-do.call(adiag,d1temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d1temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J),matrix(d11a2[[6]],ncol=J,nrow=J),matrix(d11a2[[7]],ncol=J,nrow=J),matrix(d11a2[[8]],ncol=J,nrow=J),matrix(d11a2[[9]],ncol=J,nrow=J),
matrix(d11a2[[10]],ncol=J,nrow=J),matrix(d11a2[[11]],ncol=J,nrow=J),matrix(d11a2[[12]],ncol=J,nrow=J),matrix(d11a2[[13]],ncol=J,nrow=J),matrix(d11a2[[14]],ncol=J,nrow=J),matrix(d11a2[[15]],ncol=J,nrow=J),matrix(d11a2[[16]],ncol=J,nrow=J),matrix(d11a2[[17]],ncol=J,nrow=J),
matrix(d11a2[[18]],ncol=J,nrow=J),matrix(d11a2[[19]],ncol=J,nrow=J),matrix(d11a2[[20]]))
d12<-do.call(adiag,d1temp2)


sig1<-vy*c+beta^2*eta*d1+beta^2*vx*d12

sig1i<-solve(sig1)
b1<-t(a1)%*%sig1i%*%a1
b1i<-solve(b1)
i1<-I
j1<-J
#############
# design 2 (I=10 different values, replicated J = 2 times)
I <-f[length(f)-1]
J<-n/I
th<-c(rep(a,J),rep(a+(b-a)/(I-1),J),rep(a+2*(b-a)/(I-1),J),rep(a+3*(b-a)/(I-1),J),rep(a+4*(b-a)/(I-1),J),rep(a+5*(b-a)/(I-1),J),rep(a+6*(b-a)/(I-1),J),rep(a+7*(b-a)/(I-1),J),rep(a+8*(b-a)/(I-1),J),rep(a+9*(b-a)/(I-1),J))
 
a2<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],th[20]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+4*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+5*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+6*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+7*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+8*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+9*(b-a)/(I-1))^2,nrow=J,ncol=J))

d2temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J),matrix(d11a[[6]],ncol=J,nrow=J),matrix(d11a[[7]],ncol=J,nrow=J),matrix(d11a[[8]],ncol=J,nrow=J),matrix(d11a[[9]],ncol=J,nrow=J),matrix(d11a[[10]],ncol=J,nrow=J))

d2<-do.call(adiag,d2temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d2temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J),matrix(d11a2[[6]],ncol=J,nrow=J),matrix(d11a2[[7]],ncol=J,nrow=J),matrix(d11a2[[8]],ncol=J,nrow=J),matrix(d11a2[[9]],ncol=J,nrow=J),matrix(d11a2[[10]],ncol=J,nrow=J))

d22<-do.call(adiag,d2temp2)

sig2<-vy*c+beta^2*eta*d2+beta^2*vx*d22
sig2i<-solve(sig2)
b2<-t(a2)%*%sig2i%*%a2
b2i<-solve(b2)
i2<-I
j2<-J
#############
# design 3 (I = 5 different values, replicated J=4 times)
I<-f[length(f)-2]
J<-n/I
th<-c(rep(a,J), rep(a+(b-a)/(I-1),J), rep(a+2*(b-a)/(I-1),J),rep(a+3*(b-a)/(I-1),J),rep(a+4*(b-a)/(I-1),J))
a3<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],th[20]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+4*(b-a)/(I-1))^2,nrow=J,ncol=J))

d3temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J))

d3<-do.call(adiag,d3temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d3temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J))

d32<-do.call(adiag,d3temp2)

sig3<-vy*c+beta^2*eta*d3+beta^2*vx*d32
sig3i<-solve(sig3)
b3<-t(a3)%*%sig3i%*%a3
b3i<-solve(b3)
i3<-I
j3<-J
#################
# design 4 (I = 4 different values, replicated J = 5 times)
I<-f[length(f)-3]
J<-n/I
th<-c(rep(a,J), rep(a+(b-a)/(I-1),J), rep(a+2*(b-a)/(I-1),J),rep(a+3*(b-a)/(I-1),J))
a4<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],th[20]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(I-1))^2,nrow=J,ncol=J))

d4temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J))

d4<-do.call(adiag,d4temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d4temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J))

d42<-do.call(adiag,d4temp2)


sig4<-vy*c+beta^2*eta*d4+beta^2*vx*d42
sig4i<-solve(sig4)
b4<-t(a4)%*%sig4i%*%a4
b4i<-solve(b4)
i4<-I
j4<-J
#################
# design 5 (I = 2 different values, replicated J = 10 times)
I<-f[length(f)-4]
J<-n/I
th<-c(rep(a,J), rep(a+(b-a)/(I-1),J))

a5<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],th[20]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J))

d5temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J))

d5<-do.call(adiag,d5temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d5temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J))

d52<-do.call(adiag,d5temp2)


sig5<-vy*c+beta^2*eta*d5+beta^2*vx*d52
sig5i<-solve(sig5)
b5<-t(a5)%*%sig5i%*%a5
b5i<-solve(b5)
i5<-I
j5<-J
##################


deter<-c(det(b1i),det(b2i),det(b3i),det(b4i),det(b5i))
which(deter==min(deter))

eig1<-eigen(b1i)
tr1<-sum(eig1$values)
eig2<-eigen(b2i)
tr2<-sum(eig2$values)
eig3<-eigen(b3i)
tr3<-sum(eig3$values)
eig4<-eigen(b4i)
tr4<-sum(eig4$values)
eig5<-eigen(b5i)
tr5<-sum(eig5$values)


 
trace<-c(tr1,tr2,tr3,tr4,tr5)
best<-which(trace==min(trace))


vb<-c((b1i[[2,2]])^0.5,(b2i[[2,2]])^0.5,(b3i[[2,2]])^0.5,(b4i[[2,2]])^0.5,(b5i[[2,2]])^0.5)
bestb<-which(vb==min(vb))

cq<-diag(1,nq)
thq<-c(rep(beta,nq))
th2q<-matrix(c(1,xnew),nrow=1)
a5q<-matrix(c(rep(beta,nq)), ncol=1)
varab1<-th2q%*%b1i%*%t(th2q)
d11aq1<-list(matrix(varab1, nrow=nq,ncol=nq))
dq1<-do.call(adiag,d11aq1)
sigq1<-dq1
sigq1i<-solve(sigq1)
bq1<-t(a5q)%*%sigq1i%*%a5q 
stuq1<-(1/bq1+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b1i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq1<-stuq1/xnew*100

varab2<-th2q%*%b2i%*%t(th2q)
d11aq2<-list(matrix(varab2, nrow=nq,ncol=nq))
dq2<-do.call(adiag,d11aq2)
sigq2<-dq2
sigq2i<-solve(sigq2)
bq2<-t(a5q)%*%sigq2i%*%a5q 
stuq2<-(1/bq2+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b2i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq2<-stuq2/xnew*100

varab3<-th2q%*%b3i%*%t(th2q)
d11aq3<-list(matrix(varab3, nrow=nq,ncol=nq))
dq3<-do.call(adiag,d11aq3)
sigq3<-dq3
sigq3i<-solve(sigq3)
bq3<-t(a5q)%*%sigq3i%*%a5q 
stuq3<-(1/bq3+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b3i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq3<-stuq3/xnew*100

varab4<-th2q%*%b4i%*%t(th2q)
d11aq4<-list(matrix(varab4, nrow=nq,ncol=nq))
dq4<-do.call(adiag,d11aq4)
sigq4<-dq4
sigq4i<-solve(sigq4)
bq4<-t(a5q)%*%sigq4i%*%a5q 
stuq4<-(1/bq4+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b4i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq4<-stuq4/xnew*100

varab5<-th2q%*%b5i%*%t(th2q)
d11aq5<-list(matrix(varab5, nrow=nq,ncol=nq))
dq5<-do.call(adiag,d11aq5)
sigq5<-dq5
sigq5i<-solve(sigq5)
bq5<-t(a5q)%*%sigq5i%*%a5q 
stuq5<-(1/bq5+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b5i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq5<-stuq5/xnew*100


finalrelst<-c(relstuq1,relstuq2,relstuq3,relstuq4,relstuq5)
finbest<-which(finalrelst==min(finalrelst))
finali<-c(i1,i2,i3,i4,i5)
finalj<-c(j1,j2,j3,j4,j5)
jbest<-finalj[finbest]
ibest<-finali[finbest]
totrep<-n+ny
finstd<-finalrelst[finbest]
samplrep<-ny/r
result<-list(samplrep=samplrep,nosamples=r,finstd=finstd,ibest=ibest,jbest=jbest,totrep=totrep)
return(result)
}

#########################################################################################
#########################################################################################
caldesign122<-function(n,a,b,vy,vx,beta,eta,ny,vs,xnew,vi){
f<-divisors(n)
fy<-divisors(ny)
r<-1
if(vs>0 |vi>0){r=ny}
c<-diag(1,n)
v1<-(beta^2*eta*a^2)^0.5
v2<-(beta^2*eta*b^2)^0.5
nq<-1
#########
# design 1 (all I = 22 values are different, replicated J = 1 time)

J<-1
I<-f[length(f)]
th<-c(a,a+(b-a)/(I-1),a+2*(b-a)/(I-1),a+3*(b-a)/(I-1),a+4*(b-a)/(I-1),
                  a+5*(b-a)/(I-1),a+6*(b-a)/(I-1),a+7*(b-a)/(I-1),a+8*(b-a)/(I-1),a+9*(b-a)/(I-1),a+10*(b-a)/(I-1),a+11*(b-a)/(I-1),
a+12*(b-a)/(I-1),a+13*(b-a)/(I-1),a+14*(b-a)/(I-1),a+15*(b-a)/(I-1),
                  a+16*(b-a)/(I-1),a+17*(b-a)/(I-1),a+18*(b-a)/(I-1),a+19*(b-a)/(I-1),a+20*(b-a)/(I-1),b)
 
a1<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],th[20],th[21],th[22]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+2*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+4*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+5*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+6*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+7*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+8*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+9*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+10*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+11*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+12*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+13*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+14*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+15*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+16*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+17*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+18*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+19*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+20*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+21*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J))

d1temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J),matrix(d11a[[6]],ncol=J,nrow=J),matrix(d11a[[7]],ncol=J,nrow=J),matrix(d11a[[8]],ncol=J,nrow=J),matrix(d11a[[9]],ncol=J,nrow=J),
matrix(d11a[[10]],ncol=J,nrow=J),matrix(d11a[[11]],ncol=J,nrow=J),matrix(d11a[[12]],ncol=J,nrow=J),matrix(d11a[[13]],ncol=J,nrow=J),matrix(d11a[[14]],ncol=J,nrow=J),matrix(d11a[[15]],ncol=J,nrow=J),matrix(d11a[[16]],ncol=J,nrow=J),matrix(d11a[[17]],ncol=J,nrow=J),
matrix(d11a[[18]],ncol=J,nrow=J),matrix(d11a[[19]],ncol=J,nrow=J),matrix(d11a[[20]],ncol=J,nrow=J),matrix(d11a[[21]],ncol=J,nrow=J),
matrix(d11a[[22]],ncol=J,nrow=J))
d1<-do.call(adiag,d1temp)


d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d1temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J),matrix(d11a2[[6]],ncol=J,nrow=J),matrix(d11a2[[7]],ncol=J,nrow=J),matrix(d11a2[[8]],ncol=J,nrow=J),matrix(d11a2[[9]],ncol=J,nrow=J),
matrix(d11a2[[10]],ncol=J,nrow=J),matrix(d11a2[[11]],ncol=J,nrow=J),matrix(d11a2[[12]],ncol=J,nrow=J),matrix(d11a2[[13]],ncol=J,nrow=J),matrix(d11a2[[14]],ncol=J,nrow=J),matrix(d11a2[[15]],ncol=J,nrow=J),matrix(d11a2[[16]],ncol=J,nrow=J),matrix(d11a2[[17]],ncol=J,nrow=J),
matrix(d11a2[[18]],ncol=J,nrow=J),matrix(d11a2[[19]],ncol=J,nrow=J),matrix(d11a2[[20]],ncol=J,nrow=J),matrix(d11a2[[21]],ncol=J,nrow=J),
matrix(d11a2[[22]],ncol=J,nrow=J))
d12<-do.call(adiag,d1temp2)


sig1<-vy*c+beta^2*eta*d1+beta^2*vx*d12

sig1i<-solve(sig1)
b1<-t(a1)%*%sig1i%*%a1
b1i<-solve(b1)
i1<-I
j1<-J
#############
# design 2 (I=11 different values, replicated J = 2 times)
I <-f[length(f)-1]
J<-n/I
th<-c(rep(a,J),rep(a+(b-a)/(I-1),J),rep(a+2*(b-a)/(I-1),J),rep(a+3*(b-a)/(I-1),J),rep(a+4*(b-a)/(I-1),J),rep(a+5*(b-a)/(I-1),J),rep(a+6*(b-a)/(I-1),J),rep(a+7*(b-a)/(I-1),J),rep(a+8*(b-a)/(I-1),J),rep(a+9*(b-a)/(I-1),J),rep(a+10*(b-a)/(I-1),J))
 
a2<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],th[20],th[21],th[22]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+4*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+5*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+6*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+7*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+8*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+9*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+10*(b-a)/(I-1))^2,nrow=J,ncol=J))

d2temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J),matrix(d11a[[6]],ncol=J,nrow=J),matrix(d11a[[7]],ncol=J,nrow=J),matrix(d11a[[8]],ncol=J,nrow=J),matrix(d11a[[9]],ncol=J,nrow=J),matrix(d11a[[10]],ncol=J,nrow=J),matrix(d11a[[11]],ncol=J,nrow=J))

d2<-do.call(adiag,d2temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d2temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J),matrix(d11a2[[6]],ncol=J,nrow=J),matrix(d11a2[[7]],ncol=J,nrow=J),matrix(d11a2[[8]],ncol=J,nrow=J),matrix(d11a2[[9]],ncol=J,nrow=J),matrix(d11a2[[10]],ncol=J,nrow=J),matrix(d11a2[[11]],ncol=J,nrow=J))

d22<-do.call(adiag,d2temp2)


sig2<-vy*c+beta^2*eta*d2+beta^2*vx*d22
sig2i<-solve(sig2)
b2<-t(a2)%*%sig2i%*%a2
b2i<-solve(b2)
i2<-I
j2<-J
#############
# design 3 (I = 2 different values, replicated J = 11 times)
I<-f[length(f)-2]
J<-n/I
th<-c(rep(a,J), rep(a+(b-a)/(I-1),J))

a3<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],th[20],th[21],th[22]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J))

d3temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J))

d3<-do.call(adiag,d3temp)


d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d3temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J))

d32<-do.call(adiag,d3temp2)


sig3<-vy*c+beta^2*eta*d3+beta^2*vx*d32
sig3i<-solve(sig3)
b3<-t(a3)%*%sig3i%*%a3
b3i<-solve(b3)
i3<-I
j3<-J
#################

deter<-c(det(b1i),det(b2i),det(b3i))
which(deter==min(deter))

eig1<-eigen(b1i)
tr1<-sum(eig1$values)
eig2<-eigen(b2i)
tr2<-sum(eig2$values)
eig3<-eigen(b3i)
tr3<-sum(eig3$values)
 

 
trace<-c(tr1,tr2,tr3)
best<-which(trace==min(trace))


vb<-c((b1i[[2,2]])^0.5,(b2i[[2,2]])^0.5,(b3i[[2,2]])^0.5)
bestb<-which(vb==min(vb))

cq<-diag(1,nq)
thq<-c(rep(beta,nq))
th2q<-matrix(c(1,xnew),nrow=1)
a5q<-matrix(c(rep(beta,nq)), ncol=1)
varab1<-th2q%*%b1i%*%t(th2q)
d11aq1<-list(matrix(varab1, nrow=nq,ncol=nq))
dq1<-do.call(adiag,d11aq1)
sigq1<-dq1
sigq1i<-solve(sigq1)
bq1<-t(a5q)%*%sigq1i%*%a5q
stuq1<-(1/bq1+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b1i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq1<-stuq1/xnew*100

varab2<-th2q%*%b2i%*%t(th2q)
d11aq2<-list(matrix(varab2, nrow=nq,ncol=nq))
dq2<-do.call(adiag,d11aq2)
sigq2<-dq2
sigq2i<-solve(sigq2)
bq2<-t(a5q)%*%sigq2i%*%a5q
stuq2<-(1/bq2+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b2i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq2<-stuq2/xnew*100

varab3<-th2q%*%b3i%*%t(th2q)
d11aq3<-list(matrix(varab3, nrow=nq,ncol=nq))
dq3<-do.call(adiag,d11aq3)
sigq3<-dq3
sigq3i<-solve(sigq3)
bq3<-t(a5q)%*%sigq3i%*%a5q
stuq3<-(1/bq3+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b3i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq3<-stuq3/xnew*100
finalrelst<-c(relstuq1,relstuq2,relstuq3)
finbest<-which(finalrelst==min(finalrelst))

finali<-c(i1,i2,i3)
finalj<-c(j1,j2,j3)
jbest<-finalj[finbest]
ibest<-finali[finbest]
totrep<-n+ny
finstd<-finalrelst[finbest]
samplrep<-ny/r
result<-list(samplrep=samplrep,nosamples=r,finstd=finstd,ibest=ibest,jbest=jbest,totrep=totrep)

return(result)
}
############################################
############################################
caldesign124<-function(n,a,b,vy,vx,beta,eta,ny,vs,xnew,vi){
f<-divisors(n)
fy<-divisors(ny)
r<-1
if(vs>0 |vi>0){r=ny}
c<-diag(1,n)
v1<-(beta^2*eta*a^2)^0.5
v2<-(beta^2*eta*b^2)^0.5
nq<-1
#########
# design 1 (all I = 24 values are different, replicated J = 1 time)

J<-1
I<-f[length(f)]
th<-c(a,a+(b-a)/(I-1),a+2*(b-a)/(I-1),a+3*(b-a)/(I-1),a+4*(b-a)/(I-1),
                  a+5*(b-a)/(I-1),a+6*(b-a)/(I-1),a+7*(b-a)/(I-1),a+8*(b-a)/(I-1),a+9*(b-a)/(I-1),a+10*(b-a)/(I-1),a+11*(b-a)/(I-1),
a+12*(b-a)/(I-1),a+13*(b-a)/(I-1),a+14*(b-a)/(I-1),a+15*(b-a)/(I-1),
                  a+16*(b-a)/(I-1),a+17*(b-a)/(I-1),a+18*(b-a)/(I-1),a+19*(b-a)/(I-1),a+20*(b-a)/(I-1),a+21*(b-a)/(I-1),a+22*(b-a)/(I-1),b)
 
a1<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],th[20],th[21],th[22],th[23],th[24]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+2*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+4*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+5*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+6*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+7*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+8*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+9*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+10*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+11*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+12*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+13*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+14*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+15*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+16*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+17*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+18*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+19*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+20*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+21*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+22*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+23*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J))

d1temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J),matrix(d11a[[6]],ncol=J,nrow=J),matrix(d11a[[7]],ncol=J,nrow=J),matrix(d11a[[8]],ncol=J,nrow=J),matrix(d11a[[9]],ncol=J,nrow=J),
matrix(d11a[[10]],ncol=J,nrow=J),matrix(d11a[[11]],ncol=J,nrow=J),matrix(d11a[[12]],ncol=J,nrow=J),matrix(d11a[[13]],ncol=J,nrow=J),matrix(d11a[[14]],ncol=J,nrow=J),matrix(d11a[[15]],ncol=J,nrow=J),matrix(d11a[[16]],ncol=J,nrow=J),matrix(d11a[[17]],ncol=J,nrow=J),
matrix(d11a[[18]],ncol=J,nrow=J),matrix(d11a[[19]],ncol=J,nrow=J),matrix(d11a[[20]],ncol=J,nrow=J),matrix(d11a[[21]],ncol=J,nrow=J),
matrix(d11a[[22]],ncol=J,nrow=J),matrix(d11a[[23]],ncol=J,nrow=J),matrix(d11a[[24]],ncol=J,nrow=J))
d1<-do.call(adiag,d1temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d1temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J),matrix(d11a2[[6]],ncol=J,nrow=J),matrix(d11a2[[7]],ncol=J,nrow=J),matrix(d11a2[[8]],ncol=J,nrow=J),matrix(d11a2[[9]],ncol=J,nrow=J),
matrix(d11a2[[10]],ncol=J,nrow=J),matrix(d11a2[[11]],ncol=J,nrow=J),matrix(d11a2[[12]],ncol=J,nrow=J),matrix(d11a2[[13]],ncol=J,nrow=J),matrix(d11a2[[14]],ncol=J,nrow=J),matrix(d11a2[[15]],ncol=J,nrow=J),matrix(d11a2[[16]],ncol=J,nrow=J),matrix(d11a2[[17]],ncol=J,nrow=J),
matrix(d11a2[[18]],ncol=J,nrow=J),matrix(d11a2[[19]],ncol=J,nrow=J),matrix(d11a2[[20]],ncol=J,nrow=J),matrix(d11a2[[21]],ncol=J,nrow=J),
matrix(d11a2[[22]],ncol=J,nrow=J),matrix(d11a2[[23]],ncol=J,nrow=J),matrix(d11a2[[24]],ncol=J,nrow=J))
d12<-do.call(adiag,d1temp2)


sig1<-vy*c+beta^2*eta*d1+beta^2*vx*d12

sig1i<-solve(sig1)
b1<-t(a1)%*%sig1i%*%a1
b1i<-solve(b1)
i1<-I
j1<-J
#############
# design 2 (I=12 different values, replicated J = 2 times)
I <-f[length(f)-1]
J<-n/I
th<-c(rep(a,J),rep(a+(b-a)/(I-1),J),rep(a+2*(b-a)/(I-1),J),rep(a+3*(b-a)/(I-1),J),rep(a+4*(b-a)/(I-1),J),rep(a+5*(b-a)/(I-1),J),rep(a+6*(b-a)/(I-1),J),rep(a+7*(b-a)/(I-1),J),rep(a+8*(b-a)/(I-1),J),rep(a+9*(b-a)/(I-1),J),rep(a+10*(b-a)/(I-1),J),rep(a+11*(b-a)/(I-1),J))
 
a2<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],th[20],th[21],th[22],th[23],th[24]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+4*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+5*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+6*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+7*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+8*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+9*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+10*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+11*(b-a)/(I-1))^2,nrow=J,ncol=J))

d2temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J),matrix(d11a[[6]],ncol=J,nrow=J),matrix(d11a[[7]],ncol=J,nrow=J),matrix(d11a[[8]],ncol=J,nrow=J),matrix(d11a[[9]],ncol=J,nrow=J),matrix(d11a[[10]],ncol=J,nrow=J),matrix(d11a[[11]],ncol=J,nrow=J),matrix(d11a[[12]],ncol=J,nrow=J))

d2<-do.call(adiag,d2temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d2temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J),matrix(d11a2[[6]],ncol=J,nrow=J),matrix(d11a2[[7]],ncol=J,nrow=J),matrix(d11a2[[8]],ncol=J,nrow=J),matrix(d11a2[[9]],ncol=J,nrow=J),matrix(d11a2[[10]],ncol=J,nrow=J),matrix(d11a2[[11]],ncol=J,nrow=J),matrix(d11a2[[12]],ncol=J,nrow=J))

d22<-do.call(adiag,d2temp2)


sig2<-vy*c+beta^2*eta*d2+beta^2*vx*d22
sig2i<-solve(sig2)
b2<-t(a2)%*%sig2i%*%a2
b2i<-solve(b2)
i2<-I
j2<-J
#############
# design 3 (I = 8 different values, replicated J=3 times)
I<-f[length(f)-2]
J<-n/I
th<-c(rep(a,J), rep(a+(b-a)/(I-1),J), rep(a+2*(b-a)/(I-1),J),rep(a+3*(b-a)/(I-1),J),rep(a+4*(b-a)/(I-1),J), rep(a+5*(b-a)/(I-1),J),rep(a+6*(b-a)/(I-1),J),rep(a+7*(b-a)/(I-1),J))
a3<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],th[20],th[21],th[22],th[23],th[24]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+4*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+5*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+6*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+7*(b-a)/(I-1))^2,nrow=J,ncol=J))

d3temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J),matrix(d11a[[6]],ncol=J,nrow=J),matrix(d11a[[7]],ncol=J,nrow=J),matrix(d11a[[8]],ncol=J,nrow=J))

d3<-do.call(adiag,d3temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d3temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J),matrix(d11a2[[6]],ncol=J,nrow=J),matrix(d11a2[[7]],ncol=J,nrow=J),matrix(d11a2[[8]],ncol=J,nrow=J))

d32<-do.call(adiag,d3temp2)



sig3<-vy*c+beta^2*eta*d3+beta^2*vx*d32
sig3i<-solve(sig3)
b3<-t(a3)%*%sig3i%*%a3
b3i<-solve(b3)
i3<-I
j3<-J
#################
# design 4 (I = 6 different values, replicated J = 4 times)
I<-f[length(f)-3]
J<-n/I
th<-c(rep(a,J), rep(a+(b-a)/(I-1),J), rep(a+2*(b-a)/(I-1),J),rep(a+3*(b-a)/(I-1),J),rep(a+4*(b-a)/(I-1),J), rep(a+5*(b-a)/(I-1),J))
a4<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],th[20],th[21],th[22],th[23],th[24]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+4*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+5*(b-a)/(I-1))^2,nrow=J,ncol=J))

d4temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J),matrix(d11a[[6]],ncol=J,nrow=J))

d4<-do.call(adiag,d4temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d4temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J),matrix(d11a2[[6]],ncol=J,nrow=J))

d42<-do.call(adiag,d4temp2)


sig4<-vy*c+beta^2*eta*d4+beta^2*vx*d42
sig4i<-solve(sig4)
b4<-t(a4)%*%sig4i%*%a4
b4i<-solve(b4)
i4<-I
j4<-J
#################
# design 5 (I = 4 different values, replicated J = 6 times)
I<-f[length(f)-4]
J<-n/I
th<-c(rep(a,J), rep(a+(b-a)/(I-1),J), rep(a+2*(b-a)/(I-1),J),rep(a+3*(b-a)/(I-1),J))

a5<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],th[20],th[21],th[22],th[23],th[24]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(I-1))^2,nrow=J,ncol=J))

d5temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J))

d5<-do.call(adiag,d5temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d5temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J))

d52<-do.call(adiag,d5temp2)


sig5<-vy*c+beta^2*eta*d5+beta^2*vx*d52
sig5i<-solve(sig5)
b5<-t(a5)%*%sig5i%*%a5
b5i<-solve(b5)
i5<-I
j5<-J
##################
# design 6 (I = 3 different values, replicated J = 8 times)
I<-f[length(f)-5]
J<-n/I

th<-c(rep(a,J), rep(a+(b-a)/(I-1),J), rep(a+2*(b-a)/(I-1),J))

a6<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],th[20],th[21],th[22],th[23],th[24]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(I-1))^2,nrow=J,ncol=J))

d6temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J))

d6<-do.call(adiag,d6temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d6temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J))

d62<-do.call(adiag,d6temp2)


sig6<-vy*c+beta^2*eta*d6+beta^2*vx*d62
sig6i<-solve(sig6)
b6<-t(a6)%*%sig6i%*%a6
b6i<-solve(b6)
i6<-I
j6<-J

##################
# design 7 (I = 2 different values, replicated J = 12 times)
I<-f[length(f)-6]
J<-n/I
th<-c(rep(a,J), rep(a+(b-a)/(I-1),J))

a7<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],th[20],th[21],th[22],th[23],th[24]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J))

d7temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J))

d7<-do.call(adiag,d7temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d7temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J))

d72<-do.call(adiag,d7temp2)


sig7<-vy*c+beta^2*eta*d7+beta^2*vx*d72
sig7i<-solve(sig7)
b7<-t(a7)%*%sig7i%*%a7
b7i<-solve(b7)
i7<-I
j7<-J
#################

deter<-c(det(b1i),det(b2i),det(b3i),det(b4i),det(b5i),det(b6i),det(b7i))
which(deter==min(deter))

eig1<-eigen(b1i)
tr1<-sum(eig1$values)
eig2<-eigen(b2i)
tr2<-sum(eig2$values)
eig3<-eigen(b3i)
tr3<-sum(eig3$values)
eig4<-eigen(b4i)
tr4<-sum(eig4$values)
eig5<-eigen(b5i)
tr5<-sum(eig5$values)
eig6<-eigen(b6i)
tr6<-sum(eig6$values)
eig7<-eigen(b7i)
tr7<-sum(eig7$values)

 
trace<-c(tr1,tr2,tr3,tr4,tr5,tr6,tr7)
best<-which(trace==min(trace))


vb<-c((b1i[[2,2]])^0.5,(b2i[[2,2]])^0.5,(b3i[[2,2]])^0.5,(b4i[[2,2]])^0.5,(b5i[[2,2]])^0.5,(b6i[[2,2]])^0.5,(b7i[[2,2]])^0.5)
bestb<-which(vb==min(vb))

cq<-diag(1,nq)
thq<-c(rep(beta,nq))
th2q<-matrix(c(1,xnew),nrow=1)
a5q<-matrix(c(rep(beta,nq)), ncol=1)
varab1<-th2q%*%b1i%*%t(th2q)
d11aq1<-list(matrix(varab1, nrow=nq,ncol=nq))
dq1<-do.call(adiag,d11aq1)
sigq1<-dq1
sigq1i<-solve(sigq1)
bq1<-t(a5q)%*%sigq1i%*%a5q 
stuq1<-(1/bq1+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b1i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq1<-stuq1/xnew*100

varab2<-th2q%*%b2i%*%t(th2q)
d11aq2<-list(matrix(varab2, nrow=nq,ncol=nq))
dq2<-do.call(adiag,d11aq2)
sigq2<-dq2
sigq2i<-solve(sigq2)
bq2<-t(a5q)%*%sigq2i%*%a5q 
stuq2<-(1/bq2+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b2i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq2<-stuq2/xnew*100

varab3<-th2q%*%b3i%*%t(th2q)
d11aq3<-list(matrix(varab3, nrow=nq,ncol=nq))
dq3<-do.call(adiag,d11aq3)
sigq3<-dq3
sigq3i<-solve(sigq3)
bq3<-t(a5q)%*%sigq3i%*%a5q 
stuq3<-(1/bq3+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b3i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq3<-stuq3/xnew*100

varab4<-th2q%*%b4i%*%t(th2q)
d11aq4<-list(matrix(varab4, nrow=nq,ncol=nq))
dq4<-do.call(adiag,d11aq4)
sigq4<-dq4
sigq4i<-solve(sigq4)
bq4<-t(a5q)%*%sigq4i%*%a5q 
stuq4<-(1/bq4+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b4i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq4<-stuq4/xnew*100

varab5<-th2q%*%b5i%*%t(th2q)
d11aq5<-list(matrix(varab5, nrow=nq,ncol=nq))
dq5<-do.call(adiag,d11aq5)
sigq5<-dq5
sigq5i<-solve(sigq5)
bq5<-t(a5q)%*%sigq5i%*%a5q 
stuq5<-(1/bq5+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b5i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq5<-stuq5/xnew*100

varab6<-th2q%*%b6i%*%t(th2q)
d11aq6<-list(matrix(varab6, nrow=nq,ncol=nq))
dq6<-do.call(adiag,d11aq6)
sigq6<-dq6
sigq6i<-solve(sigq6)
bq6<-t(a5q)%*%sigq6i%*%a5q 
stuq6<-(1/bq6+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b6i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq6<-stuq6/xnew*100

varab7<-th2q%*%b7i%*%t(th2q)
d11aq7<-list(matrix(varab7, nrow=nq,ncol=nq))
dq7<-do.call(adiag,d11aq7)
sigq7<-dq7
sigq7i<-solve(sigq7)
bq7<-t(a5q)%*%sigq7i%*%a5q 
stuq7<-(1/bq7+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b7i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq7<-stuq7/xnew*100
finalrelst<-c(relstuq1,relstuq2,relstuq3,relstuq4,relstuq5,relstuq6,relstuq7)
finbest<-which(finalrelst==min(finalrelst))

finali<-c(i1,i2,i3,i4,i5,i6,i7)
finalj<-c(j1,j2,j3,j4,j5,j6,j7)
jbest<-finalj[finbest]
ibest<-finali[finbest]
totrep<-n+ny
finstd<-finalrelst[finbest]
samplrep<-ny/r
result<-list(samplrep=samplrep,nosamples=r,finstd=finstd,ibest=ibest,jbest=jbest,totrep=totrep)

return(result)
}

#########################################################################################

############################################
############################################
caldesign128<-function(n,a,b,vy,vx,beta,eta,ny,vs,xnew,vi){
f<-divisors(n)
fy<-divisors(ny)
r<-1
if(vs>0 |vi>0){r=ny}
c<-diag(1,n)
v1<-(beta^2*eta*a^2)^0.5
v2<-(beta^2*eta*b^2)^0.5
nq<-1
#########
# design 1 (all I = 28 values are different, replicated J = 1 time)

J<-1
I<-f[length(f)]
th<-c(a,a+(b-a)/(I-1),a+2*(b-a)/(I-1),a+3*(b-a)/(I-1),a+4*(b-a)/(I-1),
                  a+5*(b-a)/(I-1),a+6*(b-a)/(I-1),a+7*(b-a)/(I-1),a+8*(b-a)/(I-1),a+9*(b-a)/(I-1),a+10*(b-a)/(I-1),a+11*(b-a)/(I-1),
a+12*(b-a)/(I-1),a+13*(b-a)/(I-1),a+14*(b-a)/(I-1),a+15*(b-a)/(I-1),
                  a+16*(b-a)/(I-1),a+17*(b-a)/(I-1),a+18*(b-a)/(I-1),a+19*(b-a)/(I-1),a+20*(b-a)/(I-1),a+21*(b-a)/(I-1),a+22*(b-a)/(I-1),a+23*(b-a)/(I-1),a+24*(b-a)/(I-1),a+25*(b-a)/(I-1),a+26*(b-a)/(I-1),b)
 
a1<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],th[20],th[21],th[22],th[23],th[24],th[25],th[26],th[27],th[28]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+2*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+4*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+5*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+6*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+7*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+8*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+9*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+10*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+11*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+12*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+13*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+14*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+15*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+16*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+17*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+18*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+19*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+20*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+21*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+22*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+23*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+24*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+25*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+26*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+27*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J))

d1temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J),matrix(d11a[[6]],ncol=J,nrow=J),matrix(d11a[[7]],ncol=J,nrow=J),matrix(d11a[[8]],ncol=J,nrow=J),matrix(d11a[[9]],ncol=J,nrow=J),
matrix(d11a[[10]],ncol=J,nrow=J),matrix(d11a[[11]],ncol=J,nrow=J),matrix(d11a[[12]],ncol=J,nrow=J),matrix(d11a[[13]],ncol=J,nrow=J),matrix(d11a[[14]],ncol=J,nrow=J),matrix(d11a[[15]],ncol=J,nrow=J),matrix(d11a[[16]],ncol=J,nrow=J),matrix(d11a[[17]],ncol=J,nrow=J),
matrix(d11a[[18]],ncol=J,nrow=J),matrix(d11a[[19]],ncol=J,nrow=J),matrix(d11a[[20]],ncol=J,nrow=J),matrix(d11a[[21]],ncol=J,nrow=J),
matrix(d11a[[22]],ncol=J,nrow=J),matrix(d11a[[23]],ncol=J,nrow=J),matrix(d11a[[24]],ncol=J,nrow=J),matrix(d11a[[25]],ncol=J,nrow=J),
matrix(d11a[[26]],ncol=J,nrow=J),matrix(d11a[[27]],ncol=J,nrow=J),matrix(d11a[[28]],ncol=J,nrow=J))
d1<-do.call(adiag,d1temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d1temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J),matrix(d11a2[[6]],ncol=J,nrow=J),matrix(d11a2[[7]],ncol=J,nrow=J),matrix(d11a2[[8]],ncol=J,nrow=J),matrix(d11a2[[9]],ncol=J,nrow=J),
matrix(d11a2[[10]],ncol=J,nrow=J),matrix(d11a2[[11]],ncol=J,nrow=J),matrix(d11a2[[12]],ncol=J,nrow=J),matrix(d11a2[[13]],ncol=J,nrow=J),matrix(d11a2[[14]],ncol=J,nrow=J),matrix(d11a2[[15]],ncol=J,nrow=J),matrix(d11a2[[16]],ncol=J,nrow=J),matrix(d11a2[[17]],ncol=J,nrow=J),
matrix(d11a2[[18]],ncol=J,nrow=J),matrix(d11a2[[19]],ncol=J,nrow=J),matrix(d11a2[[20]],ncol=J,nrow=J),matrix(d11a2[[21]],ncol=J,nrow=J),
matrix(d11a2[[22]],ncol=J,nrow=J),matrix(d11a2[[23]],ncol=J,nrow=J),matrix(d11a2[[24]],ncol=J,nrow=J),matrix(d11a2[[25]],ncol=J,nrow=J),
matrix(d11a2[[26]],ncol=J,nrow=J),matrix(d11a2[[27]],ncol=J,nrow=J),matrix(d11a2[[28]],ncol=J,nrow=J))
d12<-do.call(adiag,d1temp2)


sig1<-vy*c+beta^2*eta*d1+beta^2*vx*d12

sig1i<-solve(sig1)
b1<-t(a1)%*%sig1i%*%a1
b1i<-solve(b1)
i1<-I
j1<-J
#############
# design 2 (I=14 different values, replicated J = 2 times)
I <-f[length(f)-1]
J<-n/I
th<-c(rep(a,J),rep(a+(b-a)/(I-1),J),rep(a+2*(b-a)/(I-1),J),rep(a+3*(b-a)/(I-1),J),rep(a+4*(b-a)/(I-1),J),rep(a+5*(b-a)/(I-1),J),rep(a+6*(b-a)/(I-1),J),rep(a+7*(b-a)/(I-1),J),rep(a+8*(b-a)/(I-1),J),rep(a+9*(b-a)/(I-1),J),rep(a+10*(b-a)/(I-1),J),rep(a+11*(b-a)/(I-1),J),rep(a+12*(b-a)/(I-1),J),rep(a+13*(b-a)/(I-1),J))
 
a2<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],th[20],th[21],th[22],th[23],th[24],th[25],th[26],th[27],th[28]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+4*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+5*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+6*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+7*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+8*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+9*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+10*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+11*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+12*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+13*(b-a)/(I-1))^2,nrow=J,ncol=J))

d2temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J),matrix(d11a[[6]],ncol=J,nrow=J),matrix(d11a[[7]],ncol=J,nrow=J),matrix(d11a[[8]],ncol=J,nrow=J),matrix(d11a[[9]],ncol=J,nrow=J),matrix(d11a[[10]],ncol=J,nrow=J),matrix(d11a[[11]],ncol=J,nrow=J),matrix(d11a[[12]],ncol=J,nrow=J),matrix(d11a[[13]],ncol=J,nrow=J),matrix(d11a[[14]],ncol=J,nrow=J))

d2<-do.call(adiag,d2temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d2temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J),matrix(d11a2[[6]],ncol=J,nrow=J),matrix(d11a2[[7]],ncol=J,nrow=J),matrix(d11a2[[8]],ncol=J,nrow=J),matrix(d11a2[[9]],ncol=J,nrow=J),matrix(d11a2[[10]],ncol=J,nrow=J),matrix(d11a2[[11]],ncol=J,nrow=J),matrix(d11a2[[12]],ncol=J,nrow=J),matrix(d11a2[[13]],ncol=J,nrow=J),matrix(d11a2[[14]],ncol=J,nrow=J))

d22<-do.call(adiag,d2temp2)


sig2<-vy*c+beta^2*eta*d2+beta^2*vx*d22
sig2i<-solve(sig2)
b2<-t(a2)%*%sig2i%*%a2
b2i<-solve(b2)
i2<-I
j2<-J
#############
# design 3 (I = 7 different values, replicated J=4 times)
I<-f[length(f)-2]
J<-n/I
th<-c(rep(a,J), rep(a+(b-a)/(I-1),J), rep(a+2*(b-a)/(I-1),J),rep(a+3*(b-a)/(I-1),J),rep(a+4*(b-a)/(I-1),J), rep(a+5*(b-a)/(I-1),J),rep(a+6*(b-a)/(I-1),J))
a3<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],th[20],th[21],th[22],th[23],th[24],th[25],th[26],th[27],th[28]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+4*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+5*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+6*(b-a)/(I-1))^2,nrow=J,ncol=J))

d3temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J),matrix(d11a[[6]],ncol=J,nrow=J),matrix(d11a[[7]],ncol=J,nrow=J))

d3<-do.call(adiag,d3temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d3temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J),matrix(d11a2[[6]],ncol=J,nrow=J),matrix(d11a2[[7]],ncol=J,nrow=J))

d32<-do.call(adiag,d3temp2)



sig3<-vy*c+beta^2*eta*d3+beta^2*vx*d32
sig3i<-solve(sig3)
b3<-t(a3)%*%sig3i%*%a3
b3i<-solve(b3)
i3<-I
j3<-J
#################
# design 4 (I = 4 different values, replicated J = 7 times)
I<-f[length(f)-3]
J<-n/I
th<-c(rep(a,J), rep(a+(b-a)/(I-1),J), rep(a+2*(b-a)/(I-1),J),rep(a+3*(b-a)/(I-1),J))
a4<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],th[20],th[21],th[22],th[23],th[24],th[25],th[26],th[27],th[28]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(I-1))^2,nrow=J,ncol=J))

d4temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J))

d4<-do.call(adiag,d4temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))


d4temp2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))


d42<-do.call(adiag,d4temp2)


sig4<-vy*c+beta^2*eta*d4+beta^2*vx*d42
sig4i<-solve(sig4)
b4<-t(a4)%*%sig4i%*%a4
b4i<-solve(b4)
i4<-I
j4<-J
#################
# design 5 (I = 2 different values, replicated J = 14 times)
I<-f[length(f)-4]
J<-n/I
th<-c(rep(a,J), rep(a+(b-a)/(I-1),J))

a5<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],th[20],th[21],th[22],th[23],th[24],th[25],th[26],th[27],th[28]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J))

d5temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J))

d5<-do.call(adiag,d5temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d5temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J))

d52<-do.call(adiag,d5temp2)


sig5<-vy*c+beta^2*eta*d5+beta^2*vx*d52
sig5i<-solve(sig5)
b5<-t(a5)%*%sig5i%*%a5
b5i<-solve(b5)
i5<-I
j5<-J
##################


deter<-c(det(b1i),det(b2i),det(b3i),det(b4i),det(b5i))
which(deter==min(deter))

eig1<-eigen(b1i)
tr1<-sum(eig1$values)
eig2<-eigen(b2i)
tr2<-sum(eig2$values)
eig3<-eigen(b3i)
tr3<-sum(eig3$values)
eig4<-eigen(b4i)
tr4<-sum(eig4$values)
eig5<-eigen(b5i)
tr5<-sum(eig5$values)

 
trace<-c(tr1,tr2,tr3,tr4,tr5)
best<-which(trace==min(trace))


vb<-c((b1i[[2,2]])^0.5,(b2i[[2,2]])^0.5,(b3i[[2,2]])^0.5,(b4i[[2,2]])^0.5,(b5i[[2,2]])^0.5)
bestb<-which(vb==min(vb))

cq<-diag(1,nq)
thq<-c(rep(beta,nq))
th2q<-matrix(c(1,xnew),nrow=1)
a5q<-matrix(c(rep(beta,nq)), ncol=1)
varab1<-th2q%*%b1i%*%t(th2q)
d11aq1<-list(matrix(varab1, nrow=nq,ncol=nq))
dq1<-do.call(adiag,d11aq1)
sigq1<-dq1
sigq1i<-solve(sigq1)
bq1<-t(a5q)%*%sigq1i%*%a5q 
stuq1<-(1/bq1+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b1i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq1<-stuq1/xnew*100

varab2<-th2q%*%b2i%*%t(th2q)
d11aq2<-list(matrix(varab2, nrow=nq,ncol=nq))
dq2<-do.call(adiag,d11aq2)
sigq2<-dq2
sigq2i<-solve(sigq2)
bq2<-t(a5q)%*%sigq2i%*%a5q 
stuq2<-(1/bq2+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b2i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq2<-stuq2/xnew*100

varab3<-th2q%*%b3i%*%t(th2q)
d11aq3<-list(matrix(varab3, nrow=nq,ncol=nq))
dq3<-do.call(adiag,d11aq3)
sigq3<-dq3
sigq3i<-solve(sigq3)
bq3<-t(a5q)%*%sigq3i%*%a5q 
stuq3<-(1/bq3+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b3i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq3<-stuq3/xnew*100

varab4<-th2q%*%b4i%*%t(th2q)
d11aq4<-list(matrix(varab4, nrow=nq,ncol=nq))
dq4<-do.call(adiag,d11aq4)
sigq4<-dq4
sigq4i<-solve(sigq4)
bq4<-t(a5q)%*%sigq4i%*%a5q 
stuq4<-(1/bq4+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b4i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq4<-stuq4/xnew*100

varab5<-th2q%*%b5i%*%t(th2q)
d11aq5<-list(matrix(varab5, nrow=nq,ncol=nq))
dq5<-do.call(adiag,d11aq5)
sigq5<-dq5
sigq5i<-solve(sigq5)
bq5<-t(a5q)%*%sigq5i%*%a5q 
stuq5<-(1/bq5+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b5i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq5<-stuq5/xnew*100


finalrelst<-c(relstuq1,relstuq2,relstuq3,relstuq4,relstuq5)
finbest<-which(finalrelst==min(finalrelst))

finali<-c(i1,i2,i3,i4,i5)
finalj<-c(j1,j2,j3,j4,j5)
jbest<-finalj[finbest]
ibest<-finali[finbest]
totrep<-n+ny
finstd<-finalrelst[finbest]
samplrep<-ny/r
result<-list(samplrep=samplrep,nosamples=r,finstd=finstd,ibest=ibest,jbest=jbest,totrep=totrep)

return(result)
}



############################################
caldesign130<-function(n,a,b,vy,vx,beta,eta,ny,vs,xnew,vi){
f<-divisors(n)
fy<-divisors(ny)
r<-1
if(vs>0 |vi>0){r=ny}
c<-diag(1,n)
v1<-(beta^2*eta*a^2)^0.5
v2<-(beta^2*eta*b^2)^0.5
nq<-1
#########
# design 1 (all I = 30 values are different, replicated J = 1 time)

J<-1
I<-f[length(f)]
th<-c(a,a+(b-a)/(I-1),a+2*(b-a)/(I-1),a+3*(b-a)/(I-1),a+4*(b-a)/(I-1),
                  a+5*(b-a)/(I-1),a+6*(b-a)/(I-1),a+7*(b-a)/(I-1),a+8*(b-a)/(I-1),a+9*(b-a)/(I-1),a+10*(b-a)/(I-1),a+11*(b-a)/(I-1),
a+12*(b-a)/(I-1),a+13*(b-a)/(I-1),a+14*(b-a)/(I-1),a+15*(b-a)/(I-1),
                  a+16*(b-a)/(I-1),a+17*(b-a)/(I-1),a+18*(b-a)/(I-1),a+19*(b-a)/(I-1),a+20*(b-a)/(I-1),a+21*(b-a)/(I-1),a+22*(b-a)/(I-1),a+23*(b-a)/(I-1),a+24*(b-a)/(I-1),a+25*(b-a)/(I-1),a+26*(b-a)/(I-1),a+27*(b-a)/(I-1),a+28*(b-a)/(I-1),b)
 
a1<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],th[20],th[21],th[22],th[23],th[24],th[25],th[26],th[27],th[28],th[29],th[30]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+2*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+4*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+5*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+6*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+7*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+8*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+9*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+10*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+11*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+12*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+13*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+14*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+15*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+16*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+17*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+18*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+19*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+20*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+21*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+22*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+23*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+24*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+25*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+26*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+27*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+28*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+29*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J))

d1temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J),matrix(d11a[[6]],ncol=J,nrow=J),matrix(d11a[[7]],ncol=J,nrow=J),matrix(d11a[[8]],ncol=J,nrow=J),matrix(d11a[[9]],ncol=J,nrow=J),
matrix(d11a[[10]],ncol=J,nrow=J),matrix(d11a[[11]],ncol=J,nrow=J),matrix(d11a[[12]],ncol=J,nrow=J),matrix(d11a[[13]],ncol=J,nrow=J),matrix(d11a[[14]],ncol=J,nrow=J),matrix(d11a[[15]],ncol=J,nrow=J),matrix(d11a[[16]],ncol=J,nrow=J),matrix(d11a[[17]],ncol=J,nrow=J),
matrix(d11a[[18]],ncol=J,nrow=J),matrix(d11a[[19]],ncol=J,nrow=J),matrix(d11a[[20]],ncol=J,nrow=J),matrix(d11a[[21]],ncol=J,nrow=J),
matrix(d11a[[22]],ncol=J,nrow=J),matrix(d11a[[23]],ncol=J,nrow=J),matrix(d11a[[24]],ncol=J,nrow=J),matrix(d11a[[25]],ncol=J,nrow=J),matrix(d11a[[26]],ncol=J,nrow=J),matrix(d11a[[27]],ncol=J,nrow=J),
matrix(d11a[[28]],ncol=J,nrow=J),matrix(d11a[[29]],ncol=J,nrow=J),matrix(d11a[[30]],ncol=J,nrow=J))
d1<-do.call(adiag,d1temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d1temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J),matrix(d11a2[[6]],ncol=J,nrow=J),matrix(d11a2[[7]],ncol=J,nrow=J),matrix(d11a2[[8]],ncol=J,nrow=J),matrix(d11a2[[9]],ncol=J,nrow=J),
matrix(d11a2[[10]],ncol=J,nrow=J),matrix(d11a2[[11]],ncol=J,nrow=J),matrix(d11a2[[12]],ncol=J,nrow=J),matrix(d11a2[[13]],ncol=J,nrow=J),matrix(d11a2[[14]],ncol=J,nrow=J),matrix(d11a2[[15]],ncol=J,nrow=J),matrix(d11a2[[16]],ncol=J,nrow=J),matrix(d11a2[[17]],ncol=J,nrow=J),
matrix(d11a2[[18]],ncol=J,nrow=J),matrix(d11a2[[19]],ncol=J,nrow=J),matrix(d11a2[[20]],ncol=J,nrow=J),matrix(d11a2[[21]],ncol=J,nrow=J),
matrix(d11a2[[22]],ncol=J,nrow=J),matrix(d11a2[[23]],ncol=J,nrow=J),matrix(d11a2[[24]],ncol=J,nrow=J),matrix(d11a2[[25]],ncol=J,nrow=J),matrix(d11a2[[26]],ncol=J,nrow=J),matrix(d11a2[[27]],ncol=J,nrow=J),
matrix(d11a2[[28]],ncol=J,nrow=J),matrix(d11a2[[29]],ncol=J,nrow=J),matrix(d11a2[[30]],ncol=J,nrow=J))
d12<-do.call(adiag,d1temp2)


sig1<-vy*c+beta^2*eta*d1+beta^2*vx*d12

sig1i<-solve(sig1)
b1<-t(a1)%*%sig1i%*%a1
b1i<-solve(b1)
i1<-I
j1<-J
#############
# design 2 (I=15 different values, replicated J = 2 times)
I <-f[length(f)-1]
J<-n/I
th<-c(rep(a,J),rep(a+(b-a)/(I-1),J),rep(a+2*(b-a)/(I-1),J),rep(a+3*(b-a)/(I-1),J),rep(a+4*(b-a)/(I-1),J),rep(a+5*(b-a)/(I-1),J),rep(a+6*(b-a)/(I-1),J),rep(a+7*(b-a)/(I-1),J),rep(a+8*(b-a)/(I-1),J),rep(a+9*(b-a)/(I-1),J),rep(a+10*(b-a)/(I-1),J),rep(a+11*(b-a)/(I-1),J),rep(a+12*(b-a)/(I-1),J),rep(a+13*(b-a)/(I-1),J),rep(a+14*(b-a)/(I-1),J))
 
a2<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],th[20],th[21],th[22],th[23],th[24],th[25],th[26],th[27],th[28],th[29],th[30]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+4*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+5*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+6*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+7*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+8*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+9*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+10*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+11*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+12*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+13*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+14*(b-a)/(I-1))^2,nrow=J,ncol=J))

d2temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J),matrix(d11a[[6]],ncol=J,nrow=J),matrix(d11a[[7]],ncol=J,nrow=J),matrix(d11a[[8]],ncol=J,nrow=J),matrix(d11a[[9]],ncol=J,nrow=J),matrix(d11a[[10]],ncol=J,nrow=J),matrix(d11a[[11]],ncol=J,nrow=J),matrix(d11a[[12]],ncol=J,nrow=J),matrix(d11a[[13]],ncol=J,nrow=J),matrix(d11a[[14]],ncol=J,nrow=J),matrix(d11a[[15]],ncol=J,nrow=J))

d2<-do.call(adiag,d2temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d2temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J),matrix(d11a2[[6]],ncol=J,nrow=J),matrix(d11a2[[7]],ncol=J,nrow=J),matrix(d11a2[[8]],ncol=J,nrow=J),matrix(d11a2[[9]],ncol=J,nrow=J),matrix(d11a2[[10]],ncol=J,nrow=J),matrix(d11a2[[11]],ncol=J,nrow=J),matrix(d11a2[[12]],ncol=J,nrow=J),matrix(d11a2[[13]],ncol=J,nrow=J),matrix(d11a2[[14]],ncol=J,nrow=J),matrix(d11a2[[15]],ncol=J,nrow=J))

d22<-do.call(adiag,d2temp2)


sig2<-vy*c+beta^2*eta*d2+beta^2*vx*d22
sig2i<-solve(sig2)
b2<-t(a2)%*%sig2i%*%a2
b2i<-solve(b2)
i2<-I
j2<-J
#############
# design 3 (I = 10 different values, replicated J=3 times)
I<-f[length(f)-2]
J<-n/I
th<-c(rep(a,J), rep(a+(b-a)/(I-1),J), rep(a+2*(b-a)/(I-1),J),rep(a+3*(b-a)/(I-1),J),rep(a+4*(b-a)/(I-1),J), rep(a+5*(b-a)/(I-1),J),rep(a+6*(b-a)/(I-1),J),rep(a+7*(b-a)/(I-1),J),rep(a+8*(b-a)/(I-1),J),rep(a+9*(b-a)/(I-1),J))
a3<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],th[20],th[21],th[22],th[23],th[24],th[25],th[26],th[27],th[28],th[29],th[30]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+4*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+5*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+6*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+7*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+8*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+9*(b-a)/(I-1))^2,nrow=J,ncol=J))

d3temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J),matrix(d11a[[6]],ncol=J,nrow=J),matrix(d11a[[7]],ncol=J,nrow=J),matrix(d11a[[8]],ncol=J,nrow=J),matrix(d11a[[9]],ncol=J,nrow=J),matrix(d11a[[10]],ncol=J,nrow=J))

d3<-do.call(adiag,d3temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d3temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J),matrix(d11a2[[6]],ncol=J,nrow=J),matrix(d11a2[[7]],ncol=J,nrow=J),matrix(d11a2[[8]],ncol=J,nrow=J),matrix(d11a2[[9]],ncol=J,nrow=J),matrix(d11a2[[10]],ncol=J,nrow=J))

d32<-do.call(adiag,d3temp2)



sig3<-vy*c+beta^2*eta*d3+beta^2*vx*d32
sig3i<-solve(sig3)
b3<-t(a3)%*%sig3i%*%a3
b3i<-solve(b3)
i3<-I
j3<-J
#################
# design 4 (I = 6 different values, replicated J = 5 times)
I<-f[length(f)-3]
J<-n/I
th<-c(rep(a,J), rep(a+(b-a)/(I-1),J), rep(a+2*(b-a)/(I-1),J),rep(a+3*(b-a)/(I-1),J),rep(a+4*(b-a)/(I-1),J), rep(a+5*(b-a)/(I-1),J))
a4<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],th[20],th[21],th[22],th[23],th[24],th[25],th[26],th[27],th[28],th[29],th[30]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+4*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+5*(b-a)/(I-1))^2,nrow=J,ncol=J))

d4temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J),matrix(d11a[[6]],ncol=J,nrow=J))

d4<-do.call(adiag,d4temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d4temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J),matrix(d11a2[[6]],ncol=J,nrow=J))

d42<-do.call(adiag,d4temp2)


sig4<-vy*c+beta^2*eta*d4+beta^2*vx*d42
sig4i<-solve(sig4)
b4<-t(a4)%*%sig4i%*%a4
b4i<-solve(b4)
i4<-I
j4<-J
#################
# design 5 (I = 5 different values, replicated J = 6 times)
I<-f[length(f)-4]
J<-n/I
th<-c(rep(a,J), rep(a+(b-a)/(I-1),J), rep(a+2*(b-a)/(I-1),J),rep(a+3*(b-a)/(I-1),J),rep(a+4*(b-a)/(I-1),J))

a5<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],th[20],th[21],th[22],th[23],th[24],th[25],th[26],th[27],th[28],th[29],th[30]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+4*(b-a)/(I-1))^2,nrow=J,ncol=J))

d5temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),
matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J))

d5<-do.call(adiag,d5temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d5temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),
matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J))

d52<-do.call(adiag,d5temp2)


sig5<-vy*c+beta^2*eta*d5+beta^2*vx*d52
sig5i<-solve(sig5)
b5<-t(a5)%*%sig5i%*%a5
b5i<-solve(b5)
i5<-I
j5<-J
##################
# design 6 (I = 3 different values, replicated J = 10 times)
I<-f[length(f)-5]
J<-n/I

th<-c(rep(a,J), rep(a+(b-a)/(I-1),J), rep(a+2*(b-a)/(I-1),J))

a6<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],th[20],th[21],th[22],th[23],th[24],th[25],th[26],th[27],th[28],th[29],th[30]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(I-1))^2,nrow=J,ncol=J))

d6temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J))

d6<-do.call(adiag,d6temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d6temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J))

d62<-do.call(adiag,d6temp2)


sig6<-vy*c+beta^2*eta*d6+beta^2*vx*d62
sig6i<-solve(sig6)
b6<-t(a6)%*%sig6i%*%a6
b6i<-solve(b6)
i6<-I
j6<-J

##################
# design 7 (I = 2 different values, replicated J = 15 times)
I<-f[length(f)-6]
J<-n/I
th<-c(rep(a,J), rep(a+(b-a)/(I-1),J))

a7<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],th[20],th[21],th[22],th[23],th[24],th[25],th[26],th[27],th[28],th[29],th[30]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J))

d7temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J))

d7<-do.call(adiag,d7temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d7temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J))

d72<-do.call(adiag,d7temp2)


sig7<-vy*c+beta^2*eta*d7+beta^2*vx*d72
sig7i<-solve(sig7)
b7<-t(a7)%*%sig7i%*%a7
b7i<-solve(b7)
i7<-I
j7<-J
#################

deter<-c(det(b1i),det(b2i),det(b3i),det(b4i),det(b5i),det(b6i),det(b7i))
which(deter==min(deter))

eig1<-eigen(b1i)
tr1<-sum(eig1$values)
eig2<-eigen(b2i)
tr2<-sum(eig2$values)
eig3<-eigen(b3i)
tr3<-sum(eig3$values)
eig4<-eigen(b4i)
tr4<-sum(eig4$values)
eig5<-eigen(b5i)
tr5<-sum(eig5$values)
eig6<-eigen(b6i)
tr6<-sum(eig6$values)
eig7<-eigen(b7i)
tr7<-sum(eig7$values)

 
trace<-c(tr1,tr2,tr3,tr4,tr5,tr6,tr7)
best<-which(trace==min(trace))


vb<-c((b1i[[2,2]])^0.5,(b2i[[2,2]])^0.5,(b3i[[2,2]])^0.5,(b4i[[2,2]])^0.5,(b5i[[2,2]])^0.5,(b6i[[2,2]])^0.5,(b7i[[2,2]])^0.5)
bestb<-which(vb==min(vb))

cq<-diag(1,nq)
thq<-c(rep(beta,nq))
th2q<-matrix(c(1,xnew),nrow=1)
a5q<-matrix(c(rep(beta,nq)), ncol=1)
varab1<-th2q%*%b1i%*%t(th2q)
d11aq1<-list(matrix(varab1, nrow=nq,ncol=nq))
dq1<-do.call(adiag,d11aq1)
sigq1<-dq1
sigq1i<-solve(sigq1)
bq1<-t(a5q)%*%sigq1i%*%a5q 
stuq1<-(1/bq1+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b1i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq1<-stuq1/xnew*100

varab2<-th2q%*%b2i%*%t(th2q)
d11aq2<-list(matrix(varab2, nrow=nq,ncol=nq))
dq2<-do.call(adiag,d11aq2)
sigq2<-dq2
sigq2i<-solve(sigq2)
bq2<-t(a5q)%*%sigq2i%*%a5q 
stuq2<-(1/bq2+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b2i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq2<-stuq2/xnew*100

varab3<-th2q%*%b3i%*%t(th2q)
d11aq3<-list(matrix(varab3, nrow=nq,ncol=nq))
dq3<-do.call(adiag,d11aq3)
sigq3<-dq3
sigq3i<-solve(sigq3)
bq3<-t(a5q)%*%sigq3i%*%a5q 
stuq3<-(1/bq3+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b3i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq3<-stuq3/xnew*100

varab4<-th2q%*%b4i%*%t(th2q)
d11aq4<-list(matrix(varab4, nrow=nq,ncol=nq))
dq4<-do.call(adiag,d11aq4)
sigq4<-dq4
sigq4i<-solve(sigq4)
bq4<-t(a5q)%*%sigq4i%*%a5q 
stuq4<-(1/bq4+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b4i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq4<-stuq4/xnew*100

varab5<-th2q%*%b5i%*%t(th2q)
d11aq5<-list(matrix(varab5, nrow=nq,ncol=nq))
dq5<-do.call(adiag,d11aq5)
sigq5<-dq5
sigq5i<-solve(sigq5)
bq5<-t(a5q)%*%sigq5i%*%a5q 
stuq5<-(1/bq5+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b5i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq5<-stuq5/xnew*100

varab6<-th2q%*%b6i%*%t(th2q)
d11aq6<-list(matrix(varab6, nrow=nq,ncol=nq))
dq6<-do.call(adiag,d11aq6)
sigq6<-dq6
sigq6i<-solve(sigq6)
bq6<-t(a5q)%*%sigq6i%*%a5q 
stuq6<-(1/bq6+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b6i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq6<-stuq6/xnew*100

varab7<-th2q%*%b7i%*%t(th2q)
d11aq7<-list(matrix(varab7, nrow=nq,ncol=nq))
dq7<-do.call(adiag,d11aq7)
sigq7<-dq7
sigq7i<-solve(sigq7)
bq7<-t(a5q)%*%sigq7i%*%a5q 
stuq7<-(1/bq7+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b7i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq7<-stuq7/xnew*100
finalrelst<-c(relstuq1,relstuq2,relstuq3,relstuq4,relstuq5,relstuq6,relstuq7)
finbest<-which(finalrelst==min(finalrelst))

finali<-c(i1,i2,i3,i4,i5,i6,i7)
finalj<-c(j1,j2,j3,j4,j5,j6,j7)
jbest<-finalj[finbest]
ibest<-finali[finbest]
totrep<-n+ny
finstd<-finalrelst[finbest]
samplrep<-ny/r
result<-list(samplrep=samplrep,nosamples=r,finstd=finstd,ibest=ibest,jbest=jbest,totrep=totrep)

return(result)
}

#########################################################################################
############################################
caldesign136<-function(n,a,b,vy,vx,beta,eta,ny,vs,xnew,vi){
f<-divisors(n)
fy<-divisors(ny)
r<-1
if(vs>0 |vi>0){r=ny}
c<-diag(1,n)
v1<-(beta^2*eta*a^2)^0.5
v2<-(beta^2*eta*b^2)^0.5
nq<-1
#########
# #########
# design 1 (all I = 36 values are different, replicated J = 1 time)

J<-1
I<-f[length(f)]
th<-c(a,a+(b-a)/(I-1),a+2*(b-a)/(I-1),a+3*(b-a)/(I-1),a+4*(b-a)/(I-1),
                  a+5*(b-a)/(I-1),a+6*(b-a)/(I-1),a+7*(b-a)/(I-1),a+8*(b-a)/(I-1),a+9*(b-a)/(I-1),a+10*(b-a)/(I-1),a+11*(b-a)/(I-1),
a+12*(b-a)/(I-1),a+13*(b-a)/(I-1),a+14*(b-a)/(I-1),a+15*(b-a)/(I-1),
                  a+16*(b-a)/(I-1),a+17*(b-a)/(I-1),a+18*(b-a)/(I-1),a+19*(b-a)/(I-1),a+20*(b-a)/(I-1),a+21*(b-a)/(I-1),a+22*(b-a)/(I-1),
a+23*(b-a)/(I-1),a+24*(b-a)/(I-1),a+25*(b-a)/(I-1),a+26*(b-a)/(I-1),a+27*(b-a)/(I-1),a+28*(b-a)/(I-1),a+29*(b-a)/(I-1),a+30*(b-a)/(I-1),a+31*(b-a)/(I-1),a+32*(b-a)/(I-1),a+33*(b-a)/(I-1),a+34*(b-a)/(I-1),b)
 
a1<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],
th[20],th[21],th[22],th[23],th[24],th[25],th[26],th[27],th[28],th[29],th[30],th[31],th[32],th[33],th[34],th[35],th[36]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+2*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+4*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+5*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+6*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+7*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+8*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+9*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+10*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+11*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+12*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+13*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+14*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+15*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+16*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+17*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+18*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+19*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+20*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+21*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+22*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+23*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+24*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+25*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+26*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+27*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+28*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+29*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+30*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+31*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+32*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+33*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+34*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+35*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J))

d1temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J),matrix(d11a[[6]],ncol=J,nrow=J),matrix(d11a[[7]],ncol=J,nrow=J),matrix(d11a[[8]],ncol=J,nrow=J),matrix(d11a[[9]],ncol=J,nrow=J),
matrix(d11a[[10]],ncol=J,nrow=J),matrix(d11a[[11]],ncol=J,nrow=J),matrix(d11a[[12]],ncol=J,nrow=J),matrix(d11a[[13]],ncol=J,nrow=J),matrix(d11a[[14]],ncol=J,nrow=J),matrix(d11a[[15]],ncol=J,nrow=J),matrix(d11a[[16]],ncol=J,nrow=J),matrix(d11a[[17]],ncol=J,nrow=J),
matrix(d11a[[18]],ncol=J,nrow=J),matrix(d11a[[19]],ncol=J,nrow=J),matrix(d11a[[20]],ncol=J,nrow=J),matrix(d11a[[21]],ncol=J,nrow=J),
matrix(d11a[[22]],ncol=J,nrow=J),matrix(d11a[[23]],ncol=J,nrow=J),matrix(d11a[[24]],ncol=J,nrow=J),matrix(d11a[[25]],ncol=J,nrow=J),matrix(d11a[[26]],ncol=J,nrow=J),matrix(d11a[[27]],ncol=J,nrow=J),
matrix(d11a[[28]],ncol=J,nrow=J),matrix(d11a[[29]],ncol=J,nrow=J),matrix(d11a[[30]],ncol=J,nrow=J),matrix(d11a[[31]],ncol=J,nrow=J),matrix(d11a[[32]],ncol=J,nrow=J),matrix(d11a[[33]],ncol=J,nrow=J),
matrix(d11a[[34]],ncol=J,nrow=J),matrix(d11a[[35]],ncol=J,nrow=J),matrix(d11a[[36]],ncol=J,nrow=J))
d1<-do.call(adiag,d1temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d1temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J),matrix(d11a2[[6]],ncol=J,nrow=J),matrix(d11a2[[7]],ncol=J,nrow=J),matrix(d11a2[[8]],ncol=J,nrow=J),matrix(d11a2[[9]],ncol=J,nrow=J),
matrix(d11a2[[10]],ncol=J,nrow=J),matrix(d11a2[[11]],ncol=J,nrow=J),matrix(d11a2[[12]],ncol=J,nrow=J),matrix(d11a2[[13]],ncol=J,nrow=J),matrix(d11a2[[14]],ncol=J,nrow=J),matrix(d11a2[[15]],ncol=J,nrow=J),matrix(d11a2[[16]],ncol=J,nrow=J),matrix(d11a2[[17]],ncol=J,nrow=J),
matrix(d11a2[[18]],ncol=J,nrow=J),matrix(d11a2[[19]],ncol=J,nrow=J),matrix(d11a2[[20]],ncol=J,nrow=J),matrix(d11a2[[21]],ncol=J,nrow=J),
matrix(d11a2[[22]],ncol=J,nrow=J),matrix(d11a2[[23]],ncol=J,nrow=J),matrix(d11a2[[24]],ncol=J,nrow=J),matrix(d11a2[[25]],ncol=J,nrow=J),matrix(d11a2[[26]],ncol=J,nrow=J),matrix(d11a2[[27]],ncol=J,nrow=J),
matrix(d11a2[[28]],ncol=J,nrow=J),matrix(d11a2[[29]],ncol=J,nrow=J),matrix(d11a2[[30]],ncol=J,nrow=J),matrix(d11a2[[31]],ncol=J,nrow=J),matrix(d11a2[[32]],ncol=J,nrow=J),matrix(d11a2[[33]],ncol=J,nrow=J),
matrix(d11a2[[34]],ncol=J,nrow=J),matrix(d11a2[[35]],ncol=J,nrow=J),matrix(d11a2[[36]],ncol=J,nrow=J))
d12<-do.call(adiag,d1temp2)


sig1<-vy*c+beta^2*eta*d1+beta^2*vx*d12

sig1i<-solve(sig1)
b1<-t(a1)%*%sig1i%*%a1
b1i<-solve(b1)
i1<-I
j1<-J
#############
# design 2 (I=18 different values, replicated J = 2 times)
I <-f[length(f)-1]
J<-n/I
th<-c(rep(a,J),rep(a+(b-a)/(I-1),J),rep(a+2*(b-a)/(I-1),J),rep(a+3*(b-a)/(I-1),J),rep(a+4*(b-a)/(I-1),J),rep(a+5*(b-a)/(I-1),J),rep(a+6*(b-a)/(I-1),J),rep(a+7*(b-a)/(I-1),J),rep(a+8*(b-a)/(I-1),J),
rep(a+9*(b-a)/(I-1),J),rep(a+10*(b-a)/(I-1),J),rep(a+11*(b-a)/(I-1),J),rep(a+12*(b-a)/(I-1),J),rep(a+13*(b-a)/(I-1),J),rep(a+14*(b-a)/(I-1),J),rep(a+15*(b-a)/(I-1),J),rep(a+16*(b-a)/(I-1),J),rep(a+17*(b-a)/(I-1),J))
 
a2<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],
th[20],th[21],th[22],th[23],th[24],th[25],th[26],th[27],th[28],th[29],th[30],th[31],th[32],th[33],th[34],th[35],th[36]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+4*(b-a)/(I-1))^2,nrow=J,ncol=J),
matrix((a+5*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+6*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+7*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+8*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+9*(b-a)/(I-1))^2,nrow=J,ncol=J),
matrix((a+10*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+11*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+12*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+13*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+14*(b-a)/(I-1))^2,nrow=J,ncol=J),
matrix((a+15*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+16*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+17*(b-a)/(I-1))^2,nrow=J,ncol=J))

d2temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J),matrix(d11a[[6]],ncol=J,nrow=J),
matrix(d11a[[7]],ncol=J,nrow=J),matrix(d11a[[8]],ncol=J,nrow=J),matrix(d11a[[9]],ncol=J,nrow=J),matrix(d11a[[10]],ncol=J,nrow=J),matrix(d11a[[11]],ncol=J,nrow=J),matrix(d11a[[12]],ncol=J,nrow=J),
matrix(d11a[[13]],ncol=J,nrow=J),matrix(d11a[[14]],ncol=J,nrow=J),matrix(d11a[[15]],ncol=J,nrow=J),matrix(d11a[[16]],ncol=J,nrow=J),matrix(d11a[[17]],ncol=J,nrow=J),matrix(d11a[[18]],ncol=J,nrow=J))

d2<-do.call(adiag,d2temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d2temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J),matrix(d11a2[[6]],ncol=J,nrow=J),
matrix(d11a2[[7]],ncol=J,nrow=J),matrix(d11a2[[8]],ncol=J,nrow=J),matrix(d11a2[[9]],ncol=J,nrow=J),matrix(d11a2[[10]],ncol=J,nrow=J),matrix(d11a2[[11]],ncol=J,nrow=J),matrix(d11a2[[12]],ncol=J,nrow=J),
matrix(d11a2[[13]],ncol=J,nrow=J),matrix(d11a2[[14]],ncol=J,nrow=J),matrix(d11a2[[15]],ncol=J,nrow=J),matrix(d11a2[[16]],ncol=J,nrow=J),matrix(d11a2[[17]],ncol=J,nrow=J),matrix(d11a2[[18]],ncol=J,nrow=J))

d22<-do.call(adiag,d2temp2)


sig2<-vy*c+beta^2*eta*d2+beta^2*vx*d22
sig2i<-solve(sig2)
b2<-t(a2)%*%sig2i%*%a2
b2i<-solve(b2)
i2<-I
j2<-J
#############
# design 3 (I = 12 different values, replicated J=3 times)
I<-f[length(f)-2]
J<-n/I
th<-c(rep(a,J), rep(a+(b-a)/(I-1),J), rep(a+2*(b-a)/(I-1),J),rep(a+3*(b-a)/(I-1),J),rep(a+4*(b-a)/(I-1),J), rep(a+5*(b-a)/(I-1),J),rep(a+6*(b-a)/(I-1),J),rep(a+7*(b-a)/(I-1),J),rep(a+8*(b-a)/(I-1),J),
rep(a+9*(b-a)/(I-1),J),rep(a+10*(b-a)/(I-1),J),rep(a+11*(b-a)/(I-1),J))

a3<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],
th[20],th[21],th[22],th[23],th[24],th[25],th[26],th[27],th[28],th[29],th[30],th[31],th[32],th[33],th[34],th[35],th[36]), ncol=2)


d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+4*(b-a)/(I-1))^2,nrow=J,ncol=J),
matrix((a+5*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+6*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+7*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+8*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+9*(b-a)/(I-1))^2,nrow=J,ncol=J),
matrix((a+10*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+11*(b-a)/(I-1))^2,nrow=J,ncol=J))

d3temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J),matrix(d11a[[6]],ncol=J,nrow=J),
matrix(d11a[[7]],ncol=J,nrow=J),matrix(d11a[[8]],ncol=J,nrow=J),matrix(d11a[[9]],ncol=J,nrow=J),matrix(d11a[[10]],ncol=J,nrow=J),matrix(d11a[[11]],ncol=J,nrow=J),matrix(d11a[[12]],ncol=J,nrow=J))

d3<-do.call(adiag,d3temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d3temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J),
matrix(d11a2[[6]],ncol=J,nrow=J),matrix(d11a2[[7]],ncol=J,nrow=J),matrix(d11a2[[8]],ncol=J,nrow=J),matrix(d11a2[[9]],ncol=J,nrow=J),matrix(d11a2[[10]],ncol=J,nrow=J),
matrix(d11a2[[11]],ncol=J,nrow=J),matrix(d11a2[[12]],ncol=J,nrow=J))

d32<-do.call(adiag,d3temp2)



sig3<-vy*c+beta^2*eta*d3+beta^2*vx*d32
sig3i<-solve(sig3)
b3<-t(a3)%*%sig3i%*%a3
b3i<-solve(b3)
i3<-I
j3<-J
#################
# design 4 (I = 9 different values, replicated J = 4 times)
I<-f[length(f)-3]
J<-n/I
th<-c(rep(a,J), rep(a+(b-a)/(I-1),J), rep(a+2*(b-a)/(I-1),J),rep(a+3*(b-a)/(I-1),J),rep(a+4*(b-a)/(I-1),J), rep(a+5*(b-a)/(I-1),J),rep(a+6*(b-a)/(I-1),J), rep(a+7*(b-a)/(I-1),J),rep(a+8*(b-a)/(I-1),J))
a4<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],
th[20],th[21],th[22],th[23],th[24],th[25],th[26],th[27],th[28],th[29],th[30],th[31],th[32],th[33],th[34],th[35],th[36]), ncol=2)


d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+4*(b-a)/(I-1))^2,nrow=J,ncol=J),
matrix((a+5*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+6*(b-a)/(I-1))^2,nrow=J,ncol=J),
matrix((a+7*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+8*(b-a)/(I-1))^2,nrow=J,ncol=J))

d4temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J),matrix(d11a[[6]],ncol=J,nrow=J),
matrix(d11a[[7]],ncol=J,nrow=J),matrix(d11a[[8]],ncol=J,nrow=J),matrix(d11a[[9]],ncol=J,nrow=J))

d4<-do.call(adiag,d4temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d4temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J),matrix(d11a2[[6]],ncol=J,nrow=J),
matrix(d11a2[[7]],ncol=J,nrow=J),matrix(d11a2[[8]],ncol=J,nrow=J),matrix(d11a2[[9]],ncol=J,nrow=J))

d42<-do.call(adiag,d4temp2)


sig4<-vy*c+beta^2*eta*d4+beta^2*vx*d42
sig4i<-solve(sig4)
b4<-t(a4)%*%sig4i%*%a4
b4i<-solve(b4)
i4<-I
j4<-J
j4<-J
#################
# design 5 (I = 6 different values, replicated J = 6 times)
I<-f[length(f)-4]
J<-n/I
th<-c(rep(a,J), rep(a+(b-a)/(I-1),J), rep(a+2*(b-a)/(I-1),J),rep(a+3*(b-a)/(I-1),J),rep(a+4*(b-a)/(I-1),J),rep(a+5*(b-a)/(I-1),J))

a5<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],
th[20],th[21],th[22],th[23],th[24],th[25],th[26],th[27],th[28],th[29],th[30],th[31],th[32],th[33],th[34],th[35],th[36]), ncol=2)


d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+4*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+5*(b-a)/(I-1))^2,nrow=J,ncol=J))

d5temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),
matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J),matrix(d11a[[6]],ncol=J,nrow=J))

d5<-do.call(adiag,d5temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d5temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),
matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J),matrix(d11a2[[6]],ncol=J,nrow=J))

d52<-do.call(adiag,d5temp2)


sig5<-vy*c+beta^2*eta*d5+beta^2*vx*d52
sig5i<-solve(sig5)
b5<-t(a5)%*%sig5i%*%a5
b5i<-solve(b5)
i5<-I
j5<-J
##################
# design 6 (I = 4 different values, replicated J = 9 times)
I<-f[length(f)-5]
J<-n/I

th<-c(rep(a,J), rep(a+(b-a)/(I-1),J), rep(a+2*(b-a)/(I-1),J),rep(a+3*(b-a)/(I-1),J))

a6<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],
th[20],th[21],th[22],th[23],th[24],th[25],th[26],th[27],th[28],th[29],th[30],th[31],th[32],th[33],th[34],th[35],th[36]), ncol=2)


d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(I-1))^2,nrow=J,ncol=J))

d6temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J))

d6<-do.call(adiag,d6temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d6temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J))

d62<-do.call(adiag,d6temp2)


sig6<-vy*c+beta^2*eta*d6+beta^2*vx*d62
sig6i<-solve(sig6)
b6<-t(a6)%*%sig6i%*%a6
b6i<-solve(b6)
i6<-I
j6<-J

##################
# design 7 (I = 3 different values, replicated J = 12 times)
I<-f[length(f)-6]
J<-n/I
th<-c(rep(a,J), rep(a+(b-a)/(I-1),J), rep(a+2*(b-a)/(I-1),J))

a7<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],
th[20],th[21],th[22],th[23],th[24],th[25],th[26],th[27],th[28],th[29],th[30],th[31],th[32],th[33],th[34],th[35],th[36]), ncol=2)


d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(I-1))^2,nrow=J,ncol=J))

d7temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J))

d7<-do.call(adiag,d7temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d7temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J))

d72<-do.call(adiag,d7temp2)


sig7<-vy*c+beta^2*eta*d7+beta^2*vx*d72
sig7i<-solve(sig7)
b7<-t(a7)%*%sig7i%*%a7
b7i<-solve(b7)
i7<-I
j7<-J
#################
#design 8 (I = 2 different values, replicated J = 18 times)
I<-f[length(f)-7]
J<-n/I
th<-c(rep(a,J), rep(a+(b-a)/(I-1),J))

a8<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],
th[20],th[21],th[22],th[23],th[24],th[25],th[26],th[27],th[28],th[29],th[30],th[31],th[32],th[33],th[34],th[35],th[36]), ncol=2)


d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J))

d8temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J))

d8<-do.call(adiag,d8temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d8temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J))

d82<-do.call(adiag,d8temp2)


sig8<-vy*c+beta^2*eta*d8+beta^2*vx*d82
sig8i<-solve(sig8)
b8<-t(a8)%*%sig8i%*%a8
b8i<-solve(b8)
i8<-I
j8<-J

deter<-c(det(b1i),det(b2i),det(b3i),det(b4i),det(b5i),det(b6i),det(b7i),det(b8i))
which(deter==min(deter))

eig1<-eigen(b1i)
tr1<-sum(eig1$values)
eig2<-eigen(b2i)
tr2<-sum(eig2$values)
eig3<-eigen(b3i)
tr3<-sum(eig3$values)
eig4<-eigen(b4i)
tr4<-sum(eig4$values)
eig5<-eigen(b5i)
tr5<-sum(eig5$values)
eig6<-eigen(b6i)
tr6<-sum(eig6$values)
eig7<-eigen(b7i)
tr7<-sum(eig7$values)
eig8<-eigen(b8i)
tr8<-sum(eig8$values)

 
trace<-c(tr1,tr2,tr3,tr4,tr5,tr6,tr7,tr8)
best<-which(trace==min(trace))


vb<-c((b1i[[2,2]])^0.5,(b2i[[2,2]])^0.5,(b3i[[2,2]])^0.5,(b4i[[2,2]])^0.5,(b5i[[2,2]])^0.5,(b6i[[2,2]])^0.5,(b7i[[2,2]])^0.5,(b8i[[2,2]])^0.5)
bestb<-which(vb==min(vb))

cq<-diag(1,nq)
thq<-c(rep(beta,nq))
th2q<-matrix(c(1,xnew),nrow=1)
a5q<-matrix(c(rep(beta,nq)), ncol=1)
varab1<-th2q%*%b1i%*%t(th2q)
d11aq1<-list(matrix(varab1, nrow=nq,ncol=nq))
dq1<-do.call(adiag,d11aq1)
sigq1<-dq1
sigq1i<-solve(sigq1)
bq1<-t(a5q)%*%sigq1i%*%a5q 
stuq1<-(1/bq1+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b1i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq1<-stuq1/xnew*100

varab2<-th2q%*%b2i%*%t(th2q)
d11aq2<-list(matrix(varab2, nrow=nq,ncol=nq))
dq2<-do.call(adiag,d11aq2)
sigq2<-dq2
sigq2i<-solve(sigq2)
bq2<-t(a5q)%*%sigq2i%*%a5q 
stuq2<-(1/bq2+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b2i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq2<-stuq2/xnew*100

varab3<-th2q%*%b3i%*%t(th2q)
d11aq3<-list(matrix(varab3, nrow=nq,ncol=nq))
dq3<-do.call(adiag,d11aq3)
sigq3<-dq3
sigq3i<-solve(sigq3)
bq3<-t(a5q)%*%sigq3i%*%a5q 
stuq3<-(1/bq3+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b3i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq3<-stuq3/xnew*100

varab4<-th2q%*%b4i%*%t(th2q)
d11aq4<-list(matrix(varab4, nrow=nq,ncol=nq))
dq4<-do.call(adiag,d11aq4)
sigq4<-dq4
sigq4i<-solve(sigq4)
bq4<-t(a5q)%*%sigq4i%*%a5q 
stuq4<-(1/bq4+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b4i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq4<-stuq4/xnew*100

varab5<-th2q%*%b5i%*%t(th2q)
d11aq5<-list(matrix(varab5, nrow=nq,ncol=nq))
dq5<-do.call(adiag,d11aq5)
sigq5<-dq5
sigq5i<-solve(sigq5)
bq5<-t(a5q)%*%sigq5i%*%a5q 
stuq5<-(1/bq5+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b5i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq5<-stuq5/xnew*100

varab6<-th2q%*%b6i%*%t(th2q)
d11aq6<-list(matrix(varab6, nrow=nq,ncol=nq))
dq6<-do.call(adiag,d11aq6)
sigq6<-dq6
sigq6i<-solve(sigq6)
bq6<-t(a5q)%*%sigq6i%*%a5q 
stuq6<-(1/bq6+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b6i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq6<-stuq6/xnew*100

varab7<-th2q%*%b7i%*%t(th2q)
d11aq7<-list(matrix(varab7, nrow=nq,ncol=nq))
dq7<-do.call(adiag,d11aq7)
sigq7<-dq7
sigq7i<-solve(sigq7)
bq7<-t(a5q)%*%sigq7i%*%a5q 
stuq7<-(1/bq7+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b7i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq7<-stuq7/xnew*100

varab8<-th2q%*%b8i%*%t(th2q)
d11aq8<-list(matrix(varab8, nrow=nq,ncol=nq))
dq8<-do.call(adiag,d11aq8)
sigq8<-dq8
sigq8i<-solve(sigq8)
bq8<-t(a5q)%*%sigq8i%*%a5q 
stuq8<-(1/bq8+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b8i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq8<-stuq8/xnew*100



finalrelst<-c(relstuq1,relstuq2,relstuq3,relstuq4,relstuq5,relstuq6,relstuq7,relstuq8)
finbest<-which(finalrelst==min(finalrelst))

finali<-c(i1,i2,i3,i4,i5,i6,i7,i8)
finalj<-c(j1,j2,j3,j4,j5,j6,j7,j8)
jbest<-finalj[finbest]
ibest<-finali[finbest]
totrep<-n+ny
finstd<-finalrelst[finbest]
samplrep<-ny/r
result<-list(samplrep=samplrep,nosamples=r,finstd=finstd,ibest=ibest,jbest=jbest,totrep=totrep)

return(result)
}

#########################################################################################

#########################################################################################
############################################
caldesign140<-function(n,a,b,vy,vx,beta,eta,ny,vs,xnew,vi){
f<-divisors(n)
fy<-divisors(ny)
r<-1
if(vs>0 |vi>0){r=ny}
c<-diag(1,n)
v1<-(beta^2*eta*a^2)^0.5
v2<-(beta^2*eta*b^2)^0.5
nq<-1
#########
# #########
# design 1 (all I = 40 values are different, replicated J = 1 time)

J<-1
I<-f[length(f)]
th<-c(a,a+(b-a)/(I-1),a+2*(b-a)/(I-1),a+3*(b-a)/(I-1),a+4*(b-a)/(I-1),
                  a+5*(b-a)/(I-1),a+6*(b-a)/(I-1),a+7*(b-a)/(I-1),a+8*(b-a)/(I-1),a+9*(b-a)/(I-1),a+10*(b-a)/(I-1),a+11*(b-a)/(I-1),
a+12*(b-a)/(I-1),a+13*(b-a)/(I-1),a+14*(b-a)/(I-1),a+15*(b-a)/(I-1),
                  a+16*(b-a)/(I-1),a+17*(b-a)/(I-1),a+18*(b-a)/(I-1),a+19*(b-a)/(I-1),a+20*(b-a)/(I-1),a+21*(b-a)/(I-1),a+22*(b-a)/(I-1),
a+23*(b-a)/(I-1),a+24*(b-a)/(I-1),a+25*(b-a)/(I-1),a+26*(b-a)/(I-1),a+27*(b-a)/(I-1),a+28*(b-a)/(I-1),a+29*(b-a)/(I-1),a+30*(b-a)/(I-1),
a+31*(b-a)/(I-1),a+32*(b-a)/(I-1),a+33*(b-a)/(I-1),a+34*(b-a)/(I-1),a+35*(b-a)/(I-1),a+36*(b-a)/(I-1),a+37*(b-a)/(I-1),a+38*(b-a)/(I-1),b)
 
a1<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],
th[20],th[21],th[22],th[23],th[24],th[25],th[26],th[27],th[28],th[29],th[30],th[31],th[32],th[33],th[34],th[35],th[36],th[37],th[38],th[39],th[40]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+2*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+4*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+5*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+6*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+7*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+8*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+9*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+10*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+11*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+12*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+13*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+14*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+15*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+16*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+17*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+18*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+19*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+20*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+21*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+22*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+23*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+24*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+25*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+26*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+27*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+28*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+29*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+30*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+31*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+32*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+33*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+34*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+35*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+36*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+37*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+38*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+39*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J))

d1temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J),matrix(d11a[[6]],ncol=J,nrow=J),matrix(d11a[[7]],ncol=J,nrow=J),matrix(d11a[[8]],ncol=J,nrow=J),matrix(d11a[[9]],ncol=J,nrow=J),
matrix(d11a[[10]],ncol=J,nrow=J),matrix(d11a[[11]],ncol=J,nrow=J),matrix(d11a[[12]],ncol=J,nrow=J),matrix(d11a[[13]],ncol=J,nrow=J),matrix(d11a[[14]],ncol=J,nrow=J),matrix(d11a[[15]],ncol=J,nrow=J),matrix(d11a[[16]],ncol=J,nrow=J),matrix(d11a[[17]],ncol=J,nrow=J),
matrix(d11a[[18]],ncol=J,nrow=J),matrix(d11a[[19]],ncol=J,nrow=J),matrix(d11a[[20]],ncol=J,nrow=J),matrix(d11a[[21]],ncol=J,nrow=J),
matrix(d11a[[22]],ncol=J,nrow=J),matrix(d11a[[23]],ncol=J,nrow=J),matrix(d11a[[24]],ncol=J,nrow=J),matrix(d11a[[25]],ncol=J,nrow=J),matrix(d11a[[26]],ncol=J,nrow=J),matrix(d11a[[27]],ncol=J,nrow=J),
matrix(d11a[[28]],ncol=J,nrow=J),matrix(d11a[[29]],ncol=J,nrow=J),matrix(d11a[[30]],ncol=J,nrow=J),matrix(d11a[[31]],ncol=J,nrow=J),matrix(d11a[[32]],ncol=J,nrow=J),matrix(d11a[[33]],ncol=J,nrow=J),
matrix(d11a[[34]],ncol=J,nrow=J),matrix(d11a[[35]],ncol=J,nrow=J),matrix(d11a[[36]],ncol=J,nrow=J),matrix(d11a[[37]],ncol=J,nrow=J),
matrix(d11a[[38]],ncol=J,nrow=J),matrix(d11a[[39]],ncol=J,nrow=J),matrix(d11a[[40]],ncol=J,nrow=J))
d1<-do.call(adiag,d1temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d1temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J),matrix(d11a2[[6]],ncol=J,nrow=J),matrix(d11a2[[7]],ncol=J,nrow=J),matrix(d11a2[[8]],ncol=J,nrow=J),matrix(d11a2[[9]],ncol=J,nrow=J),
matrix(d11a2[[10]],ncol=J,nrow=J),matrix(d11a2[[11]],ncol=J,nrow=J),matrix(d11a2[[12]],ncol=J,nrow=J),matrix(d11a2[[13]],ncol=J,nrow=J),matrix(d11a2[[14]],ncol=J,nrow=J),matrix(d11a2[[15]],ncol=J,nrow=J),matrix(d11a2[[16]],ncol=J,nrow=J),matrix(d11a2[[17]],ncol=J,nrow=J),
matrix(d11a2[[18]],ncol=J,nrow=J),matrix(d11a2[[19]],ncol=J,nrow=J),matrix(d11a2[[20]],ncol=J,nrow=J),matrix(d11a2[[21]],ncol=J,nrow=J),
matrix(d11a2[[22]],ncol=J,nrow=J),matrix(d11a2[[23]],ncol=J,nrow=J),matrix(d11a2[[24]],ncol=J,nrow=J),matrix(d11a2[[25]],ncol=J,nrow=J),matrix(d11a2[[26]],ncol=J,nrow=J),matrix(d11a2[[27]],ncol=J,nrow=J),
matrix(d11a2[[28]],ncol=J,nrow=J),matrix(d11a2[[29]],ncol=J,nrow=J),matrix(d11a2[[30]],ncol=J,nrow=J),matrix(d11a2[[31]],ncol=J,nrow=J),matrix(d11a2[[32]],ncol=J,nrow=J),matrix(d11a2[[33]],ncol=J,nrow=J),
matrix(d11a2[[34]],ncol=J,nrow=J),matrix(d11a2[[35]],ncol=J,nrow=J),matrix(d11a2[[36]],ncol=J,nrow=J),matrix(d11a2[[37]],ncol=J,nrow=J),
matrix(d11a2[[38]],ncol=J,nrow=J),matrix(d11a2[[39]],ncol=J,nrow=J),matrix(d11a2[[40]],ncol=J,nrow=J))
d12<-do.call(adiag,d1temp2)


sig1<-vy*c+beta^2*eta*d1+beta^2*vx*d12

sig1i<-solve(sig1)
b1<-t(a1)%*%sig1i%*%a1
b1i<-solve(b1)
i1<-I
j1<-J
#############
# design 2 (I=20 different values, replicated J = 2 times)
I <-f[length(f)-1]
J<-n/I
th<-c(rep(a,J),rep(a+(b-a)/(I-1),J),rep(a+2*(b-a)/(I-1),J),rep(a+3*(b-a)/(I-1),J),rep(a+4*(b-a)/(I-1),J),rep(a+5*(b-a)/(I-1),J),rep(a+6*(b-a)/(I-1),J),rep(a+7*(b-a)/(I-1),J),rep(a+8*(b-a)/(I-1),J),
rep(a+9*(b-a)/(I-1),J),rep(a+10*(b-a)/(I-1),J),rep(a+11*(b-a)/(I-1),J),rep(a+12*(b-a)/(I-1),J),rep(a+13*(b-a)/(I-1),J),rep(a+14*(b-a)/(I-1),J),rep(a+15*(b-a)/(I-1),J),rep(a+16*(b-a)/(I-1),J),
rep(a+17*(b-a)/(I-1),J),rep(a+18*(b-a)/(I-1),J),
rep(a+19*(b-a)/(I-1),J))
 
a2<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],
th[20],th[21],th[22],th[23],th[24],th[25],th[26],th[27],th[28],th[29],th[30],th[31],th[32],th[33],th[34],th[35],th[36],th[37],th[38],th[39],th[40]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+4*(b-a)/(I-1))^2,nrow=J,ncol=J),
matrix((a+5*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+6*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+7*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+8*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+9*(b-a)/(I-1))^2,nrow=J,ncol=J),
matrix((a+10*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+11*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+12*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+13*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+14*(b-a)/(I-1))^2,nrow=J,ncol=J),
matrix((a+15*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+16*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+17*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+18*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+19*(b-a)/(I-1))^2,nrow=J,ncol=J))

d2temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J),matrix(d11a[[6]],ncol=J,nrow=J),
matrix(d11a[[7]],ncol=J,nrow=J),matrix(d11a[[8]],ncol=J,nrow=J),matrix(d11a[[9]],ncol=J,nrow=J),matrix(d11a[[10]],ncol=J,nrow=J),matrix(d11a[[11]],ncol=J,nrow=J),matrix(d11a[[12]],ncol=J,nrow=J),
matrix(d11a[[13]],ncol=J,nrow=J),matrix(d11a[[14]],ncol=J,nrow=J),matrix(d11a[[15]],ncol=J,nrow=J),matrix(d11a[[16]],ncol=J,nrow=J),matrix(d11a[[17]],ncol=J,nrow=J),matrix(d11a[[18]],ncol=J,nrow=J),
matrix(d11a[[19]],ncol=J,nrow=J),matrix(d11a[[20]],ncol=J,nrow=J))

d2<-do.call(adiag,d2temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d2temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J),matrix(d11a2[[6]],ncol=J,nrow=J),
matrix(d11a2[[7]],ncol=J,nrow=J),matrix(d11a2[[8]],ncol=J,nrow=J),matrix(d11a2[[9]],ncol=J,nrow=J),matrix(d11a2[[10]],ncol=J,nrow=J),matrix(d11a2[[11]],ncol=J,nrow=J),matrix(d11a2[[12]],ncol=J,nrow=J),
matrix(d11a2[[13]],ncol=J,nrow=J),matrix(d11a2[[14]],ncol=J,nrow=J),matrix(d11a2[[15]],ncol=J,nrow=J),matrix(d11a2[[16]],ncol=J,nrow=J),matrix(d11a2[[17]],ncol=J,nrow=J),matrix(d11a2[[18]],ncol=J,nrow=J),
matrix(d11a2[[19]],ncol=J,nrow=J),matrix(d11a2[[20]],ncol=J,nrow=J))

d22<-do.call(adiag,d2temp2)


sig2<-vy*c+beta^2*eta*d2+beta^2*vx*d22
sig2i<-solve(sig2)
b2<-t(a2)%*%sig2i%*%a2
b2i<-solve(b2)
i2<-I
j2<-J
#############
# design 3 (I = 10 different values, replicated J=4 times)
I<-f[length(f)-2]
J<-n/I
th<-c(rep(a,J), rep(a+(b-a)/(I-1),J), rep(a+2*(b-a)/(I-1),J),rep(a+3*(b-a)/(I-1),J),rep(a+4*(b-a)/(I-1),J), rep(a+5*(b-a)/(I-1),J),rep(a+6*(b-a)/(I-1),J),rep(a+7*(b-a)/(I-1),J),rep(a+8*(b-a)/(I-1),J),
rep(a+9*(b-a)/(I-1),J))

a3<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],
th[20],th[21],th[22],th[23],th[24],th[25],th[26],th[27],th[28],th[29],th[30],th[31],th[32],th[33],th[34],th[35],th[36],th[37],th[38],th[39],th[40]), ncol=2)


d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+4*(b-a)/(I-1))^2,nrow=J,ncol=J),
matrix((a+5*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+6*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+7*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+8*(b-a)/(I-1))^2,nrow=J,ncol=J),
matrix((a+9*(b-a)/(I-1))^2,nrow=J,ncol=J))

d3temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J),matrix(d11a[[6]],ncol=J,nrow=J),
matrix(d11a[[7]],ncol=J,nrow=J),matrix(d11a[[8]],ncol=J,nrow=J),matrix(d11a[[9]],ncol=J,nrow=J),matrix(d11a[[10]],ncol=J,nrow=J))

d3<-do.call(adiag,d3temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d3temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J),
matrix(d11a2[[6]],ncol=J,nrow=J),matrix(d11a2[[7]],ncol=J,nrow=J),matrix(d11a2[[8]],ncol=J,nrow=J),matrix(d11a2[[9]],ncol=J,nrow=J),matrix(d11a2[[10]],ncol=J,nrow=J))

d32<-do.call(adiag,d3temp2)



sig3<-vy*c+beta^2*eta*d3+beta^2*vx*d32
sig3i<-solve(sig3)
b3<-t(a3)%*%sig3i%*%a3
b3i<-solve(b3)
i3<-I
j3<-J
#################
# design 4 (I = 8 different values, replicated J = 5 times)
I<-f[length(f)-3]
J<-n/I
th<-c(rep(a,J), rep(a+(b-a)/(I-1),J), rep(a+2*(b-a)/(I-1),J),rep(a+3*(b-a)/(I-1),J),rep(a+4*(b-a)/(I-1),J), rep(a+5*(b-a)/(I-1),J),rep(a+6*(b-a)/(I-1),J), rep(a+7*(b-a)/(I-1),J))
a4<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],
th[20],th[21],th[22],th[23],th[24],th[25],th[26],th[27],th[28],th[29],th[30],th[31],th[32],th[33],th[34],th[35],th[36],th[37],th[38],th[39],th[40]), ncol=2)


d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+4*(b-a)/(I-1))^2,nrow=J,ncol=J),
matrix((a+5*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+6*(b-a)/(I-1))^2,nrow=J,ncol=J),
matrix((a+7*(b-a)/(I-1))^2,nrow=J,ncol=J))

d4temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J),matrix(d11a[[6]],ncol=J,nrow=J),
matrix(d11a[[7]],ncol=J,nrow=J),matrix(d11a[[8]],ncol=J,nrow=J))

d4<-do.call(adiag,d4temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d4temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J),matrix(d11a2[[6]],ncol=J,nrow=J),
matrix(d11a2[[7]],ncol=J,nrow=J),matrix(d11a2[[8]],ncol=J,nrow=J))

d42<-do.call(adiag,d4temp2)


sig4<-vy*c+beta^2*eta*d4+beta^2*vx*d42
sig4i<-solve(sig4)
b4<-t(a4)%*%sig4i%*%a4
b4i<-solve(b4)
i4<-I
j4<-J
j4<-J
#################
# design 5 (I = 5 different values, replicated J = 8 times)
I<-f[length(f)-4]
J<-n/I
th<-c(rep(a,J), rep(a+(b-a)/(I-1),J), rep(a+2*(b-a)/(I-1),J),rep(a+3*(b-a)/(I-1),J),rep(a+4*(b-a)/(I-1),J))

a5<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],
th[20],th[21],th[22],th[23],th[24],th[25],th[26],th[27],th[28],th[29],th[30],th[31],th[32],th[33],th[34],th[35],th[36],th[37],th[38],th[39],th[40]), ncol=2)


d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(I-1))^2,nrow=J,ncol=J),
matrix((a+4*(b-a)/(I-1))^2,nrow=J,ncol=J))

d5temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),
matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J))

d5<-do.call(adiag,d5temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d5temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),
matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J))

d52<-do.call(adiag,d5temp2)


sig5<-vy*c+beta^2*eta*d5+beta^2*vx*d52
sig5i<-solve(sig5)
b5<-t(a5)%*%sig5i%*%a5
b5i<-solve(b5)
i5<-I
j5<-J
##################
# design 6 (I = 4 different values, replicated J = 10 times)
I<-f[length(f)-5]
J<-n/I

th<-c(rep(a,J), rep(a+(b-a)/(I-1),J), rep(a+2*(b-a)/(I-1),J),rep(a+3*(b-a)/(I-1),J))

a6<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],
th[20],th[21],th[22],th[23],th[24],th[25],th[26],th[27],th[28],th[29],th[30],th[31],th[32],th[33],th[34],th[35],th[36],th[37],th[38],th[39],th[40]), ncol=2)


d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(I-1))^2,nrow=J,ncol=J))

d6temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J))

d6<-do.call(adiag,d6temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d6temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J))

d62<-do.call(adiag,d6temp2)


sig6<-vy*c+beta^2*eta*d6+beta^2*vx*d62
sig6i<-solve(sig6)
b6<-t(a6)%*%sig6i%*%a6
b6i<-solve(b6)
i6<-I
j6<-J

##################
# design 7 (I = 2 different values, replicated J = 20 times)
I<-f[length(f)-6]
J<-n/I
th<-c(rep(a,J), rep(a+(b-a)/(I-1),J))

a7<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],
th[20],th[21],th[22],th[23],th[24],th[25],th[26],th[27],th[28],th[29],th[30],th[31],th[32],th[33],th[34],th[35],th[36],th[37],th[38],th[39],th[40]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J))

d7temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J))

d7<-do.call(adiag,d7temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d7temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J))

d72<-do.call(adiag,d7temp2)


sig7<-vy*c+beta^2*eta*d7+beta^2*vx*d72
sig7i<-solve(sig7)
b7<-t(a7)%*%sig7i%*%a7
b7i<-solve(b7)
i7<-I
j7<-J
#################


deter<-c(det(b1i),det(b2i),det(b3i),det(b4i),det(b5i),det(b6i),det(b7i))
which(deter==min(deter))

eig1<-eigen(b1i)
tr1<-sum(eig1$values)
eig2<-eigen(b2i)
tr2<-sum(eig2$values)
eig3<-eigen(b3i)
tr3<-sum(eig3$values)
eig4<-eigen(b4i)
tr4<-sum(eig4$values)
eig5<-eigen(b5i)
tr5<-sum(eig5$values)
eig6<-eigen(b6i)
tr6<-sum(eig6$values)
eig7<-eigen(b7i)
tr7<-sum(eig7$values)

 
trace<-c(tr1,tr2,tr3,tr4,tr5,tr6,tr7)
best<-which(trace==min(trace))


vb<-c((b1i[[2,2]])^0.5,(b2i[[2,2]])^0.5,(b3i[[2,2]])^0.5,(b4i[[2,2]])^0.5,(b5i[[2,2]])^0.5,(b6i[[2,2]])^0.5,(b7i[[2,2]])^0.5)
bestb<-which(vb==min(vb))

cq<-diag(1,nq)
thq<-c(rep(beta,nq))
th2q<-matrix(c(1,xnew),nrow=1)
a5q<-matrix(c(rep(beta,nq)), ncol=1)
varab1<-th2q%*%b1i%*%t(th2q)
d11aq1<-list(matrix(varab1, nrow=nq,ncol=nq))
dq1<-do.call(adiag,d11aq1)
sigq1<-dq1
sigq1i<-solve(sigq1)
bq1<-t(a5q)%*%sigq1i%*%a5q 
stuq1<-(1/bq1+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b1i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq1<-stuq1/xnew*100

varab2<-th2q%*%b2i%*%t(th2q)
d11aq2<-list(matrix(varab2, nrow=nq,ncol=nq))
dq2<-do.call(adiag,d11aq2)
sigq2<-dq2
sigq2i<-solve(sigq2)
bq2<-t(a5q)%*%sigq2i%*%a5q 
stuq2<-(1/bq2+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b2i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq2<-stuq2/xnew*100

varab3<-th2q%*%b3i%*%t(th2q)
d11aq3<-list(matrix(varab3, nrow=nq,ncol=nq))
dq3<-do.call(adiag,d11aq3)
sigq3<-dq3
sigq3i<-solve(sigq3)
bq3<-t(a5q)%*%sigq3i%*%a5q 
stuq3<-(1/bq3+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b3i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq3<-stuq3/xnew*100

varab4<-th2q%*%b4i%*%t(th2q)
d11aq4<-list(matrix(varab4, nrow=nq,ncol=nq))
dq4<-do.call(adiag,d11aq4)
sigq4<-dq4
sigq4i<-solve(sigq4)
bq4<-t(a5q)%*%sigq4i%*%a5q 
stuq4<-(1/bq4+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b4i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq4<-stuq4/xnew*100

varab5<-th2q%*%b5i%*%t(th2q)
d11aq5<-list(matrix(varab5, nrow=nq,ncol=nq))
dq5<-do.call(adiag,d11aq5)
sigq5<-dq5
sigq5i<-solve(sigq5)
bq5<-t(a5q)%*%sigq5i%*%a5q 
stuq5<-(1/bq5+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b5i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq5<-stuq5/xnew*100

varab6<-th2q%*%b6i%*%t(th2q)
d11aq6<-list(matrix(varab6, nrow=nq,ncol=nq))
dq6<-do.call(adiag,d11aq6)
sigq6<-dq6
sigq6i<-solve(sigq6)
bq6<-t(a5q)%*%sigq6i%*%a5q 
stuq6<-(1/bq6+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b6i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq6<-stuq6/xnew*100

varab7<-th2q%*%b7i%*%t(th2q)
d11aq7<-list(matrix(varab7, nrow=nq,ncol=nq))
dq7<-do.call(adiag,d11aq7)
sigq7<-dq7
sigq7i<-solve(sigq7)
bq7<-t(a5q)%*%sigq7i%*%a5q 
stuq7<-(1/bq7+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b7i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq7<-stuq7/xnew*100




finalrelst<-c(relstuq1,relstuq2,relstuq3,relstuq4,relstuq5,relstuq6,relstuq7)
finbest<-which(finalrelst==min(finalrelst))

finali<-c(i1,i2,i3,i4,i5,i6,i7)
finalj<-c(j1,j2,j3,j4,j5,j6,j7)
jbest<-finalj[finbest]
ibest<-finali[finbest]
totrep<-n+ny
finstd<-finalrelst[finbest]
samplrep<-ny/r
result<-list(samplrep=samplrep,nosamples=r,finstd=finstd,ibest=ibest,jbest=jbest,totrep=totrep)

return(result)
}

#########################################################################################


########################################################

#########################################################################################
############################################
caldesign144<-function(n,a,b,vy,vx,beta,eta,ny,vs,xnew,vi){
f<-divisors(n)
fy<-divisors(ny)
r<-1
if(vs>0 |vi>0){r=ny}
c<-diag(1,n)
v1<-(beta^2*eta*a^2)^0.5
v2<-(beta^2*eta*b^2)^0.5
nq<-1
#########
# #########
# design 1 (all I = 44 values are different, replicated J = 1 time)

J<-1
I<-f[length(f)]
th<-c(a,a+(b-a)/(I-1),a+2*(b-a)/(I-1),a+3*(b-a)/(I-1),a+4*(b-a)/(I-1),
                  a+5*(b-a)/(I-1),a+6*(b-a)/(I-1),a+7*(b-a)/(I-1),a+8*(b-a)/(I-1),a+9*(b-a)/(I-1),a+10*(b-a)/(I-1),a+11*(b-a)/(I-1),
a+12*(b-a)/(I-1),a+13*(b-a)/(I-1),a+14*(b-a)/(I-1),a+15*(b-a)/(I-1),
                  a+16*(b-a)/(I-1),a+17*(b-a)/(I-1),a+18*(b-a)/(I-1),a+19*(b-a)/(I-1),a+20*(b-a)/(I-1),a+21*(b-a)/(I-1),a+22*(b-a)/(I-1),
a+23*(b-a)/(I-1),a+24*(b-a)/(I-1),a+25*(b-a)/(I-1),a+26*(b-a)/(I-1),a+27*(b-a)/(I-1),a+28*(b-a)/(I-1),a+29*(b-a)/(I-1),a+30*(b-a)/(I-1),a+31*(b-a)/(I-1),a+32*(b-a)/(I-1),
a+33*(b-a)/(I-1),a+34*(b-a)/(I-1),a+35*(b-a)/(I-1),a+36*(b-a)/(I-1),a+37*(b-a)/(I-1),a+38*(b-a)/(I-1),a+39*(b-a)/(I-1),a+40*(b-a)/(I-1),a+41*(b-a)/(I-1),a+42*(b-a)/(I-1),b)
 
a1<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],
th[20],th[21],th[22],th[23],th[24],th[25],th[26],th[27],th[28],th[29],th[30],th[31],th[32],th[33],th[34],th[35],th[36],th[37],th[38],th[39],th[40],th[41],th[42],th[43],th[44]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+2*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+4*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+5*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+6*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+7*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+8*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+9*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+10*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+11*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+12*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+13*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+14*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+15*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+16*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+17*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+18*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+19*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+20*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+21*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+22*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+23*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+24*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+25*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+26*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+27*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+28*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+29*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+30*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+31*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+32*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+33*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+34*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+35*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+36*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+37*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+38*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+39*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+40*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+41*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+42*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+43*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J))

d1temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J),matrix(d11a[[6]],ncol=J,nrow=J),matrix(d11a[[7]],ncol=J,nrow=J),matrix(d11a[[8]],ncol=J,nrow=J),matrix(d11a[[9]],ncol=J,nrow=J),
matrix(d11a[[10]],ncol=J,nrow=J),matrix(d11a[[11]],ncol=J,nrow=J),matrix(d11a[[12]],ncol=J,nrow=J),matrix(d11a[[13]],ncol=J,nrow=J),matrix(d11a[[14]],ncol=J,nrow=J),matrix(d11a[[15]],ncol=J,nrow=J),matrix(d11a[[16]],ncol=J,nrow=J),matrix(d11a[[17]],ncol=J,nrow=J),
matrix(d11a[[18]],ncol=J,nrow=J),matrix(d11a[[19]],ncol=J,nrow=J),matrix(d11a[[20]],ncol=J,nrow=J),matrix(d11a[[21]],ncol=J,nrow=J),
matrix(d11a[[22]],ncol=J,nrow=J),matrix(d11a[[23]],ncol=J,nrow=J),matrix(d11a[[24]],ncol=J,nrow=J),matrix(d11a[[25]],ncol=J,nrow=J),matrix(d11a[[26]],ncol=J,nrow=J),matrix(d11a[[27]],ncol=J,nrow=J),
matrix(d11a[[28]],ncol=J,nrow=J),matrix(d11a[[29]],ncol=J,nrow=J),matrix(d11a[[30]],ncol=J,nrow=J),matrix(d11a[[31]],ncol=J,nrow=J),matrix(d11a[[32]],ncol=J,nrow=J),matrix(d11a[[33]],ncol=J,nrow=J),
matrix(d11a[[34]],ncol=J,nrow=J),matrix(d11a[[35]],ncol=J,nrow=J),matrix(d11a[[36]],ncol=J,nrow=J),matrix(d11a[[37]],ncol=J,nrow=J),matrix(d11a[[38]],ncol=J,nrow=J),matrix(d11a[[39]],ncol=J,nrow=J),
matrix(d11a[[40]],ncol=J,nrow=J),matrix(d11a[[41]],ncol=J,nrow=J),matrix(d11a[[42]],ncol=J,nrow=J),matrix(d11a[[43]],ncol=J,nrow=J),matrix(d11a[[44]],ncol=J,nrow=J))
d1<-do.call(adiag,d1temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d1temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J),matrix(d11a2[[6]],ncol=J,nrow=J),matrix(d11a2[[7]],ncol=J,nrow=J),matrix(d11a2[[8]],ncol=J,nrow=J),matrix(d11a2[[9]],ncol=J,nrow=J),
matrix(d11a2[[10]],ncol=J,nrow=J),matrix(d11a2[[11]],ncol=J,nrow=J),matrix(d11a2[[12]],ncol=J,nrow=J),matrix(d11a2[[13]],ncol=J,nrow=J),matrix(d11a2[[14]],ncol=J,nrow=J),matrix(d11a2[[15]],ncol=J,nrow=J),matrix(d11a2[[16]],ncol=J,nrow=J),matrix(d11a2[[17]],ncol=J,nrow=J),
matrix(d11a2[[18]],ncol=J,nrow=J),matrix(d11a2[[19]],ncol=J,nrow=J),matrix(d11a2[[20]],ncol=J,nrow=J),matrix(d11a2[[21]],ncol=J,nrow=J),
matrix(d11a2[[22]],ncol=J,nrow=J),matrix(d11a2[[23]],ncol=J,nrow=J),matrix(d11a2[[24]],ncol=J,nrow=J),matrix(d11a2[[25]],ncol=J,nrow=J),matrix(d11a2[[26]],ncol=J,nrow=J),matrix(d11a2[[27]],ncol=J,nrow=J),
matrix(d11a2[[28]],ncol=J,nrow=J),matrix(d11a2[[29]],ncol=J,nrow=J),matrix(d11a2[[30]],ncol=J,nrow=J),matrix(d11a2[[31]],ncol=J,nrow=J),matrix(d11a2[[32]],ncol=J,nrow=J),matrix(d11a2[[33]],ncol=J,nrow=J),
matrix(d11a2[[34]],ncol=J,nrow=J),matrix(d11a2[[35]],ncol=J,nrow=J),matrix(d11a2[[36]],ncol=J,nrow=J),matrix(d11a2[[37]],ncol=J,nrow=J),matrix(d11a2[[38]],ncol=J,nrow=J),matrix(d11a2[[39]],ncol=J,nrow=J),
matrix(d11a2[[40]],ncol=J,nrow=J),matrix(d11a2[[41]],ncol=J,nrow=J),matrix(d11a2[[42]],ncol=J,nrow=J),matrix(d11a2[[43]],ncol=J,nrow=J),matrix(d11a2[[44]],ncol=J,nrow=J))
d12<-do.call(adiag,d1temp2)


sig1<-vy*c+beta^2*eta*d1+beta^2*vx*d12

sig1i<-solve(sig1)
b1<-t(a1)%*%sig1i%*%a1
b1i<-solve(b1)
i1<-I
j1<-J
#############
# design 2 (I=22 different values, replicated J = 2 times)
I <-f[length(f)-1]
J<-n/I
th<-c(rep(a,J),rep(a+(b-a)/(I-1),J),rep(a+2*(b-a)/(I-1),J),rep(a+3*(b-a)/(I-1),J),rep(a+4*(b-a)/(I-1),J),rep(a+5*(b-a)/(I-1),J),rep(a+6*(b-a)/(I-1),J),rep(a+7*(b-a)/(I-1),J),rep(a+8*(b-a)/(I-1),J),
rep(a+9*(b-a)/(I-1),J),rep(a+10*(b-a)/(I-1),J),rep(a+11*(b-a)/(I-1),J),rep(a+12*(b-a)/(I-1),J),rep(a+13*(b-a)/(I-1),J),rep(a+14*(b-a)/(I-1),J),rep(a+15*(b-a)/(I-1),J),rep(a+16*(b-a)/(I-1),J),
rep(a+17*(b-a)/(I-1),J),rep(a+18*(b-a)/(I-1),J),rep(a+19*(b-a)/(I-1),J),rep(a+20*(b-a)/(I-1),J),rep(a+21*(b-a)/(I-1),J))
 
a2<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],
th[20],th[21],th[22],th[23],th[24],th[25],th[26],th[27],th[28],th[29],th[30],th[31],th[32],th[33],th[34],th[35],th[36],th[37],th[38],th[39],th[40],th[41],th[42],th[43],th[44]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+4*(b-a)/(I-1))^2,nrow=J,ncol=J),
matrix((a+5*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+6*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+7*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+8*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+9*(b-a)/(I-1))^2,nrow=J,ncol=J),
matrix((a+10*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+11*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+12*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+13*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+14*(b-a)/(I-1))^2,nrow=J,ncol=J),
matrix((a+15*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+16*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+17*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+18*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+19*(b-a)/(I-1))^2,nrow=J,ncol=J),
matrix((a+20*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+21*(b-a)/(I-1))^2,nrow=J,ncol=J))

d2temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J),matrix(d11a[[6]],ncol=J,nrow=J),
matrix(d11a[[7]],ncol=J,nrow=J),matrix(d11a[[8]],ncol=J,nrow=J),matrix(d11a[[9]],ncol=J,nrow=J),matrix(d11a[[10]],ncol=J,nrow=J),matrix(d11a[[11]],ncol=J,nrow=J),matrix(d11a[[12]],ncol=J,nrow=J),
matrix(d11a[[13]],ncol=J,nrow=J),matrix(d11a[[14]],ncol=J,nrow=J),matrix(d11a[[15]],ncol=J,nrow=J),matrix(d11a[[16]],ncol=J,nrow=J),matrix(d11a[[17]],ncol=J,nrow=J),matrix(d11a[[18]],ncol=J,nrow=J),
matrix(d11a[[19]],ncol=J,nrow=J),matrix(d11a[[20]],ncol=J,nrow=J),matrix(d11a[[21]],ncol=J,nrow=J),matrix(d11a[[22]],ncol=J,nrow=J))

d2<-do.call(adiag,d2temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J))

d2temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J),matrix(d11a2[[6]],ncol=J,nrow=J),
matrix(d11a2[[7]],ncol=J,nrow=J),matrix(d11a2[[8]],ncol=J,nrow=J),matrix(d11a2[[9]],ncol=J,nrow=J),matrix(d11a2[[10]],ncol=J,nrow=J),matrix(d11a2[[11]],ncol=J,nrow=J),matrix(d11a2[[12]],ncol=J,nrow=J),
matrix(d11a2[[13]],ncol=J,nrow=J),matrix(d11a2[[14]],ncol=J,nrow=J),matrix(d11a2[[15]],ncol=J,nrow=J),matrix(d11a2[[16]],ncol=J,nrow=J),matrix(d11a2[[17]],ncol=J,nrow=J),matrix(d11a2[[18]],ncol=J,nrow=J),
matrix(d11a2[[19]],ncol=J,nrow=J),matrix(d11a2[[20]],ncol=J,nrow=J),matrix(d11a2[[21]],ncol=J,nrow=J),matrix(d11a2[[22]],ncol=J,nrow=J))

d22<-do.call(adiag,d2temp2)


sig2<-vy*c+beta^2*eta*d2+beta^2*vx*d22
sig2i<-solve(sig2)
b2<-t(a2)%*%sig2i%*%a2
b2i<-solve(b2)
i2<-I
j2<-J
#############
# design 3 (I = 11 different values, replicated J=4 times)
I<-f[length(f)-2]
J<-n/I
th<-c(rep(a,J), rep(a+(b-a)/(I-1),J), rep(a+2*(b-a)/(I-1),J),rep(a+3*(b-a)/(I-1),J),rep(a+4*(b-a)/(I-1),J), rep(a+5*(b-a)/(I-1),J),rep(a+6*(b-a)/(I-1),J),rep(a+7*(b-a)/(I-1),J),rep(a+8*(b-a)/(I-1),J),
rep(a+9*(b-a)/(I-1),J),rep(a+10*(b-a)/(I-1),J))

a3<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],
th[20],th[21],th[22],th[23],th[24],th[25],th[26],th[27],th[28],th[29],th[30],th[31],th[32],th[33],th[34],th[35],th[36],th[37],th[38],th[39],th[40],th[41],th[42],th[43],th[44]), ncol=2)



d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+4*(b-a)/(I-1))^2,nrow=J,ncol=J),
matrix((a+5*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+6*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+7*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+8*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+9*(b-a)/(I-1))^2,nrow=J,ncol=J),
matrix((a+10*(b-a)/(I-1))^2,nrow=J,ncol=J))

d3temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J),matrix(d11a[[6]],ncol=J,nrow=J),
matrix(d11a[[7]],ncol=J,nrow=J),matrix(d11a[[8]],ncol=J,nrow=J),matrix(d11a[[9]],ncol=J,nrow=J),matrix(d11a[[10]],ncol=J,nrow=J),matrix(d11a[[11]],ncol=J,nrow=J))

d3<-do.call(adiag,d3temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d3temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J),
matrix(d11a2[[6]],ncol=J,nrow=J),matrix(d11a2[[7]],ncol=J,nrow=J),matrix(d11a2[[8]],ncol=J,nrow=J),matrix(d11a2[[9]],ncol=J,nrow=J),matrix(d11a2[[10]],ncol=J,nrow=J),
matrix(d11a2[[11]],ncol=J,nrow=J))

d32<-do.call(adiag,d3temp2)



sig3<-vy*c+beta^2*eta*d3+beta^2*vx*d32
sig3i<-solve(sig3)
b3<-t(a3)%*%sig3i%*%a3
b3i<-solve(b3)
i3<-I
j3<-J
#################
# design 4 (I = 4 different values, replicated J = 11 times)
I<-f[length(f)-3]
J<-n/I
th<-c(rep(a,J), rep(a+(b-a)/(I-1),J), rep(a+2*(b-a)/(I-1),J),rep(a+3*(b-a)/(I-1),J))
a4<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],
th[20],th[21],th[22],th[23],th[24],th[25],th[26],th[27],th[28],th[29],th[30],th[31],th[32],th[33],th[34],th[35],th[36],th[37],th[38],th[39],th[40],th[41],th[42],th[43],th[44]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(I-1))^2,nrow=J,ncol=J))

d4temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J))

d4<-do.call(adiag,d4temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d4temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J))

d42<-do.call(adiag,d4temp2)


sig4<-vy*c+beta^2*eta*d4+beta^2*vx*d42
sig4i<-solve(sig4)
b4<-t(a4)%*%sig4i%*%a4
b4i<-solve(b4)
i4<-I
j4<-J
j4<-J
#################
# design 5 (I = 2 different values, replicated J = 22 times)
I<-f[length(f)-4]
J<-n/I
th<-c(rep(a,J), rep(a+(b-a)/(I-1),J))

a5<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],
th[20],th[21],th[22],th[23],th[24],th[25],th[26],th[27],th[28],th[29],th[30],th[31],th[32],th[33],th[34],th[35],th[36],th[37],th[38],th[39],th[40],th[41],th[42],th[43],th[44]), ncol=2)


d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J))

d5temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J))

d5<-do.call(adiag,d5temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))
d5temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J))

d52<-do.call(adiag,d5temp2)


sig5<-vy*c+beta^2*eta*d5+beta^2*vx*d52
sig5i<-solve(sig5)
b5<-t(a5)%*%sig5i%*%a5
b5i<-solve(b5)
i5<-I
j5<-J
##################


deter<-c(det(b1i),det(b2i),det(b3i),det(b4i),det(b5i))
which(deter==min(deter))

eig1<-eigen(b1i)
tr1<-sum(eig1$values)
eig2<-eigen(b2i)
tr2<-sum(eig2$values)
eig3<-eigen(b3i)
tr3<-sum(eig3$values)
eig4<-eigen(b4i)
tr4<-sum(eig4$values)
eig5<-eigen(b5i)
tr5<-sum(eig5$values)
 
trace<-c(tr1,tr2,tr3,tr4,tr5)
best<-which(trace==min(trace))


vb<-c((b1i[[2,2]])^0.5,(b2i[[2,2]])^0.5,(b3i[[2,2]])^0.5,(b4i[[2,2]])^0.5,(b5i[[2,2]])^0.5)
bestb<-which(vb==min(vb))

cq<-diag(1,nq)
thq<-c(rep(beta,nq))
th2q<-matrix(c(1,xnew),nrow=1)
a5q<-matrix(c(rep(beta,nq)), ncol=1)
varab1<-th2q%*%b1i%*%t(th2q)
d11aq1<-list(matrix(varab1, nrow=nq,ncol=nq))
dq1<-do.call(adiag,d11aq1)
sigq1<-dq1
sigq1i<-solve(sigq1)
bq1<-t(a5q)%*%sigq1i%*%a5q 
stuq1<-(1/bq1+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b1i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq1<-stuq1/xnew*100

varab2<-th2q%*%b2i%*%t(th2q)
d11aq2<-list(matrix(varab2, nrow=nq,ncol=nq))
dq2<-do.call(adiag,d11aq2)
sigq2<-dq2
sigq2i<-solve(sigq2)
bq2<-t(a5q)%*%sigq2i%*%a5q 
stuq2<-(1/bq2+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b2i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq2<-stuq2/xnew*100

varab3<-th2q%*%b3i%*%t(th2q)
d11aq3<-list(matrix(varab3, nrow=nq,ncol=nq))
dq3<-do.call(adiag,d11aq3)
sigq3<-dq3
sigq3i<-solve(sigq3)
bq3<-t(a5q)%*%sigq3i%*%a5q 
stuq3<-(1/bq3+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b3i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq3<-stuq3/xnew*100

varab4<-th2q%*%b4i%*%t(th2q)
d11aq4<-list(matrix(varab4, nrow=nq,ncol=nq))
dq4<-do.call(adiag,d11aq4)
sigq4<-dq4
sigq4i<-solve(sigq4)
bq4<-t(a5q)%*%sigq4i%*%a5q 
stuq4<-(1/bq4+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b4i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq4<-stuq4/xnew*100

varab5<-th2q%*%b5i%*%t(th2q)
d11aq5<-list(matrix(varab5, nrow=nq,ncol=nq))
dq5<-do.call(adiag,d11aq5)
sigq5<-dq5
sigq5i<-solve(sigq5)
bq5<-t(a5q)%*%sigq5i%*%a5q 
stuq5<-(1/bq5+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b5i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq5<-stuq5/xnew*100



finalrelst<-c(relstuq1,relstuq2,relstuq3,relstuq4,relstuq5)
finbest<-which(finalrelst==min(finalrelst))

finali<-c(i1,i2,i3,i4,i5)
finalj<-c(j1,j2,j3,j4,j5)
jbest<-finalj[finbest]
ibest<-finali[finbest]
totrep<-n+ny
finstd<-finalrelst[finbest]
samplrep<-ny/r
result<-list(samplrep=samplrep,nosamples=r,finstd=finstd,ibest=ibest,jbest=jbest,totrep=totrep)

return(result)
}

#########################################################################################

#########################################################################################
############################################
caldesign148<-function(n,a,b,vy,vx,beta,eta,ny,vs,xnew,vi){
f<-divisors(n)
fy<-divisors(ny)
r<-1
if(vs>0 |vi>0){r=ny}
c<-diag(1,n)
v1<-(beta^2*eta*a^2)^0.5
v2<-(beta^2*eta*b^2)^0.5
nq<-1
#########
# #########
# design 1 (all I = 48 values are different, replicated J = 1 time)

J<-1
I<-f[length(f)]
th<-c(a,a+(b-a)/(I-1),a+2*(b-a)/(I-1),a+3*(b-a)/(I-1),a+4*(b-a)/(I-1),
                  a+5*(b-a)/(I-1),a+6*(b-a)/(I-1),a+7*(b-a)/(I-1),a+8*(b-a)/(I-1),a+9*(b-a)/(I-1),a+10*(b-a)/(I-1),a+11*(b-a)/(I-1),
a+12*(b-a)/(I-1),a+13*(b-a)/(I-1),a+14*(b-a)/(I-1),a+15*(b-a)/(I-1),
                  a+16*(b-a)/(I-1),a+17*(b-a)/(I-1),a+18*(b-a)/(I-1),a+19*(b-a)/(I-1),a+20*(b-a)/(I-1),a+21*(b-a)/(I-1),a+22*(b-a)/(I-1),
a+23*(b-a)/(I-1),a+24*(b-a)/(I-1),a+25*(b-a)/(I-1),a+26*(b-a)/(I-1),a+27*(b-a)/(I-1),a+28*(b-a)/(I-1),a+29*(b-a)/(I-1),a+30*(b-a)/(I-1),a+31*(b-a)/(I-1),a+32*(b-a)/(I-1),
a+33*(b-a)/(I-1),a+34*(b-a)/(I-1),a+35*(b-a)/(I-1),a+36*(b-a)/(I-1),a+37*(b-a)/(I-1),a+38*(b-a)/(I-1),a+39*(b-a)/(I-1),a+40*(b-a)/(I-1),a+41*(b-a)/(I-1),a+42*(b-a)/(I-1),a+43*(b-a)/(I-1),a+44*(b-a)/(I-1),
a+45*(b-a)/(I-1),a+46*(b-a)/(I-1),b)
 
a1<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],
th[20],th[21],th[22],th[23],th[24],th[25],th[26],th[27],th[28],th[29],th[30],th[31],th[32],th[33],th[34],th[35],th[36],th[37],th[38],th[39],th[40],th[41],th[42],th[43],th[44],th[45],th[46],
th[47],th[48]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+2*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+4*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+5*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+6*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+7*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+8*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+9*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+10*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+11*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+12*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+13*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+14*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+15*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+16*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+17*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+18*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+19*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+20*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+21*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+22*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+23*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+24*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+25*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+26*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+27*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+28*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+29*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+30*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+31*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+32*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+33*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+34*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+35*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+36*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+37*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+38*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+39*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+40*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+41*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+42*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+43*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+44*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+45*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),
matrix((a+46*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J),matrix((a+47*(b-a)/(f[length(f)-1]-1))^2,nrow=J,ncol=J))

d1temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J),matrix(d11a[[6]],ncol=J,nrow=J),matrix(d11a[[7]],ncol=J,nrow=J),matrix(d11a[[8]],ncol=J,nrow=J),matrix(d11a[[9]],ncol=J,nrow=J),
matrix(d11a[[10]],ncol=J,nrow=J),matrix(d11a[[11]],ncol=J,nrow=J),matrix(d11a[[12]],ncol=J,nrow=J),matrix(d11a[[13]],ncol=J,nrow=J),matrix(d11a[[14]],ncol=J,nrow=J),matrix(d11a[[15]],ncol=J,nrow=J),matrix(d11a[[16]],ncol=J,nrow=J),matrix(d11a[[17]],ncol=J,nrow=J),
matrix(d11a[[18]],ncol=J,nrow=J),matrix(d11a[[19]],ncol=J,nrow=J),matrix(d11a[[20]],ncol=J,nrow=J),matrix(d11a[[21]],ncol=J,nrow=J),
matrix(d11a[[22]],ncol=J,nrow=J),matrix(d11a[[23]],ncol=J,nrow=J),matrix(d11a[[24]],ncol=J,nrow=J),matrix(d11a[[25]],ncol=J,nrow=J),matrix(d11a[[26]],ncol=J,nrow=J),matrix(d11a[[27]],ncol=J,nrow=J),
matrix(d11a[[28]],ncol=J,nrow=J),matrix(d11a[[29]],ncol=J,nrow=J),matrix(d11a[[30]],ncol=J,nrow=J),matrix(d11a[[31]],ncol=J,nrow=J),matrix(d11a[[32]],ncol=J,nrow=J),matrix(d11a[[33]],ncol=J,nrow=J),
matrix(d11a[[34]],ncol=J,nrow=J),matrix(d11a[[35]],ncol=J,nrow=J),matrix(d11a[[36]],ncol=J,nrow=J),matrix(d11a[[37]],ncol=J,nrow=J),matrix(d11a[[38]],ncol=J,nrow=J),matrix(d11a[[39]],ncol=J,nrow=J),
matrix(d11a[[40]],ncol=J,nrow=J),matrix(d11a[[41]],ncol=J,nrow=J),matrix(d11a[[42]],ncol=J,nrow=J),matrix(d11a[[43]],ncol=J,nrow=J),matrix(d11a[[44]],ncol=J,nrow=J),matrix(d11a[[45]],ncol=J,nrow=J),
matrix(d11a[[46]],ncol=J,nrow=J),matrix(d11a[[47]],ncol=J,nrow=J),matrix(d11a[[48]],ncol=J,nrow=J))
d1<-do.call(adiag,d1temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d1temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J),matrix(d11a2[[6]],ncol=J,nrow=J),matrix(d11a2[[7]],ncol=J,nrow=J),matrix(d11a2[[8]],ncol=J,nrow=J),matrix(d11a2[[9]],ncol=J,nrow=J),
matrix(d11a2[[10]],ncol=J,nrow=J),matrix(d11a2[[11]],ncol=J,nrow=J),matrix(d11a2[[12]],ncol=J,nrow=J),matrix(d11a2[[13]],ncol=J,nrow=J),matrix(d11a2[[14]],ncol=J,nrow=J),matrix(d11a2[[15]],ncol=J,nrow=J),matrix(d11a2[[16]],ncol=J,nrow=J),matrix(d11a2[[17]],ncol=J,nrow=J),
matrix(d11a2[[18]],ncol=J,nrow=J),matrix(d11a2[[19]],ncol=J,nrow=J),matrix(d11a2[[20]],ncol=J,nrow=J),matrix(d11a2[[21]],ncol=J,nrow=J),
matrix(d11a2[[22]],ncol=J,nrow=J),matrix(d11a2[[23]],ncol=J,nrow=J),matrix(d11a2[[24]],ncol=J,nrow=J),matrix(d11a2[[25]],ncol=J,nrow=J),matrix(d11a2[[26]],ncol=J,nrow=J),matrix(d11a2[[27]],ncol=J,nrow=J),
matrix(d11a2[[28]],ncol=J,nrow=J),matrix(d11a2[[29]],ncol=J,nrow=J),matrix(d11a2[[30]],ncol=J,nrow=J),matrix(d11a2[[31]],ncol=J,nrow=J),matrix(d11a2[[32]],ncol=J,nrow=J),matrix(d11a2[[33]],ncol=J,nrow=J),
matrix(d11a2[[34]],ncol=J,nrow=J),matrix(d11a2[[35]],ncol=J,nrow=J),matrix(d11a2[[36]],ncol=J,nrow=J),matrix(d11a2[[37]],ncol=J,nrow=J),matrix(d11a2[[38]],ncol=J,nrow=J),matrix(d11a2[[39]],ncol=J,nrow=J),
matrix(d11a2[[40]],ncol=J,nrow=J),matrix(d11a2[[41]],ncol=J,nrow=J),matrix(d11a2[[42]],ncol=J,nrow=J),matrix(d11a2[[43]],ncol=J,nrow=J),matrix(d11a2[[44]],ncol=J,nrow=J),matrix(d11a2[[45]],ncol=J,nrow=J),
matrix(d11a2[[46]],ncol=J,nrow=J),matrix(d11a2[[47]],ncol=J,nrow=J),matrix(d11a2[[48]],ncol=J,nrow=J))
d12<-do.call(adiag,d1temp2)


sig1<-vy*c+beta^2*eta*d1+beta^2*vx*d12

sig1i<-solve(sig1)
b1<-t(a1)%*%sig1i%*%a1
b1i<-solve(b1)
i1<-I
j1<-J
#############
# design 2 (I=24 different values, replicated J = 2 times)
I <-f[length(f)-1]
J<-n/I
th<-c(rep(a,J),rep(a+(b-a)/(I-1),J),rep(a+2*(b-a)/(I-1),J),rep(a+3*(b-a)/(I-1),J),rep(a+4*(b-a)/(I-1),J),rep(a+5*(b-a)/(I-1),J),rep(a+6*(b-a)/(I-1),J),rep(a+7*(b-a)/(I-1),J),rep(a+8*(b-a)/(I-1),J),
rep(a+9*(b-a)/(I-1),J),rep(a+10*(b-a)/(I-1),J),rep(a+11*(b-a)/(I-1),J),rep(a+12*(b-a)/(I-1),J),rep(a+13*(b-a)/(I-1),J),rep(a+14*(b-a)/(I-1),J),rep(a+15*(b-a)/(I-1),J),rep(a+16*(b-a)/(I-1),J),
rep(a+17*(b-a)/(I-1),J),rep(a+18*(b-a)/(I-1),J),rep(a+19*(b-a)/(I-1),J),rep(a+20*(b-a)/(I-1),J),rep(a+21*(b-a)/(I-1),J),rep(a+22*(b-a)/(I-1),J),
rep(a+23*(b-a)/(I-1),J))
 
a2<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],
th[20],th[21],th[22],th[23],th[24],th[25],th[26],th[27],th[28],th[29],th[30],th[31],th[32],th[33],th[34],th[35],th[36],th[37],th[38],th[39],th[40],th[41],th[42],th[43],th[44],th[45],th[46],
th[47],th[48]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+4*(b-a)/(I-1))^2,nrow=J,ncol=J),
matrix((a+5*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+6*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+7*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+8*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+9*(b-a)/(I-1))^2,nrow=J,ncol=J),
matrix((a+10*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+11*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+12*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+13*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+14*(b-a)/(I-1))^2,nrow=J,ncol=J),
matrix((a+15*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+16*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+17*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+18*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+19*(b-a)/(I-1))^2,nrow=J,ncol=J),
matrix((a+20*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+21*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+22*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+23*(b-a)/(I-1))^2,nrow=J,ncol=J))

d2temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J),matrix(d11a[[6]],ncol=J,nrow=J),
matrix(d11a[[7]],ncol=J,nrow=J),matrix(d11a[[8]],ncol=J,nrow=J),matrix(d11a[[9]],ncol=J,nrow=J),matrix(d11a[[10]],ncol=J,nrow=J),matrix(d11a[[11]],ncol=J,nrow=J),matrix(d11a[[12]],ncol=J,nrow=J),
matrix(d11a[[13]],ncol=J,nrow=J),matrix(d11a[[14]],ncol=J,nrow=J),matrix(d11a[[15]],ncol=J,nrow=J),matrix(d11a[[16]],ncol=J,nrow=J),matrix(d11a[[17]],ncol=J,nrow=J),matrix(d11a[[18]],ncol=J,nrow=J),
matrix(d11a[[19]],ncol=J,nrow=J),matrix(d11a[[20]],ncol=J,nrow=J),matrix(d11a[[21]],ncol=J,nrow=J),matrix(d11a[[22]],ncol=J,nrow=J),matrix(d11a[[23]],ncol=J,nrow=J),matrix(d11a[[24]],ncol=J,nrow=J))

d2<-do.call(adiag,d2temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d2temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J),matrix(d11a2[[6]],ncol=J,nrow=J),
matrix(d11a2[[7]],ncol=J,nrow=J),matrix(d11a2[[8]],ncol=J,nrow=J),matrix(d11a2[[9]],ncol=J,nrow=J),matrix(d11a2[[10]],ncol=J,nrow=J),matrix(d11a2[[11]],ncol=J,nrow=J),matrix(d11a2[[12]],ncol=J,nrow=J),
matrix(d11a2[[13]],ncol=J,nrow=J),matrix(d11a2[[14]],ncol=J,nrow=J),matrix(d11a2[[15]],ncol=J,nrow=J),matrix(d11a2[[16]],ncol=J,nrow=J),matrix(d11a2[[17]],ncol=J,nrow=J),matrix(d11a2[[18]],ncol=J,nrow=J),
matrix(d11a2[[19]],ncol=J,nrow=J),matrix(d11a2[[20]],ncol=J,nrow=J),matrix(d11a2[[21]],ncol=J,nrow=J),matrix(d11a2[[22]],ncol=J,nrow=J),matrix(d11a2[[23]],ncol=J,nrow=J),matrix(d11a2[[24]],ncol=J,nrow=J))

d22<-do.call(adiag,d2temp2)


sig2<-vy*c+beta^2*eta*d2+beta^2*vx*d22
sig2i<-solve(sig2)
b2<-t(a2)%*%sig2i%*%a2
b2i<-solve(b2)
i2<-I
j2<-J
#############
# design 3 (I = 16 different values, replicated J=3 times)
I<-f[length(f)-2]
J<-n/I
th<-c(rep(a,J), rep(a+(b-a)/(I-1),J), rep(a+2*(b-a)/(I-1),J),rep(a+3*(b-a)/(I-1),J),rep(a+4*(b-a)/(I-1),J), rep(a+5*(b-a)/(I-1),J),rep(a+6*(b-a)/(I-1),J),rep(a+7*(b-a)/(I-1),J),rep(a+8*(b-a)/(I-1),J),
rep(a+9*(b-a)/(I-1),J),rep(a+10*(b-a)/(I-1),J),rep(a+11*(b-a)/(I-1),J),rep(a+12*(b-a)/(I-1),J),rep(a+13*(b-a)/(I-1),J),
rep(a+14*(b-a)/(I-1),J),rep(a+15*(b-a)/(I-1),J))

a3<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],
th[20],th[21],th[22],th[23],th[24],th[25],th[26],th[27],th[28],th[29],th[30],th[31],th[32],th[33],th[34],th[35],th[36],th[37],th[38],th[39],th[40],th[41],th[42],th[43],th[44],th[45],th[46],
th[47],th[48]), ncol=2)



d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+4*(b-a)/(I-1))^2,nrow=J,ncol=J),
matrix((a+5*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+6*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+7*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+8*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+9*(b-a)/(I-1))^2,nrow=J,ncol=J),
matrix((a+10*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+11*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+12*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+13*(b-a)/(I-1))^2,nrow=J,ncol=J),
matrix((a+14*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+15*(b-a)/(I-1))^2,nrow=J,ncol=J))

d3temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J),matrix(d11a[[6]],ncol=J,nrow=J),
matrix(d11a[[7]],ncol=J,nrow=J),matrix(d11a[[8]],ncol=J,nrow=J),matrix(d11a[[9]],ncol=J,nrow=J),matrix(d11a[[10]],ncol=J,nrow=J),matrix(d11a[[11]],ncol=J,nrow=J),matrix(d11a[[12]],ncol=J,nrow=J),
matrix(d11a[[13]],ncol=J,nrow=J),matrix(d11a[[14]],ncol=J,nrow=J),matrix(d11a[[15]],ncol=J,nrow=J),matrix(d11a[[16]],ncol=J,nrow=J))

d3<-do.call(adiag,d3temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J))

d3temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J),
matrix(d11a2[[6]],ncol=J,nrow=J),matrix(d11a2[[7]],ncol=J,nrow=J),matrix(d11a2[[8]],ncol=J,nrow=J),matrix(d11a2[[9]],ncol=J,nrow=J),matrix(d11a2[[10]],ncol=J,nrow=J),
matrix(d11a2[[11]],ncol=J,nrow=J),matrix(d11a2[[12]],ncol=J,nrow=J),matrix(d11a2[[13]],ncol=J,nrow=J),matrix(d11a2[[14]],ncol=J,nrow=J),
matrix(d11a2[[15]],ncol=J,nrow=J),matrix(d11a2[[16]],ncol=J,nrow=J))

d32<-do.call(adiag,d3temp2)



sig3<-vy*c+beta^2*eta*d3+beta^2*vx*d32
sig3i<-solve(sig3)
b3<-t(a3)%*%sig3i%*%a3
b3i<-solve(b3)
i3<-I
j3<-J
#################
# design 4 (I = 12 different values, replicated J = 4 times)
I<-f[length(f)-3]
J<-n/I
th<-c(rep(a,J), rep(a+(b-a)/(I-1),J), rep(a+2*(b-a)/(I-1),J),rep(a+3*(b-a)/(I-1),J),rep(a+4*(b-a)/(I-1),J), rep(a+5*(b-a)/(I-1),J),rep(a+6*(b-a)/(I-1),J), rep(a+7*(b-a)/(I-1),J),rep(a+8*(b-a)/(I-1),J),
rep(a+9*(b-a)/(I-1),J),rep(a+10*(b-a)/(I-1),J), rep(a+11*(b-a)/(I-1),J),rep(a+12*(b-a)/(I-1),J))
a4<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],
th[20],th[21],th[22],th[23],th[24],th[25],th[26],th[27],th[28],th[29],th[30],th[31],th[32],th[33],th[34],th[35],th[36],th[37],th[38],th[39],th[40],th[41],th[42],th[43],th[44],th[45],th[46],
th[47],th[48]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+4*(b-a)/(I-1))^2,nrow=J,ncol=J),
matrix((a+5*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+6*(b-a)/(I-1))^2,nrow=J,ncol=J),
matrix((a+7*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+8*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+9*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+10*(b-a)/(I-1))^2,nrow=J,ncol=J),
matrix((a+11*(b-a)/(I-1))^2,nrow=J,ncol=J))

d4temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J),matrix(d11a[[6]],ncol=J,nrow=J),
matrix(d11a[[7]],ncol=J,nrow=J),matrix(d11a[[8]],ncol=J,nrow=J),matrix(d11a[[9]],ncol=J,nrow=J),matrix(d11a[[10]],ncol=J,nrow=J),
matrix(d11a[[11]],ncol=J,nrow=J),matrix(d11a[[12]],ncol=J,nrow=J))

d4<-do.call(adiag,d4temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),
matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d4temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J),matrix(d11a2[[6]],ncol=J,nrow=J),
matrix(d11a2[[7]],ncol=J,nrow=J),matrix(d11a2[[8]],ncol=J,nrow=J),matrix(d11a2[[9]],ncol=J,nrow=J),matrix(d11a2[[10]],ncol=J,nrow=J),matrix(d11a2[[11]],ncol=J,nrow=J),matrix(d11a2[[12]],ncol=J,nrow=J))

d42<-do.call(adiag,d4temp2)


sig4<-vy*c+beta^2*eta*d4+beta^2*vx*d42
sig4i<-solve(sig4)
b4<-t(a4)%*%sig4i%*%a4
b4i<-solve(b4)
i4<-I
j4<-J
j4<-J
#################
# design 5 (I = 8 different values, replicated J = 6 times)
I<-f[length(f)-4]
J<-n/I
th<-c(rep(a,J), rep(a+(b-a)/(I-1),J), rep(a+2*(b-a)/(I-1),J),rep(a+3*(b-a)/(I-1),J),rep(a+4*(b-a)/(I-1),J),rep(a+5*(b-a)/(I-1),J),rep(a+6*(b-a)/(I-1),J),rep(a+7*(b-a)/(I-1),J))

a5<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],
th[20],th[21],th[22],th[23],th[24],th[25],th[26],th[27],th[28],th[29],th[30],th[31],th[32],th[33],th[34],th[35],th[36],th[37],th[38],th[39],th[40],th[41],th[42],th[43],th[44],th[45],th[46],
th[47],th[48]), ncol=2)


d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(I-1))^2,nrow=J,ncol=J),
matrix((a+4*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+5*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+6*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+7*(b-a)/(I-1))^2,nrow=J,ncol=J))

d5temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),
matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J),matrix(d11a[[6]],ncol=J,nrow=J),matrix(d11a[[7]],ncol=J,nrow=J),matrix(d11a[[8]],ncol=J,nrow=J))

d5<-do.call(adiag,d5temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d5temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),
matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J),matrix(d11a2[[6]],ncol=J,nrow=J),matrix(d11a2[[7]],ncol=J,nrow=J),matrix(d11a2[[8]],ncol=J,nrow=J))

d52<-do.call(adiag,d5temp2)


sig5<-vy*c+beta^2*eta*d5+beta^2*vx*d52
sig5i<-solve(sig5)
b5<-t(a5)%*%sig5i%*%a5
b5i<-solve(b5)
i5<-I
j5<-J
##################
# design 6 (I = 6 different values, replicated J = 8 times)
I<-f[length(f)-5]
J<-n/I

th<-c(rep(a,J), rep(a+(b-a)/(I-1),J), rep(a+2*(b-a)/(I-1),J),rep(a+3*(b-a)/(I-1),J),rep(a+4*(b-a)/(I-1),J),rep(a+5*(b-a)/(I-1),J))

a6<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],
th[20],th[21],th[22],th[23],th[24],th[25],th[26],th[27],th[28],th[29],th[30],th[31],th[32],th[33],th[34],th[35],th[36],th[37],th[38],th[39],th[40],th[41],th[42],th[43],th[44],th[45],th[46],
th[47],th[48]), ncol=2)

d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(I-1))^2,nrow=J,ncol=J),
matrix((a+3*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+4*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+5*(b-a)/(I-1))^2,nrow=J,ncol=J))

d6temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J),matrix(d11a[[5]],ncol=J,nrow=J),matrix(d11a[[6]],ncol=J,nrow=J))

d6<-do.call(adiag,d6temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d6temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J),matrix(d11a2[[5]],ncol=J,nrow=J),matrix(d11a2[[6]],ncol=J,nrow=J))

d62<-do.call(adiag,d6temp2)


sig6<-vy*c+beta^2*eta*d6+beta^2*vx*d62
sig6i<-solve(sig6)
b6<-t(a6)%*%sig6i%*%a6
b6i<-solve(b6)
i6<-I
j6<-J

##################
# design 7 (I = 4 different values, replicated J = 12 times)
I<-f[length(f)-6]
J<-n/I
th<-c(rep(a,J), rep(a+(b-a)/(I-1),J), rep(a+2*(b-a)/(I-1),J),rep(a+3*(b-a)/(I-1),J))

a7<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],
th[20],th[21],th[22],th[23],th[24],th[25],th[26],th[27],th[28],th[29],th[30],th[31],th[32],th[33],th[34],th[35],th[36],th[37],th[38],th[39],th[40],th[41],th[42],th[43],th[44],th[45],th[46],
th[47],th[48]), ncol=2)



d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+3*(b-a)/(I-1))^2,nrow=J,ncol=J))

d7temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J),matrix(d11a[[4]],ncol=J,nrow=J))

d7<-do.call(adiag,d7temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d7temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J),matrix(d11a2[[4]],ncol=J,nrow=J))

d72<-do.call(adiag,d7temp2)


sig7<-vy*c+beta^2*eta*d7+beta^2*vx*d72
sig7i<-solve(sig7)
b7<-t(a7)%*%sig7i%*%a7
b7i<-solve(b7)
i7<-I
j7<-J
#################
#design 8 (I = 3 different values, replicated J = 16 times)
I<-f[length(f)-7]
J<-n/I
th<-c(rep(a,J), rep(a+(b-a)/(I-1),J),rep(a+2*(b-a)/(I-1),J))

a8<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],
th[20],th[21],th[22],th[23],th[24],th[25],th[26],th[27],th[28],th[29],th[30],th[31],th[32],th[33],th[34],th[35],th[36],th[37],th[38],th[39],th[40],th[41],th[42],th[43],th[44],th[45],th[46],
th[47],th[48]), ncol=2)


d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J),matrix((a+2*(b-a)/(I-1))^2,nrow=J,ncol=J))

d8temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J),matrix(d11a[[3]],ncol=J,nrow=J))

d8<-do.call(adiag,d8temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d8temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J),matrix(d11a2[[3]],ncol=J,nrow=J))

d82<-do.call(adiag,d8temp2)


sig8<-vy*c+beta^2*eta*d8+beta^2*vx*d82
sig8i<-solve(sig8)
b8<-t(a8)%*%sig8i%*%a8
b8i<-solve(b8)
i8<-I
j8<-J
########################################
#design 9 (I=2 different values, replicated J=24 times)
I<-f[length(f)-8]
J<-n/I
th<-c(rep(a,J), rep(a+(b-a)/(I-1),J))

a9<-matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10],th[11],th[12],th[13],th[14],th[15],th[16],th[17],th[18],th[19],
th[20],th[21],th[22],th[23],th[24],th[25],th[26],th[27],th[28],th[29],th[30],th[31],th[32],th[33],th[34],th[35],th[36],th[37],th[38],th[39],th[40],th[41],th[42],th[43],th[44],th[45],th[46],
th[47],th[48]), ncol=2)



d11a<-list(matrix(a^2, nrow=J,ncol=J),matrix((a+(b-a)/(I-1))^2,nrow=J,ncol=J))

d9temp<-list(matrix(d11a[[1]],ncol=J,nrow=J),matrix(d11a[[2]],ncol=J,nrow=J))

d9<-do.call(adiag,d9temp)

d11a2<-list(matrix(nq, nrow=J,ncol=J),matrix(nq,nrow=J,ncol=J))

d9temp2<-list(matrix(d11a2[[1]],ncol=J,nrow=J),matrix(d11a2[[2]],ncol=J,nrow=J))

d92<-do.call(adiag,d9temp2)


sig9<-vy*c+beta^2*eta*d9+beta^2*vx*d92
sig9i<-solve(sig9)
b9<-t(a9)%*%sig9i%*%a9
b9i<-solve(b9)
i9<-I
j9<-J

deter<-c(det(b1i),det(b2i),det(b3i),det(b4i),det(b5i),det(b6i),det(b7i),det(b8i),det(b9i))
which(deter==min(deter))

eig1<-eigen(b1i)
tr1<-sum(eig1$values)
eig2<-eigen(b2i)
tr2<-sum(eig2$values)
eig3<-eigen(b3i)
tr3<-sum(eig3$values)
eig4<-eigen(b4i)
tr4<-sum(eig4$values)
eig5<-eigen(b5i)
tr5<-sum(eig5$values)
eig6<-eigen(b6i)
tr6<-sum(eig6$values)
eig7<-eigen(b7i)
tr7<-sum(eig7$values)
eig8<-eigen(b8i)
tr8<-sum(eig8$values)
eig9<-eigen(b9i)
tr9<-sum(eig9$values)
 
trace<-c(tr1,tr2,tr3,tr4,tr5,tr6,tr7,tr8,tr9)
best<-which(trace==min(trace))


vb<-c((b1i[[2,2]])^0.5,(b2i[[2,2]])^0.5,(b3i[[2,2]])^0.5,(b4i[[2,2]])^0.5,(b5i[[2,2]])^0.5,(b6i[[2,2]])^0.5,(b7i[[2,2]])^0.5,(b8i[[2,2]])^0.5,(b9i[[2,2]])^0.5)
bestb<-which(vb==min(vb))

cq<-diag(1,nq)
thq<-c(rep(beta,nq))
th2q<-matrix(c(1,xnew),nrow=1)
a5q<-matrix(c(rep(beta,nq)), ncol=1)
varab1<-th2q%*%b1i%*%t(th2q)
d11aq1<-list(matrix(varab1, nrow=nq,ncol=nq))
dq1<-do.call(adiag,d11aq1)
sigq1<-dq1
sigq1i<-solve(sigq1)
bq1<-t(a5q)%*%sigq1i%*%a5q 
stuq1<-(1/bq1+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b1i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq1<-stuq1/xnew*100

varab2<-th2q%*%b2i%*%t(th2q)
d11aq2<-list(matrix(varab2, nrow=nq,ncol=nq))
dq2<-do.call(adiag,d11aq2)
sigq2<-dq2
sigq2i<-solve(sigq2)
bq2<-t(a5q)%*%sigq2i%*%a5q 
stuq2<-(1/bq2+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b2i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq2<-stuq2/xnew*100

varab3<-th2q%*%b3i%*%t(th2q)
d11aq3<-list(matrix(varab3, nrow=nq,ncol=nq))
dq3<-do.call(adiag,d11aq3)
sigq3<-dq3
sigq3i<-solve(sigq3)
bq3<-t(a5q)%*%sigq3i%*%a5q 
stuq3<-(1/bq3+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b3i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq3<-stuq3/xnew*100

varab4<-th2q%*%b4i%*%t(th2q)
d11aq4<-list(matrix(varab4, nrow=nq,ncol=nq))
dq4<-do.call(adiag,d11aq4)
sigq4<-dq4
sigq4i<-solve(sigq4)
bq4<-t(a5q)%*%sigq4i%*%a5q 
stuq4<-(1/bq4+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b4i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq4<-stuq4/xnew*100

varab5<-th2q%*%b5i%*%t(th2q)
d11aq5<-list(matrix(varab5, nrow=nq,ncol=nq))
dq5<-do.call(adiag,d11aq5)
sigq5<-dq5
sigq5i<-solve(sigq5)
bq5<-t(a5q)%*%sigq5i%*%a5q 
stuq5<-(1/bq5+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b5i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq5<-stuq5/xnew*100

varab6<-th2q%*%b6i%*%t(th2q)
d11aq6<-list(matrix(varab6, nrow=nq,ncol=nq))
dq6<-do.call(adiag,d11aq6)
sigq6<-dq6
sigq6i<-solve(sigq6)
bq6<-t(a5q)%*%sigq6i%*%a5q 
stuq6<-(1/bq6+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b6i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq6<-stuq6/xnew*100

varab7<-th2q%*%b7i%*%t(th2q)
d11aq7<-list(matrix(varab7, nrow=nq,ncol=nq))
dq7<-do.call(adiag,d11aq7)
sigq7<-dq7
sigq7i<-solve(sigq7)
bq7<-t(a5q)%*%sigq7i%*%a5q 
stuq7<-(1/bq7+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b7i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq7<-stuq7/xnew*100

varab8<-th2q%*%b8i%*%t(th2q)
d11aq8<-list(matrix(varab8, nrow=nq,ncol=nq))
dq8<-do.call(adiag,d11aq8)
sigq8<-dq8
sigq8i<-solve(sigq8)
bq8<-t(a5q)%*%sigq8i%*%a5q 
stuq8<-(1/bq8+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b8i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq8<-stuq8/xnew*100

varab9<-th2q%*%b9i%*%t(th2q)
d11aq9<-list(matrix(varab9, nrow=nq,ncol=nq))
dq9<-do.call(adiag,d11aq9)
sigq9<-dq9
sigq9i<-solve(sigq9)
bq9<-t(a5q)%*%sigq9i%*%a5q 
stuq9<-(1/bq9+vy/ny/beta/beta+vs/r/beta/beta+(beta^2+b9i[2,2])*xnew^2*vi/r/beta^2)^0.5
relstuq9<-stuq9/xnew*100

finalrelst<-c(relstuq1,relstuq2,relstuq3,relstuq4,relstuq5,relstuq6,relstuq7,relstuq8,relstuq9)
finbest<-which(finalrelst==min(finalrelst))

finali<-c(i1,i2,i3,i4,i5,i6,i7,i8,i9)
finalj<-c(j1,j2,j3,j4,j5,j6,j7,j8,j9)
jbest<-finalj[finbest]
ibest<-finali[finbest]
totrep<-n+ny
finstd<-finalrelst[finbest]
samplrep<-ny/r
result<-list(samplrep=samplrep,nosamples=r,finstd=finstd,ibest=ibest,jbest=jbest,totrep=totrep)

return(result)
}

#########################################################################################
#########################################################################################
#########################################################################################

optx<-NULL
A<-NULL
quantn<-totn-4
if(quantn>0){
for(i in 1:quantn){
des104<-caldesign104(4,a,b,vy,vx,beta,eta,i,vs,xnew,vi)
if(caln>3){A<-rbind(A,c(des104$ibest,des104$jbest,des104$nosamples,des104$samplrep,des104$finstd))}
}}
quantn<-totn-6
if(quantn>0){
for(i in 1:quantn){
des106<-caldesign106(6,a,b,vy,vx,beta,eta,i,vs,xnew,vi)
if(caln>5){A<-rbind(A,c(des106$ibest,des106$jbest,des106$nosamples,des106$samplrep,des106$finstd))}
}}
quantn<-totn-8
if(quantn>0){
for(i in 1:quantn){
des108<-caldesign108(8,a,b,vy,vx,beta,eta,i,vs,xnew,vi)
if(caln>7){A<-rbind(A,c(des108$ibest,des108$jbest,des108$nosamples,des108$samplrep,des108$finstd))}
}}
quantn<-totn-10
if(quantn>0){
for(i in 1:quantn){
des110<-caldesign110(10,a,b,vy,vx,beta,eta,i,vs,xnew,vi)
if(caln>9){A<-rbind(A,c(des110$ibest,des110$jbest,des110$nosamples,des110$samplrep,des110$finstd))}
}}
quantn<-totn-12
if(quantn>0){
for(i in 1:quantn){
des112<-caldesign112(12,a,b,vy,vx,beta,eta,i,vs,xnew,vi)
if(caln>11){A<-rbind(A,c(des112$ibest,des112$jbest,des112$nosamples,des112$samplrep,des112$finstd))}
}}

quantn<-totn-14
if(quantn>0){
for(i in 1:quantn){
des114<-caldesign114(14,a,b,vy,vx,beta,eta,i,vs,xnew,vi)
if(caln>13){A<-rbind(A,c(des114$ibest,des114$jbest,des114$nosamples,des114$samplrep,des114$finstd))}
}}

quantn<-totn-15
if(quantn>0){
for(i in 1:quantn){
des115<-caldesign115(15,a,b,vy,vx,beta,eta,i,vs,xnew,vi)
if(caln>14){A<-rbind(A,c(des115$ibest,des115$jbest,des115$nosamples,des115$samplrep,des115$finstd))}
}}
quantn<-totn-16
if(quantn>0){
for(i in 1:quantn){
des116<-caldesign116(16,a,b,vy,vx,beta,eta,i,vs,xnew,vi)
if(caln>15){A<-rbind(A,c(des116$ibest,des116$jbest,des116$nosamples,des116$samplrep,des116$finstd))}
}}
quantn<-totn-18
if(quantn>0){
for(i in 1:quantn){
des118<-caldesign118(18,a,b,vy,vx,beta,eta,i,vs,xnew,vi)
if(caln>17){A<-rbind(A,c(des118$ibest,des118$jbest,des118$nosamples,des118$samplrep,des118$finstd))}
}}
quantn<-totn-20
if(quantn>0){
for(i in 1:quantn){
des120<-caldesign120(20,a,b,vy,vx,beta,eta,i,vs,xnew,vi)
if(caln>19){A<-rbind(A,c(des120$ibest,des120$jbest,des120$nosamples,des120$samplrep,des120$finstd))}
}}

quantn<-totn-22
if(quantn>0){
for(i in 1:quantn){
des122<-caldesign122(22,a,b,vy,vx,beta,eta,i,vs,xnew,vi)
if(caln>21){A<-rbind(A,c(des122$ibest,des122$jbest,des122$nosamples,des122$samplrep,des122$finstd))}
}}

quantn<-totn-24
if(quantn>0){
for(i in 1:quantn){
des124<-caldesign124(24,a,b,vy,vx,beta,eta,i,vs,xnew,vi)

if(caln>23){A<-rbind(A,c(des124$ibest,des124$jbest,des124$nosamples,des124$samplrep,des124$finstd))}
}}

quantn<-totn-28
if(quantn>0){
for(i in 1:quantn){
des128<-caldesign128(28,a,b,vy,vx,beta,eta,i,vs,xnew,vi)

if(caln>27){A<-rbind(A,c(des128$ibest,des128$jbest,des128$nosamples,des128$samplrep,des128$finstd))}
}}

quantn<-totn-30
if(quantn>0){
for(i in 1:quantn){
des130<-caldesign130(30,a,b,vy,vx,beta,eta,i,vs,xnew,vi)

if(caln>29){A<-rbind(A,c(des130$ibest,des130$jbest,des130$nosamples,des130$samplrep,des130$finstd))}
}}

quantn<-totn-36
if(quantn>0){
for(i in 1:quantn){
des136<-caldesign136(36,a,b,vy,vx,beta,eta,i,vs,xnew,vi)

if(caln>35){A<-rbind(A,c(des136$ibest,des136$jbest,des136$nosamples,des136$samplrep,des136$finstd))}
}}
 
quantn<-totn-40
if(quantn>0){
for(i in 1:quantn){
des140<-caldesign140(40,a,b,vy,vx,beta,eta,i,vs,xnew,vi)

if(caln>39){A<-rbind(A,c(des140$ibest,des140$jbest,des140$nosamples,des140$samplrep,des140$finstd))}
}}

quantn<-totn-44
if(quantn>0){
for(i in 1:quantn){
des144<-caldesign144(44,a,b,vy,vx,beta,eta,i,vs,xnew,vi)

if(caln>43){A<-rbind(A,c(des144$ibest,des144$jbest,des144$nosamples,des144$samplrep,des144$finstd))}
}}


quantn<-totn-48
if(quantn>0){
for(i in 1:quantn){
des148<-caldesign148(48,a,b,vy,vx,beta,eta,i,vs,xnew,vi)

if(caln>47){A<-rbind(A,c(des148$ibest,des148$jbest,des148$nosamples,des148$samplrep,des148$finstd))}
}}

rowmin<-apply( A, 2, which.min)
optI<-A[rowmin[5],1]
if(optI==2){optx[1]<-a
optx[2]<-b}

if(optI==3){optx[1]<-a
optx[2]<-(a+(b-a))/2
optx[3]<-b}

if(optI==4){optx[1]<-a
optx[2]<-(a+(b-a))/3
optx[3]<-(a+2*(b-a))/3
optx[4]<-b}

if(optI==5){optx[1]<-a
for(i in 1:4){
optx[i+1]<-a+i*(b-a)/4
}}

if(optI==6){optx[1]<-a
for(i in 1:5){
optx[i+1]<-a+i*(b-a)/5
}}

if(optI==7){optx[1]<-a
for(i in 1:6){
optx[i+1]<-a+i*(b-a)/6
}}

if(optI==8){optx[1]<-a
for(i in 1:7){
optx[i+1]<-a+i*(b-a)/7
}}


if(optI==9){optx[1]<-a
for(i in 1:8){
optx[i+1]<-a+i*(b-a)/8
}}

if(optI==10){optx[1]<-a
for(i in 1:9){
optx[i+1]<-a+i*(b-a)/9
}}

if(optI==11){optx[1]<-a
for(i in 1:10){
optx[i+1]<-a+i*(b-a)/10
}}

if(optI==12){optx[1]<-a
for(i in 1:11){
optx[i+1]<-a+i*(b-a)/11
}}

if(optI==13){optx[1]<-a
for(i in 1:12){
optx[i+1]<-a+i*(b-a)/12
}}

if(optI==14){optx[1]<-a
for(i in 1:13){
optx[i+1]<-a+i*(b-a)/13
}}

if(optI==15){optx[1]<-a
for(i in 1:14){
optx[i+1]<-a+i*(b-a)/14
}}

if(optI==16){optx[1]<-a
for(i in 1:15){
optx[i+1]<-a+i*(b-a)/15
}}

if(optI==18){optx[1]<-a
for(i in 1:17){
optx[i+1]<-a+i*(b-a)/17
}}

if(optI==20){optx[1]<-a
for(i in 1:19){
optx[i+1]<-a+i*(b-a)/19
}}
 
if(optI==22){optx[1]<-a
for(i in 1:21){
optx[i+1]<-a+i*(b-a)/21
}}

if(optI==24){optx[1]<-a
for(i in 1:23){
optx[i+1]<-a+i*(b-a)/23
}}

if(optI==28){optx[1]<-a
for(i in 1:27){
optx[i+1]<-a+i*(b-a)/27
}}

if(optI==30){optx[1]<-a
for(i in 1:29){
optx[i+1]<-a+i*(b-a)/29
}}

if(optI==36){optx[1]<-a
for(i in 1:35){
optx[i+1]<-a+i*(b-a)/35
}}

if(optI==40){optx[1]<-a
for(i in 1:39){
optx[i+1]<-a+i*(b-a)/39
}}

if(optI==44){optx[1]<-a
for(i in 1:43){
optx[i+1]<-a+i*(b-a)/43
}}

if(optI==48){optx[1]<-a
for(i in 1:47){
optx[i+1]<-a+i*(b-a)/47
}}

optJ<-A[rowmin[5],2]
optr<-A[rowmin[5],3]
optnr<-A[rowmin[5],4]
 
optstd<-A[rowmin[5],5]
result2<-list(optI=optI,optJ=optJ,optx=unlist(optx),optr=optr,optnr=optnr,optstd=optstd)
return(result2)}


