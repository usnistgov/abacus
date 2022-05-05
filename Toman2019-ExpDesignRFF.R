### This function calculates the optimal design for the RRF experiment 
 
expdesign2<-function(caln,totn,theta,b,bsigh,vyh,vxh,lxh,vsh,xnew,lnewh){
 

caldes<-function(n,ny,theta,b,bsigh,vyh,vxh,lxh,vsh,xnew,lnewh){
r<-1
if(vsh>0 |lnewh>0){r=ny}
vbhat<-1/theta^2/n*(bsigh*theta^2+vyh+(lxh*theta^2+vxh)*(b^2+bsigh))
sy<-(xnew^2*vbhat/b^2+vyh/ny/b^2+vsh/r/b^2+(b^2+vbhat)*xnew^2*lnewh/r/b^2)^0.5
relsq<-sy/xnew*100
totrep<-n+ny
samplrep<-ny/r
result<-list(samplrep=samplrep,nosamples=r,relsq=relsq,totrep=totrep)
return(result)
}


optx<-NULL
A<-NULL
for(j in 4:40){
quantn<-totn-j
if(quantn>0){
for(i in 1:quantn){
des<-caldes(j,i,theta,b,bsigh,vyh,vxh,lxh,vsh,xnew,lnewh)
if(caln>j-1){A<-rbind(A,c(j,des$relsq,des$nosamples,des$samplrep,des$totrep))}
}}}

rowmin<-apply( A, 2, which.min)
optI<-A[rowmin[2],1]
opts<-A[rowmin[2],3]
optq<-A[rowmin[2],4]
optstd<-A[rowmin[2],2]

result2<-list(optI=optI,opts=opts,optq=optq,optstd=optstd)
return(result2)}
