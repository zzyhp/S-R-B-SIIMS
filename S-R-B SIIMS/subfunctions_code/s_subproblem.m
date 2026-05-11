function ss=s_subproblem(eta3,D,g,kappa0,kappa)
%TTk:TT(x+v^{k+1})
%f=T;
ebslong=0.01;sigma=1;
[N1,N2]=size(D);
Nt=5;
s=zeros(N1,N2,Nt+1);
s(1:N1,1:N2,1)=zeros(N1,N2);
ss=zeros(N1,N2);
dt=1;



for n=1:Nt



for i=2:N1-1
for j=2:N2-1
DXP=s(i+1,j,n)-s(i,j,n);
DXM=s(i,j,n)-s(i-1,j,n);
DX0=0.5*(s(i+1,j,n)-s(i-1,j,n));
DYP=s(i,j+1,n)-s(i,j,n);
DYM=s(i,j,n)-s(i,j-1,n);
DY0=0.5*(s(i,j+1,n)-s(i,j-1,n));
CE=1/sqrt(ebslong^2+DXP^2+DY0^2);
CW=1/sqrt(ebslong^2+DXM^2+DY0^2);
CS=1/sqrt(ebslong^2+DX0^2+DYP^2);
CN=1/sqrt(ebslong^2+DX0^2+DYM^2);
FM=1+eta3*dt*(CE+CW+CS+CN);
FZ=s(i,j,n)+eta3*dt*(CE*s(i+1,j,n)+CW*s(i-1,j,n)+CS*s(i,j+1,n)+CN*s(i,j-1,n))+dt*g(i,j)+dt*sigma*ss(i,j);
s(i,j,n+1)=FZ/FM;
if s(i,j,n+1)<=kappa0
    s(i,j,n+1)=kappa0;
elseif s(i,j,n+1)>=kappa
    s(i,j,n+1)=kappa;
end
end
end

end
ss(1:N1,1:N2)=s(1:N1,1:N2,Nt);
end



