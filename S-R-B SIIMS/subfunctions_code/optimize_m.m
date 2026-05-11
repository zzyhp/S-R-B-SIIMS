function m=optimize_m(beta2,zeta,g,eta1,ITE)
N=128;
dt=0.1;


m=zeros(N+1);
s=zeros(N+1);
p1=zeros(N+1);
p2=zeros(N+1);

n=1;
while n<=ITE 

div_p=divergence(p1,p2);


for i=2:N 
for j=2:N 
m(i,j)=g(i,j)-(eta1/2*(beta2+zeta))*div_p(i,j);
end
end


for i=2:N 
for j=2:N 
s(i,j)=(eta1^2/(4*(beta2+zeta)))*div_p(i,j)-(eta1^2/2)*g(i,j);
end
end

[s_x,s_y]=gradient(s);


for i=2:N 
for j=2:N 
p1(i,j)=(p1(i,j)+dt*s_x(i,j))/(1+dt*sqrt(s_x(i,j).^2+s_y(i,j).^2));
p2(i,j)=(p2(i,j)+dt*s_y(i,j))/(1+dt*sqrt(s_x(i,j).^2+s_y(i,j).^2));
end
end


n=n+1;
end

%result
div_p=divergence(p1,p2);

for i=2:N 
for j=2:N 
    m(i,j)=g(i,j)-(eta1/2*(beta2+zeta))*div_p(i,j);
end
end
end