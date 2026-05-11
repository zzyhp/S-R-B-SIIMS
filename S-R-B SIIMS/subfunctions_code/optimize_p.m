function ps = optimize_p(lambda1,p,p1,beta1,gs,kexi,q)
mu3=2*lambda1;
[M,N]=size(p);
epslong=0.000001;
phy_x=zeros(129,129);
phy_y=zeros(129,129);
for i=1:129
for j=1:129
phy_x(i,j)=i;phy_y(i,j)=j;
end
end
    d = zeros(M, N, 2);  
    b = zeros(M, N, 2); 
   term=zeros(129,129);

    for n=1:3


         % Step 1: v = d - b
        v1 = d(:,:,1) - b(:,:,1);
        v2 = d(:,:,2) - b(:,:,2);
        fi=divergence(v1,v2);

 for i=2:128
     for j=2:128
         p(i,j)=(mu3*(p(i+1,j)+p(i-1,j)+p(i,j+1)+p(i,j-1))-mu3*fi(i,j)-(beta1*gs(i,j))/2+2*kexi*(q(i,j)-p1(i,j)))/4*mu3;
     end
 end

[ppx,ppy]=gradient(p);
    for i=1:129
        for j=1:129
            term(i,j)=sqrt((ppx(i,j)+b(i,j,1))^2+(ppy(i,j)+b(i,j,2))^2);
        end
    end

    for i=1:129
        for j=1:129
    if term(i,j)-lambda1/mu3>0
        term(i,j)=term(i,j)-lambda1/mu3;
    else
        term(i,j)=0;
    end
        end
    end

    for i=1:129
        for j=1:129
            d(i,j,1)=(ppx(i,j)+b(i,j,1))/(sqrt((ppx(i,j)+b(i,j,1))^2+(ppy(i,j)+b(i,j,2))^2)+epslong^2)*term(i,j);
            d(i,j,2)=(ppy(i,j)+b(i,j,2))/(sqrt((ppx(i,j)+b(i,j,1))^2+(ppy(i,j)+b(i,j,2))^2)+epslong^2)*term(i,j);
        end
    end

    for i=1:129
        for j=1:129
            b(i,j,1)=b(i,j,1)+ppx(i,j)-d(i,j,1);
            b(i,j,2)=b(i,j,2)+ppy(i,j)-d(i,j,2);
        end
    end

    end

ps=p;
end