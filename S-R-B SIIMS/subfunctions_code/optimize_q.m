function qs = optimize_q(p,lambda2,beta2,gs,kexi,q)
mu3=2*lambda2;
[M,N]=size(q);
epslong=0.0001;
    d = zeros(M, N, 2); 
    b = zeros(M, N, 2); 
term=zeros(129,129);


    for n=1:3
        v1 = d(:,:,1) - b(:,:,1);
        v2 = d(:,:,2) - b(:,:,2);
        fi=divergence(v1,v2);

 for i=2:128
     for j=2:128
         q(i,j)=(mu3*(q(i+1,j)+q(i-1,j)+q(i,j+1)+q(i,j-1))-mu3*fi(i,j)-(beta2*gs(i,j))/2+2*kexi*p(i,j))/(4*mu3+2*kexi);
     end
 end

[qqx,qqy]=gradient(q);
    for i=1:129
        for j=1:129
            term(i,j)=sqrt((qqx(i,j)+b(i,j,1))^2+(qqy(i,j)+b(i,j,2))^2);
        end
    end

    for i=1:129
        for j=1:129
    if term(i,j)-lambda2/mu3>0
        term(i,j)=term(i,j)-lambda2/mu3;
    else
        term(i,j)=0;
    end
        end
    end



    for i=1:129
        for j=1:129
            d(i,j,1)=(qqx(i,j)+b(i,j,1))/(sqrt((qqx(i,j)+b(i,j,1))^2+(qqy(i,j)+b(i,j,2))^2)+epslong^2)*term(i,j);
            d(i,j,2)=(qqy(i,j)+b(i,j,2))/(sqrt((qqx(i,j)+b(i,j,1))^2+(qqy(i,j)+b(i,j,2))^2)+epslong^2)*term(i,j);
        end
    end

    for i=1:129
        for j=1:129
            b(i,j,1)=b(i,j,1)+qqx(i,j)-d(i,j,1);
            b(i,j,2)=b(i,j,2)+qqy(i,j)-d(i,j,2);
        end
    end

    end

qs=q;
end