% Huan Han, Yuanhao Zha, Daoping Zhang and Yimin Zhang. A novel framework for
% joint segmentation-registration-bias correction in low-contrast images.
% SIAM Journal on Imaging Sciences.

% The Code is partially created based on the method described in the following paper 
%   [1] Peng Chen, Ke Chen, Huan Han*,Daoping Zhang.Multiscale approach for variational problem joint diffeomorphic
%       image registration and intensity correction: theory and application.
%       SIAM Journal on Multiscale Modeling and Simulation,accepted.
%   [2] Huan Han, Zhengping Wang, Yimin Zhang. Multiscale approach for Two-Dimensional diffeomorphic image registration. 
%      SIAM Journal on Multiscale Modeling and Simulation,19(4) (2021)1538-1572.
%   [3] Huan Han, Andong Wang. A fast multi grid algorithm for 2D diffeomorphic image 
%       registration model. Journal of Computational and Applied Mathematics, 394 (2021) :113576
%   [4] Huan Han, Zhengping Wang.A diffeomorphic image registration model with fractional-order regularization
%       and Cauchy-Riemann constraint.SIAM Journal on Imaging Sciences 13(3)(2020) 1240-1271.
clc; clear; close all; eval(['help ' mfilename])
addpath(genpath('data'));
addpath(genpath('subfunctions_code'));
N1=129; N2=129; 
data=load('Brain_M3.mat');

%%%%%%%%%%%%%%%%%%%%parameter%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K=3;beta2=0.2;K1=0;K2=5;zeta=3;kexi=0.05;gamma=0.1;ITE=800;mu2=0.00000000000001;beta1=0.15 ;mu1=0.0000000000000;
T=data.T;D=data.D;kappa0=-0.001;kappa=0.001;lambda1=0.15;lambda2=0.01;eta1=0.01;eta2=0.00000001;
%%%%%%%%%%%%%%%%%%%%parameter%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%processing%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmin=min(min(T));dmin=min(min(D));
for i=1:129
for j=1:129
T(i,j)=T(i,j)-tmin+0.01;D(i,j)=D(i,j)-dmin+0.01;
end
end 
%%%%%%%%%%%%%%%%%%%%processing%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%Initial value%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[idx, C1] = kmeans1(T, K);
[idx1, C2] = kmeans1(T, K);
c=zeros(K,1);cc=zeros(K,1);
for l=1:K
    c(l,1)=log(C1(l,1));
    cc(l,1)=log(C2(l,1));
end
ccc=c;
u1=zeros(N1,N2);
u2=zeros(N1,N2);
phy_x=zeros(N1,N2);
phy_y=zeros(N1,N2);
filter1=fspecial('gaussian',3,1000000.1);
filter2=fspecial('gaussian',3,1.1);
qi=(1/K)*ones(N1,N2,K);
p_i=(1/K)*ones(N1,N2,K);
qqi=(1/K)*ones(N1,N2,K);

ppi=p_i;
s=zeros(N1,N2);
TT=T;TTk=T;
m=zeros(N1,N2);

for i=1:129
for j=1:129
phy_x(i,j)=i;phy_y(i,j)=j;
end
end
%%%%%%%%%%%%%%%%%%%%Initial value%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
%%%%%%%%%%%%%%%%%%%%AMM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=1;
 while n<ITE

%%%%%%%%%%u-subproblem%%%%%%%%%%
ss1_1=zeros(N1,N2);ss1_2=zeros(N1,N2);lnTTk=zeros(N1,N2);
for i=2:128
    for j=2:128
        lnTTk(i,j)=log(abs(TTk(i,j))+0.000001);
    end
end
[lnTTkx,lnTTky]=gradient(lnTTk);
[TTkx,TTky]=gradient(TTk);
[ppix,ppiy]=gradient(ppi);ssCE_1=zeros(N1,N2);ssCE_2=zeros(N1,N2);
ss2_1=zeros(N1,N2);ss2_2=zeros(N1,N2);
f1=zeros(N1,N2);f2=zeros(N1,N2);

for l=1:K
for i=2:128
    for j=2:128
        ss1_1(i,j)=ss1_1(i,j)+(qi(i,j,l)*(log(abs(TTk(i,j)-s(i,j))+0.0001)-m(i,j)-cc(l,1))*lnTTkx(i,j));
        ss1_2(i,j)=ss1_2(i,j)+(qi(i,j,l)*(log(abs(TTk(i,j)-s(i,j))+0.0001)-m(i,j)-cc(l,1))*lnTTky(i,j));
    end
end
end

for l=1:K
    for i=2:128
        for j=2:128
            ssCE_1(i,j)=ssCE_1(i,j)+(ppi(i,j,l)-qi(i,j,l))*ppix(i,j,l);
            ssCE_2(i,j)=ssCE_2(i,j)+(ppi(i,j,l)-qi(i,j,l))*ppiy(i,j,l);
        end
    end
end

for i=2:128
    for j=2:128
        ss2_1(i,j)=-(m(i,j)+log(D(i,j))-log(abs(TTk(i,j)-s(i,j))+0.000001))*lnTTkx(i,j);
        ss2_2(i,j)=-(m(i,j)+log(D(i,j))-log(abs(TTk(i,j)-s(i,j))+0.000001))*lnTTky(i,j);
    end
end

for i=2:128
    for j=2:128
        f1(i,j)=beta2*ss1_1(i,j)+2*kexi*ssCE_1(i,j)+zeta*ss2_1(i,j);
        f2(i,j)=beta2*ss1_2(i,j)+2*kexi*ssCE_2(i,j)+zeta*ss2_2(i,j);
    end
end
uu1=MCG2D(5,u1-zeta*ss2_1,u1,gamma);uu2=MCG2D(5,u2-zeta*ss2_2,u2,gamma);
u1=imfilter(uu1,filter2,'conv');u2=imfilter(uu2,filter2,'conv');
%%%%%%%%%%u-subproblem%%%%%%%%%%%%%

for i=2:128
for j=2:128
phy_x(i,j)=i+u2(i,j);phy_y(i,j)=j+u1(i,j);deltx=(phy_x(i,j)-fix(phy_x(i,j)));delty=(phy_y(i,j)-fix(phy_y(i,j)));m1=fix(phy_x(i,j));m2=fix(phy_y(i,j));
if m1<=1
m1=1;
elseif m1>=128
m1=128;
end
if m2<=1
m2=1;
elseif m2>=128
m2=128;
end
TTk(i,j)=(1-deltx)*(1-delty)*double(TT(m1,m2))+deltx*(1-delty)*double(TT(m1+1,m2))+delty*(1-deltx)*double(TT(m1,m2+1))+deltx*delty*double(TT(m1+1,m2+1));
end
end

%%%%%%%%%%%%%%m-subproblem%%%%%%%%%%%%%%%%%
gm1=zeros(N1,N2,K);gm2=zeros(N1,N2);term1=zeros(N1,N2);
for l=1:K
for i=2:128
for j=2:128
term1(i,j) = log(abs(TTk(i,j) - s(i,j)) + 0.00001);
gm1(i,j,l) = term1(i,j) - cc(l,1);
gm2(i,j) = -log(abs(D(i,j)))+ term1(i,j) ;
end
end
end

gm=zeros(N1,N2);  qq_i=zeros(N1,N2);
for l=1:K
qq_i=qq_i+qi(:,:,l).*gm1(:,:,l);
end
for i=2:128
for j=2:128
gm(i,j)=(beta2*qq_i(i,j)+zeta*gm2(i,j))/(beta2+zeta);
end
end

ITE1=10;m=optimize_m(beta2,zeta,gm,eta1,ITE1);
Min=min(min(m));Max=max(max(m));
if Min<K1
for i=2:128
for j=2:128
if m(i,j)<0
m(i,j)=abs(m(i,j))/abs(Min)*K1;
end
end
end 
end
if Max>K2
for i=2:128
for j=2:128
if m(i,j)>0
m(i,j)=abs(m(i,j))/abs(Max)*K2;
end
end
end 
end
%%%%%%%%%%m-subproblem%%%%%%%%%%

%%%%%%%%%%s-subproblem%%%%%%%%%%
Tk=zeros(N1,N2,K);t=zeros(N1,N2);t1=zeros(N1,N2);g1=zeros(N1,N2);
for l=1:K
for i=2:128
    for j=2:128
        Tk(i,j,l)=qi(i,j,l)*(log(abs(TTk(i,j)-s(i,j))+0.000001)-m(i,j)-cc(l,1));
    end
end
end
for l=1:K
        t=t+Tk(:,:,l);
end

for i=2:128
   for j=2:128
    g1(i,j)=-beta2*t(i,j)+zeta*(m(i,j)+log(abs(D(i,j)))-log(abs(TTk(i,j)-s(i,j))+0.000001));
   end
end
s=s_subproblem(eta2,D,g1,kappa0,kappa);s=imfilter(s,filter1,'conv');
%%%%%%%%%%s-subproblem%%%%%%%%%%

%%%%%%%%%%p-subproblem%%%%%%%%%%
ggs=zeros(129,129,K);
for l=1:K
for i=2:128
for j=2:128
    ggs(i,j,l)=(log(T(i,j))-c(l,1))^2+8*log(2);
end
end
end

for l=1:K
    p_i(:,:,l)=optimize_p(lambda1,p_i(:,:,l),ppi(:,:,l),beta1,ggs(:,:,l),kexi,qi(:,:,l));
end

pmin = min(p_i, [], 3); 
for l=1:K
    p_i(:,:,l)=p_i(:,:,l)-pmin+0.01;
end
denominator = sum(p_i, 3);
for l=1:K
for i=2:128
for j=2:128
    p_i(i,j,l)=p_i(i,j,l)/denominator(i,j);
end
end
end
output1 = zeros(129, 129, 3);
[~, max_idx] = max(p_i, [], 3);  
output1(:,:,1) = (max_idx == 1);  
output1(:,:,2) = (max_idx == 2); 
output1(:,:,3) = (max_idx == 3);  

for l=1:K
for i=2:128
for j=2:128
phy_x(i,j)=i+u2(i,j);phy_y(i,j)=j+u1(i,j);deltx=(phy_x(i,j)-fix(phy_x(i,j)));delty=(phy_y(i,j)-fix(phy_y(i,j)));m1=fix(phy_x(i,j));m2=fix(phy_y(i,j));
if m1<=1
m1=1;
elseif m1>=128
m1=128;
end
if m2<=1
m2=1;
elseif m2>=128
m2=128;
end
ppi(i,j,l)=(1-deltx)*(1-delty)*double(p_i(m1,m2,l))+deltx*(1-delty)*double(p_i(m1+1,m2,l))+delty*(1-deltx)*double(p_i(m1,m2+1,l))+deltx*delty*double(p_i(m1+1,m2+1,l));
end
end
end

gg=zeros(129,129,K);
for l=1:K
for i=2:128
for j=2:128
    gg(i,j,l)=(log(abs(TTk(i,j)-s(i,j))+0.000001)-m(i,j)-cc(l,1))^2+8*log(2);
end
end
end

for l=1:K
qi(:,:,l)=optimize_q(ppi(:,:,l),lambda2,beta2,gg(:,:,l),kexi,qi(:,:,l));
end

qmin = min(qi, [], 3);  
for l=1:K
    qi(:,:,l)=qi(:,:,l)-qmin+0.01;
end
denominator1 = sum(qi, 3);
for l=1:K
    for i=2:128
        for j=2:128
    qi(i,j,l)=qi(i,j,l)/denominator1(i,j);
        end
    end
end

output = zeros(129, 129, 3);
[~, max_i1dx] = max(qi, [], 3);  
output(:,:,1) = (max_i1dx == 1);  
output(:,:,2) = (max_i1dx == 2);  
output(:,:,3) = (max_i1dx == 3); 
%%%%%%%%%%q-subproblem%%%%%%%%%%

%%%%%%%%%%%%%%Theta-subproblem%%%%%%%%%%%%

for l=1:K
    term21=sum(sum(p_i(:,:,l).*log(T(:,:))));
    term22=sum(sum(p_i(:,:,l)));
    c(l,1)=term21/term22;
end

for l=1:K
    term31=sum(sum(qi(:,:,l).*(log(abs(TTk(:,:)-s(:,:))+0.00001)-m(:,:))));
    term32=sum(sum(qi(:,:,l)));
    cc(l,1)=term31/term32;
end

%%%%%%%%%%%%%%Theta-subproblem%%%%%%%%%%%%
n=n+1;
 end

v=1;ITE1=100;
while v<ITE1
 beta=5;ggg=zeros(N1,N1,K);kexi1=0.2;
 for l=1:K
for i=2:128
for j=2:128
    ggg(i,j,l)=(log(abs(T(i,j)-s(i,j))+0.000001)-m(i,j)-ccc(l,1))^2+8*log(2);
end
end
end
 for l=1:K
 qqi(:,:,l)=optimize_p(lambda1,qqi(:,:,l),ppi(:,:,l),beta,ggg(:,:,l),kexi1,qi(:,:,l));
 end
qqmin = min(qqi, [], 3);
for l=1:K
    qqi(:,:,l)=qqi(:,:,l)-qqmin+0.01;
end
denominator12 = sum(qqi, 3);
for l=1:K
    for i=2:128
        for j=2:128
    qqi(i,j,l)=qqi(i,j,l)/denominator12(i,j);
        end
    end
end

output3 = zeros(129, 129, 3);
[~, max_i3dx] = max(qqi, [], 3);  
output3(:,:,1) = (max_i3dx == 1); 
output3(:,:,2) = (max_i3dx == 2); 
output3(:,:,3) = (max_i3dx == 3); 

for l=1:K
    term31=sum(sum(qqi(:,:,l).*(log(abs(TT(:,:)-s(:,:))+0.00001)-m(:,:))));
    term32=sum(sum(qqi(:,:,l)));
    ccc(l,1)=term31/term32;
end
v=v+1;
end
%%%%%%%%%%%%%%%%%%%%AMM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 elapsedTime = toc;
fprintf('runtime: %.4f sec\n', elapsedTime);
%%%%%%%%%%%%%%%%%%%%plot%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TTS=zeros(N1,N2);
  for i=2:128
    for j=2:128
 TTS(i,j)=(abs(TTk(i,j)-s(i,j)))/exp(m(i,j));
    end
  end
figure
imshow(TTS,[0,255]);
title('reg result')

 SSm=0;
SSm0=0.0;
for i=1:129
for j=1:129
SSm=SSm+(TTS(i,j)-D(i,j))^2;
SSm0=SSm0+(T(i,j)-D(i,j))^2;
end
end
Re_ssd=SSm/SSm0;
fprintf('Re_ssd： %.4f\n',Re_ssd );
I255 =draw_region_boundaries((max_i1dx == 1), (max_i1dx == 2),(max_i1dx == 3), TTS, ...
     'interface_red2.png', 1);
I256 =draw_region_boundaries((max_idx == 1), (max_idx == 2),(max_idx == 3), T, ...
     'interface_red.png', 1);
I254 =draw_region_boundaries((max_i3dx == 1), (max_i3dx == 2),(max_i3dx == 3), T, ...
     'interface_red3.png', 1);

figure;
subplot(1, 2, 1); 
imshow(T, [0,255]);  
title('T'); 

subplot(1, 2, 2);  
imshow(D, []);
title('R'); 

figure;
subplot(1, 2, 1); 
imshow(I255, []);  
title('Seg result of T_c(\phi)'); 


subplot(1, 2, 2);  
imshow(I254, [0,255]);
title('Seg result of T'); 

figure;
subplot(1, 2, 1);    
imshow(output, []);   
title('T_c(\phi)');

subplot(1, 2, 2);   
imshow(output3, [0,255]);
title('T');
