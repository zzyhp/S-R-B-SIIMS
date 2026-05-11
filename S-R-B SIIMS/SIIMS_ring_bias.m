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
data=load('O-Oms.mat');

%%%%%%%%%%%%%%%Parameter%%%%%%%%%%%%
% All other parameters remain unchanged; only parameter kexi=0;0.04;0.06;0.08;0.1 is adjusted.
K=2;beta2=0.01;K1=-10.5;K2=16;zeta=1;kexi=0.0;gamma=5;ITE=100;mu2=0.00000000000001;beta1=0.03;mu1=0.0000000000000;
T=data.T;D=data.D;kappa0=-0.001;kappa=0.001;lambda1=0.15;lambda2=0.01;eta1=0.01;eta2=0.00000000001;

%%%%%%%%%%%%%%%%%parameter%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%processing%%%%%%%%%%%
tmin=min(min(T));dmin=min(min(D));
for i=1:129
for j=1:129
T(i,j)=T(i,j)-tmin+0.01;D(i,j)=D(i,j)-dmin+0.01;
end
end 
%%%%%%%%%%%%%%%%%%%processing%%%%%%%%%%%

%%%%%%%%%%%%%%Initial value%%%%%%%%
[idx, C1] = kmeans1(T, K);
[idx1, C2] = kmeans1(T, K);
c=zeros(K,1);cc=zeros(K,1);
for l=1:K
    c(l,1)=log(C1(l,1));
    cc(l,1)=log(C2(l,1));
end
u1=zeros(N1,N2);
u2=zeros(N1,N2);
phy_x=zeros(N1,N2);
phy_y=zeros(N1,N2);
filter1=fspecial('gaussian',3,1000000.1);
filter2=fspecial('gaussian',3,1.1);
qi=(1/K)*ones(N1,N2,K);
p_i=(1/K)*ones(N1,N2,K);
ppi=p_i;
s=zeros(N1,N2);
TT=T;TTk=T;
m=zeros(N1,N2);

for i=1:129
for j=1:129
phy_x(i,j)=i;phy_y(i,j)=j;
end
end
%%%%%%%%%%%%%%Initial value%%%%%%%%


tic
%%%%%%%%%%%%%%%%%%%%AMM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=1;
 while n<ITE
%%%%%%%%%%%%%%%u—subproblem%%%%%%%%%
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


uu1=MCG2D(5,u1-zeta*ss2_1,u1,gamma);
uu2=MCG2D(5,u2-zeta*ss2_2,u2,gamma);
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

ITE1=10;
m=optimize_m(beta2,zeta,gm,eta1,ITE1);
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
% m=imfilter(m,filter1,'conv');
%%%%%%%%%%%%m-subproblem%%%%%%%%
E=m-gm;


%%%%%%%%%%%%%%s—subproblem%%%%%%
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

%%%%%%%%%%%%%%s—subproblem%%%%%%

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
mask_l11 = p_i(:,:,1) > p_i(:,:,2); 
mask_l22 = p_i(:,:,2) > p_i(:,:,1); 
output1(:,:,1) = mask_l11; 
output1(:,:,2) = mask_l22; 

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
%%%%%%%%%%q-subproblem%%%%%%%%%%
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
mask_l1 = qi(:,:,1) > qi(:,:,2);
mask_l2 = qi(:,:,2) > qi(:,:,1); 
output(:,:,1) = mask_l1; 
output(:,:,2) = mask_l2; 

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
elapsedTime = toc;
fprintf('runtime: %.4f sec\n', elapsedTime);
%%%%%%%%%%%%%%%%AMM%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%plot%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 TTS=zeros(N1,N2);
  for i=2:128
    for j=2:128
 TTS(i,j)=(abs(TTk(i,j)-s(i,j)))/exp(m(i,j));
    end
  end
  figure
 imshow(TTS,[0,255]);
 title('bias correction result')
  figure
I255 = single_interface_curve_rgb(mask_l1, mask_l2, TTS, ...
    'interface_red2.png', 1.2); 
imshow(I255);
title('Segmentation result of T_c(\phi)')

 figure
I256 = single_interface_curve_rgb(mask_l11, mask_l22, T, ...
    'interface_red.png', 1.2); 
imshow(I256);
title('Segmentation result of T')

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




