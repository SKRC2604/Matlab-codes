%% Stiffness matrix method
% Input
clc;
clear;
n=3;   % Number of members
I=[1.5,2,1];   %Moment of inertia in m4
L=[3,4,4];   %length in m
m=[1,2,3];   %Member in number 
uu=2;   %Number of unrestrained degree of freedom
ur=6;   %Number of restrained degree of freedom
dof=uu+ur;   %Degree of freedom
uul=[1,2];   %Global labels of unrestrained degree of freedom
url=[3,4,5,6,7,8];   %Global labels of restrained degree of freedom
l1=[3,1,4,5];   %Global labels of member 1
l2=[1,2,5,6];   %Global labels of member 2
l3=[2,7,6,8];   %Global labels of member 3
l=[l1;l2;l3];
Ktotal=zeros(dof);
fem1=[35.56,-17.78,53.33,26.67];   %Local fixed end moments of member 1
fem2=[66.67,-66.67,100,100];  %Local fixed end moments of member 2
fem3=[20,-20,20,20];    %Local fixed en3 moments of member 3

%% Rotation co-efficient for each member
rc1=4.*I./L;
rc2=2.*I./L;

%% Stiffness matrix 4 by 4 (Axial deformation neglected)
for i=1:n
    Knew=zeros(dof);
    k1=[rc1(i);rc2(i);(rc1(i)+rc2(i))/L(i); (-(rc1(i)+rc2(i))/L(i))];
    k2=[rc2(i);rc1(i);(rc1(i)+rc2(i))/L(i);(-(rc1(i)+rc2(i))/L(i))];
    k3=[(rc1(i)+rc2(i))/L(i);(rc1(i)+rc2(i))/L(i); (2*(rc1(i)+rc2(i))/(L(i)^2));(-2*(rc1(i)+rc2(i))/(L(i)^2))];
    k4=-k3;
    K=[k1,k2,k3,k4];
    fprintf('Member Number=');
    disp(i);
    fprintf('Local Stiffness matrix of member, [K]=');
    disp(K);
        for p=1:4
            for q=1:4
            Knew((l(i,p)),(l(i,q)))=K(p,q);
         end
    end
    Ktotal=Ktotal + Knew;
    if i==1
        Kg1=K;
    elseif i==2
        Kg2=K;
    elseif i==3
        Kg3=K;
    end
end
fprintf('Stiffness matrix of complete structure,[Ktotal]=\n');
disp(Ktotal);
Kunr=zeros(uu);
for x=1:uu
    for y=1:uu
        Kunr(x,y)=Ktotal(x,y);
    end
end
fprintf('unrestrained stiffness sub-matrix,[kuu]=\n');
disp(Kunr);
KuuInv=inv(Kunr);
fprintf('Inverse of unrestrained stiffness sub-matrix,[kuuInverse]=\n');
disp(KuuInv);

%% Creation of joint load vector
jl=[-48.89;46.67;-35.56;-53.33;-126.67;-120;20;-20];
jlu=[-48.89;46.67];
delu=KuuInv*jlu;
fprintf('Joint Load Vector,[jl]=\n');
disp(jl);
delr=[0;0;0;0;0;0;0;0];
del=zeros(dof,1);
del=[delu;delr];
deli=zeros(4,1);
for i=1:n
    for p=1:4
        deli(p,1)=del((l(i,p)),1);
    end
    if i==1  
        delbar1=deli;
        mbar1=(Kg1 * delbar1)+fem1';
        fprintf('Member Number=');
        disp(i);
        fprintf('Global displacement matrix [DeltaBar]=\n');
        disp(delbar1)
         fprintf('Global End moment matrix [MBar]=\n');
         disp(mbar1);
    elseif i==2
        delbar2=deli;
        mbar2=(Kg2*delbar2)+fem2';
         fprintf('Member Number');
         disp(i);
          fprintf('Global displacement matrix [DeltaBar]=\n');
          disp(delbar2);
           fprintf('Global global end moment matrix [MBar]=\n');
           disp(mbar2);
           
            elseif i==3
        delbar3=deli;
        mbar3=(Kg3*delbar3)+fem3';
         fprintf('Member Number');
         disp(i);
          fprintf('Global displacement matrix [DeltaBar]=\n');
          disp(delbar3);
           fprintf('Global global end moment matrix [MBar]=\n');
           disp(mbar3);
    end
end

    
    
    