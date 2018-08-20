%% PLANER NON-ORTHOGONAL PROBLEMS
% Input
clc;
clear;
n=2;   % Number of members
I=[2,3];   %Moment of inertia in m4
L=[4, 5];   %length in m
A=[1,1];   %Member in number 
theta=[0,126.87];   %Angle in degrees
uu=3;   %Number of unrestrained degree of freedom
ur=6;   %Number of restrained degree of freedom
uul=[1,2,3];   %Global labels of unrestrained degree of freedom
url=[4,5,6,7,8,9];   %Global labels of restrained degree of freedom
l1=[4,1,7,3,6,2];   %Global labels of member 1
l2=[5,1,9,3,8,2];   %Global labels of member 2
   %Global labels of member 2
l=[l1;l2];
dof=uu+ur;   %Degree of freedom
Ktotal=zeros(dof);
Tt1=zeros(6);   %Transformation matrix for member 1
Tt2=zeros(6);   %Transformation matrix for member 2
  %Transformation matrix for member 3

fem1=[0;0;0;0;0;0];   %Local fixed end moments of member 1
fem2=[0;0;0;0;0;0];  %Local fixed end moments of member 2
  %Local fixed end moments of member 3


%% Rotation co-efficient for each member
rc1=4.*I./L;
rc2=2.*I./L;
rc3=A./L;
cx=cosd(theta);
cy=sind(theta);

%% Stiffness matrix 4 by 4 (Axial deformation neglected)
for i=1:n
    Knew=zeros(dof);
    k1=[rc1(i);rc2(i);(rc1(i)+rc2(i))/L(i); (-(rc1(i)+rc2(i))/L(i));0;0];
    k2=[rc2(i);rc1(i);(rc1(i)+rc2(i))/L(i);(-(rc1(i)+rc2(i))/L(i));0;0];
    k3=[(rc1(i)+rc2(i))/L(i);(rc1(i)+rc2(i))/L(i); (2*(rc1(i)+rc2(i))/(L(i)^2));(-2*(rc1(i)+rc2(i))/(L(i)^2));0;0];
    k4=-k3;
    k5=[0;0;0;0;rc3(i);-rc3(i)];
    k6=[0;0;0;0;-rc3(i);rc3(i)];
    K=[k1,k2,k3,k4,k5,k6];
    fprintf('Member Number=');
    disp(i);
    fprintf('Local Stiffness matrix of member, [K]=');
    disp(K);
   T1=[1;0;0;0;0;0];
   T2=[0;1;0;0;0;0];
   T3=[0;0;cx(i);0;cy(i);0];
   T4=[0;0;0;cx(i);0;cy(i)];
   T5=[0;0;-cy(i);0;cx(i);0];
   T6=[0;0;0;-cy(i);0;cx(i)];
   T=[T1,T2,T3,T4,T5,T6];
   fprintf('Transformation matrix of member, [T]=\n');
    disp(T);
    Ttr=T';
    fprintf('Transformation matrix Transpose, [T]=\n');
    disp(Ttr);
    Kg=Ttr*K*T;
    fprintf('Global Matrix, [K global]=\n');
    disp(Kg);
   
        for p=1:6
            for q=1:6
            Knew((l(i,p)),(l(i,q)))=Kg(p,q);
         end
    end
    Ktotal=Ktotal + Knew;
    if i==1
        Tt1=T;
        Kg1=Kg;
        fembar1=Tt1'*fem1;
    elseif i==2
        Tt2=T;
        Kg2=Kg;
        fembar2=Tt2'*fem2;
   
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
fprintf('unrestrained stiffness sub-matrix,[Kuu]=\n');
disp(Kunr);
KuuInv=inv(Kunr);
fprintf('Inverse of unrestrained stiffness sub-matrix,[KuuInverse]=\n');
disp(KuuInv);

%% Creation of joint load vector
jl=[0 0 -50 0 0 0 0 0 0]';   %Values given in kN or kNm
jlu=[0;0;-50];   %Load vector in unrestrained dof
delu=KuuInv*jlu;
fprintf('Joint Load Vector,[jl]=\n');
disp(jl);
fprintf('Unrestrained displacement,[DelU]=\n');
disp(delu);
delr=[0;0;0;0;0;0;0;0;0;0];
del=zeros(dof,1);
del=[delu;delr];
deli=zeros(6,1);
for i=1:n
    for p=1:6
        deli(p,1)=del((l(i,p)),1);
    end
    if i==1
        delbar1=deli;
        mbar1=(Kg1 * delbar1)+fem1;
        fprintf('Member Number=');
        disp(i);
        fprintf('Global displacement matrix [DeltaBar]=\n');
        disp(delbar1)
         fprintf('Global End moment matrix [MBar]=\n');
         disp(mbar1);
         
    elseif i==2
        delbar2=deli;
        mbar2=(Kg2*delbar2)+fem2;
         fprintf('Member Number');
         disp(i);
          fprintf('Global displacement matrix [DeltaBar]=\n');
          disp(delbar2);
           fprintf('Global global end moment matrix [MBar]=\n');
           disp(mbar2);
           
         
    end
end

    
    
    