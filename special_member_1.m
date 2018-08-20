%%stiffness matrix method for special members
%input
clc;
clear;
n=2; %number of members
I=[1.5 1]*0.0031; %moment of inertia in m4
A=[1.25 1]*0.15;
theta=[0 0];
L=[4 2]; %length in m
m=[1 2]; %member number
uu=3; %Number of unrestrained degrees of freedom
ur=6; %Number of restrained degrees of freedom
uul=[1 2 3]; %Global labels of unrestrained dof
url=[4 5 6 7 8 9]; %Global labels for restrained dof
l1=[4 1 6 2 8 3]; %Global labels for member 1
l2=[1 5 2 7 3 9]; %Global labels for member 2

l=[l1; l2];
Ktotal=zeros(uu+ur);
fem1=[26.667; -26.667; 40; 40; 0; 0]; %Local fixed end moments of member 1
fem2=[0; 0; 0; 0; 0; 0]; %Local fixed end moments of member 2

rc1=4.*I./L;
rc2=2.*I./L;
rc3=A./L;
cx=cosd(theta);
cy=sind(theta);
%%Stiffness matrix 4 by 4 (axial deformation neglected)
for i=1:n
    Knew=zeros (uu+ur);
    k1=[rc1(i); rc2(i); (rc1(i)+rc2(i))/(L(i)); (-(rc1(i)+rc2(i))/L(i)); 0; 0];
    k2=[rc2(i);rc1(i);(rc1(i)+rc2(i))/(L(i));(-(rc1(i)+rc2(i))/L(i)); 0; 0];
    k3=[(rc1(i)+rc2(i))/L(i);(rc1(i)+rc2(i))/(L(i));(2*(rc1(i)+rc2(i))/(L(i)^2));(-2*(rc1(i)+rc2(i))/(L(i)^2)); 0; 0];
    k4=-k3;
    k5=[0; 0; 0; 0; rc3(i); -rc3(i)];
    k6=-k5;
    K=[k1 k2 k3 k4 k5 k6];
    fprintf('member number=');
    disp(i);
    fprintf('stiffness matrix=\n');
    disp(K);
    T1=[1; 0; 0; 0; 0; 0];
    T2=[0; 1; 0; 0; 0; 0];
    T3=[0; 0; cx(i); 0; cy(i); 0];
    T4=[0; 0; 0; cx(i); 0; cy(i)];
    T5=[0; 0; -cy(i); 0; cx(i); 0];
    T6=[0; 0; 0; -cy(i); 0; cx(i)];
    T=[T1 T2 T3 T4 T5 T6];
    fprintf('Transformation matrix of member, [T]=\n');
    disp(T);
    Ttr=T';
    Kg=Ttr*K*T;
    fprintf('Global Matrix, [K global]=\n');
    disp(Kg);
    for p=1:6
        for q=1:6
            Knew((l(i,p)),(l(i,q)))=Kg(p,q);
        end
    end
    Ktotal=Ktotal+Knew;
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
fprintf('Stiffness Matrix of compete structure, [Ktotal]=\n');
  disp(Ktotal);
  Kunr=zeros (uu);
  for x=1:uu
      for y=1:uu
          Kunr(x,y)=Ktotal(x,y);
      end
  end
  Kur=zeros(uu,ur);
  Kru=zeros(ur,uu);
  Krr=zeros(ur,ur);
  for x=1:uu
      for y=uu+1:uu+ur
          Kur(x,y-uu)=Ktotal(x,y);
      end
  end
  for x=uu+1:uu+ur
      for y=1:uu
          Kru(x-uu,y)=Ktotal(x,y);
      end
  end
  for x=uu+1:uu+ur
      for y=uu+1:uu+ur
          Krr(x-uu,y-uu)=Ktotal(x,y);
      end
  end
  fprintf('Unrestrained Stiffness sub-matrix, [Kuu]=\n');
  disp(Kunr);
  KuuInv=inv(Kunr);
  fprintf('Inverse of Unrestrained sub-matrix, [KuuInverse]=\n');
  disp(KuuInv);
  Kmod=Krr-Kru*KuuInv*Kur;
  fprintf('Modified stiffness matrix for Special member AC');
  disp(Kmod);
  jl=[26.667 -40 0 -26.667 0 -40 0 0 0]';
  jlu=jl(1:uu,1);
  jlr=jl(uu+1:uu+ur,1);
  delu=KuuInv*jlu;
  fprintf('unrestrained displacement=\n');
  disp(delu);
  % Restrained end reactions
  Rr=(Kru*KuuInv*jlu)-jlr;
  fprintf('Reactions at restrained dof=\n');
  disp(Rr);
  fprintf('Joint Load vector, [Jl]=\n');
  disp(jl);
  delr=zeros(ur,1);
  del=zeros(uu+ur,1);
  del=[delu; delr];
  deli=zeros (6,1);
  for i=1:n

      for p=1:6
          deli(p,1)=del((l(i,p)),1);
      end
      if i==1
          delbar1=deli;
          mbar1=(Kg1*delbar1)+fembar1;
          fprintf('Member Number=');
          disp(i);
          fprintf('Global displacement matrix [DeltaBar]=\n');
          disp(delbar1);
          fprintf('Global End moment matrix [MBar]=\n');
          disp(mbar1);
      elseif i==2
          delbar2=deli;
          mbar2=(Kg2*delbar2)+fembar2;
          fprintf('Member Number =');
          disp(i);
          fprintf('Global displacement matrix [DeltaBar]=\n');
          disp(delbar2);
          fprintf('Global End moment matrix [MBar]=\n');
          disp(mbar2);

      end
  end
  
           
  