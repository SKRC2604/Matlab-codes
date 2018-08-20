%% TRUSS PROBLEMS
% Input
clc;
clear;
n=4;   % Number of members
I=[1,1,1,1,1];   %Moment of inertia in m4
L=[4,4,4,5.656,5.656];   %length in m
A=[0.5,0.4,0.5,0.6,0.6];   %Member in number 
theta=[90,0,90,-45,45];   %Angle in degrees
uu=4;   %Number of unrestrained degree of freedom
ur=4;   %Number of restrained degree of freedom
uul=[1,2,3,4];   %Global labels of unrestrained degree of freedom
url=[4,5,6,7,8];   %Global labels of restrained degree of freedom
l1=[6,2,5,1];   %Global labels of member 1
l2=[2,4,1,3];   %Global labels of member 2
l3=[8,4,7,3];   %Global labels of member 3
l4=[2,8,1,7];   %Global labels of member 4

l=[l1;l2;l3;l4];
dof=uu+ur;   %Degree of freedom
Ktotal=zeros(dof);
Tt1=zeros(4);   %Transformation matrix for member 1
Tt2=zeros(4);   %Transformation matrix for member 2
Tt3=zeros(4);   %Transformation matrix for member 3
Tt4=zeros(4);   %Transformation matrix for member 4

fem1=[0;0;0;0];   %Local fixed end moments of member 1
fem2=[0;0;0;0];  %Local fixed end moments of member 2
fem3=[0;0;0;0];   %Local fixed end moments of member 3
fem4=[0;0;0;0];  %Local fixed end moments of member 4

rc1 = A./L;
cx = cosd(theta);
cy = sind(theta);
for i =1:n
        Knew = zeros(dof);
        k1 = [0;0;0;0];
        k2 = k1;
        k3 = [0;0;rc1(i);-rc1(i)];
        k4 = -k3;
        K = [k1 k2 k3 k4];
        fprintf('Member Number =');
        disp(i);
        fprintf('Local stiffness matrix of member ,[K] = \n');
        disp (K);
        T1 = [cx(i);0;cy(i);0];
        T2 = [0;cx(i);0;cy(i)];
        T3 = [-cy(i);0;cx(i);0];
        T4 = [0;-cy(i);0;cx(i)];
        T = [T1 T2 T3 T4];
        fprintf('Transformation matrix of member,[T] = \n');
        disp(T);
        Ttr = T';
        fprintf('Transformation matrix Transpose,[T] = \n');
        disp(Ttr);
        Kg = Ttr*K*T;
        fprintf('Global Matrix, [K global] = \n');
        disp(Kg);
        for p = 1:4
            for q = 1:4
                Knew((l(i,p)),(l(i,q))) = Kg(p,q);
            end
        end
        Ktotal = Ktotal + Knew;
        if i == 1
            Tt1 = T;
            Kg1 = Kg;
            fembar1 = Tt1'*fem1;
        elseif i ==2
            Tt2 = T;
            Kg2= Kg;
            fembar2 = Tt2'*fem2;
        elseif i==3
            Tt3 = T;
            Kg3 = Kg;
            fembar3 = Tt3'*fem3;
        elseif i==4
            Tt4 = T;
            Kg4 = Kg;
            fembar4 = Tt4'*fem4;
     
        end
end
fprintf('Stiffness Matrix of complete structure, [Ktotal] = \n');
disp(Ktotal);
Kunr = zeros(uu);
for x = 1:uu
    for y = 1:uu
        Kunr(x,y) = Ktotal(x,y);
    end
end
fprintf('Unrestrained Stiffness sub-matrix,[Kuu] = \n');
disp(Kunr);
KuuInv = inv(Kunr);
fprintf('Inverse of Unrestrained Stiffness sub-matrix ,[KunInverse] = \n');
disp(KuuInv);
jl= [40;0;0;-20;0;0;0;0];
jlu = [40;0;0;-20];
delu = KuuInv*jlu;
fprintf('Joint Load vector, [J1] = \n');
disp(jl);
fprintf('Unrestrained displacements, [DelU] = \n');
disp(delu);
delr = [0;0;0;0];
del = zeros(dof ,1);
del = [delu;delr];
deli = zeros(4 ,1);
for i = 1:n
    for p = 1:4
        deli(p,1) = del((l(i,p)),1);
    end
    if i==1
          delbar1 = deli;
          mbar1 = (Kg1 * delbar1) + fembar1;
          fprintf('Member Number = ');
          disp(i);
          fprintf('Global displacement of matrix [Deltabar] = \n');
          disp(delbar1);
          fprintf('Global End moment matrix [MBar] = \n');
          disp(mbar1);
    elseif i==2 
              delbar2 = deli;
              mbar2 = (Kg2*delbar2) + fembar2;
              fprintf('Member Number = ');
              disp(i);
              fprintf('Global displacement of matrix [Deltabar] = \n');
              disp(delbar2);
              fprintf('Global End moment matrix [MBar] = \n');
              disp(mbar2);
              
       elseif i==3 
              delbar3= deli;
              mbar3 = (Kg3*delbar3) + fembar3;
              fprintf('Member Number = ');
              disp(i);
              fprintf('Global displacement of matrix [Deltabar] = \n');
              disp(delbar3);
              fprintf('Global End moment matrix [MBar] = \n');
              disp(mbar3);
      
       elseif i==4
              delbar4= deli;
              mbar4 = (Kg4*delbar4) + fembar4;
              fprintf('Member Number = ');
              disp(i);
              fprintf('Global displacement of matrix [Deltabar] = \n');
              disp(delbar4);
              fprintf('Global End moment matrix [MBar] = \n');
              disp(mbar4);
              
   
              
              
            
              
              
    end
end