%% 3D analysis of space frame
clc;
clear;

%% Input
n = 5; % number of members
EI = [1 1 1 1 1]; %Flexural rigidity
EIy = EI;
EIz = EI;
GI = [0.25 0.25 0.25 0.25 0.25].*EI; %Torsional constant
EA = [0.25 0.25 0.25 0.25 0.25].*EI; %Axial rigidity
L = [3 3 3 3 5]; % length in m 
nj = n+1; % Number of Joints
codm = [0 5 0;3 5 0; 6 5 0; 0 0 0; 3 0 0; 6 0 0]; %Coordinate wrt X,Y.Z: 
size=nj,3
dc = [1 0 0; 1 0 0; 1 0 0; 1 0 0; 0 1 0]; % Direction cosines for each member
tytr = [1 1 1 1 2]; % Type of transformation fo each member
psi = [0 0 0 0 90]; % Psi angle in degrees for each member

% C matrix
c1 =[1 0 0; 0 1 0; 0 0 1]; % C	matrix for	member 1
c2 =[1 0 0; 0 1 0; 0 0 1]; % C	matrix for	member 2
c3 =[1 0 0; 0 1 0; 0 0 1]; % C	matrix for	member 3
c4 =[1 0 0; 0 1 0; 0 0 1]; % C	matrix for	member 4
c5 =[0 1 0; 0 0 1; 1 0 0]; % C	matrix for	member 5

uu = 12; % Number of unrestrained degrees of freedom 
ur = 24; % Number of restrained degrees of freedom
uul = [1 2 3 4 5 6 7 8 9 10 11 12]; % global labels of unrestained dof
url = [13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36];	% global labels	of	restained dof

l1 =[13 14 15 16 17	18 1 2 3 4 5 6]; % Global labels for member 1
l2 =[1 2 3 4 5 6 19	20 21 22 23 24]; % Global labels for member 2
l3 =[25 26 27 28 29	30 7 8 9 10 11 12]; % Global labels for member 3
l4 =[7 8 9 10 11 12	31 32 33 34 35 36]; % Global labels for member 4
l5 =[7 8 9 10 11 12 1 2 3 4 5 6]; % Global labels for member 5
l= [l1; l2; l3; l4; l5];
dof = uu + ur; % Degrees of freedom 
Ktotal = zeros (dof);

fem1= [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]; % Local Fixed end moments of member 1
fem2= [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]; % Local Fixed end moments of member 2
fem3= [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]; % Local Fixed end moments of member 3
fem4= [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]; % Local Fixed end moments ofmember 4
fem5= [0; 0; 75; 62.5; 0; 0; 0; 0; 75; -62.5; 0; 0]; % Local Fixed end moments of member 5

%% Transformation matrix 
T1 = zeros(12);
T2 = zeros(12);
T3 = zeros(12); 
T4 = zeros(12);
T5 = zeros(12);

for i = 1:3
    for j = 1:3
        T1(i,j)=c1(i,j);
        T1(i+3,j+3)=c1(i,j);
        T1(i+6,j+6)=c1(i,j);
        T1(i+9,j+9)=c1(i,j);
        T2(i,j)=c2(i,j);
        T2(i+3,j+3)=c2(i,j);
        T2(i+6,j+6)=c2(i,j);
        T2(i+9,j+9)=c2(i,j);
        T3(i,j)=c3(i,j);
        T3(i+3,j+3)=c3(i,j);
        T3(i+6,j+6)=c3(i,j);
        T3(i+9,j+9)=c3(i,j);
        T4(i,j)=c4(i,j);
        T4(i+3,j+3)=c4(i,j);
        T4(i+6,j+6)=c4(i,j);
        T4(i+9,j+9)=c4(i,j);
        T5(i,j)=c5(i,j);
        T5(i+3,j+3)=c5(i,j);
        T5(i+6,j+6)=c5(i,j);
        T5(i+9,j+9)=c5(i,j);
    end
end 

%% Getting Type of transformation and Psi angle
for i = 1:n
    if tytr(i) ==1
        fprintf ('Member Number =');
        disp (i);
        fprintf ('Type of transformation is Y-Z-X \n');
    else
        fprintf ('Member Number ='); 
        disp (i);
        fprintf ('Type of transformation is Z-Y-X \n');
    end
 
end

fprintf ('Psi angle='); 
disp (psi(i));
 
%% Stiffness coefficients for each member 
sc1 = EA./L;
sc2 = 6*EIz./(L.^2); 
sc3 = 6*EIy./(L.^2);
sc4 = GI./L;
sc5 = 2*EIy./L;
sc6 = 12*EIz./(L.^3); 
sc7 = 12*EIy./(L.^3);
sc8 = 2*EIz./L;

%% stiffness matrix 6 by 6
for i = 1:n
    Knew = zeros (dof);
    k1 = [sc1(i); 0; 0; 0; 0; 0; -sc3(i); 0; 0; 0; 0; 0];
    k2 = [0; sc6(i); 0; 0; 0; sc2(i); 0; -sc6(i); 0; 0; 0; sc2(i)];
    k3 = [0; 0; sc7(i); 0; -sc3(i); 0; 0; 0; -sc7(i); 0; -sc3(i); 0];
    k4 = [0; 0; 0; sc4(i); 0; 0; 0; 0; 0; -sc4(i); 0; 0];
    k5 = [0; 0; -sc3(i); 0; (2*sc5(i)); 0; 0; 0; sc3(i); 0; sc5(i); 0];
    k6 = [0; sc2(i); 0; 0; 0; (2*sc8(i)); 0; -sc2(i); 0; 0; 0; sc8(i)];
    k7 = -k1;
    k8 = -k2;
    k9 = -k3;
    k10 = -k4;
    k11 = [0; 0; -sc3(i); 0; sc5(i); 0; 0; 0; sc3(i); 0;(2*sc5(i));0];
    k12 = [0; sc2(i); 0; 0; 0; sc5(i); 0; -sc2(i); 0; 0; 0;(2*sc8(i))];
    K = [k1 k2 k3 k4 k5 k6 k7 k8 k9 k10 k11 k12]; 
    fprintf ('Member Number =');
    disp (i);
    fprintf ('Local Stiffness matrix of member, [K] = \n'); 
    disp (K);
    if i == 1
        T = T1;
    elseif i == 2
        T = T2;
    else
        T = T3;
    end
    Ttr = T';
    Kg = Ttr*K*T;
    fprintf ('Global Matrix, [K global] = \n');
    disp (Kg);

    for p = 1:12
        for q = 1:12
            Knew((l(i,p)),(l(i,q))) =Kg(p,q);
        end
    end

    Ktotal = Ktotal + Knew;
    if i == 1
        Tt1= T;
        Kg1=Kg;
        fembar1= Tt1'*fem1;
    elseif i == 2
        Tt2 = T;
        Kg2 = Kg;
        fembar2= Tt2'*fem2;
    elseif i ==3
        Tt3 = T;
        Kg3 = Kg;
        fembar3= Tt3'*fem3;
    elseif i ==4
        Tt4 = T;
        Kg4 = Kg;
        fembar4= Tt4'*fem4;
    else
        Tt5 = T;
        Kg5 = Kg;
        fembar5= Tt5'*fem5;
    end
end

fprintf ('Stiffness Matrix of complete structure, [Ktotal] = \n');
disp (Ktotal);
Kunr = zeros(12);

for x=1:uu
    for y=1:uu
        Kunr(x,y)= Ktotal(x,y);
    end
end

fprintf ('Unrestrained Stiffness sub-matix, [Kuu] = \n');
disp (Kunr);
KuuInv= inv(Kunr);
fprintf ('Inverse of Unrestrained Stiffness sub-matix, [KuuInverse] = \n');
disp (KuuInv);

%% Creation of joint load vector
jl= [0; 0; -75; 62.5; 0; 0; 0; 0; -75; -62.5; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0]; % values given in kN or kNm
jlu = jl(1:12,1); % load vector in unrestrained dof 
delu = KuuInv*jlu;
fprintf ('Joint Load vector, [Jl] = \n');
disp (jl');
fprintf ('Unrestrained displacements, [DelU] = \n');
disp (delu');
delr = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
del = zeros (dof,1); 
del = [delu; delr]; 
deli= zeros (12,1);

for i = 1:n
    for p = 1:12
        deli(p,1) = del((l(i,p)),1) ;
    end
    
    if i == 1
        delbar1 = deli;
        mbar1= (Kg1 * delbar1)+fembar1;
        fprintf ('Member Number =');
        disp (i);
        fprintf ('Global displacement matrix [DeltaBar] = \n'); 
        disp (delbar1');
        fprintf ('Global End moment matrix [MBar] = \n');
        disp (mbar1');
    elseif i == 2
        delbar2 = deli;
        mbar2= (Kg2 * delbar2)+fembar2;
        fprintf ('Member Number =');
        disp (i);
        fprintf ('Global displacement matrix [DeltaBar] = \n');
        disp (delbar2');
        fprintf ('Global End moment matrix [MBar] = \n');
        disp (mbar2');
    elseif i ==3
        delbar3 = deli;
        mbar3= (Kg3 * delbar3)+fembar3; 
        fprintf ('Member Number ='); 
        disp (i);
        fprintf ('Global displacement matrix [DeltaBar] = \n');
        disp (delbar3');
        fprintf ('Global End moment matrix [MBar] = \n');
        disp (mbar3');
    elseif i ==4
        delbar4 = deli;
        mbar4= (Kg4 * delbar4)+fembar4; 
        fprintf ('Member Number ='); 
        disp (i);
        fprintf ('Global displacement matrix [DeltaBar] = \n');
        disp (delbar4');
        fprintf ('Global End moment matrix [MBar] = \n'); 
        disp (mbar4');
    else
        delbar5 = deli;
        mbar5= (Kg5 * delbar5)+fembar5;
        fprintf ('Member Number ='); 
        disp (i);
        fprintf ('Global displacement matrix [DeltaBar] = \n');
        disp (delbar5');
        fprintf ('Global End moment matrix [MBar] = \n');
        disp (mbar5');
    end
end
