
clear;
n = 3; % number of members
uu=12;
ur=12;
dof=uu+ur;
EI = 1; %Flexural rigidity
EIy = EI;
EIz = EI;
GI = (0.25).*EI; %Torsional constant
EA = (0.25).*EI; %Axial rigidity
L = 3; % length in m  
%codm = [0 5 0;3 5 0; 6 5 0; 0 0 0; 3 0 0; 6 0 0]; %Coordinate wrt X,Y.Z: 
%size=nj,3
%dc = [1 0 0; 1 0 0; 1 0 0; 1 0 0; 0 1 0]; % Direction cosines for each member
%tytr = [1 1 1 1 2]; % Type of transformation fo each member
%psi = [0 0 0 0 90]; % Psi angle in degrees for each member

T=[0 0 -1 0 0 0 0 0 0 0 0 0;
    0 1 0 0 0 0 0 0 0 0 0 0;
    1 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 -1 0 0 0 0 0 0;
    0 0 0 0 1 0 0 0 0 0 0 0;
    0 0 0 1 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 -1 0 0 0;
    0 0 0 0 0 0 0 1 0 0 0 0;
    0 0 0 0 0 0 1 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 -1;
    0 0 0 0 0 0 0 0 0 0 1 0;
    0 0 0 0 0 0 0 0 0 1 0 0];
fprintf('T of member= \n');
disp(T);
fprintf('Ttr of member= \n');
disp(T');
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

    Knew = zeros (dof);
    k1 = [sc1; 0; 0; 0; 0; 0; -sc3; 0; 0; 0; 0; 0];
    k2 = [0; sc6; 0; 0; 0; sc2; 0; -sc6; 0; 0; 0; sc2];
    k3 = [0; 0; sc7; 0; -sc3; 0; 0; 0; -sc7; 0; -sc3; 0];
    k4 = [0; 0; 0; sc4; 0; 0; 0; 0; 0; -sc4; 0; 0];
    k5 = [0; 0; -sc3; 0; (2*sc5); 0; 0; 0; sc3; 0; sc5; 0];
    k6 = [0; sc2; 0; 0; 0; (2*sc8); 0; -sc2; 0; 0; 0; sc8];
    k7 = -k1;
    k8 = -k2;
    k9 = -k3;
    k10 = -k4;
    k11 = [0; 0; -sc3; 0; sc5; 0; 0; 0; sc3; 0;(2*sc5);0];
    k12 = [0; sc2; 0; 0; 0; sc5; 0; -sc2; 0; 0; 0;(2*sc8)];
    K = [k1 k2 k3 k4 k5 k6 k7 k8 k9 k10 k11 k12]; 
   % fprintf ('Member Number =');
    fprintf ('Local Stiffness matrix of member, [K] = \n'); 
    disp (K);
    
    Ttr=T';
    Kg=Ttr*K*T;
    fprintf('Global Stiffness Matrix, [Kg]= \n');
    disp(Kg);

