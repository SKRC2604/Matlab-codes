 %% TRUSS PROBLEMS
% Input
clear;
clc;
n=3;
I=[1 1 1];
L=[3.16,4.242,4];
A=[.1575,.2,.2];
theta=[71.56,-45,0];
uu=3; % unrestrained degrees of freedom 
ur=3; % nrestrained degrees of freedom
uul=[1 2 3]; %Global labels of unrestraind dof
url=[ 4 5 6]; %Global labels of restraind dof
l1=[6 1 5 2]; % Global labels of degrees of freedom for member 1 
l2=[1 4 2 3];  % Global labels of degrees of freedom for member 2
l3=[6 4 5 3];  % Global labels of degrees of freedom for member 3
l=[l1;l2;l3];
dof=uu+ur;
Ktotal=zeros(dof); %For global K matrix for all members
Tt1=zeros(4); %Transformation matrix
Tt2=zeros(4); %Transformation matrix
Tt3=zeros(4); %Transformation matrix

fem1=[0 0 0 0 ]'; %Local fixed end moments
fem2=[0 0 0 0 ]';
fem3=[0 0 0 0 ]';
%% Elements of local K matrix

r = A./L;
cx=cosd(theta);
cy=sind(theta);
%% Finding Global K matrix for each member and final global K matrix 
for i=1:n;
    Knew=zeros(dof);
    k1=[0;0;0;0];
    k2=[0;0;0;0];
    k3=[0;0;r(i) ; -r(i)];
    k4=-k3;

    K=[k1 k2 k3 k4];
    fprintf('member number=');
    disp(i);
    fprintf('stiffness matrix=\n');
    disp(K);
   
    T1=[ cx(i); 0; cy(i); 0];
    T2=[ 0; cx(i); 0; cy(i)];
    T3=[ -cy(i); 0; cx(i); 0];
    T4=[ 0; -cy(i); 0; cx(i)];
    T=[T1 T2 T3 T4];
    fprintf('Transformation matrix of member, [T]=\n');
    disp(T);
    Ttr=T';
    Kg=Ttr*K*T;
    fprintf('Global Matrix, [K global]=\n');
    disp(Kg);
    for p=1:4 
        for q=1:4
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
    elseif i==3
        Tt3=T;
        Kg3=Kg;
        fembar3=Tt3'*fem3;
    end
end
fprintf('Stiffness Matrix of complete structure, [Ktotal]=\n');
disp(Ktotal);
Kunr=zeros(uu);
for x=1:uu
    for y=1:uu
        Kunr(x,y)=Ktotal(x,y);
    end
end
fprintf('unrestrained Stiffness Matrix, [Kuu]=\n');
disp(Kunr);
KuuInv=inv(Kunr);
fprintf('Inverse of unrestrained stiffness matrix, [KuuInv]=\n');
disp(KuuInv);
%% creation of joint load vector
jl=[-30; 20; 0; 0; 0; 0];
jlu=[-30; 20; 0];
delu=KuuInv*jlu;
fprintf('joint load vector, [Jl]=\n');
disp(jl);
fprintf('unrestrained displacement, [DelU]=\n');
disp(delu);
delr=zeros(ur,1);
del=zeros(dof,1);
del=[delu; delr];
deli=zeros(4,1);
for i=1:n
    for p=1:4
        deli(p,1)=del((l(i,p)), 1);
    end
    if i==1
        delbar1=deli;
        mbar=(Kg1*delbar1)+fembar1;
        fprintf('member number=');
        disp(i);
        fprintf('Global displacement matrix=\n');
        disp(delbar1);
        fprintf('Global end moment matrix=\n');
        disp(mbar);
    elseif i==2
        delbar2=deli;
        mbar=(Kg2*delbar2)+fembar2;
        fprintf('member number=');
        disp(i);
        fprintf('Global displacement matrix=\n');
        disp(delbar2);
        fprintf('Global end moment matrix=\n');
        disp(mbar);
        elseif i==3
        delbar3=deli;
        mbar=(Kg3*delbar3)+fembar3;
        fprintf('member number=');
        disp(i);
        fprintf('Global displacement matrix=\n');
        disp(delbar3);
        fprintf('Global end moment matrix=\n');
        disp(mbar);
    end
end
        
        

