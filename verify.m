% =======================================================================
% Part of the Research for: "On McDanell’s Conjecture"
% If you use this code, please cite the corresponding paper.
% Author: Mirhan Urkmez
% Date: 24/11/2024
% Description: 
% This script is a companion to the paper "Junction Points in Singular 
% Optimal Control" and verifies the counterexample values provided therein. 
% =======================================================================
%% The values generated with this code can alternatively be loaded  from "results.mat".

%%
clc; clear;

%% Check q(Check if beta is nonzero in the singular control interval.),r,m. GLC condition, range of u_s and u_n
%% System Definitions
n=3;
syms x [n 1] real
syms l [n 1] real
syms d [n 1] real
syms t u  real 
syms xi(t) eta(t) zeta(t) u_s(t) 

%Input on the nonsingular side
u_n=1;

% Dynamics
f=[2+3*x1+(19/12)*x1^4;1;1+x3-2*x3^2];
g=[-1-x1-(1/4)*x1^2-(8/3)*x1^3-(197/120)*x1^5+(7871/720)*x1^6;-2*x2+x2^2+(37/3)*x2^3;1+2*x3-3*x3^3];
xdot=f+g*u;
%cost
J=(1-u)*((x1-xi)+(x2-eta)+(x3-zeta))^2;

%Derivative of \xi, \eta, \zeta
ddot=subs(subs(xdot,u,u_s),x,d);
udot_n=diff(u_n,t);
udot_s=diff(u_s,t);

assumeAlso(xi(t),'real')
assumeAlso(eta(t),'real')
assumeAlso(zeta(t),'real')
assumeAlso(u_s(t),'real')

%Values at the junction point
lj=[1 1 1]'; %\lambda values at the junction
xj=[0 0 0]'; % state values at the junction
tj=0; % time at the junction
u_sj=1; %singular input at the junction
%% Hamiltonian calculations
H=l'*xdot+J(t);
ldot=-jacobian(H,x)';
H_u=diff(H,u); %Switching Function
q=1;
k=6; %Number of derivatives of \alpha and \beta to compute

%Taking time derivatives of the switching function
h_u_dot{1}=simplify(subs(diff(H_u,t),[diff(xi(t),t) diff(eta(t),t) diff(zeta(t),t)],ddot')+jacobian(H_u,x)*xdot+jacobian(H_u,l)*ldot+jacobian(H_u,d)*ddot);
for i=2:2*q
    h_u_dot{i}=simplify(subs(diff(h_u_dot{i-1},t),[diff(xi(t),t) diff(eta(t),t) diff(zeta(t),t)],ddot')+jacobian(h_u_dot{i-1},x)*xdot+jacobian(h_u_dot{i-1},l)*ldot+jacobian(h_u_dot{i-1},d)*ddot);
end
%\beta, \alpha and their derivatives on both nonsingular and singular arc
beta_n{1}=diff(h_u_dot{2*q},u);
alpha_n{1}=h_u_dot{2*q}-beta_n{1}*u;
alpha_s{1}=alpha_n{1};
beta_s{1}=beta_n{1};
for i=2:k+1
    alpha_n{i}=subs(diff(alpha_n{i-1},t),[diff(xi(t),t) diff(eta(t),t) diff(zeta(t),t)],ddot')+jacobian(alpha_n{i-1},x)*xdot+jacobian(alpha_n{i-1},l)*ldot+jacobian(alpha_n{i-1},u)*udot_n+jacobian(alpha_n{i-1},d)*ddot;
    beta_n{i}=subs(diff(beta_n{i-1},t),[diff(xi(t),t) diff(eta(t),t) diff(zeta(t),t)],ddot')+jacobian(beta_n{i-1},x)*xdot+jacobian(beta_n{i-1},l)*ldot+jacobian(beta_n{i-1},u)*udot_n+jacobian(beta_n{i-1},d)*ddot;
    alpha_s{i}=subs(diff(alpha_s{i-1},t),[diff(xi(t),t) diff(eta(t),t) diff(zeta(t),t)],ddot')+jacobian(alpha_s{i-1},x)*xdot+jacobian(alpha_s{i-1},l)*ldot+jacobian(alpha_s{i-1},u)*udot_s(t)+jacobian(alpha_s{i-1},d)*ddot;
    beta_s{i}=subs(diff(beta_s{i-1},t),[diff(xi(t),t) diff(eta(t),t) diff(zeta(t),t)],ddot')+jacobian(beta_s{i-1},x)*xdot+jacobian(beta_s{i-1},l)*ldot+jacobian(beta_s{i-1},u)*udot_s(t)+jacobian(beta_s{i-1},d)*ddot;
end

for i=1:k+1
    alpha0_n{i}=subs(subs(subs(alpha_n{i},[u d'],[ u_n x']),[t  ],tj),[ x' xi(0) eta(0) zeta(0) l'], [xj' xj' lj']);
    alpha0_s{i}=subs(subs(subs(alpha_s{i},[u d'],[ u_s x']),[t  ],tj),[ x' xi(0) eta(0) zeta(0) l'], [xj' xj' lj']);

    beta0_n{i}=subs(subs(subs(beta_n{i},[u d'],[ u_n x']),[t ],tj),[ x' xi(0) eta(0) zeta(0) l'], [xj' xj' lj']);
    beta0_s{i}=subs(subs(subs(beta_s{i},[u d'],[ u_s x']),[t  ],tj),[ x' xi(0) eta(0) zeta(0) l'], [xj' xj' lj']);
end

%% Finding u_s(0) and u_s'(0)
fprintf("The way we find u_s(0) and u_s'(0) is straightforward. \n We apply L'Hopital to u_s(t)=-alpha(t)/beta(t) at t=0. \n At each L'Hopital step we consider both the inconclusive and conclusive cases. \n If u_s(0) and u_s'(0) values can be found such tha L'Hopital step is conclusive and compatible with the u_s(0) value,then that constitutes one solution. \n Alternatively, if there exist values for u_s(0) and u_s'(0) that make the L'Hôpital step inconclusive (i.e., we must continue applying L'Hôpital), \n we verify that the first conclusive step ultimately confirms the initially assumed value of u_s(0).\n\n");

syms alpha(t) beta(t)

usClosed=-alpha(t)/beta(t); %Closed form equation of u_s

%See that both alpha(0) and beta(0) are 0
fprintf("alpha(0)=%.2f, beta(0)=%.2f. Therefore the limit of u_s(t)=-alpha(t)/beta(t) at t=0 \n will be determined using l'Hopital's rule. \n \n ",double(alpha0_s{1}),double(beta0_s{1}));


%Apply l'hopital

[numer,denom]=numden(usClosed);

% Cells to hold derivatives of numerator and denomiator of u_s
dernumer{1}=diff(numer,t);
derdenom{1}=diff(denom,t);

%Derivatives of numerator and denomiator with their values substitued
dernumerSubbed{1}=subs(dernumer{1},t,0);
derdenomSubbed{1}=subs(derdenom{1},t,0);
dernumerSubbed{1}=simplify(subs(dernumerSubbed{1},[alpha(0) beta(0) subs(diff(alpha(t), t), t, 0) subs(diff(beta(t), t), t, 0)],[alpha0_s{1} beta0_s{1} alpha0_s{2} beta0_s{2}]));
derdenomSubbed{1}=simplify(subs(derdenomSubbed{1},[alpha(0) beta(0) subs(diff(alpha(t), t), t, 0) subs(diff(beta(t), t), t, 0)],[alpha0_s{1} beta0_s{1} alpha0_s{2} beta0_s{2}]));
fprintf("First Iteration of L'Hopital: alpha'(0)/beta'(0)=(%s)/(%s)=%s. \n If the test is conclusive u_s(0)=1. For this to be conclusive \n at least one of the  alpha'(0) and beta'(0) should be nonzero. \n However, for u_s(0)=1, both are zero. Therefore, the test is inconclusive. \n Note that u_s(0) must still be 1 to make both alpha'(0) and beta'(0) 0. \n Then, the first conclusive l'Hopital step must make u_s(0)=1. \n\n",char(dernumerSubbed{1}),char(derdenomSubbed{1}),char(dernumerSubbed{1}/derdenomSubbed{1}));

u0value=double(dernumerSubbed{1}/derdenomSubbed{1});

%L'Hopital continues
for i=2:3
    dernumer{i}=diff(dernumer{i-1},t);
    derdenom{i}=diff(derdenom{i-1},t);
    dernumerSubbed{i}=subs(dernumer{i},t,0);
    derdenomSubbed{i}=subs(derdenom{i},t,0);
    %Put alpha(0) and beta(0) values
    dernumerSubbed{i}=simplify(subs(dernumerSubbed{i},[alpha(0) beta(0) subs(diff(alpha(t), t), t, 0) subs(diff(beta(t), t), t, 0) subs(diff(alpha(t), t, t), t, 0) subs(diff(beta(t), t, t), t, 0) subs(diff(alpha(t), t, t, t), t, 0) subs(diff(beta(t), t, t, t), t, 0)],[alpha0_s{1} beta0_s{1} alpha0_s{2} beta0_s{2} alpha0_s{3} beta0_s{3}  alpha0_s{4} beta0_s{4}]));
    derdenomSubbed{i}=simplify(subs(derdenomSubbed{i},[alpha(0) beta(0) subs(diff(alpha(t), t), t, 0) subs(diff(beta(t), t), t, 0) subs(diff(alpha(t), t, t), t, 0) subs(diff(beta(t), t, t), t, 0) subs(diff(alpha(t), t, t, t), t, 0) subs(diff(beta(t), t, t, t), t, 0)],[alpha0_s{1} beta0_s{1} alpha0_s{2} beta0_s{2} alpha0_s{3} beta0_s{3} alpha0_s{4} beta0_s{4}]));
end

%Check L'Hopital steps

syms udummy uderdummy
%Find values that make both alpha'' and beta'' 0.
i=2;
u_zero_solutions=solve([subs(dernumerSubbed{i},[u_s(0),subs(diff(u_s(t), t), t, 0)],[udummy, uderdummy])==0, subs(derdenomSubbed{i},[u_s(0), subs(diff(u_s(t), t), t, 0)],[udummy, uderdummy])==0],[udummy,uderdummy]);
% There is only one case where u_s(0)=1. 
%Then, we plug these values into alpha'''(0)/beta'''(0) below.



%Find the solution where u_s(0)=1
ind=find(double(u_zero_solutions.udummy)==1);

%check alpha'''(0)/beta'''(0) with that solution, check if it makes
%u_s(0)=1. It can be seen that u_s(0) is not equal to 1. 
%  Therefore, the values that make alpha''(0)/beta''(0)=0/0 are not
%  compatible with u_s(0)=0.
%Therefore, the  L'Hoptial has to be conclusive in this step.
us_check=subs(dernumerSubbed{3},[u_s(0),subs(diff(u_s(t), t), t, 0)],[u_zero_solutions.udummy(ind), u_zero_solutions.uderdummy(ind)])/subs(derdenomSubbed{3},[u_s(0),subs(diff(u_s(t), t), t, 0)],[u_zero_solutions.udummy(ind), u_zero_solutions.uderdummy(ind)]);
fprintf("Second Iteration of L'Hopital: There is only one u_s'(0) value \n such that u_s(0)= 0 and this step is inconclusive for alpha(0)''/beta''(0). \n We plug this value to alpha(0)'''/beta'''(0) which gives alpha(0)'''/beta'''(0)=%s. \n This can not be 1 or 0/0. Therefore this u_s(0)' value is not feasible. \n Then L'Hopital must be conclusive for alpha(0)''/beta''(0). \n\n",us_check);

% Find values that make L'Hopital conclusive such that
% u_s(0)==alpha''(0)/beta''(0)=1.
u_conclusive_solutions=solve([subs(dernumerSubbed{i},[u_s(0),subs(diff(u_s(t), t), t, 0)],[udummy, uderdummy])/ subs(derdenomSubbed{i},[u_s(0), subs(diff(u_s(t), t), t, 0)],[udummy, uderdummy])==udummy],[udummy,uderdummy]);
ind_conc=find(double(u_zero_solutions.udummy)==1);
fprintf("Now we check the  u_s'(0) values that makes alpha''(0)/beta''(0)=u_s(0)=1 conclusive. \n There is only one such u'(0)=%s \n\n",u_conclusive_solutions.uderdummy(ind_conc));
fprintf("Therefore u_s(0)=1, u_s'(0)=0. \n\n")

udot0value= double(u_conclusive_solutions.uderdummy(ind_conc));

%% Verify that the found u_s(0), u_s'(0) values comply when L'Hopital is applied to u_s'(0).
udot = diff(usClosed,t);
fprintf("Now, we have found u_s(0) and u_s'(0).Let's apply L'Hopital to u_s'(t) at t=0 \n and see we get the u_s'(0) value with the found u_s(0) and u_s'(0) values.\n\n")

%Apply l'hopital
[numudot,denudot]=numden(udot);

% Cells to hold derivatives of numerator and denomiator of u_s'
dernumudot{1}=diff(numudot,t);
derdenudot{1}=diff(denudot,t);

%Derivatives of numerator and denomiator with their values substitued
dernumudotSubbed{1}=subs(dernumudot{1},t,0);
derdenudotSubbed{1}=subs(derdenudot{1},t,0);
dernumudotSubbed{1}=simplify(subs(dernumudotSubbed{1},[u_s(0), subs(diff(u_s(t), t), t, 0) alpha(0) beta(0) subs(diff(alpha(t), t), t, 0) subs(diff(beta(t), t), t, 0)],[u0value udot0value alpha0_s{1} beta0_s{1} alpha0_s{2} beta0_s{2}]));
derdenudotSubbed{1}=simplify(subs(derdenudotSubbed{1},[u_s(0), subs(diff(u_s(t), t), t, 0) alpha(0) beta(0) subs(diff(alpha(t), t), t, 0) subs(diff(beta(t), t), t, 0)],[u0value udot0value alpha0_s{1} beta0_s{1} alpha0_s{2} beta0_s{2}]));

udotlhopitalnumerator=dernumudotSubbed{1};
udotlhopitaldenominator=derdenudotSubbed{1};
i=2;
while udotlhopitalnumerator==0 && udotlhopitaldenominator==0
    dernumudot{i}=diff(dernumudot{i-1},t);
    derdenudot{i}=diff(derdenudot{i-1},t);
    dernumudotSubbed{i}=subs(dernumudot{i},t,0);
    derdenudotSubbed{i}=subs(derdenudot{i},t,0);
    dernumudotSubbed{i}=simplify(subs(dernumudotSubbed{i},[alpha(0) beta(0) subs(diff(alpha(t), t), t, 0) subs(diff(beta(t), t), t, 0) subs(diff(alpha(t), t, t), t, 0) subs(diff(beta(t), t, t), t, 0) subs(diff(alpha(t), t, t, t), t, 0) subs(diff(beta(t), t, t, t), t, 0) subs(diff(alpha(t), t, t, t, t), t, 0) subs(diff(beta(t), t, t, t, t), t, 0)],[alpha0_s{1} beta0_s{1} alpha0_s{2} beta0_s{2} alpha0_s{3} beta0_s{3} alpha0_s{4} beta0_s{4} alpha0_s{5} beta0_s{5}]));
    derdenudotSubbed{i}=simplify(subs(derdenudotSubbed{i},[alpha(0) beta(0) subs(diff(alpha(t), t), t, 0) subs(diff(beta(t), t), t, 0) subs(diff(alpha(t), t, t), t, 0) subs(diff(beta(t), t, t), t, 0) subs(diff(alpha(t), t, t, t), t, 0) subs(diff(beta(t), t, t, t), t, 0) subs(diff(alpha(t), t, t, t,t), t, 0) subs(diff(beta(t), t, t, t, t), t, 0)],[alpha0_s{1} beta0_s{1} alpha0_s{2} beta0_s{2} alpha0_s{3} beta0_s{3} alpha0_s{4} beta0_s{4} alpha0_s{5} beta0_s{5}]));
    dernumudotSubbed{i}=subs(dernumudotSubbed{i},[u_s(0) subs(diff(u_s(t), t), t, 0)],[u0value udot0value ]);
    derdenudotSubbed{i}=subs(derdenudotSubbed{i},[u_s(0) subs(diff(u_s(t), t), t, 0)],[u0value udot0value ]);
    udotlhopitalnumerator=dernumudotSubbed{i};
    udotlhopitaldenominator=derdenudotSubbed{i};
    i=i+1;
end
fprintf("We have found that u_s(0)=1 and u_s'(0)=0. And when  apply L'Hopital to u_s'(t) at t=0 \n and we put these values, we get u_s'(0)=%s. \n\n",udotlhopitalnumerator/udotlhopitaldenominator);

%% Find u_s''(0)
fprintf("Now, we have found u'(0). The way we find u''(0) \n is by applying L'Hopital to u''(t) at t=0 \n and find the u''(0) by putting the found u_s(0) and u_s'(0) values.\n\n")

uddot = diff(udot,t);

%Apply l'hopital
[numuddot,denuddot]=numden(uddot);

% Cells to hold derivatives of numerator and denomiator of u_s'
dernumuddot{1}=diff(numuddot,t);
derdenuddot{1}=diff(denuddot,t);

%Derivatives of numerator and denomiator with their values substitued
dernumuddotSubbed{1}=subs(dernumuddot{1},t,0);
derdenuddotSubbed{1}=subs(derdenuddot{1},t,0);
dernumuddotSubbed{1}=simplify(subs(dernumuddotSubbed{1},[u_s(0), subs(diff(u_s(t), t), t, 0) alpha(0) beta(0) subs(diff(alpha(t), t), t, 0) subs(diff(beta(t), t), t, 0)],[u0value udot0value alpha0_s{1} beta0_s{1} alpha0_s{2} beta0_s{2}]));
derdenuddotSubbed{1}=simplify(subs(derdenuddotSubbed{1},[u_s(0), subs(diff(u_s(t), t), t, 0) alpha(0) beta(0) subs(diff(alpha(t), t), t, 0) subs(diff(beta(t), t), t, 0)],[u0value udot0value alpha0_s{1} beta0_s{1} alpha0_s{2} beta0_s{2}]));

uddotlhopitalnumerator=dernumuddotSubbed{1};
uddotlhopitaldenominator=derdenuddotSubbed{1};
i=2;
while uddotlhopitalnumerator==0 && uddotlhopitaldenominator==0
    dernumuddot{i}=diff(dernumuddot{i-1},t);
    derdenuddot{i}=diff(derdenuddot{i-1},t);
    dernumuddotSubbed{i}=subs(dernumuddot{i},t,0);
    derdenuddotSubbed{i}=subs(derdenuddot{i},t,0);
    dernumuddotSubbed{i}=simplify(subs(dernumuddotSubbed{i},[alpha(0) beta(0) subs(diff(alpha(t), t), t, 0) subs(diff(beta(t), t), t, 0) subs(diff(alpha(t), t, t), t, 0) subs(diff(beta(t), t, t), t, 0) subs(diff(alpha(t), t, t, t), t, 0) subs(diff(beta(t), t, t, t), t, 0) subs(diff(alpha(t), t, t, t, t), t, 0) subs(diff(beta(t), t, t, t, t), t, 0)],[alpha0_s{1} beta0_s{1} alpha0_s{2} beta0_s{2} alpha0_s{3} beta0_s{3} alpha0_s{4} beta0_s{4} alpha0_s{5} beta0_s{5}]));
    derdenuddotSubbed{i}=simplify(subs(derdenuddotSubbed{i},[alpha(0) beta(0) subs(diff(alpha(t), t), t, 0) subs(diff(beta(t), t), t, 0) subs(diff(alpha(t), t, t), t, 0) subs(diff(beta(t), t, t), t, 0) subs(diff(alpha(t), t, t, t), t, 0) subs(diff(beta(t), t, t, t), t, 0) subs(diff(alpha(t), t, t, t,t), t, 0) subs(diff(beta(t), t, t, t, t), t, 0)],[alpha0_s{1} beta0_s{1} alpha0_s{2} beta0_s{2} alpha0_s{3} beta0_s{3} alpha0_s{4} beta0_s{4} alpha0_s{5} beta0_s{5}]));
    dernumuddotSubbed{i}=subs(dernumuddotSubbed{i},[u_s(0) subs(diff(u_s(t), t), t, 0)],[u0value udot0value ]);
    derdenuddotSubbed{i}=subs(derdenuddotSubbed{i},[u_s(0) subs(diff(u_s(t), t), t, 0)],[u0value udot0value ]);
    uddotlhopitalnumerator=dernumuddotSubbed{i};
    uddotlhopitaldenominator=derdenuddotSubbed{i};
    i=i+1;
end
fprintf("After L'Hopital, we get \n u_s''(0)=(%s)/(%s), \n where subs(diff(u_s(t), t, t), t, 0)=u_s''(0).  That means u_s''(0)=0.\n\n",uddotlhopitalnumerator,uddotlhopitaldenominator);
syms uddotdummy
uddot0value=solve(subs(uddotlhopitalnumerator/uddotlhopitaldenominator,subs(diff(u_s(t), t, t), t, 0),uddotdummy)==uddotdummy);

%% Find u_s'''(0)
fprintf("Now we find u_s'''(0) in the same way as before.  We apply  L'Hopital \n to u_s'''(t) and put the found u_s(0), u_s'(0), u_s''(0) values. \n\n ");

udddot = diff(uddot,t);

%Apply l'hopital
[numudddot,denudddot]=numden(udddot);

% Cells to hold derivatives of numerator and denomiator of u_s'
dernumudddot{1}=diff(numudddot,t);
derdenudddot{1}=diff(denudddot,t);

%Derivatives of numerator and denomiator with their values substitued
dernumudddotSubbed{1}=subs(dernumudddot{1},t,0);
derdenudddotSubbed{1}=subs(derdenudddot{1},t,0);
dernumudddotSubbed{1}=simplify(subs(dernumudddotSubbed{1},[u_s(0), subs(diff(u_s(t), t), t, 0) subs(diff(u_s(t), t, t), t, 0) alpha(0) beta(0) subs(diff(alpha(t), t), t, 0) subs(diff(beta(t), t), t, 0)],[u0value udot0value uddot0value alpha0_s{1} beta0_s{1} alpha0_s{2} beta0_s{2}]));
derdenudddotSubbed{1}=simplify(subs(derdenudddotSubbed{1},[u_s(0), subs(diff(u_s(t), t), t, 0) subs(diff(u_s(t), t, t), t, 0) alpha(0) beta(0) subs(diff(alpha(t), t), t, 0) subs(diff(beta(t), t), t, 0)],[u0value udot0value uddot0value alpha0_s{1} beta0_s{1} alpha0_s{2} beta0_s{2}]));

udddotlhopitalnumerator=dernumudddotSubbed{1};
udddotlhopitaldenominator=derdenudddotSubbed{1};
i=2;
while udddotlhopitalnumerator==0 && udddotlhopitaldenominator==0
    dernumudddot{i}=diff(dernumudddot{i-1},t);
    derdenudddot{i}=diff(derdenudddot{i-1},t);
    dernumudddotSubbed{i}=subs(dernumudddot{i},t,0);
    derdenudddotSubbed{i}=subs(derdenudddot{i},t,0);
    dernumudddotSubbed{i}=simplify(subs(dernumudddotSubbed{i},[alpha(0) beta(0) subs(diff(alpha(t), t), t, 0) subs(diff(beta(t), t), t, 0) subs(diff(alpha(t), t, t), t, 0) subs(diff(beta(t), t, t), t, 0) subs(diff(alpha(t), t, t, t), t, 0) subs(diff(beta(t), t, t, t), t, 0) subs(diff(alpha(t), t, t, t, t), t, 0) subs(diff(beta(t), t, t, t, t), t, 0) subs(diff(alpha(t), t, t, t, t, t), t, 0) subs(diff(beta(t), t, t, t, t, t), t, 0) subs(diff(alpha(t), t, t, t, t, t, t), t, 0) subs(diff(beta(t), t, t, t, t, t, t), t, 0)],[alpha0_s{1} beta0_s{1} alpha0_s{2} beta0_s{2} alpha0_s{3} beta0_s{3} alpha0_s{4} beta0_s{4} alpha0_s{5} beta0_s{5} alpha0_s{6} beta0_s{6} alpha0_s{7} beta0_s{7}]));
    derdenudddotSubbed{i}=simplify(subs(derdenudddotSubbed{i},[alpha(0) beta(0) subs(diff(alpha(t), t), t, 0) subs(diff(beta(t), t), t, 0) subs(diff(alpha(t), t, t), t, 0) subs(diff(beta(t), t, t), t, 0) subs(diff(alpha(t), t, t, t), t, 0) subs(diff(beta(t), t, t, t), t, 0) subs(diff(alpha(t), t, t, t, t), t, 0) subs(diff(beta(t), t, t, t, t), t, 0) subs(diff(alpha(t), t, t, t, t, t), t, 0) subs(diff(beta(t), t, t, t, t, t), t, 0) subs(diff(alpha(t), t, t, t, t, t, t), t, 0) subs(diff(beta(t), t, t, t, t, t, t), t, 0)],[alpha0_s{1} beta0_s{1} alpha0_s{2} beta0_s{2} alpha0_s{3} beta0_s{3} alpha0_s{4} beta0_s{4} alpha0_s{5} beta0_s{5} alpha0_s{6} beta0_s{6} alpha0_s{7} beta0_s{7}]));
    dernumudddotSubbed{i}=subs(dernumudddotSubbed{i},[u_s(0) subs(diff(u_s(t), t), t, 0) subs(diff(u_s(t), t, t), t, 0)],[u0value udot0value uddot0value]);
    derdenudddotSubbed{i}=subs(derdenudddotSubbed{i},[u_s(0) subs(diff(u_s(t), t), t, 0) subs(diff(u_s(t), t, t), t, 0)],[u0value udot0value uddot0value]);
    udddotlhopitalnumerator=dernumudddotSubbed{i};
    udddotlhopitaldenominator=derdenudddotSubbed{i};
    i=i+1;
end

% Find u_s'''(0) value that complies with L'Hopital.
syms udddotdummy
udddot0value=solve(subs(udddotlhopitalnumerator/udddotlhopitaldenominator,subs(diff(u_s(t), t, t, t), t, 0),udddotdummy)==udddotdummy,udddotdummy);
fprintf("The value of u_s'''(0) that complies with L'Hopital is %s, \n which is the same as given in the paper. \n\n",udddot0value);

%% Verify alpha and beta values
fprintf("Now, we verify the values of alpha(t),alpha'(t),alpha''(t) and \n beta(t),beta'(t),beta''(t) using the found values of u_s(0), u_s'(0), u_s''(0).\n\n");

%Find values of alpha, beta
alpha0value=subs(alpha0_s{1},[u_s(0) subs(diff(u_s(t), t), t, 0) subs(diff(u_s(t), t, t), t, 0)],[u0value udot0value uddot0value]);
dalpha0value=subs(alpha0_s{2},[u_s(0) subs(diff(u_s(t), t), t, 0) subs(diff(u_s(t), t, t), t, 0)],[u0value udot0value uddot0value]);
ddalpha0value=subs(alpha0_s{3},[u_s(0) subs(diff(u_s(t), t), t, 0) subs(diff(u_s(t), t, t), t, 0)],[u0value udot0value uddot0value]);

beta0value=subs(beta0_s{1},[u_s(0) subs(diff(u_s(t), t), t, 0) subs(diff(u_s(t), t, t), t, 0)],[u0value udot0value uddot0value]);
dbeta0value=subs(beta0_s{2},[u_s(0) subs(diff(u_s(t), t), t, 0) subs(diff(u_s(t), t, t), t, 0)],[u0value udot0value uddot0value]);
ddbeta0value=subs(beta0_s{3},[u_s(0) subs(diff(u_s(t), t), t, 0) subs(diff(u_s(t), t, t), t, 0)],[u0value udot0value uddot0value]);

fprintf("alpha values: alpha(0)=%s, alpha'(0)=%s, alpha''(0)=%s. \n\n",alpha0value,dalpha0value,ddalpha0value);
fprintf("beta values: beta(0)=%s, beta'(0)=%s, beta''(0)=%s. \n\n",beta0value,dbeta0value,ddbeta0value);
fprintf("All the found values are the same as the ones in the paper. \n\n");