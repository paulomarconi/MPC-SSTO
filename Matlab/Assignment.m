% ACS6116 
% Paulo Loma Marconi
% University of Sheffield, UK

clear variables; close all; clc;
path('ACS6116_mcode',path); % temporal
%%
mode = 'reg'; %regulation=reg, tracking=trk

%% Global variables
x0 = [0; 0; 0];
r = 0.3; % reference

model_d = 1; % model with disturbance, 1=yes, 0=no
d = 0; % init value of the perturbation value

N = 3;  % Horizon increase N and cost value J(k) is reduced
nk = 80; % simulation steps

rho = 1; % SSTO optimization

Qdesign = 1; Rdesign = .08; % Q and R for tunning
%% System continuos+discrete

Msys = 10; Dsys = 0.8; Rsys = 0.1; Tt = 0.5; Tg = 0.2;

Ac = [-Dsys/Msys       1/Msys        0   ;
          0            -1/Tt        1/Tt ;
     -1/(Rsys*Tg)        0         -1/Tg  ];
Bc = [0; 0; 1/Tg];

if model_d==1
    Bd = [-1/Msys; 0; 0;]; % Bd=disturbance matrix
else
    Bd = [0; 0; 0];
end

C = [50, 0, 0];
D = 0;
Ts = 0.1; %sampling time

% Continuos + disturbance
sys = idss(Ac,Bc,C,D,Bd,x0,0); 

% Converting to discrete model
sysd = c2d(sys,Ts);
A = sysd.A; 
B = sysd.B; 
C = sysd.C; 
D = sysd.D; 
Bd = -sysd.K; 

if model_d == 1
    Dd = 1;
else 
    Dd = 0;
end



%% System Dimensions
% n  = size(A, 1); % Number of State variables
% m  = size(B, 2); % Number of Input variables
% q = size(Bd, 2); % Number of Disturbance Inputs
% p  = size(C, 1); % Number of Output measurments 
% nPx = size(Px);   % Dimensions of State Constrains
% mPu = size(Pu);   % Dimensions of Input Constrains

%% Constraints
% Input constraints
u_max = 0.5; u_min = 0.5;
Pu = [1; -1];
qu = [u_max; u_min]; 
% State constraints
y_max = 0.5; y_min = 0.5;
Px = [eye(3).*C; -eye(3).*C];
qx = [ones(3,1)*y_max; ones(3,1)*y_min];
% Terminal State constraints
Pxf = Px;
qxf = qx;


%% Part II: Target equilibrium pairs and optimization under constraints
% Target
T = [eye(3)-A, -B; C 0];
Tinv = T\[Bd*d; r-Dd*d];   
xss = Tinv(1:3);
uss = Tinv(4);

% Optimization problem
x_optVect = [xss; uss]; % vector x for optimization problem A*x<b

fun = @(x_optVect)  ( C(1)*x_optVect(1)+Dd*d-r ) + rho*x_optVect(4) ;  
x0_opt = [r, 0, 0, 0]; % init values

% Equality constraints
Aeq = T; 
beq = [Bd; r-Dd*d];
% Inequality constraints
Aopt = [Px, zeros(6,1); zeros(2,3), Pu]; 
bopt = [qx; qu];

% Minimization 
[f,fval] = fmincon(fun,x0_opt,Aopt,bopt,Aeq,beq);
xss = f(1:3)';
uss = f(4);
yss = C*xss+Dd*d;

%% Part II: Constraints for SSTO model
% Input constraints
Pu_ssto = Pu;
qu_ssto = qu-Pu*uss; 
% State constraints
Px_ssto = Px; %
qx_ssto = qx-Px*xss;

% Cost matrices
Q = (C'*C)*Qdesign; % Q^(1/2)
R = 1*Rdesign;

% Deadbeat Terminal Inequality State constraints
% K  = -acker(A,B,[0,0,0]);  
K = -dlqr(A,B,Q,R);
% K = [K(1) 0 0];

Maux = [];
for i =0:N-1
    Maux = [Maux;(A+B*K)^(i)];
end 
Mm = kron(eye(N),[Px; Pu*K]);

% regulation
PxN = Mm*Maux; 
qxN = kron(ones(N,1),[qx; qu]);% deadbeat terminal inequality constraints
% tracking
PxN_ssto = Mm*Maux; 
qxN_ssto = kron(ones(N,1),[qx-Px*xss; qu-Pu*uss]);% deadbeat terminal inequality constraints


%% Prediction, Constraint and Cost matrices
% Prediction matrices
[F,G] = predict_mats(A,B,N);

% Constructing Constraint matrices
[Pc,qc,Sc] = constraint_mats(F,G,Pu,qu,Px,qx,PxN,qxN);
[Pc_ssto,qc_ssto,Sc_ssto] = constraint_mats(F,G,Pu_ssto,qu_ssto,Px_ssto,qx_ssto,PxN_ssto,qxN_ssto);

% % Cost matrices
% Q = C'*C; % Q^(1/2)
% R = 1;

%% Dead-beat controller K for mode-2 
% K = -acker(A,B,[0,0,0]);

% Closed-loop
Acl = A+B*K;
% Check if it's stable
eigVal = eig(Acl);
if abs(real(eigVal)) <= 1; disp('Mode-2 CL system x(k+1) = (A+B*K)*x(k) is stable!');    
else; disp('Mode-2 CL system x(k+1) = (A+B*K)*x(k) is unstable!'); end
disp('eigen values:'); disp(eigVal);

%% Calculate the terminal cost matrix P
Phi = A+B*K;
S = Q+K'*R*K;
P = dlyap(Phi',S);

% Cost matrices
[H,L,M] = cost_mats(F,G,Q,R,P);

%% Constrained - Closed-loop states (xs) and inputs (us) predictions
% initialize
x = x0;

n = 3; m = 3;
xs = zeros(n,nk+1);
us = zeros(m,nk+1);

for k = 1:nk+1
    
    % Function to call perturbation at specific time steps    
    if k == 30
        d = 0.6;
    elseif k == 60
        d = -0.2;
    else 
        d = 0; 
    end
    ref(1,k) = r; 
    
    
    if mode == 'trk'       
        % Optimization problem
        fun = @(x_optVect2)  ( C(1)*x_optVect2(1)+Dd*d-r ) + (rho*x_optVect2(4)) ;  
        x0_opt = [x;us(1,k)];

        % Equality constraints
        Aeq = [eye(3)-A, -B ;C 0]; 
        beq = [Bd*d; r-Dd*d];

        % Minimization 
        [f,fval] = fmincon(fun,x0_opt,Aopt,bopt,Aeq,beq);
        xss = f(1:3);
        uss = f(4);
        yss = C*xss+Dd*d;
        x_optVect2 = [xss;uss]; % vector x for optimization problem A*x<b 

        % SSTO states
        epsilon(:,k) = x-xss;
        nu(:,k) = us(:,k)-uss;

        % input constraints
        qu_ssto = qu-Pu*uss;
        % state constraints
        qx_ssto = qx-Px*xss;
        % deadbeat terminal inequality constraints
        qxN_ssto = kron(ones(N,1),[qx-Px*xss; qu-Pu*uss]);
        % Constraint matrices
        [Pc_ssto,qc_ssto,Sc_ssto] = constraint_mats(F,G,Pu_ssto,qu_ssto,Px_ssto,qx_ssto,PxN_ssto,qxN_ssto);
    end
   
    % States
    xs(:,k) = x;
    % Measurement
    y = C*xs+Dd*d;
%     y1(k,:) = C*x+Dd*d; % violates y constraint |y|<0.5
    
    
    % Constrained MPC control law (RH-FH) LQ-MPC at every step k        
    if mode == 'reg'
        [Ustar,fval,flag] = quadprog(H,L*x,Pc,qc+Sc*x); 
    elseif mode == 'trk'
        [NUstar,fval,flag] = quadprog(H,L*epsilon(:,k),Pc_ssto,qc_ssto+Sc_ssto*epsilon(:,k)); 
    end
    % check feasibility
    if flag < 1 
        disp(['Optimization is infeasible at k = ',num2str(k)]);
        break;    
    end
    
    % store input
    if mode == 'reg'
        %us_C(:,k) = Ustar(1:3,1);
        us(:,k) = Ustar(1);
    elseif mode == 'trk'
        us(:,k) = NUstar(1)+uss;    
    end   
    
    % Feedback x(k+1)=A*x(k)+B*u(k)
    x = A*x+B*us(1,k)+Bd*d;      
    
    
    
    % Cost Value 
    if mode == 'reg'
        J(:,k) = 0.5*(Ustar(1:N,:)')*H*Ustar(1:N,:)+(xs(:,k)')*(L')*Ustar(1:N,:)+...
                  (xs(:,k)')*M*xs(:,k);
    elseif mode == 'trk'
        J(:,k) = 0.5*(NUstar(1:N,:)')*H*NUstar(1:N,:)+(epsilon(:,k)')*(L')*NUstar(1:N,:)+...
                  (epsilon(:,k)')*M*epsilon(:,k);
    end
    
end
% Predicted states for ploting and for value function
xs_aux = xs';
xs1 = xs_aux(:,1)';
xs2 = xs_aux(:,2)';
xs3 = xs_aux(:,3)';

%% Plotting
t = Ts*(0:nk); % time simulation

figure(1); 
    subplot(3,3,1);
    hold on;
%     stairs(0:N-1,us_U); % unconstrained open-loop    
%     stairs(0:nk-1,us(1,1:end-1)'); % unconstrained closed-loop
    stairs(t,us(1,1:end)'); % constrained closed-loop    
    title(['Input u for N=',num2str(N)]);
    xlabel('Time step k [s]'); 
    ylabel('u(k)');
%     ylabel('u(0+j|0)/u(k)');
%     legend('u(0+j|0)','us(k)','us C(k)');
    grid on;

    subplot(3,3,2); 
    hold on;
%     plot(1:N,x1,1:N,x2); % unconstrained open-loop  
%     plot(1:nk,xs1,1:nk,xs2); % unconstrained closed-loop
    plot(t,xs1,t,xs2,t,xs3); % constrained closed-loop      
    title(['States for N=',num2str(N)]);
    xlabel('Time step k [s]'); 
    ylabel('States'); 
    legend('x1(k)','x2(k)','x3(k)','Location','Southeast');
    grid on;

    subplot(3,3,3); 
    hold on;
    plot(xs1,xs2,'-o');
%     plot3(xs1_C,xs2_C,xs3_C);
    title(['phase-plot for N=',num2str(N)]);
    xlabel('x1(k)'); 
    ylabel('x2(k)');
    grid on;
    
    
    
    subplot(3,3,4); 
    hold on;
    plot(xs1,xs3,'-o');
%     plot3(xs1_C,xs2_C,xs3_C);
    title(['phase-plot for N=',num2str(N)]);
    xlabel('x1(k)'); 
    ylabel('x3(k)');
    grid on;
    
    
    subplot(3,3,5); 
    hold on;
    plot(xs2,xs3,'-o');
%     plot3(xs1_C,xs2_C,xs3_C);
    title(['phase-plot for N=',num2str(N)]);
    xlabel('x2(k)'); 
    ylabel('x3(k)');
    grid on;
    
    subplot(3,3,6); 
    plot(t,y,t,ref);
    title(['Output y for N=',num2str(N)]);
    xlabel('Time step k [s]'); 
    ylabel('y(k)');    
    legend('y(k)','r(k)');
    grid on;

    subplot(3,3,7); 
    hold on;
    plot(t,J);
    title(['Value function for N=',num2str(N)]);
    xlabel('Time step k [s]'); 
    ylabel('J(k)');
    grid on;
   

%% Plot box state constraints
% x1_min=-y_min; x1_max=y_max;
% x2_min=-y_min; x2_max=y_max;
% x1 = [x1_min, x1_max, x1_max, x1_min, x1_min];
% x2 = [x2_min, x2_min, x2_max, x2_max, x2_min];
% plot(x1,x2,'b--');


%% Plot region
% x11 = -1:.1:2; x22 = -0.5:.1:3;
% [x1,x2] = meshgrid(x11,x22);  % create a grid
% scatter(x1(:),x2(:),2,'filled')

%% Plot feasible region
% v = -1:.01:1;
% [x1,x2] = meshgrid(v);  % create a grid
% ineq = x1 + x2 <= 0;    % some inequality
% f = double(ineq);
% surf(x1,x2,f);
% view(0,90)            % rotate surface plot to top view


% [x1,x2] = meshgrid( x(1),xs(2) );  % create a grid
% ineq = Pc*Ustar(1) <= qc+Sc*x;    % some inequality
% f = double(ineq);
% surf(x1,x2,f);
% view(0,90)      
