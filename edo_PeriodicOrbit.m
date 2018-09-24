%% COCO Periodic orbit toolbox applied on the excited harmonic oscillator

%% Parameters

alpha = pi;         % Natural angular frequency [rad/s]
f_0 = 1;            % Static deformation [m]
omega = 2*pi;       % Angular frequency of the excitation [rad/s]
t_0 = 0;            % The beginning of the time interval [s]
T = 2*pi/omega;     % Time period [s]
x_10 = 0;           % The initial position [m]
x_20 = 0;           % The initial velocity [m/s]
D = 0;              % The relative damping [-]
   
%% Initial solution guess

pnames = {'D' 'alpha' 'f_0' 'omega'};
p0 = [D; alpha; f_0; omega];
x_0 = [x_10;x_20];
t_interval = [t_0;T];
[t0,x0] = ode45(@(t,x) edo(t,x,p0),t_interval,x_0);

%% COCO

% Initialize an empty continuation problem structure named prob
prob = coco_prob();

% Setting the ode to autonomous to indicate the explicit time dependence
prob = coco_set(prob, 'ode', 'autonomous', false);
% Setting the number of discretization interval (default value: 10)
prob = coco_set(prob, 'coll', 'NTST', 50);
% Setting the degree of the interpolating polynomials (default value: 4)
prob = coco_set(prob, 'coll', 'NCOL', 4);
% Number of grid points: N = NCOL*NTST+1

% Collecting the ode_isol2po() arguments
coll_args = {@edo, @edo_DFDX, @edo_DFDP, @edo_DFDT, t0, x0, pnames, p0};
% Encoding the periodic orbit continuation
prob = ode_isol2po(prob, '', coll_args{:});
% Assigning integer 1 to ensure that adaptive changes are made to the orbit discretization after each succesful step
prob = coco_set(prob, 'cont', 'NAdapt', 1);

% Calculation
coco(prob, 'run', [], 1, {'D' 'po.period'}, [-0.75 0.75]);

%% Evaluation

bd = coco_bd_read('run');  % Extract bifurcation data
labs = coco_bd_labs(bd);	% Extract labels

fig1 = figure(1);
hold on
grid on
for i_lab = 1:labs
    sol = po_read_solution('','run',i_lab);
    plot(sol.xbp(:,1), sol.xbp(:,2),'LineWidth',1)
end
hold off

%% Functions

function dx_dt = edo(t,x,p)
    
    % Parameters
    D = p(1,:);
    alpha = p(2,:);
    f_0 = p(3,:);
    omega = p(4,:);
    
    % State variables
    x1 = x(1,:);
    x2 = x(2,:);
    
    % The ODE
    dx_dt(1,:) = x2;
    dx_dt(2,:) = -alpha.^2.*x1 - 2*D.*alpha.*x2 + f_0.*alpha.^2.*cos(omega.*t);
    
end

function J = edo_DFDX(t,x,p)

    % State variables
    x1 = x(1,:);
    % x2 = x(2,:);
    
    % Parameters
    D = p(1,:);
    alpha = p(2,:);
    % f_0 = p(3,:);
    % omega = p(4,:);
    
    % The Jacobian matrix
    J = zeros(2,2,numel(x1));
    J(1,2,:) = 1;
    J(2,1,:) = -alpha.^2;
    J(2,2,:) = -2*D.*alpha;
    
end

function J = edo_DFDT(t,x,p)

    % State variables
    % x1 = x(1,:);
    % x2 = x(2,:);
    
    % Parameters
    % D = p(1,:);
    alpha = p(2,:);
    f_0 = p(3,:);
    omega = p(4,:);

    J = zeros(2,numel(t));
    J(2,:) = -f_0.*alpha.^2.*omega.*sin(omega.*t);
    
end

function J = edo_DFDP(t,x,p)

    % State variables
    x1 = x(1,:);
    x2 = x(2,:);
    
    % Parameters
    % D = p(1,:);
    alpha = p(2,:);
    f_0 = p(3,:);
    omega = p(4,:);
    
    % The Jacobian matrix
    J = zeros(2,4,numel(x1));
    J(2,1,:) = -2*alpha.*x2;
    J(2,2,:) = -2*alpha.*x1;
    J(2,3,:) = alpha.^2.*cos(omega.*t);
    J(2,4,:) = -f_0.*alpha.^2.*omega.*sin(omega.*t);
    
end