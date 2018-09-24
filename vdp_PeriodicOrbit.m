%% COCO Periodic orbit toolbox applied on the van der Pol oscillator

%% Parameters

mu = 2;
t_0 = 0;    % The beginning of the time interval [1]
t_f = 7;    % Time period [1]
x_10 = 1;   % The initial position [1]
x_20 = 0;   % The initial velocity [1]
   
%% Initial solution guess

pnames = 'mu';
p0 = mu;
x_0 = [x_10;x_20];
t_interval = [t_0;t_f];
[t0,x0] = ode45(@(t,x) vdp(x,p0),t_interval,x_0);

%% COCO

% Initialize an empty continuation problem structure named prob
prob = coco_prob();

% Setting the number of discretization interval (default value: 10)
prob = coco_set(prob, 'coll', 'NTST', 100);
% Setting the degree of the interpolating polynomials (default value: 4)
prob = coco_set(prob, 'coll', 'NCOL', 4);
% Number of grid points: N = NCOL*NTST+1

% Collecting the ode_isol2po() arguments
coll_args = {@vdp, t0, x0, pnames, p0};
% Encoding the periodic orbit continuation
prob = ode_isol2po(prob, '', coll_args{:});
% Assigning integer 1 to ensure that adaptive changes are made to the orbit discretization after each succesful step
prob = coco_set(prob, 'cont', 'NAdapt', 1);

% Calculation
coco(prob, 'run', [], 1, {'mu' 'po.period'}, [0 5]);

%% Evaluation

bd = coco_bd_read('run');  % Extract bifurcation data
labs = coco_bd_labs(bd);	% Extract labels

fig1 = figure(1);
hold on
grid on
for i_lab = labs
    sol = po_read_solution('','run',i_lab);
    plot(sol.xbp(:,1), sol.xbp(:,2),'LineWidth',1)
end
hold off

%% Functions

function dx_dt = vdp(x,p)
    
    % State variables
    x1 = x(1,:);
    x2 = x(2,:);
    
    % Parameters
    mu = p(1,:);
    
    % The ODE
    dx_dt = [x2; mu.*(1-x1.^2).*x2-x1];
    
end