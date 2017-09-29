syms tx ty tz bgx bgy bgz ux uy uz nx ny nz wx wy wz w1
syms dtx dty dyz
global A C B W F V G ip op ux_real uy_real uz_real ux uy uz
P = zeros(5);
xend = [0]
xend1 = [0]
xend2 = [0]
xend3 = [0]
xend4 = [0]
int_ux = 0;
int_uy = 0;
traj_ux = [0];
traj_uy = [0];


vec1 = [(ux - bgx - nx);(uy - bgy - ny);(uz - bgz - nz)]    %THis is the input
vec2 = [1 sin(tx)*tan(ty) cos(tx)*tan(ty); 0 cos(tx) -sin(tx); 0 sin(tx)/cos(ty) cos(tx)/cos(ty)]  %This is the state equation
outputvec = [-9.8*sin(ty); 9.8*cos(ty)*sin(tx); 9.8*cos(ty)*cos(tx); tz]
vec3 = vec2*vec1

%vec3 = transpose(vec3); %We are using transpose here as the Jacobian formula later will use a row rather than column


func = [vec3;0;0;0]; %This is the function whose jacobian we find
state = [tx;ty;tz;bgx;bgy;bgz];
input = [ux;uy;uz];  %Remember everything is in ROW not column form, INPUT
distur = [nx;ny;nz];
%distur = [nx];

NonLinearFuncA = jacobian(func,state) %The jacobian is the derivative that we evaluate at a point to find A
NonLinearFuncB = jacobian(func,input) %This is B
NonLinearFuncF = jacobian(func,distur) %This is F
NonLinearFuncC = jacobian(outputvec,state) %This is C
V = 0.00015*eye(3)
W = 0.001*eye(4)
N = zeros(3)
D = zeros(3)
H = zeros(3)


%creating vec4 to measure real states
bgx= 0;
bgy = 0;
bgz = 0;
nx = 0;
ny = 0;
nz = 0;


tx = 0;
ty = 0;
tz = 0;
ux = 0 + 0.0002 + random('norm',0,0.00015);
uy = 0 + 0.0002 + random('norm',0,0.00015);
uz = 0.2 + random('norm',0,0.00015);
bgx = 0;
bgy = 0;
bgz = 0;
nx = 0; 
ny = 0;
nz = 0;
X = [tx;ty;tz;bgx;bgy;bgz]
X_ideal = [0;0;0]
ip = [ux;uy;uz]

wx = eval(outputvec(1)) + random('norm',0,0.001)
wy = eval(outputvec(2)) + random('norm',0,0.001)
wz = eval(outputvec(3)) + random('norm',0,0.001)
wz1 = eval(outputvec(4)) + random('norm',0,0.001)
op = [wx;wy;wz;wz1]
%Insert for loop here

%initial value of P
IniP(1:36) = 0
IniP1 = 0.1*eye(6)
IniP(1:36) = [IniP1(1,1:6) IniP1(2,1:6) IniP1(3,1:6) IniP1(4,1:6) IniP1(5,1:6) IniP1(6,1:6)]




for i = 1:0.01:50
A = eval(NonLinearFuncA); %This is A
B = eval(NonLinearFuncB); %This is B
F = eval(NonLinearFuncF); %This is F
C = eval(NonLinearFuncC); %This is C


[T1,P_tf] = ode45(@fare_b,[0 0.01],IniP);

IniP(1:36) = P_tf(end,1:36);


%Getting back the Matrix P

P(1,1:6) = P_tf(end,1:6);
P(2,1:6) = P_tf(end,7:12);
P(3,1:6) = P_tf(end,13:18);
P(4,1:6) = P_tf(end,19:24);
P(5,1:6) = P_tf(end,25:30);
P(6,1:6) = P_tf(end,31:36);



%Getting G
G = P*(C')/W;


%Getting back the states required
[T2,X_tf] = ode45(@p_estimate_b,[0 0.01],(X'));
X(1:6) = (X_tf(end,1:6))';

%Getting back the all the past parameters
tx = X(1);
ty = X(2);
tz = X(3);
bgx = X(4);
bgy = X(5);
bgz = X(6);


%plotting past parameters
xend = [xend X(1)];
xend1 = [xend1 X(2)];
xend2 = [xend1 X(3)];
xend3 = [xend1 X(5)];

%plotting ideal trajectories
ux_real = 10*cos(2*pi*2*i); %for plotting real trajectories
uy_real = 10*cos(2*pi*2*i);
uz_real = 0;
[T3,Xideal_tf] = ode45(@realp_estimate_b,[0 0.01],(X_ideal'));
X_ideal(1:3) = (Xideal_tf(end,1:3))';
%X_ideal(1) = X_ideal(1) -0.2;
%X_ideal(2) = X_ideal(2) -0.2;



%input trajectories
ux = 10*cos(2*pi*2*i) + 0.2 + random('norm',0,0.00015);
uy = 10*cos(2*pi*2*i) + 0.2 + random('norm',0,0.00015);
uz = 0.2 + random('norm',0,0.00015);

%int_ux = int_ux + ux;
%int_uy = int_uy + uy;
traj_ux = [traj_ux X_ideal(1)];
traj_uy = [traj_uy X_ideal(2)];
traj_uz = [traj_uy X_ideal(3)];
ip = [ux;uy;uz];





%output trajectories
wx = eval(outputvec(1)) + random('norm',0,0.001);
wy = eval(outputvec(2)) + random('norm',0,0.001);
wz = eval(outputvec(3)) + random('norm',0,0.001);
op = [wx;wy;wz;wz1];
end