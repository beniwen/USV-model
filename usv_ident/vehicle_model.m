function [dstate,y] =vehicle_model(t,state, tau ,Xudot, Yvdot, Kpdot, Nrdot, varargin)

% 线性化船舶模型
% Inputs:
% eta= [u v p r phi psi x y]'
% xdot= [du dv dp dr dphi dpsi]'
% pos_dot=[dx dy]'
% tau = [Xe Ye Ke Ne]' 
%

% where
% u     = surge velocity          (m/s)
% v     = sway velocity           (m/s)
% p     = roll velocity           (rad/s)
% r     = yaw velocity            (rad/s)
% phi   = roll angle              (rad)
% psi   = yaw angle               (rad)
%
% Xe is the surge external force (e.g. rudder and fin forces)
% Ye is the sway external force  
% Ke is the sway external force  
% Ne is the sway external force   

% Vessel Data
rho_water     =	1014.0;	        %	water density	[kg/m^3]	
rho_air		=	1.225	;	    %	air density		[kg/m^3]	
g				=	9.81;	        %	gravity constant	[m/s^2]	
deg2rad 		=	pi/180;	        %	degrees to radians	
rad2deg 		=	180/pi;	        %	rad to degrees		
ms2kt			=	3600/1852;	    % 	m/s to kt 			
kt2ms 		=	1852/3600;	    %	kt to m/s			
RPM2rads		=	2*pi/60;	    %	RPM to rad/s		
rads2RPM		=   60/(2*pi);	    %	rad/s to RPM		
HP2W			=	745.700;	    %	HP to Watt			


Lpp    =  51.5 ;                  % Length between perpendiculars垂直间长 [m]
B      =  8.6  ;                  % Beam over all 船宽 [m]
D	   =  2.3  ;                  % Draught吃水深度 [m]     

%Load condition (Modified by T.Perez)
disp   =  350.0;                   % Displacement 排水量 [m^3]
m      =  disp*rho_water;          % Mass [Kg]
Izz    =  47.934*10^6 ;            % Yaw Inertia偏航转动惯量
Ixx    =  2.3763*10^6 ;            % Roll Inertia横滚转动惯量
U_nom  =  8.0   ;	               % Speed nominal [m/sec]速度 (app 15kts) 
KM		=  4.47;	             %  [m] Transverse metacentre above keel 横倾复原力臂
KB		=  1.53;	             %  [m] Transverse centre of bouancy 浮心高度
gm 		=  1.1;	                 %  [m]	Transverse Metacenter 
bm 		=  KM - KB;
LCG       = 20.41 ;                % [m]	Longitudinal CG (from AP considered at the rudder stock)
VCG       = 3.36  ;                % [m]	Vertical  CG  above baseline
xG        = -3.38  ;               % coordinate of CG from the body fixed frame adopted for the PMM test  
zG  	    = -(VCG-D);          % coordinate of CG from the body fixed frame adopted for the PMM test  
m_xg	    = m * xG;
m_zg	    = m * zG;
Dp        = 1.6 ;                   % Propeller diameter [m]


% Data for surge equation 纵向（均为线性流体动力导数）
% Xudot  	= -17400.0 ;  %-mx
Xuau     	= -1.96e+003 ;
Xvr    	=  0.33 * m ;
   
% Hydrodynamic coefficients in sway equation横向
% Yvdot = -393000 ;  %-my
Ypdot = -296000 ; 
Yrdot = -1400000 ; %附加质量静距
Yauv  = -11800 ; 
Yur   =  131000 ; 
Yvav  = -3700 ; 
Yrar   =  0 ;
Yvar  = -794000 ; 
Yrav  = -182000 ; 
Ybauv =  10800 ; % Y_{\phi |v u|}
Ybaur =  251000 ; 
Ybuu  = -74 ; 


% Hydrodynamic coefficients in roll equation横摇水动力系数
Kvdot =  296000 ;
% Kpdot = -774000 ;  %-Jxx
Krdot =  0 ;
Kauv  =  9260 ;
Kur   = -102000 ;
Kvav  =  29300 ;
Krar  =  0 ;
Kvar  =  621000 ;
Krav  =  142000 ;
Kbauv =  -8400 ;
Kbaur =  -196000 ;
Kbuu  =  -1180 ;
Kaup  =  -15500 ;
Kpap  =  -416000 ;
Kp    =  -500000 ;
Kb    =  0.776*m*g;
Kbbb  =  -0.325*m*g ;

% Hydrodynamic coefficients in yaw equation偏航水动力系数
Nvdot =  538000 ;   %-
Npdot =  0 ;
% Nrdot = -38.7e6;    %jzz
Nauv  = -92000 ;
Naur  = -4710000 ;
Nvav  =  0 ;
Nrar  = -202000000 ;
Nvar  =  0 ;
Nrav  = -15600000 ;
Nbauv = -214000 ;
Nbuar = -4980000 ;
Nbuau = -8000 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 定义输入
u   = state(1);
v   = state(2);
p  	= state(3);   	
r  	= state(4); 
b  	= state(5); %phi roll 		
psi = state(6);
x_pos  	= state(7);   	
y_pos  	= state(8); 
% 外部力与力矩
Xe  = tau(1);
Ye  = tau(2);
Ke  = tau(3);
Ne  = tau(4);

% Auxiliary variables辅助变量
au = abs(u);
av = abs(v); 
ar = abs(r); 
ap = abs(p); 
ab = abs(b);
L2 = Lpp^2; %参考面积？

% Total Mass Matrix 刚体质量矩阵
% x= [u v p r phi psi ]'
M =[ (m-Xudot)  0   0   0   0   0;
   0 (m-Yvdot) -(m*zG+Ypdot) (m*xG-Yrdot) 0 0;
   0 -(m*zG+Kvdot) (Ixx-Kpdot) -Krdot 0 0;
   0 (m*xG-Nvdot) -Npdot (Izz-Nrdot) 0 0;
   0 0 0 0 1 0; 
   0 0 0 0 0 1] ;
%转换矩阵
J=[cos(psi),-cos(b)*sin(psi);sin(psi),cos(b)*cos(psi)];
X=[u,v]';   
% Hydrodynamic forces without added mass terms(considered in the M matrix)
% 不考虑附加质量的水动力作用力 粘性流动动力（贵岛模型）
Xh  = Xuau*u*au+Xvr*v*r;

Yh = Yauv*au*v +  Yur*u*r + Yvav*v*av + Yvar*v*ar + Yrav*r*av + Ybauv*b*abs(u*v) + Ybaur*b*abs(u*r) + Ybuu*b*u^2;

Kh = Kauv*au*v +Kur*u*r + Kvav*v*av + Kvar*v*ar + Krav*r*av + Kbauv*b*abs(u*v) + Kbaur*b*abs(u*r) + Kbuu*b*u^2 +Kaup*au*p + Kpap*p*ap +Kp*p +Kbbb*b^3-(rho_water*g*gm*disp)*b;

Nh =Nauv*au*v + Naur*au*r +Nrar*r*ar + Nrav*r*av+Nbauv*b*abs(u*b) + Nbuar*b*u*ar + Nbuau*b*u*au;
 
% Rigid-body centripetal accelerations刚体运动产生的离心加速度影响
Xc =  m*(r*v+xG*r^2-zG*p*r);  
Yc = - m*u*r;
Kc =   m*zG*u*r;
Nc = - m*xG*u*r;

% Total forces
F1 = Xh+Xc+Xe;
F2 = Yh+Yc+Ye;
F4 = Kh+Kc+Ke;
F6 = Nh+Nc+Ne;

%运动学及动力学方程
% xdot= [du dv dp dr dphi dpsi]'
% pos_dot=[x y]'
% dstate=[xdot',pos_dot']';
%y=[u;v;p;r;b;psi;x_pos;y_pos];
pos_dot=J*X;
xdot = M\[F1; F2; F4; F6; p; r*cos(b)];
dstate=[xdot',pos_dot']';
y=[u;v;p;r;b;psi;x_pos;y_pos];
end