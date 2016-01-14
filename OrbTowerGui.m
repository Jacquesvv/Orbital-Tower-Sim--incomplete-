%Jacques Janse Van Vuuren
%12092369
%

% wolfram - (a-2*b+c)/(h^2)=(Z^2)*0.5*(((d-2*a+e)/(k^2))+((f-2*b+g)/(k^2)))
%
function varargout = OrbTowerGui(varargin)
% ORBTOWERGUI MATLAB code for OrbTowerGui.fig
%      ORBTOWERGUI, by itself, creates a new ORBTOWERGUI or raises the existing
%      singleton*.
%
%      H = ORBTOWERGUI returns the handle to a new ORBTOWERGUI or the handle to
%      the existing singleton*.
%
%      ORBTOWERGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ORBTOWERGUI.M with the given input arguments.
%
%      ORBTOWERGUI('Property','Value',...) creates a new ORBTOWERGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before OrbTowerGui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to OrbTowerGui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help OrbTowerGui

% Last Modified by GUIDE v2.5 02-Oct-2015 17:33:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @OrbTowerGui_OpeningFcn, ...
                   'gui_OutputFcn',  @OrbTowerGui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before OrbTowerGui is made visible.
function OrbTowerGui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to OrbTowerGui (see VARARGIN)

% Choose default command line output for OrbTowerGui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


%Constants
M = 6.4185*10^23;       %mass Mars
G = 6.6738480*10^-11;   %gravitational constant
f = 88643;              %rotational frequency - seconds
w = (2*pi)/f;           %angular velocity mars - rad/s
R = 3390000;            %radius mars - m

densCarbonNano = 1600;  %kg/m^3
diamTower = 10;         %m


geoStat =((G*M)/(w^2))^(1/3) - R

%syms x

%x=linspace(0,10000,100);

r = R:1000:80000000;
Aup = (w^2)*r;
Adown = (G*M)./(r.^2);

axes(handles.axes1);
plot(r, Aup,r,Adown);
legend('a-centripetal','a-gravity');


axis([0 80000000 0 1])
ylabel('Acceleration - m/s^2') % y-axis label
xlabel('Distance - m') % y-axis label
title('Normal','FontSize',12)

x = R:1000:80000000;
Fup =((w^2)*(x.^2))/2 - ((w^2)*(R^2))/2;
Fdown = -(G*M)./x + (G*M)/R;
Fdiff = (-(w^2)*(x.^2))/2 + ((w^2)*(R^2))/2 -(G*M)./x + (G*M)/R;

distTower = (sqrt(8*G*M*R*(w^2) + (R^4)*(w^4)))/(2*R*(w^2)) -R/2

axes(handles.axes2);
plot(x, Fdown,x,Fup,x,Fdiff);
legend('dF-gravity','dF-centripetal','dF-diff');
xlabel('Distance - m') % y-axis label
title('Integrated','FontSize',12)
axis([0, 80000000, 0, 15000000])


axes(handles.axes3);

mTower = distTower*(2*pi*diamTower)*densCarbonNano;

k = R:1000:80000000;
TensionAtPoint = mTower*(-((w^2)*(k.^2))/2 + ((w^2)*(R^2))/2 -(G*M)./k + (G*M)/R);

plot(k, TensionAtPoint);
ylabel('Force - N') % y-axis label
xlabel('Distance - m') % y-axis label
title('Tension','FontSize',12)
legend('Tension');
axis([0, 80000000, 0, 70000000000000000000])



nRows = 12;
nColumns = 11;



%u = zeros
%uI(1:nColumns)=1:nColumns
uI = zeros(nColumns,1)

u(1) =0;
u(2) =0;
for j = 2:10
%uI = zeros(nColumns,1);
uI(1:nColumns)=R:1000:(R+((nColumns-1)*1000))
 
%uI(1:(80000-3390))=3390001:1000:80000000;
%length(uI)

%uI(R:1000:80000000)=1

%uI(5)= uI(5)+500;

Tension = mTower*(-((w^2)*(uI.^2))/2 + ((w^2)*(R^2))/2 -(G*M)./uI + (G*M)/R)
    
%Tension(5)= Tension(5)+5;



A1 = zeros(nRows,nColumns);
A1(3:(nRows+1):(nRows-1)*nColumns)=(-Tension(1:10));
A1(1:(nRows+1):nRows*nColumns)=uI(1:11) + Tension(1:11);
A1(2:(nRows+1):(nRows)*nColumns)=(Tension(1:11))/2;
A1(1,:)=[]


B1 = zeros(nRows,nColumns);
B1(3:(nRows+1):(nRows-1)*nColumns)=(Tension(1:10))/2;
B1(1:(nRows+1):nRows*nColumns)=2*(uI(1:11)) + Tension(1:11);
B1(2:(nRows+1):(nRows)*nColumns)=(Tension(1:11))/2;
B1(1,:)=[];

C1 = zeros(nRows,nColumns);
C1(2:(nRows+1):(nRows)*nColumns)=-1*(uI(1:11));
C1(1,:)=[];



u(2) = u(2)+5;
u(j+1) = (((B1*uI)*(u(j)) + (C1*uI*(u(j-1)))) \ (A1*uI))

end

axes(handles.axes4);

plot(1:11,u);
ylabel('Displacement') % y-axis label
xlabel('Time') % y-axis label


%eq = (((w^2)*(x^2))/2) - (((w^2)*(R^2))/2) == ((G*M)/x) - ((G*M)/R);
%solx = solve(eq,x)

%eq = (sqrt(8*G*M*R*(w^2) + (R^4)*(w^4)))/(2*R*(w^2)) -R/2


%p = [(w^2)/2 0 (((G*M)/R)-((w^2)*(R^2))/2) G*M]

%root = roots(p);

%{

dist =700;
time = 1000;

i=(0:dist);
j =(0:time);


axes(handles.axes1);
n = 100;
X= (0:n)*0 +700;
Y= (0:n)*0 +700;
Z = (0:n)*20;


X(50:70) = 4.8;


plot3(X,Y,Z);
grid on
shading interp







%=========================================  Draw Disc   ================
hold on 
brown = [0 0 0.9];
x = 700;
y = 700;
z = 0;
h = 2000;

r = 10;
[Xcyl,Ycyl,Zcyl] = cylinder(r,30);
Xcyl = Xcyl + x;
Ycyl = Ycyl + y;
Zcyl = Zcyl*h + z;
Disc =surf(Xcyl,Ycyl,Zcyl,'FaceColor', brown,'EdgeColor','none');
set(Disc,'FaceAlpha',0.5)

%=======================================================================
axis([0 1000 0 1000 -50 2100]);
%set(gca,'Color',[0.18 0.18 0.18]);

%=========================================  Draw Sun    ================
%axes(handles.axes3);
hold on

rSun = 500; 
sunX =700;
sunY=700;
sunZ =-500;

phi=linspace(0,1.3,30);
theta=linspace(0,2*pi,50);
[phi,theta]=meshgrid(phi,theta);

x=rSun*sin(phi).*cos(theta);
y=rSun*sin(phi).*sin(theta);
z=rSun*cos(phi); 

faceC = [0.55 0.40 0.05];
edgeC = [0.8 0.5 0.1];

surf(x+sunX, y+sunY, z+sunZ,'FaceColor', faceC,'EdgeColor',edgeC);
%axis([-15000 25000 -15000 25000 -25000 25000]);
%set(gca,'Color',[0.18 0.18 0.18]);
%=======================================================================

%}

% UIWAIT makes OrbTowerGui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = OrbTowerGui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
