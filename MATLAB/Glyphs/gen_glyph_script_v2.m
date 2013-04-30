%% Make a Reynolds stress tensor glyph
clear all; clc; close all;
% Rotation
theta = 45;

%glyph to plot
gType=1;

switch gType
    case 1
        % 1-comp
        uiuj=[ 1   0   0;...
               0   0   0;...
               0   0   0];
    case 2
        % 2-comp
        uiuj=[ 0.5   0   0;...
               0   0   0;...
               0   0   0.5];
    case 3
        % 3-comp
        uiuj=[ 1/3   0   0;...
               0   1/3   0;...
               0   0   1/3];
    case 4
        % plane-strain
        uiuj=[ .4333   0   0;...
               0   .2337   0;...
               0   0   1/3];  
    otherwise
        uiuj=[ 0.9309  -0.4127   0.0799;...
              -0.4127   0.6827   0.0032;...
               0.0799   0.0032   0.3864];
end

k = 0.5*trace(uiuj);

[evecs,evals] = eigs(uiuj);

% rotate eigenvectors
rotz = [cosd(theta) -sind(theta) 0;...
        sind(theta) cosd(theta) 0;...
        0 0 1];
evecs = rotz*evecs;
uiuj = evecs*evals*evecs';

% plot inertial coordinates
scale = 1.1;

P1c = [1 0 0];
P2c = [0 1 0];
P3c = [0 0 1];
x1c = scale*[P1c(1) 0 -P1c(1)];
y1c = scale*[P1c(2) 0 -P1c(2)];
z1c = scale*[P1c(3) 0 -P1c(3)];

x2c = scale*[P2c(1) 0 -P2c(1)];
y2c = scale*[P2c(2) 0 -P2c(2)];
z2c = scale*[P2c(3) 0 -P2c(3)];

x3c = scale*[P3c(1) 0 -P3c(1)];
y3c = scale*[P3c(2) 0 -P3c(2)];
z3c = scale*[P3c(3) 0 -P3c(3)];

% plot principal axes
    
P1 = evecs(1:3,1)';
x1 = scale*[P1(1) 0 -P1(1)];
y1 = scale*[P1(2) 0 -P1(2)];
z1 = scale*[P1(3) 0 -P1(3)];
P2 = evecs(1:3,2)';
x2 = scale*[P2(1) 0 -P2(1)];
y2 = scale*[P2(2) 0 -P2(2)];
z2 = scale*[P2(3) 0 -P2(3)];
P3 = evecs(1:3,3)';
x3 = scale*[P3(1) 0 -P3(1)];
y3 = scale*[P3(2) 0 -P3(2)];
z3 = scale*[P3(3) 0 -P3(3)];





%% Plot 2D Glyph Slice
scatSize = 20;

% generate resolution
dim = 200;

cb = 'no';

% setup figures
scrsz = get(0,'ScreenSize');
screenH = scrsz(4);
screenW = scrsz(3);

% Z-axis slice
angle = linspace(0,2*pi,dim);
X = cos(angle);
Y = sin(angle);
Z = zeros(1,length(angle));

[X1,Y1,Z1,C1] = gen_glyph_magnitude(uiuj,X,Y,Z);

figure('Position',[100 screenH/2 600 600],'paperpositionmode','auto', 'color','white','InvertHardcopy','off');

hold on;
set(gca,'FontName','Times');
title('XY Plane');
surf([X1; X1],[Y1; Y1],[Z1; Z1],[C1; C1],'EdgeColor','flat','LineWidth',2.5);
% plot3(x1,y1,z1,'b',x2,y2,z2,'g',x3,y3,z3,'r','LineWidth',1.5);
plot3(x1c,y1c,z1c,'--k',x2c,y2c,z2c,'--k',x3c,y3c,z3c,'--k','LineWidth',1.5);
xlabel('X');ylabel('Y');zlabel('Z');
colormap('jet');caxis([0 1]);
if strcmp(cb,'yes')
    colorbar('FontName','Times');
end
axis_length=1;
axis equal;axis([-axis_length axis_length -axis_length axis_length -axis_length axis_length]);
grid on;box on;
view([0 0 1])
light('Position',[3 2 0],'Style','local');

hold off;


% Y-axis slice
X = cos(angle);
Z = sin(angle);
Y = zeros(1,length(angle));

[X2,Y2,Z2,C2] = gen_glyph_magnitude(uiuj,X,Y,Z);

figure('Position',[700 screenH/2 600 600],'paperpositionmode','auto', 'color','white','InvertHardcopy','off');

hold on;
set(gca,'FontName','Times');
title('XZ Plane');
surf([X2; X2],[Y2; Y2],[Z2; Z2],[C2; C2],'EdgeColor','flat','LineWidth',2.5);
% plot3(x1,y1,z1,'b',x2,y2,z2,'g',x3,y3,z3,'r','LineWidth',1.5);
plot3(x1c,y1c,z1c,'--k',x2c,y2c,z2c,'--k',x3c,y3c,z3c,'--k','LineWidth',1.5);
xlabel('X');ylabel('Y');zlabel('Z');
colormap('jet');caxis([0 1]);
if strcmp(cb,'yes')
    colorbar('FontName','Times');
end
axis_length=1;
axis equal;axis([-axis_length axis_length -axis_length axis_length -axis_length axis_length]);
grid on;box on;
view([0 1 0])
light('Position',[3 2 0],'Style','local');
hold off;

% X-axis slice
Y = cos(angle);
Z = sin(angle);
X = zeros(1,length(angle));

[X3,Y3,Z3,C3] = gen_glyph_magnitude(uiuj,X,Y,Z);

figure('Position',[1300 screenH/2 600 600],'paperpositionmode','auto', 'color','white','InvertHardcopy','off');

hold on;
set(gca,'FontName','Times');
title('YZ Plane');
surf([X3; X3],[Y3; Y3],[Z3; Z3],[C3; C3],'EdgeColor','flat','LineWidth',2.5);
% plot3(x1,y1,z1,'b',x2,y2,z2,'g',x3,y3,z3,'r','LineWidth',1.5);
plot3(x1c,y1c,z1c,'--k',x2c,y2c,z2c,'--k',x3c,y3c,z3c,'--k','LineWidth',1.5);
xlabel('X');ylabel('Y');zlabel('Z');
colormap('jet');caxis([0 1]);
if strcmp(cb,'yes')
    colorbar('FontName','Times');
end
axis_length=1;
axis equal;axis([-axis_length axis_length -axis_length axis_length -axis_length axis_length]);
grid on;box on;
view([1 0 0])
light('Position',[3 2 0],'Style','local');
hold off;

%% Export
if (0)
    image_file = sprintf('/home/memory/Documents/Latex-Docs/Thesis/Chapter41/images/chan_glyph_XY.eps');
    print('-f1','-depsc2','-noui','-r200',image_file);
    image_file = sprintf('/home/memory/Documents/Latex-Docs/Thesis/Chapter41/images/chan_glyph_XZ.eps');
    print('-f2','-depsc2','-noui','-r200',image_file);
    image_file = sprintf('/home/memory/Documents/Latex-Docs/Thesis/Chapter41/images/chan_glyph_YZ.eps');
    print('-f3','-depsc2','-noui','-r200',image_file);
end

%% Plot 3D Glyph

% build unit sphere
[X,Y,Z]=sphere(200);
[X4,Y4,Z4,C4] = gen_glyph_magnitude(uiuj,X,Y,Z);

% setup figures
scrsz = get(0,'ScreenSize');
screenH = scrsz(4);
screenW = scrsz(3);

figure('Position',[200 screenH/4 700 600],'paperpositionmode','auto', 'color','white','InvertHardcopy','off');

hold on;
set(gca,'FontName','Times');
glyphsurf=surf(X4,Y4,Z4,C4,'EdgeColor','k','EdgeAlpha',0,'FaceColor','interp','FaceAlpha',0.5,'FaceLighting','phong');
colormap('jet');caxis([0 1]);colorbar('FontName','Times');
axis_length=1;
axis equal;axis([-axis_length axis_length -axis_length axis_length -axis_length axis_length]);

plot3(x1,y1,z1,'b',x2,y2,z2,'g',x3,y3,z3,'r','LineWidth',1.5);
plot3(x1c,y1c,z1c,'--k',x2c,y2c,z2c,'--k',x3c,y3c,z3c,'--k','LineWidth',1.5);
grid on;box on;
view(-15,55);
xlabel('X');ylabel('Y');zlabel('Z','Rotation',0);

light('Position',[3 2 0],'Style','local');

hold off

%% Export
if (0)
    image_file = sprintf('/home/memory/Documents/Latex-Docs/Thesis/Chapter41/images/chan_glyph_3D.eps');
    print('-f4','-depsc2','-noui','-r200',image_file);
end

