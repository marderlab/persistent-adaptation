clear all; close all; clc;
burstIdx = [];

pth1 = './Simulations/EKSeq';
folderNames = strcat(pth1,'/*p1*');
checkFiles = dir(folderNames);
checkDirs = checkFiles([checkFiles.isdir]);
checkDirs = checkDirs(~ismember({checkDirs.name},{'.','..'}));
checkFiles = checkDirs;

i = 11; %Simulation `IKHP4T`

dmy = dir(strcat(pth1,'/',checkFiles(i).name,'/*.mat'));
windowMinutes = [30.0 60.0 90.0 120.0 150.0 180.0 205.0];
load(strcat(dmy(end).folder,'/',dmy(end).name), 'tme')
delta = (tme(end)-tme(1))/1000;
% which block to load
wndowPlots(1,:) = ceil((windowMinutes.*60)./delta);
% how far into the block to proceed
wndowPlots(2,:) = mod(windowMinutes.*60,delta)./delta;

out = strcat(dmy(end).folder,'/',compose('IKHP4T_%03d.mat', wndowPlots(1,:)));

for j = 1:7
	ld= load(out{j});
	idxMark = floor(length(ld.tme)*wndowPlots(2,j));
	condPlot(:,j) = ld.conds(5:7,idxMark); %% each column is a point in conductance space look at: Kd, KCa, A -- 5:7
	halfAxPlot(:,j) = ld.halfAcs(8:10,idxMark); %% each column is a point in halfAx space -- 8:10
	clear ld
end

%%%%%%%%%%%%%%%

x = condPlot(1,:);
y = condPlot(2,:);
z = condPlot(3,:);
N = numel(x);

figure; hold on;

% Points with hollow circles so labels stand out
%plot3(x, y, z, 'o', ...
%      'MarkerSize', 10, ...
%      'MarkerFaceColor', 'none', ...
%      'LineWidth', 1.5);

% ---- CIRCLED NUMBER LABELS ----
for k = 1:N
    if k <= 20
        % Unicode circled digits: ① is decimal 9312
        lbl = char(9311 + k);   % 1→①, 2→②, ...
    else
        % Fallback if you ever go beyond 20
        lbl = sprintf('%d', k);
    end
    
    % Slight vertical offset so label sits above the marker
    text(x(k), y(k), z(k), lbl, ...
        'FontSize', 18, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom');   % put label above the point
end

% ---- ARROWS from point k → k+1 ----
dx = diff(x);
dy = diff(y);
dz = diff(z);

quiver3(x(1:end-1), y(1:end-1), z(1:end-1), ...
        dx, dy, dz, 0, 'LineWidth', 1.5, 'MaxHeadSize', 0.1, 'Color', [0.6 0.6 0.6]);   % grey arrows

xlabel('$\bar{g}_{Kd} (\mu S)$','Interpreter', 'latex','FontSize',25);
ylabel('$\bar{g}_{KCa} (\mu S)$','Interpreter', 'latex','FontSize',25);
zlabel('$\bar{g}_{A} (\mu S)$','Interpreter', 'latex','FontSize',25);
set(gca,'FontSize',25);
title('Conductance Space','fontsize',20)
view(3);
xlim([6 10]);
ylim([10 25]);
zlim([8 24]);
grid on;

%%%%% Draw sphere

sx = condPlot(1,2:end);
sy = condPlot(2,2:end);
sz = condPlot(3,2:end);
pts = [sx(:) sy(:) sz(:)];

% ----- 1. Center of sphere (mean of points) -----
ctr = mean(pts,1);   % [cx cy cz]

% ----- 2. Radius: max distance from center -----
dists = sqrt(sum((pts - ctr).^2, 2));
R = max(dists);

% ----- 3. Generate sphere mesh -----
[phi,theta] = meshgrid(linspace(0,2*pi,50), linspace(0,pi,30));
Xs = ctr(1) + R .* sin(theta).*cos(phi);
Ys = ctr(2) + R .* sin(theta).*sin(phi);
Zs = ctr(3) + R .* cos(theta);

% Plot on top of your existing figure
hold on;
hs = surf(Xs, Ys, Zs);

set(hs, 'FaceAlpha', 0.08, ...   % very transparent
        'EdgeAlpha', 0.15, ...
        'FaceColor', [1, 177, 80]./256, ... 
        'EdgeColor', [1, 177, 80]./256);  

xlim([4 15])

%% orange line
p1 = [x(1), y(1), z(1)];   % point 1
v = p1 - ctr;              % direction from sphere center
dist = norm(v);

% closest point on sphere
u = v / dist;                  % unit direction
spherePoint = ctr + R * u;     % projection onto sphere surface

% vector from p1 to the sphere surface
vec = spherePoint - p1;

% draw arrow
quiver3(p1(1), p1(2), p1(3), ...
        vec(1), vec(2), vec(3), 0, ...
        'LineWidth', 2, ...
        'MaxHeadSize', 0.1, ...
        'Color', [244, 127, 32]./256);   

axis equal

%%%%%%%%%%%%


figure;

% Extract your 4 points
X = -12.3 - halfAxPlot(1,end-3:end);
Y = -28.3 - halfAxPlot(2,end-3:end);
Z = -27.2 - halfAxPlot(3,end-3:end);

% Plot hollow markers
plot3(X, Y, Z, 'o', 'MarkerSize', .1, ...
      'LineWidth', .1, 'MarkerFaceColor','none');
hold on;

circ = [];
% ---- LABELS ④⑤⑥⑦ ----
for k = 1:4
    circ = char(9311 + (k+3));   % 4→④, 5→⑤, 6→⑥, 7→⑦
    text(X(k), Y(k), Z(k), circ, ...
        'FontSize', 22, ...
        'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle');
end

% ---- QUIVERS BETWEEN THEM ----
dX = diff(X);
dY = diff(Y);
dZ = diff(Z);

quiver3(X(1:end-1), Y(1:end-1), Z(1:end-1), ...
        dX, dY, dZ, 0, ...
        'LineWidth', 1.8, ...
        'MaxHeadSize', 0.1, ...
        'Color', [0.5 0.5 0.5]);   % gray arrows

% ---- AXES & LABELING ----
grid on;
xlabel('Kd_M (mV)','FontSize',20);
ylabel('KCa_M (mV)','FontSize',20);
zlabel('A_M (mV)','FontSize',20);

title('Half-(In)Acitvation Space','FontSize',20);
set(gca,'FontSize',20);

view(3);


%ylim([-25 -20]);
%xlim([-8 -5.5]);
%zlim([-25 -19]);
