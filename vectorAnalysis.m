clear all; close all; clc

% Color Library
clrs(1,:)  = [31 120 180]/256; %dark blue
clrs(2,:)  = [115, 156, 112]/256; %light green
clrs(3,:)  = [11, 163, 0]/256; %dark green
clrs(4,:)  = [166 206 227]/256; %light blue
clrs(5,:)  = [173, 101, 101]/256; 
clrs(6,:)  = [179 50 50]/256; %red
clrs(7,:)  = [256 128 0]/256; %orange
clrs(8,:)  = [75,0,130]/256; %indigo
clrs(9,:)  = [243, 204, 255]/256; %light purp
clrs(10,:) = [255,20,147]/256; %deep pink
clrs(11,:) = [255,182,193]/256; %light pink
clrs(12,:) = [87, 179, 176]/256; %dark cyan
clrs(13,:) = [72,209,204]/256; %light cyan
clrs(14,:) = [18, 70, 105]/256; %v dark blue
clrs(15,:) = [36, 102, 100]/256; %v dark cyan
clrs(16,:) = [171, 107, 44]/256; %v dark orange
clrs(17,:) = [102, 14, 14]/256; %v dark red
clrs(18,:) = [24, 97, 20]/256; %v dark green
clrs(19,:) = [250, 178, 105]/256; %light orange


pth1 = './Simulations/EKSeq'; clrId = [1 5 3 7];
%pth1 = './simulations/EKSeq_FinalGInitH'; clrId = [8 8 8 8];
%pth1 = './simulations/EKSeq_InitGFinalHalf'; clrId = [12 12 12 12];

folderNames = strcat(pth1,'/*p1*');
checkFiles = dir(folderNames);

paramValsCell = {};
alfaThresh = .05;
simCount = 0;

for i = 1:length(checkFiles)
	skipTrigger = 0;

	if checkFiles(i).isdir == 1
		dmy = dir(strcat(pth1,'/',checkFiles(i).name,'/*.mat'));
		simCount = simCount + 1;

		% perturbation on and off times
		pertOnAndOff = [30.0 60.0 90.0 120.0 150.0 180.0 205.0];
		load(strcat(dmy(end).folder,'/',dmy(end).name), 'tme')
		delta = (tme(end)-tme(1))/1000;
		% which block to load
		wndowPlots(1,:) = ceil((pertOnAndOff.*60)./delta);
		% how far into the block to proceed
		wndowPlots(2,:) = mod(pertOnAndOff.*60,delta)./delta;

		for j = 1:length(pertOnAndOff)

			lastTme = wndowPlots(2,j)*(tme(end)-tme(1))/1000 - 10;
			firstTme = lastTme - 30;
			firstProportion = firstTme*1000/(tme(end)-tme(1));

			lastIdx = floor(wndowPlots(2,j)*length(tme));
			firstIdx = floor(firstProportion*length(tme));

			fld = strcat(pth1,'/',checkFiles(i).name,'/',checkFiles(i).name(end-5:end),'_',sprintf('%03d', wndowPlots(1,j)),'.mat');
			ld = load(fld,'conds','halfAcs','alfa','tme');

			condsFull(simCount,j,:) = mean(ld.conds(:,firstIdx:lastIdx),2);
			halfAcsFull(simCount,j,:) = mean(ld.halfAcs(:,firstIdx:lastIdx),2);
			tmeFull(simCount,j,:) = ld.tme(firstIdx:lastIdx);
			alfaFull(simCount,j,:) = ld.alfa(firstIdx:lastIdx);

		end
	end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Max Conductances
% Step 1: Compute differences along 2nd dimension
diffs = condsFull(:,2:end,:) - condsFull(:,1:end-1,:); % size 20 x 6 x 7
vector_lengths = sqrt(sum(diffs.^2, 3));  % size 20 x 6

f = figure;
hold on;
for i = 1:20
	plot([1:6],vector_lengths(i,:),'linewidth',2,'color',[.75 .75 .75]); hold on;
end

% For each of the 6 categories (x = 1:6)
for i = 1:6
    % Extract 20 values for this category
    data = vector_lengths(:,i);
    
    % Plot as scatter: x-axis = step index, y-axis = length
    if mod(i,2) == 0
    	clx = clrs(clrId(1),:);
    else
    	clx = clrs(clrId(2),:);
    end

	scatter(i*ones(size(data)), data, 250, 'filled','MarkerFaceColor', clx);
	plot([i*ones(size(data))-.25 i*ones(size(data))+.25],[median(data) median(data)], 'linewidth', 2, 'color', clx);

end

ylabel('Vector Length (uS)');
ax = gca;
xlim([0 7])
ax.FontSize = 20;
ax.XTickLabelRotation = 45;  
ax.XTickLabel = {'','After Pert 1','After Wash 1','After Pert 2','After Wash 2','After Pert 3','After Wash 3',''};
f.Position = [584 481 1024 520];

hold off;

%%%%%%%%%%%%%%%%

% Initialize storage
sphere_radii = zeros(20,1);
sphere_distances = zeros(20,1);

for n = 1:20
    % Extract last 6 points
    subset = squeeze(condsFull(n,2:7,:));  % size 6 x 7
    
    % Compute approximate center (mean)
    center = mean(subset, 1);
    
    % Compute distances from center to all 6 points
    dists = sqrt(sum((subset - center).^2, 2));
    
    % Radius: maximal distance from center
    radius = max(dists);
    
    % Get first vector
    v1 = squeeze(condsFull(n,1,:));
    
    % Distance of v1 to center
    d_v1_center = norm(center' - v1);

    % Distance to sphere surface
    dist_outside = d_v1_center - radius;
    % Store results
    sphere_radii(n) = radius;
    sphere_distances(n) = dist_outside;  % >0: outside, <0: inside
end

% Suppose M is an n-by-k matrix (rows = subjects, cols = conditions)
conditionLabels = {'After_Pert1','After_Wash1','After_Pert2', ...
                   'After_Wash2','After_Pert3','After_Wash3'};
T               = array2table([vector_lengths], ...
                  'VariableNames', [conditionLabels]);

writetable(T,'vector_lengths_maxConds.csv');

%{
% significant pairs
                  Comparison        P_adj
1  After_Pert1 - After_Pert3 4.991188e-02
2  After_Pert1 - After_Wash1 4.515942e-12
5  After_Pert1 - After_Wash2 2.257238e-11
8  After_Pert1 - After_Wash3 2.990541e-10

3  After_Pert2 - After_Wash1 1.668883e-05
4  After_Pert3 - After_Wash1 5.433679e-05

6  After_Pert2 - After_Wash2 4.882951e-05

7  After_Pert3 - After_Wash2 1.508450e-04

9  After_Pert2 - After_Wash3 2.646397e-04

10 After_Pert3 - After_Wash3 7.490201e-04
%}

%%%%%%%%%%%

f2 = figure;
tmp = [sphere_distances,sphere_radii];
for i = 1:2
    % Extract 20 values for this category
    data = tmp(:,i);
    
    % Plot as scatter: x-axis = step index, y-axis = length
    if mod(i,2) == 0
    	scatter(i*ones(size(data)), data, 250, 'filled','MarkerFaceColor',clrs(clrId(3),:))
    else
    	scatter(i*ones(size(data)), data, 250, 'filled','MarkerFaceColor',clrs(clrId(4),:))
    end

    hold on;

end

xlim([0 3])

ylabel('Vector Length (uS)');
ax = gca;
xlim([0 3])
ax.FontSize = 20;
ax.XTick = [0:3];

ax.XTickLabelRotation = 45;  
ax.XTickLabel = {'','Min Distance From Sphere','Sphere Radius',''};
f2.Position = [584 224 660 777];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Half Acs

% Step 1: Compute differences along 2nd dimension
diffs = halfAcsFull(:,2:end,:) - halfAcsFull(:,1:end-1,:); % size 20 x 6 x 11

% Step 2: Compute vector lengths (norm along 3rd dimension)
vector_lengths = sqrt(sum(diffs.^2, 3));  % size 20 x 11

% Now plotting with swapped axes:
f = figure;
hold on;
for i = 1:20
	plot([1:6],vector_lengths(i,:),'linewidth',2,'color',[.75 .75 .75]); hold on;
end
% For each of the 6 categories (x = 1:6)
for i = 1:6
    % Extract 20 values for this category
    data = vector_lengths(:,i);
    
    % Plot as scatter: x-axis = step index, y-axis = length
    if mod(i,2) == 0
    	clx = clrs(clrId(1),:);
    else
    	clx = clrs(clrId(2),:);
    end

	scatter(i*ones(size(data)), data, 250, 'filled','MarkerFaceColor', clx);
	plot([i*ones(size(data))-.25 i*ones(size(data))+.25],[median(data) median(data)], 'linewidth', 2, 'color', clx);


end

ylabel('Vector Length (mV)');
ax = gca;
xlim([0 7])
ax.FontSize = 20;
ax.XTickLabelRotation = 45;  
ax.XTickLabel = {'','After Pert 1','After Wash 1','After Pert 2','After Wash 2','After Pert 3','After Wash 3',''};
f.Position = [584 481 1024 520];

hold off;

% Suppose M is an n-by-k matrix (rows = subjects, cols = conditions)
conditionLabels = {'After_Pert1','After_Wash1','After_Pert2', ...
                   'After_Wash2','After_Pert3','After_Wash3'};
T               = array2table([vector_lengths], ...
                  'VariableNames', [conditionLabels]);

writetable(T,'vector_lengths_HalfAcs.csv');   % one tidy CSV is easiest
%{
% significant pairs

                 Comparison        P_adj
1 After_Pert1 - After_Pert2 7.775879e-04
2 After_Pert1 - After_Pert3 5.566011e-06
6 After_Pert1 - After_Wash3 6.795411e-03

4 After_Pert3 - After_Wash1 8.263003e-07
3 After_Pert2 - After_Wash1 1.669622e-04
5 After_Wash1 - After_Wash2 4.111287e-02
7 After_Wash1 - After_Wash3 1.768057e-03

%}


for n = 1:20 
    % Example for one model n:
    half = squeeze(halfAcsFull(n,:,:));  % Size: [timepoints x half-(in)ac values]

    % Wash–Pert–Wash segment
    v1 = half(6,:) - half(5,:);  % From W2 to P3
    v2 = half(7,:) - half(6,:);  % From P3 to W3
    angle1 = acosd(dot(v1,v2)/(norm(v1)*norm(v2)));

    % Pert–Wash–Pert segment
    v3 = half(4,:) - half(3,:);  % From P2 to W2
    v4 = half(5,:) - half(4,:);  % From W2 to P3
    angle2 = acosd(dot(v3,v4)/(norm(v3)*norm(v4)));

    angle_store(n,:) = [angle1, angle2];
end


f2 = figure;
for i = 1:2
    % Extract 20 values for this category
    data = angle_store(:,i);
    
    % Plot as scatter: x-axis = step index, y-axis = length
    if mod(i,2) == 0
    	scatter(i*ones(size(data)), data, 250, 'filled','MarkerFaceColor',clrs(clrId(3),:))
    else
    	scatter(i*ones(size(data)), data, 250, 'filled','MarkerFaceColor',clrs(clrId(4),:))
    end

    hold on;

end

xlim([0 3])

ylabel('Angle (deg)');
ax = gca;
xlim([0 3])
ylim([90 180])
ax.FontSize = 20;
ax.XTick = [.80 1 1.80 2];
ax.XTickLabel = {'Angle between', 'After Wash 2 & After Pert 3','Angle between','After Pert 3 & After Wash 3'};
ax.XTickLabelRotation = 45;  
f2.Position = [584 224 660 777];

