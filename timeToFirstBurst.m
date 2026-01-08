clear all; close all;

pth1 = './Simulations/EKSeq'; clrId = [1 6 3]; 
%pth1 = './simulations/EKSeq_FinalGInitH'; clrId = [8 8 8]; 
%pth1 = './simulations/EKSeq_InitGFinalHalf'; clrId = [12 12 12];

folderNames = strcat(pth1,'/*p1*');
checkFiles = dir(folderNames);

simCount = 0;

for i = 1:length(checkFiles)
	i
	skipTrigger = 0;
	if checkFiles(i).isdir == 1
		dmy = dir(strcat(pth1,'/',checkFiles(i).name,'/*.mat'));
		dmy(1).name
		simCount = simCount + 1;

		% perturbation on and off times
		pertOnAndOff = [30.0 60.0 90.0 120.0 150.0 180.0];
		load(strcat(dmy(end).folder,'/',dmy(end).name), 'tme')
		delta = (tme(end)-tme(1))/1000;
		% which block to load
		wndowPlots(1,:) = ceil((pertOnAndOff.*60)./delta);
		% how far into the block to proceed
		wndowPlots(2,:) = mod(pertOnAndOff.*60,delta)./delta;

        for jj = 1:3
			key = 2*jj
			loadBlocks = [wndowPlots(1,key-1):wndowPlots(1,key)];
			tCourse = []; VCourse = []; ELCourse = [];
			for jjj = 1:length(loadBlocks)
				ld = load(strcat(dmy(loadBlocks(jjj)).folder,'/',dmy(loadBlocks(jjj)).name));
				tCourse = [tCourse; ld.tme];
				VCourse = [VCourse, ld.V];
				ELCourse = [ELCourse; ld.ELvals];
			end

			%% Kill off all indxes outside the perturbation.
			perturbIdx = find(ELCourse > -50);
			VCourse = VCourse(perturbIdx);
			tCourse = tCourse(perturbIdx);
			% Shift
			tCourse = tCourse - tCourse(1);
            
           firstBurst = double(py.timeToFirstBurst.computeTimeToFirstBurst(VCourse,tCourse))/1000;

           %{
           %plot to check the line is in the right place
           plot(tCourse/1000,VCourse); hold on;
           xline(firstBurst,'linewidth',2,'color','r')
           pause; close all
           %}
           
		   timeStore(simCount,jj) = firstBurst;
		   tCourse = []; VCourse = []; ELCourse = []; 
		end
	end

end


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


groupA1 = timeStore(:,1);
groupA2 = timeStore(:,2);
groupA3 = timeStore(:,3);

meanA1 = median(timeStore(:,1));
meanA2 = median(timeStore(:,2));
meanA3 = median(timeStore(:,3));

% Data
dmy1 = [1*ones(size(groupA1)),2*ones(size(groupA1)),3*ones(size(groupA1))];
dmy2 = [groupA1,groupA2, groupA3];
for kk = 1:size(dmy1,1)
	plot(dmy1(kk,:),dmy2(kk,:),'linewidth',1.5,'color',[.75 .75 .75]); hold on;
end

scatter(1*ones(size(groupA1)), groupA1, 250,'LineWidth',1.3,'MarkerEdgeColor',clrs(clrId(1),:),'MarkerFaceColor', clrs(clrId(1),:)); 
scatter(2*ones(size(groupA2)), groupA2, 250,'LineWidth',1.3,'MarkerEdgeColor',clrs(clrId(2),:),'MarkerFaceColor', clrs(clrId(2),:)); 
scatter(3*ones(size(groupA3)), groupA3, 250,'LineWidth',1.3,'MarkerEdgeColor',clrs(clrId(3),:),'MarkerFaceColor', clrs(clrId(3),:)); 

% Plot median lines
plot([0.8 1.2], [meanA1 meanA1], 'color', clrs(clrId(1),:), 'LineWidth', 2);
plot([1.8 2.2], [meanA2 meanA2], 'color', clrs(clrId(2),:), 'LineWidth', 2);
plot([2.8 3.2], [meanA3 meanA3], 'color', clrs(clrId(3),:), 'LineWidth', 2);

% Enhance the plot
xticks([1 2 3]);
xticklabels({'First', 'Second', 'Third'});
xtickangle(45); 
xlim([.3 3.6]);
ylabel('Time (sec)');
xlabel('Perturbation')
%ylim([0 35]);
title('Time to First Burst')
set(gca,'fontsize',20,'Box','off','linewidth',2.5);
hold off; % Release the figure for other plots
