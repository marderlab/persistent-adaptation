clear all; close all; clc;
burstIdx = [];

pth1 = './simulations/EKSeq';
%pth1 = './simulations/EKSeq_swap';
%pth1 = './simulations/EKSeq_swap_hInf';
%pth1 = './simulations/EKSeq_FinalGInitH';
%pth1 = './simulations/EKSeq_InitGFinalHalf';

% set to 0 if youre plotting one and want to examine
% set to 1 if you're plotting a bunch
closeFigAfterRun = 0;
folderNames = strcat(pth1,'/*p1*');
checkFiles = dir(folderNames);
svdParams = [];

for i = 1:length(checkFiles)
	skipTrigger = 0;
	if checkFiles(i).isdir == 1
		for ii = 1:length(checkFiles)
			% just create plots
			if strcmp(checkFiles(ii).name,strcat(checkFiles(i).name,'.png'))
				skipTrigger = 1;
			end
		end
		tic
		if skipTrigger == 1
			sprintf(strcat(checkFiles(i).name,' appears to be plotted.'));
			continue
		end
		dmy = dir(strcat(pth1,'/',checkFiles(i).name,'/*.mat'));
		if sum(dmy(1).name(1:9) == 'xxxNEWxxx') == 9
		%	sprintf(strcat(checkFiles(i).name,' appears to be incomplete. skipping plot...'))
		%	continue
		else
			totSecs = length(dmy);
			'setting up'
			% compute perturbation
			shades.totEK = []; shades.totI = []; shades.totTme = [];
			for j = length(dmy):-1:1
				load(strcat(dmy(j).folder,'/',dmy(j).name),'EKvals','Ivals','tme')
				shades.totEK = [EKvals; shades.totEK];
				shades.totI = [Ivals; shades.totI];
				shades.totTme = [tme; shades.totTme];
				clear EKvals Ivals tme
			end

			% compute windows
			%20, 60+20, 180+20, 240+20,  20+360, 20+420
			%windowMinutes = [20.0 80.0 200.0 260.0 380.0 440.0];
			windowMinutes = [30.0 60.0 90.0 120.0 150.0 180.0 205.0];
			load(strcat(dmy(end).folder,'/',dmy(end).name), 'tme')
			delta = (tme(end)-tme(1))/1000;
			% which block to load
			wndowPlots(1,:) = ceil((windowMinutes.*60)./delta);
			% how far into the block to proceed
			wndowPlots(2,:) = mod(windowMinutes.*60,delta)./delta;

			for j = length(dmy):-1:1

			dmy(j).name
    		params = plotFigCompute(dmy(j),j,totSecs,wndowPlots,0,shades,svdParams,closeFigAfterRun);

			end
		end
		toc
	end
end
