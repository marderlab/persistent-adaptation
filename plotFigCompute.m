 function params = plotFigCompute(filePointer,currentSec,totSecs,wndwSecs,onlyStat,shades,svdParams,closeFigAfterRun)

	%fontsize
	fs = 20;
	windw = 2e3; windwPct = .95;
	resl = 1;
	timeUnit.conv = 60; timeUnit.nme = ' [mins]';
	%timeUnit.conv = 3600; timeUnit.nme = ' [hrs]';

	ld = load(strcat(filePointer.folder,'/',filePointer.name));

	chunkSize = (ld.tme(end)-ld.tme(1))/1000;

	if onlyStat ~= 1

		clrs(1,:)  = [31 120 180]/256; %[0 0.4470 0.7410]; dark blue
		clrs(2,:)  = [115, 156, 112]/256; %[0.8500 0.3250 0.0980]; light green
		clrs(3,:)  = [11, 163, 0]/256; %[0.9290 0.6940 0.1250]; dark green
		clrs(4,:)  = [166 206 227]/256; %[0.4940 0.1840 0.5560]; light blue
		clrs(5,:)  = [173, 101, 101]/256; %[240,128,128]/256; %new light red to not conflict with pink %[251 154 153]/256; %[0.4660 0.6740 0.1880]; light red
		clrs(6,:)  = [179 50 50]/256; %[227 26 28]/256; %[0.3010 0.7450 0.9330]; red
		clrs(7,:)  = [256 128 0]/256; %[253 141 1.11]/256; %[0.6350 0.0780 0.1840]; organge

		clrs(8,:)  = [75,0,130]/256; %indigo
		clrs(9,:)  = [243, 204, 255]/256; %[147,112,219]/256; %light purp
		clrs(10,:) = [255,20,147]/256; %deep pink
		clrs(11,:) = [255,182,193]/256; %light pink
		clrs(12,:) = [87, 179, 176]/256; %dark cyan
		clrs(13,:) = [72,209,204]/256; %light cyan

		clrs(14,:) = [18, 70, 105]/256; %v dark blue
		clrs(15,:) = [36, 102, 100]/256; %v dark cyan
		clrs(16,:) = [171, 107, 44]/256; %v dark orange
		clrs(17,:) = [102, 14, 14]/256; %v dark red
		clrs(18,:) = [24, 97, 20]/256; %v dark green

		clrs(19,: ) = [250, 178, 105]/256; %light orange


		clrOrder = [1,3,13,8,6,10,7];
		clrOrderShifts = [1,4,3,2,13,12,8,6,10,7,19];

		% design shadings and plots
		if totSecs == currentSec
	
			fig = figure();
			set(fig,'Position',[0 0 1400 2000]);

			ax1 = subplot(5,1,2);
			%ax1.Position(2) = ax1.Position(2) -.04;
			%xlabel(ax1,strcat('Time ',timeUnit.nme),'fontweight','bold');
			yLabStuff = ylabel(ax1,'V [mV]','fontsize',fs);
			set(ax1,'fontsize',fs);
			ylim(ax1,[-70 85]);
			xLimSet = [0 ld.tme(end)/(timeUnit.conv*1000)]; 
			xlim(ax1,xLimSet);
			shader(1,shades);

			ax3_1 = subplot(10,1,5);
			shader(1,shades);

			ax3_2 = subplot(10,1,6);
			shader(1,shades);

			ax4 = subplot(5,1,4);
			shader(1,shades);

			ax5 = subplot(5,1,5);
			shader(1,shades);

		end

				
		% inset traces
		redLnStrtIx = [];
		totTimeOfBlocks = (ld.tme(end) - ld.tme(1))/(timeUnit.conv*1000);

		trgPrm = 0;
		for jj = 1:7
			if currentSec == wndwSecs(1,jj)
				lblOnXY = 0;
				if jj == 1
					posMat = [0.02 0.875 .2 .1];
					trgPrm = 1;
					colorMat = [1 6];
				elseif jj == 2
					posMat = [0.27 0.875 .2 .1];
					trgPrm = 1;
					colorMat = [1 6];
					lblOnXY = 1;
				elseif jj == 3
					posMat = [0.52 0.875 .2 .1];
					colorMat = [1 6];
				elseif jj == 4
					posMat = [0.77 0.875 .2 .1];
					colorMat = [1 6];
				elseif jj == 5
					posMat = [0.145 0.765 .2 .1];
					colorMat = [1 6];
				elseif jj == 6
					posMat = [0.395 0.765 .2 .1];
					colorMat = [1 6];
				elseif jj == 7
					posMat = [0.645 0.765 .2 .1];
					colorMat = [1 6];
					trgTmeScl = 1;
				%{
				elseif jj == 8
					posMat = [0.36 0.6280 .3 .12];
					trgTmeScl = 1;
				elseif jj == 9
					posMat = [0.69 0.6280 .3 .12];
					trgTmeScl = 1;
				%}
				end
				ax{jj} = subplot('Position',posMat);
				lxn = totTimeOfBlocks*(wndwSecs(1,jj) - 1) + wndwSecs(2,jj)*totTimeOfBlocks;
				tracePlotSection(ax{jj},jj,lblOnXY,colorMat,ld,wndwSecs,totTimeOfBlocks);
			end
		end
		% rest of the plots.
		ax1 = subplot(5,1,2);
		set(ax1,'LineWidth',2.5)
		hold on;
		q1 = plot(ld.tme(1:resl:end)/(timeUnit.conv*1000),ld.V(1:resl:end),'LineWidth',2,'Color',[0 0.4470 0.7410],'LineStyle','-');
		set(gca, 'XTickLabel', {});

		% Conductance Plot
		ax3_2 = subplot(10,1,6); %bottom graph
		set(gca, 'XTickLabel', {});
		for jj = 1:7
			h2 = plot(ax3_2,ld.tme(1:resl:end)/(timeUnit.conv*1000),ld.conds(jj,1:resl:end),'linewidth',2.5,'color',clrs(clrOrder(jj),:));
			hold on;
		end
		box off;
		ax3_1 = subplot(10,1,5); %top graph
		for jj = 1:7
			h1 = plot(ax3_1,ld.tme(1:resl:end)/(timeUnit.conv*1000),ld.conds(jj,1:resl:end),'linewidth',2.5,'color',clrs(clrOrder(jj),:));
			hold on;
		end
		box off;
		ax3_1.XAxis.Visible = 'off'; set(ax3_1, 'XTick', []);

		if totSecs == currentSec
			xlim(ax3_1,[0 ld.tme(end)/(timeUnit.conv*1000)]); 
			xlim(ax3_2,[0 ld.tme(end)/(timeUnit.conv*1000)]);
		end

		if currentSec == 1
			% obtain max and min values for each line bleh.
			maxVal = []; minVal = [];
			lines = findobj(ax3_1, 'Type', 'Line');
			for rr = 1:7
				subsetLines = lines(rr:7:size(lines));
				ttDmy = [subsetLines(:).YData];
				minVal(rr) = min(ttDmy);
				% find maxVal but cap at 100
				dddmy = max(ttDmy);
				if dddmy < 100
					maxVal(rr) = max(ttDmy);
				else
					%maxVal(rr) = 100;
					maxVal(rr) = 200;
				end
			end
			maxVal = fliplr(maxVal); minVal = fliplr(minVal);
			clear ttDmy
			% Assuming 'ymax' is the highest value on the y-axis
			ylabmax = floor(max(maxVal) / 5) * 5;  % Round up to the nearest multiple of 5
			midpt = round(ylabmax / 2 / 5) * 5;  % Calculate a middle point that's also a multiple of 5
			% Set the y-axis ticks
			if midpt == 5
				midpt = 10;
				ylabmax = 15;
			end
			yticks([5 midpt ylabmax]);
			%yticks([5 15 25]);
			%ylim(ax3_1,[1.2 max(maxVal)*1.1]); set(ax3_1, 'XColor', 'none','linewidth', 2.5,'fontsize',fs);
			%ylim(ax3_2,[0 1.2]); set(ax3_2,'linewidth', 2.5,'fontsize',fs);
			ylim(ax3_1,[2.5 max(maxVal)*1.1]); set(ax3_1, 'XColor', 'none','linewidth', 2.5,'fontsize',fs);
			ylim(ax3_2,[0 2.5]); set(ax3_2,'linewidth', 2.5,'fontsize',fs);

			ylab = ylabel(ax3_2,'g̅_i (μS)','fontsize',fs);
			ylab.Position(1) = ax1.YLabel.Position(1);
			ylab.Position(2) = 2;
			%ax1.YLabel.Position(1)/diff(ax1.XLim)*ax1.Position(3)+ax1.Position(1)
			%dmy5 = get(ylab,'Position'); set(ylab,'Position',[dmy5(1) dmy5(2) + .05])

			%ax3_2 - bottom, ax3_1 - top
			dmy2 = get(ax3_2, 'Position');
			botLeftCornerOfTopGraph = dmy2(2) + dmy2(4) + .01;
			dmy3 = get(ax3_1, 'Position');
			topLeftCornerOfTopGraph = dmy3(2) + dmy3(4);
			newHeight = topLeftCornerOfTopGraph - botLeftCornerOfTopGraph;
			set(ax3_1,'Position',[dmy3(1) botLeftCornerOfTopGraph dmy3(3) newHeight])
			dmy3 = get(ax3_1, 'Position'); dmy2 = get(ax3_2, 'Position'); 
			dmy2_b = get(ax3_2, 'Position');
			dmy3_b = get(ax3_1, 'Position');

			%annotation('line', [dmy2_b(1)-0.005, dmy2_b(1)+0.005], [dmy2_b(2)+dmy2_b(4)-0.005-.07 dmy2_b(2)+dmy2_b(4)+0.005-.07], 'Color', 'k','linewidth',2);
			%annotation('line', [dmy2_b(1)-0.005, dmy2_b(1)+0.005], [dmy2_b(2)+dmy2_b(4)-0.005+.01-.07 dmy2_b(2)+dmy2_b(4)+0.005+.01-.07], 'Color', 'k','linewidth',2);

			dmyP2 = get(ax3_2,'Position');
			z = legend('Na','CaT','CaS','H','Kd','KCa','A');

			z.Position(1) = dmyP2(1) + dmyP2(3) + .01;
			yAvg = (z.Position(2) + dmyP2(2))/2;
			z.Position(2) = yAvg - .05;

		end

		% Half Acitvation Plots
		ax4 = subplot(5,1,4);
		set(ax4,'LineWidth',2.5);
		hold on;
		lnType{1} = '-'; lnType{2} = '--';
		ix = [1 1;1 1; 3 1;3 1;13 2;13 2; 8 1; 6 2; 10 1; 7 2; 7 1];

		for jj = [1 3 5 7 8 9 10] % ix for half acts
			plot(ax4,ld.tme(1:resl:end)/(timeUnit.conv*1000),ld.halfAcs(jj,1:resl:end),'linewidth',2.5,'linestyle',lnType{ix(jj,2)},'color',clrs(ix(jj,1),:));
			hold on;
		end

		if currentSec == 1
			ylim([-25 25]);
			ylabel('ΔV_{1/2,m} (mV)');
			%title('Half-Activation Shifts','fontsize',fs);
			set(gca,'fontsize',fs,'XTickLabel', {});
		end

		if totSecs == currentSec
			xlim(ax4,[0 ld.tme(end)/(timeUnit.conv*1000)]);
		end

		% Half Inacitvation Plots
		ax5 = subplot(5,1,5);
		set(ax5,'LineWidth',2.5);
		for jj = [2 4 6 11] % ix for half inax
			plot(ax5,ld.tme(1:resl:end)/(timeUnit.conv*1000),ld.halfAcs(jj,1:resl:end),'linewidth',2.5,'linestyle',lnType{ix(jj,2)},'color',clrs(ix(jj,1),:));
			hold on;
		end

		if currentSec == 1
			ylim([-25 25]);
			ylabel('ΔV_{1/2,h} (mV)');
			xlabel('Time (min)')
			%title('Half-Inactivation Shifts','fontsize',fs);
			set(gca,'fontsize',fs);
		end

		if totSecs == currentSec
			xlim(ax5,[0 ld.tme(end)/(timeUnit.conv*1000)]);
		end

		params = [];
		if totSecs == currentSec || trgPrm == 1
			trgPrm = 0;
			params = [mean(ld.conds(:,end-2000:end),2); mean(ld.halfAcs(:,end-2000:end),2)]; %% last section, 2 and 1
		end

		clearvars ld.V ld.alfa ld.conds ld.halfAcs ld.tme ld.EKvals

		if currentSec == 1
			%% TO BE RUN LAST!
			ax5.Position(2) = ax5.Position(2) - .04;
			ax4.Position(2) = ax4.Position(2) - .05;
			ax3_1.Position(2) = ax3_1.Position(2) -.05;
			ax3_2.Position(2) = ax3_2.Position(2) -.05;
			ax1.Position(2) = ax1.Position(2) - .065;
			ax1.Position(4) = ax1.Position(4) + .03;
			G = findall(gcf, 'Type', 'axes');
			for qq = 5:7
				G(qq).Position(2) = G(qq).Position(2) - .03;
			end

			%% Additional labelling on the trace plot
			axes(ax1);

			dmyP = get(ax1,'Position');
			dmyPy = ylim(ax1);
			botYLoc = dmyP(4)*(-50-dmyPy(1))/(diff(dmyPy))+dmyP(2);
			topYLoc = dmyP(4)*(50-dmyPy(1))/(diff(dmyPy))+dmyP(2);			
			botXLoc = dmyP(1)+dmyP(3)-.01;
			annotation('line', [dmyP(1) dmyP(1)],[(1-((dmyPy(2)-60)/diff(dmyPy)))*dmyP(4)+dmyP(2),dmyP(4)+dmyP(2)+.001],'Color','white','linewidth',2.52);

			for ii = 1:7
				redtme = totTimeOfBlocks*(wndwSecs(1,ii) - 1) + wndwSecs(2,ii)*totTimeOfBlocks;
				h = xline(redtme,'b',num2str(ii),'linewidth',4);
				%h.LabelVerticalAlignment = 'top';
				h.LabelHorizontalAlignment = 'center';
				h.FontSize = fs;
				h.LabelOrientation = 'horizontal';
				uistack(h, 'bottom');
			end
			saveas(gcf,strcat(filePointer(1).folder,'.png'),'png');
			if closeFigAfterRun == 1
		 		close all;
		 	end
	 	end

 	
 	end

	function lxn = tracePlotSection(axs,i,lblOnXY,firstSecondColor,ld,wndwSecs,totTimeOfBlocks)
		
		checkTime = ((ld.tme(end)-ld.tme(1))*wndwSecs(2,i))/1000;

		if checkTime < 3.0 %if you're less than 3 seconds, then you need to load previous 
			dmy11 = strcat(filePointer.folder,'/',filePointer.name);
			dmy12 = strcat(dmy11(1:end-7),sprintf('%03d',wndwSecs(1,i)-1),'.mat');
			ld2 = load(dmy12);

			%mdIdx = floor(wndwSecs(2,i)*length(ld.tme))+1;
			forwardOfMiddle = find(ld.tme <= ld.tme(1)+2000);
			dmyAt = [ld.tme(forwardOfMiddle)/1000];
			dmyAV = [ld.V(forwardOfMiddle)];
			dmyAEK = [ld.EKvals(forwardOfMiddle)];

			%endIdx = floor(length(ld2.tme)*(3 - (ld.tme(end)-ld.tme(stIdx))/1000)/totTimeOfBlocks);
			backwardOfMiddle = find(ld.tme >= ld.tme(end)-1000);
			dmyBt = [ld2.tme(backwardOfMiddle)/1000];
			dmyBV = [ld2.V(backwardOfMiddle)];
			dmyBEK = [ld2.EKvals(backwardOfMiddle)];

			totPlot_T = [dmyBt' dmyAt']; totPlot_V = [dmyBV dmyAV];

			plot(axs,totPlot_T,totPlot_V,'linewidth',1.5,'color',clrs(firstSecondColor(1),:)); hold on;
			perturbSelect = find([dmyBEK' dmyAEK']>-80.0);
			plot(axs,totPlot_T(perturbSelect),totPlot_V(perturbSelect),'linewidth',1.5,'color',clrs(firstSecondColor(2),:));
		else

			%dmy11 = strcat(filePointer.folder,'/',filePointer.name);

			mdIdx = floor(wndwSecs(2,i)*length(ld.tme))+1;
			dmyA = (ld.tme >= ld.tme(mdIdx)-1000); dmyB = (ld.tme <= ld.tme(mdIdx)+2000);
			midIdx = find( dmyA.*dmyB );
			dmyAt = [ld.tme(midIdx)/1000];
			dmyAV = [ld.V(midIdx)];
			dmyAEK = [ld.EKvals(midIdx)];

			totPlot_T = [dmyAt']; totPlot_V = [dmyAV];

			plot(axs,totPlot_T,totPlot_V,'linewidth',1.5,'color',clrs(firstSecondColor(1),:)); hold on;
			perturbSelect = find([dmyAEK]>-80.0);
			plot(axs,totPlot_T(perturbSelect),totPlot_V(perturbSelect),'linewidth',1.5,'color',clrs(firstSecondColor(2),:));
			xlim([totPlot_T(1)-.5,totPlot_T(end)+.5])
			
		end


		strtIx = round(length(ld.tme)*wndwSecs(2,i));

		ttl = title(axs,num2str(i),'fontsize',fs);
		ttl.Position(2) = 10;

		ylim([-70 20]); yline(-50);

		hold on;
		axis off;

		xScaleLength = 0.5; % Length of the x-axis scale bar
		yScaleLength = 40; % Length of the y-axis scale bar

		% Get current axis limits
		xlims = [totPlot_T(1),totPlot_T(end)+.5]; %xlim();
		ylims = ylim();

		% Calculate positions for scale bars
		xLblStart = xlims(2) - 0.3 * (xlims(2) - xlims(1)); % 60% from the right
		xLblHeight = -70;
		yLblStart = ylims(2) - 0.6 * (ylims(2) - ylims(1)); % 60% from the top
		yLblLeftAdjustment = xlims(1)-0.05;

		% Draw scale bars
		plot([xLblStart, xLblStart + xScaleLength], [xLblHeight, xLblHeight], 'k', 'LineWidth', 2); % X scale bar
		plot([yLblLeftAdjustment, yLblLeftAdjustment], [yLblStart, yLblStart + yScaleLength], 'k', 'LineWidth', 2); % Y scale bar

		% Labels
		if lblOnXY == 1
			text(xLblStart + xScaleLength/2, xLblHeight - 0.05*yScaleLength, [num2str(xScaleLength) ' (sec)'], 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top','fontweight','bold', 'fontsize',fs);
			text(yLblLeftAdjustment - 0.05*xScaleLength, yLblStart + yScaleLength/2, [num2str(yScaleLength) ' (mV)'], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle','fontweight','bold', 'fontsize',fs);
		end	

	end

	function fillerLxn = shader(runPatch,shades)
		a = [shades.totEK ~= -80];
		b = [diff(a);0];
		xShade1 = shades.totTme(find(b))/(60*1000);
		dmy7 = [];
		yShade1 = [-9999 -9999 9999 9999];
		fillerLxn{1} = {}; fillerLxn{2} = {};
		if ~isempty(xShade1)
			for kk = 1:length(xShade1)/2
				if runPatch == 1
					patch([xShade1(2*kk-1) xShade1(2*kk) xShade1(2*kk) xShade1(2*kk-1)], yShade1, [208 226 159]/256, 'EdgeColor', 'none','HandleVisibility','off','FaceAlpha',.5);
					hold on;
				end
				dmy = xlim;
				dmyX = get(gca,'Position');
				dmy7(kk,:) = (xShade1(2*kk-1:2*kk)-dmy(1))./(dmy(2)-dmy(1))*dmyX(3)+dmyX(1);
				%((1+((xShade1(2*kk-1:2*kk)-dmy(2))./(dmy(2)-dmy(1)))).*dmyX(3))+dmyX(1)
			end
			fillerLxn{1} = dmy7; 
		end
		a = [shades.totI ~= 0];
		b = [diff(a);0];
		xShade2 = shades.totTme(find(b))/(60*1000);
		if ~isempty(xShade2)
			for kk = 1:length(xShade2)/2
				if runPatch == 1
					patch([xShade2(2*kk-1) xShade2(2*kk) xShade2(2*kk) xShade2(2*kk-1)], yShade1,  [1.0 0.9 0.9], 'EdgeColor', 'none','HandleVisibility','off','FaceAlpha',.5);
					hold on;
				end
				dmy = xlim;
				dmyX = get(gca,'Position');
				dmy7(kk,:) = (xShade2(2*kk-1:2*kk)-dmy(1))./(dmy(2)-dmy(1))*dmyX(3)+dmyX(1);
			end
			fillerLxn{2} = dmy7;
		end
	end

end
