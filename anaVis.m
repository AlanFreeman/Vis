function anaVis

% anaVis: analyse the visual signal processing model
%
% The model computes the generator potential, p, and action potential rate, a,
% for both subcortical and cortical neurons. anaVis analyses and displays
% this activity.
% Dimensions:
%   p(t, c, z), a(t, c, z) for subcortical neurons
%   p(t, k, z), a(t, k, z) for cortical neurons
%   where:
%     t = time
%     c = channel number
%     k = location number
%     z = stage number
% Cortical stages:
%   1: inhibitory neuron soma
%   2: excitatory neuron
%   3: inhibitory neuron axon terminal
%   4: summed inhibitory input to excitatory neuron
%   5: summed subcortical input to excitatory neuron

	% Initialise
	m = setVis; % set the structural, stimulus and response parameters
	%	m.loc.loc = [-1.9, -1.2; -1.8, .9; -1.7, 1.7; -.9, -.4; -.5, -2.3; ...
	%	-.2, 2.4; .3, -2.9; .8, .8; 1.6, -2; 2.3, 2.4]; % high DSI
	%	m.loc.loc = [-2.4, -1.5; -2.4, 1.2; -.8, -1.7; -.2, -.2; .1, -1.4; ...
	%	.2, .7; .8, -2; 1.2, .1; 1.7, -1.1; 1.5, 1.7]; % low DSI
	m.loc.loc = [-.8, 0];
	w = .5 * (m.p.wid - 2); % half-width of cropped visual field (deg)
	lim = [- w, w]; tick = [- w, 0, w]; % axis limits and ticks (deg)

	% Specify tasks to be performed
	switch 'resp.dir'
		case 'cor.ind' % plot correlation between DSI maps
			m.tasks = 'corInd plot set print';
			m.corInd.file = 'Offset'; % Covariance or Offset
			m.x = 'ind'; m.y = 'resp'; m.plot.group = {'cor', 'prob'};
			m.plot.funFun = @(d, m)plot(d.ind(:), d.resp(:), 'ok');
			x = {'xLim', [0, 1], 'xTick', [0, .5, 1]}; % set x-axis
			switch m.corInd.file % set y-axis
				case 'Covariance', y = {'yLim', [-.1, 1], 'yTick', [0, .5, 1]};
				case 'Offset', y = {'yLim', [-10, 90], 'yTick', [0, 45, 90]};
			end, m.set.axes = [x, y];
			%	m.set.axes = {'xScale', 'log', 'xLim', [10 ^ -2, .5], ...
			%	'yScale', 'log', 'yLim', [3 * 10 ^ -5, .012]};
		case 'cor.x.y' % plot correlation of summed input with nearby inhibition
			switch 'mean'
				case 'cor' % raw correlations
					m.tasks = 'select resp unpack cor selDir loc desc plot add set';
					m.axes = {'x', 'y', 'dirR'};
					m.set.caxis = [-.35, .35]; m.set.colorbar = {};
				case 'dif' % correlation differences
					m.tasks = 'select resp unpack cor selDir dif loc desc plot add set';
					m.dif.mean = 0; % no mean
					m.axes = {'x', 'y'};
					m.set.caxis = [-1, 1]; m.set.colorbar = {};
				case 'mean' % mean of correlation differences
					m.tasks = 'select resp unpack cor selDir dif desc plot set';
					%	m.tasks = 'select resp unpack cor selDir dif save';
					m.dif.mean = 1; % mean difference across locations
					m.set.caxis = [0, .012]; m.set.caxis = [-.1, 1];
					m.set.colorbar = {'ticks', [0, .5, 1]};
			end
			m.select.cycle = m.p.cycles(end);
			m.resp.par = 'dir';
			dir = linspace(0, 360, 16 + 1); m.resp.dir = dir(1: end - 1);
			%	m.resp.dir = 90 + [0, 180];
			m.p.act = 0; m.p.domain = 'time';
			m.cor.group = 'dir';
			m.loc.group = 'dirR';
			m.x = 'x'; m.y = 'y';
			m.plot.funFun = @(d, m)imagesc(m.p.x, m.p.x, shiftdim(d.resp, 1)');
			m.add.funFun = @(d, m)plot(d.x, d.y, 'ok');
			m.set.axes = {'xLim', lim, 'xTick', tick, 'yLim', lim, 'yTick', tick};
			m.save.name = 'CovDif.mat';
		case 'cor.x.y.spat' % plot correlation of receptive fields
			m.tasks = 'select load corSpat loc desc plot add set';
			m.select.cycle = m.p.cycles(end);
			m.load.name = 'Field.mat';
			%	m.p.stage = 3; % correlate excitatory neuron with inhibitory neurons
			m.x = 'x'; m.y = 'y'; m.axes = {'x', 'y'};
			m.plot.funFun = @(d, m)imagesc(m.p.x, m.p.x, shiftdim(d.resp, 1)');
			m.add.funFun = @(d, m)plot(d.x, d.y, 'ok');
			m.set.axes = {'xLim', lim, 'xTick', tick, 'yLim', lim, 'yTick', tick};
			m.set.colorbar = {}; m.set.caxis = [-1, 1];
		case 'count.cor' % plot histogram of excitation/inhibition correlation
			m.tasks = ['select resp unpack corSum selDir ', ...
				'loc change hist plot set'];
			%	m.tasks = 'select resp unpack corSum selDir loc kstest list';
			m.select.cycle = m.p.cycles(end);
			m.resp.par = 'dir';
			dir = linspace(0, 360, 16 + 1); m.resp.dir = dir(1: end - 1);
			m.p.act = 0; m.p.domain = 'time';
			m.corSum.group = 'dir';
			m.loc.group = 'dirR'; m.loc.range = [-2.2, -.8; -.5, .1];
			m.kstest.samp = 'dirR';
			m.change.arg = {'m.x', 'resp', 'm.y', 'count'};
			m.hist.group = 'dirR'; m.hist.edges = linspace(-1, 1, 21);
			m.line = 'dirR';
			m.plot.funFun = @(d, m)histogram('binEdges', d.edge, ...
				'binCounts', d.count, 'displayStyle', 'stairs');
		case 'count.ind' % plot direction selectivity index histogram
			m.tasks = 'select resp amp loc fit change hist plot stop set';
			m.select.cycle = m.p.cycles(end);
			m.resp.par = 'dir';
			dir = linspace(-180, 180, 20 + 1); m.resp.dir = dir(1: end - 1);
			m.p.stage = 3; % inhibitory neurons
			m.loc.range = [-2.2, -.8; -.5, .1]; m.loc.range = [-3, -3; 3, 3];
			m.x = 'dir'; m.y = 'resp';
			m.fit.group = {'y', 'x'}; m.fit.dim = 2;
			m.change.arg = {'m.x', 'ind', 'm.y', 'count'};
			m.plot.funFun = @(d, m)histogram('binEdges', d.edge, ...
				'binCounts', d.count, 'displayStyle', 'stairs');
			m.set.axes = {'xLim', [0, 1], 'xTick', [0, .5, 1], ...
				'yLim', [0, 32], 'yTick', [0, 16, 32], 'clipping', 'off'};
		case 'count.offset' % plot histogram of exc./inh. spatial phase offset
			m.tasks = 'load unpack selDir select offset hist desc plot set print';
			m.load.name = 'Phase';
			m.unpack.par = 'dir';
			dir = linspace(0, 360, 16 + 1); m.resp.dir = dir(1: end - 1);
			m.select.dirR = 'pref';
			m.offset.crop = 61; % cropped size (x-locations)
			m.x = 'resp'; m.y = 'count';
			m.hist.edges = -10: 5: 90;
			m.plot.funFun = @(d, m)histogram('binEdges', d.edge, ...
				'binCounts', d.count, 'displayStyle', 'stairs');
			m.set.axes = {'xLim', [- 10, 90], 'xTick', [0, 45, 90], ...
				'yLim', [0, 1500], 'yTick', [0, 700, 1400]};
		case 'dir.x.y' % plot map of preferred direction
			m.tasks = 'load dir desc plot set';
			m.load.name = 'Index';
			m.x = 'x'; m.y = 'y';
			m.plot.funFun =	@(d, m)quiver(d.x, d.y, real(d.dirP), imag(d.dirP), 'k');
			%	'k', 'lineWidth', 1, 'maxHeadSize', .25);
			m.set.axes = {'xLim', lim, 'xTick', tick, 'yLim', lim, 'yTick', tick};
		case 'ind.x.y' % plot map of direction selectivity index
			m.tasks = 'select resp amp loc fit pack plot set';
			m.tasks = 'select resp amp loc fit pack save plot set';
			m.tasks = 'load plot set';
			m.select.cycle = m.p.cycles(end);
			m.resp.par = 'dir';
			dir = linspace(-180, 180, 20 + 1); m.resp.dir = dir(1: end - 1);
			%	m.p.act = 0; % m.p.stage = 3;
			m.loc.range = [-3, -3; 3, 3]; % m.loc.all = 1;
			m.fit.group = {'y', 'x'}; m.fit.dim = 2;
			m.x = 'x'; m.y = 'y';
			m.plot.funFun = @(d, m)imagesc(d.(m.x), d.(m.y), shiftdim(d.ind, 1)');
			m.set.axes = {'xLim', lim, 'xTick', tick, 'yLim', lim, 'yTick', tick};
			m.set.colorbar = {}; m.set.caxis = [0, 1];
			m.save.name = 'Index.mat'; m.load.name = m.save.name;
		case 'ind.x.y.pred' % plot DSI map predicted from receptive field correl.
			m.tasks = 'select load corSpat ind desc plot set';
			m.select.cycle = m.p.cycles(end);
			m.load.name = 'Field.mat';
			m.x = 'x'; m.y = 'y';
			m.plot.funFun = @(d, m)imagesc(m.p.x, m.p.x, shiftdim(d.resp, 1)');
			m.set.axes = {'xLim', lim, 'xTick', tick, 'yLim', lim, 'yTick', tick};
			m.set.colorbar = {}; % m.set.caxis = [-1, 1];
		case 'offset.x.y' % plot map of exc./inh. spatial phase offset
			m.tasks = 'load unpack selDir select offset desc plot set';
			m.load.name = 'Phase'; m.save.name = 'Offset';
			m.unpack.par = 'dir';
			dir = linspace(0, 360, 16 + 1); m.resp.dir = dir(1: end - 1);
			m.select.dirR = 'pref';
			m.x = 'x'; m.y = 'x';
			m.plot.funFun = @(d, m)imagesc(d.x, d.x, shiftdim(d.resp, 1)');
			m.set.axes = {'xLim', lim, 'xTick', tick, 'yLim', lim, 'yTick', tick};
			m.set.colorbar = {'xLim', [-10, 50], 'xTick', [0, 25, 50]};
			m.set.caxis = [-10, 50];
		case 'orient.x.y' % plot map of preferred orientation or direction
			m.tasks = 'select orient desc plot set';
			m.select.cycle = m.p.cycles(end);
			m.orient.orient = 1; % 0 for direction, 1 for orientation
			m.x = 'x'; m.y = 'y';
			m.plot.funFun = @(d, m)imagesc(d.x, d.y, shiftdim(d.dirP, 1)');
			m.set.axes = {'xLim', lim, 'xTick', tick, 'yLim', lim, 'yTick', tick};
			m.set.colormap = 'hsv';
			switch m.orient.orient % direction or orientation?
				case 0 % direction
					m.set.caxis = [0, 360];
					m.set.colorbar = {'xLim', [0, 360], 'xTick', [0, 180, 360]};
				case 1 % orientation
					m.set.caxis = [0, 180];
					m.set.colorbar = {'xLim', [0, 180], 'xTick', [0, 90, 180]};
			end
		case 'resp.dir' % plot response amplitude vs direction
			m.tasks = 'select resp amp loc desc plot fit refine pred add set';
			m.select.cycle = m.p.cycles(end);
			m.resp.par = 'dir';
			dir = linspace(-180, 180, 20 + 1); m.resp.dir = dir(1: end - 1);
			%	m.p.stage = 3; % m.p.act = 0; 
			m.x = 'dir'; m.y = 'resp';
			g = {'x', 'y'}; m.axes = g;
			m.plot.arg = {'lineStyle', 'none', 'marker', 'o', ...
				'markerEdgeColor','k', 'markerFaceColor','k', 'markerSize', 6};
			m.set.axes = {'xLim', [-180, 180], 'xTick', [-180, 0, 180], ...
				'yLim', [0, 20], 'yTick', [0, 10, 20]};
			%	'yLim', [0, 60], 'yTick', [0, 30, 60]}; % inhibitory stage
			m.fit.group = g; m.refine.group = g; m.pred.group = g;
			m.plot.colourOrder = [0, 0, 0]; % black line
		case 'resp.dir.half' % plot resp. amp. vs direction for half-field inhib.
			m.tasks = ...
				'select resp inhib ampNew desc plot fit refine pred add set print';
			m.select.cycle = m.p.cycles(end);
			m.resp.par = {'dir', 'tauI'};
			dir = linspace(-180, 180, 20 + 1); m.resp.dir = dir(1: end - 1);
			m.resp.tauI = [.1, 1000];
			m.p.act = 0; % generator potential
			m.x = 'dir'; m.y = 'resp';
			g = {'side'}; m.plot.group = g;
			m.plot.arg = {'lineStyle', 'none', 'marker', 'o', ...
				'markerEdgeColor','k', 'markerFaceColor','k', 'markerSize', 6};
			m.set.axes = {'xLim', [-180, 180], 'xTick', [-180, 0, 180], ...
				'yLim', [0, 20], 'yTick', [0, 10, 20]};
			m.fit.group = g; m.refine.group = g; m.pred.group = g;
		case 'resp.dir.inh' % plot response amplitude vs direction vs inhibition
			m.tasks = 'select resp amp loc unpack desc plot fit refine pred add set';
			m.select.cycle = m.p.cycles(end);
			m.resp.par = {'dir', 'tauI'};
			dir = linspace(-180, 180, 20 + 1); m.resp.dir = dir(1: end - 1);
			m.resp.tauI = [.1, 1000];
			m.unpack.par = 'tauI';
			m.x = 'dir'; m.y = 'resp';
			g = {'x', 'y', 'tauI'}; m.axes = g;
			m.plot.arg = {'lineStyle', 'none', 'marker', 'o', ...
				'markerEdgeColor','k', 'markerFaceColor','k', 'markerSize', 6};
			m.set.axes = {'xLim', [-180, 180], 'xTick', [-180, 0, 180], ...
				'yLim', [0, 20], 'yTick', [0, 10, 20]};
			%	'yLim', [0, 60], 'yTick', [0, 30, 60]};
			m.fit.group = g; m.refine.group = g; m.pred.group = g;
		case 'resp.phase' % plot response amplitude or phase vs spatial phase
			m.tasks = 'select resp amp save';
			m.tasks = 'load unpack selDir select loc phase desc plot set';
			%	m.tasks = 'load unpack selDir select loc phase desc plot set add';
			m.select.cycle = m.p.cycles(end);
			m.resp.par = {'phaseS', 'dir'};
			p = linspace(0, 360, 16 + 1); m.resp.phaseS = p(1: end - 1);
			dir = linspace(0, 360, 16 + 1); m.resp.dir = dir(1: end - 1);
			m.p.stim = 'reverse'; m.p.act = 0;
			m.p.comp = 'complex'; m.p.stage = [];
			m.save.name = 'Phase.mat'; m.load.name = m.save.name;
			m.unpack.par = 'dir';
			m.select.dirR = 'pref'; m.loc.group = 'dirR';
			m.phase.group = {'x', 'y', 'dirR'};
			m.phase.comp = 'amp'; m.phase.stage = 4; m.phase.shift = 1;
			m.plot.group = {'dirR'}; m.line = {'x', 'y'};
			m.x = 'phaseS'; m.y = 'resp';
			a = {'xLim', [0, 360], 'xTick', [0, 180, 360]};
			switch m.phase.comp % which component?
				case 'amp' % amplitude vs spatial phase
					m.set.axes = [a, {'yLim', [0, 16], 'yTick', [0, 8, 16], ...
						'clipping', 'off'}];
				case 'complex' % polar plot
					m.x = 'resp'; m.y = 'resp';
					m.plot.funFun = @(d, m)plot(real(d.resp), imag(d.resp), ...
						'lineStyle', 'none', 'marker', 'o', ...
						'markerEdgeColor','k', 'markerFaceColor','k', 'markerSize', 6);
					%	m.plot.funFun = @(d, m)plot([0, 0; real(d.respComp)], ...
					%	[0, 0; imag(d.respComp)]);
					m.set.axes = {'xLim', [-15, 15], 'xTick', [-15, 0, 15], ...
						'yLim', [-15, 15], 'yTick', [-15, 0, 15]};
					m.add.funFun = @(d, m)plot(real(d.respSim), imag(d.respSim));
				case 'phase' % response phase vs spatial phase
					m.set.axes = [a, {'yLim', [0, 360], 'yTick', [0, 180, 360]}];
					%	m.set.axes = [a, {'yLim', [-360, 360], 'yTick', [-360, 0, 360]}, ...
					%	'clipping', 'off'];
					%	m.set.axes = [a, {'yLim', [150, 540], 'yTick', [180, 360, 540]}];
			end
		case 'resp.t' % plot response time course
			m.tasks = 'select resp loc unpack select desc plot set';
			m.select.cycle = m.p.cycles(end);
			m.resp.par = 'dir'; m.resp.dir = -54 + [0, 180]; m.resp.dir = [-56, 103];
			%	m.p.stim = 'reverse';
			%	m.p.solver = 'solveT'; m.p.time = 1.5;
			m.p.act = 0; m.p.domain = 'time'; % m.p.tauI = 1000;
			m.unpack.group = {'x', 'y'};
			m.unpack.par = {'dir', 'stage'};
			m.select.stage = [2, 4, 5]; % m.select.stage = 1: 3;
			m.x = 'time'; m.y = 'resp';
			m.plot.group = {'x', 'y', 'dir'}; m.plot.line = {'stage'};
			m.plot.funFun = @(d, m)plot(d.time, shiftdim(d.resp, 1));
			m.set.axes = {'xLim', [0, .5], 'xTick', [0, .25, .5], ...
				'yLim', [-30, 30], 'yTick', [-30, 0, 30], 'clipping', 'off'}; % pot
			%	'yLim', [0, 60], 'yTick', [0, 30, 60], 'clipping', 'off'}; % act
			%	m.set.axes = {'xLim', [1, 1.5], 'xTick', [1, 1.25, 1.5], ...
			%	'yLim', [-30, 30], 'yTick', [-30, 0, 30]}; % pot
		case 'resp.t.comp' % plot response time course compiled from reversing grat.
			m.tasks = 'select resp comp loc unpack select desc plot set';
			m.select.cycle = m.p.cycles(end);
			m.p.stim = 'reverse'; m.p.act = 0;
			m.resp.par = {'dir', 'phaseS'};
			m.resp.dir = [-90, 90]; m.resp.phaseS = [0, 90];
			m.unpack.par = {'dir', 'stage'}; m.resp.stage = 1: 5;
			m.select.stage = [2, 4, 5]; % m.select.stage = 1: 3;
			m.x = 'time'; m.y = 'resp';
			m.axes = {'x', 'y', 'dir'}; m.line = 'stage';
			m.plot.funFun = @(d, m)plot(d.time, shiftdim(d.resp, 1));
			m.set.axes = {'xLim', [0, .5], 'xTick', [0, .25, .5], ...
				'yLim', [-30, 30], 'yTick', [-30, 0, 30]};
		case 'resp.t.half' % plot response time course for half-field inhibition
			m.tasks = 'select resp inhib unpack select desc plot set';
			m.select.cycle = m.p.cycles(end);
			m.resp.group = 'cycle';
			m.resp.par = {'dir', 'tauI'};
			m.resp.dir = [-90, 90];
			m.resp.tauI = [.1, 1000];
			m.p.act = 0; % generator potential
			m.unpack.group = 'side';
			m.unpack.par = {'dir', 'stage'};
			m.resp.stage = 1: 5;
			m.select.stage = [2, 4, 5];
			switch 'resp'
				case 'field' % plot preferred half-field
					m.x = 'x'; m.y = 'x'; m.axes = {'cycle', 'dir', 'side'};
					m.plot.funFun = @(d, m)imagesc(d.(m.x), d.(m.y), ...
						shiftdim(d.field(1, :, :), 1)');
					m.set.colorbar = {}; m.set.colormap = 'gray';
					w = .5 * m.p.wid; % half-width of visual field
					m.set.axes = {'xLim', [- w, w], 'yLim', [- w, w], 'clipping', 'off'};
				case 'resp'
					m.x = 'time'; m.y = 'resp';
					m.axes = {'cycle', 'dir', 'side'}; m.line = 'stage';
					m.plot.funFun = @(d, m)plot(d.time, shiftdim(d.resp, 1));
					m.set.axes = {'xLim', [0, .5], 'xTick', [0, .25, .5], ...
						'yLim', [-30, 30], 'yTick', [-30, 0, 30], 'clipping', 'off'}; % pot
			end
		case 'resp.t.inh' % plot response time course: dynamic and static inhibition
			m.tasks = 'select resp loc unpack select desc plot set';
			m.select.cycle = m.p.cycles(end);
			m.resp.group = 'cycle';
			m.resp.par = {'dir', 'tauI'};
			m.resp.dir = [-90, 90]; m.resp.tauI = [.1, 1000];
			% m.p.act = 0; m.p.cort = 0; m.p.solver = 'solveT'; m.p.stim = 'reverse';
			m.p.act = 0; m.p.domain = 'time';
			m.unpack.par = {'dir', 'stage', 'tauI'}; m.resp.stage = 1: 5;
			m.select.stage = [2, 4, 5]; % m.select.stage = 1: 3;
			m.x = 'time'; m.y = 'resp';
			m.axes = {'cycle', 'x', 'y', 'dir', 'tauI'}; m.line = {'stage'};
			m.plot.funFun = @(d, m)plot(d.time, shiftdim(d.resp, 1));
			m.set.axes = {'xLim', [0, .5], 'xTick', [0, .25, .5], ...
				'yLim', [-30, 30], 'yTick', [-30, 0, 30], 'clipping', 'off'}; % pot
			%	'yLim', [0, 60], 'yTick', [0, 30, 60], 'clipping', 'off'}; % act
			%	m.set.axes = {'xLim', [1, 1.5], 'xTick', [1, 1.25, 1.5], ...
			%	'yLim', [-30, 30], 'yTick', [-30, 0, 30]}; % pot
		case 'resp.t.phase' % plot response time course vs spatial phase
			m.tasks = 'select resp loc unpack select desc plot set';
			m.select.cycle = m.p.cycles(end);
			m.resp.par = 'phaseS';
			p = linspace(0, 180, 8 + 1); m.resp.phaseS = p(1: end - 1);
			m.p.dir = 90;
			m.p.act = 0; m.p.domain = 'time'; m.p.stim = 'reverse';
			m.unpack.group = {'x', 'y'};
			m.unpack.par = {'phaseS', 'stage'}; %	m.resp.stage = 1: 5;
			m.select.stage = [2, 4, 5]; % m.select.stage = 1: 3;
			m.x = 'time'; m.y = 'resp';
			m.axes = {'x', 'y', 'phaseS'}; m.line = {'stage'};
			m.plot.funFun = @(d, m)plot(d.time, shiftdim(d.resp, 1));
			m.set.axes = {'xLim', [0, .5], 'xTick', [0, .25, .5], ...
				'yLim', [-30, 30], 'yTick', [-30, 0, 30], 'clipping', 'off'}; % pot
			%	'yLim', [0, 60], 'yTick', [0, 30, 60], 'clipping', 'off'}; % act
		case 'resp.t.rev' % plot time course for half-field inhib., reversing grat.
			m.tasks = 'select resp inhib unpack select desc plot set';
			m.select.cycle = m.p.cycles(end);
			m.resp.par = {'phaseS', 'tauI'};
			m.resp.phaseS = [-90, -45, 0, 45]; m.resp.tauI = [.1, 1000];
			m.p.dir = 90; m.p.stim = 'reverse'; m.p.act = 0;
			m.unpack.group = 'side';
			m.unpack.par = {'phaseS', 'stage'};
			m.resp.stage = 1: 5;
			m.select.stage = 4;
			m.x = 'time'; m.y = 'resp';
			m.axes = {'side'}; m.line = 'phaseS';
			m.plot.funFun = @(d, m)plot(d.time, shiftdim(d.resp, 1));
			m.set.axes = {'xLim', [0, .5], 'xTick', [0, .25, .5], ...
				'yLim', [-30, 30], 'yTick', [-30, 0, 30], 'clipping', 'off'}; % pot
		case 'resp.t.sparse' % plot response time course for sparse noise
			m.tasks = 'select resp loc unpack select desc plot set';
			m.select.cycle = m.p.cycles(end);
			m.resp.par = 'cont'; m.resp.cont = [-1, 1];
			m.p.stim = 'sparse';
			m.p.locE = m.loc.loc;
			m.p.act = 0; m.p.domain = 'time';
			%	m.p.solver = 'solveF'; m.p.time = .5; %	m.p.error = inf;
			m.p.solver = 'solveT'; m.p.time = .2;
			m.unpack.group = {'x', 'y'}; m.unpack.par = {'cont', 'stage'};
			m.select.stage = [2, 4, 5]; % m.select.stage = 1: 3;
			m.x = 'time'; m.y = 'resp';
			m.plot.group = {'x', 'y', 'cont'}; m.plot.line = 'stage';
			m.set.axes = {'xLim', [0, .2], 'xTick', [0, .1, .2], ...
				'yLim', [-32, 32], 'yTick', [-30, 0, 30], 'clipping', 'off'}; % pot
		case 'resp.x.y' % plot receptive field
			%	m.tasks = 'select field desc save';
			m.tasks = 'load loc map desc plot add set';
			m.select.cycle = m.p.cycles(end);
			m.p.stim = 'sparse';
			%	m.p.solver = 'solveF'; m.p.time = .5; m.p.error = inf;
			m.p.solver = 'solveT'; m.p.time = .2;
			m.save.name = 'Field.mat'; m.load.name = m.save.name;
			m.loc.group = 'cont';
			m.p.stage = 4; m.p.act = 0;
			m.map.dif = 1; % 0 for separate on/off responses, 1 for on/off difference
			m.map.contour = 0; % 0 for image, 1 for contour plot
			m.x = 'xStim'; m.y = 'yStim'; m.plot.group = {'x', 'y'};
			if m.map.contour % contour plot
				m.map.group = {'cont', 'x', 'y'};
				m.plot.line = {'cont', 'max'}; % show on and off resp. separately
				m.plot.funFun = @(d, m)contour(d.xStim, d.xStim, ...
					shiftdim(d.resp, 1)', d.contour, d.colour); % contours
			else % image plot
				m.map.group = {'x', 'y'};
				m.plot.funFun = @(d, m)imagesc(d.xStim, d.xStim, shiftdim(d.resp, 1)');
				m.set.colorbar = {}; m.set.colormap = 'jet'; m.set.caxis = [-20, 20];
			end
			m.set.axes = {'xLim', lim, 'xTick', tick, 'yLim', lim, 'yTick', tick};
			m.add.funFun = @(d, m)plot(d.x, d.y, 'ok');
		case 'stim.x.y' % show movie of stimulus
			m.tasks = 'stim select plot set';
			%	m.tasks = 'stim plot movie set';
			m.p.dir = 270; m.p.phaseS = 90; m.p.stim = 'reverse'; m.p.time = .5;
			m.select.t = 0;
			m.x = 'x'; m.y = 'y';
			m.plot.funFun = @(d, m)imagesc(d.x(1, :), d.y(1, :), ...
				shiftdim(d.stim(1, :, :), 1)');
			m.set.axes = {'xLim', lim, 'xTick', tick, 'yLim', lim, 'yTick', tick};
			m.set.colormap = 'gray'; m.set.colorbar = {}; m.set.caxis = [-1, 1];
		case 'weight' % calculate weights remotely
			m.tasks = 'weight';
			%	m.p.gpu = 1; % use GPU
			%	m.p.parallel = 1;	% parallel processing
		case 'weight.x.y' % plot channel weight versus location
			m.tasks = 'select chan desc plot add set';
			m.select.cycle = [0, m.p.cycles(end)]; % m.select.cycle = 0;
			m.chan.group = 'cycle';
			m.x = 'x'; m.y = 'y';
			m.axes = 'cycle'; m.line = 'sign'; % m.plot.subplot = [1, 1];
			m.plot.funFun = @(d, m)scatter(d.x, d.y, max(1, 12 * d.weight), 'filled');
			m.set.axes = {'xLim', lim, 'xTick', tick, 'yLim', lim, 'yTick', tick};
			m.set.legend = {'off'};
			m.add.funFun = @(d, m)plot(m.loc.loc(1), m.loc.loc(2), 'ok');
	end
	
	% Set up plotting for export
	%	m.set.axesLength = [5, 5]; m.set.lineWidth = 1.5; % for uniform plotting
	
	% Obtain handles for functions in calVis
	m = calVis(m);
	
	% Define handles for task functions in this file
	switch 1, case 1 % hide with code folding
		m.amp.fun = @doAmp;
		m.ampNew.fun = @doAmpNew;
		m.chan.fun = @doChan;
		m.comp.fun = m.calComp;
		m.corInd.fun = @doCorInd;
		m.desc.fun = @doDesc;
		m.dir.fun = @doDir;
		m.field.fun = m.calField;
		m.fit.fun = @doFit;
		m.hist.fun = @doHist;
		m.ind.fun = m.calInd;
		m.loc.fun = @doLoc;
		m.kstest.fun = @doKstest;
		m.map.fun = @doMap;
		m.mod.fun = @doMod;
		m.movie.fun = @doMovie;
		m.orient.fun = @doOrient;
		m.pack.fun = @doPack;
		m.print.fun = @doPrint;
		m.refine.fun = @doRefine;
		m.resp.fun = m.calResp;
		m.stim.fun = @doStim;
		m.unpack.fun = @doUnpack;
		m.weight.fun = @doWeight;
	end
	
	% Analyse
	m = m.calChan(m); % calculate channel locations
	m = m.calGain(m); % calculate geniculocortical gain
	m = m.calTemp(m); % calculate temporal parameters
	[d, m] = m.calWeight(m); % calculate geniculocortical weight
	anaTab(d, m); % run the analysis tasks

function [d, m] = doAmp(d, m)

% Calculate response fundamental component at specified stage

	% Select response component and, possibly, stage
	r = d.resp; % response: 1 x fs x ks x zs x ps, ps can be multi-dimensional
	s = size(r); % size
	dim = d.Properties.CustomProperties.RespDim; % response dimensions
	if isempty(m.p.stage) % select fundamental component
		r = r(:, 2, :, :, :); % fundamental component: 1 x 1 x ks x zs x ps
		r = reshape(r, [1, s(3: end)]); % prepare for storage: 1 x ks x zs x ps
		dim(2) = []; % remove dimension
	else % also select stage
		r = r(:, 2, :, m.p.stage, :); % fundamental component: 1 x 1 x ks x 1 x ps
		r = reshape(r, [1, s(3), s(5: end)]); % prepare for storage: 1 x ks x ps
		dim([2, 4]) = [];
	end
	d.Properties.CustomProperties.RespDim = dim; % update
	
	% Calculate response metric
	switch m.p.comp % response metric?
		case 'amp' % amplitude
			d.resp = 2 * abs(r) / m.p.fs; % amplitude (mV or Hz): 1 x ks x ps
		case 'complex'
			d.resp = 2 * r / m.p.fs; % leave as complex (mV or Hz): 1 x ks x ps
		case 'phase'
			d.resp = 180 * angle(r) / pi; % phase (deg): 1 x ks x ps
			d.Properties.VariableDescriptions{'resp'} = 'Response phase (deg)';
	end
	
function [d, m] = doAmpNew(d, m)

% Calculate response amplitude for cortical cells at specified stage

	% Convert response to impulse rate frequency response
	r = d.resp; % time course of potential (mV)
	r = m.p.gainGen * max(0, r); % convert to impulse rate (Hz)
	m.p.act = 1; % update the flag
	r = fft(r, [], 2); % convert to frequency response
	
	% Shift frequency and stage parameters to left of response array
	dim = d.Properties.CustomProperties.RespDim; % response dimensions
	[~, i] = ismember({'freq', 'stage'}, dim); % dimensions of interest
	iO = 1: ndims(r); iO(i) = []; % other dimensions
	s = size(r); % response array size
	s(i) = []; % size of other dimensions
	r = permute(r, [i, iO]); % shift dimensions of interest to left
	
	% Calculate response amplitude and store
	r = abs(r(2, m.p.stage, :)); % amplitude of excitatory fundamental component
	r = 2 * r / m.p.fs; % convert amplitude from fft to physical units
	r = reshape(r, [s, 1]); % prepare for storage
	d.resp = r; % store
	d.Properties.CustomProperties.RespDim = dim(iO); % update
	
function [d, m] = doChan(d, m)

% Store channel data in data table

	% Make one row for each channel
	l = m.p.chan; % channel data
	cs = size(l, 1); % number of channels
	w = d.weight{1}; % geniculocortical weights
	d.weight = []; % don't want to replicate this: too large
	d = repmat(d, [cs, 1]); % one row for each channel
	d.x = l(:, 1); d.y = l(:, 2); d.sign = l(:, 3);
	d.Properties.VariableDescriptions{'x'} = 'Horizontal location (deg)';
	d.Properties.VariableDescriptions{'y'} = 'Vertical location (deg)';
	
	% Select channel closest to visual field location m.loc.loc
	if m.p.cort == 0 % obtain list of locations of channels or cortical neurons
		l = m.p.chan; % subcortical channel location
	else
		l = m.p.loc; % cortical location
	end
	i = knnsearch(l(:, 1: 2), m.loc.loc); % location closest to m.loc.loc
	d.weight = w(:, i); % select
	d.Properties.VariableDescriptions{'weight'} = 'Synaptic weight';

function [d, m] = doCorInd(~, m)

% Calculate correlation between DSI maps

	load(m.corInd.file, 'd'); rC = d.resp; % load covariance difference or
		%	spatial phase offset: 1 x xsC x ysC
	[~, xsC, ~] = size(rC); % number of x values
	load('Index', 'd'); rT = d.ind; % tuning DSI: 1 x xs x xs
	[~, xs, ~] = size(rT); % number of x values
	x = round(.5 * (xsC - xs)); % number of x values to crop at each edge
	rC = rC(1, 1 + x: xsC - x, 1 + x: xsC - x); % cropped image: 1 x xs x xs
	[c, p] = corrcoef(rT, rC, 'rows', 'complete'); % corr. coeff. and prob.
	i = 1: 2: xs; rC = rC(1, i, i); rT = rT(1, i, i); % subsample
	d.resp = rC; d.ind = rT; % store for plotting
	d.cor = c(2); d.prob = p(2); % store

function [d, m] = doDesc(d, m)

% Describe the variables in the data table

	desc = d.Properties.VariableDescriptions; % variable descriptions
	undesc = find(strcmp('', desc)); % indices of undescribed variables
	for i = undesc % loop over undescribed variables
		name = d.Properties.VariableNames{i}; % name of current variable
		switch name % describe this variable
			case 'cont', desc{i} = 'Contrast';	
			case 'dir', desc{i} = 'Stimulus direction (deg ACW from rightward)';
			case 'phaseS', desc{i} = 'Spatial phase (deg)';
			case 'resp'
				if m.p.act, desc{i} = 'Impulse rate (Hz)';
				else, desc{i} = 'Generator potential (mV)'; end
			case 'time', desc{i} = 'Time (s)';
			case 'x', desc{i} = 'Horizontal location (deg)';
			case 'y', desc{i} = 'Vertical location (deg)';
		end
	end
	d.Properties.VariableDescriptions = desc;

function [d, m] = doDir(d, m)

% Prepare map of preferred direction

	% Calculate map of preferred direction
	s = 21; % map size (elements)
	x = linspace(d.x(1), d.x(end), s); % x values (deg)
	[x, y] = ndgrid(x); % x and y values at which to plot (deg): xs x ys
	ind = d.ind; % direction selectivity index: 1 x xs x ys
	dirP = pi * d.dirP / 180; % preferred direction (radians): 1 x xs x ys
	dirP = exp(1i * dirP); % make it a unit vector: 1 x xs x ys
	dirP = ind .* dirP; % set the magnitude to the direction selectivity index
	dirP = shiftdim(dirP, 1); % convert to image: xs x ys
	dirP = imresize(dirP, [s, s]); % resize to reduce the number of elements
	d = repmat(d, [s, 1]); % one row for each element
	d.x = x; d.y = y; d.dirP = dirP; % store

function [d, m] = doFit(d, m)

% Fit model to direction tuning curve

	% Prepare for regression
	x = d.dir'; % predictor: use column vector
	y = d.resp'; % response
	fun = m.modelTuning; % model to fit
	[yMax, i] = max(y); % estimate model coefficients
	dirP = x(i); % preferred direction (deg)
	i = i - 2: i + 2; i = mod(i - 1, length(x)) + 1; % indices around maximum
	yA = y; yA(i) = 0; % response without maximum
	[~, i] = max(yA); % maximum of antipreferred direction
	dirA = x(i); % antipreferred direction (deg)
	b0 = [0, yMax, .5 * yMax, 10, dirP, dirA]; % initial estimate of coefficients
	cName = {'r0', 'r1', 'r2', 'band', 'dirP', 'dirA'}; % coefficient names
	arg = {'coefficientNames', cName}; % fitnlm arguments
	warning('off', 'MATLAB:rankDeficientMatrix'); % hide warnings
	warning('off', 'stats:nlinfit:ModelConstantWRTParam');
	warning('off', 'stats:nlinfit:IllConditionedJacobian');
	warning('off', 'stats:nlinfit:IterationLimitExceeded');
	lastwarn(''); % remove warning in preparation for next
	
	% Prepare output file
	c = nan(1, length(cName)); % coefficient values are missing
	dC = array2table(c, 'variableNames', cName);
	dC.rSq = nan; % adjusted r-squared
	dC.ind = nan; % direction selectivity index
	d.dirP = []; % remove preferred direction used for calculating weights
	dC = [d, dC]; % concatenate with input file
	dC.Properties.VariableDescriptions{'r1'} = 'Response amplitude (Hz)';
	dC.Properties.VariableDescriptions{'band'} = 'Direction bandwidth (deg)';
	dC.Properties.VariableDescriptions{'dirP'} = 'Preferred direction (deg)';
	dC.Properties.VariableDescriptions{'rSq'} = 'Adjusted r-squared';
	dC.Properties.VariableDescriptions{'ind'} = 'Direction selectivity index';
	switch 'nan' % make default file in case regression doesn't work
		case 'empty', d = []; % old method
		case 'nan', d = dC;
	end
	
	% Fit
	try % try fitting
		model = fitnlm(x, y, fun, b0, arg{:});
	catch % there was an error
		return % use default file
	end
	m.model{m.group} = model;
	
	% Process warnings
	[~, w] = lastwarn; % find the warning ID
	if ~ isempty(w) % there was a warning
		%	return % use default file
	end
	
	% Store
	if isfield(m.fit, 'dim') && m.fit.dim == 2 % store coefficients
		
		% Make table of regression coefficients
		dC(:, cName) = array2table(model.Coefficients.Estimate'); % make it a table
		k = dC.band; % coefficient of exponent
		dC.band = acosd(1 - log(2) / dC.band); % half-width at half-height (deg)
		dC.rSq = model.Rsquared.Adjusted;
		
		% Fix coefficient values
		r1 = dC.r1; r2 = dC.r2; dirP = dC.dirP; dirA = dC.dirA;
		if r1 < r2 % r1 is assumed to be preferred direction
			r1 = r2; r2 = dC.r1; % switch amplitudes
			dirP = dirA; dirA = dC.dirP; % switch directions
		end
		if r1 < 0 % can't fix this
			return % use default file
		end
		if r2 < 0 % presumably results from noise
			r2 = 0;
		end
		if ~ isreal(dC.band) % can't fix this
			return % use default file
		end
		switch 'anti' % choose response for direction selectivity index
			case 'anti' % use antipreferred response
				r = r1 * exp(k * (-2)) + ...
					r2 * exp(k * (cos((pi / 180) * (dirP + 180 - dirA)) - 1));
			case 'sub', r = r2; % use suboptimal peak
		end
		dC.ind = (r1 - r) / (r1 + r); % direction selectivity index
		dC.r1 = r1; dC.r2 = r2; % amplitudes
		dC.dirP = mod(dirP + 180, 360) - 180; % range is [-180, 180)
		dC.dirA = mod(dirA + 180, 360) - 180;
		
	else % store x and y as columns, for model prediction
		dC = repmat(dC(1, :), size(x)); % one row per observation
		dC.(m.x) = x; dC.(m.y) = y;
	end
	d = dC; % success

function [d, m] = doHist(d, m)

% Calculate histogram

	if isfield(m.hist, 'edges')
		[n, e] = histcounts(d.(m.x), m.hist.edges); % histogram counts
	else
		[n, e] = histcounts(d.(m.x), 10);
	end
	w = e(2) - e(1); % bin width
	c = .5 * w + e(1: end - 1); % bin centres
	d = d(1, :); % keep only the first row
	d.(m.x) = c; % bin centres
	d.edge = e; % bin edges
	d.(m.y) = n; % counts
	d.Properties.VariableDescriptions{'edge'} = 'Edge';
	d.Properties.VariableDescriptions{m.y} = 'Count';

function [d, m] = doKstest(d, m)

% Perform two-sample Kolmogorov-Smirnov test

	g = d.(m.kstest.samp); % grouping variable
	i = g == g(1); % i is true for first sample, false for second
	x1 = d.resp(i); x2 = d.resp(~ i); % two samples
	[~, p] = kstest2(x1, x2); % two-sample Kolmogorov-Smirnov test
	d = d(1, :); % keep only first line
	d.prob = p; % store

function [d, m] = doLoc(d, m)

% Select visual field locations and unpack them

	% Prepare response variable
	y = d.resp; % response variable
	y = shiftdim(y, 1); % remove leading single-element dimension
	dim = d.Properties.CustomProperties.RespDim; % response dimensions
	[~, i] = ismember('loc', dim); % index of location in response
	i = i - 1; % dimension representing location
	y = permute(y, [i, 1: i - 1, i + 1: ndims(y)]); % shift it to first
	s = size(y); % size
	
	% Specify locations to select
	if m.p.cort, l = m.p.loc; % cortical locations
	else, l = m.p.chan; % subcortical locations
	end % l is list of all locations
	if isfield(m.loc, 'all') % find locations
		loc = l; % all locations
	elseif isfield(m.loc, 'range') % range of locations
		loc = m.loc.range; % lower left and upper right corners
		locs = round(1 + m.p.densC * (loc(2, :) - loc(1, :))); % number of locations
		[locX, locY] = ndgrid(linspace(loc(1, 1), loc(2, 1), locs(1)), ...
			linspace(loc(1, 2), loc(2, 2), locs(2))); % locations as grid
		loc = [locX(:), locY(:)]; % locations as list
	else % list of individual locations
		loc = m.loc.loc; % locations to select
	end
	
	% Select locations
	if ~ isfield(m.loc, 'all') % subset of locations
		i = knnsearch(l(:, 1: 2), loc); % location closest to loc
		loc = l(i, :); % use model locations, not user estimates
		y = y(i, :); % keep only selected locations
	end
	locs = size(y, 1); % number of locations
	s(1) = locs; y = reshape(y, s); % restore shape
	
	% Unpack
	desc = d.Properties.VariableDescriptions{'resp'}; % response description
	d.resp = []; % remove during replication, restore below
	d = repmat(d, [locs, 1]); % one row for each location
	d.resp = y; % store response
	d.Properties.VariableDescriptions{'resp'} = desc; % restore description
	d.x = loc(:, 1); d.y = loc(:, 2); % store location
	[~, i] = ismember('loc', dim); % index of location dimension
	dim(i) = []; d.Properties.CustomProperties.RespDim = dim; % update dimensions
	
function [d, m] = doMap(d, m)

% Prepare for receptive field map

	r = d.resp; % generator potential (mV): cs x xs x ys x zs
	r = r(:, :, :, m.p.stage); % select required stage: cs x xs x ys
	if m.p.act % convert from generator potential to impulse rate
		r = m.p.gainGen * max(0, r); % impulse rate (Hz)
	end
	if m.map.dif % take difference of on- and off-responses
		r = r(2, :, :) - r(1, :, :); % difference
		d = d(1, :); % keep only first line
	end
	rMax = max(r(:)); % store maximum
	if m.p.act && rMax == 0 % impulse rate is zero
		d = []; return % remove row
	end
	if m.map.contour % contour map
		r = r / rMax; % normalise response
		d.contour = [.05, .5, .95]; % contour levels
		if d.cont == -1, d.colour = 'b'; % dark stimulus
		else, d.colour = 'r'; % light stimulus
		end
	else % image
		r = shiftdim(r, 1); % xs x ys
		%	r = imgaussfilt(r); % smooth
		r = shiftdim(r, -1); % 1 x xs x ys
	end
	d.resp = r; d.max = rMax; % store
	
function [d, m] = doMod(d, m)

% Calculate relative modulation = fundamental response amplitude / mean resp.

	% Calculate
	[d, m] = m.calResp(d, m); % calculate response: 1 x fs x ks x zs x dirs
	r = d.resp; % response
	[~, fs, ks, ~, ~] = size(r); % size
	r1 = 2 * abs(squeeze(r(:, 2, :, 2, :))) / fs;
		% fundamental amplitude of excitatory cells: ks x dirs
	[r1, i] = max(r1, [], 2); % choose direction in which resp. is maximum: ks x 1
	r0 = squeeze(r(:, 1, :, 2, :)) / fs; % mean response: ks x dirs
	i = sub2ind(size(r0), (1: ks)', i); % indices of maxima
	r0 = r0(i); % keep mean amplitude where fundamental is maximum: ks x 1
	r = r1 ./ r0; % relative modulation: ks x 1
	
	% Store
	d.resp = r'; % relative modulation
	d.Properties.VariableDescriptions{'resp'} = 'Relative modulation';

function [d, m] = doMovie(d, m)
		
	% Compile and run stimulus movie
	
	% Initialise
	s = d.stim; % stimulus
	fs = size(s, 1); % number of movie frames
	x = d.x(1, :); y = x;
	colormap(m.set.colormap); % colour map

	% Compile movie
	figure(gcf); % make the figure visible
	for i = 1: fs % loop over frames
		sC = shiftdim(s(i, :, :), 1)'; % current stimulus, x is horizontal
		imagesc(x, y, sC); % draw image
		colorbar; caxis(m.set.caxis); % colour axis limits
		if i == 1, pause; % wait at first frame until a key is pressed
		else, pause(.1); % pause between frames
		end
	end
	
function [d, m] = doOrient(d, m)

% Prepare map of preferred orientation

	dirP = d.dirP; % preferred direction (deg): 1 x ks
	if m.orient.orient, dirP = mod(dirP, 180); end  % convert direction to orient.
	xs = length(m.p.x); % number of cortical x locations
	dirP = reshape(dirP, [1, xs, xs]); % convert to image: xs x ys
	d.y = d.x; d.dirP = dirP; % store

function [d, m] = doPack(d, m)

% Pack cortical locations

	% Process inputs
	x = unique(d.x); xs = length(x); % xs, number of xs
	y = unique(d.y); ys = length(y); % ys, number of ys
 	band = d.band; % bandwidth
	dirP = d.dirP; % preferred direction
	ind = d.ind; % direction selectivity index
	band = reshape(band, [xs, ys]); % reverse of unpack
	dirP = reshape(dirP, [xs, ys]);
	ind = reshape(ind, [xs, ys]);
	
	% Store
	d = d(1, :); % keep only first line
	d.x = x'; d.y = y'; % x, y
	d.band = shiftdim(band, -1); % 1 x xs x ys
	d.dirP = shiftdim(dirP, -1);
	d.ind = shiftdim(ind, -1);

function [d, m] = doPrint(d, m)

% Print the current figure to a PDF file

	set(gcf, 'paperOrientation', 'landscape');
	print(gcf, '-dpdf', 'Vis');

function [d, m] = doRefine(d, m)

% Calculate finely-spaced range of xs

	if isfield(m.refine, 'n') % set number of xs
		n = m.refine.n; % user-supplied
	else
		n = 100; % default
	end
	x = d.(m.x); % original xs
	xR = linspace(min(x), max(x), n); % refined xs
	if size(x, 2) == 1 % column vector
		xR = xR'; % keep it that way
		d = d(1, :); d = repmat(d, [n, 1]); % expand data table
	end
	d.(m.x) = xR; % store

function [d, m] = doStim(~, m)

% Calculate the stimulus

	% Calculate
	t = 0: m.p.dt: m.p.time; % times
	ts = length(t); % number of times
	x = m.p.x; y = x; % x and y values
	s = m.calStim(t, x, y, m); % stimulus: ts x xs x ys
	
	% Store
	d = table(t, x, y, 'variableNames', {'t', 'x', 'y'});
	d = repmat(d, [ts, 1]); % one row for each time
	d.t = t'; % time
	d.stim = s; % stimulus

function [d, m] = doUnpack(d, m)

% Unpack parameters in the response array: assume one-row input

	% Obtain values of parameters to be unpacked
	r = d.resp; % response: first dimension is assumed to have a single element
	if isfield(m.unpack, 'par')
		par = m.unpack.par; % names of parameters to unpack
	else
		par = m.resp.par; % default
	end
	if ~ iscell(par), par = {par}; end % make sure that it's a cell array
	pars = length(par); % number of parameters to be unpacked
	val = cell(1, pars); % values of parameters to be unpacked: allocate storage
	valE = val; % expanded values, to be stored in d
	s = zeros(1, pars); % size of parameter array: allocate storage
	for i = 1: pars % loop over parameters
		parC = par{i}; % name of current parameter
		%	val{i} = m.resp.(parC); % values of this parameter
		val{i} = d.(parC); % values of this parameter
		s(i) = length(val{i}); % number of values
	end
	[valE{:}] = ndgrid(val{:}); % expand values
	ss = prod(s); % number of values
	
	% Store
	d = repmat(d, [ss, 1]); % new number of rows equals number of values
	for i = 1: pars % loop over parameters
		parC = par{i}; % name of current parameter
		d.(parC) = valE{i}(:); % store as column vector
	end
	dim = d.Properties.CustomProperties.RespDim; % response dimensions
	[~, iU] = ismember(par, dim); % resp dimensions of parameters
	i = 1: ndims(r); i(iU) = []; % indices of remaining parameters
	i = circshift(i, -1); % remove leading singleton
	r = permute(r, [iU, i]); % move unpacked parameters to left
	s = size(r); s(1: pars) = []; s = [ss, s]; % new size of response array
	r = reshape(r, s); % reshape response
	d.resp = r; % store
	dim(iU) = []; d.Properties.CustomProperties.RespDim = dim; % update

function [d, m] = doWeight(~, m)

% Obtain geniculocortical weights for selected cortical cell

	m = m.calGain(m); % calculate geniculocortical gain
	[d, m] = m.calWeight(m); % calculate geniculocortical weight

function m = setVis
	
% setVis: set the parameters for the model of visual signal processing model

	% Set the structural parameters
	m.p.densS = [26.6, 24.4]; % ganglion cell density (cells/deg^2): off, on
	m.p.densC = 10; % (nominal) cortical cell density (cells/deg)
	m.p.dev = .166; % s.d. of same-sign ganglion cell separation: s.d/mean
	m.p.stages = [4, 5]; % number of stages: subcortical, cortical
	m.p.wid = 8; % visual field width (deg)
	
	% Set the stimulus parameters
	m.p.cont = .3; % contrast
	m.p.dir = 0; % direction relative to rightward (deg anticlockwise)
	m.p.dur = .05; % stimulus duration (s)
	m.p.freqS = .5; % spatial frequency (cycles/deg)
	m.p.freqT = 2; % temporal frequency (Hz)
	m.p.locE = [0, 0]; % sparse noise element location (deg)
	m.p.phaseS = 0; % grating spatial phase (deg)
	m.p.stim = 'grating'; % grating, reverse or sparse
	m.p.widE = 1; % sparse noise element width (deg)
	
	% Set the functional parameters
	m.p.actI = 21; % inhibitory neuron resting impulse rate (Hz)
	m.p.actS = 14; % subcortical resting impulse rate (Hz)
	m.p.cycles = floor(((0: 3) / 3) * 16000); % development cycles to store
	m.p.gainGen = 7.2; % generator function gain (Hz/mV)
  m.p.gainS = 450; % contrast sensitivity of centre mechanism (Hz/c.u.)
	m.p.gainSC0 = 3.5; % gain of subcortical-cortical convergence
	m.p.gainSE0 = 1; % gain of geniculo-excitatory convergence
	m.p.gainSIMax = 1; % gain of geniculo-inhibitory convergence
	m.p.gainIEMax = 2.25; % gain of inhibitory-excitatory convergence
	m.p.radCen = .4; % centre mechanism radius (deg)
	m.p.radGE = .95; % radius of geniculo-excitatory convergence (deg)
	m.p.radIE = .95; % radius of inhibitory-excitatory convergence (deg)
	m.p.tau = .01; % neuron time constant (s)
	m.p.tauDev = .0005; % deviation from m.p.tau (s)
	m.p.tauI = .1; % inhibitory time constant (s)
	
	% Set the computational and display parameters
	m.p.act = 1; % response type: 0 for potential, 1 for impulse rate
	m.p.calIE = 'img'; % inh. to exc. calculation: 'dot', 'imf' or 'img'
	m.p.comp = 'amp'; % response component: amp, complex or phase
	m.p.cort = 1; % stages to analyse or display: 0 for subcortex, 1 for cortex
	m.p.domain = 'freq'; % response domain: freq or time
	m.p.dt = .01; % temporal sampling interval (s)
	m.p.error = 1e-3; % maximum magnitude of imaginary components in ifft
	m.p.file = 'Vis8.mat'; % file containing weights
	m.p.fs = 64; % number of temporal frequencies: make it even
	m.p.gpu = 0; % 0 for no GPU, 1 for GPU
	m.p.implicit = 1; % calculate stimulus implicitly in response simulation
	m.p.parallel = 0; % 0 for single-stream processing, 1 for parallel
	m.p.solver = 'solveF'; % function to solve model: solveF or solveT
	m.p.stage = 2; % cortical processing stage to display, excitatory by default
	m.p.time = 1.5; % simulation time (s)
	
	% Set special cases
	if 1 % prevent trailing comments when code folding
	end
