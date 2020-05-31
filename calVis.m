function m = calVis(m)

% calVis: calculate structure, stimuli and responses for the visual signal
% processing model

	% Define function handles for functions in this file
	m.calChan = @calChan;
	m.calComp = @calComp;
	m.calField = @calField;
	m.calGain = @calGain;
	m.calInd = @calInd;
	m.calResp = @calResp;
	m.calStim = @calStim;
	m.calTemp = @calTemp;
	m.calWeight = @calWeight;
	m.cor.fun = @calCor;
	m.corSpat.fun = @calCorSpat;
	m.corSum.fun = @calCorSum;
	m.dif.fun = @calDif;
	m.inhib.fun = @calInhib;
	m.modelTuning = @modelTuning;
	m.offset.fun = @calOffset;
	m.phase.fun = @calPhase;
	m.selDir.fun = @selDir;
	m.solveF = @solveF;
	m.solveT = @solveT;

function m = calChan(m)
	
% Calculate the location of the subcortical channels and cortical locations

	% Calculate subcortical channels
	rng('default'); % make the results reproducible
	loc = cell(2, 1); % allocate storage
	for i = 1: 2 % off- then on-cells
		
		% Calculate locations, assuming square array
		sep = 1 / sqrt(m.p.densS(i)); % mean separation between cells
		if i == 1 % off-cells centred on (0, 0)
			xP = sep: sep: .5 * m.p.wid; % cells with positive x
			x = [- fliplr(xP), 0, xP]; % all cells
		else % on-cells offset
			xP = .5 * sep: sep: .5 * m.p.wid; % offset from centre
			x = [- fliplr(xP), xP]; % no central cell
		end
		[x, y] = ndgrid(x); % two-dimensional locations

		% Perturb the array
		dev = m.p.dev * sep; % standard deviation of x location
		x = x + normrnd(0, dev, size(x)); % perturb with Gaussian deviates
		y = y + normrnd(0, dev, size(y));

		% Store
		s = repmat(2 * (i - 1) - 1, [numel(x), 1]); % channel sign
		loc{i} = [x(:), y(:), s]; % one row for each channel
		
	end
	m.p.chan = vertcat(loc{:}); % concatenate off- and on-locations

	% Calculate cortical locations
	n = floor(.5 * m.p.wid * m.p.densC); % number of cells with positive x
	x = n / m.p.densC; % maximum x (deg)
	x = linspace(- x, x, 2 * n + 1); % x locations (deg)
	m.p.x = x; % store
	[x, y] = ndgrid(x); % two-dimensional locations (deg)
	m.p.loc = [x(:), y(:)]; % cortical locations (deg): ks x 2
	
function [d, m] = calComp(d, m)

% Compile responses to contrast-reversing grating into resp. to moving grating

	% Initialise
	r = d.resp; % response (mV): 1 x fs x ks x zs x ds x ps
	r0 = r(1, 1, :, :, :, 1); % time-invariant component: 1 x 1 x ks x zs x ds
	f = m.p.f; % frequencies (Hz): 1 x fs
	fs = length(f); % number of frequencies
	f0 = f(2); % fundamental frequency
	dir = d.dir; ds = length(dir); % directions (deg), number of directions
	p = m.resp.phaseS; % spatial phase (deg): 1 x ps
	p = pi * p / 180; % convert to radians
	ps = length(p); % number of phases
	
	% Time-shift the contrast-reversing responses and sum them
	s = (f' / f0); % frequency as a multiple of fundamental: fs x 1
	s = s .* sign(dir); % include direction: fs x ds
	s = s .* reshape (p, [1, 1, ps]); % phase shift (radians): fs x ds x ps
	s = exp(1i * s); % phase shift multiplier
	s = reshape(s, [1, fs, 1, 1, ds, ps]); % 1 x fs x 1 x 1 x ds x ps
	r = s .* r; % time-shift the response
	r = (2 / ps) * sum(r, 6); % sum over temporal phases: 1 x fs x ks x zs x ds
	r(1, 1, :, :, :) = r0; % restore time-invariant component
	
	% Store
	d.resp = ifftReal(r); % convert frequency response to temporal response
	d.time = m.p.t; % time (s)

function [d, m] = calCor(d, m)

% Correlate summed excitation at reference cell with each inhibitory input

	% Calculate correlation
	r = d.resp; % generator potential time course (mV): 1 x ts x ks x zs
	[~, ~, ks, ~] = size(r); % number of elements along each dimension
	rE = r(:, :, :, 5); % summed excitation: 1 x ts x ks
	rE = shiftdim(rE, 1); % shift time to the first dimension: ts x ks
	rI = r(:, :, :, 3); % inhibition from individual cells: 1 x ts x ks
	rI = shiftdim(rI, 1); % shift time to first dimension: ts x ks
	r = [rE, rI]; % prepare for correlation: ts x 2 * ks
	r = cov(r); % covariance: 2 * ks x 2 * ks
	v = diag(r); v = v(1: ks); % variance of excitatory cells: ks x 1
	r = r(1: ks, ks + 1: end); % cov. of excitation and inhibition: ks x ks
		% each row represents one excitatory cell
	r = r ./ v; % normalise covariance by variance: ks x ks
	
	% Weight correlation by distance from excitatory cell
	g = m.p.loc; % cortical cell locations (deg): ks x 2
	g = reshape(g, [1, ks, 2]) - reshape(g, [ks, 1, 2]);
		% displacement of neighbour location from reference location: ks x ks x 2
	g = exp(- dot(g, g, 3) / m.p.radIE ^ 2); % weighting function: ks x ks
	r = g .* r; % weight orientation vector by separation: ks x ks
	
	%	Store
	d.y = d.x; % vertical location (deg)
	xs = length(m.p.x); % number of cells in horizontal direction
	d.resp = reshape(r, [1, ks, xs, xs]); % exc./inh. correl.: 1 x ks x xs x ys
	d.Properties.CustomProperties.RespDim = {'dir', 'loc', 'x', 'y'}; % resp. dim.

function [d, m] = calCorSum(d, m)

% Correlate summed excitation with summed inhibition

	% Calculate correlation
	r = d.resp; % generator potential time course (mV): 1 x ts x ks x zs
	[~, ~, ks, ~] = size(r); % number of elements along each dimension
	rE = r(:, :, :, 5); % summed excitation: 1 x ts x ks
	rE = shiftdim(rE, 1); % shift time to the first dimension: ts x ks
	rI = r(:, :, :, 4); % summed inhibition: 1 x ts x ks
	rI = shiftdim(rI, 1); % shift time to first dimension: ts x ks
	r = [rE, rI]; % prepare for correlation: ts x 2 * ks
	r = corrcoef(r); % correlation coefficient: 2 * ks x 2 * ks
	r = diag(r, ks); % corr. of excitation and inhibition: ks x 1
	
	% Store
	d.resp = r'; % correlations: 1 x ks
	d.Properties.CustomProperties.RespDim = {'dir', 'loc'}; % response dimensions
	
function [d, m] = calCorSpat(d, m)

% Correlate receptive fields

	% Calculate correlation
	r = d.resp; % potential time course (mV): cs x xs x ys x ks x zs
	[~, xs, ys, ks, zs] = size(r); % number of stimulus locations, neurons, stages
	r = m.p.gainGen * max(0, r); % convert to impulse rate (Hz)
	r = r(2, :, :, :, :) - r(1, :, :, :, :); % on - off: 1 x xs x ys x ks x zs
	r = reshape(r, [xs * ys, ks, zs]); % prepare for correl'n: xs * ys x ks x zs
	switch m.p.stage % excitatory or inhibitory?
		case 2 % correlate excitatory receptive fields
			r = r(:, :, 2); % select required stage: xs * ys x ks
			r = cov(r); % covariance of reference and surrounding r.f.s: ks x ks
			v = diag(r); % variance of reference r.f.s: ks x 1
			r = r ./ v; % normalise: ks x ks
		case 3 % correlate excitatory r.f. with r.f. of surrounding inhibitory cells
			rE = r(:, :, 2); % reference cells are excitatory: xs * ys x ks
			rI = r(:, :, 3); % surrounding cells are inhibitory: xs * ys x ks
			r = [rE, rI]; % prepare for correlation: xs * ys x 2 * ks
			r = cov(r); % covariance of reference and surround r.r.s: 2 * ks x 2 * ks
			r = r(1: ks, ks + 1: end); % cov. of excitation and inhibition: ks x ks
				% each row represents one excitatory cell
			v = diag(r); v = v(1: ks); % variance of excitatory cells: ks x 1
			r = r ./ v; % normalise covariance by variance: ks x ks
	end
	
	%	Store
	d = d(1, :); % keep only the first row
	xs = length(m.p.x); % number of cortical neurons in x-direction
	d.resp = reshape(r, [1, ks, xs, xs]); % turn into image: 1 x ks x xs x ys
	d.Properties.CustomProperties.RespDim = {'', 'loc', 'x', 'y'};

function [d, m] = calDif(d, m)

% Take difference between antipreferred and preferred exc./inh. correlations

	% Calculate difference and its mean
	r = d.resp; % correlations: 2 x ks x xs x xs
	[~, ks, xs, ~] = size(r); % number of elements along each dimension
	r = r(2, :, :, :) - r(1, :, :, :); % antipref. - preferred: 1 x ks x xs x xs
	rM = reshape(r, [1, xs, xs, ks]); % prepare for mean: 1 x xs x xs x ks
	rM = mean(rM, 4); % mean over all inhibitory inputs: 1 x xs x xs

	% Store
	d = d(1, :); % keep only first line
	d.resp = r; % difference of correlations: 1 x ks x xs x xs
	if m.dif.mean % response is mean difference
		d.resp = rM / .009; % mean diff. norm. by 95th percentile: 1 x xs x xs
	end

function dp = calDeriv(t, p, m)

%	Calculate the derivative of the generator potential
%
%	Internal variables:
%		q(c, z) = generator potential for subcortical channel c and stage z
%		dq(c, z) = time derivative of q
%		r(k, z) = generator potential for cortical location k, stage z
%		dr(k, z) = time derivative of r

	% Set the computed parameters
	dir = 2 * pi * m.p.dir / 360; % motion direction (radians from rightward)
	fS = 2 * pi * m.p.freqS; % spatial frequency (radians/deg)
	fT = 2 * pi * m.p.freqT; % temporal frequency (radians/s)
	tauS = m.p.tau + m.p.chan(:, 3) * m.p.tauDev; % t/c depends on ch. sign
	
	% Sort input into channels and stages
	cs = size(m.p.chan, 1); % number of channels
	zs = m.p.stages(1); % number of subcortical stages
	ns = cs * zs; % number of subcortical neurons
	q = reshape(p(1: ns), [cs, zs]); % cs x zs
	dq = zeros(cs, zs); % allocate storage
	ks = size(m.p.loc, 1); % number of cortical locations
	xs = length(m.p.x); % number of cortical locations along a single dimension
	zs = m.p.stages(2); % number of cortical stages
	r = reshape(p(ns + 1: end), [ks, zs]); % ks x zs
	dr = zeros(ks, zs); % allocate storage

	% Calculate subcortical input: dot product of stimulus and receptive field
	switch m.p.stim % stimulus type
		case 'grating' % drifting grating
			if m.p.implicit % stimulus is calculated implicitly
				inp = m.p.cont * exp(-.25 * (m.p.radCen * fS) ^ 2); % dot product
				u = m.p.chan(:, 1: 2) * [cos(dir); sin(dir)];
					% channel locations (deg from centre in direction of motion)
				inp = inp * cos(fS * u - fT * t); % attenuation due to channel location
			else % stimulus is calculated explicity
				inp = calInp(t, m); % doesn't work: don't know why
			end
		case 'sparse' % sparse noise
			inp = calInp(t, m)'; % dot product
	end
	inp = - (m.p.gainS / m.p.gainGen) * inp; % convert from contrast units to mV

	% Calculate subcortical potentials
	dq(:, 1) = (1 / m.p.tau) * (inp - q(:, 1));
	dq(:, 2) = (1 ./ tauS) .* (- m.p.chan(:, 3) .* q(:, 1) - q(:, 2));
	dq(:, 3) = (1 ./ tauS) .* (q(:, 2) + m.p.actS / m.p.gainGen - q(:, 3));
	dq(:, 4) = (1 ./ tauS) .* (q(:, 3) - q(:, 4));

	% Use GPU
	if m.p.gpu % use GPU
		g = gpuArray(m.p.gainGE); w = gpuArray(m.p.weight);
	else % no GPU
		g = m.p.gainGE; w = m.p.weight;
	end
	
	% Input from subcortex
	inpS = q(:, 4); % input from subcortex: cs x 1
	inpS = sum(g .* w .* inpS)'; % sum over subcortical channels
	
	% Retreive from GPU
	if m.p.gpu, inpS = gather(inpS); end
	
	% Inhibitory dendrites and soma
	inp = m.p.gainSI0 * inpS; % include gain: ks x 1
	dr(:, 1) = (1 ./ m.p.tau) .* (inp - r(:, 1));
	
	% Inhibitory output
	inpI = max(0, r(:, 1)); % rectified input from inhib. cells
	dr(:, 3) = (1 ./ m.p.tauI) .* (inpI - r(:, 3)); % filter
	
	% Inhibitory input to excitatory stage
	inpI = reshape(r(:, 3), [xs, xs]); % xs x xs
	switch 'im' % convolution or imgaussfilt
		case 'conv' % convolution
			inpI = conv2(m.p.gainIE, inpI, 'same'); % convolution with Gaussian
		case 'im' % filter in spatial or frequency domain
			dev = m.p.radIE / sqrt(2); % standard deviation of Gaussian (deg)
			dev = m.p.densC * dev; % standard deviation of Gaussian (samples)
			inpI = imgaussfilt(inpI, dev); % filter
	end
	inpI = inpI(:); % ks x 1

	% Excitatory cortical stage
	inp = m.p.gainSE0 * inpS; % input from subcortex
	inpI = m.p.gainIE0 * inpI; % inhibitory input
	dr(:, 2) = (1 ./ m.p.tau) .* (inp - inpI - r(:, 2));

	% Vectorise the derivatives
	dp = [dq(:); dr(:)];

function [d, m] = calField(d, m)

% Calculate receptive field

	% Initialise
	m.p.weight = d.weight{1}; % current weights
	m.p.gainSI0 = d.gainSI0; m.p.gainIE0 = d.gainIE0; % current gains
	cont = [-1; 1]; % dark, light
	wid = m.p.wid; % width of area in which to present elements (deg)
	locs = 4 * wid + 1; % number of element locations in each direction
	x = linspace(-.5 * wid, .5 * wid, locs); % element x-locations (deg)
	y = x; % element y-locations
	t = m.p.t; % times (s)
	ks = size(m.p.loc, 1); % number of cortical locations
	zs = m.p.stages(2); % number of cortical stages
	r = zeros(2, locs, locs, ks, zs); % allocate storage: cs x xs x ys x ks x zs
	
	% Calculate response
	for i = 1: 2 % dark then light element
		m.p.cont = cont(i); % contrast
		for j = 1: locs % x locations
			for k = 1: locs % y locations
				m.p.locE = [x(j), y(k)]; % element location
				switch m.p.solver % calculate cortical potential
					case 'solveF'
						[~, ~, p] = solveF(m); % Fourier transform: fs x ks x zs
						p = ifftReal(p, [], 1, inf); % inverse transform: ts x ks x zs
					case 'solveT', [~, ~, p] = solveT(m); % time course: ts x ks x zs
				end
				r(i, j, k, :, :) = max(p); % maximum response: cs x xs x ys x ks x zs
				p = p(t < .1, :, :); % remove inhibitory rebound: time < .1 s
				r(i, j, k, :, [1, 2, 5]) = max(p(:, :, [1, 2, 5])); % fast stages
			end
		end
	end
	
	% Store
	x = repmat(x, [2, 1]); y = x; % one for each contrast
	dirP = repmat(d.dirP, [2, 1]); % preferred direction (deg): 1 x ks
	d = table(cont, x, y, r, dirP, ...
		'variableNames', {'cont', 'xStim', 'yStim', 'resp', 'dirP'});
	d.Properties.VariableDescriptions{'xStim'} = 'Stimulus x-location (deg)';
	d.Properties.VariableDescriptions{'yStim'} = 'Stimulus y-location (deg)';
	d.Properties.VariableDescriptions{'resp'} = 'Generator potential (mV)';
	d = addprop(d, 'RespDim', 'table'); % add property: response dimensions
	d.Properties.CustomProperties.RespDim = ...
		{'cont', 'xStim', 'yStim', 'loc', 'stage'};

function m = calGain(m)

% Calculate gains

	% Initialise
	cs = size(m.p.chan, 1); % number of channels
	ks = size(m.p.loc, 1); % number of cortical locations

	% Calculate location array: cs x ks x 2
	locS = m.p.chan(:, 1: 2); % subcortical locations: cs x 2
	locS = reshape(locS, [cs, 1, 2]); % create index for cortex: cs x 1 x 2
	locC = m.p.loc; % cortical locations: ks x 2
	locC = reshape(locC, [1, ks, 2]); % create index for subcortex: 1 x ks x 2
	
	% Calculate geniculocortical gain
	dSq = sum((locC - locS) .^ 2, 3); % squares of geniculocortical distances
	g = exp(- dSq / m.p.radGE ^ 2); % weight: cs x ks
	g = g ./ sum(g, 1); % normalise to unity-sum
	m.p.gainGE = m.p.gainSC0 * g; % store
	
	% Calculate inhibitory to excitatory gain
	switch m.p.calIE % dot product?
		case 'dot' % yes, using inhibitory to excitatory weights
			locI = m.p.loc; locC = locI; % cortical locations (deg): ks x 2
			locI = reshape(locI, [ks, 1, 2]); % create inh. neuron index: ks x 1 x 2
			locC = reshape(locC, [1, ks, 2]); % create exc. neuron index: 1 x ks x 2
			dSq = sum((locC - locI) .^ 2, 3); % squares of intracort. dist.: ks x ks
			g = exp(- dSq / m.p.radIE ^ 2); % weight: ks x ks
			g = g ./ sum(g, 1); % normalise to unity-sum
			m.p.gainIE = g; % store
		otherwise % no, not using inhibitory to excitatory weights
			x = m.p.x; % cortical locations
			g = exp(- (x / m.p.radIE) .^ 2); % Gaussian profile: 1 x xs
			g = g' * g; % gain: xs x xs
			g = g / sum(g(:)); % normalise to unity-sum
			m.p.gainIE = g; % store
	end

function [d, m] = calInd(d, m)

% Calculate direction selectivity index from receptive field correlations

	% Obtain receptive field correlations, probabilities
	r = d.resp; % receptive field correlations: 1 x ks x xs x ys
	[~, ks, xs, ys] = size(r); % ks is the number of excitatory cortical cells
	r = reshape(r, [ks, ks]); % correlations between their receptive fields
	p = reshape(d.prob, [ks, ks]); % probability of null hypothesis: ks x ks
	
	% Find locations in preferred half-field and significant locations
	dis = m.p.loc; % cortical cell locations (deg): ks x 2
	dis = reshape(dis, [1, ks, 2]) - reshape(dis, [ks, 1, 2]);
		% displacement of neighbour location from reference location: ks x ks x 2
	dir = pi * d.dirP' / 180; % preferred direction (radians): ks x 1
	dir = [cos(dir), sin(dir)]; % direction as vector: ks x 2
	dir = reshape(dir, [ks, 1, 2]); % reshape for dot product: ks x 1 x 2
	i = sum(dis .* dir, 3); % projection of displacement onto direction: ks x ks
	i = i <= 0; % 1 for locations in preferred half-field: ks x ks
	j = p <= 0.05; % 1 for locations in which correlation is significant: ks x ks
	
	% Calculate index and store
	rP = r; rP(~ i | ~ j) = nan; rP = mean(rP, 2, 'omitnan');
		% mean of significant correlations in preferred half-field: ks x 1
	rN = r; rN(i | ~ j) = nan; rN = mean(rN, 2, 'omitnan'); % nonpreferred: ks x 1
	i = rP - rN; % direction selectivity index: ks x 1
	d.resp = reshape(i, [1, xs, ys]); % make it an image: 1 x xs x ys
	d.x = m.p.x; d.y = d.x; % store locations

function [d, m] = calInhib(d, m)

% Calculate response for static inhibition in one half-field
%
% Input: response with inhibitory time constant as last dimension:
%   first element is for dynamic inhibition, second for static.
% Output: first row contains response for dynamic inhibition in
%   nonpreferred half-field, second row in preferred.

	% Obtain grating response and select components of interest
	r = d.resp; % response (mV): 1 x fs x ks x zs x ss x is
	[~, fs, ks, zs, ss, ~] = size(r); % size of response array
	rD = r(:, :, :, :, :, 1); % resp. with dynamic inhib.: 1 x fs x ks x zs x ss
	rS = r(:, :, :, :, :, 2); % resp. with static inhib.: 1 x fs x ks x zs x ss
	
	% Change inhibition to static in half-field
	dis = m.p.loc; % cortical locations (deg): ks x 2
	dis = dis - m.loc.loc; % distance from selected location
	k = knnsearch(m.p.loc(:, 1: 2), m.loc.loc); % location closest to m.loc.loc
	dir = pi * d.dirP(k) / 180; % convert preferred direction to radians
	dir = [cos(dir); sin(dir)]; % direction as vector: 2 x 1
	i = dis * dir < 0; % truth of locations in preferred half-field: ks x 1
	rN = rD; % response for dynamic inhibition in nonpreferred half-field
	rN(:, :, i, :, :) = rS(:, :, i, :, :); % change preferred half-field
	rP = rD; % response for dynamic inhibition in preferred half-field
	rP(:, :, ~ i, :, :) = rS(:, :, ~ i, :, :); % change nonpreferred half-field
	rS = [rN; rP]; % store: 2 x fs x ks x zs x ss

	% Weight response by distance from selected location
	g = exp(- dot(dis, dis, 2) / m.p.radIE ^ 2); % weighting function: ks x 1
	g = g / sum(g); % normalise to unity-sum
	g = reshape(g, [1, 1, ks]); % prepare for weighting: 1 x 1 x ks 
	rS = g .* rS; % weight response: 2 x fs x ks x zs x ss
	rS = sum(rS, 3); % sum responses in half-field: 2 x fs x 1 x zs x ss
	rS = d.gainIE0(1) * rS; % include inhibitory-to-excitatory gain
	
	% Recalculate response in excitatory neuron
	inpS = rD(:, :, k, 5, :); % summed subcortical input: 1 x fs x 1 x 1 x ss
	inpI = rS(:, :, :, 3, :); % summed inhibitory input: 2 x fs x 1 x 1 x ss
	f = 2 * pi * m.p.f; % frequency (radians/s)
	r = r(:, :, k, :, :, 1); % original response: 1 x fs x 1 x zs x ss
	r = [r; r]; % make it 2 rows
	r(:, :, :, 2, :) = (inpS - inpI) ./ (1 + 1i * m.p.tau * f);
		% excitatory response: 2 x fs x 1 x 1 x ss
	r(:, :, :, 4, :) = inpI; % summed inhibition: 2 x fs x 1 x 1 x ss
	r = ifftReal(r, [], 2); % convert to temporal response
	
	% Store
	%	d = d(1, :); % keep only first row
	xs = length(m.p.x); % prepare to reshape i
	d.field = reshape(i, [1, xs, xs]);
		% 1 if location is in preferred half-field, 0 otherwise: 1 x xs x xs
	d.Properties.VariableDescriptions{'field'} = 'Preferred half-field';
	d.time = m.p.t; % times (s)
	d = repmat(d, [2, 1]); % 2 rows
	d.side = categorical(["non"; "pref"]); % side of dynamic inhibition
	d.resp = reshape(r, [2, fs, zs, ss]); % response (mV): 2 x fs x zs x ss
	dim = d.Properties.CustomProperties.RespDim; % response dimensions
	[~, i] = ismember('loc', dim); % response dimensions of parameter
	dim(i) = []; d.Properties.CustomProperties.RespDim = dim; % update
	
function inp = calInp(t, m)

% Calculate subcortical input: dot product of stimulus and convergence function
%	Input: one or more stimulus times (s)
% Output: subcortical input at time(s) t for each channel: ts x js

	switch m.p.stim % stimulus type
		
		case 'grating' % drifting grating
			
			x = m.p.x; y = x; % x and y values
			s = squeeze(calStim(t, x, y, m)); % stimulus
			g = exp(- (x / m.p.radCen) .^ 2); % Gaussian profile: 1 x xs
			g = g' * g; % gain: xs x xs
			g = g / sum(g(:)); % normalise to unity-sum
			inp = conv2(s, g, 'same'); % convolution: xs x xs
			inp = interp2(x, y, inp, m.p.chan(:, 1), m.p.chan(:, 2));
				% interpolate at channel locations: js x 1

		case 'sparse' % sparse noise: each element is a square of light or dark

			% Calculate stimulus:
			%		m.p.locE = location of element centre (deg): 1 x 2
			%		m.p.widE = element width (deg): 1 x 1
			sig = m.p.radCen / sqrt(2); % s.d. of convergence function (deg)
			locL = (m.p.locE - .5 * m.p.widE - m.p.chan(:, 1: 2)) / sig;
				% location of element lower edge (standard deviations): js x 2
			locU = (m.p.locE + .5 * m.p.widE - m.p.chan(:, 1: 2)) / sig;
				% location of element upper edge (standard deviations): js x 2
			locL = locL'; locU = locU'; % 2 x js

			% Calculate subcortical input
			t = t(:); % column vector
			i = t <= m.p.dur; % indices of time(s) when element is on: ts x 1
			inp = m.p.cont * i .* (normcdf(locU(1, :)) - normcdf(locL(1, :))) .* ...
				(normcdf(locU(2, :)) - normcdf(locL(2, :))); % dot product
			
	end

function [d, m] = calOffset(d, m)

% Calculate spatial offset between exc. and inh. phase sensitivity functions

	% Calculate
	r = d.resp; % response phase: 1 x ks x zs x ps
	[~, ks, ~, ~] = size(r); % number of cells
	rIn = r(:, :, [4, 5], :); % keep only summed responses: 1 x ks x zs x ps
	p = pi * d.phaseS' / 180; % spatial phase (radians): ps x 1
	u = [cos(p), sin(p)]; % unit vectors, cosine and sine phase: ps x 2
	r = zeros(ks, 2, 2); % allocate storage: ks x zs x 2
	for j = 1: 2 % inhibitory then excitatory
		for i = 1: ks % loop over cells
			rC = shiftdim(rIn(1, i, j, :), 3); % ps x 1
			r(i, j, :) = u \ rC; % regression of response on unit vector
		end
	end
	r = atan(real(r(:, :, 2) ./ r(:, :, 1))); % spatial phase (radians): ks x 2
	r = 180 * r / pi; % spatial phase (deg): ks x 2
	r = r(:, 1) - r(:, 2); % inh. phase - exc. phase (deg): ks x 1
	r = mod(r + 10, 180) - 10; % range is [-10, 170) (deg): ks x 1
	
	%	Store
	xs = length(m.p.x); % number of cortical locations along a single dimension
	r = reshape(r, [1, xs, xs]); % offset (deg): 1 x xs x ys
	if isfield(m.offset, 'crop') % crop the image
		xsC = m.offset.crop; % cropped size
		x = round(.5 * (xs - xsC)); % number of x values to crop at each edge
		r = r(1, 1 + x: xs - x, 1 + x: xs - x); % cropped image: 1 x xs x xs
	end
	%	[p, h, s] = signtest(r(:), 0, 'tail', 'right')
	d.resp = r; % offset (deg): 1 x xs x ys

function [d, m] = calPhase(d, m)

% Prepare to plot phase sensitivity function

	rIn = d.resp; % response: 1 x zs x ps
	r = rIn(:, m.phase.stage, :); r = shiftdim(r, 1); % select stage: 1 x ps
	switch m.phase.comp % response component
		case 'amp'
			r = abs(r); % response amplitude
		case 'complex' % calculate response and simulated response
			rIn = rIn(1, [4, 5], :); % keep only summed responses: 1 x 2 x ps
			p = pi * d.phaseS' / 180; % spatial phase (radians): ps x 1
			u = [cos(p), sin(p)]; % unit vectors, cosine and sine phase: ps x 2
			rS = zeros(2, 2); % simulated response: stages x components, zs x cs
			for i = 1: 2 % inhibitory then excitatory
				rC = shiftdim(rIn(1, i, :), 2); % ps x 1
				rS(i, :) = u \ rC; % regression of response on unit vector
			end
			u = sum(rS, 2); % inh. and exc. vectors: zs x 1
			u = u ./ abs(u); % make them unit vectors: zs x 1
			rM = sqrt(sum(abs(rS) .^ 2, 2)); % response magnitude (mV): zs x 1
			u = rM .* u; % asssign magnitude: zs x 1
			d.respComp = u.'; % excitatory and inhibitory components: 1 x zs
			rA = atan(real(rS(:, 2) ./ rS(:, 1))); % response angle (radians): zs x 1
			rS = u .* cos(p' - rA); % modulate by spatial phase: zs x ps
			rS = diff(rS); % simulated response is exc. - inh.: 1 x ps
			rS = rS / (1 + 1i * m.p.tau * 2 * pi * m.p.freqT); % first-order delay
			rS(end + 1) = rS(1); % complete the ellipse
			d.respSim = rS; % store
		case 'phase'
			r = angle(r); % phase (radians): 1 x ps
			r = unwrap(r); % unwrapped phase (radians): 1 x ps
			if m.phase.shift % shift laterally
				i = diff(r); % differences between successive phases
				[~, i] = max(abs(i)); % index of maximum slope magnitude
				r = circshift(r, - i + 1); % shift so that the maximum slope is first
			end
			r = r - r(1); % set the origin to (0, 0)
			r = unwrap(r); % unwrap again
			r = 180 * r / pi; % phase (deg): 1 x ps
			d.Properties.VariableDescriptions{'resp'} = 'Response phase (deg)';
	end
	d.resp = r;
	
function [d, m] = calResp(d, m)

% Calculate cortical frequency response for specified model parameter values
%
%	Inputs:
%   m.resp.par = names of parameters with specified values: cell
%   m.resp.pari = values for parameter pari: double
% Output:
%   d.resp = response(freq, cell, stage, par1, par2, ...): double (complex)

	% Initialise model parameters
	m.p.weight = d.weight{1}; % current weights
	m.p.gainSI0 = d.gainSI0; m.p.gainIE0 = d.gainIE0; % current gains
	
	% Initialise response array
	%	fs = m.p.fs; % number of frequencies
	fs = length(m.p.f); % number of frequencies
	ks = size(m.p.loc, 1); % number of cortical cells
	zs = m.p.stages(2); % number of cortical stages
	vs = 1; % default number of parameter values
	if isfield(m.resp, 'par') % there are multiple parameter values
		par = m.resp.par; % parameter names
		if ~ iscell(par), par = {par}; end % make it cell if it isn't already
		pars = length(par); % number of parameters
		s = zeros(1, pars); % number of values for each parameter
		for i = 1: pars % loop over parameters
			s(i) = length(m.resp.(par{i})); % number of values for this parameter
			d.(par{i}) = m.resp.(par{i}); % store values
		end
		vs = prod(s); % number of values
		sub = cell(1, pars); % for subscripts of a value
	end
	resp = zeros(fs, ks, zs, vs); % allocate storage: fs x ks x zs x vs
	
	% Calculate response
	for i = 1: vs % loop over parameter values
		if exist('par', 'var') % parameter values are specified
			[sub{:}] = ind2sub(s, i); % subscripts for current value
			for j = 1: pars % loop over parameters
				v = m.resp.(par{j}); % values for this parameter
				m.p.(par{j}) = v(sub{j}); % set value for this parameter
			end
		end % set current stimulus value
		
		% Calculate response
		switch m.p.solver % frequency or temporal domain?
			case 'solveF' % frequency domain
				[~, pS, pC] = solveF(m); % fft of generator potential (mV)
				dom = "freq"; % domain: frequency response
			case 'solveT' % temporal domain
				[~, pS, pC] = solveT(m); % time course of generator potential (mV)
				dom = "time"; % temporal response
		end
		switch m.p.cort % subcortex or cortex?
			case 0, r = pS; % subcortex: fs x js x zs
			case 1, r = pC; % cortex: fs x ks x zs
		end
		if m.p.act % convert to impulse rate
			if dom == "freq" % this is a frequency response
				r = ifftReal(r); % convert to temporal response
				dom = "time";
			end
			r = m.p.gainGen * max(0, r); % impulse rate (Hz)
		end
		switch m.p.domain % output frequency or temporal response?
			case 'freq' % frequency response
				if dom == "time", r = fft(r); end % convert temporal response
			case 'time' % temporal response
				if dom == "freq", r = ifftReal(r); end % convert frequency resp.
		end
		resp(:, :, :, i) = r; % store the response for this set of parameter values
				
	end
	
	% Store
	d.weight = []; % don't need it anymore
	d.resp = reshape(resp, [1, fs, ks, zs, s]);
		% 1 x fs x ks x zs x par1s x par2s x ...
	d.Properties.CustomProperties.RespDim = ...
		[{'', m.p.domain, 'loc', 'stage'}, par]; % response dimensions
	d.stage = 1: m.p.stages(2); % store stage numbers
	switch m.p.domain % store domain values
		case 'freq', d.freq = m.p.f;
		case 'time', d.time = m.p.t;
	end

function s = calStim(t, x, y, m)

% Calculate the stimulus
% Inputs:
%		t = times (s)
%		x = x locations (deg)
%		y = y locations (deg)
%		m = metadata
%	Ouput:
%		s(t, x, y) (contrast-units): ts x xs x ys

	% Initialise
	xs = length(x); % number of visual field locations
	s = zeros(xs); % zero stimulus
	
	% Calculate grating
	[t, x, y] = ndgrid(t, x, y); % grid times and locations
	phaseT = 2 * pi * m.p.freqT * t; % temporal phase (radians)
	dir = 2 * pi * m.p.dir / 360; % direction (radians)
	u = cos(dir) * x + sin(dir) * y; % location in stimulus direction
	freqS = 2 * pi * m.p.freqS; % spat. freq. (radians / deg)
	phaseS = (pi / 180) * m.p.phaseS; % spatial phase (radians)
	switch m.p.stim
		case 'grating' % drifting grating
			s = cos(u * freqS - phaseT);
		case 'flash' % flashed stationary grating
			s = cos((cos(dir) * x + sin(dir) * y) * freqS + phaseS);
		case 'reverse' % contrast-reversing stationary grating
			%	s = cos(phaseT) .* cos(u * freqS + phaseS);
			s = cos(phaseT) .* cos(u * freqS - phaseS);
	end
	s = m.p.cont * s; % set contrast

function m = calTemp(m)

% Calculate temporal parameters: frequencies and times

	% Calculate fundamental frequency and times
	switch m.p.stim % stimulus type
		case {'grating', 'reverse'} % cyclic
			fs = m.p.fs; % number of frequencies: assume that it's even
			f = m.p.freqT; % fundamental frequency (Hz)
			t = 1 / f; % duration (s)
			t = linspace(0, t, fs + 1); % times (s)
			t = t(1: end - 1); % open-ended interval
		case 'sparse' % transient
			t = 0: m.p.dt: m.p.time; % simulation times
			t = t(1: end - 1); % open-ended interval: 1 x ts
			fs = length(t); % number of frequencies: assume that it's even
			f = 1 / m.p.time; % fundamental frequency (Hz)
	end
	
	% Calculate frequencies
	i = .5 * fs; % index of highest positive frequency
	i = fftshift(linspace(- i, i - 1, fs)); % indices of fft frequencies
	f = i * f; % frequencies (Hz)
	
	% Store
	m.p.f = f; % frequencies (Hz)
	m.p.t = t; % times (s)

function [gainSI0, gainIE0] = calVarGain(m, cycle)

% Calculate gains that change during development

	% Increase gains linearly with development cycle
	cycles = max(m.p.cycles); % number of development cycles
	gainSI0 = 1 + (cycle / cycles) * (m.p.gainSIMax - 1); % subcortex to inh.
	gainIE0 = 1 + (cycle / cycles) * (m.p.gainIEMax - 1); % inhibitory to exc.

function [d, m] = calWeight(m)

% Calculate geniculocortical synaptic weights with Hebbian process

	% Calculate weights
	if exist(m.p.file, 'file') % retrieve weights from file
		load(m.p.file, 'd'); % d = s.d;
		m.p.cycles = d.cycle'; % development cycles
		[d.gainSI0, d.gainIE0] = calVarGain(m, d.cycle);
			% store variable gains in file
	else
		
		% Initialise for development
		fs = m.p.fs; % number of frequencies
		cs = size(m.p.chan, 1); % number of subcortical channels
		ks = size(m.p.loc, 1); % number of cortical locations
		cycles = max(m.p.cycles); % number of development cycles
		w = ones(cs, ks); % synaptic weights, 0 - 2: cs x ks
		wInc = .2; % weight increment on a cycle
		dirs = 16; % number of directions (deg)
		dir = linspace(0, 360, dirs + 1); % directions (deg)
		dir = dir(1: dirs); % open-ended interval
		dirOrig = m.p.dir; % original setting

		% Calculate weights
		%	rng('default'); % test whether changing channel choice changes results
		timeS = clock; timeR = timeS; % starting time, report time
		dC = cell(cycles, 1); % allocate storage
		for cycle = 0: cycles % loop over development cycles

			% Vary weights and gains
			if cycle > 0 % vary weights
				c = randi(cs); % choose a random channel
				w(c, :) = min(2, w(c, :) + wInc); % increase the weight of that channel
			end
			m.p.weight = w; % current weights
			[m.p.gainSI0, m.p.gainIE0] = calVarGain(m, cycle);
			
			% Calculate cortical response as a maximum over orientations
			q = zeros(2, ks, dirs); % allocate storage
			if m.p.parallel % parallel processing
				parfor i = 1: dirs % loop over directions, using parallel processors
					[~, ~, p] = solveF(m, dir(i)); % solve model for cortical potential
					q(:, :, i) = p(1: 2, :, 2); % DC and fundamental for excitatory stage
				end
			else % single-stream processing
				for i = 1: dirs % loop over directions
					m.p.dir = dir(i); % current direction
					[~, ~, p] = solveF(m); % solve model for cortical potential
					q(:, :, i) = p(1: 2, :, 2); % DC and fundamental for excitatory stage
				end
			end
			p0 = q(1, :, :) / fs; % DC potential (mV)
			p1 = abs(q(2, :, :)) * 2 / fs; % fundamental potential (mV)
			a = p0 + p1; % amplitude of potential (mV): 1 x ks x dirs
			a = m.p.gainGen * max(0, a); % impulse rate (Hz): 1 x ks x dirs
			[a, i] = max(a, [], 3); % maximum over orientations: 1 x ks

			% Adjust weights
			if cycle > 0 % adjust weights relative to previous values
				j = a == aP; % neurons for which the response is unchanged
				w(c, j) = wP(c, j); % leave weight unchanged from previous value
				j = a < aP; % neurons for which the response has decreased
				w(c, j) = max(0, wP(c, j) - wInc); % decrease weight relative to previous
			end
			wP = w; aP = a; % save comparison values for next cycle
			
			% Store
			[t, j] = ismember(cycle, m.p.cycles); % store this cycle?
			if t % yes: store
				dirP = dir(i); % preferred direction: 1 x ks
				dC{j} = table(m.p.wid, cycle, m.p.x, dir, {w}, m.p.gainSI0, ...
					m.p.gainIE0, dirP, 'variableNames', {'width', 'cycle', 'x', ...
					'dir', 'weight', 'gainSI0', 'gainIE0', 'dirP'});
			end
			
			% Report progress
			if etime(clock, timeR) > 120 % 2 minutes after last report
				timeC = etime(clock, timeS); % time since start
				fprintf('Time (s), cycle: %4g, %4g\n', floor(timeC), cycle); % report
				timeR = clock; % report time
			end				

		end
		d = vertcat(dC{:}); % combine data tables
		m.p.dir = dirOrig; % restore original value
		
		% Save weights
		save(m.p.file, 'd'); % save file

	end
	
	% Initialise variable descriptions
	d.Properties.VariableDescriptions{1} = ''; % set all descriptions to empty
	d = addprop(d, 'RespDim', 'table'); % add property: response dimensions

function y = ifftReal(y, n, dim, err)

% Take the inverse transform, assuming that the result should be real

	% Calculate the inverse
	if exist('dim', 'var')
		y = ifft(y, n, dim); % inverse transform
	else
		y = ifft(y); % default dimension: first non-singleton
	end
	yR = real(y); yI = imag(y); % real and imaginary parts
	
	% Check that the imaginary components are not too large and return real
	if ~ exist('err', 'var'), err = 1e-3; end
	e = abs(yI); % magnitude of imaginary component
	e = max(e(:)); % maximum error
	if e > err % it has substantial imaginary components
		error('Inverse transform has substantial imaginary components');
	else
		y = yR; % make it real
	end

function y = modelTuning(b, x)

% Calculate the model for the time course in response to moving stimuli

	r0 = b(1); % resting impulse rate (Hz)
	r1 = b(2); % peak impulse rate 1 (Hz)
	r2 = b(3); % peak impulse rate 2 (Hz)
	k = b(4); % coefficient of exponent
	dirP = b(5); % preferred direction (deg)
	dirA = b(6); % antipreferred direction (deg)
	y = r0 + r1 * exp(k * (cos((pi / 180) * (x - dirP)) - 1)) + ...
		r2 * exp(k * (cos((pi / 180) * (x - dirA)) - 1));

function [d, m] = selDir(d, m)

% Select preferred and antipreferred directions

	% Initialise
	rIn = d.resp; % correlations: ds x ks x ps: ps can be multidimensional
	s = size(rIn); % number of elements along each dimension
	ds = s(1); % number of directions
	s(1) = 2; % we are reducing ds directions to 2
	dir = d.dir; % stimulus directions: ds x 1
	dirP{1} = d.dirP(1, :); % preferred direction (deg): 1 x ks
	dirP{2} = mod(dirP{1} + 180, 360); % antipreferred direction (deg): 1 x ks
	r = zeros(s); % response, keeping only pref. and anti. directions
	
	% Select and combine
	for i = 1: 2 % loop over relative directions: preferred then antipreferred
		for j = 1: ds % loop over stimulus directions
			k = dirP{i} == dir(j); % cells with relative direction in stimulus dir.
			r(i, k, :) = rIn(j, k, :); % transfer from input response
		end
	end
	
	% Store
	d = d(1: 2, :); % keep only first two rows
	d.dirR = categorical({'pref'; 'anti'}); % label relative direction
	d.resp = r; % response: 2 x ks x ps
	d.Properties.CustomProperties.RespDim(1) = {'dirR'}; % response dimensions
		
function [fD, pS, pC] = solveF(m, dir)

%	Solve the model equations in the frequency domain.
%	The stimulus is a drifting or contrast-reversing grating, or sparse noise.
%	The function calculates the discrete Fourier transform of generator potential:
%		pS = subcortical potential with size fs x cs x zs
%		pC = cortical potential, with size fs x ks x zs
%	where fs is the number of frequencies, cs the number of subcortical channels,
%	ks the number of cortical locations, and zs the number of stages.

	% Initialise stimulus parameters and neuronal numbers
	if ~ exist('dir', 'var') % direction is not supplied as argument
		dir = m.p.dir;
	end
	dir = 2 * pi * dir / 360; % direction (radians)
	fS = 2 * pi * m.p.freqS; % spatial frequency (radians / deg)
	fT = 2 * pi * m.p.freqT; % temporal frequency (radians / s)
	phaseS = pi * m.p.phaseS / 180; % spatial phase (radians)
	cs = size(m.p.chan, 1); % number of subcortical channels
	ks = size(m.p.loc, 1); % number of cortical locations
	xs = length(m.p.x); % number of cortical locations along a single dimension
	tauS = m.p.tau + m.p.chan(:, 3)' * m.p.tauDev; % t/c depends on ch. sign
	
	% Calculate frequencies
	fD = m.p.f; fs = length(m.p.f); % frequencies (Hz); number of frequencies
	f = 2 * pi * fD; % radians/s
	f = f'; % column vector
	
	% Zero potentials
	pS = zeros(fs, cs, m.p.stages(1)); % subcortical potential: fs x cs x zs
	pC = zeros(fs, ks, m.p.stages(2)); % cortical potential: fs x ks x zs

	% First subcortical stage: photoreceptors
	switch m.p.stim % stimulus type
		case {'grating', 'reverse'} % grating
			r = m.p.cont * exp(-.25 * (m.p.radCen * fS) ^ 2);
				% dot product of stimulus and receptive field
			u = m.p.chan(:, 1: 2) * [cos(dir); sin(dir)];
				% channel locations (deg from centre in direction of motion)
			switch m.p.stim
				case 'grating' % drifting grating
					r = r * exp(- 1i * u * fS); % phase shift for channel location: cs x 1
				case 'reverse' % contrast-reversing grating
					%	u = abs(u); % direction is irrelevant
					%	r = r * cos(u * fS + phaseS); % amp. shift for channel loc.: cs x 1
					r = r * cos(u * fS - phaseS); % amp. shift for channel loc.: cs x 1
			end
			r = r / (1 + 1i * m.p.tau * fT); % output of first stage
			pS(2, :, 1) = .5 * fs * r; % fundamental fft component
			pS(end, :, 1) = conj(pS(2, :, 1)); % because time course is real
		case 'sparse' % sparse noise
			inp = calInp(m.p.t, m); % dot product: fs x cs
			inp = fft(inp); % Fourier transform: fs x cs
			pS(:, :, 1) = inp ./ (1 + 1i * m.p.tau * f); % output of first stage
	end
	pS(:, :, 1) = - (m.p.gainS / m.p.gainGen) * pS(:, :, 1); % convert to mV
	
	% Later subcortical stages
	pS(:, :, 2) = - m.p.chan(:, 3)' .* pS(:, :, 1) ./ (1 + 1i * tauS .* f);
		% bipolar cells
	pS(:, :, 3) = pS(:, :, 2) ./ (1 + 1i * tauS .* f); % ganglion cells
	pS(1, :, 3) = fs * m.p.actS / m.p.gainGen; % DC fft component
	pS(:, :, 4) = pS(:, :, 3) ./ (1 + 1i * tauS .* f); % lateral geniculate
	
	% Inputs to cortex
	g = reshape(m.p.gainGE, [1, cs, ks]); w = reshape(m.p.weight, [1, cs, ks]);
		% prepare to sum: 1 x cs x ks
	switch m.p.stim % stimulus type
		case {'grating', 'reverse'} % cyclic
			p = pS(1: 2, :, 4); % first two components only, to save time
			inpS = sum(g .* w .* p, 2); % sum over channels: 2 x 1 x ks
			inpS = squeeze(inpS); % remove channel dimension: 2 x ks
			inpS(3: fs, :) = zeros(fs - 2, ks); % restore missing components
			inpS(end, :) = conj(inpS(2, :)); % fs x ks
		case 'sparse' % sparse noise
			%	inpS = tall(pS(:, :, 4)); % geniculate output: fs x cs
			%	write('Datastore', inpS); % create datastore
			inpS = pS(:, :, 4); % geniculate output: fs x cs
			inpS = g .* w .* inpS; % sum over channels: fs x cs x ks
			inpS = sum(inpS, 2); % sum over channels: fs x 1 x ks
			%	inpS = gather(inpS); % return to standard array
			inpS = reshape(inpS, [fs, ks]); % remove channel dimension: fs x ks
	end
	
	% Inhibitory cortical stage
	pC(:, :, 1) = m.p.gainSI0 * inpS ./ (1 + 1i * m.p.tau * f); % fs x ks
	
	% Post-somal inhibitory stage
	inpI = ifft(pC(:, :, 1)); % time course: fs x ks
	inpI = max(0, inpI); % rectified time course: fs x ks
	inpI = fft(inpI); % back to frequency domain: fs x ks
	inpI = inpI ./ (1 + 1i * m.p.tauI * f); % filter: fs x ks
	pC(:, :, 3) = inpI; % store
	
	% Inhibitory input to excitatory stage
	switch m.p.calIE % method
		case 'dot' % dot product
			g = shiftdim(m.p.gainIE, -1); % prepare to sum: 1 x ks x ks
			m = calWeightI(m); % calculate excitatory-to-inhibitory weight
			w = shiftdim(m.p.weightI, -1); % 1 x ks x ks
			inpI = sum(g .* w .* inpI, 2); % sum over inhibitory neurons: fs x 1 x ks
			inpI = squeeze(inpI); % remove inhibitory neuron dimension: fs x ks
		case 'imf' % imfilter: this takes longer than img, don't know why
			inpI = reshape(inpI, [fs, xs, xs]); % fs x xs x xs
			g = shiftdim(m.p.gainIE, -1); % 1 x xs x xs
			inpI = imfilter(inpI, g); % cross-correlate inhib.
			inpI = reshape(inpI, [fs, ks]); % reverse reshape: fs x ks
		case 'img' % imgaussfilt
			inpI = inpI.'; % transpose: ks x fs
			inpI = reshape(inpI, [xs, xs, fs]); % xs x xs x fs
			dev = m.p.radIE / sqrt(2); % standard deviation of Gaussian (deg)
			dev = m.p.densC * dev; % standard deviation of Gaussian (samples)
			inpI = imgaussfilt(real(inpI), dev) + ...
				1i * imgaussfilt(imag(inpI), dev); % filter
			inpI = reshape(inpI, [ks, fs]); % reverse reshape
			inpI = inpI.'; % transpose: fs x ks
	end

	% Excitatory cortical stage
	inpS = m.p.gainSE0 * inpS; % include gains
	inpI = m.p.gainIE0 * inpI;
	pC(:, :, 2) = (inpS - inpI) ./ (1 + 1i * m.p.tau * f); % filter and store
	
	% Store summed inputs
	pC(:, :, 4) = inpI; % store summed inhibitory input as fourth stage
	pC(:, :, 5) = inpS; % store summed excitatory input as fifth stage

function [t, pS, pC] = solveT(m)

% Solve the model equations in the temporal domain.
%	The function calculates the time course of generator potential p(t, c, z),
%	where t is time, c is channel number and z is stage number.

	% Initialise
	cs = size(m.p.chan, 1); % number of subcortical channels
	zsS = m.p.stages(1); % number of subcortical stages
	ks = size(m.p.loc, 1); % number of cortical locations
	zsC = m.p.stages(2); % number of cortical stages

	% Calculate potentials
	fun = @(t, p)calDeriv(t, p, m); % function to calculate derivatives
	tSpan = m.p.t; % simulation times
	pS = zeros(cs, zsS); % initial values for subcortex
	pI = m.p.actS / m.p.gainGen; % subcortical resting potential (mV)
	pS(:, 3: 4) = pI; % ganglion cells, geniculate cells
	pC = zeros(ks, zsC); % initial values for cortex
	pI = m.p.gainSC0 * pI; % inhibitory resting potential
	pC(:, [1, 3, 4, 5]) = pI; % inhibitory soma, axon
	pI = (1 - m.p.gainIEMax) * pI; % excitatory resting potential
	pC(:, 2) = pI; % excitatory neurons
	p0 = [pS(:); pC(:)]; % vectorise
	[t, p] = ode45(fun, tSpan, p0); % solve
	
	% Sort potentials into subcortical, cortical
	t = t'; % row vector
	ts = length(t); % number of times
	i = ts * cs * zsS; % number of subcortical values
	pS = reshape(p(1: i), [ts, cs, zsS]); % subcortical potentials: ts x cs x zs
	pC = reshape(p(i + 1: end), [ts, ks, zsC]); % cortical pot.: ts x ks x zs
	
	% Calculate summed inhibitory input to cortical neurons
	p = pC(:, :, 3); % individual inhibitory responses: ts x ks
	xs = length(m.p.x); % number of cortical neurons in one dimension
	p = reshape(p', [xs, xs, ts]); % prepare for filter: xs x xs x ts
	dev = m.p.densC * m.p.radIE / sqrt(2); % standard dev. of Gaussian (samples)
	p = imgaussfilt(p, dev); % filter with Gaussian: xs x xs x ts
	p = reshape(p, [ks, ts])'; % revert to original shape: ts x ks
	pC(:, :, 4) = m.p.gainIE0 * p; % summed inhibition: ts x ks

	% Calculate summed excitatory input to cortical neurons
	g = shiftdim(m.p.gainGE, -1); % gain from geniculate to cortex: 1 x cs x ks
	w = shiftdim(m.p.weight, -1); % synaptic weights: 1 x cs x ks
	p = pS(:, :, 4); % geniculate output: ts x cs x 1
	p = m.p.gainSE0 * sum(g .* w .* p, 2); % sum over channels: ts x 1 x ks
	pC(:, :, 5) = reshape(p, [ts, ks]); % summed excitation: ts x ks
	