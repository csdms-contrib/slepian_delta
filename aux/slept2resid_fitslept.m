%% SLEPT2RESID_FITSLEPT
% This function is a helper function for SLEPT2RESID.
% This allows us to fit the Slepian coefficients to a polynomial in parallel

% Last modified by
%   williameclee-at-arizona.edu, 7/16/2024

function varargout = slept2resid_fitslept(sleptCoeffs, fitwhat, givenerrors, isSpecialterm, ...
	freq, phase, freqSpec, phaseSpec, dateNormalised, dateExtraNormalised, moredates, ...
	G1, G2, G3, GSpec1, GSpec2, GSpec3, nMonth)
P2ftest = 0;
P3ftest = 0;
% If we have a priori error information, create a weighting matrix,
% and change the G and d matrices to reflect this. Since each
% coefficient has its own weighting, we have to invert them
% separately.
Weight = diag(1 ./ givenerrors);
data = sleptCoeffs;
G1weighted = Weight * G1;
G2weighted = Weight * G2;
G3weighted = Weight * G3;
dataWeighted = Weight * data;

extravalues = zeros([length(dateExtraNormalised), 1]);

% This is in case you request a single special term to be looked at
% in a special way
if isSpecialterm
	GSpec1weighted = Weight * GSpec1;
	GSpec2weighted = Weight * GSpec2;
	GSpec3weighted = Weight * GSpec3;
	nOmega = length(freqSpec);
	th = phaseSpec;
	myomega = freqSpec;
else
	nOmega = length(fitwhat(2:end));

	if nOmega > 0
		myomega = freq;
		th = phase;
	end

end

%% First order polynomial
% Do the fitting by minimizing least squares
if ~isSpecialterm
	% The linear fit with periodics
	coeff1 = (G1weighted' * G1weighted) ...
		\ (G1weighted' * dataWeighted);
else
	% If there was a special one, substitute
	coeff1 = (GSpec1weighted' * GSpec1weighted) ...
		\ (GSpec1weighted' * dataWeighted);
end

data1fit = coeff1(1) + coeff1(2) * dateNormalised;

% Add the periodic components (if they exist)
if nOmega > 0
	startP = length(coeff1) - 2 * nOmega + 1;
	% The amplitude of the periodic components
	amp1 = [coeff1(startP:startP + nOmega - 1), ...
				coeff1(startP + nOmega:end)];
	amp1 = sqrt(amp1(:, 1) .^ 2 + amp1(:, 2) .^ 2);

	% The phase in time of the periodic components
	phase1 = [coeff1(startP:startP + nOmega - 1), ...
				  coeff1(startP + nOmega:end)];
	phase1 = atan2(phase1(:, 1), phase1(:, 2));

	% Adding the periodic components to the fit
	data1fit = data1fit ...
		+ sum(repmat(amp1, 1, nMonth) ...
		.* sin(th' + repmat(phase1, 1, nMonth)), 1);
end

% Compute the residual time series for this coefficient
resid1 = data - data1fit';

% Here's the definition of the signal at lm vs time
sleptCoeffsSignal = data1fit;
% Here's the definition of the residual at lm vs time
sleptCoeffsResid = resid1;

% Do extra dates if you have them
if moredates
	extravaluesfn1 = coeff1(1) + coeff1(2) * (dateExtraNormalised);

	if nOmega > 0
		th_extras = repmat(myomega, length(dateExtraNormalised), 1) ...
			* 2 * pi .* repmat((dateExtraNormalised)', 1, nOmega);
		% Evaluate at the missing dates
		extravaluesfn1 = extravaluesfn1 ...
			+ sum(repmat(amp1, 1, 1) .* sin(th_extras' + repmat(phase1, 1, 1)), 1);
	end

	extravalues = extravaluesfn1;
end

% Get the residual sum of squares for later F tests
rss1 = sum(resid1 .^ 2);

%% Second order polynomial
% Now repeat that all again with second order polynomial, if we want
if fitwhat(1) >= 2

	if ~isSpecialterm
		coeff2 = (G2weighted' * G2weighted) ...
			\ (G2weighted' * dataWeighted);
	else
		coeff2 = (GSpec2weighted' * GSpec2weighted) ...
			\ (GSpec2weighted' * dataWeighted);
	end

	data2fit = coeff2(1) + coeff2(2) * (dateNormalised) ...
		+ coeff2(3) * (dateNormalised) .^ 2;

	if nOmega > 0
		startP = length(coeff2) - 2 * nOmega + 1;
		amp2 = [coeff2(startP:startP + nOmega - 1), ...
					coeff2(startP + nOmega:end)];
		amp2 = sqrt(amp2(:, 1) .^ 2 + amp2(:, 2) .^ 2);

		phase2 = [coeff2(startP:startP + nOmega - 1), ...
					  coeff2(startP + nOmega:end)];
		phase2 = atan2(phase2(:, 1), phase2(:, 2));

		data2fit = data2fit + ...
			sum(repmat(amp2, 1, nMonth) ...
			.* sin(th' + repmat(phase2, 1, nMonth)), 1);
	end

	% Compute the residual time series for this coefficient
	resid2 = data - data2fit';

	% Do extra dates if you have them
	if moredates
		extravaluesfn2 = coeff2(1) + coeff2(2) * (dateExtraNormalised) ...
			+ coeff2(3) * (dateExtraNormalised) .^ 2;

		if nOmega > 0
			th_extras = repmat(myomega, length(dateExtraNormalised), 1) ...
				* 2 * pi .* repmat((dateExtraNormalised)', 1, nOmega);
			extravaluesfn2 = extravaluesfn2 ...
				+ sum(repmat(amp1, 1, 1) ...
				.* sin(th_extras' + repmat(phase1, 1, 1)), 1);
		end

	end

	% Get the residual sum of squares
	rss2 = sum(resid2 .^ 2);
	% Calculate an F-score for this new fit
	fratioP2 = (rss1 - rss2) / 1 / (rss2 / (length(sleptCoeffs) - length(coeff2)));
	fscore = finv(.95, 1, length(sleptCoeffs) - length(coeff2));

	if fratioP2 > fscore
		P2ftest = 1;
		% We pass, so update the signal and residuals with this new fit
		sleptCoeffsResid = resid2;
		sleptCoeffsSignal = data2fit;

		if moredates
			extravalues = extravaluesfn2;
		end

	else
		P2ftest = 0;
	end

end

%% Third order polynomial
% Now repeat that all again with third order polynomial, if we want
if fitwhat(1) >= 3

	if ~isSpecialterm
		coeff3 = (G3weighted' * G3weighted) ...
			\ (G3weighted' * dataWeighted);
	else
		coeff3 = (GSpec3weighted' * GSpec3weighted) ...
			\ (GSpec3weighted' * dataWeighted);
	end

	data3fit = coeff3(1) + coeff3(2) * (dateNormalised) ...
		+ coeff3(3) * (dateNormalised) .^ 2 ...
		+ coeff3(4) * (dateNormalised) .^ 3;

	if nOmega > 0
		startP = length(coeff3) - 2 * nOmega + 1;
		amp3 = [coeff3(startP:startP + nOmega - 1), ...
					coeff3(startP + nOmega:end)];
		amp3 = sqrt(amp3(:, 1) .^ 2 + amp3(:, 2) .^ 2);

		phase3 = [coeff3(startP:startP + nOmega - 1), ...
					  coeff3(startP + nOmega:end)];
		phase3 = atan2(phase3(:, 1), phase3(:, 2));

		data3fit = data3fit + ...
			sum(repmat(amp3, 1, nMonth) ...
			.* sin(th' + repmat(phase3, 1, nMonth)), 1);
	end

	resid3 = data - data3fit';

	% Do extra dates if you have them
	if moredates
		extravaluesfn3 = coeff3(1) + coeff3(2) * (dateExtraNormalised) ...
			+ coeff3(3) * (dateExtraNormalised) .^ 2 ...
			+ coeff3(4) * (dateExtraNormalised) .^ 3;

		if nOmega > 0
			th_extras = repmat(myomega, length(dateExtraNormalised), 1) ...
				* 2 * pi .* repmat((dateExtraNormalised)', 1, nOmega);
			extravaluesfn3 = extravaluesfn3 ...
				+ sum(repmat(amp1, 1, 1) ...
				.* sin(th_extras' + repmat(phase1, 1, 1)), 1);
		end

	end

	% Get the residual sum of squares
	rss3 = sum(resid3 .^ 2);
	% Calculate an F-score for this new fit
	fratioP3 = (rss1 - rss3) / 2 / (rss3 / (length(sleptCoeffs) - length(coeff3)));
	fscore = finv(.95, 1, length(sleptCoeffs) - length(coeff3));

	if fratioP3 > fscore
		P3ftest = 1;
		% We pass, so update the signal and residuals with this new fit
		sleptCoeffsResid = resid3;
		sleptCoeffsSignal = data3fit';

		if moredates
			extravalues = extravaluesfn3';
		end

	else
		P3ftest = 0;
	end

end

% Make the matrix ftests
ftests = [0, P2ftest, P3ftest];

varargout = {sleptCoeffsSignal, sleptCoeffsResid, ftests, extravalues};
end