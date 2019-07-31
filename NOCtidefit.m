function [amp,phase,fit,ctl_names,p] = NOCtidefit(mjd,hts,varargin)
% [amp,phase,fit] = NOCtidefit(mjd,hts,varargin)
%
%=======================================================================
%
% Given a time series (hts) length nsize, with corresponding times mjd
% performs a least squares fit of ndc tidal constants with speeds given in the list
% sig, returning corresponding amplitude and phase.
%
% Input:
%
% mjd (double array length nsize)    - modified julian day of each
%                                      point in time series
% hts (double array length nsize)    - input time series
%
% If no other parameters are supplied then the program will assume an
% appropriate set of tides based on the range of data
% If you wish to tell the program what harmonics to use then you can either
% (1) specify a "ctl" file using 'filename'
% (2) input a structure, S, consisting of parts sig, sigrel, relamp, and relpha
%
%      sig (double array len ndc+nrel)  - speeds of primary (1:ndc) and
%                                         related (ndc+1,ndc+nrel)constituents
%    sigrel(double array length nrel)   - (n)=speed of primary constituent
%                                         related to related constituent
%                                         with speed sig(ndc+n)
%    relamp(double array length nrel)   - related/primary amplitude
%    relpha(double array length nrel)   - related-primary phase
%
% If sigrel is not empty then an 
% extra set of tides are included without increasing the number of
% fitted parameters, by assuming a constant amplitude factor (relamp)
% and relative phase (relpha) in comparison with one of the fitted
% constituents.  Most common is for relamp to reflect relative
% amplitudes of the equilibrium constituents, and relpha=0.
%
% sigrel contains 2 lists of speeds: one for the speed of the
% constituent to be calculated by association, the other for the speed of
% the constituent it is related to.
%
% (3) enter a choice for the harmonics out of the following list using 'force'
%      '18yr' suitable for 18 years of data
%      '1yr' suitable for 1 year of data
%      '1yr_nolong' suitable for 1 year without long period harmonics
%      '6mo' suitable for 6 months of data
%      '6mo_nolong' suitable for 6 months of data without long period harmonics
%      '1mo' suitable for 1 month of data
%
% OPTIONAL ARGUMENTS :
%
% These come in pairs with the name and then the variable e.g. 'eldot',F
% or 'verbose','long'
%
% eldot (double array length nsize)     - This is used when we are estimating
%                                         tidal parameters from GNSS-MR data
%                                         where we are taking into account the Hdot
%                                         effect.
% flag (logical array length nsize)    - This allows you to input all the data (including flagged data)
%                                        so that you can output the fit to all the results and not
%                                        just the flagged values - useful if including eldot
%
% verbose (char)                        - Output some extra information to screen
%                                         Can be either :
%                                         on/yes/short : brief output
%                                         no/off       : no output (DEFAULT)  
%                                         long         : Outputs the constituents and their 
%                                                        estimated amplitude and phase
%
%=======================================================================
%
% Output:
%
% amp (double array, 1:ndc+nrel+1)      - (1)                = constant offset
%                                         (2:ndc+1)          =  amplitudes of primary constituents
%                                         (ndc+2:ndc+nrel+1) =  amplitudes of relative constituents
%
% phase (double array, length ndc+nrel) - (1:ndc)            = phases of primary constituents
%                                         (ndc+1:ndc+nrel)   = phases of relative constituents
% fit (double array length nsize)       - least squares fit to the data
%
%
% Edited 2019-May-08 to output ctl_names instead of names - ie the list
% actually used rather than the complete 115 list! 

% Set up some parameters

ncmax  = 120;
idmin  = -1409059;
rad    = pi/180;
deg    = 1.0 ./ rad;

sig0 = [0.0410686;   0.0821373;   0.5443747;   1.0158958;
        1.0980330;  12.8542862;  12.9271398;  13.3986609;
       13.4715145;  13.9430356;  14.0251729;  14.4920521;
       14.5695476;  14.9178647;  14.9589314;  15.0000000;
       15.0410686;  15.0821353;  15.1232059;  15.5125897;
       15.5854433;  16.0569644;  16.1391017;  27.3416965;
       27.4238338;  27.8953548;  27.9682085;  28.4397295;
       28.5125832;  28.9019670;  28.9841042;  29.0662415;
       29.4556253;  29.5284789;  29.9589333;  30.0000000;
       30.0410667;  30.0821373;  30.5443747;  30.6265120;
       31.0158958;  42.9271398;  43.4761564;  43.9430356;
       44.0251729;  45.0410686;  57.4238338;  57.9682085;
       58.4397295;  58.9841042;  59.0662415;  60.0000000;
       60.0821373;  86.4079380;  86.9523127;  87.4238338;
       87.9682085;  88.0503458;  88.9841042;  89.0662415;
       26.4079380;  26.8701754;  26.9523127;  27.5059711;
       28.3575923;  29.9178627;  31.0887494;  42.3827651;
       43.0092771;  44.5695476;  56.8701754;  56.9523127;
       57.8860712;  71.9112441;  72.4602606;  73.0092771;
       84.8476675;  85.3920423;  85.8542797;  85.9364170;
       86.3258007;  86.4807917;  86.8701754;  87.4966874;
       88.5125832;  88.5947205; 114.8476675; 115.3920423;
      115.9364170; 116.4079380; 116.9523127; 117.0344500;
      117.5059711; 117.9682085; 118.0503458; 145.9364170;
      146.9523127; 174.3761465; 174.9205212; 175.9364170;
       27.4966874;  27.8860712;  28.9430356;  29.0251729;
       30.4715211;  31.0980330;  56.4079380;  57.4966874;
       58.5125832;  59.5284789;  28.3986629;  28.4807962;
       72.9271398;  74.0251729;  29.5284789;   0.0000000;
        0.0000000;   0.0000000;   0.0000000;   0.0000000];


names = {'SA';   'SSA'; 'MM'; 'MSF'; 'MF'; '2Q1'; 'SIG1'; 'Q1'; 'RO1'; 
         'O1';   'MP1'; 'M1'; 'CHI1'; 'PI1'; 'P1'; 'S1'; 'K1'; 'PSI1'; 'PHI1'; 
         'TH1';  'J1'; 'SO1'; 'OO1'; 'OQ2'; 'MNS2'; '2N2'; 'MU2'; 'N2'; 'NU2'; 
         'OP2';  'M2'; 'MKS2'; 'LAM2'; 'L2'; 'T2'; 'S2'; 'R2'; 'K2'; 'MSN2'; 
         'KJ2';  '2SM2'; 'MO3'; 'M3'; 'SO3'; 'MK3'; 'SK3'; 'MN4'; 'M4'; 'SN4'; 
         'MS4';  'MK4'; 'S4'; 'SK4'; '2MN6'; 'M6'; 'MSN6'; '2MS6'; '2MK6'; '2SM6'; 
         'MSK6'; '2MN2S2'; '3MSK2'; '3M2S2'; 'MNK2S2'; 'SNK2'; '2SK2'; '2MS2N2'; 
         'MQ3';  '2MP3'; '2MQ3'; '3MK4'; '3MS4'; '2MSK4'; '3MK5'; 'M5'; '3MO5'; 
         '2MNS6';'3MNS6'; '4MK6'; '4MS6'; '2MSNK6'; '2MV6'; '3MSK6'; '4MN6'; '3MSN6'; 
         'MKL6'; '2MN8'; '3MN8'; 'M8'; '2MSN8'; '3MS8'; '3MK8'; 'MSNK8'; '2MS8'; 
         '2MSK8';'4MS10'; '3M2S10'; '4MSN12'; '5MS12'; '4M2S12'; 'MVS2'; '2MK2'; 
         'MA2';  'MB2'; 'MSV2'; 'SKM2'; '2MNS4'; 'MV4'; '3MN4'; '2MSN4'; 'NA2'; 
         'NB2';  'MSO5'; 'MSK5'; '2MN2'}; 

if mod(nargin,2) == 1
	error('Incorrect number of input arguments : arguments need to be in pairs.');
end

% first check the mjd and hts to be common in size and make them column vectors

if ~isvector(mjd)
	error('variable mjd must be a vector');
end
if ~isvector(hts)
	error('variable hts must be a vector');
end

if xor(isrow(mjd),isrow(hts))
	warning('mjd and hts should both be column or row vectors : Converting them to column vectors');
	if isrow(mjd)
		mjd = mjd';
	end
	if isrow(hts) 
		hts = hts';
	end
end

rowcheck = 0;
if isrow(mjd)
	% Both mjd and hts must be row vectors so convert to column vectors for calculation but set tg_fit at end back to row vector to match
	mjd = mjd';
	hts = hts';
	rowcheck = 1;
end

if size(mjd,1) ~= size(hts,1)
	error('mjd and hts matrix dimensions must agree.');
end

nsize = length(mjd);

mjdns = floor(mjd);
hrs = (mjd - mjdns) .* 24;

gotCtl  = 0;
gotBB   = 0;
gotFlag = 0;
verbose = 0;

for j = 1:2:nargin-2
	switch lower(varargin{j})
		case 'filename'
			filename = varargin{j+1};
                        validateattributes(filename,{'char'},{'nonempty'});
			[sig,sigrel,relamp,relpha,ctl_names] = readctl(filename);
			gotCtl = 1;
		case 'struct'
			S = varargin{j+1};
                        if ~isstruct(S)
                                error(sprintf('struct %s is not a structure', inputname(j+1)));
                        end
			sig = S.sig;
			sigrel = S.sigrel;
			relamp = S.relamp;
			relpha = S.relpha;
			ctl_names = S.names;
			gotCtl = 1;
		case 'force'
			force_ctl = varargin{j+1};
                        validateattributes(force_ctl,{'char'},{'nonempty'});
			gotCtl = -1;
		case 'eldot'
			bb = varargin{j+1};
			validateattributes(bb,{'numeric'},{'size',[nsize 1]});
			gotBB = 1;
		case 'flag'
			I = varargin{j+1};
			validateattributes(I,{'numeric','logical'},{'size',[nsize 1]});
			if isnumeric(I)
				I = logical(I);
			end
			gotFlag = 1;
		case 'verbose'
			verbose_cmd = varargin{j+1};
                        validateattributes(verbose_cmd,{'char'},{'nonempty'});
			switch lower(verbose_cmd)
			 	case 'on'
					verbose = 1;
				case 'yes'
					verbose = 1;
				case 'no'
					verbose = 0;
				case 'off'
					verbose = 0;
				case 'long';
					verbose = 2;
				case 'short'
					verbose = 1;
				otherwise 
					error(sprintf('Unknown verbose option %s', verbose_cmd))
			end
		otherwise
			error(sprintf('Unknown input argument : %s',varargin{j}))
	end
end
			
T = range(mjd);

if (gotCtl == 0)
	if T >= 3000
		force_ctl = '18yr';
	elseif T >= 365
		force_ctl = '1yr';
	elseif T >= 180
		force_ctl = '1yr_nolong';
	elseif T >= 90
		force_ctl = '6mo_nolong';
	else
		force_ctl = '1mo';
	end
end

	
if gotCtl <= 0
	switch lower(force_ctl)
		case '18yr'
			p = (1:115)';
			sig = sig0(p);
			sigrel = [];
			relamp = [];
			relpha = [];
			ctl_names = names(p);
		case '1yr'
			p = find(logical([ones(60,1);zeros(42,1);1;1;zeros(11,1)]));
			sig = sig0(p);
			sigrel = [];
			relamp = [];
			relpha = [];
			ctl_names = names(p);
		case '1yr_nolong'
			p = find(logical([zeros(5,1);ones(55,1);zeros(42,1);1;1;zeros(11,1)]));
			sig = sig0(p);
			sigrel = [];
			relamp = [];
			relpha = [];
			ctl_names = names(p);
		case '6mo'
			p       = [2;3;4;5;6;7;8;9;10;11;12;13;15;17;19;20;21;22;23;24;25;26;27;28;29;30;31;32;33;34;36;38;39;40;41;42;43;44;45;46;47;48;49;50;51;52;53;54;55;56;57;58;59;60;14;35];
			sig = sig0(p);
			sigrel  = sig0([17;36]);
			relamp  = [0.019;0.059];
			relpha  = [0.0;0.0];
			ctl_names = names(p);
		case '6mo_nolong'
			p       = [6;7;8;9;10;11;12;13;15;17;19;20;21;22;23;24;25;26;27;28;29;30;31;32;33;34;36;38;39;40;41;42;43;44;45;46;47;48;49;50;51;52;53;54;55;56;57;58;59;60;14;35];
			sig = sig0(p);
			sigrel  = sig0([17;36]);
			relamp  = [0.019;0.059];
			relpha  = [0.0;0.0];
			ctl_names = names(p);
		case '1mo'
			p = [3;4;8;10;12;17;21;23;27;28;31;34;36;41;42;43;45;47;48;49;50;54;55;56;57;59;14;15;18;19;26;29;35;38];
	
			sig    = sig0(p);
			sigrel = sig0([17;17;17;17;28;28;36;36]);
			relamp = [0.019;0.331;0.008;0.014;0.133;0.194;0.059;0.272];
			relpha =  zeros(8,1);
			ctl_names = names(p);
		otherwise
			error(sprintf('unknown choice of harmonics : %s', force_ctl))

	end
end


%
% check validity of input parameters (time range of data, number of constituents)
%

nrel = length(sigrel);
ndc  = length(sig)-nrel;

nfit=2*ndc+1;

amp   = zeros(ndc+nrel+1,1);
phase = zeros(ndc+nrel,1);
fit   = zeros(length(mjdns),1);

if length(relamp) ~= nrel | length(relpha) ~= nrel
	error('relative amplitude and or phase are not the right size')
end


if any(mjdns < idmin)
	error('Modified Julian Date is outside of valid range');
end

if (ndc+nrel >= ncmax |  ndc <= 0 | nrel < 0) 
	error('ndc + nrel > ncmax - too many constituents (or zero)');
end

%
% compile a table matching speeds in sig with those in sig0
%

if any(sig == 0)
	error('Input speeds == 0')
end

ntab = zeros(ndc+nrel,1);

j = 1;
for i = 1:ndc
	while abs(sig(i)-sig0(j)) > 0.00001
		j = j + 1;
		if (j > ncmax)
			error('speeds not compatable or in wrong order compared with full list')
		end
	end
	ntab(i) = j;
end

j = 1;
for i = ndc+1:ndc+nrel
	while abs(sig(i)-sig0(j)) > 0.00001
		j = j + 1;
		if (j > ncmax)
			error('speeds not compatable or in wrong order compared with full list')
		end
	end
	ntab(i) = j;
end

%
% compile a table matching speeds in sigrel with those in sig(1:ndc),
%

ntabrel = zeros(nrel,1);

j = 1;
for i = 1:nrel
	while abs(sigrel(i)-sig(j)) > 0.00001
		j = j + 1;
		if (j > ndc)
			nerr = 2;
			return
		end
	end
	ntabrel(i) = j;
end

%
% calculate mean hmn and standard deviation hsd of (valid values in) hts.
% Fit is performed on a rescaled timeseries, with mean 0 and standard
% deviation 1, and rescaled at the end.
%


if gotFlag

	nvalid = length(hts);

	if sum(I) < 2*ndc+1
		error('Not enough valid data points for fit')
	end

	hmn = sum(hts(I)) / sum(I);
	hsd = sqrt( sum(hts(I).^2) / sum(I) - hmn.^2);
	hts = (hts - hmn) ./ hsd;

else
	nvalid = length(hts);

	if nvalid < 2*ndc+1
		error('Not enough valid data points for fit')
	end

	hmn = sum(hts) / nvalid;
	hsd = sqrt( sum(hts.^2) / nvalid - hmn.^2);
	hts = (hts - hmn) ./ hsd;
end

% Loop through valid data, calculating tides at each time, and thus
% contribution to covar and b. A call to phamp0 causes reference
% phases and nodal amplitude factors to be calculated at 00:00
% on day mjdn0. mjdn0 is chosen to be within ndayfu/2 days of
% prediction time, on a day with mjdn0 divisible by ndayfu.
% mjdn0 is changed only when prediction time is out of range
% mjdn0(+-)ndayfu/2. Change of mjdn0 forces a call to phamp0.
%
% Time measured relative to reference time is tprime (days)
%

b       = zeros(nfit,1);
covar   = zeros(nfit,nfit);

phamp_tic = tic; [f,v] = phamp0fast(mjdns); timespent = toc(phamp_tic);

if verbose == 1
	disp(sprintf('Time taken to estimate phases and nodal amplitude factors : %.1f seconds', timespent))
end

A       = zeros(nvalid,nfit);

A(:,1) = 1.0;
if gotBB == 0
	for j = 1:ndc
		k = ntab(j);
		A(:,2*j)   = f(:,k) .* cos(rad.*(sig0(k).*hrs + v(:,k) ));
		A(:,2*j+1) = f(:,k) .* sin(rad.*(sig0(k).*hrs + v(:,k) ));
	end
else
	for j = 1:ndc
   	     k = ntab(j);

        	% This is the "normal" equation except negative because the heights are distance to sea level from antenna
        	A(:,2*j)   = -f(:,k) .* cos(rad.*(sig0(k).*hrs + v(:,k) ));
        	A(:,2*j+1) = -f(:,k) .* sin(rad.*(sig0(k).*hrs + v(:,k) ));

        	% This is the additonal hdot * tan(elevation) ./ elevation_dot effect
        	A(:,2*j)   = A(:,2*j)   + bb .* 24 .* f(:,k) .* rad .* sig0(k) .* sin(rad.*(sig0(k).*hrs + v(:,k)));
        	A(:,2*j+1) = A(:,j*2+1) - bb .* 24 .* f(:,k) .* rad .* sig0(k) .* cos(rad.*(sig0(k).*hrs + v(:,k)));
	end
end

%
% 		then add in effects of relative constituents
%
if gotBB == 0
	for j = 1:nrel
		k=ntab(j+ndc);
		j2=ntabrel(j);
       		A(:,2*j2)   = A(:,2*j2)   + relamp(j).*f(:,k) .* cos(rad.*( sig0(k).* hrs + v(:,k) + relpha(j) ));
       		A(:,2*j2+1) = A(:,2*j2+1) + relamp(j).*f(:,k) .* sin(rad.*( sig0(k).* hrs + v(:,k) + relpha(j) ));
	end
else
	for j = 1:nrel
	        k=ntab(j+ndc);
        	j2=ntabrel(j);

        	% This is the "normal" equation except negative because the heights are distance to sea level from antenna
        	A(:,2*j2)   = A(:,2*j2)   - relamp(j) .* f(:,k) .* cos(rad.*( sig0(k).* hrs + v(:,k) + relpha(j) ));
        	A(:,2*j2+1) = A(:,2*j2+1) - relamp(j) .* f(:,k) .* sin(rad.*( sig0(k).* hrs + v(:,k) + relpha(j) ));

        	% This is the additonal hdot * tan(elevation) ./ elevation_dot effect
        	A(:,2*j2)   = A(:,2*j2)   + relamp(j) .* f(:,k) .* bb .* 24 .* rad .* sig0(k) .* sin(rad.*( sig0(k).* hrs + v(:,k) + relpha(j) ));
        	A(:,2*j2+1) = A(:,2*j2+1) - relamp(j) .* f(:,k) .* bb .* 24 .* rad .* sig0(k) .* cos(rad.*( sig0(k).* hrs + v(:,k) + relpha(j) ));
	end
end

% Now solve the least squares

if gotFlag 
	b = A(I,:)'*hts(I);
	covar = A(I,:)'*A(I,:);
else
	b = A'*hts;
	covar = A'*A;
end

%
% solve normal equations
%

b = covar \ b;

m = b;
m(1) = hsd * b(1) + hmn;
m(2:end) = m(2:end) * hsd;

fit = A*m;

if rowcheck
	fit = fit';
end

%c b now contains the solution vector x.
%c
%c separate out related constituents: fitted functions are (ignoring
%c nodal terms)
%c
%c C.cos(at) + CB.cos(bt+phi)
%c D.sin(at) + DB.sin(bt+phi)
%c
%c where C,D are elements of x - amplitudes of fit, B = relative amplitude
%c of related constituent, and phi = relative phase of related constituent.
%c
%c The sum of these can be rewritten as
%c
%c   C.cos(at) + D.sin(at) + 
%c   B[D.sin(phi)+C.cos(phi)].cos(bt) + B[D.cos(phi)-C.sin(phi)].sin(bt)
%c
%c so C,D are amplitudes for cos and sin parts of primary constituent,
%c B[D.sin(phi)+C.cos(phi)], B[D.cos(phi)-C.sin(phi)]
%c are amplitudes for cos and sin parts of related constituent.
%c
%c Note, nodal variations of primary and related constituents are NOT assumed
%c to be the same.
%c

for i = 1:nrel
	j  = 2*(i+ndc);
	j2 = 2*ntabrel(i);
	b(j)  = relamp(i) * (  b(j2+1)*sin(rad*relpha(i)) +  b(j2)*cos(rad*relpha(i))  );
	b(j+1)= relamp(i) * (  b(j2+1)*cos(rad*relpha(i)) -  b(j2)*sin(rad*relpha(i))  );
end

%
% convert from sine and cosine amplitudes to amplitude/phase, and
% rescale back to original scaling
%
phase = deg * atan2( b(3:2:end), b(2:2:end-1));
J = phase < 0;
phase(J) = phase(J) + 360;

amp = hsd * sqrt( b(3:2:end).^2 + b(2:2:end-1).^2); 
amp =  [hsd * b(1) + hmn; amp];

if verbose == 2
	disp(sprintf('Z0 = %.3f', amp(1)))
	for j = 1:length(phase)
		disp(sprintf('%6s %9.0f %8.3f %7.2f', ctl_names{j}, sig(j)*1e7, amp(j+1), phase(j) ))
	end
end
		
end
