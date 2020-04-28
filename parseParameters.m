function parser = parseParameters(parser,varargin)

addParameter(parser,'test','exact',...
                @(x) any(validatestring(x,{'exact','standard','asymptotic'})));
              
addParameter(parser,'varianceEstimator','bartlett',...
                @(x) any(validatestring(x,{'none','bartlett','roy'})));
              
addParameter(parser,'surrogates',5000,@isnumeric);

addParameter(parser,'taperMethod','none',...
                @(x) any(validatestring(x,{'none','tukey','parzen','bartlett'})));
              
addParameter(parser,'multivariateBartlett',false,...
                @(x) (islogical(x) && isscalar(x)));

parse(parser,varargin{:});