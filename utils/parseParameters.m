function parser = parseParameters(parser,varargin)

addParameter(parser,'test','modified',...
                @(x) any(validatestring(x,{'asymptotic','finite','modified'})));
              
addParameter(parser,'varianceEstimator','bartlett',...
                @(x) any(validatestring(x,{'none','bartlett','roy'})));
              
addParameter(parser,'surrogates',5000,@isnumeric);

addParameter(parser,'taperMethod','none',...
                @(x) any(validatestring(x,{'none','tukey','parzen','bartlett'})));
              
addParameter(parser,'multivariateBartlett',false,...
                @(x) (islogical(x) && isscalar(x)));
            
addParameter(parser,'p','auto',@ischar);
addParameter(parser,'q','auto',@ischar);

parse(parser,varargin{:});