function f = flag_synthesis(flmn, L, N, varargin)


p = inputParser;
p.addRequired('flmn', @isnumeric);          
p.addRequired('L', @isnumeric);          
p.addRequired('N', @isnumeric);   
p.addParamValue('Reality', false, @islogical);
p.parse(flmn, L, N, varargin{:});
args = p.Results;

% Computing inverse transform.
f = flag_synthesis_mex(flmn, L, N, args.Reality);
     

  