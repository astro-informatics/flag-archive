function flmn = flag_analysis(f, L, N, varargin)


p = inputParser;
p.addRequired('f', @isnumeric);          
p.addRequired('L', @isnumeric);          
p.addRequired('N', @isnumeric);   
p.addParamValue('Reality', false, @islogical);
p.parse(f, L, N, varargin{:});
args = p.Results;

% Computing inverse transform.
flmn = flag_analysis_mex(f, L, N, args.Reality);
     

  