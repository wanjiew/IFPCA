function options = ifpcaSet(Data, opt) 
% ifpcaSet is to set the options for ifpca_paper.

% Set default options
  options.KSrep = 2000;
  options.kmeansrep = 30;
  options.Nrep = 30;
  options.per = 1;
  options.KS = 'To simulate';
  flag = 0;
% End default options.

% Error when no input data
 if nargin == 0 || isempty (Data)
	error('Error occurred: no data file is found')
 end

% Set default options with input data
  [n, p] = size(Data);
  options.pvalcut = log(p)/p;
% End default options

% When there is input, read the input options
if nargin > 1 && ~isempty(opt)
% List of valid field names
  vfields = fieldnames( options );

% Grab valid fields from user's opts
  for i = 1:length(vfields)
    field = vfields{i};
    if isfield( opt, field );
      options.(field) = opt.(field);
    end
  end
  flag = isfield( opt, 'KS' );
% End grabs
end

% Quick return and display the options when no outputs
 if nargout == 0   
   disp('options:')
   disp( options )
   return
 end

% Simulate KSvalue if it is not put in
 if(flag == 1)
    rep = options.KSrep*p;
    KSvalue = KSsimulate(n, rep);
 end
% End of simulation
