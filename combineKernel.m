function ComK = combineKernel(varargin)

[row, col] = size(varargin{1}{1});
ComK = zeros(row,col);
for i=1: nargin-1;
    ComK = ComK + varargin{1}{i}*varargin{i+1};
end
