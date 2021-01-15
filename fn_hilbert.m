function h = fn_hilbert(s)
%SUMMARY
%   Performs Hilbert transform
%USAGE
%   h = fn_hilbert(s)
%AUTHOR
%   Paul Wilcox (2007)
%INPUTS
%   s = vector or matrix to transform. For a matrix, the Hilbert transform 
%   is performed on each column
%OUTPUTS
%   h = Hilbert transform of s

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[s, n] = shiftdim(s);
s1 = size(s,1);
p = 2 ^ nextpow2(size(s, 1));
sp = fft(s, p);
sp = sp(1:floor(p/2), :);
h = ifft(sp, p) * 2;
h = h(1:s1,:);
h = shiftdim(h, -n);

return;