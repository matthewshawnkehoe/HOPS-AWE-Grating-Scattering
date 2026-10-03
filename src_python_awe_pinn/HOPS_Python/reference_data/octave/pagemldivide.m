function X = pagemldivide(A,B)
% Octave shim for MATLAB's pagemldivide
n = size(A,3);
X = zeros(size(A,2),size(B,2),n);
for k=1:n
  X(:,:,k) = A(:,:,k)\B(:,:,k);
end
end
