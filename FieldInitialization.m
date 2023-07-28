function [Initialization] = FieldInitialization(d1, d1sub, bwImg, narrowBand, alignment_factor, smoothness_factor,...
    alignment_wt_init, smoothness_wt_init)
%% Docstring
%{
  Takes in a grayscale image (bwImg), its black and white filtered
  image (narrowBand), and the d1 operator on the narrowBand and returns 
  a vector field over the image. 
  Let n = # of black pixels. Then, Initialization is a 2n x 1 array where
  the first n components are the x components of the vector field and the
  last n components are the y-components of the vector field.
%}
%% Initialize Field
%D = d1 * d1';
D = d1sub * d1sub';
n = size(D,1);
DTensored = sparse(2*n, 2*n);
DTensored(1:n, 1:n) = D;
DTensored(n+1:end,n+1:end) = D;

Smoothness = @(f) f' * DTensored * f;
[Gmag, Gdir] = imgradient(bwImg);

% Floor all gradients that are very small/ not on the narrowBand
% Make all other gradients 1
Gmag(Gmag < 0.01) = 0;
Gmag(~narrowBand) = 0;
Gmag(Gmag ~= 0) = 1;

%convert deg to radians
%ComplexG = Gmag .* (exp(1i .* Gdir .* (pi / 180)));
ComplexG = (Gmag.^0.25) .* (exp(1i .* Gdir .* (pi / 180)));
ComplexG = ComplexG(narrowBand);

F = @(f) f(1:length(f) / 2) + 1i * f(length(f) / 2 + 1 : end);

alignmentTerm = @(f) sum(abs(ComplexToLong(f, Gmag, narrowBand) - ComplexG(Gmag(narrowBand) ~= 0).^4).^2);

ObjectiveFunction = @(f) alignmentTerm(f) + Smoothness(f);

%% Linear Solve

cG4 = ComplexG.^4;
g = [real(cG4); imag(cG4)];
A = DTensored;
maskRow = ones(nnz(narrowBand), 1);
I = 1:2 * nnz(narrowBand);
zeroGradInd = find(ComplexG == 0);
maskRow(zeroGradInd) = 0;
maskRow = [maskRow; maskRow];
Mask = sparse(I, I, maskRow);

%lambda is the alignment weight 
lambda = 0.2;
lambda_normalized = lambda / alignment_factor; 

Initialization = ((smoothness_wt_init/smoothness_factor)*A +...
    (alignment_wt_init/alignment_factor)* Mask) \ ((alignment_wt_init/alignment_factor)* Mask * g);

end

