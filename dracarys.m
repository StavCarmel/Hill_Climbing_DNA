function Y = dracarys(X)

X = upper(X);
assert(isa(X, 'char'),'Input should be a char vector');

if mod(length(X), 3)
    error('Coding sequence not divisible by 3');
end

if any(X ~= 'A' & X ~= 'C' & X ~= 'G' & X ~= 'T')
    error('Input is not a valid DNA sequence');
end

if ~strcmp(X(1:3), 'ATG')
    error('Missing a start codon');
end

if ~(strcmp(X(end-2:end), 'TAG') || strcmp(X(end-2:end), 'TGA') || strcmp(X(end-2:end), 'TAA'))
    error('Missing a stop codon');
end

Y = (sum(X(3:3:end-1) == X(4:3:end)))^2 / sqrt(length(X)) * 28 * pi;