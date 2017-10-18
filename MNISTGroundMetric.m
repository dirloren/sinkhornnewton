function M = MNISTGroundMetric(gridSize, h)
%MNISTGroundMetric - generate euclidean distance matrix for 2D grid
%
% Inputs:
%    gridSize - total number of grid  points
%    h - discretization step length
%
% Outputs:
%    M - gridSize x gridSize distance matrix

% Author: Christoph Brauer
% email: ch.brauer@tu-bs.de
% Website: https://www.tu-braunschweig.de/iaa/personal/brauer
% October 2017; Last revision: 5-October-2017

sqrtGridSize = sqrt(gridSize);

% assert that the number of grid points is square
assert(mod(sqrtGridSize, 1) == 0, 'gridSize must be square');

% preallocate distance matrix
M = zeros(gridSize);

% generate distance matrix
for i0 = 1:sqrtGridSize
    for j0 = 1:sqrtGridSize
        for i = 1:sqrtGridSize
            for j = 1:sqrtGridSize
                M((j0 - 1) * sqrtGridSize + i0, (j - 1) * sqrtGridSize + i) ...
                    = h * sqrt((i - i0)^2 + (j - j0)^2);
            end
        end
    end
end