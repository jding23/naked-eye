function downsampledPoints = randomDownsample(points, percentage)
    % points: N-by-3 matrix where each row represents the coordinates of a point
    % percentage: Downsampling ratio, e.g., 0.1 means retaining 10% of the points
    N = size(points, 1);
    numPointsToKeep = round(N * percentage);
    shuffledIndices = randperm(N);
    indicesToKeep = shuffledIndices(1:numPointsToKeep);
    downsampledPoints = points(indicesToKeep, :);
end