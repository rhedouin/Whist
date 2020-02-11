% test_averageOrientationsLogEuclidean
close all

orientations = randn(3, 3);
orientations = orientations ./ vecnorm(orientations')';

avg_orientation = averageOrientationsLogEuclidean(orientations);

visualize3dArrow(orientations)
visualize3dArrow(avg_orientation, 0, 'r')

