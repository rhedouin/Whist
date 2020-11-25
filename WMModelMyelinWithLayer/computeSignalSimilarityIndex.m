function RMSE = computeSignalSimilarityIndex(S1, S2)
    RMSE = sqrt(abs(dot(S1 - S2, S1 - S2)) / length(S1));
end

