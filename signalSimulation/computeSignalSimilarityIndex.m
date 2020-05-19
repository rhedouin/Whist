% Compare signals
function similarity = computeSignalSimilarityIndex(signal_1, signal_2)

    diff_signal = signal_1 - signal_2;
    similarity = sqrt(dot(diff_signal, diff_signal) / length(diff_signal));
end