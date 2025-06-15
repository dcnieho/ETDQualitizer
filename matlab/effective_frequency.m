function freq = effective_frequency(x, y, duration)
N_valid = sum(~(isnan(x) | isnan(y)));
freq    = N_valid/duration;
