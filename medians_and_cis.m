function A = medians_and_cis(upper_ci_level,lower_ci_level,sequence)

% this function calculates the median and upper and lower
% cis of a time sequence generated by a GEM. This
% only works for sequences that occur at a set of standardized times.

ci_up = prctile(sequence,upper_ci_level);
ci_down = prctile(sequence,lower_ci_level);
percentile_50 = prctile(sequence,50);

A = [percentile_50; ci_up; ci_down];