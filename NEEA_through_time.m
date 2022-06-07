function [LRSmax] = NEEA_through_time(bs,ds,slope, stand_times,R_data_out)

LRSmax = zeros(1,length(stand_times)); % set up storage
bseq = 2:0.01:6; % set up a series of b strategies to test

for j = 1:length(stand_times)
    N = R_data_out(j:length(stand_times)); % grab the time series of population sizes from time t to the end
    for k = 1:length(bseq)
        LRS = -(bseq(k)-bs.*N) .* exp(-cumsum((slope*bseq(k)^2+ds.*N).*0.01)).*0.01;
        % check equation and negatives!!
    end
    max(LRS)
    find(max(LRS))
    LRSmax(j) =  bseq(LRS == max(LRS)); % find the maximum
end