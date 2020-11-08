function [hitRate, p] = binomialSigThresh(n,prob)

hits = 0:n;
p = zeros(size(hits));
for i = 1:length(hits)
    p(i) = binopdf(hits(i),n,prob);
end
half1hits = hits(n/2+1:end);
half1p = p(n/2+1:end);

[a,b] = sort(abs(half1p-0.05));
p1 = half1p(b(1:2));
hits1 = half1hits(b(1:2));


half2hits = hits(1:n/2+1);
half2p = p(1:n/2+1);

[a,b] = sort(abs(half2p-0.05));
p2 = half2p(b(1:2));
hits2 = half2hits(b(1:2));

p = [p1 p2];
hits = [hits1 hits2];

hitRate = hits/n;
