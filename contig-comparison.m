#!/usr/local/bin/octave -q

[k n50 b] = textread('contig-comparison.mat', '%s %d %d');

plot(b, n50, 'bx')
axis([0 5000 0 25000])

for i = 1:length(k)
  text(b(i), n50(i), k(i))
end

title(sprintf('Aligned contig N50 vs breakpoints of NA12878'))
xlabel 'Contig breakpoints'
ylabel 'Aligned contig N50 (bp)'

print -F:16 "contig-comparison.png"
