#!/usr/local/bin/octave -q

[k n50 b] = textread('scaffold-comparison.mat', '%s %d %d');

plot(b, n50, 'bx')
axis([0 4000 0 2500000])

for i = 1:length(k)
  text(b(i), n50(i), k(i))
end

title(sprintf('Aligned scaffold N50 vs breakpoints of NA12878'))
xlabel 'Scaffold breakpoints'
ylabel 'Aligned scaffold N50 (bp)'

print -F:16 "scaffold-comparison.png"
