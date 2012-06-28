figures= \
	fork.png tip.png repeat.png transitive.png bubble.png \
	weakedge.png ambiguousedge.png \
	contig-comparison.png scaffold-comparison.png

all: scaffold.html

clean:
	rm -f scaffold.html $(figures)

.PHONY: all clean

scaffold.html: $(figures)

%.html: %.md
	multimarkdown $< >$@

%.pdf: %.html
	cupsfilter $< >$@ 2>/dev/null

%.png: %.dot
	dot -Tpng $< >$@

%.png: %.m %.mat
	octave -q $<
