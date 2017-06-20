JAXO = *.xml
EPS = $(patsubst %.xml,%.eps,$(wildcard $(JAXO)))
PDF = $(patsubst %.xml,%.pdf,$(wildcard $(JAXO)))

all:	$(PDF)

clean:
	rm -f $(EPS)
	rm -f *.pdf

%.pdf:	%.eps
	epstopdf $<

%.eps:	%.xml
	jaxodraw --convert `readlink -f $<` `readlink -f $*.tex`
	tex2correct $*.tex
	latex $*.tex
	dvips -x 400 $*.dvi
	ps2eps -f $*.ps
	rm $*.tex $*.tex.bak $*.tex.xml $*.dvi $*.aux $*.log $*.ps

