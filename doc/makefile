LC = pdflatex
PDFV = evince

all : TrigenicInteractionPredictor.pdf

TrigenicInteractionPredictor.pdf : TrigenicInteractionPredictor.tex
	$(LC) TrigenicInteractionPredictor.tex

run : TrigenicInteractionPredictor.pdf
	$(PDFV) TrigenicInteractionPredictor.pdf

.PHONY : clean
clean :
	rm -f TrigenicInteractionPredictor.pdf TrigenicInteractionPredictor.toc TrigenicInteractionPredictor.log TrigenicInteractionPredictor.dvi TrigenicInteractionPredictor.aux TrigenicInteractionPredictor.synctex.gz
