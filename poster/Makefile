# Makefile for xray periodic table

.PHONY: pdf clean distclean

pdf:
	make ptable_xray.pdf

atoms.tex:  gen_atoms.py
	python gen_atoms.py > atoms.tex

ptable_xray.pdf: atoms.tex ptable_xray.tex elines.tex
	pdflatex ptable_xray

clean:
	rm -f *.log *.aux *.pyc *.synctex.gz

distclean:
	make clean
	rm -f ptable_xray.pdf atoms.tex
