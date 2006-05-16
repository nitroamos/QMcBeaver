WEBSITE = /www/QMcBeaver
REMOTE_DOXY_PDF_NAME = qmc_doc.pdf

all: html/index.html beaver.pdf

html/index.html : ../include/* ../src/*
	@echo "Updating help system..."
	@doxygen

beaver.pdf : beaver.tex
	@echo "Updating beaver.pdf..."
	@pdflatex beaver.tex
	@bibtex beaver
	@pdflatex beaver.tex
	@pdflatex beaver.tex

.PHONY:website
website: $(WEBSITE)/$(REMOTE_DOXY_PDF_NAME) $(WEBSITE)/html $(WEBSITE)/beaver.pdf

$(WEBSITE)/$(REMOTE_DOXY_PDF_NAME) : latex/refman.pdf
	@cp latex/refman.pdf $(WEBSITE)/$(REMOTE_DOXY_PDF_NAME)
	@chmod a+r $(WEBSITE)/$(REMOTE_DOXY_PDF_NAME)

$(WEBSITE)/beaver.pdf : beaver.pdf
	@cp beaver.pdf $(WEBSITE)/.
	@chmod a+r $(WEBSITE)/beaver.pdf

$(WEBSITE)/html : html/index.html
	@cp -r html $(WEBSITE)/.
	@chmod -R a+rx $(WEBSITE)/html

