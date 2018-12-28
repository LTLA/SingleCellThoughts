REPORTS = ./index.html \
    general/bootstrapping.html general/cdr.html general/clustering.html general/linearity.html general/standardization.html \
    software/correlations/corsim.html software/marker_detection/comments.html software/doublet_detection/bycell.html software/doublet_detection/bycluster.html \
    workflows/qc.html 

all: $(REPORTS)

$(REPORTS): %.html: %.Rmd
	cd $(shell dirname $<) && R --no-save --slave -e 'rmarkdown::render(basename("$<"))'

clean:
	rm -f $(REPORTS)
