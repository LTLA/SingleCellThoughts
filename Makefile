GENERAL_REPORTS = general/bootstrapping.html general/cdr.html general/clustering.html general/linearity.html general/standardization.html
SOFTWARE_REPORTS = software/correlations/corsim.html software/marker_detection/comments.html
WORKFLOW_REPORTS = workflows/qc.html 

all: index.html $(GENERAL_REPORTS) $(SOFTWARE_REPORTS) $(WORKFLOW_REPORTS)

index.html: index.Rmd general.Rmd software.Rmd workflows.Rmd
	echo "rmarkdown::render_site()" | R --no-save --slave
	mv _site/*.html .
	mv _site/site_libs .

$(GENERAL_REPORTS): %.html: %.Rmd
	cd general && echo 'rmarkdown::render(basename("$<"))' | R --no-save --slave

$(SOFTWARE_REPORTS): %.html: %.Rmd
	cd `dirname $<` && echo 'rmarkdown::render(basename("$<"))' | R --no-save --slave

$(WORKFLOW_REPORTS): %.html: %.Rmd
	cd workflows && echo 'rmarkdown::render(basename("$<"))' | R --no-save ----slave

clean:
	rm -f *.html $(GENERAL_REPORTS) $(SOFTWARE_REPORTS) $(WORKFLOW_REPORTS)
	rm -rf site_libs/
