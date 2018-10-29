all: index.html

index.html: index.Rmd general.Rmd software.Rmd workflows.Rmd
	echo "rmarkdown::render_site()" | R --no-save --slave
	mv _site/*.html .
	mv _site/site_libs .
