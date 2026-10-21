gitbook:
	R_BOOKDOWN_RENDER=1 Rscript -e 'library(bookdown); render_book("index.Rmd", "gitbook")'

pdfbook:
	R_BOOKDOWN_RENDER=1 Rscript -e 'library(bookdown); render_book("index.Rmd", "pdf_book")'

epub:
	R_BOOKDOWN_RENDER=1 Rscript -e 'library(bookdown); render_book("index.Rmd", "epub_book")'

clean:
	Rscript -e 'bookdown::clean_book()';\
	rm -rf _book _bookdown_files _main.Rmd docs/sclet_cache

cover:
	Rscript -e 'source("book-cover.R")'
