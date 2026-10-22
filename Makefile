# The book executes ~110 chunks unconditionally and ~83 more behind switches,
# against real data. It is rendered locally on purpose: the Python-backed
# steps make a CI run far too heavy. That means the coverage the book provides
# for 69 exported functions is only as good as someone remembering to render --
# so every render leaves a record (`render-record.md`) naming the
# sclet commit and how many chunks each chapter actually ran. Diff two records
# to see a chapter that quietly stopped executing.
.PHONY: gitbook pdfbook epub clean cover record

gitbook:
	SCLET_BOOKDOWN_EVAL=true R_BOOKDOWN_RENDER=1 \
	  Rscript -e 'library(bookdown); render_book("index.Rmd", "gitbook", output_dir = "_book_full")'
	$(MAKE) record

pdfbook:
	SCLET_BOOKDOWN_EVAL=true R_BOOKDOWN_RENDER=1 \
	  Rscript -e 'library(bookdown); render_book("index.Rmd", "pdf_book")'
	$(MAKE) record

epub:
	SCLET_BOOKDOWN_EVAL=true R_BOOKDOWN_RENDER=1 \
	  Rscript -e 'library(bookdown); render_book("index.Rmd", "epub_book")'

# What the last render executed. Cheap (reads the .Rmd sources; no knitting).
record:
	Rscript record-render.R . render-record.md

clean:
	Rscript -e 'bookdown::clean_book()';\
	rm -rf _book _bookdown_files _main.Rmd docs/sclet_cache

cover:
	Rscript -e 'source("book-cover.R")'
