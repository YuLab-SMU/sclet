# The book executes ~110 chunks unconditionally and ~83 more behind switches,
# against real data. It is rendered locally on purpose: the Python-backed
# steps make a CI run far too heavy. That means the coverage the book provides
# for 69 exported functions is only as good as someone remembering to render --
# so every render leaves a record (`render-record.md`) naming the
# sclet commit and how many chunks each chapter actually ran. Diff two records
# to see a chapter that quietly stopped executing.
.PHONY: gitbook gitbook-ai pdfbook epub clean cover record

gitbook:
	SCLET_BOOKDOWN_EVAL=true R_BOOKDOWN_RENDER=1 \
	  Rscript -e 'library(bookdown); render_book("index.Rmd", "gitbook", output_dir = "_book_full")'
	$(MAKE) record

# The same render with the AI chapters executing. Five chunks in
# `ai-in-package` call a provider; the rest of the AI material runs on every
# render. This costs real credits and the prose changes with the model, which
# is why it is a separate target and not the default -- but it is the only way
# to check that the AI path still works end to end, so run it when the AI layer
# changes and diff `render-record.md` afterwards.
gitbook-ai:
	SCLET_BOOKDOWN_EVAL=true SCLET_BOOKDOWN_EVAL_AI=1 R_BOOKDOWN_RENDER=1 \
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
