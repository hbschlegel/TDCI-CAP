# docs/Makefile

####################
# Variables
####################

# List all markdown sources in the order they should appear
SOURCES = src/manual.md

# Output filenames
PDF_OUTPUT  = build/manual.pdf
HTML_OUTPUT = build/index.html

# Template files
LATEX_HEADER = src/header.tex
HTML_HEADER  = src/header.html

# Pandoc & Engine
PANDOC     = pandoc
PDF_ENGINE = xelatex

####################
# Build Targets
####################

# Default target: build both PDF and HTML
all: pdf html

pdf: $(PDF_OUTPUT)
html: $(HTML_OUTPUT)

# Build the PDF by merging all sources in order
$(PDF_OUTPUT): $(SOURCES) $(LATEX_HEADER)
	$(PANDOC) $(SOURCES) \
      --verbose \
      -o debug.tex \
	    -o $(PDF_OUTPUT) \
	    --pdf-engine=$(PDF_ENGINE) \
	    --resource-path=.:src \
	    -H $(LATEX_HEADER) \
      -V tables=yes -V longtable \
      --lua-filter=src/landscape-tables.lua \
	    --toc --toc-depth=2

# Build single-page HTML by merging all sources
$(HTML_OUTPUT): $(SOURCES) $(HTML_HEADER)
	$(PANDOC) $(SOURCES) \
      --verbose \
	    -o $(HTML_OUTPUT) \
	    --mathjax \
	    --resource-path=.:src \
	    --include-in-header=$(HTML_HEADER) \
	    --toc --toc-depth=2

# Clean up generated files
clean:
	rm -f $(PDF_OUTPUT) $(HTML_OUTPUT)

.PHONY: all pdf html clean
