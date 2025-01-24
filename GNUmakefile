#
# Example makefile using mjotex and a BibTeX references database.
#
# After defining $(PN) to be the name of your document, this will
# compile the source file $(PN).tex to $(PN).pdf.
#
# To support workflows where the document is kept open in a PDF reader
# during recompilation, the build process takes place in a separate
# $(BUILDDIR) directory that defaults to "build". Only when the build
# has completed do we move the resulting $(PN).pdf to the current
# directory. This avoids renaming or modifying the open document in
# place, which in turn helps prevent flickering, loss of bookmarks, et
# cetera, in several PDF readers.

# The build directory. Not broken? Don't fix.
BUILDDIR=build

# The latex compiler. The SOURCE_DATE_EPOCH=0 prevents the creation
# and modification dates from being embedded as metadata into the
# output file; that in turn is important because it allows us to tell
# when the output stops changing (that is, when we are done). The
# variable is supported in pdftex v1.40.17 and later.
LATEX = SOURCE_DATE_EPOCH=0 pdflatex -file-line-error -halt-on-error --output-directory $(BUILDDIR)

# The name of this document.
#
# For example, to use the name of our parent directory:
#
# PN = $(notdir $(realpath .))
#
PN = the_uniqueness_of_lyapunov_rank_among_symmetric_cones

# A space-separated list of bib files. These must all belong to paths
# contained in your $BIBINPUTS environment variable.
#
# Leave commented if you don't use a bibliography database.
#
BIBS = local-references.bib

# A space-separated list of the mjotex files that you use. The path to
# mjotex must be contain in your $TEXINPUTS environment variable.
#
# MJOTEX  = mjotex.sty
#
MJOTEX  = mjo-algebra.tex mjo-common.tex mjo-cone.tex mjo-font.tex
MJOTEX += mjo-hurwitz.tex mjo-proof_by_cases.tex mjo-set.tex
MJOTEX += mjo-theorem.tex mjo.bst

# Use kpsewhich (from the kpathsea suite) to find the absolute paths
# of the bibtex/mjotex files listed in in $(BIBS)/$(MJOTEX). The SRCS
# variable should contain all (Bib)TeX source files for the document.
SRCS = $(PN).tex
ifdef BIBS
BIBPATHS = $(shell kpsewhich $(BIBS))
SRCS += $(BIBPATHS)
endif
ifdef MJOTEX
MJOTEXPATHS = $(shell kpsewhich $(MJOTEX))
SRCS += $(MJOTEXPATHS)
endif

# The first target is the default, so put the PDF document first.
#
# This voodoo is all designed to find a "fixed point" of calling
# $(LATEX). When you build a LaTeX document, it requires an unknown
# number of compilation passes. How do you know when to stop? Easy,
# stop when the output file stops changing! But how to encode that
# in a makefile?
#
# At the start of this rule, we call $(LATEX) to compile $(PN).tex.
# It stores the resulting PDF in a separate "build" directory, so this
# will not actually create or overwrite the target $(PN).pdf. If there
# is no $(PN).pdf already, then we rename the new one to $(PN).pdf and
# we are done. But if there was a previous version, then we compare
# the two (new & old) versions. In either case, we rename the new one
# over the old one. But if the two differ, then we repeat this process
# in a loop until the just-built PDF is identical to the one from the
# previous iteration.
$(PN).pdf: $(SRCS) $(BUILDDIR)/$(PN).bbl $(INDEX_DSTS)
	$(LATEX) $(PN).tex

	if [ -f $@ ]; then \
	  while ! cmp -s $@ $(BUILDDIR)/$(PN).pdf; do \
	    mv $(BUILDDIR)/$(PN).pdf $@; \
	    $(LATEX) $(PN).tex; \
	  done; \
	fi;

	if grep -q 'Rerun to get' $(BUILDDIR)/$(PN).log; then \
	  $(LATEX) $(PN).tex; \
	fi;

	mv $(BUILDDIR)/$(PN).pdf $@

$(BUILDDIR):
	mkdir $@

$(BUILDDIR)/$(PN).aux: $(SRCS) | $(BUILDDIR)
	$(LATEX) $(PN).tex


# The pipe below indicates an "order-only dependency" on the aux file.
# Without it, every compilation of $(PN).tex would produce a new
# $(PN).aux, and thus $(PN).bbl would be rebuilt. This in turn causes
# $(PN).pdf to appear out-of-date, which leads to a recompilation of
# $(PN).tex... and so on. The order-only dependency means we won't
# rebuild $(PN).bbl if $(PN).aux changes.
#
# As a side effect, we now need to depend on $(SRCS) here, since we
# won't pick it up transitively from $(PN).aux.
#
# If the $BIBS variable is undefined, we presume that there are no
# references and create an empty bbl file. Otherwise, we risk trying
# to run biblatex on an aux file containing no citations. If you do
# define $BIBS but don't cite anything, you'll run into a similar
# problem. Don't do that.
#
$(BUILDDIR)/$(PN).bbl: $(SRCS) | $(BUILDDIR)/$(PN).aux
ifdef BIBS
	bibtex $(BUILDDIR)/$(PN).aux
else
	printf '' > $@
endif

# If the output PDF exists but the log file does not, then an attempt
# to "build the log file" (i.e. build the PDF) would do nothing. Thus
# whenever the log file does not exist, we do a fresh build.
$(BUILDDIR)/$(PN).log: $(SRCS)
	$(MAKE) clean
	$(MAKE)

# Ensure that there are no overfull or underfull boxes in the output
# document by parsing the log for said warnings.
.PHONY: check-boxes
check-boxes: $(BUILDDIR)/$(PN).log
	@! grep -i 'overfull\|underfull' $<

# Run chktex to find silly mistakes. There is some exit code weirdness
# (Savannah bug 53129), so we just look for empty output.
.PHONY: check-chktex
CHKTEX = chktex --localrc .chktexrc --quiet --inputfiles=0
check-chktex:
	@chktexout=$$($(CHKTEX) $(PN).tex); \
	  test -z "$${chktexout}" || { echo "$${chktexout}" 1>&2; exit 1; }

# Ensure that there are no undefined references in the document by
# parsing the log file for said warnings.
.PHONY: check-undefined
check-undefined: $(BUILDDIR)/$(PN).log
	@! grep -i 'undefined' $<

# Run python doctests
.PHONY: check-python
check-python:
	cd tests && python -m doctest ./*.py && cd ../

# Run a suite of checks.
.PHONY: check
check: check-boxes check-chktex check-undefined check-python

# Clean up leftover junk. This only looks overcomplicated because
# the *.{foo,bar} syntax supported by Bash is not POSIX, and Make
# will execute these commands using /bin/sh (which should be POSIX).
JUNK_EXTENSIONS  = aux bbl bcf blg glo ilg ist listing lof log nav out pdf
JUNK_EXTENSIONS += snm spl toc xml
.PHONY: clean
clean:
	for ext in $(JUNK_EXTENSIONS); do rm -f *.$$ext; done;
	rm -rf dist/ $(BUILDDIR)/
	rm -f $(SAGE_LISTING_DSTS) $(INDEX_SRCS) $(INDEX_DSTS)

# If this document will be published, the publisher isn't going to
# have your BibTeX database or your mjotex files. So, you need to
# package them up along with the code for your document. This target
# will create a "dist" directory and copy the necessary stuff there.
#
.PHONY: dist
dist: $(BUILDDIR)/$(PN).bbl
	mkdir -p dist
	cp $(SRCS) dist/

