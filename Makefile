## This is Adam's thesisProject repo

current: target
-include target.mk
Ignore = target.mk

vim_session:
	bash -ic "vmt README.md TODO.md notes.md"

######################################################################

Sources += README.md TODO.md
Sources += $(wildcard *md)
Ignore += *_files

######################################################################

## Naive exploration

## This may be how we should do in-outing of P values
binomTest.Rout: binomTest.R

## Meeting looks at naive location statistic
nlsMeet.Rout: nlsMeet.R

## AC calculates order statistics 2026 Jul 03 (Fri)
unifCpStats.Rout: unifCpStats.R trianglePlot_Poisson.R

######################################################################

Ignore += *.html
Sources += *.R

autoknit = defined
autopipeR = defined

TrianglePlots.html: TrianglePlots.Rmd

## Conservative.Rmd
## ProbabilityIntegralTransform.Rmd
TestingFunctions.html: TestingFunctions.Rmd

## PoissonTests.html: PoissonTests.Rmd

week3.html: week3.Rmd

poisson.Rout: poisson.R

## These are just functions
## trianglePlot_Poisson.Rout: trianglePlot_Poisson.R

jdTriangle.Rout: jdTriangle.R trianglePlot_Poisson.rda

## jdTest.Rout: jdTest.R trianglePlot_Poisson.R
jdTest.Rout: jdTest.R trianglePlot_Poisson.rda

### Makestuff

Sources += Makefile

Ignore += makestuff
msrepo = https://github.com/dushoff

## ln -s ../makestuff . ## Do this first if you want a linked makestuff
Makefile: makestuff/00.stamp
makestuff/%.stamp: | makestuff
	- $(RM) makestuff/*.stamp
	cd makestuff && $(MAKE) pull
	touch $@
makestuff:
	git clone --depth 1 $(msrepo)/makestuff

-include makestuff/os.mk

-include makestuff/pipeR.mk

-include makestuff/git.mk
-include makestuff/visual.mk
