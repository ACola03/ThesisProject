## This is Adam's thesisProject repo

current: target
-include target.mk
Ignore = target.mk

vim_session:
	bash -cl "vmt README.md"

######################################################################

Sources += README.md TODO.md

######################################################################

Ignore += *.html
Sources += *.R

autoknit = defined
autopipeR = defined
TrianglePlots.html: TrianglePlots.Rmd

## Conservative.Rmd
## ProbabilityIntegralTransform.Rmd
TestingFunctions.html: TestingFunctions.Rmd

poisson.Rout: poisson.R

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
