# not a GNU package. You can remove this line, if
# have all needed files, that a GNU package needs
AUTOMAKE_OPTIONS = foreign 1.4

SUBDIRS = src tests 

# test rule
TSTDIR = tests

check:	all
	@echo "-------------------------------------------------------------------"
	@echo "Running tests"
	@make test -s -C $(TSTDIR)
	@echo "-------------------------------------------------------------------"

# test gmm rule
TSTDIRGMM = tests/gmm

checkgmm:	all
	@echo "-------------------------------------------------------------------"
	@echo "Running gmm tests"
	@make test -s -C $(TSTDIRGMM)
	@echo "-------------------------------------------------------------------"
