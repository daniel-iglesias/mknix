# this is example-file: tests/Makefile.am

TESTSOURCE = $(wildcard test*.cpp)

TESTOBJS=$(TESTSOURCE:%.cpp=%.o)

TESTEXECS = $(TESTOBJS:%.o=%)

EXECS=$(TESTOBJS:.o=.exec)

DIFFS=$(TESTOBJS:.o=.diff)

TESTDEPENDFILE = .depend

# adding source and verified output files to the distribution
EXTRA_DIST = $(TESTSOURCE)

# additional linker and compiler flags
LDFLAGS = $(GLOBALLDFLAGS) -lblas -llapack
CPPFLAGS = $(GLOBALCPPFLAGS) -I$(top_srcdir)/src
CPPFLAG_LAPACK = -DHAVE_LAPACK

#################################################################

test: $(TESTEXECS) $(EXECS)

$(TESTEXECS): $(TESTDEPENDFILE) $(TESTOBJS) 
	@$(CXX) $@.o $(LDFLAGS) -o $@

%.exec: %
	@name=$<; \
	    echo "----------------------------------------" ;\
			printf " Testing %-40s \n" $$name; \
	    echo "----------------------------------------" ;\
	./$< ;\
	    echo "----------------------------------------"

%.diff: %
	@./$< > $<.lastout
	@name=$<; \
	name=$${name/test/}  ;\
	if [ -e $<.verified ] ; then \
		if diff -w $<.verified $<.lastout >$@ ; then \
			printf " Testing %-40s OK\n" $$name; \
		else \
			printf " Testing %-40s FAILED\n" $$name;\
		    echo "----------------------------------------" ;\
		    echo "Diff for $$name (verified '<', current '>')" ;\
	  		cat $@ ;\
			echo "----------------------------------------" ;\
			echo "[ If new output is correct: mv $<.lastout $<.verified ]" ;\
		fi ;\
		$(RM) $@ ;\
	else \
	    echo "----------------------------------------" ;\
	    echo "Output of $$name" ;\
		cat $<.lastout ;\
	    echo "----------------------------------------" ;\
		echo "WARNING: verified output of $$name is missing." ;\
	    echo "Please generate a file $<.verified that" ;\
	    echo "contains the verfied output of $$name."; \
		echo "[ If output is correct: mv $<.lastout $<.verified ]" ;\
	fi

$(TESTDEPENDFILE):
		@$(CXX) $(INCLUDES) -MM   $(TESTSOURCES) > $(TESTDEPENDFILE)

-include $(TESTDEPENDFILE)

%.o: %.cpp
	@name=$<; \
	name=$${name/test/}  ;\
	printf " Compiling %-38s OK\n" $$name; \
	$(CXX) $(INCLUDES) $(CPPFLAGS) $(CPPFLAG_LAPACK) -c $<
	
clean-local:
	$(RM) -rf $(TESTOBJS) $(TESTEXECS) $(TESTDEPENDFILE) $(TESTEXECS:=.lastout)\
		$(TESTEXECS:=.diff) *~ $(TESTEXECS:=.exe)
