# Compares the accuracy against MPFR and the runtime of the variadic
# AccurateNorm and AccurateRNorm implemented by several revisions.
#
#   make accuracy   # writes build/accuracy_<variant>.csv
#   make benchmark  # writes build/benchmark_<compiler>_<variant>_<round>.csv
#   make report     # prints the results as Markdown tables
#
# Requires GCC, Clang, MPFR, GMP, Python 3, and the revisions below to be
# present in the Git repository.

# Compensated accumulation of all squares, the base of both compositions
REVISION_baseline := 74fadb55
# Recursive composition of 2-argument norms
REVISION_recursive := 44268ea1
# Sequential composition of 2-argument norms
REVISION_sequential := 048d4d69

VARIANTS := baseline recursive sequential
COMPILERS := gcc clang
ROUNDS := 1 2

TOP := $(shell git rev-parse --show-toplevel)
BUILD := build
CXXFLAGS := -std=c++17 -O3 -march=native
CXX_gcc := g++
CXX_clang := clang++

ACCURACY := $(VARIANTS:%=$(BUILD)/accuracy_%.csv)
BENCHMARKS := $(foreach compiler,$(COMPILERS),\
                $(VARIANTS:%=$(BUILD)/benchmark_$(compiler)_%))

.PHONY: accuracy benchmark report clean
.PRECIOUS: $(BUILD)/%/include

accuracy: $(ACCURACY)

# Runs the benchmarks of all variants one after another and repeats the
# sequence to reduce the influence of other processes.
benchmark: $(BENCHMARKS)
	@for round in $(ROUNDS); do \
	  for compiler in $(COMPILERS); do \
	    for variant in $(VARIANTS); do \
	      $(BUILD)/benchmark_$${compiler}_$${variant} \
	        > $(BUILD)/benchmark_$${compiler}_$${variant}_$${round}.csv \
	        2> /dev/null; \
	    done; \
	  done; \
	done

report:
	@python3 report.py $(BUILD) $(VARIANTS)

$(BUILD)/%/include:
	@mkdir -p $(BUILD)/$*
	git -C $(TOP) archive $(REVISION_$*) include | tar -x -C $(BUILD)/$*
	@test -f $@/ceres/accurate_norm.h

$(BUILD)/accuracy_%: accuracy.cc | $(BUILD)/%/include
	$(CXX_gcc) $(CXXFLAGS) -I$(BUILD)/$*/include $< -o $@ -lmpfr -lgmp

$(BUILD)/accuracy_%.csv: $(BUILD)/accuracy_%
	./$< > $@

$(BUILD)/benchmark_gcc_%: benchmark.cc | $(BUILD)/%/include
	$(CXX_gcc) $(CXXFLAGS) -I$(BUILD)/$*/include $< -o $@

$(BUILD)/benchmark_clang_%: benchmark.cc | $(BUILD)/%/include
	$(CXX_clang) $(CXXFLAGS) -I$(BUILD)/$*/include $< -o $@

clean:
	rm -rf $(BUILD)
