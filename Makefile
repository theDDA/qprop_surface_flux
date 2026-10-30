.DEFAULT_GOAL := all

QPROP_SRC ?= ../qprop/qprop_v3.2/src
BUILD ?= build
CASES := attoclock long-wavelength vortex lissajous
GSL_CFLAGS ?= $(shell gsl-config --cflags)
GSL_LIBDIR ?= $(shell gsl-config --prefix)/lib
FINUFFT_CFLAGS ?=
FINUFFT_LIBS ?= -lfinufft

ifeq ($(shell uname -s),Darwin)
NATIVE ?= -mcpu=native
LINALG_LIBS ?= -framework Accelerate
else
NATIVE ?= -march=native
LINALG_LIBS ?= -llapacke -lopenblas
endif
CXXFLAGS ?= -O3 -ffast-math $(NATIVE) -std=c++11 -Wno-deprecated
ifneq ($(strip $(FINUFFT_LIBS)),)
A1_FLAGS := -DHAVE_FINUFFT $(FINUFFT_CFLAGS)
endif
LDLIBS := -L$(GSL_LIBDIR) -lgsl $(LINALG_LIBS) -lm
OVERLAY := $(BUILD)/overlay
INCLUDES := -Isrc -I$(OVERLAY)/base -I$(QPROP_SRC)/base $(GSL_CFLAGS)
BASE := bar ylm wavefunction hamop grid fluid factorial parameter gfunc
OBJECTS := $(addprefix $(BUILD)/obj/,$(addsuffix .o,$(BASE)))
LIB := $(BUILD)/libqprop.a
PATCHED := base/tsurffSpectrum.hh $(addsuffix /isurfv.cc,$(CASES))
TARGETS := $(foreach c,$(CASES),$(BUILD)/bin/$(c)/tsurff $(BUILD)/bin/$(c)/isurfv)
B1_TARGETS := $(foreach c,$(CASES),$(BUILD)/bin/$(c)/tsurff_b1)

.PHONY: all b1
all: $(TARGETS)
b1: $(B1_TARGETS)

$(OVERLAY)/.stamp: qprop-3.2.patch $(addprefix $(QPROP_SRC)/,$(PATCHED))
	@mkdir -p $(OVERLAY)/base $(addprefix $(OVERLAY)/,$(CASES))
	@for file in $(PATCHED); do cp "$(QPROP_SRC)/$$file" "$(OVERLAY)/$$file" || exit 1; done
	patch --batch --fuzz=0 -p1 -d $(OVERLAY) < qprop-3.2.patch
	@touch $@

$(BUILD)/obj/%.o: $(QPROP_SRC)/base/%.cc
	@mkdir -p $(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(GSL_CFLAGS) -I$(QPROP_SRC)/base -MMD -MP -c $< -o $@

$(LIB): $(OBJECTS)
	$(AR) rcs $@ $^

define case_rules
$(BUILD)/bin/$(1)/tsurff: $(QPROP_SRC)/$(1)/tsurff.cc $(LIB) $(OVERLAY)/.stamp $(wildcard src/*)
	@mkdir -p $$(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(A1_FLAGS) $(INCLUDES) -I$(QPROP_SRC)/$(1) -MMD -MP $$< $(LIB) $(LDLIBS) $(FINUFFT_LIBS) -o $$@
$(BUILD)/bin/$(1)/tsurff_b1: $(QPROP_SRC)/$(1)/tsurff.cc $(LIB) $(OVERLAY)/.stamp $(wildcard src/*)
	@mkdir -p $$(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -DQPROP_DIRECT_ANGULAR $(INCLUDES) -I$(QPROP_SRC)/$(1) -MMD -MP $$< $(LIB) $(LDLIBS) -o $$@
$(BUILD)/bin/$(1)/isurfv: $(QPROP_SRC)/$(1)/isurfv.cc $(LIB) $(OVERLAY)/.stamp $(wildcard src/*.hh)
	@mkdir -p $$(@D)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) -I$(QPROP_SRC)/$(1) -MMD -MP $(OVERLAY)/$(1)/isurfv.cc $(LIB) $(LDLIBS) -o $$@
endef
$(foreach c,$(CASES),$(eval $(call case_rules,$(c))))

-include $(OBJECTS:.o=.d) $(addsuffix .d,$(TARGETS) $(B1_TARGETS))
