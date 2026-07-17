FULL_DIAG ?= OFF

DIAG_SUFFIX := $(if $(filter ON,$(FULL_DIAG)),_ed,)

all: PRESET = default$(DIAG_SUFFIX)
cascadelake: PRESET = cascadelake$(DIAG_SUFFIX)
icelake: PRESET = icelake$(DIAG_SUFFIX)
debug: PRESET = debug$(DIAG_SUFFIX)
test: PRESET = test

all cascadelake icelake debug:
	@cmake --preset $(PRESET)
	@cmake --build --preset $(PRESET)

test:
	@cmake --preset $(PRESET)
	@cmake --build --preset $(PRESET)
	@ctest --preset $(PRESET)

clean:
	@rm -rf build
	@rm -rf auto_generated*

.PHONY: all clean icelake cascadelake debug test