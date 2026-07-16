FULL_DIAG ?= OFF

DIAG_SUFFIX := $(if $(filter ON,$(FULL_DIAG)),_ed,)

all: PRESET = default$(DIAG_SUFFIX)
cascadelake: PRESET = cascadelake$(DIAG_SUFFIX)
icelake: PRESET = icelake$(DIAG_SUFFIX)
debug: PRESET = debug$(DIAG_SUFFIX)

all cascadelake icelake debug:
	@cmake --preset $(PRESET)
	@cmake --build --preset $(PRESET)

clean:
	@rm -rf build
	@rm -rf auto_generated*

.PHONY: all clean icelake cascadelake debug