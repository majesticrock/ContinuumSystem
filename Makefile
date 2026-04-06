FULL_DIAG ?= OFF

BUILD_SUFFIX := $(if $(filter ON,$(FULL_DIAG)),_ed,)

BUILD_DIR               := build$(BUILD_SUFFIX)
CASCADELAKE_BUILD_DIR   := build_CascadeLake$(BUILD_SUFFIX)
ICELAKE_BUILD_DIR       := build_IceLake$(BUILD_SUFFIX)
DEBUG_BUILD_DIR         := build_debug$(BUILD_SUFFIX)

# Default target to build the project
all: $(BUILD_DIR)/Makefile
	@$(MAKE) -C $(BUILD_DIR)

$(BUILD_DIR)/Makefile: CMakeLists.txt
	@mkdir -p $(BUILD_DIR)
	@cd $(BUILD_DIR) && cmake \
		-DCMAKE_C_COMPILER=gcc \
		-DCMAKE_CXX_COMPILER=g++ \
		-DCONTINUUM_FULL_DIAG=$(FULL_DIAG) \
		..

cascadelake: $(CASCADELAKE_BUILD_DIR)/Makefile
	@$(MAKE) -C $(CASCADELAKE_BUILD_DIR)

$(CASCADELAKE_BUILD_DIR)/Makefile: CMakeLists.txt
	@mkdir -p $(CASCADELAKE_BUILD_DIR)
	@cd $(CASCADELAKE_BUILD_DIR) && cmake \
		-DCMAKE_C_COMPILER=gcc \
		-DCMAKE_CXX_COMPILER=g++ \
 		-DCLUSTER_BUILD=cascadelake \
		-DCONTINUUM_FULL_DIAG=$(FULL_DIAG) \
		..

icelake: $(ICELAKE_BUILD_DIR)/Makefile
	@$(MAKE) -C $(ICELAKE_BUILD_DIR)

$(ICELAKE_BUILD_DIR)/Makefile: CMakeLists.txt
	@mkdir -p $(ICELAKE_BUILD_DIR)
	@cd $(ICELAKE_BUILD_DIR) && cmake \
		-DCMAKE_C_COMPILER=gcc \
		-DCMAKE_CXX_COMPILER=g++ \
 		-DCLUSTER_BUILD=icelake \
		-DCONTINUUM_FULL_DIAG=$(FULL_DIAG) \
		..

debug: $(DEBUG_BUILD_DIR)/Makefile
	@$(MAKE) -C $(DEBUG_BUILD_DIR)

$(DEBUG_BUILD_DIR)/Makefile: CMakeLists.txt
	@mkdir -p $(DEBUG_BUILD_DIR)
	@cd $(DEBUG_BUILD_DIR) && cmake \
		-DCMAKE_C_COMPILER=gcc \
		-DCMAKE_CXX_COMPILER=g++ \
 		-DCMAKE_BUILD_TYPE=Debug \
		-DCONTINUUM_FULL_DIAG=$(FULL_DIAG) \
		..

clean:
	@rm -rf $(BUILD_DIR) $(CASCADELAKE_BUILD_DIR) $(ICELAKE_BUILD_DIR) $(DEBUG_BUILD_DIR) build_header
	@rm -rf $(BUILD_DIR)_ed $(CASCADELAKE_BUILD_DIR)_ed $(ICELAKE_BUILD_DIR)_ed $(DEBUG_BUILD_DIR)_ed
	@rm -rf auto_generated*

.PHONY: all clean icelake cascadelake debug