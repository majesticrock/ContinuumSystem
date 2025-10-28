BUILD_DIR = build
CASCADELAKE_BUILD_DIR = build_CascadeLake
ICELAKE_BUILD_DIR = build_IceLake
DEBUG_BUILD_DIR = build_debug

# Default target to build the project
all: $(BUILD_DIR)/Makefile
	@$(MAKE) -C $(BUILD_DIR)

$(BUILD_DIR)/Makefile: CMakeLists.txt
	@mkdir -p $(BUILD_DIR)
	@cd $(BUILD_DIR) && cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ ..

cascadelake: $(CASCADELAKE_BUILD_DIR)/Makefile
	@$(MAKE) -C $(CASCADELAKE_BUILD_DIR)

$(CASCADELAKE_BUILD_DIR)/Makefile: CMakeLists.txt
	@mkdir -p $(CASCADELAKE_BUILD_DIR)
	@cd $(CASCADELAKE_BUILD_DIR) && cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCLUSTER_BUILD=cascadelake ..

icelake: $(ICELAKE_BUILD_DIR)/Makefile
	@$(MAKE) -C $(ICELAKE_BUILD_DIR)

$(ICELAKE_BUILD_DIR)/Makefile: CMakeLists.txt
	@mkdir -p $(ICELAKE_BUILD_DIR)
	@cd $(ICELAKE_BUILD_DIR) && cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCLUSTER_BUILD=icelake ..

debug: $(DEBUG_BUILD_DIR)/Makefile
	@$(MAKE) -C $(DEBUG_BUILD_DIR)

$(DEBUG_BUILD_DIR)/Makefile: CMakeLists.txt
	@mkdir -p $(DEBUG_BUILD_DIR)
	@cd $(DEBUG_BUILD_DIR) && cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_BUILD_TYPE=Debug ..

clean:
	@rm -rf $(BUILD_DIR) $(CASCADELAKE_BUILD_DIR) $(ICELAKE_BUILD_DIR) $(DEBUG_BUILD_DIR) $(NDEBUG_BUILD_DIR) $(NO_MPI_BUILD_DIR) build_header
	@rm -rf auto_generated*

.PHONY: all clean icelake cascadelake debug