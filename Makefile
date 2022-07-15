BUILD_FILE_NAME=fft_accel
all: create_proj


create_proj: $(BUILD_FILE_NAME)

$(BUILD_FILE_NAME):
	@echo "Create build folder"
	mkdir -p build
	@echo $(PWD)
	@echo "Copy files"
	cp -f $(PWD)/tcl/run_hls.tcl $(PWD)/build
	@echo "go to build dir"
	cd $(PWD)/build	&&	vitis_hls -f run_hls.tcl	&&	vitis_hls -p $(BUILD_FILE_NAME)

clean:
	rm -rf build/ sim_files/

.PHONY: all create_proj clean