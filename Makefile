BUILD_FILE_NAME=build/board

all: create_proj


create_proj: $(BUILD_FILE_NAME)

$(BUILD_FILE_NAME):
	echo create build exists file
	mkdir -p build
	echo $PWD
	vivado_hls -f build_hls.tcl
	# vivado -nolog -nojournal -mode batch -source ./tcl/build.tcl

clean:
	rm -rf build

.PHONY: all create_proj clean