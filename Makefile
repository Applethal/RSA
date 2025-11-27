GCC			= gcc

INPUT		= Instances/model_9.csv

SRC_DIR		= src
BUILD_DIR	= build

SRC			= $(wildcard $(SRC_DIR)/*.c)
OBJ			= $(SRC:$(SRC_DIR)/%.c=$(BUILD_DIR)/%.o)


TARGET = $(BUILD_DIR)/rsa

$(TARGET): $(OBJ)
	$(GCC) -g -Iinclude -o $@ $^

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	$(GCC) -g -Iinclude -c $< -o $@

run: $(TARGET)
	$(TARGET) $(INPUT)

debug: $(TARGET)
	$(TARGET) $(INPUT) -Debug

clean:
	rm -rf build/*

.PHONY: run debug
