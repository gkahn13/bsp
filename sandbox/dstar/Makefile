BUILD = build

FLAGS = -optimize
DLL = $(BUILD)/ClassLibrary1.dll

# timesteps
T = 15

# types of sym functions to create
SYM_FUNCS=Cost CostGradDiagHess

####### POINT ############
PT_DIR = point
PT_SYM_DIR = $(PT_DIR)/sym
PT_MASKS_DIR = $(PT_DIR)/masks

# Point.cs
PT_CS = $(PT_DIR)/Point.cs

# make point-state
PT_STATE_CS = $(PT_DIR)/PointState.cs
PT_STATE_EXEC = $(BUILD)/PointState.exe
point-state: $(PT_STATE_CS) $(PT_CS)
	make clean
	mono-csc $(FLAGS) -r:$(DLL) -out:$(PT_STATE_EXEC) $^
	bash modifyCompiledSource.sh 'state' $(PT_STATE_EXEC) $(T) $(PT_SYM_DIR) $(PT_MASKS_DIR) '$(SYM_FUNCS)'
	rm -f $(PT_STATE_EXEC)
	
# make point-control
PT_CONTROL_CS = $(PT_DIR)/PointControl.cs
PT_CONTROL_EXEC = $(BUILD)/PointControl.exe
point-control: $(PT_CONTROL_CS) $(PT_CS)
	make clean
	mono-csc $(FLAGS) -r:$(DLL) -out:$(PT_CONTROL_EXEC) $^
	bash modifyCompiledSource.sh 'control' $(PT_CONTROL_EXEC) $(T) $(PT_SYM_DIR) $(PT_MASKS_DIR) '$(SYM_FUNCS)'
	rm -f $(PT_CONTROL_EXEC)
	
	
####### arm ##################
ARM_DIR = ARM
ARM_SYM_DIR = $(ARM_DIR)/sym
ARM_MASKS_DIR = $(ARM_DIR)/masks

# Arm.cs
# PT_CS = $(PT_DIR)/Point.cs

# make point-state
ARM_STATE_CS = $(ARM_DIR)/ArmState.cs
ARM_STATE_EXEC = $(BUILD)/ArmState.exe
arm-state: $(ARM_STATE_CS)
	make clean
	mono-csc $(FLAGS) -r:$(DLL) -out:$(ARM_STATE_EXEC) $^
	bash modifyCompiledSource.sh 'state' $(ARM_STATE_EXEC) $(T) $(ARM_SYM_DIR) $(ARM_MASKS_DIR) '$(SYM_FUNCS)'
	rm -f $(ARM_STATE_EXEC)
	
clean:
	rm -f *.txt
	rm -f *.c
	rm -f $(BUILD)/*.exe
	rm -f tmp*