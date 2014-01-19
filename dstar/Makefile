BUILD = build

FLAGS = -optimize
DLL = $(BUILD)/ClassLibrary1.dll

# timesteps
T = 15

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
	mono-csc $(FLAGS) -r:$(DLL) -out:$(PT_STATE_EXEC) $^
	mono $(PT_STATE_EXEC) $(T)
	mv *.c $(PT_SYM_DIR)/state-symeval-$(T).c
	mv *.txt $(PT_MASKS_DIR)/state-masks-$(T).txt
	rm -f $(PT_STATE_EXEC)
	
# make point-control
PT_CONTROL_CS = $(PT_DIR)/PointControl.cs
PT_CONTROL_EXEC = $(BUILD)/PointControl.exe
point-control: $(PT_CONTROL_CS) $(PT_CS)
	mono-csc $(FLAGS) -r:$(DLL) -out:$(PT_CONTROL_EXEC) $^
	mono $(PT_CONTROL_EXEC) $(T)
	mv *.c $(PT_SYM_DIR)/control-symeval-$(T).c
	mv *.txt $(PT_MASKS_DIR)/control-masks-$(T).txt
	rm -f $(PT_CONTROL_EXEC)
	
clean:
	rm -f *.txt
	rm -f *.c
	rm -f $(BUILD)/*.exe