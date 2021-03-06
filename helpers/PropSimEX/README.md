# PropSimEX

PropSim integrated executable source code. Runs a liquid rocket engine simulator.

## MATLAB Compiler

To generate an executable (.exe file for Windows and .sh file for macOS/Linux) from this source code, run

```matlab
mcc -m PropSimIntegrated.m
```

in the MATALB command window. Once the compilation is complete, a `PropSimIntegrated.exe` or `PropSimIntegrated.sh` is generated, which can be run via a command prompt. These executables are also called from DesignLiquidEngine.

## Outputs

`PropSimIntegrated` outputs a MATLAB data file called `PropSimOutput.mat` in the `./helpers` directory. A copy is made and moved to the current case folder (under `./helpers/case-files/[your-case-name]`). There, the data file is read by DesignLiquidEngine and plotted given that the runtime option (specified in `./config.yaml`) is set to `true`.
