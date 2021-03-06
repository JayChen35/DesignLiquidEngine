# DesignLiquidEngine

## Quick Start

To run:

```properties
pip3 install -r requirements.txt
```

```properties
python main.py
```

### Usage Options

Purely running `python main.py` takes the user through a series of user keyboard inputs. If this is undesirable, there exists
ways to execute the script using command line arguments.

To run using default settings:

```properties
python main.py --default
```

or

```properties
python main.py -d
```

By default, the configuration file is `config.yaml` and output file names indicate the date and time in which they were created (e.g. `05-Mar-2021_212414EST`).

To run while specifying the configuration file path:

```properties
python main.py -c config.yaml
```

Here, `-c` obviously abbreviates "configuration." The user may replace `config.yaml` in the example above with the path to his/her desired configuration file. Note that DesignLiquidEngine only accepts `.yaml` configuration files. `.json` file support is a work-in-progress; YAML files were chosen so that the user can see comments next to a parameter, which describes that parameter and its unit.

To run while specifying the desired name for all output files:

```properties
python main.py -n project_caelus
```

Here, `-n` obviously abbreviates "name." Any name that can be used to name files and directories in your operating system are valid here, and DesignLiquidEngine will correct any invalid names.

### Customizing the Configuration File (`config.yaml`)

Allowable data types: `ints`, `floats`, `null`.

## Documentation

The following section details the theory behind DesignLiquidEngine, as well as further usage information for any subdirectories.

### `./helpers/PropSimEX`

Contains the PropSim integrated executable source code. Runs a liquid rocket engine simulator that is optimized via the `mcc` MATLAB Compiler.

#### MATLAB Compiler

To generate an executable (.exe file for Windows and .sh file for macOS/Linux) from this source code, run

```properties
mcc -m PropSimIntegrated.m
```

in the MATALB command window. Once the compilation is complete, a `PropSimIntegrated.exe` or `PropSimIntegrated.sh` is generated, which can be run via a command prompt. These executables are also called from DesignLiquidEngine.

#### Outputs

`PropSimIntegrated` outputs a MATLAB data file called `PropSimOutput.mat` in the `./helpers` directory. A copy is made and moved to the current case folder (under `./helpers/case-files/[your-case-name]`). There, the data file is read by DesignLiquidEngine and plotted given that the runtime option (specified in `./config.yaml`) is set to `true`.
