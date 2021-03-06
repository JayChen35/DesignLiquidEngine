# DesignLiquidEngine

## Quick Start

To run:

```bash
pip3 install -r requirements.txt
```

```bash
python main.py
```

### Usage Options

Purely running `python main.py` takes the user through a series of user keyboard inputs. If this is undesirable, there exists
ways to execute the script using command line arguments.

To run using default settings:

```bash
python main.py --default
```

or

```bash
python main.py -d
```

By default, the configuration file is `config.yaml` and output file names indicate the date and time in which they were created (e.g. `05-Mar-2021_212414EST`).

To run while specifying the configuration file path:

```bash
python main.py -c config.yaml
```

Here, `-c` obviously abbreviates "configuration." The user may replace `config.yaml` in the example above with the path to his/her desired configuration file. Note that DesignLiquidEngine only accepts `.yaml` configuration files. `.json` file support is a work-in-progress; YAML files were chosen so that the user can see comments next to a parameter, which describes that parameter and its unit.

To run while specifying the desired name for all output files:

```bash
python main.py -n project_caelus
```

Here, `-n` obviously abbreviates "name." Any name that can be used to name files and directories in your operating system are valid here, and DesignLiquidEngine will correct any invalid names.

### Customizing the Configuration File (`config.yaml`)

Allowable data types: `ints`, `floats`, `null`.
