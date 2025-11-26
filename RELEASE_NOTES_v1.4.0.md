# CoastalME v1.4.0 Release Notes

**Release Date:** November 26, 2025

## Overview

CoastalME v1.4.0 is a significant release that introduces YAML configuration file support, providing a modern alternative to the traditional INI format. This release maintains backward compatibility with existing INI files while offering enhanced configuration capabilities through YAML.

## Major Features

### YAML Configuration Support
- **New Default Format**: CoastalME now accepts `cme.yaml` as the default input configuration file
- **Backward Compatible**: Falls back to `cme.ini` if YAML file is not found
- **Full Feature Migration**: All input features have been migrated to the YAML parsing workflow
- **Time-Varying Wave Data**: Complete support for time-varying wave data in YAML format

### Key Benefits of YAML Format
- More human-readable configuration files
- Better structure for complex configurations
- Easier to maintain and version control
- Support for comments and documentation within config files

## Improvements and Bug Fixes

### Configuration Parsing
- Fixed YAML cliff collapse flag reading
- Fixed reading of time-variable wave data in YAML input
- Full migration of input features to YAML parsing workflow

### Model Improvements
- Depth of closure improvements
- CShore bug fixes
- Bug fixes related to talus handling

### Performance Optimizations
- Fixed multiple instances of passing strings by value instead of by reference
- Enhanced debug output capabilities

### Code Quality
- Extensive code tidying and refactoring
- Improved code maintainability

## Installation and Usage

### Building from Source
CoastalME builds easily using Linux. For Windows users, we recommend using Windows Subsystem for Linux (WSL).

1. Clone the repository:
   ```bash
   git clone -b feature/ini-to-yaml https://github.com/CoastalModellingEnvironment/coastalme.git
   ```

2. Navigate to the src folder:
   ```bash
   cd coastalme/src
   ```

3. Run the build script:
   ```bash
   chmod a+x run_cmake.sh
   ./run_cmake.sh
   ```

4. Install:
   ```bash
   make install
   ```

### Using YAML Configuration
Create a `cme.yaml` file in the CoastalME root directory. The software will automatically use YAML format if available, or fall back to `cme.ini` if not found.

Example YAML configuration files are provided in the repository:
- `cme.yaml` - Standard configuration
- `cme.yaml.CLIFF` - Cliff-specific configuration
- `cme.yaml.DUNE` - Dune-specific configuration
- `cme.yaml.THORPNESS` - Thorpness case study configuration

## Migration Guide

### From INI to YAML
If you have existing `cme.ini` files, you can continue using them without modification. However, we encourage migrating to YAML for the benefits mentioned above.

To migrate:
1. Use the provided `cme.yaml` examples as templates
2. Transfer your existing configuration parameters to the YAML structure
3. Take advantage of YAML's improved readability and commenting features

## Compatibility

- **Backward Compatible**: Existing INI files continue to work
- **Recommended**: GDAL >= 2.1
- **Required**: CMake >= 3.16
- **Required**: C++11 or later
- **Required**: Fortran 2008 compiler

## Known Issues

- This is a development/testing release from the `feature/ini-to-yaml` branch
- Some experimental features may require additional testing
- Please report any issues via the [GitHub issue tracker](https://github.com/coastalme/coastalme/issues)

## Future Development

Version 1.4.0 establishes the foundation for YAML-based configuration. Future releases will:
- Expand YAML configuration options
- Add more comprehensive validation
- Provide migration tools for complex configurations

## Credits

- **Lead Developers**: Andres Payo (British Geological Survey), David Favis-Mortlock (British Geological Survey)
- **Contributors**: See [COMMITERS.md](COMMITERS.md) for full list

## References

- Main site: https://www.osgeo.org/projects/coastalme/
- Wiki: https://earthwise.bgs.ac.uk/index.php/Category:Coastal_Modeling_Environment
- Repository: https://github.com/coastalme/coastalme
- Documentation: https://coastalme.github.io/coastalme/

## Citation

Payo, A., Favis-Mortlock, D., Dickson, M., Hall, J. W., Hurst, M. D., Walkden, M. J. A., Townend, I., Ives, M. C., Nicholls, R. J., and Ellis, M. A.: Coastal Modelling Environment version 1.0: a framework for integrating landform-specific component models in order to simulate decadal to centennial morphological changes on complex coasts, Geosci. Model Dev., 10, 2715–2740, https://doi.org/10.5194/gmd-10-2715-2017, 2017.

## License

CoastalME is free software distributed under the GNU General Public License v3.0. See [LICENSE.md](LICENSE.md) for details.
