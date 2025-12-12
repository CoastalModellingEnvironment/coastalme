<p align="left">
  <img src="https://www.osgeo.org/wp-content/uploads/preview_CoastalME_logo_white_740x412_acf_cropped.png" alt="CoastalME logo white background" width="250">
</p>

# CoastalME
The [latest release](https://github.com/coastalme/coastalme/releases) is version 1.4.0 (26 Nov 2025).

The [main branch][https://github.com/CoastalModellingEnvironment/coastalme] is also 1.4.0.

The [development branch][https://github.com/CoastalModellingEnvironment/coastalme/tree/develop] is version 1.4.1 (12 Dec 2025). It should work well, but use at your own risk.

There are other experimental versions which may not even compile. 
=======

By <a href="https://www.bgs.ac.uk/people/payo-garcia-andres/">Andres Payo</a>, <a href="https://en.wikipedia.org/wiki/David_Favis-Mortlock">David Favis-Mortlock</a>, <a href="https://www.linkedin.com/in/wilf-chun-a35a3a148"> Wilf Chun, and [others](COMMITERS.md).

## Table of contents

- [About](#about)
- [Quick Start](#quick-start)
- [References](#references)
- [How to build from source](#how-to-build-from-source)
- [How to contribute](#how-to-contribute)
- [Code of Conduct](#code-of-conduct)
- [Credits](#credits)

## About
CoastalME (Coastal Modelling Environment) is a Free and Open Source software for geospatial modelling to simulate decadal and longer coastal morphological changes.

It is an engineering tool for advanced modellers seeking to simulate the interaction of multiple coastal landforms and different types of human interventions.

[Payo et al. (2017)](https://doi.org/10.5194/gmd-10-2715-2017) described in detail the rationale behind CoastalME and demonstrated how it can be used to integrate; the Soft Cliff and Platform Erosion model [SCAPE](http://www.bioone.org/doi/abs/10.2112/JCOASTRES-D-10-00099.1), the Coastal Vector Evolution Model [COVE](http://onlinelibrary.wiley.com/doi/10.1002/2015JF003704/full) and the Cross Shore model [CSHORE](http://ascelibrary.org/doi/10.1061/(ASCE)WW.1943-5460.0000347).

The software is written in C++ using an object-oriented approach. Source code documentation using Doxygen is available [here](https://coastalme.github.io/coastalme/).

The C++ source code is available for [download](https://github.com/coastalme/coastalme) under [GNU](https://github.com/coastalme/CoastalME/tree/master?tab=GPL-3.0-1-ov-file) open source license.

* Main site: https://www.osgeo.org/projects/coastalme/ - Developer docs, links to other resources
* Wiki site: https://earthwise.bgs.ac.uk/index.php/Category:Coastal_Modeling_Environment - User docs, links to other resources
* GIT repository: https://github.com/coastalme/coastalme
* Bug tracker: https://github.com/coastalme/coastalme/issues
* [CoastalME Discussions forum](https://github.com/coastalme/coastalme/discussions/)
* [Which GitHub-discussion tool should I use?](https://docs.github.com/en/get-started/using-github/communicating-on-github#which-discussion-tool-should-i-use)
* [Doxygen documentation](https://coastalme.github.io/coastalme/)

## Quick Start

CoastalME builds easily using Linux. If you wish to run CoastalME on Windows, then we currently recommend using the Windows Subsystem Linux (WSL) software to do this.

Create a local copy of the github repository, for example by downloading a zipfile, then unpacking it or cloning. We suggest unpacking it to something like "/home/YOUR NAME/Projects/CoastalME/", this is then your CoastalME folder.

```
git clone https://github.com/coastalme/coastalme
```

In a terminal window (i.e. at a command-line prompt) move to the CoastalME folder. 

Then move to the the src folder Then run run_cmake.sh. 
```
cd CoastalME/src
./run_cmake.sh
```
If you get a "Permission denied" message `-bash: ./run_cmake.sh: Permission denied` you will have to grant permission using `chmod a+x run_cmake.sh`, `chmod a+x cshore/make_cshore.sh` and then `./run_cmake.sh`

If you see error messages about missing software (for example, telling you that CMake cannot be found or is too old, or GDAL cannot be found or is too old) then you need to install or update the software that is causing the problem.

Run make install `make install`. This will create an executable file called cme in the CoastalME folder.

Edit cme.ini to tell CoastalME which input file to read (for example, in/test_suite/minimal_wave_angle_230/minimal.dat). You'll need to make sure that the "Path for output" listed in cme.ini (for example, out/test_suite/minimal_wave_angle_230) exists. If it doesn't exist, then create it: `mkdir -p out/test_suite/minimal_wave_angle_230/`.

Run cme `./cme`. Output will appear in the "Path for output" folder.

To test that your installation is running correctly, you can run a suite of pre-defined tests by running the following commands:

	chmod a+x run_test_suite.sh
	./run_test_suite.sh

The `chmod` comand will ensure that you have permission to execute the run_test_suite.sh file.

Once you have CoastalME (CME) up and running, you can reduce the quantity of output (it can be overwhelming!) in several ways. 

* Change "Content of log file" in the main CME input file for any of the test suite runs (the name of this main input file is listed in cme.ini, both are simple text files). If you set "Content of log file" to zero, then CME won't output a log file; setting it to 4 (all output) is really only useful to developers.

* Change "GIS vector files to output" and "GIS vector files to output" in the main CME input file. These are both set to "all" in the test suite files on GitHub. Instead of "all" you can list the space-separated codes for only the output that you want to see. A list of CME GIS output codes is in codes.txt"

Enjoy!

## References
  
  See <a href="https://doi.org/10.5281/zenodo.1418810"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.1418810.svg" alt="DOI"></a> for the first release of the source code.

## How to cite

Payo, A., Favis-Mortlock, D., Dickson, M., Hall, J. W., Hurst, M. D., Walkden, M. J. A., Townend, I., Ives, M. C., Nicholls, R. J., and Ellis, M. A.: Coastal Modelling Environment version 1.0: a framework for integrating landform-specific component models in order to simulate decadal to centennial morphological changes on complex coasts, Geosci. Model Dev., 10, 2715–2740, https://doi.org/10.5194/gmd-10-2715-2017, 2017.
	


## How to build from source

See [BUILDING.md](BUILDING.md)


## How to contribute

See [CONTRIBUTING.md](CONTRIBUTING.md)

## Code of Conduct

See [CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md)

## Credits

Devised by [Andres Payo](https://github.com/coastalme) and [David Favis-Mortlock](https://github.com/davefavismortlock/). Lead developer and Version Control Administrator is [David Favis-Mortlock](https://github.com/davefavismortlock/) 
