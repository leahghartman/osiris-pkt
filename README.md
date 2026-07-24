![OSIRIS dev logo](https://github.com/GoLP-IST/osiris/wiki/files/osiris-dev-logo.png)

# OSIRIS

[![latest version](https://img.shields.io/badge/osiris-v4.4.4-blue.svg?style=flat-square)](https://github.com/GoLP-IST/osiris/releases)

OSIRIS is a fully relativistic, massively parallel particle-in-cell (PIC) code
developed and maintained by the OSIRIS consortium. The consortium consists of
the [Extreme Plasma Physics Group](http://epp.ist.utl.pt/) and the [UCLA Plasma
Simulation Group](https://plasmasim.physics.ucla.edu/).

## Notes for using OSIRIS

If you want to use OSIRIS to perform large scale simulations, we recommend to
use latest stable version. You can obtain a specific version by downloading it
[here](https://github.com/GoLP-IST/osiris/releases) or using a more
sophisticate way via git. If you run

```shell
# clones OSIRIS repository into the folder osiris
git clone https://github.com/GoLP-IST/osiris.git
# move to the osiris repository
cd osiris
```

then you will have the latest stable version of OSIRIS. If you wish to update
OSIRIS to include the latest changes, run the command

```shell
git checkout master && git pull
```

to update the local master branch which includes the latest stable release.
**Note, this approach is design for users who just want to run simulations using
OSIRIS.** If you wish to contribute to OSIRIS or use the development version,
please read a guide on [contributing to OSIRIS](https://github.com/GoLP-IST/osiris/wiki).

An guide on running simulations and configuring your inputdeck can be found
[here](https://osirisdoc.wimpzilla.ist.utl.pt/).

