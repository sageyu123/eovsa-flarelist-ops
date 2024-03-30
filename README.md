# EOVSA Flarelist Operations (`eovsa-flarelist-ops`)

The `eovsa-flarelist-ops` repository is for the maintenance and data enrichment of the [EOVSA Flarelist website](http://www.ovsa.njit.edu/flarelist). It comprises essential scripts for updating the MySQL database with the latest solar flare observations by EOVSA and for generating spectrogram data and figures. This repository functions as a submodule within the broader [EOVSA Flarelist main repository](https://github.com/ovro-eovsa/eovsa-flarelist), ensuring the flare list remains up-to-date with recent observations.

## Prerequisites

To utilize the scripts in this repository effectively, you should have Python installed on your system along with several key dependencies. The exact requirements might differ based on your specific setup and the scripts you intend to use, but generally, you'll need:

- MySQL connectors for database operations.
- Data processing and visualization libraries, such as NumPy and Matplotlib.

**Environment Setup:** The scripts are currently deployed on the OVSA server, accessible under the `sjyu` account. To begin working with the scripts, log in to the server (`sjyu@ovsa`), and attach to the pre-existing screen session dedicated to flare list updates:

```bash
screen -rd flarelist
```

**Repository Update:** Ensure you're working with the latest version of the `eovsa-flarelist-ops` repository. It's synced with the main repository at [https://github.com/ovro-eovsa/eovsa-flarelist.git](https://github.com/ovro-eovsa/eovsa-flarelist.git). To update:

```bash
cd $HOME/eovsa-flarelist-ops
sudo git fetch origin
sudo git merge origin/main
```

## Usage

Manual updates to the MySQL database might be necessary at times. To access the MySQL console on the OVSA server with the relevant credentials, use:

```bash
mysql -h"${FLARE_DB_HOST}" -u"${FLARE_DB_USER}" -p"${FLARE_DB_PASSWORD}"
```

**Setting Permissions:** Before running the `run_flarelist2sql.sh` script, switch to the `user` group to ensure you have the appropriate file access permissions:

```bash
newgrp user
```

**Script Execution:** Execute the script with the following syntax, providing optional arguments as needed:

```bash
./run_flarelist2sql.sh -t "YYYY-MM-DD HH:MM:SS" "YYYY-MM-DD HH:MM:SS"
```

Options include:
- `-t` or `--timerange` to define the time range for fetching and processing flare data.
- `--do_manu` for manual adjustments of the start/end times for radio bursts.
- `-h` or `--help` to show this help message and exit.

## License

This project is under the [MIT License](LICENSE.md). For more details, refer to the LICENSE file.