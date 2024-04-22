# cochlea-nerve

Code of the model of mammalian auditory periphery by Tichacek et al. (2023)

## Usage


### Download
Please download a GitHub release or clone directly
```
git clone https://github.com/ondrejtichacek/cochlea-nerve
```

### Dependencies
Several packages from FileExchange are used, some of them are modified. The packages are shipped with the code within the `externals` directory.

### First time-config

The following commands are to be executed in MATLAB:

```
cd path_to_repo
START
ex1_1_outer_ear
```

-> this should result in an error:

```
Error using globalOpt
The system ondrej@aviendha is not configured. Add it to globalOpt to configure.
```

Please execute
```
edit globalOpt
```

Now add the options for your system, which will be used based on your user and computername.

