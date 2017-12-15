import platform
import os

from libconfig import *

# Reset before starting is good, it avoids issues when working with jupyter and %autoreload.
# Otherwise an error will be raised trying to register again already registered options.
reset_options()

# Register the options for the library.
# Options to run RosettaScripts
register_option("system", "overwrite",  False, "bool", "Allow overwriting already existing files")
register_option("system", "output", "./", "path_out", "Default folder to output generated files")

register_option("rosetta", "path",  os.path.expanduser('~'), "path_in", "Path to the rosetta binaries")
if platform.linux_distribution()[0] != "":
    register_option("rosetta", "compilation", "linuxgccrelease", "string", "Target binaries of rosetta")
elif platform.mac_ver()[0] != "":
    register_option("rosetta", "compilation", "macosclangrelease", "string", "Target binaries of rosetta")
else:
    register_option("rosetta", "compilation", "winccrelease", "string", "Target binaries of rosetta")

# Let's assume one wants to give the option to generate a user's configuration file.
# First time the library is called, it can create this config file with the current defaults.
# The user can change those for further runs. Ideally, the file would be something such as:
config_file = os.path.join(os.getenv("HOME"), ".rstoolbox.cfg")

# Either make or read from the file.
if not os.path.isfile(config_file):
    write_options_to_YAML( config_file )
else:
    set_options_from_YAML( config_file )

# Finally, register_option and reset_option are taken out from the global view so that they are not
# imported with the rest of the functions. This way the user can not access to them when importing
# the library and has to work through the rest of the available functions.
for name in ["register_option", "reset_options"]:
    del globals()[name]
