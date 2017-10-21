import platform
import os

from libconfig import *

reset_options()

register_option("system", "overwrite",  False, "bool", "Allow overwriting already existing files")
register_option("system", "output", "./", "path_out", "Default folder to output generated files")

register_option("rosetta", "path",  os.path.expanduser('~'), "path_in", "Path to the rosetta binaries")
if platform.linux_distribution()[0] != "":
    register_option("rosetta", "compilation", "linuxgccrelease", "string", "Target binaries of rosetta")
elif platform.mac_ver()[0] != "":
    register_option("rosetta", "compilation", "macosclangrelease", "string", "Target binaries of rosetta")
else:
    register_option("rosetta", "compilation", "winccrelease", "string", "Target binaries of rosetta")
