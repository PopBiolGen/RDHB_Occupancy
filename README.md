# RDHB_Occupancy

Analyses of occupancy of Red Dwarf Honeybees (RDHB).

## Data access

This repo contains code for analysing data from the RDHB response.  The data are not publicly available and are not included in the repo.  For those with access to the data, you will need to set up a new variable in your .Renviron, with something like:

`usethis::edit_r_environ()`

And then define a new variable in this file,

`DATA_PATH="your/local/path/to/data/directory"`.

If you then restart R, this variable will be available, and scripts like `a-setup.R` will now know where to find the data.

## Repo structure

It is currently a bit organic.  All scripts are in the /src folder.  Model descriptions are in the /ms folder.

Scripts a-f fit occupancy models of increasing complexity to the data.  Script f is probably the most useful of these.

Scripts g-h develop, test, and fit a complicated dynamic occupancy model to the data.  The intent of this model was to attempt to thread control effort into the occupancy framework.  It kind of works, but needs more development.

Scripts i-j switches gears entirely and fits a model designed to identify the location of colonies, given (imperfect) detection of foraging bees (and data on a foraging kernel).  This model is given the shorthand name, "point-process-model".






