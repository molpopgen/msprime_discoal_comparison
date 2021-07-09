#!/bin/bash

bash make_sim_commands.sh | parallel
bash extractS.sh
R --no-save < compare.R
