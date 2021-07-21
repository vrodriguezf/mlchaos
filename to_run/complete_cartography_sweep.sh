#!/bin/bash
for i in $(seq 2 10); do
    papermill -p config.reduction_factor $i /home/victor/work/nbs_experiments/complete_cartography.ipynb /home/victor/pm_output.ipynb
done