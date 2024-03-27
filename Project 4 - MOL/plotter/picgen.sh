#!/bin/bash
mkdir ../report/figures
bash plot_heat_true.sh
bash plot_heat_lp.sh heat-a-1 -0.001 0.5 -0.001 0.5 -0.001 0.5
bash plot_heat_lp.sh heat-a-2 -0.2 0.4 -0.001 0.5 -0.001 0.5
bash plot_heat_lp.sh heat-a-3 -0.001 0.5 -0.001 0.5 -0.001 0.5
bash plot_heat_lp.sh heat-a-4 -0.001 0.5 -0.001 0.5 -0.001 0.5
bash plot_heat_lp.sh heat-b-1 -0.001 0.5 -0.001 0.5 -0.001 0.5
bash plot_heat_lp.sh heat-b-2 -0.001 0.5 -0.001 0.5 -0.001 0.5
bash plot_heat_lp.sh heat-c-1 -0.001 0.5 -0.001 0.5 -0.001 0.5
bash plot_heat_lp.sh heat-c-2 -1.5 1.5 -3 4 -10000 10000
bash plot_heat_lp.sh heat-d-1 -0.001 0.5 -0.001 0.5 -0.001 0.5
bash plot_heat_lp.sh heat-d-2 -0.6 0.4 -0.3 0.6 -0.3 0.6
bash plot_heat_l.sh heat-e-1
bash plot_advection_true.sh
bash plot_advection.sh advection-1
bash plot_advection.sh advection-2
bash plot_advection.sh advection-3
bash plot_advection.sh advection-4
bash plot_advection.sh advection-5
bash plot_advection.sh advection-6
bash plot_advection.sh advection-7
rm *.txt