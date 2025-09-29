#!/usr/bin/bash
module load intel
ifort -o patch_bulk_eps.x patch_bulk_eps.f get_ALI.f gaussj.f
