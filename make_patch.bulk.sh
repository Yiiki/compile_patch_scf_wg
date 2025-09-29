#!/usr/bin/bash
module load intel
mpiifort -o patch_bulk_chi.x mod_mixing.f90 mch_pulay.f mch_kerk.f patch_bulk_chi.f  getpot2.f90  gaussj.f  get_ALI.f cfft.f90 cfftd.f90 UxcCA.f90
#cp patch_scf_list2.x ../bin/
