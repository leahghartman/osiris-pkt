write(output_unit,'(A)') 'build_info = {'
write(output_unit,'(A)') '  "version": "",'
write(output_unit,'(A)') '  "branch": "",'
write(output_unit,'(A)') '  "dimensions": "1",'
write(output_unit,'(A)') '  "precision": "DOUBLE",'
write(output_unit,'(A)') '  "build_type": "debug",'
write(output_unit,'(A)') '  "system": "macosx-m1.gnu",'
write(output_unit,'(A)') '  "build_tool": "make",'
write(output_unit,'(A)') '  "build_flags": {'
write(output_unit,'(A)') '    "FPP"  :"/opt/homebrew//bin/mpicc -C -E -x as&
   &sembler-with-cpp -D__HAS_MPI_IN_PLACE__ -D_OPENMP -DHDF5",'
write(output_unit,'(A)') '    "FPPF" :"-DP_X_DIM=1 -DOS_REV=\"\" -DFORTRANS&
   &INGLEUNDERSCORE -DPRECISION_DOUBLE -DENABLE_RAD -DENABLE_TILES -DENABLE&
   &_PGC -DENABLE_QED -DENABLE_SHEAR -DENABLE_CYLMODES -DENABLE_QEDCYL -DEN&
   &ABLE_RADCYL -DENABLE_OVERDENSE -DENABLE_OVERDENSE_CYL -DENABLE_NEUTRAL_&
   &SPIN -DENABLE_XXFEL -DENABLE_GR -DENABLE_PKT",'
write(output_unit,'(A)') '    "F90"  :"/opt/homebrew//bin/mpif90 -Wa,-q",'
write(output_unit,'(A)') '    "F90F" :"-pipe -ffree-line-length-none -fno-r&
   &ange-check -lgcc_s.1.1 --openmp -g -Og -fbacktrace -fbounds-check -fimp&
   &licit-none -Wimplicit-interface -Wconversion  -Wsurprising -Wunderflow &
   & -I/opt/homebrew/Cellar/open-mpi/5.0.9/include -Wl,-flat_namespace -I/o&
   &pt/homebrew/Cellar/open-mpi/5.0.9/lib -I/opt/homebrew//include",'
write(output_unit,'(A)') '    "CF"   :"-Og -g -Wall -pedantic -mcpu=apple-m&
   &1 -std=c99  -DFORTRANSINGLEUNDERSCORE -DPRECISION_DOUBLE -D__MACH_TIMER&
   &__",'
write(output_unit,'(A)') '    "cc"   :"/opt/homebrew//bin/mpicc",'
write(output_unit,'(A)') '    "LDF"  :" -L/opt/homebrew//lib -lhdf5_fortran&
   & -lhdf5 -lz -lm -I/opt/homebrew/Cellar/open-mpi/5.0.9/include -Wl,-flat&
   &_namespace -I/opt/homebrew/Cellar/open-mpi/5.0.9/lib -L/opt/homebrew/Ce&
   &llar/open-mpi/5.0.9/lib -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_m&
   &pifh -lmpi  ",'
write(output_unit,'(A)') '   }'
write(output_unit,'(A)') '}'

