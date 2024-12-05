# Pedantic Makefile testing

All testing is done manually, and each test is accomplished in this sequence. Each command is run, and then it's expected behavior is written below.

## 1) `make pristine`

Removes all data from the canon_maps directory, except on data where `k=8`. If ran on an initially empty repository, it has no effect. This ensures that the repository is as close to a just-initialized state as possible, allowing the following tests to occur.

## 2) `make`

Makes the `base` target, which creates executables for libwayne, blant, and the canon map data. A lot of output is shown here, but end result should produce a full `canon_maps` directory.

## 3) `make`

Runs `make` again, but because everything has already been made and nothing has changed, no new files are generated and everything stays the same.

## 4) `touch ./src/fast-canon-map.c && make`

Simulates changing of the fast-canon-map code. Remakes canon map data. Since all canon maps depend on the `fast-canon-map`, which in turn depends on `./src/fast-canon-map.c` (or wherever the source directory is), this should regenerate all canon maps. In a real situation, if `fast-canon-map.c` was changed, then of course, all canon maps would also change. Note, however, that canon map data where `k=8` is not changed whatsoever, as directed.

## 5) `touch ./src/compute-alphas-NBE.c && make`

Simulates the changing of the compute-alphas-NBE code. Since all alpha list NBE data depends on `compute-alphas-NBE`, which in turn depends on `./src/compute-alphas-NBE.c`, this should regenerate all alpha list NBE files.

## 6) `touch ./src/compute-alphas-EBE.c && make`

Simulates the changing of the compute-alphas-EBE code, similar to (5).

## 7) `touch ./src/compute-alphas-MCMC.c && make`

Simulates the changing of the compute-alphas-MCMC code, similar to (5). However, these are simply just copied from `canon_maps.correct` due to their long compute time.

## 8) `touch canon_maps/canon_list7.txt && make`

Simulates the changing of canon list data, which should never happen, since these data files are generated. However, if it does happen somehow, all alpha list data targets (which depend on the canon_list data) are re-made.