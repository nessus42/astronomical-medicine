# -*- Mode: CMake; fill-column: 79 -*-

add_library(libwcs
  imsetwcs.c imgetwcs.c matchstar.c findstar.c daoread.c wcscon.c
  fitswcs.c wcsinit.c wcs.c ty2read.c webread.c tmcread.c
  gscread.c gsc2read.c ujcread.c uacread.c ubcread.c ucacread.c
  sdssread.c tabread.c binread.c ctgread.c actread.c catutil.c
  imrotate.c fitsfile.c imhfile.c
  hget.c hput.c imio.c dateutil.c imutil.c
  worldpos.c tnxpos.c dsspos.c platepos.c
  sortstar.c platefit.c iget.c fileutil.c polfit.c
  wcslib.c lin.c cel.c proj.c sph.c wcstrig.c distort.c
  )
