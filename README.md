Tools and techniques for using medical imaging software to visualize astronomy data
-----------------------------------------------------------------------------------

More information on the Astronomical Medicine project can be found [here](http://am.iic.harvard.edu/).

Unfortunately, the project is largely defunct at the moment, and you may have a very hard time compiling this code due to its
dependence on third-party libraries that have changed significantly since this code was originally developed.

I (Douglas Alan, doug AT alum.mit.edu) have two pre-built binaries, which you can get by contacting me. One of them will only run on rather
old versions of OS X. The other one is for Linux, but I do not know whether it will run properly on recent versions of Linux.

**NOTE:** The software developed specifically for Astronomical Medicine is distributed with the permissive MIT license,
but we use some third-party libraries that carry GPL licenses. Due to the viral nature of GPL, this means that any
executable generated from this code is bound by the terms of the GPL.

There are three branches in the respository:

**stable** -- This is the version that we released to the world. It is not known to have any huge glaring bugs, but it is missing some features that were developed later.

**internal-releases** -- This branch includes more recent versions that were used internally at the IIC and have some additional features, but these versions were never beat on quite thoroughly enough to release to the world.

**dev** -- This branch represents a significant rewrite of the code base to be much more elegant and maintainable, and to have some additional features. Unfortunately, the project ran out of funding before any scientists ever had a chance to use it and provide feedback, so it is probably not in a usable state. It does compile and passes some very minimal regression tests.
