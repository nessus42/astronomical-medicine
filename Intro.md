# Introduction #

The Astronomical Medicine source code repository has three branches in it. These branches are visible via the `hg branches` command. To check out a specific branch from a Mecurial clone into the working directory, you would use the following Mercurial command:
```
hg co stable
```
for instance, if you want the "stable" branch. If, on the other hand, you want the "dev" branch, you would replace "stable" with "dev" in the above Mercurial command.

There are three branches at the moment:

  * **stable** -- This is the version that we released to the world. It is not known to have any huge glaring bugs, but it is missing some features that were developed later.

  * **internal-releases** -- This branch includes more recent versions that were used internally at the IIC and have some additional features, but these versions were never beat on quite thoroughly enough to release to the world.

  * **dev** -- This branch represents a significant rewrite of the code base to be much more elegant and maintainable, and to have some additional features. Unfortunately, the project ran out of funding before any scientists ever had a chance to use it and provide feedback, so it is probably not in a usable state. It _does_ compile and passes some very minimal regression tests.

# More information #

More information on the Astronomical Medicine project can be found [here](http://am.iic.harvard.edu).
