# banks

This library provides an interface to evio banks.

To build, the env GEMC should point to the gemc directory containin "sources", and
the env BANKS should point to the location of this README.


To compile:

```
scons [-jN]  [OPT=1]
```

where N is the number of processors on your system (optional, for faster compilation) and OPT will optimize the build.


Two utilities are provided:

- bin/evioDump will print on screen the content of an evio file
- bin/evio2root will convert evio banks into ROOT trees. For this conversion you'll need the bank definitions
to be in the current directory (note: this should improve).


