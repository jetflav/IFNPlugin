# To Do

## if the code moves to fjcontrib

- make sure data/... -> ../data/...
- make check script in Makefile points to fjcontrib one
- decide what to do about the -std compile flag

## General things

- make sure MassFlavRecombiner also gets the relevant arguments?

- integrate that in terms of the modulo_2 handling elsewhere in the code
  (so that all previous modulo_2 functionality goes into FlavRecombiner?
  Watch out if that complicates the interface and potentially retain
  the current modulo_2 arguments)