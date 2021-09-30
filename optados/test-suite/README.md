 Optados test suite

## Dependencies

The code needs the `configparser` module, that can be installed e.g. via
`pip install --user configparser`.

## How to write a new test

### Writing a new test for optados.x 

1. Create a new folder for the test inside `test-suite/tests`, 
   with a short but meaningful name. 
   The name *must* start with the prefix `testopt_` if it is a Optados
   test. This is needed to properly group tests in categories.

2. modify the file `test-suite/tests/jobconfig` adding a new section, 
   following the template of existing tests. E.g.:
   ```
   # Testing preconditioner
   [testopt_precond_1]
   program = OPTADOS_ODO_OK
   inputs_args = ('gaas1.odi', '')
   output = gaas1.odo
   ```
   where:
   - Line 1: comment on what the test is supposed to do or test; 
     feel free to write a long description spanning multiple lines;
   - Line 2: test folder name created above, in square brackets
   - Line 3: name of a program that you want to run to test.
      The possible program names are defined in `test-suite/tests/userconfig`, discussed
      below. The "program" defines which executable to use to run the program, 
      which parser needs to be used, and which custom tolerances should be used to
      compare results of the test with reference results. It also defines if the test
      is expected to have the code fail (e.g. to check the code stops for unexpected 
      inputs).
   - Line 4: comma-separated list of tuples of length 2, containing the
     `('inputfile', 'cmdline-params')` for each subtest to run. Typically, you 
     will have only one subtest. The first string is the name of the input file within
     the test folder, the second are the command line parameters to pass to the 
     executable (use an empty string `''` if you don't have any custom parameter).
   - Line 5: defines the name of the file to parse.
  
   Additional parameters can be specified: check the 
   [testcode documentation](http://testcode.readthedocs.io/).

   One additional parameter is worth mentioning here:
   
   - `max_nprocs = 0`: this specifies the max number of CPUs that this test
     is allowed to run on. Zero means to never call the test with `mpirun`, 
     a larger number indicates the max number of MPI processes to use. Setting to
     zero is useful for tests that can only run in serial (e.g. for the gamma-only
     runs).

3. Put all needed input files in the folder (`.win`, `.amn`, `.mmn`, ...), making sure
   that the input file has the name you specified above in the `jobconfig` file.
   Add the files to the git repository. 
   Also, add a `.gitignore` file in your test folder for those files that are dynamically 
   created at runtime and should be ignored. The content of this file is also used by
   the `clean_tests` code (described later) to decide if a file can be safely deleted.

4. Compile the code in its most recent version.
5. Run the code (with the same parameters as defined in `jobconfig`) 
   in the test folder to create a reference output.

   `TODO: MAKE THIS AND THE NEXT POINT EASIER TO DO`

6. Open the file in point 5 above, verify that it is the actual expected output,
   and copy/move the output file to a file named `benchmark.out.default.inp=<inputfilename>` 
   (or `benchmark.out.default.inp=<inputfilename>.args=<cmdline_args>` if you have
   some command line parameters in your `input_args`). Add this to the git repo.

7. If you have chosen to use one of the existing "programs" already present in `userconfig`,
   you are probably already ok: just run the test again to check that the test now passes
   without errors. 

   **STRONG SUGGESTION**: To make sure the tests is working, try to change in the 
   reference/benchmark file one of the values that should be checked, to verify that 
   the tests actually fails if the value is unexpected. Remember to put back the correct
   value afterwards!

8. If instead you need to parse a different file, or parse additional data:
   * for an additional value, edit/improve the python parsing functions;
   * to parse a different file, create a new section in the `userconfig` file.
   In both cases, see the description below in the `userconfig` section.


# Parsing new files: the `userconfig` file
The `userconfig` file defines the "programs" to run. Each section
defines (at least) an executable to run and a function to parse the output.
Moreover, additional options can be provided like some custom tolerances.

An example program section looks like this:
```
[OPTADOS_ODO_OK]
exe = ../../optados.x
extract_fn = tools parsers.parse_odo.parse
tolerance = ( (1.0e-3, 5.0e-3, 'bandenergy'),
	             (1.0e-6, 1.0e-6, 'bandidx')))
```
Each line defines the following:
* Line 1: the name of the program in square brackets (that will be used in the `jobconfig` file)
  Try to comply to the following syntax: `<CODENAME>_<FILETOPARSE>_<SHOULDFAIL>` where:
  * `<CODENAME>` is `OPTADOS` 
  * `<FILETOPARSE>` is a short string defining the type of file that is expected to be
    parsed/checked (e.g. `ODO` or `DAT`)
  * `<SHOULDFAIL>` is `OK` if this is a standard run that should end with error code zero,
    or `FAIL` if you expect the code to fail.
* Line 2: specify the executable to run. Typically this is either `../../optados.x` 
   (the location is with respect to the folder in which `userconfig` is
  located)
* Line 3: define the (python) function to parse the output files. In the example above, 
  the two parameters indicate that the python modules live in the folder `tools` and, 
  within it, the function `parse` will be called, defined inside
  `parsers/parse_geninterp_dat.py`. 
  To parse a new file, define a new file inside `tools/parsers`, and within it define
  a `parse(fname)` python function that accepts a python function and returns
  a dictionary in the form `{'testvaluekey': [value1, value2, ...]}` for the values
  to test. We strongly suggest that you take inspiration from existing parsers (also
  for the logic to manage verbose output).
* Line 4-...: custom tolerances, using the format specified in the
  [testcode documentation](http://testcode.readthedocs.io/). In particular, for each
  value, the first number is a absolute tolarance, the second a relative tolerance, and
  the third is the `testvaluekey` returned by the parser.
  *Note*: if you don't specify a `testvaluekey` here, this is still checked with a fairly 
  strict tolerance. We still suggest to define explicitly all the `testvaluekey`s returned
  by the parser for clarity.

* the case `<SHOULDFAIL>=FAIL` is useful if you want to test an expected failure: e.g. 
  if you want to check that the code fails if you provide an unexpected input.
  In this case, you have to add an option to the section:
  ```
  can_fail = true
  ```
  to make sure that `testcode` does not mark the test as failed because the error code is
  non-zero.

Additional parameters can be specified: check the 
[testcode documentation](http://testcode.readthedocs.io/).

# How to run the tests
In the `test-suite` folder, run `./run_tests`. It will prompt you for the test
you want to run, and then run them.
The code has a number of command-line options to run in non-interactive mode, or
to specify some options (number of MPI processors for parallel runs, verbose mode).
Run `./run_tests -h` for further info.

**Note**: you will need to have compiled `optados.x` and to be able
to run the tests.

## Cleaning up
When running tests, a number of temporary files are created. While these should not
create problems when the test is run multiple times, sometimes it is better to remove
these files.
To do this, run `./clean_tests` in the `test-code` folder. This will delete files
that are for sure an output of the test code. If you want a deeper clean-up 
of files that are ignored by git (see also the output of the code for further info), 
you can run `./clean_tests -i`. This is typically safe, but double check (especially if
you are creating a new test, to avoid the unexpected deletion of files).
For this reason, we strongly suggest that you commit a `.gitignore` file
inside your testfolder for files you know might be generated during the run.

# Acknowledgements
The core of the test suite uses `testcode` by J. Spencer, hosted
on [this GitHub repository](https://github.com/jsspencer/testcode).

We also acknowledge S. Poncé for the first implementation of the test-suite
in Wannier90. G. Pizzi wrote the current Wannier90 test-suite. Cloned with
permission for Optados in Sept 2019.
