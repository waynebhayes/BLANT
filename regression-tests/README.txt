To add a regression test, create a directory named to mkae it obvious what you're testing (eg CompressedFiles)
Then, put some example shell scripts in there with lines designed to test blant.
Any filename that ends in .sh (eg, regression-tests/you/test.sh) will be run nightly. It should exit(0) if the
test was successful, or some number > 0 if not. (This is standard Unix convention, 0 =success, nonzero=error.)

Take a look at the CompressedFiles/ directory for an example.

Be sure to 'chmod 755 *.sh' your scripts before doing a pull request.

