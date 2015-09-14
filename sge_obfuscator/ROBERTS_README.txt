Setup and run the main program just by following the advice in the other readme files (the originals from Malozemoff)

BUT FIRST:
-In code/src, run: g++ -o zen2 zencode.cpp circuit.cpp clt_mlm.cpp utils.cpp -lgmpxx -lgomp -lgmp -fopenmp

The -fopenmp part is optional. It may cause the program not to work.

-In code/src/launch_zen.sh, replace the filepath with whatever the absolute filepath of zen1 is. Then, give launch_zen.sh permission to run.

-Before running code, run the commands in code/setenv.sh. Don't actually run the file, it does nothing. Run the commands in it by copy+pasting to terminal.

In code/src/_zobfuscator.cpp, in line 252 (or search for "launch"), replace the filepath with the absolute filepath of launch_zen.sh.

When running, can specify number of slots to use with --nslots

